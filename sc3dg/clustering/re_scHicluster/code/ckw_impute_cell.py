import os
import time
import h5py
import cv2
import json
cv2.useOptimized()
import argparse
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix, save_npz, diags, eye
from scipy.sparse.linalg import norm
from scipy import ndimage as nd

def impute_cell(config, 
                cell,
                chrom,
                logscale=False, pad=1, std=1, rp=0.5, tol=0.01, window_size=500000000, step_size=10000000,
                output_dist=500000000, output_format='hdf5', mode=None):
    
    config = json.load(open(config))
    rp = config.get('rp', rp)
    indir = config['temp_dir'] + '/cell_matrix'
    outdir = config['temp_dir'] + '/imputed_matrix'
    res = int(config['res'])
    chrom_file = config['genome_reference_path']
    
    chromsize = pd.read_csv(chrom_file, sep='\t', header=None, index_col=0).to_dict()[1]
    
    
    def random_walk_cpu(P, rp, tol):
        if rp==1:
            return P
        ngene = P.shape[0]
        I = eye(ngene)
        Q = P.copy()
        start_time = time.time()
        for i in range(30):
            Q_new = P.dot(Q * (1 - rp) + rp * I)
            delta = norm(Q - Q_new)
            Q = Q_new.copy()
            sparsity = Q.nnz / ngene / ngene
            end_time = time.time()
            # print('Iter', i+1, 'takes', end_time-start_time, 'seconds, loss', delta, 'sparsity', sparsity)
            if delta < tol:
                break
        return Q

    def output():
        if output_format=='npz':
            os.makedirs(f'{outdir}/{cell}', exist_ok=True)
            save_npz(f'{outdir}/{cell}/{cell}_{c}_{mode}.npz', E.astype(np.float32))
        else:
            os.makedirs(f'{outdir}/{cell}', exist_ok=True)
            f = h5py.File(f'{outdir}/{cell}/{cell}_{c}_{mode}.hdf5', 'w')
            g = f.create_group('Matrix')
            g.create_dataset('data', data=E.data, dtype='float32', compression='gzip')
            g.create_dataset('indices', data=E.indices, dtype=int, compression='gzip') 
            g.create_dataset('indptr', data=E.indptr, dtype=int, compression='gzip')
            g.attrs['shape'] = E.shape
            f.close()
        return

    if not os.path.exists(outdir):
        os.makedirs(outdir, exist_ok=True)
     
    if chrom[:3]=='chr':
        c = chrom
    else:
        c = 'chr' + chrom
    if not mode:
        mode = f'pad{str(pad)}_std{str(std)}_rp{str(rp)}_sqrtvc'

    ws = window_size // res
    ss = step_size // res

    chromsize = pd.read_csv(chrom_file, sep='\t', header=None, index_col=0).to_dict()[1]

    start_time = time.time()
    ngene = int(chromsize[c] // res) + 1
    D = np.loadtxt(f'{indir}/{c}/{cell}_{c}.txt')
    
    # 处理空数据
    if len(D)==0:
        E = eye(ngene).tocsr()
        output()
        return

    elif len(D.shape)==1:
        D = D.reshape(1,-1)

    # 简单粗暴的解决方案：直接过滤掉超出范围的索引
    original_count = len(D)
    valid_mask = (D[:, 0] < ngene) & (D[:, 1] < ngene) & (D[:, 0] >= 0) & (D[:, 1] >= 0)
    D = D[valid_mask]
    
    # 输出过滤信息
    filtered_count = original_count - len(D)
    # if filtered_count > 0:
    #     print(f"{cell}_{c}: 过滤掉 {filtered_count} 个超出范围的数据点 (总数: {original_count})")
    
    # 如果所有数据都被过滤掉了，创建空矩阵
    if len(D) == 0:
        print(f"{cell}_{c}: 所有数据都被过滤掉了，创建空矩阵")
        E = eye(ngene).tocsr()
        output()
        return

    # 确保索引是整数类型
    D = D.astype(int)
    
    A = csr_matrix((D[:, 2], (D[:, 0], D[:, 1])), shape = (ngene, ngene))
    if logscale:
        A.data = np.log2(A.data + 1)

    end_time = time.time()
    # print('Loading takes', end_time-start_time, 'seconds')

    start_time = time.time()
    if pad > 0:
        A = cv2.GaussianBlur((A + A.T).astype(np.float32).toarray(), (pad*2+1, pad*2+1), std)
    else:
        A = (A + A.T).astype(np.float32)

    end_time = time.time()
    # print('Convolution takes', end_time-start_time, 'seconds')

    start_time = time.time()
    # remove diagonal before rwr
    A = csr_matrix(A)
    A = A - diags(A.diagonal())
    if ws>=ngene or rp==1:
        B = A + diags((A.sum(axis=0).A.ravel()==0).astype(int))
        d = diags(1 / B.sum(axis=0).A.ravel())
        P = d.dot(B)
        E = random_walk_cpu(P, rp, tol)
        E = E.multiply(E > 0.01/ngene)
    else:
        idx = (np.repeat(np.arange(ws), ws), np.tile(np.arange(ws), ws))
        idxfilter = (np.abs(idx[1] - idx[0]) < (output_dist // res + 1))
        idx = (idx[0][idxfilter], idx[1][idxfilter])
        # first filter
        idxfilter = ((idx[0] + idx[1]) < (ws + ss))
        idx1 = (idx[0][idxfilter], idx[1][idxfilter])
        mask1 = csr_matrix((np.ones(len(idx1[0])), (idx1[0], idx1[1])), (ws, ws))
        # last filter
        idxfilter = ((idx[0] + idx[1]) >= ((ngene-ws) // ss * 2 + 1) * ss + 3 * ws - 2 * ngene)
        idx2 = (idx[0][idxfilter], idx[1][idxfilter])
        mask2 = csr_matrix((np.ones(len(idx2[0])), (idx2[0], idx2[1])), (ws, ws))
        # center filter
        idxfilter = np.logical_and((idx[0] + idx[1]) < (ws + ss), (idx[0] + idx[1]) >= (ws - ss))
        idx0 = (idx[0][idxfilter], idx[1][idxfilter])
        mask0 = csr_matrix((np.ones(len(idx0[0])), (idx0[0], idx0[1])), (ws, ws))
        start_time = time.time()
        E = csr_matrix(A.shape)
        for ll in [x for x in range(0, ngene - ws, ss)] + [ngene - ws]:
            B = A[ll:(ll+ws), ll:(ll+ws)]
            B = B + diags((B.sum(axis=0).A.ravel()==0).astype(int))
            d = diags(1 / B.sum(axis=0).A.ravel())
            P = d.dot(B)
            Etmp = random_walk_cpu(P, rp, tol)
            if ll==0:
                E[ll:(ll+ws), ll:(ll+ws)] += Etmp.multiply(mask1)
            elif ll==(ngene-ws):
                E[ll:(ll+ws), ll:(ll+ws)] += Etmp.multiply(mask2)
            else:
                E[ll:(ll+ws), ll:(ll+ws)] += Etmp.multiply(mask0)
            print('Window', ll)
        E = E.multiply(E > 0.01/ws)
    # print('RWR takes', time.time() - start_time, 'seconds')

    start_time = time.time()
    E = E + E.T
    d = E.sum(axis=0).A.ravel()
    d[d==0] = 1
    b = diags(1 / np.sqrt(d))
    E = b.dot(E).dot(b)
    # print('SQRTVC takes', time.time() - start_time, 'seconds')

    # longest distance filter mask
    start_time = time.time()
    idx = np.triu_indices(E.shape[0], 0)
    if (output_dist // res + 1) < ngene:
        idxfilter = ((idx[1] - idx[0]) < (output_dist // res + 1))
        idx = (idx[0][idxfilter], idx[1][idxfilter])
    mask = csr_matrix((np.ones(len(idx[0])), (idx[0], idx[1])), E.shape)
    E = E.tocsr().multiply(mask)
    # print('Filter takes', time.time() - start_time, 'seconds')

    output()
    return

if __name__ == '__main__':
    impute_cell(
        config='config.JSON',
        cell='GSM7682104_BJ',
        chrom='chr3'
    )