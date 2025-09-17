# command time python /gale/ddn/snm3C/humanPFC/code/embed_concatcell_chr.py 
# --cell_list /gale/ddn/snm3C/humanPFC/smoothed_matrix/filelist/imputelist_pad1_std1_rp0.5_sqrtvc_chr${c}.txt 
# --outprefix /gale/ddn/snm3C/humanPFC/smoothed_matrix/${res0}b_resolution/merged/pad1_std1_rp0.5_sqrtvc_chr${c} 
# --res ${res}
import os
import time
import h5py
import argparse
import numpy as np
from scipy.sparse import load_npz, save_npz, csr_matrix, vstack
from sklearn.decomposition import TruncatedSVD

def embed_concatcell_chr(config, celllist, outprefix, dist=10000000, save_raw=True, dim=50):
     
    res = config['res']


    with h5py.File(celllist[0], 'r') as f:
        ngene = f['Matrix'].attrs['shape'][0]
    idx = np.triu_indices(ngene, k=1)
    idxfilter = np.array([(yy - xx) < (dist / res + 1) for xx,yy in zip(idx[0], idx[1])])
    idx = (idx[0][idxfilter], idx[1][idxfilter])

    if os.path.exists(f'{outprefix}.svd{dim}.npy'):
        print(f'{outprefix}.svd{dim}.npy already exists, skip embedding')
        return
    
    start_time = time.time()
    # matrix = np.zeros((len(celllist), np.sum(idxfilter)))
    matrix = []
    for i,cell in enumerate(celllist):
        with h5py.File(cell, 'r') as f:
            g = f['Matrix']
            A = csr_matrix((g['data'][()], g['indices'][()], g['indptr'][()]), g.attrs['shape'])
        # matrix[i] = A[idx]
        matrix.append(csr_matrix(A[idx]))
        # if i%100==0:
        #     print(i, 'cells loaded', time.time() - start_time, 'seconds')

    matrix = vstack(matrix)

    if save_raw:
        save_npz(f'{outprefix}.npz', matrix)

    scalefactor = 100000
    matrix.data = matrix.data * scalefactor
    svd = TruncatedSVD(n_components=dim, algorithm='arpack')
    matrix_reduce = svd.fit_transform(matrix)
    matrix_reduce = np.concatenate((svd.singular_values_[None,:], matrix_reduce), axis=0)
    np.save(f'{outprefix}.svd{dim}.npy', matrix_reduce)
    return


