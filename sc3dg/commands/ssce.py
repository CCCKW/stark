import os
import cooler
import sc3dg.cooltools.cooltools as cooltools
import numpy as np
import pandas as pd
from sc3dg.cooltools.cooltools import insulation
import bioframe
import sys
import statsmodels.api as sm
# linear
from sklearn.linear_model import Ridge
from joblib import Parallel, delayed
import click

pre_chrom = ['chr' + str((i + 1)) for i in range(22)] + ['chrX', 'chrY']


def rg(prepare):
    prepare = prepare.loc[~(prepare == 0).all(axis=1)]
    prepare = np.log2(prepare + 1)
    X = prepare[["dup", "trans", "tad", "ab"]]
    y = prepare["contacts"]
    # 使用statsmodels进行线性回归
    # 使用statsmodels进行线性回归
    X = sm.add_constant(X)

    # 使用岭回归进行拟合
    ridge = Ridge(alpha=1.0, fit_intercept=False)  # 设置fit_intercept=False,让Ridge不再估计截距项
    ridge.fit(X, y)

    # 提取标准化系数的绝对值
    beta_abs = np.abs(ridge.coef_)

    # 根据变量的相关性分为正负两组,并将常数项视为负相关指标
    pos_vars = ["tad", "ab"]
    neg_vars = ["dup", "trans", "const"]

    # 计算正负变量标准化系数的和
    pos_sum = beta_abs[X.columns.isin(pos_vars)].sum()
    neg_sum = beta_abs[X.columns.isin(neg_vars)].sum()

    # 计算有效信息指标
    valid_info = pos_sum / neg_sum

    return valid_info




def rg(prepare):
    prepare = prepare.loc[~(prepare == 0).all(axis=1)]
    prepare = np.log2(prepare + 1)
    X = prepare[["trans", "tad", "ab"]]
    y = prepare["contacts"]
    # 使用statsmodels进行线性回归
    # 使用statsmodels进行线性回归
    X = sm.add_constant(X)

    # 使用岭回归进行拟合
    ridge = Ridge(alpha=1.0, fit_intercept=False)  # 设置fit_intercept=False,让Ridge不再估计截距项
    ridge.fit(X, y)

    # 提取标准化系数的绝对值
    beta_abs = np.abs(ridge.coef_)

    # 根据变量的相关性分为正负两组,并将常数项视为负相关指标
    pos_vars = ["tad", "ab"]
    neg_vars = ["trans", "const"]

    # 计算正负变量标准化系数的和
    pos_sum = beta_abs[X.columns.isin(pos_vars)].sum()
    neg_sum = beta_abs[X.columns.isin(neg_vars)].sum()

    # 计算有效信息指标
    valid_info = pos_sum / neg_sum

    return valid_info


def deal(mcool, g, name, output):
    cache_path = './cache'
    os.makedirs(cache_path, exist_ok=True)
    if os.path.exists(output + '/' + name + '.txt'):
        return
    print(mcool)

    try:
        clr = cooler.Cooler('%s::/resolutions/40000' % mcool)
    except:
        print('err load, re-zoomify')
        cooler.zoomify_cooler(base_uris='%s::/resolutions/10000' % mcool,
                              outfile=cache_path + '/' + os.path.basename(mcool),
                              resolutions=[10000, 40000, 500000, 1000000],
                              chunksize=100000000, nproc=10)
        os.system('mv %s %s' % (cache_path + '/' + os.path.basename(mcool), mcool))
        clr = cooler.Cooler('%s::/resolutions/40000' % mcool)
    if not cooltools.lib.is_cooler_balanced(clr):
        cooler.balance_cooler(clr, store=True)
    win = [10 * 40000]
    insulation_table = insulation(clr, win, verbose=False)
    insulation_table = insulation_table.groupby('chrom')['is_boundary_400000'].sum()
    insulation_table = pd.DataFrame(insulation_table[insulation_table.index.isin(pre_chrom)])
    insulation_table.columns = ['tad']
    print('done IS')
    clr = cooler.Cooler('%s::/resolutions/1000000' % mcool)
    if not cooltools.lib.is_cooler_balanced(clr):
        cooler.balance_cooler(clr, store=True)
    genome = bioframe.load_fasta(g)
    bins = clr.bins()[:]
    gc_cov = bioframe.frac_gc(bins[['chrom', 'start', 'end']], genome)
    view_df = pd.DataFrame({'chrom': clr.chromnames,
                            'start': 0,
                            'end': clr.chromsizes.values,
                            'name': clr.chromnames}
                           )

    try:
        cis_eigs = cooltools.eigs_cis(
            clr,
            gc_cov,
            view_df=view_df,
            n_eigs=3,
        )

        # cis_eigs[0] returns eigenvalues, here we focus on eigenvectors
        eigenvector_track = cis_eigs[1][['chrom', 'start', 'end', 'E1']]

        eigenvector_track['E1'] = np.nan_to_num(eigenvector_track['E1'])
        eigenvector_track['sign'] = '$'
        eigenvector_track['sign'][eigenvector_track['E1'] > 0] = 'A'
        eigenvector_track['sign'][eigenvector_track['E1'] < 0] = 'B'
        a_list = []
        b_list = []
        ab_switch = {}
        for c in pre_chrom:
            ab_switch[c] = 0
        for c in pre_chrom:
            if not c in set(eigenvector_track['chrom']):
                continue
            seq = list(eigenvector_track[eigenvector_track['chrom'] == c]['sign'])
            base = eigenvector_track[eigenvector_track['chrom'] == c].index[0]

            seq = ''.join(seq)
            #     print(seq)
            new_seq = ''
            idx = []
            tmp_idx = []
            for i in range(len(seq)):
                sym = seq[i]
                if new_seq == '':
                    new_seq += sym
                    tmp_idx.append(i + base)
                else:
                    if sym == seq[i - 1]:
                        tmp_idx.append(i + base)
                    else:
                        idx.append(tmp_idx)
                        tmp_idx = [i + base]
                        new_seq += sym
            for i in range(len(new_seq) - 3):
                if new_seq[i:i + 3] == 'ABA':
                    tmp_A = idx[i] + idx[i + 2]
                    tmp_B = idx[i + 1]
                    ab_switch[c] += 1

                elif new_seq[i:i + 3] == 'BAB':
                    tmp_A = idx[i + 1]
                    tmp_B = idx[i] + idx[i + 2]
                    ab_switch[c] += 1
        ab_switch = pd.DataFrame(ab_switch, index=[0]).T
        ab_switch.columns = ['ab']
        print(ab_switch)
    except:
        ab_switch = pd.DataFrame({c: 0 for c in pre_chrom}, index=[0]).T
        ab_switch.columns = ['ab']

    clr = cooler.Cooler('%s::/resolutions/40000' % mcool)
    bins = clr.bins()[:]
    pixels = clr.pixels()[:]
    pixels['chr1'] = list(bins.iloc[list(pixels['bin1_id']), 0])
    pixels['chr2'] = list(bins.iloc[list(pixels['bin2_id']), 0])
    pixels['start'] = list(bins.iloc[list(pixels['bin1_id']), 1])
    pixels['end'] = list(bins.iloc[list(pixels['bin1_id']), 2])
    pixels = pixels[pixels['chr1'].isin(pre_chrom) & pixels['chr2'].isin(pre_chrom)]
    trans_result = pixels[pixels['chr1'] != pixels['chr2']]
    trans1 = pd.DataFrame(trans_result.groupby('chr1')['count'].sum())
    trans2 = pd.DataFrame(trans_result.groupby('chr2')['count'].sum())
    trans_result = pd.concat([trans1, trans2], axis=1)
    trans_result.fillna(0, inplace=True)
    trans_result = pd.DataFrame(trans_result.sum(1))
    trans_result.columns = ['trans']
    intra = pixels[pixels['chr1'] == pixels['chr2']]
    intra = pd.DataFrame(intra.groupby('chr1')['count'].sum())

    raw = {}
    for c in set(pixels['chr1']):
        raw[c] = pixels[(pixels['chr1'] == c) | (pixels['chr2'] == c)]['count'].sum()

    raw = pd.DataFrame(raw, index=[0]).T
    raw.columns = ['contacts']

    result = pd.concat([raw, trans_result, ab_switch, insulation_table], axis=1)
    result.fillna(0, inplace=True)
    result.to_csv(output + '/' + name + '.csv')
    result = pd.read_csv(output + '/' + name + '.csv')
    result = result.set_index('Unnamed: 0')
    try:
        valid = rg(result)
    except:
        valid = 0

    with open(output + '/' + name + '.txt', 'w+') as f:
        f.write(str(valid))
    f.close()

    return valid



@click.command('ssce', short_help='cal Spatial Structure Capture Efficiency')
@click.option('--mcool', help='mcool file path')
@click.option('--genome', help='genome file path')
@click.option('--output', help='save path')
@click.option('--nproc', help='number of process',default=5)
def ssce(mcool, genome,  output, nproc):
    os.makedirs(output, exist_ok=True)
    task = []
    for mc in os.listdir(mcool):
        if mc.endswith('.mcool'):
            name = mc.split('.')[0]
            mc = mcool + '/' +mc
            task.append((mc, genome, name, output))
    print(task)
    Parallel(n_jobs=nproc)(delayed(deal)(*t) for t in task)
