from sc3dg.utils.cell_calling import find_nonambient_barcodes
from sc3dg.utils.matrix import CountMatrix
import numpy as np
import pandas as pd
import cooler
import os
import scanpy as sc
import click


@click.command('emptycells', short_help='emptycells for QC')
@click.option('--mcool-path', help='mcool file path', type=str)
@click.option('--support_chrom', help='support chrom file path', type=str, default=None)
@click.option('--output', help='output csv file path', type=str)
def emptycells(mcool_path, support_chrom, output):
    os.makedirs(output, exist_ok=True)
    res = {}
    
    if os.path.exists(output + '/umi_1000_parallel.csv'):
        print('Results already exist, skipping computation.')
        return

    # 设置进程数，默认使用CPU核心数
    if n_processes is None:
        n_processes = min(cpu_count(), 20)  # 限制最大进程数避免内存溢出

    mat = pd.DataFrame()
    if not support_chrom:
        chrom = ['chr' + str(i) for i in range(1, 22)] + ['chrX', 'chrY']
    else:
        chrom = support_chrom.split(',') if isinstance(support_chrom, str) else support_chrom
    
    tmp_chrom = []
    for i in range(len(chrom)):
        for j in range(i, len(chrom)):
            tmp_chrom.append(chrom[i] + '-' + chrom[j])

    mat['type'] = tmp_chrom
    mat = mat.set_index('type')

    # 获取所有mcool文件
    cell_name = []
    mcool_files = []
    for val in os.listdir(mcool_path + '/'):
        if val.endswith('mcool'):
            cell_name.append(val)
            mcool_files.append((mcool_path, val, chrom))

    print(f'Found {len(mcool_files)} mcool files, using {n_processes} processes')

    # 并行处理mcool文件
    with Pool(processes=n_processes) as pool:
        results = []
        for i, result in enumerate(pool.imap(process_single_mcool, mcool_files)):
            val, pixels = result
            mat = pd.concat([mat, pixels], axis=1)
            if (i + 1) % 100 == 0:
                print(f'Processed {i + 1}/{len(mcool_files)} files')
        
    mat = mat.fillna(0)
    mat.columns = cell_name
    mat.to_csv(output + '/matrix.csv')
    
    # 构建cell calling的配置
    cell_calling_config = {
        'MIN_UMIS': config.get('MIN_UMIS', 1000),
        'MIN_UMI_FRAC_OF_MEDIAN': config.get('MIN_UMI_FRAC_OF_MEDIAN', 0.01),
        'MAX_ADJ_PVALUE': config.get('MAX_ADJ_PVALUE', 0.01),
        'N_CANDIDATE_BARCODES': config.get('N_CANDIDATE_BARCODES', 10000),
        'NUM_SIMS': config.get('NUM_SIMS', 2000),
        'N_CORES': config.get('N_CORES', min(cpu_count(), 20)),
        'MAX_MEM_GB': config.get('MAX_MEM_GB', 0.1),
        'EMPTY_BC_FRAC': config.get('EMPTY_BC_FRAC', 0.03),
    }
    
    matrix = CountMatrix.from_dataFrame(mat)
    orig_cell_bcs = mat.columns.values.astype(bytes)
    
    result = find_nonambient_barcodes_ultra_fast(
        matrix,
        orig_cell_bcs,
        custom_config=cell_calling_config
    )
    
    if result is None:
        print("❌ Cell calling failed")
        return
    
    eval_bcs, log_likelihood, pvalues, pvalues_adj, is_nonambient = result
    
    new_mat = pd.DataFrame(mat.sum(0))
    new_mat['ambient'] = True
    new_mat['ambient'].iloc[eval_bcs[is_nonambient]] = False
    new_mat.columns = ['count', 'ambient']
    new_mat = new_mat.sort_values('count', ascending=False)
    new_mat['rank'] = list(range(new_mat.shape[0]))
    print('#######shape#######', new_mat.shape[0])
    
    # 使用config中的参数生成文件名
    tag = config.get('tag', 'parallel')
    new_mat.to_csv(output + f'/umi_{umi}_{tag}.csv')
    new_mat_filtered = new_mat[new_mat['ambient'] == False]
    res[umi] = new_mat_filtered.shape[0]
    
    import matplotlib.pyplot as plt
    df = new_mat
    plt.figure(figsize=(6, 5), dpi=200)
    plt.xscale('log')
    
    plt.plot(df[df['ambient']]['rank'], df[df['ambient']]['count'], color='grey', alpha=0.7, label='Background')
    plt.plot(df[~df['ambient']]['rank'], df[~df['ambient']]['count'], color='blue', alpha=0.7, label='Cells')
    
    plt.xlabel('Barcodes')
    plt.ylabel('UMI counts')
    plt.title('Barcode Rank Plot')
    plt.legend()
    plt.yscale('log')
    plt.ylim(bottom=1)
    
    plt.savefig(output + f'/rankplot_{umi}_{tag}.png', dpi=200)
    pd.DataFrame(res, index=[0]).to_csv(output + '/res.csv')
