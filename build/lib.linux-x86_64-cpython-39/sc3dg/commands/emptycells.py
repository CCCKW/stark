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

    mat = pd.DataFrame()
    if not support_chrom:
        chrom = ['chr' + str(i) for i in range(1, 22)] + ['chrX', 'chrY']
    else:
        chrom = support_chrom.split(',')
    tmp_chrom = []
    for i in range(len(chrom)):
        for j in range(i, len(chrom)):
            tmp_chrom.append(chrom[i] + '-' + chrom[j])

    mat['type'] = tmp_chrom
    mat = mat.set_index('type')
    c = 0

    cell_name = []
    for val in os.listdir(mcool_path + '/'):
        if val.endswith('mcool'):
            cell_name.append(val)
            cr = cooler.Cooler(mcool_path + '/' + val + '::resolutions/1000000')
            bins = cr.bins()[:]
            pixels = cr.pixels()[:]
            # print(bins)
            pixels['chr1'] = list(bins.iloc[list(pixels['bin1_id']), 0])
            pixels['chr2'] = list(bins.iloc[list(pixels['bin2_id']), 0])
            pixels['type'] = pixels['chr1'] + '-' + pixels['chr2']
            pixels = pixels.loc[:, ['type', 'count']]

            pixels = pixels.groupby('type').sum('count')
            mat = pd.concat([mat, pixels], axis=1)
            c += 1
            if c % 1000 == 0:
                print('Dealing ', c)
    mat = mat.fillna(0)

    mat.columns = cell_name
    mat.to_csv(output + '/matrix.csv')
    umi = 500
    matrix = CountMatrix.from_dataFrame(mat)
    orig_cell_bcs = mat.columns.values.astype(bytes)
    eval_bcs, log_likelihood, pvalues, pvalues_adj, is_nonambient = find_nonambient_barcodes(
        matrix,  # Full expression matrix
        orig_cell_bcs,  # (iterable of str): Strings of initially-called cell barcodes
        min_umi_frac_of_median=0.01,
        min_umis_nonambient=umi,
        max_adj_pvalue=0.01
    )
    new_mat = pd.DataFrame(mat.sum(0))
    new_mat['ambient'] = True
    new_mat['ambient'][list(eval_bcs[list(is_nonambient)])] = False
    new_mat.columns = ['count', 'ambient']
    new_mat = new_mat.sort_values('count', ascending=False)
    new_mat['rank'] = list(range(new_mat.shape[0]))
    print('#######shape#######', new_mat.shape[0])
    new_mat.to_csv('umi_%s_%s.csv' % (umi, tag))
    new_mat[new_mat['ambient'] == False]
    res[umi] = new_mat[new_mat['ambient'] == False].shape[0]
    # 为了模仿上面的图表，我们需要创建一个对数刻度的散点图。
    # 由于我们没有相同的数据，我们将使用之前创建的数据并进行一些调整。
    import matplotlib.pyplot as plt
    df = new_mat
    # 创建一个对数刻度的散点图
    plt.figure(figsize=(6, 5), dpi=200)

    # 我们将使用对数刻度来更好地模仿给定的图表
    plt.xscale('log')

    # 绘制散点图，使用不同的颜色区分ambient值
    plt.plot(df[df['ambient']]['rank'], df[df['ambient']]['count'], color='grey', alpha=0.7, label='Background')
    plt.plot(df[~df['ambient']]['rank'], df[~df['ambient']]['count'], color='blue', alpha=0.7, label='Cells')

    # 设置坐标轴的标题和图表的标题
    plt.xlabel('Barcodes')
    plt.ylabel('UMI counts')
    plt.title('Barcode Rank Plot')

    # 添加图例
    plt.legend()

    # 由于原始图像的y轴是对数刻度，我们也设置y轴为对数刻度
    plt.yscale('log')

    # 设置y轴的范围，使其更加匹配原始图表的范围
    plt.ylim(bottom=1)

    plt.savefig(output + '/rankplit.png' % (umi, tag), dpi=200)
    # 显示图表
    # plt.show()
    pd.DataFrame(res, index=[0]).to_csv(output + '/res.csv')