import cooler
import bioframe
import pandas as pd
import click
import os
import sc3dg.cooltools.cooltools as cooltools
from joblib import Parallel, delayed

def analyze_chromatin_interactions(resolution, mcool_file_path, chromsizes_file, centromere_file, output_csv, nthread=4):
    """
    分析染色体相互作用数据，计算预期的cis相互作用，并识别'dots'。

    参数:
    resolution: int - 分辨率，用于确定'dots'的清晰度。
    mcool_file_path: str - mcool文件的路径。
    chromsizes_file: str - 染色体大小文件的路径。
    centromere_file: str - 中心粒位置文件的路径。
    output_csv: str - 输出结果CSV文件的路径。
    """

    # 打开mcool文件
    clr = cooler.Cooler(mcool_file_path)

    # 获取染色体大小信息
    hg38_chromsizes = pd.read_csv(chromsizes_file, sep="\t", header=None)
    hg38_chromsizes.columns = ['name', 'length']
    hg38_chromsizes = hg38_chromsizes.set_index("name")["length"]
    # 读取中心粒位置信息
    hg38_cens = pd.read_csv(centromere_file, sep="\t")

    # 制作染色体臂信息
    hg38_arms = bioframe.make_chromarms(hg38_chromsizes, hg38_cens)

    # 选择只包含在cooler中的染色体
    hg38_arms = hg38_arms.set_index("chrom").loc[clr.chromnames].reset_index()

    # 计算预期的cis相互作用
    expected = cooltools.expected_cis(
        clr,
        view_df=hg38_arms,
        nproc=nthread,
    )

    # 识别'dots'
    dots_df = cooltools.dots(
        clr,
        expected=expected,
        view_df=hg38_arms,
        max_loci_separation=10000000,  # 指定主对角线两侧的分离距离
        nproc=nthread,
    )

    # 将结果保存到CSV文件
    dots_df.to_csv(output_csv)

# # 使用示例
# resolution = 10000
# mcool_file_path = '/Users/jiangwenjie/D/ScHIC_pipeline/scHIC_process/CallLoop/Brain_HIC_mm10.mcool::/resolutions/'+str(resolution)
# chromsizes_file = '/Users/jiangwenjie/D/ScHIC_pipeline/scHIC_process/CallLoop/mm10.chromsize'
# centromere_file = '/Users/jiangwenjie/D/ScHIC_pipeline/scHIC_process/CallLoop/mm10_centromere.txt'
# output_csv = '/Users/jiangwenjie/D/ScHIC_pipeline/scHIC_process/CallLoop/Loop_summary.csv'


@click.command('loop', short_help='detect loop')
@click.option('--resolution', help='resolution', type=int)
@click.option('--mcool', help='mcool file path')
@click.option('--chromsize', help='chromsizes file path')
@click.option('--centromere', help='centromere file path')
@click.option('--output', help='output csv file path')
@click.option('--nthread', help='number of threads', default=4)
def loop(resolution, mcool, chromsize, centromere, output, nthread):
    task = []
    for r, d, fs in os.walk(mcool):
        for f in fs:
            if f.endswith('.mcool'):
                mcc = mcool + '/' + f + '::/resolutions/' + str(resolution)
                output = output + '/' + f.split('.')[0] + '.csv'
                task.append((resolution, mcc, chromsize, centromere, output, nthread))
            

    print(task)

    Parallel(n_jobs=nthread)(delayed(analyze_chromatin_interactions)(*t) for t in task)
    

