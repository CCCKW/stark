import os
import json
import shutil

data_path = '/cluster2/home/Kangwen/0_benchmark/0_benchmark_dataset'
template_path = '/cluster2/home/Kangwen/0_benchmark/1_method/re_higashi/template.JSON'

reslutions_list = ['10k', '50k','100k', '500k', '1Mb']


result_path = '/cluster2/home/Kangwen/0_benchmark/2_result/0_higashi'
os.makedirs(result_path, exist_ok=True)


config = json.load(open(template_path, 'r'))


# 根据 data_path目录得到 dataset_list
dataset_list = []
info_list = []
outdir_base = []
filelist_path_list = []
for dr in os.listdir(data_path):
    dataset_list.append(dr)
    info_path = os.path.join(data_path, dr, 'info.JSON')
    info_list.append(os.path.join(data_path, dr, 'info.JSON'))
    outdir_base.append(os.path.join(result_path, dr))
    filelist_path_list.append(os.path.join(data_path, dr, 'filelist.txt'))
    print(f'dataset: {dr}')
    print(f'info_path: {info_path}')
    print(f'outdir_base: {outdir_base[-1]}')
    print(f'filelist_path: {filelist_path_list[-1]}')
    print(f'\n')



# 根据dataset和分辨率 创建结果目录
for res in reslutions_list:
    for dr_base,filelist,info in zip(outdir_base, filelist_path_list,info_list ):
        outdir = dr_base + '_' + res
        if res == '50k':
            intres = 50000
        elif res == '500k':
            intres = 500000
        elif res == '100k':
            intres = 100000
        elif res == '1Mb':
            intres = 1000000
        elif res == '10k':
            intres = 10000
        os.makedirs(outdir, exist_ok=True)
        # copy filelist
        shutil.copy(filelist, outdir)
        # zhunbei config
        print(info)
        info = json.load(open(info, 'r'))
        config['genome_reference_path'] = info['genome_reference_path']
        config['temp_dir'] = outdir
        config['resolution'] =  intres
        config['resolution_cell'] = intres
        config['minimum_distance'] = intres
        config['resolution_fh'] = [intres]
        config['data_dir'] = outdir
        config['temp_dir'] = outdir
        config['cytoband_path'] = info['cytoband_path']
        if info['wuzhong'] == 'hg38':
            chromList = ['chr' + str(i) for i in range(1, 23)] 
        else:
            chromList = ['chr' + str(i) for i in range(1, 20)] 
        config['chrom_list'] = chromList
        config['impute_list'] = chromList
        # save config to outdir
        config_path = os.path.join(outdir, 'config.JSON')
        print(config)
        with open(config_path, 'w') as f:
            json.dump(config, f, indent=4) 
        
        
# 创建slurm脚本
code = '/cluster2/home/Kangwen/0_benchmark/3_workflow/higashi.sh'
python_code = '/cluster2/home/Kangwen/0_benchmark/1_method/re_higashi/0_higashi_pipeline.py'




slurm_script = '''#!/bin/bash
#SBATCH -J higashi
#SBATCH --nodelist=node1
#SBATCH --output=log.%x.out
#SBATCH --error=log.%x.err #报错
#SBATCH --mail-type=all # 指定状态发生时，发送邮件通知，有效种类为（NONE, BEGIN, END, FAIL, REQUEUE, ALL）；
#SBATCH --mail-user=18516291316@163.com # 发送给指定邮箱；
#SBATCH --cpus-per-task=10  # 每个任务所需要的核心数，默认为1；
#SBATCH --open-mode=append
ulimit -s unlimited
ulimit -l unlimited

source activate /cluster/home/Kangwen/anaconda3/envs/ckw



'''


with open(code, 'w') as f:
    f.write(slurm_script)
    for res in reslutions_list:
        for dr_base in outdir_base :
    
            cmd = f'python {python_code} {dr_base}_{res}/config.JSON\n'
            f.write(cmd)
            f.write('sleep 1\n')
