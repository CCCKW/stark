#!/bin/bash
#SBATCH -J fasthigashi
#SBATCH --nodelist=node1
#SBATCH --output=log.%x.out
#SBATCH --error=log.%x.err #报错
#SBATCH --mail-type=all # 指定状态发生时，发送邮件通知，有效种类为（NONE, BEGIN, END, FAIL, REQUEUE, ALL）；
#SBATCH --mail-user=18516291316@163.com # 发送给指定邮箱；
#SBATCH --cpus-per-task=20  # 每个任务所需要的核心数，默认为1；
#SBATCH --open-mode=append
ulimit -s unlimited
ulimit -l unlimited

source activate /cluster/home/Kangwen/anaconda3/envs/ckw
config=''
script='/cluster/home/Kangwen/dream/1_method/re_fasthigashi/0_fasthigashi_pipeline.py'
python -u $script $config
###end