import os
import gzip
import sys
import argparse
import re
import numpy as np
import pandas as pd
import subprocess
import time
import logging
from pathlib import Path
from .scNano_barcode import generate_barcode
from .ReformSAM import reSAM
from .Generate_scNanoHIC_pairs import WriteCellPairs
from multiprocessing import Pool
barcodesh = '/cluster/home/Kangwen/Hic/ckw/utils/barcode.sh'

def run_command(command):
    print(command)
    subprocess.run(command, shell=True, timeout=None)
    

def pp(opt, fastq,log_out):
    print('==========run scNano================')

 
    T = time.time()
    # 移动到工作目录
    os.chdir(opt['output'])

    # 创建文件夹tmp,并移动
    if not os.path.exists(opt['type'] + '_' + fastq + '_tmp'): # scHic_sample_tmp
        os.makedirs(opt['type'] + '_' + fastq + '_tmp')
    output = opt['output'] + '/' +opt['type'] + '_' + fastq + '_tmp'
    os.chdir(opt['type'] + '_' + fastq + '_tmp')

    fq = opt['fastq_dir'][opt['fastq_log'].index(fastq)]
    PCR = os.path.dirname(fq) + '/' + 'PCR.txt'
    TN5 = os.path.dirname(fq) + '/' + 'TN5.txt'
    bc = os.path.dirname(fq) + '/' + 'index.txt'

    print(PCR,TN5,bc)
    if not os.path.exists(PCR) or not os.path.exists(TN5) or not os.path.exists(bc):
        print('not exit pcr, tn5, bc file')
        return
    else:
        print('exit pcr, tn5, bc file')

    # 创建log
    logging.basicConfig(filename=fastq + '_logging.log', level=logging.DEBUG)
    for k in opt.keys():

     
        if (k == 'fastq_dir') or (k == 'fastq_log') or (k == 'run_sample'):
            continue
        else:
            logging.debug('# ' + k + ': ' + str(opt[k]))

    
    
    # 设置文件夹路径
    raw_folder = './raw_data'
    trim_folder = './trim_data'
    pairs_folder = './pairs_data'
    mapping_folder = './mapping_data'
    digestpairs_folder = './digestpairs_folder'
    filtedpairs_folder = './filted_pairs_data'
    sortedpairs_folder = './sorted_pairs_data'
    cool_folder = './cool_data'
    KRcool_folder = './KRcool_data'
    mcool_folder = './Mcool_data'

    # 创建必要的文件夹
    os.makedirs(trim_folder, exist_ok=True)
    os.makedirs(raw_folder, exist_ok=True)
    os.makedirs(pairs_folder, exist_ok=True)
    os.makedirs(mapping_folder, exist_ok=True)
    os.makedirs(digestpairs_folder, exist_ok=True)
    os.makedirs(filtedpairs_folder, exist_ok=True)
    os.makedirs(sortedpairs_folder, exist_ok=True)
    os.makedirs(cool_folder, exist_ok=True)
    os.makedirs(KRcool_folder, exist_ok=True)
    os.makedirs(mcool_folder, exist_ok=True)

    for val in os.listdir(opt['output'] + '/' + opt['type'] + '_' + fastq + '_tmp/Mcool_data'):
        if val.endswith('.mcool'):
            log_out.info('The result of ' + fastq + ' has been generated')
            print('The result of ' + fastq + ' has been generated')
            return

    with open(bc , 'r') as f:
        lines = f.readlines()
    lines = [x.strip() for x  in lines]

    divide_list = []
    bcname = []
    for i in range(0, len(lines), 3):
        divide_list.append(lines[i:i + 3])
        bcname.append('barcode_%s.fa' % i)


    for i,out in enumerate(bcname):
        print(out, divide_list[i])
        generate_barcode(divide_list[i],PCR, TN5,out)
    print('cjw')
    for fa in bcname:
        node = time.time()
        print(fa)
        cmd = 'nanoplexer -b %s -B 200M -t %s -p %s %s' % (fa, opt['thread'],raw_folder, fq)
        run_command(cmd)
        print('nanoplexer:%s'% (time.time() - node))

    
    # 第一步: 运行 cutadapt 处理所有文件
    for file_name in os.listdir(raw_folder):
        if file_name.endswith(".fastq"):
            input_file = os.path.join(raw_folder, file_name)
            output_file = os.path.join(trim_folder, f"{file_name}_Trim.fastq")
            log_file = os.path.join(trim_folder, f"{file_name}_Trim.log")
            cutadapt_command = f"cutadapt -g AATGATACGGCGACCACCGAGATCT -g TCGTCGGCAGCGTC -g AGATGTGTATAAGAGACAG -a AGATCTCGGTGGTCGCCGTATCATT -a GACGCTGCCGACGA -a CTGTCTCTTATACACATCT --times 2 -j {opt['thread']} --minimum-length 500 -e 0.2 -O 12 -o {output_file} {input_file} > {log_file}"
            run_command(cutadapt_command)
    #         break

    
    for input_file in os.listdir(trim_folder):
        if input_file.endswith('.fastq'):
            # 提取文件名（不包括路径和扩展名）
            file_name = os.path.splitext(input_file)[0]
            
            # 设置输出文件的路径
            output_sam = os.path.join(mapping_folder, f'{file_name}.sam')
            reformed_sam = os.path.join(mapping_folder, f'{file_name}_reformed.sam')

            # 运行miniMap2命令
            reference_genome = opt['index']
            print('minimap2')
            t_=opt['thread']
            minimap2_cmd = f'minimap2 -ax map-ont -t {t_} {reference_genome} {os.path.join(trim_folder, input_file)} > {output_sam}'
            subprocess.run(minimap2_cmd, shell=True)

            # 运行Python脚本
            print('reSAM')
            reSAM(input_sam=output_sam,
                   output_sam=reformed_sam,
                   input_is_bam=False,
                     output_is_bam=False, 
                     set_bx_flag=False)

            WriteCellPairs(opt['genomesize'], reformed_sam,opt['species'],pairs_folder)
       

    # 查找指定目录中的所有 .pairs 文件
  
    for input_file in os.listdir(pairs_folder):
        if input_file.endswith(".pairs"):
            print(input_file)
            input_file_path = os.path.join(pairs_folder, input_file)
            file_name = os.path.splitext(input_file)[0].replace("_Trim_reformed", "")
            print(file_name)
            # 定义文件路径
            digest_pairs = os.path.join(digestpairs_folder, f"{file_name}_digest.pairs")
            filted_pairs = os.path.join(filtedpairs_folder, f"{file_name}_filted.pairs")
            sorted_pairs = os.path.join(sortedpairs_folder, f"{file_name}_sorted.pairs")

            # 执行 pairtools 命令
            print('restricting....')
            subprocess.run(["pairtools", "restrict", "-f", opt['enzyme_bed'], input_file_path, "-o", digest_pairs])
            print('selecting......')
            subprocess.run(["pairtools", "select", "mapq1>=30 and mapq2>=30", digest_pairs, "-o", filted_pairs])
            print('sorting.....')
            subprocess.run(["pairtools", "sort", "--nproc", str(opt['thread']), filted_pairs, "-o", sorted_pairs])

            resolution = int(opt['resolution'])
            cool_file = os.path.join(cool_folder, f"{file_name}.{resolution}.cool")
            KR_cool_file = os.path.join(KRcool_folder, f"{file_name}.{resolution}.KR.cool")
            mcool_file = os.path.join(mcool_folder, f"{file_name}.mcool")
            # 执行 cooler 和 hicCorrectMatrix 命令
            chrom_size = opt['genomesize']
            subprocess.run(["cooler", "cload", "pairs", "-c1", "2", "-p1", "3", "-c2", "4", "-p2", "5", f"{chrom_size}:{resolution}", sorted_pairs, cool_file])
            # subprocess.run(["hicCorrectMatrix", "correct", "--matrix", cool_file, "--correctionMethod", "KR", "--outFileName", KR_cool_file, "--filterThreshold", "-1.5", "5.0", "--chromosomes", "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19"])

            subprocess.run(['cooler', 'zoomify',cool_file,'-r', opt['zoomify_res'], '-o', mcool_file])