import argparse
import re
import numpy as np
import subprocess
import time
import os
import logging
from .scaffold import *
import pandas as pd
from joblib import Parallel, delayed


def assembly(type_,opt, fastq,log_out):

    print('********run %s || %s*********' % (type_, fastq))
    os.chdir(opt['output'])
    start = time.time()

    # 创建文件夹tmp,并移动
    if not os.path.exists(opt['type'] + '_' + fastq + '_tmp'):
        os.makedirs(opt['type'] + '_' + fastq + '_tmp')
    os.chdir(opt['type'] + '_' + fastq + '_tmp')
    
    # 创建logging文件
    logging.basicConfig(filename=fastq + '_logging.log', level=logging.DEBUG)
    for k in opt.keys():
        if (k == 'fastq_dir') or (k == 'fastq_log') or (k == 'run_sample'):
            continue
        else:
            logging.debug('# ' + k + ': ' + str(opt[k]))
        
    if os.path.exists(opt['output'] + '/' + opt['type'] + '_' + fastq + '_tmp/' + fastq +'.bam'):
        log_out.info('The result of ' + fastq + ' has been generated')
        pass
    if os.path.exists(opt['output'] + '/' + opt['type'] + '_' + fastq + '_tmp/Result/SCpair'):
        for val in os.listdir(opt['output'] + '/' + opt['type'] + '_' + fastq + '_tmp/Result/SCpair'):
            if val.endswith('.mcool'):
                return

    fq_r1 = opt['fastq_dir'][opt['fastq_log'].index(fastq)]
    fq_r2 = fq_r1.replace('_1.fastq', '_2.fastq')
    
    if type_ == 'scHic' and opt['exist_barcode']:
        fq_r3 = fq_r1.replace('_1.fastq', '_3.fastq')
        fq_r4 = fq_r1.replace('_1.fastq', '_4.fastq')
        if not os.path.exists(fq_r3) or not os.path.exists(fq_r4):
            print(fastq , fq_r3, fq_r4)
            print('The barcode fastq file is not exist')
            log_out.error('The barcode fastq file is not exist')
            return
    
    threads = opt['thread']
    
    if not os.path.exists('Result'):
        os.makedirs('Result')

    # 默认inner和outer barcode和fastq文件放在一个目录下面
    if type_ == 'sciHic':
        opt['inner-barcode'] = os.path.dirname(opt['fastq_dir'][opt['fastq_log'].index(fastq)]) + '/' + fastq + '_inner_barcodes.txt'
        opt['outer-barcode'] = os.path.dirname(opt['fastq_dir'][opt['fastq_log'].index(fastq)]) + '/' + fastq + '_outer_barcodes.txt'

    # ==================== FASTP 处理步骤 ====================
    # 根据不同类型进行不同的前处理
    
    if type_ == 'droplet':
        # DROPLET 特殊处理流程
        fastq_dir = os.path.dirname(fq_r1)
        
        # 1. 使用hictools合并数据
        t = hictools_combine(fastq_dir, fastq, 'atac')
        logging.info(f'hictools_combine: {t}')
        
        # 2-5. droplet特有的barcode处理步骤
        t = process_droplet_barcode_steps(opt, fastq, fastq_dir, threads)
        logging.info(f'droplet_barcode_processing: {t}')
        
        # 6. fastp质控 - 使用droplet处理后的文件
        fq_r1_processed = f'{fastq_dir}/{fastq}_R1_BC_cov.fq'
        fq_r2_processed = f'{fastq_dir}/{fastq}_R3_BC_cov.fq'
        
        t = fastp(fq_r1_processed, fq_r2_processed, 
                  'trimmed-pair1.fastq.gz', 'trimmed-pair2.fastq.gz',
                  threads, max=False)
        logging.info('fastp::: %s '%t)
        
        # 7. 清理droplet特有的中间文件
        cleanup_droplet_intermediate_files(fastq_dir, fastq)
        
    elif opt['exist_barcode']:
        # 其他有barcode的类型处理
        if type_ == 'GAGE-seq':
            from sc3dg.utils.muti_demultiplex import main as muti_demultiplex
            adaptor_file_bc2 = opt['adaptor-file-bc2']
            adaptor_file_l2 = opt['adaptor-file-l2']
            adaptor_file_bc1 = opt['adaptor-file-bc1']
            outfile1 = 'demulti_R1.fastq.gz'
            outfile2 = 'demulti_R2.fastq.gz'
            barcode = [adaptor_file_bc2, adaptor_file_l2, adaptor_file_bc1]
            pos1 = [[], [], []]
            pos2 = [slice(0, 8), slice(8, 23), slice(23, 31)]
            mode = [1, 5, 1]
            unknown1 = None
            unknown2 = None
            barcode2keep = None
            chunk_size = 100000
            nproc = threads
            
            muti_demultiplex(fq_r1, fq_r2, outfile1, outfile2, barcode, pos1, pos2, mode, 
                           unknown1, unknown2, barcode2keep, chunk_size, nproc)
            
            t = fastp('demulti_R1.fastq.gz', 'demulti_R2.fastq.gz',
                      'trimmed-pair1.fastq.gz', 'trimmed-pair2.fastq.gz',
                      threads, max=False)
             
        elif type_ == 'snHic':
            t = fastp(fq_r1, fq_r2, 
                      'trimmed-pair1.fastq.gz', 'trimmed-pair2.fastq.gz',
                      threads, max=False)
            snHic_barcode('trimmed-pair1.fastq.gz', 'trimmed-pair2.fastq.gz')

        elif type_ == 'scHic':
            print('making index..............')
            t = GetFastqIndex(fastq, fq_r1, fq_r4, fq_r2, fq_r3)
            logging.info('make_index::: %s '%t)

            t = fastp('%s._1_barcode.fastq.gz'%fastq,
                    '%s._2_barcode.fastq.gz'%fastq, 
                    'trimmed-pair1.fastq.gz', 
                    'trimmed-pair2.fastq.gz', threads, max=True)
            logging.info('fastp::: %s '%t)

        elif type_ == 'sciHic':
            t = fastp(fq_r1, fq_r2, 'trimmed-pair1_before.fastq.gz',
                        'trimmed-pair2_before.fastq.gz',
                        threads, max=False,
                        adapter=['AGATCGGAAGAGCGATCGG', 'AGATCGGAAGAGCGTCGTG'])
            logging.info('fastp::: %s '%t)

            t = time.time()
            inline_splitter('trimmed-pair1_before.fastq.gz', 
                            'trimmed-pair2_before.fastq.gz',
                            opt['outer-barcode'], 
                            'trimmed-pair1_before2.fastq.gz', 
                            'trimmed-pair2_before2.fastq.gz',
                            fastq + '.splitting_stats.html')
            logging.info('inline_splitter::: %s '% (time.time()- t))

            t = time.time()
            analyze_scDHC_V2design(opt['inner-barcode'],
                            'trimmed-pair1_before2.fastq.gz', 
                            'trimmed-pair2_before2.fastq.gz',
                            'trimmed-pair1.fastq.gz', 
                            'trimmed-pair2.fastq.gz')
            logging.info('analyze_scDHC_V2design::: %s ' % (time.time()- t))

            t = time.time()
            awk()
            logging.info('awk::: %s ' % (time.time()- t))
    else:
        # 无barcode的常规处理
        t = fastp(fq_r1, fq_r2, 'trimmed-pair1.fastq.gz', 'trimmed-pair2.fastq.gz', threads, max=False)
        logging.info('fastp::: %s '%t)

    # ==================== 通用处理步骤 ====================
    # 从这里开始所有类型都使用相同的流程
    
    # bwa比对
    if opt['aligner'] == 'bwa':
        if type_ == 'droplet':
            t = bwa_mem_droplet(opt['index'], threads, fastq)
        else:
            t = bwa(opt['index'], threads, fastq)
        logging.info('bwa::: %s '%t)
    else:
        t = bowtie2(opt['index'], threads, fastq)
        logging.info('bowtie2::: %s '%t)

    # 比对的sam转bam
    t = sam_to_bam(fastq, threads)
    logging.info('sam_to_bam::: %s '%t)

    # 质控bam文件
    t = samtools_qc(threads, opt['qc'], fastq)
    logging.info('samtools_qc::: %s '%t)
    
    # pairtools parse
    t = pairtools_parse(fastq, 
                        opt['genomesize'],
                        opt['species'],
                         opt['add_columns'],
                         threads)
    logging.info('pairtools_parse::: %s '%t)

    # pairtools restrict
    t = pairtools_restrict(fastq, opt['enzyme_bed'])
    logging.info('pairtools_restrict::: %s '%t)

    # pairtools select
    t = pairtools_select(fastq, opt['select'])
    logging.info('pairtools_select::: %s '%t)

    # pairtools sort 
    t = pairtools_sort(fastq)
    logging.info('pairtools_sort::: %s '%t)

    # pairtools dedup
    t = pairtools_dedup(fastq, str(opt['max_mismatch']))
    logging.info('pairtools_dedup::: %s '%t)

    # 整体的mcool
    t = cooler_cload_pairs(opt['genomesize'],
                            str(opt['resolution']),
                            fastq, prefix='dedup')
    logging.info('cooler_cload_pairs::: %s '%t)

    if isinstance(opt['zoomify_res'], list):
        opt['zoomify_res'] = ','.join(opt['zoomify_res'])

    t = cooelr_zoomify(fastq, opt['resolution'], opt['zoomify_res'])
    logging.info('cooelr_zoomify::: %s '%t)

    # 单细胞分析处理
    if opt['exist_barcode']:
        if type_ == 'GAGE-seq':
            from sc3dg.utils.parseSCpairs import main as split_cells
            split_cells(fastq, './Result/SCpair')
        else:
            split_cells(fastq)
        t = time.time()
        pair_list = []
        for p in os.listdir('./Result/SCpair'):
            if p.endswith('pairs.gz'):
                pair_list.append(p)
    
        Parallel(n_jobs=6)(delayed(cload_correct_zoomify)('./Result/SCpair/',
                                      p,
                                      opt['genomesize'],
                                      str(opt['resolution']),
                                       opt['zoomify_res']
                            ) for p in pair_list)
        logging.info('sub_cload_correct_zoomify::: %s '% (time.time() - t))
 
    mv_Result()
    logging.info('total time: %s' % (time.time() - start))