import os
import gzip
import sys
import argparse
import re
import numpy as np
import subprocess
import time
from .digest_genome import digest_genomes
from .inline_splitter import inline_splitter
import logging
from .analyze_scDHC_V2design import analyze_scDHC_V2design


def pp(opt, fastq,log_out):
    T = time.time()
    print('********run scHic*********')
    # 移动到工作目录
    os.chdir(opt['output'])

 

    # 创建文件夹tmp,并移动
    if not os.path.exists(opt['type'] + '_' + fastq + '_tmp'): # scHic_sample_tmp
        os.makedirs(opt['type'] + '_' + fastq + '_tmp')
    os.chdir(opt['type'] + '_' + fastq + '_tmp')

    # 创建logging文件
    logging.basicConfig(filename=fastq + '_logging.log', level=logging.DEBUG)
    
    
    # 创建logging文件
    logging.basicConfig(filename=fastq + '_logging.log', level=logging.DEBUG)
    if os.path.exists(opt['output'] + '/' + opt['type'] + '_' + fastq + '_tmp/Result/' + fastq +'.mcool'):
        log_out.info('The result of ' + fastq + ' has been generated')
        
        return



    fq_r1 = opt['fastq_dir'][opt['fastq_log'].index(fastq)]
    fq_r2 = fq_r1.replace('_1.fastq', '_2.fastq')



    if not os.path.exists('Result'):
        os.makedirs('Result')
    if not os.path.exists('Result/SCPair'):
        os.makedirs('Result/SCPair')

    # 默认inner和outer barcode和fastq文件放在一个目录下面
    opt['inner-barcode'] = os.path.dirname(opt['fastq_dir'][opt['fastq_log'].index(fastq)]) + '/' + fastq + '_inner_barcode.txt'
    opt['outer-barcode'] = os.path.dirname(opt['fastq_dir'][opt['fastq_log'].index(fastq)]) + '/' + fastq + '_outer_barcode.txt'

    # 开始运行pipeline
    # Get paird fastq
    cmd = 'fastp --adapter_sequence AGATCGGAAGAGCGATCGG '
    cmd += ' --adapter_sequence_r2 AGATCGGAAGAGCGTCGTG '
    cmd += ' --overrepresentation_analysis'
    cmd += ' --correction --thread ' + str(opt['thread'])
    cmd += ' --html trimmed.fastp.html '
    cmd += ' --json trimmed.fastp.json '
    cmd += ' -i ' + fq_r1 + ' -I ' + fq_r2
    cmd += ' --out1 trimmed-pair1.fastq.gz --out2 trimmed-pair2.fastq.gz'
    cmd += ' --failed_out failed.fq.gz'
    print('1++++++\n', cmd)
    t = time.time()
    logging.info(cmd)
    os.system(cmd)
    logging.info('fastp time: ' + str(time.time() - t))

    ####Get inner split pair File
    print('2+++++++++++++\n', '#Get inner split pair File')
    t = time.time()
    logging.info(cmd)
    inline_splitter('trimmed-pair1.fastq.gz', 'trimmed-pair2.fastq.gz', opt['outer-barcode'], 'trimmed-pair1.fastq.split', 'trimmed-pair2.fastq.split' ,fastq + '.splitting_stats.html')
    logging.info('inline_splitter time: ' + str(time.time() - t))

    #Run analyze_scDHC_V2design.py, which searches for the adaptor and clips it out
    print('3+++++++++++++\n#Run analyze_scDHC_V2design.py')
    t = time.time()
    logging.info(cmd)
    analyze_scDHC_V2design(opt['inner-barcode'],
                           'trimmed-pair1.fastq.split',
                           'trimmed-pair2.fastq.split',
                            'trimmed-pair1.bc_clipped',
                            'trimmed-pair2.bc_clipped',
                            'read_barcode.index')
    
    logging.info('analyze_scDHC_V2design time: ' + str(time.time() - t))
 
    

    cmd = 'awk -F\'\\t\' \'{print $1 \"\\t\"  $3 \"_\" $4}\''
    cmd +=  ' '  'read_barcode.index > ' + fastq + '.index'
    print('4+++++++\n' , cmd)
    t = time.time()
    logging.info(cmd)
    os.system(cmd)

    logging.info('awk time: ' + str(time.time() - t))

    # Align reads to combo-reference using bwa/bowtie2
    t = time.time()
    if opt['aligner'] == 'bwa':
        cmd = 'bwa mem ' + opt['index']
        cmd += ' -t ' + str(opt['thread'])
        cmd += ' trimmed-pair1.fastq.gz trimmed-pair2.fastq.gz > '  + fastq + '.sam'
        # print(os.getcwd())
        print('5+++++++\n', cmd)
        os.system(cmd)
    else:
        cmd = 'bowtie2 -p ' + str(opt['thread'])
        cmd += ' -x ' + opt['index']
        cmd += ' -1 trimmed-pair1.fastq.gz -2 trimmed-pair2.fastq.gz -S ' + fastq + '.sam'
        print('5+++++++\n', cmd)
        os.system(cmd)
    logging.info(cmd)
    logging.info(opt['aligner'] + ' time: ' + str(time.time() - t))

    # sam to bam
    cmd = 'samtools view -@ 8 -S ' + fastq + '.sam -1b -o ' + fastq + '.bam'
    print('6+++++++\n',cmd)
    t = time.time()
    logging.info(cmd)
    os.system(cmd)

    logging.info('samtools view time: ' + str(time.time() - t))

        # Generate pairs
    cmd = 'pairtools parse ' + fastq + '.bam ' + ' -c ' + opt['genomesize']
    cmd += ' --drop-sam --drop-seq --output-stats ' + fastq + '.stats '
    cmd += ' --assembly ' + opt['species']  + ' --no-flip '
    cmd += ' --add-columns ' + opt['add_columns']
    cmd += ' --walks-policy all '
    cmd += '-o ' + fastq + '.pairs.gz'
    cmd += ' --nproc-in ' + str(opt['thread']) +  ' --nproc-out ' + str(opt['thread'])
    print('7++++++++++++\n',cmd) 
    t = time.time()
    logging.info(cmd)
    os.system(cmd)
 
    logging.info('pairtools parse time: ' + str(time.time() - t))

    # #QC select
    cmd = 'pairtools select ' + '\"' + opt['select'] + '\"'
    cmd += ' ' + fastq + '.pairs.gz -o ' + fastq + '.filtered.pairs.gz'
    print('8++++++++++\n', cmd)
    t = time.time()
    logging.info(cmd)
    os.system(cmd)
    logging.info('pairtools select time: ' + str(time.time() - t))

    # Generate restriction Site
    print('9++++++++++\n', 'Generate restriction Site')
    t = time.time()
    logging.info(cmd)
    digest_genomes(opt['fa_index'], ['dpnii'], 'Speces_restrect.bed')
    logging.info('digest_genomes time: ' + str(time.time() - t))

    ###Calculate restricted pairs
    cmd = 'pairtools restrict -f Speces_restrect.bed '
    cmd += fastq + '.filtered.pairs.gz -o ' + fastq + '.filtered.restricted.pairs.gz'
    print('10++++++++++\n', cmd)
    t = time.time()
    logging.info(cmd)
    os.system(cmd)
    logging.info('pairtools restrict time: ' + str(time.time() - t))


    #Sort pairs
    cmd = 'pairtools sort --nproc ' + str(opt['thread']) + ' ' + fastq + '.filtered.restricted.pairs.gz -o ' + fastq + '.filtered.restricted.sorted.pairs.gz'
    print('11++++++++++\n', cmd)
    t = time.time()
    logging.info(cmd)
    os.system(cmd)
    logging.info('pairtools sort time: ' + str(time.time() - t))


    # Split Cell
    print('Split Cell')
    t = time.time()
    ReadDict, SClist = GenerateSCDict(fastq+".index")
    GenerateSCPairs(ReadDict, SClist, fastq + '.filtered.restricted.sorted.pairs.gz', "./Result/SCPair")
    logging.info('Split Cell time: ' + str(time.time() - t))


     #dedup
    cmd = 'pairtools dedup '
    cmd += '--max-mismatch ' + str(opt['max_mismatch'])
    cmd += ' --mark-dups --output ' + '>( pairtools split --output-pairs ' + fastq + \
        '.nodups.pairs.gz --output-sam ' + fastq + \
        '.nodups.bam ) --output-unmapped >( pairtools split --output-pairs  ' + fastq + \
        '.unmapped.pairs.gz --output-sam  ' + fastq + \
        '.unmapped.bam ) --output-dups >( pairtools split --output-pairs  ' + fastq + \
        '.dups.pairs.gz --output-sam ' + fastq + \
        '.dups.bam ) --output-stats  ' + fastq + \
        '.dedup.stats ' + fastq + \
        '.filtered.restricted.sorted.pairs.gz'
    print('12+++++++++\n', cmd)
    t = time.time()
    logging.info(cmd)
    subprocess.run(cmd, shell=True, executable="/bin/bash")
    logging.info('pairtools dedup time: ' + str(time.time() - t))

        # 
    cmd = 'gunzip ' + fastq + '.nodups.pairs.gz'
    print('13+++++++++++\n',cmd)
    t = time.time()
    logging.info(cmd)
    os.system(cmd)
    logging.info('gunzip time: ' + str(time.time() - t))

    cmd = 'cooler cload pairs -c1 2 -p1 3 -c2 4 -p2 5 '
    cmd += opt['genomesize'] + ':' + str(opt['resolution']) + ' ' + fastq + '.nodups.pairs ' + fastq + str(opt['resolution']) + '.cool'
    print('14+++++++++++\n', cmd)
    t = time.time()
    logging.info(cmd)
    os.system(cmd)
    logging.info('cooler cload time: ' + str(time.time() - t))

    cmd = 'gzip ' + fastq + '.nodups.pairs'
    print('15+++++++++++\n', cmd)
    t = time.time()
    os.system(cmd)
    logging.info('gzip time: ' + str(time.time() - t))
   
    # KR correction
    cmd = 'hicCorrectMatrix correct --matrix  ' + fastq + str(opt['resolution']) + '.cool '
    cmd += ' --correctionMethod KR --outFileName '
    cmd += ' ' + fastq + str(opt['resolution']) + '.KR.cool'
    cmd += ' --chromosomes chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr 20 chr21 chr22 chrX chrY '
    print('16+++++++++++\n', cmd)
    t = time.time()
    logging.info(cmd)
    os.system(cmd)
    logging.info('hicCorrectMatrix time: ' + str(time.time() - t))


    cmd = 'cooler zoomify '  + fastq + str(opt['resolution']) + '.KR.cool -r 10000,50000,100000,500000 -o ' + fastq + '.mcool'
    print('17+++++++++++\n', cmd)
    t = time.time()
    logging.info(cmd)
    os.system(cmd)
    logging.info('cooler zoomify time: ' + str(time.time() - t))
    logging.info('total time: ' + str(time.time() - T))
   





def GetReadsIndex(file1_path, file2_path,indexName):

    #####记录barcode对应的reads
    fp = open(indexName, "w")
    DirtyFile = open("DirtyReads.txt", "w")

    # 打开两个文件并同时读取内容
    with gzip.open(file1_path, 'rt') as file1, gzip.open(file2_path, 'rt') as file2:
        # 同时逐行读取两个文件的内容
        for line1, line2 in zip(file1, file2):
            # 处理文件内容
            # 例如，打印每一行
            line1 = str(line1)
            line2 = str(line2)
            if '@' == line1[0] and '@' == line2[0]:
                if line1.split(" ")[0] == line2.split(" ")[0]:
                    if line1.split(" ")[1][-18:-1] == line2.split(" ")[1][-18:-1]:
                        fp.write(line1.split(" ")[0][1:]+"\t"+line2.split(" ")[1][-18:-1]+"\n")
                    else:
                        DirtyFile.write("Reads1\t"+line1[0:-1]+"\tReads2\t"+line2)
                else:
                    DirtyFile.write("Reads1\t" + line1[0:-1] + "\tReads2\t" + line2)
    fp.close()




def GenerateSCDict(ReadIndex):
    ReadDict = {}
    SCset = set()
    count = 0
    startTime = time.time()
    with open(ReadIndex, "r") as file:
        indexInfo = file.readlines()
        print("index has loaded")
        for line in indexInfo:
            IndexComponent = line.rstrip().split("\t")
            barcode = IndexComponent[1]
            ReadDict[IndexComponent[0]] = barcode
            SCset.add(barcode)
            count+=1
            if count%1000000 == 0:
                print("%d index has processed in %f s"%(count, time.time() - startTime))
                startTime = time.time()

    SClist = list(SCset)
    return ReadDict, SClist


def GenerateSCPairs(ReadDict, SClist, pairFile, ResultDir):
    SCPairs = {barcode: [] for barcode in SClist}
    HeadInfo = []
    count = 0
    print("*****Start make single cell pairs*****")
    with gzip.open(pairFile, "rt") as pf:
        startTime = time.time()
        for line in pf:
            if "#" in line:
                HeadInfo.append(line)
                continue
            PairComponent = line.split("\t")
            barcode = ReadDict[PairComponent[0]]
            SCPairs[barcode].append(line)
            count += 1
            if count % 100000 == 0:
                print("%d pairs has processed in %f s" % (count, time.time() - startTime))
                startTime = time.time()

    print("*****Start generate single cell pairs files*****")
    cell_count = 0
    scPairsSummary = open("singleCell_HIC_pairs_summary.txt", "w")
    scPairsSummary.write("Cell Num" + "\t" + "Cell Barcode" + "\t" + "Pairs Count" + "\n")
    for i in range(len(SClist)):
        if len(SCPairs[SClist[i]])>1:
            scPairsSummary.write(str(i+1)+"\t"+str(SClist[i])+"\t"+str(len(SCPairs[SClist[i]]))+"\n")
        if len(SCPairs[SClist[i]])<1000:
            continue
        cell_count+=1
        print("Process cell %d,\t"%(cell_count),end='')
        fp = open(ResultDir+"/cell_"+str(cell_count)+"_"+str(i)+".pairs","w")
        for singel_head in HeadInfo:
            fp.write(singel_head)
        count=0

        for pair in SCPairs[SClist[i]]:
            fp.write(pair)
            count+=1
        print("%d pairs"%(count))
        fp.close()
    scPairsSummary.close()