import os
import gzip
import sys
import argparse
import re
import numpy as np
import subprocess
import time
from .digest_genome import digest_genomes
import logging


def pp(opt, fastq):
    print('********run scHic*********')
    T = time.time()
    # 移动到工作目录
    os.chdir(opt['output'])

    # 创建文件夹tmp,并移动
    if not os.path.exists(opt['type'] + '_' + fastq + '_tmp'): # scHic_sample_tmp
        os.makedirs(opt['type'] + '_' + fastq + '_tmp')
    os.chdir(opt['type'] + '_' + fastq + '_tmp')

    fq_r1 = opt['fastq'] + '/' + fastq + '_1.fastq' + opt['fastq_prefix']
    fq_r2 = opt['fastq'] + '/' + fastq + '_2.fastq' + opt['fastq_prefix']

    # 创建log
    logging.basicConfig(filename=fastq + '_logging.log', level=logging.DEBUG)
    for k in opt.keys():
        logging.debug('# ' + k + ': ' + str(opt[k]))

    
    # fastp
        # Get paird fastq
    cmd = 'fastp --overrepresentation_analysis'
    cmd += ' --correction --thread ' + str(opt['thread'])
    cmd += ' --html trimmed.fastp.html '
    cmd += ' --json trimmed.fastp.json '
    cmd += ' -i ' + fq_r1 + ' -I ' + fq_r2
    cmd += ' --out1 trimmed-pair1.fastq.gz --out2 trimmed-pair2.fastq.gz'
    cmd += ' --failed_out failed.fq.gz'
    t = time.time()
    logging.info(cmd)
    os.system(cmd)
    logging.info('fastp time: ' + str(time.time() - t))


    # Align reads to combo-reference using bwa/bowtie2
    t = time.time()
    if opt['aligner'] == 'bwa':
        cmd = 'bwa mem ' + opt['index']
        cmd += ' -t ' + str(opt['thread'])
        cmd += ' trimmed-pair1.fastq.gz trimmed-pair2.fastq.gz > '  + fastq + '.sam'
        os.system(cmd)
    else:
        cmd = 'bowtie2 -p ' + str(opt['thread'])
        cmd += ' -x ' + opt['index']
        cmd += ' -1 trimmed-pair1.fastq.gz -2 trimmed-pair2.fastq.gz -S ' + fastq + '.sam'
        os.system(cmd)
    logging.info(cmd)
    logging.info('align time: ' + str(time.time() - t))

    # sam to bam
    cmd = 'samtools view -@ 8 -S ' + fastq + '.sam -1b -o ' + fastq + '.bam'
    t = time.time()
    logging.info(cmd)
    os.system(cmd)
    logging.info('sam to bam time: ' + str(time.time() - t))

    #Get Reads Index
    print('Get Reads Index')
    t = time.time()
    logging.info('Get Reads Index')
    GetReadsIndex("trimmed-pair1.fastq.gz","trimmed-pair2.fastq.gz",fastq+".index")
    logging.info('Get Reads Index time: ' + str(time.time() - t))

    # Generate pairs
    cmd = 'pairtools parse ' + fastq + '.bam ' + ' -c ' + opt['genomesize']
    cmd += ' --drop-sam --drop-seq --output-stats ' + fastq + 'stats '
    cmd += ' --assembly ' + opt['species']  + ' --no-flip '
    cmd += ' --add-columns ' + opt['add_columns']
    cmd += ' --walks-policy all '
    cmd += '-o ' + fastq + '.pairs.gz'
    cmd += ' --nproc-in ' + str(opt['thread']) +  ' --nproc-out ' + str(opt['thread'])
    logging.info(cmd)
    os.system(cmd)

    # #QC select
    cmd = 'pairtools select ' + '\"' + opt['select'] + '\"'
    cmd += ' ' + fastq + '.pairs.gz -o ' + fastq + '.filtered.pairs.gz'
    print('5++++++++++\n', cmd)
    t = time.time()
    logging.info(cmd)
    os.system(cmd)
    logging.info('pairtools select time: ' + str(time.time() - t))



    # Split Cell
    # 创建结果目录
    print('Split Cell')
    if not os.path.exists('Result'):
        os.makedirs('Result')
    if not os.path.exists('Result/SCPair'):
        os.makedirs('Result/SCPair')
    t = time.time()
    logging.info('Split Cell')
    ReadDict, SClist = GenerateSCDict(fastq+".index")
    GenerateSCPairs(ReadDict, SClist, fastq+".filtered.pairs.gz", "./Result/SCPair")
    logging.info('Split Cell time: ' + str(time.time() - t))

    # Sort pairs
    cmd ='pairtools sort --nproc ' + str(opt['thread'])
    cmd += ' ' + fastq + '.filtered.pairs.gz -o ' + fastq + '.filtered.sorted.pairs.gz'
    print('6++++++\n', cmd)
    t = time.time()
    logging.info(cmd)
    os.system(cmd)
    logging.info('pairtools sort time: ' + str(time.time() - t))



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
        '.filtered.sorted.pairs.gz'
    print('7+++++++++\n', cmd)
    t = time.time()
    logging.info(cmd)
    subprocess.run(cmd, shell=True, executable="/bin/bash")
    logging.info('pairtools dedup time: ' + str(time.time() - t))

    # QC select
    cmd = 'pairtools select ' + '\"' + opt['select']  + '\"'
    cmd += ' ' + fastq + '.nodups.pairs.gz -o ' + fastq + '.nodups.UU.pairs.gz'
    print('8++++++++++++\n', cmd)
    t = time.time()
    logging.info(cmd)
    os.system(cmd)
    logging.info('pairtools select time: ' + str(time.time() - t))

    # Generate restriction Site
    print('9++++++++++\n', 'Generate restriction Site')
    t = time.time()
    logging.info('digest_genomes')
    digest_genomes(opt['fa_index'], opt['enzyme'], 'Speces_restrect.bed')
    logging.info('digest_genomes time: ' + str(time.time() - t))

        ###Calculate restricted pairs
    cmd = 'pairtools restrict -f Speces_restrect.bed '
    cmd += fastq + '.nodups.UU.pairs.gz -o ' + fastq + '.nodups.restricted.pairs.gz'
    print('10++++++++++\n', cmd)
    t = time.time()
    logging.info(cmd)
    os.system(cmd)
    logging.info('pairtools restrict time: ' + str(time.time() - t))


    # 
    cmd = 'gunzip ' + fastq + '.nodups.restricted.pairs.gz'
    print('11+++++++++++\n',cmd)
    t = time.time()
    logging.info(cmd)
    os.system(cmd)
    logging.info('gunzip time: ' + str(time.time() - t))
    
    cmd = 'cooler cload pairs -c1 2 -p1 3 -c2 4 -p2 5 '
    cmd += opt['genomesize'] + ':' + str(opt['resolution']) + ' ' + fastq + '.nodups.restricted.pairs ' + fastq + str(opt['resolution']) + '.cool'
    print('12+++++++++++\n', cmd)
    t = time.time()
    logging.info(cmd)
    os.system(cmd)
    logging.info('cooler cload time: ' + str(time.time() - t))
   
    cmd = 'gzip ' + fastq + '.nodups.restricted.pairs'
    print('13+++++++++++\n', cmd)
    t = time.time()
    logging.info(cmd)
    os.system(cmd)
    logging.info('gzip time: ' + str(time.time() - t))

    # KR correction
    cmd = 'hicCorrectMatrix correct --matrix  ' + fastq + str(opt['resolution']) + '.cool '
    cmd += ' --correctionMethod KR --outFileName '
    cmd += ' ' + fastq + str(opt['resolution']) + '.KR.cool'
    cmd += ' --filterThreshold -1.5 5.0 --chromosomes chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX '
    print('14+++++++++++\n', cmd)
    t = time.time()
    logging.info(cmd)
    os.system(cmd)
    logging.info('hicCorrectMatrix time: ' + str(time.time() - t))
   
    

    cmd = 'cooler zoomify '  + fastq + str(opt['resolution']) + '.KR.cool -r 10000,50000,100000,500000 -o ' + fastq + '.mcool'
    print('15+++++++++++\n', cmd)
    t = time.time()
    logging.info(cmd)
    os.system(cmd)
    logging.info('cooler zoomify time: ' + str(time.time() - t))






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
            if '@' == line1[0] and '@' == line2[0] and len(line1.split("\t")) == 2 and len(line2.split("\t")) == 2:
                if line1.split(" ")[0] == line2.split(" ")[0]:
                    if line1.split("\t")[1][0:-1] == line2.split("\t")[1][0:-1]:
                        fp.write(line1.split(" ")[0][1:]+"\t"+line2.split("\t")[1][0:-1]+"\n")
                    else:
                        DirtyFile.write("Reads1\t"+line1[0:-1]+"\tReads2\t"+line2)
                else:
                    DirtyFile.write("Reads1\t" + line1[0:-1] + "\tReads2\t" + line2)
    fp.close()
    DirtyFile.close()


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