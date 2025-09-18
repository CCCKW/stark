import gzip
import sys
import argparse
import re
import numpy as np
import subprocess
import time
import logging
import os
from Bio import SeqIO
from .scaffold import *
from joblib import Parallel, delayed

def pp(opt, fastq,log_out):

    
    print('********run snMethyl*********')
    T = time.time()
    # 移动到工作目录
    os.chdir(opt['output'])
    opt['zoomify_res'] = ','.join(opt['zoomify_res'])
        # 创建文件夹tmp,并移动
    if not os.path.exists(opt['type'] + '_' + fastq + '_tmp'): # scHic_sample_tmp
        os.makedirs(opt['type'] + '_' + fastq + '_tmp')
    

    os.chdir(opt['type'] + '_' + fastq + '_tmp')
    if os.path.exists(opt['output'] + '/' + opt['type'] + '_' + fastq + '_tmp/Result/' + fastq +'.mcool'):
        log_out.info('The result of ' + fastq + ' has been generated')
        # return
 
    else:
        for val in os.listdir(opt['output'] + '/' + opt['type'] + '_' + fastq + '_tmp'):
            print('rm -r ' + val)
            os.system('rm -r ' + val)
            print('clearing done')

    print('********run scMethyl after clearing*********')

    
    fq_r1 = opt['fastq_dir'][opt['fastq_log'].index(fastq)]
    fq_r2 = fq_r1.replace('_1.fastq', '_2.fastq')


    # 创建logging文件
    logging.basicConfig(filename=fastq + '_logging.log', level=logging.DEBUG)
 


    for k in opt.keys():

     
        if (k == 'fastq_dir') or (k == 'fastq_log') or (k == 'run_sample'):
            continue
        else:
            logging.debug('# ' + k + ': ' + str(opt[k]))

    if not os.path.exists('Result'):
        os.makedirs('Result')
    if not os.path.exists('Result/SCPair'):
        os.makedirs('Result/SCPair')
    

    # get index
    t = time.time()
    read_fastq_gz(fq_r1, 'index.txt')
    logging.info('read_fastq_gz time: ' + str(time.time() - t))



    
    # Get paird fastq
    cmd = 'fastp '
    cmd += '--trim_front1 25 --trim_front2 25 --trim_tail1 3 --trim_tail2 3 '
    cmd += ' -i ' + fq_r1 + ' -I ' + fq_r2
    cmd += ' --out1 trimmed-pair1.fastq.gz --out2 trimmed-pair2.fastq.gz'
    cmd += ' --failed_out failed.fq.gz'
    cmd += ' --html trimmed.fastp.html '
    cmd += ' --json trimmed.fastp.json '
    print('1++++++\n', cmd)
    t = time.time()
    logging.info(cmd)
    os.system(cmd)
    logging.info('fastp time: ' + str(time.time() - t))
    
    # Align reads to combo-reference using bwa/bowtie2
    opt['index'] = opt['index'].rsplit('/', 1)[0]

    t = time.time()
    if os.path.exists(opt['index'] + '/Bisulfite_Genome' ) :
        pass
    else:
        logging.info('create bismark index')
        cmd = 'bismark_genome_preparation --bowtie2 ' + opt['index']
        logging.info(cmd)
        os.system(cmd)
        logging.info('bismark_genome_preparation time: ' + str(time.time() - t))
   


    # ####第一次mapping
    cmd = 'bismark --bowtie2 ' + opt['index']
    cmd += ' --fastq trimmed-pair1.fastq.gz --pbat -un --parallel ' + str(opt['thread'])  
    t = time.time()
    logging.info(cmd)
    os.system(cmd)
    logging.info('bismark1 time: ' + str(time.time() - t))

    cmd = 'bismark --bowtie2 ' + opt['index']
    cmd += ' --fastq trimmed-pair2.fastq.gz -un --parallel ' + str(opt['thread']) 
    t = time.time()
    logging.info(cmd)
    os.system(cmd)
    logging.info('bismark2 time: ' + str(time.time() - t))

    # 对umapped的进行剪切继续map
    t = time.time()
    logging.info('split unmapped reads')
    split( 'trimmed-pair1.fastq.gz_unmapped_reads.fq.gz', 40 ,40)
    split( 'trimmed-pair2.fastq.gz_unmapped_reads.fq.gz', 40 ,40)
    logging.info('split time: ' + str(time.time() - t))


    # 对umapped进行mapped
    cmd = 'bismark --bowtie2 ' + opt['index']
    cmd += ' --fastq trimmed-pair1.fastq.gz_unmapped_reads.fq.gz_r1.fq --pbat --parallel ' + str(opt['thread']) 
    t = time.time()
    logging.info(cmd)
    os.system(cmd)
    logging.info('bismark3 time: ' + str(time.time() - t))

    cmd = 'bismark --bowtie2 ' + opt['index']
    cmd += ' --fastq trimmed-pair2.fastq.gz_unmapped_reads.fq.gz_r1.fq  --parallel ' + str(opt['thread']) 
    t = time.time()
    logging.info(cmd)
    os.system(cmd)
    logging.info('bismark4 time: ' + str(time.time() - t))

    # rm
    cmd = 'rm *_unmapped_reads.fq.*'
    logging.info(cmd)
    t = time.time()
    os.system(cmd)
    logging.info('rm time: ' + str(time.time() - t))
    
    # merge
    cmd = 'samtools merge *.bam -o ' + fastq + '.merged.bam'
    logging.info(cmd)
    t = time.time()
    os.system(cmd)
    logging.info('samtools merge time: ' + str(time.time() - t))


    t = time.time()
    qc = opt['qc']
    cmd = 'samtools view -q %s -@ %s -b ' % (qc,str(opt['thread'])) + fastq + '.merged.bam -o ' + fastq + '.filtered.bam'
    logging.info(cmd)
    os.system(cmd)

    cmd = 'samtools sort -@ %s -n ' % str(opt['thread']) + fastq + '.filtered.bam -o ' + fastq + '.sorted.bam'
    logging.info(cmd)
    os.system(cmd)

    cmd = 'samtools view -@ %s -h ' % str(opt['thread']) + fastq + '.sorted.bam -o ' + fastq  + '.sorted.sam'
    logging.info(cmd)
    os.system(cmd)
    logging.info('samtools view time: ' + str(time.time() - t))

    logging.info('filtSam')
    t = time.time()
    filtSam(fastq + '.sorted.sam', fastq + '.sam')
    logging.info('filtSam time: ' + str(time.time() - t))


    # Generate pairs
    cmd = 'pairtools parse ' + fastq + '.sam ' + ' -c ' + opt['genomesize']
    cmd += ' --drop-sam --drop-seq --output-stats ' + fastq + '.stats '
    cmd += ' --assembly ' + opt['species']  + ' --no-flip '
    cmd += ' --add-columns ' + opt['add_columns']
    cmd += ' --walks-policy all '
    cmd += '-o ' + fastq + '.pairs.gz'
    cmd += ' --nproc-in ' + str(opt['thread']) +  ' --nproc-out ' + str(opt['thread'])
    # print('4++++++++++++\n',cmd)
    t = time.time()
    logging.info(cmd)
    os.system(cmd)
    logging.info('pairtools parse time: ' + str(time.time() - t))

    cmd = 'pairtools restrict -f ' + opt['enzyme_bed'] + ' ' + fastq + '.pairs.gz -o ' + fastq + '.restrict.pairs.gz'
    # print('4++++++++++++\n',cmd)
    t = time.time()   
    logging.info(cmd)
    os.system(cmd)
    logging.info('pairtools restrict time: ' + str(time.time() - t))

    os.makedirs('./Result/SCPair', exist_ok=True)
    generate(fastq + '.restrict.pairs.gz', 'index.txt', './Result/SCPair')


    # #QC select
    cmd = 'pairtools select ' + '\"' + opt['select'] + '\"'
    cmd += ' ' + fastq + '.restrict.pairs.gz -o ' + fastq + '.filtered.pairs.gz'
    # print('5++++++++++\n', cmd)
    t = time.time()
    logging.info(cmd)
    os.system(cmd)
    logging.info('pairtools select time: ' + str(time.time() - t))

    # Sort pairs
    cmd ='pairtools sort --nproc ' + str(opt['thread'])
    cmd += ' ' + fastq + '.filtered.pairs.gz -o ' + fastq + '.sorted.pairs.gz'
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
        '.sorted.pairs.gz'
    # print('7+++++++++\n', cmd)
    t = time.time()
    logging.info(cmd)
    subprocess.run(cmd, shell=True, executable="/bin/bash")
    logging.info('pairtools dedup time: ' + str(time.time() - t))
    

    # #QC select
    cmd = 'pairtools select ' + '\"' + opt['select'] + '\"'
    cmd += ' ' + fastq + '.nodups.pairs.gz -o ' + fastq + '.nodups.UU.pairs.gz'
    # print('8++++++++++\n', cmd)
    t = time.time()
    logging.info(cmd)
    os.system(cmd)
    logging.info('pairtools select time: ' + str(time.time() - t))


        # 
    cmd = 'gunzip ' + fastq + '.nodups.UU.pairs.gz'
    # print('9+++++++++++\n',cmd)
    t = time.time()
    logging.info(cmd)
    os.system(cmd)
    logging.info('gunzip time: ' + str(time.time() - t))
    
    cmd = 'cooler cload pairs -c1 2 -p1 3 -c2 4 -p2 5 '
    cmd += opt['genomesize'] + ':' + str(opt['resolution']) + ' ' + fastq + '.nodups.UU.pairs ' + fastq + str(opt['resolution']) + '.cool'
    # print('10+++++++++++\n', cmd)
    t = time.time()
    logging.info(cmd)
    os.system(cmd)
    logging.info('cooler cload time: ' + str(time.time() - t))
   
    cmd = 'gzip ' + fastq + '.nodups.UU.pairs'
    # print('11+++++++++++\n', cmd)
    t = time.time()
    logging.info(cmd)
    os.system(cmd)
    logging.info('gzip time: ' + str(time.time() - t))

    

    cmd = 'cooler zoomify '  + fastq + str(opt['resolution']) + '.cool -r %s -o ' % opt['zoomify_res'] + fastq + '.mcool'
    # print('13+++++++++++\n', cmd)
    t = time.time()
    logging.info(cmd)
    os.system(cmd)
    logging.info('cooler zoomify time: ' + str(time.time() - t))


    os.system('mv ' + fastq + '.mcool ./Result/')
    os.system('mv ' + fastq + str(opt['resolution']) + '.KR.cool ./Result/')
    os.system('mv ' + fastq + str(opt['resolution']) + '.cool ./Result/')
    os.system('mv ' + fastq + '.nodups.pairs.gz ./Result/')

    pair_list = []
    for p in os.listdir('./Result/SCPair'):
        if p.endswith('pairs.gz'):
            pair_list.append(p)
    Parallel(n_jobs=6)(delayed(cload_correct_zoomify)('./Result/SCPair/',
                                      p,
                                      opt['genomesize'],
                                      str(opt['resolution']),
                                       opt['zoomify_res']
                            ) for p in pair_list)



    logging.info('total time: ' + str(time.time() - T))




def read_fastq_gz(file_path, index_file):
    fp = open(index_file, "w")
    with gzip.open(file_path, "rt") as handle:
        for record in SeqIO.parse(handle, "fastq"):
            fp.write(str(record.id) + "\t" + str(record.seq[0:6]) + "\n")
    fp.close()
    return

def split(fn1, size1, size2):
    if 'gz' in fn1:
            dfh1=gzip.open(fn1,'r')

    if 'gz' not in fn1:
        dfh1=open(fn1,'r')
    if sys.version_info[0]==3:
        ID1=str(dfh1.readline().rstrip())[2:-1]
        line1=str(dfh1.readline().rstrip())[2:-1]
        plus1=str(dfh1.readline().rstrip())[2:-1]
        QS1=str(dfh1.readline().rstrip())[2:-1]
    if sys.version_info[0]==2:
        ID1=dfh1.readline().rstrip()
        line1=dfh1.readline().rstrip()
        plus1=dfh1.readline().rstrip()
        QS1=dfh1.readline().rstrip()
    rfhs=[]


    trim1=5
    trim2=5

    for i in range(1,4):
        rfhs.append(open(fn1+'_r'+str(i)+'.fq','w'))


    while line1:
        if len(line1[trim1:size1+trim1])>=30:
            rfhs[0].write(str(ID1)+'-1'+'\n'+str(line1[trim1:size1+trim1])+'\n'+str(plus1[0])+'\n'+str(QS1[trim1:size1+trim1])+'\n')
        if len(line1[trim1+size1:(-1*size2)-trim2])>=30:
            rfhs[0].write(str(ID1)+'-2'+'\n'+str(line1[trim1+size1:(-1*size2)-trim2])+'\n'+str(plus1[0])+'\n'+str(QS1[trim1+size1:(-1*size2)-trim2])+'\n')
        if len(line1[(-1*size2)-trim2:-1*trim2])>=30:
            rfhs[0].write(str(ID1)+'-3'+'\n'+str(line1[(-1*size2)-trim2:-1*trim2])+'\n'+str(plus1[0])+'\n'+str(QS1[(-1*size2)-trim2:-1*trim2])+'\n')
        if sys.version_info[0]==3:
            ID1=str(dfh1.readline().rstrip())[2:-1]
            line1=str(dfh1.readline().rstrip())[2:-1]
            plus1=str(dfh1.readline().rstrip())[2:-1]
            QS1=str(dfh1.readline().rstrip())[2:-1]
        if sys.version_info[0]==2:
            ID1=dfh1.readline().rstrip()
            line1=dfh1.readline().rstrip()
            plus1=dfh1.readline().rstrip()
            QS1=dfh1.readline().rstrip()  

    map(lambda x:x.close(),rfhs)


def filtSam(SamFile, CleanFile ):
  
    fp = open(CleanFile, "w")
    Read_list = []
    ReadID_list = []
    with open(SamFile, 'r') as file:
        lines = file.readlines()

    for i in range(len(lines)-1):
        if "@" in lines[i]:
            #print(lines[i])
            #print(i)
            fp.write(lines[i])

        else:
            Read_ID1 = lines[i].split("\t")[0].split("_")[0]
            Read_ID2 = lines[i+1].split("\t")[0].split("_")[0]
            components_1 = lines[i].split("\t")
            components_2 = lines[i+1].split("\t")

            if Read_ID1 == Read_ID2:
                components_1[0] = Read_ID1
                components_2[0] = Read_ID2
                lines1 = "\t".join(components_1)
                lines2 = "\t".join(components_2)

                fp.write(lines1+lines2)
                lines[i+1] = ''
    fp.close()

def generate(pairFile,  ReadIndex,ResultDir):
    ReadDict, SClist = GenerateSCDict(ReadIndex)
    GenerateSCPairs(ReadDict, SClist, pairFile, ResultDir)


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
        # if SClist[i] not in ['ACTTGA', 'CGATGT', 'GCCAAT', 'TAGCTT']:
        #     continue
        if len(SCPairs[SClist[i]])>1:
            scPairsSummary.write(str(i+1)+"\t"+str(SClist[i])+"\t"+str(len(SCPairs[SClist[i]]))+"\n")
        # if len(SCPairs[SClist[i]])<1000:
        #     continue
        cell_count+=1
        print("Process cell %d,\t"%(cell_count),end='')
        if not os.path.exists(ResultDir):
            os.mkdir(ResultDir)
        fp = gzip.open(ResultDir+"/cell_"+str(cell_count)+"_"+str(i)+".pairs.gz","w+")
        for singel_head in HeadInfo:
            fp.write(singel_head.encode('utf-8'))
        count=0

        for pair in SCPairs[SClist[i]]:
            fp.write(pair.encode('utf-8'))
            count+=1
        print("%d pairs"%(count))
        fp.close()
    scPairsSummary.close()