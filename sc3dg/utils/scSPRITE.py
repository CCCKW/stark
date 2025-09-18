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

from multiprocessing import Pool

# jar_dir = sys.prefix + '/inHic'
# jar_path = jar_dir + '/utils/BarcodeIdentification_v1.2.0.jar'
jar_path = '/cluster/home/Kangwen/Hic/seq/4_scSPRITE/Code/BarcodeIdentification_v1.2.0.jar'


def pp(opt, fastq, log_out):
    print('====================scSPRITE====================')
    T = time.time()
    # 移动到工作目录
    os.chdir(opt['output'])

    # 创建文件夹tmp,并移动
    if not os.path.exists(opt['type'] + '_' + fastq + '_tmp'):  # scHic_sample_tmp
        os.makedirs(opt['type'] + '_' + fastq + '_tmp')
    output = opt['output'] + '/' + opt['type'] + '_' + fastq + '_tmp'
    os.chdir(opt['type'] + '_' + fastq + '_tmp')

    fq_r1 = opt['fastq_dir'][opt['fastq_log'].index(fastq)]
    fq_r2 = fq_r1.replace('_1.fastq', '_2.fastq')

    # 创建log
    logging.basicConfig(filename=fastq + '_logging.log', level=logging.DEBUG)
    for k in opt.keys():

        if (k == 'fastq_dir') or (k == 'fastq_log') or (k == 'run_sample'):
            continue
        else:
            logging.debug('# ' + k + ': ' + str(opt[k]))

    os.makedirs('./Result', exist_ok=True)

    os.makedirs('./Result/SCpair', exist_ok=True)
    os.makedirs('./Result/mcool_folder/', exist_ok=True)
    os.makedirs('./Result/cool_folder/', exist_ok=True)
    for val in os.listdir(opt['output'] + '/' + opt['type'] + '_' + fastq + '_tmp/Result/mcool_folder/'):
        if val.endswith('mcool'):
            pass
            # print('==========generated+++++++')
            # return

    # Generate Barcode Reads
    cmd = 'java -jar ' + jar_path
    cmd += ' ' + '--input1 ' + fq_r1 + ' --input2 ' + fq_r2
    cmd += ' ' + '--output1 ' + fastq + '_1_barcode.fastq.gz'
    cmd += ' ' + '--output2  ' + fastq + '_2_barcode.fastq.gz'
    cmd += ' --config ' + opt['sprite_config']
    t = time.time()
    logging.info(cmd)
    # os.system(cmd)
    logging.info('BarcodeIdentification time: ' + str(time.time() - t))

    # Get paird fastq
    cmd = 'fastp --overrepresentation_analysis'
    cmd += ' --correction --thread ' + str(opt['thread'])
    cmd += ' --html trimmed.fastp.html '
    cmd += ' --json trimmed.fastp.json '
    cmd += ' -i ' + fastq + '_1_barcode.fastq.gz' + ' -I ' + fastq + '_2_barcode.fastq.gz'
    cmd += ' --out1 trimmed-pair1.fastq.gz --out2 trimmed-pair2.fastq.gz'
    cmd += ' --failed_out failed.fq.gz'
    cmd += ' --max_len1 100'
    cmd += ' --max_len2 100'
    t = time.time()
    logging.info(cmd)
    os.system(cmd)
    logging.info('fastp time: ' + str(time.time() - t))

    t = time.time()
    logging.info('filter_reads')

    filter_reads('trimmed-pair1.fastq.gz', 'trimmed-pair1.clean.fastq.gz')
    filter_reads('trimmed-pair2.fastq.gz', 'trimmed-pair2.clean.fastq.gz')
    logging.info('filter_reads time: ' + str(time.time() - t))

    # Align reads to combo-reference using bwa/bowtie2
    t = time.time()
    if opt['aligner'] == 'bwa':
        cmd = 'bwa mem ' + opt['index']
        cmd += ' -t ' + str(opt['thread'])
        cmd += ' trimmed-pair1.clean.fastq.gz trimmed-pair2.clean.fastq.gz > ' + fastq + '.sam'
        print('6+++++++\n', cmd)
        os.system(cmd)
    else:
        cmd = 'bowtie2 -p ' + str(opt['thread'])
        cmd += ' -x ' + opt['index']
        cmd += ' -1 trimmed-pair1.clean.fastq.gz -2 trimmed-pair2.clean.fastq.gz -S ' + fastq + '.sam'
        print('6+++++++\n', cmd)
        os.system(cmd)
    logging.info(cmd)
    logging.info('align time: ' + str(time.time() - t))

    # sam to bam
    cmd = 'samtools view -@ 8 -S ' + fastq + '.sam -b -o ' + fastq + '.bam'
    t = time.time()
    logging.info(cmd)
    os.system(cmd)
    logging.info('sam to bam time: ' + str(time.time() - t))

    qc = opt['qc']

    print(123)

    cmd = 'samtools view -h -@ ' + str(opt['thread'])
    if qc == 0:
        cmd += ' -b ' + fastq + '.bam > ' + fastq + '.filtered.bam'
    else:
        cmd += ' -b -q  %s ' % qc + fastq + '.bam > ' + fastq + '.filtered.bam'
    t = time.time()
    logging.info(cmd)
    os.system(cmd)
    logging.info('samtools view time: ' + str(time.time() - t))

    cmd = 'bedtools intersect -v -a ' + fastq + '.filtered.bam'
    cmd += ' -b ' + opt['repeat_masker'] + ' > ' + fastq + '.filtered.masked.bam'
    t = time.time()
    logging.info(cmd)
    os.system(cmd)
    logging.info('bedtools intersect time: ' + str(time.time() - t))

    cmd = 'samtools view -h -@ ' + str(opt['thread'])
    cmd += ' -o ' + fastq + '.filtered.masked.sam ' + fastq + '.filtered.masked.bam'
    t = time.time()
    logging.info(cmd)
    os.system(cmd)
    logging.info('samtools view time: ' + str(time.time() - t))

    t = time.time()
    logging.info('GenerateInfo')
    GenerateInfo(fastq + '.filtered.masked.sam')
    logging.info('GenerateInfo time: ' + str(time.time() - t))

    cmd = 'sort --parallel=16 -k2,2 -k3,3 '
    cmd += ' ScSpriteInfo.txt > ScSpriteInfo.sorted.txt'
    t = time.time()
    logging.info(cmd)
    os.system(cmd)
    logging.info('sort time: ' + str(time.time() - t))

    t = time.time()
    WriteSCPairCluster('ScSpriteInfo.sorted.txt')
    logging.info('WriteSCPairCluster time: ' + str(time.time() - t))
    t = time.time()
    WriteCellPairs(opt['genomesize'],
                   'HeadInfo.txt',
                   opt['species'], 'SCPairCellSummary.txt',
                   'SCPairClusterSummary.txt', output)
    logging.info('WriteCellPairs time: ' + str(time.time() - t))

    mcool_folfer = './Result/mcool_folder'
    os.makedirs('./Result/mcool_folder', exist_ok=True)

    # for val in os.listdir('./Result/SCpair'):
    #     if val.endswith('.pairs.gz') and '[' in val and ']' in val:
    #         val = val.split('.')[0]
    #         cmd = 'cooler cload pairs -c1 2 -p1 3 -c2 4 -p2 5 '
    #         cmd += opt['genomesize'] + ':' + str(opt['resolution']) + ' ./Result/SCpair/' + val + '.pairs.gz ./Result/SCpair/' + val  + '_'+str(opt['resolution']) + '.cool'
    #         t = time.time()
    #         logging.info(cmd)
    #         print(cmd)
    #         os.system(cmd)
    #         logging.info('cooler cload time: ' + str(time.time() - t))

    #         cmd = 'hicCorrectMatrix correct --matrix  ./Result/SCpair/' + val + '_' +   str(opt['resolution']) + '.cool '
    #         cmd += ' --correctionMethod KR --outFileName '
    #         cmd += ' ./Result/SCpair/' + val + str(opt['resolution']) + '.KR.cool'
    #         cmd += ' --filterThreshold -1.5 5.0 --chromosomes chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19'
    #         print('14+++++++++++\n', cmd)
    #         t = time.time()
    #         logging.info(cmd)
    #         os.system(cmd)
    #         logging.info('hicCorrectMatrix time: ' + str(time.time() - t))

    #         cmd = 'cooler zoomify ./Result/SCpair/'  + val + str(opt['resolution']) + '.KR.cool -r 10000,40000,100000,500000 -o ./Result/mcool_folder/' + val + '.mcool'
    #         print('15+++++++++++\n', cmd)
    #         t = time.time()
    #         logging.info(cmd)
    #         os.system(cmd)
    #         logging.info('cooler zoomify time: ' + str(time.time() - t))

    os.system('mv ./Resul/SCpair/*.coowl ./Result/cool_folder/')
    os.system('mv ./Resul/SCpair/*.mcool ./Result/mcool_folder/')

    logging.info('Total time: ' + str(time.time() - T))


def WriteSCPairCluster(ScSpriteInfo):
    PairDir = {}
    currentCell = 'None'
    SCPairClusterSummary = open("SCPairClusterSummary.txt", "wt")
    SCPairCellSummary = open("SCPairCellSummary.txt", "wt")

    count = 0
    CellPairCount = {}
    lineCount = 0
    StartTime = time.time()
    for line in open("ScSpriteInfo.sorted.txt"):

        lineCount += 1
        if lineCount % 1000000 == 0:
            print(f"Process {lineCount} Reads in {time.time() - StartTime}s ")
            StartTime = time.time()

        ReadComponents = line.split("\t")

        if currentCell != ReadComponents[1]:
            CellPairCount[currentCell] = [count, 0]
            count = 0
            currentCell = ReadComponents[1]

            if len(PairDir) > 2000:
                for key in PairDir.keys():

                    split_result = key.split('][')

                    # 获取最后符合要求的部分
                    last_three_brackets = '[' + ']['.join(split_result[-3:])

                    if CellPairCount[last_three_brackets][1] == 0:
                        SCPairCellSummary.write(
                            last_three_brackets + "\t" + str(CellPairCount[last_three_brackets][0]) + "\n")
                        CellPairCount[last_three_brackets][1] = 1

                    # if CellPairCount[last_three_brackets][0] < 1000:
                    #     continue

                    Locs = list(PairDir[key])
                    if len(Locs) < 2 or len(Locs) > 10000:
                        continue

                    flag = 0
                    for SingleTerm in Locs:
                        if flag == 0:
                            flag += 1
                            SCPairClusterSummary.write(key)
                        SCPairClusterSummary.write("\t" + SingleTerm)
                    SCPairClusterSummary.write("\n")

                PairDir = {}
                CellPairCount = {}

        #####判断正负链
        if bool(int(ReadComponents[3]) & 0x10):
            ReadComponents[3] = '-'
        else:
            ReadComponents[3] = '+'

        if ReadComponents[2] not in PairDir:
            PairDir[ReadComponents[2]] = set()
            PairDir[ReadComponents[2]].add(
                ReadComponents[4] + ":" + ReadComponents[5] + ":" + ReadComponents[3] + ":" + ReadComponents[6])
            count += 1
        else:
            PairDir[ReadComponents[2]].add(
                ReadComponents[4] + ":" + ReadComponents[5] + ":" + ReadComponents[3] + ":" + ReadComponents[6])
            count += 1

    SCPairClusterSummary.close()
    SCPairCellSummary.close()


def WriteCellPairs(ChromSize, HeadInfo, GenomeInfo, SCPairCellSummary, SCPairClusterSummary, output):
    df = pd.read_csv(SCPairCellSummary, delimiter='\t', header=None)

    sorted_df = df.sort_values(by=df.columns[1], ascending=False)

    CellList = list(sorted_df.iloc[:, 0])
    ProcessedCell = None
    for line in open(SCPairClusterSummary):

        ReadComponents = line[0:-1].split("\t")

        split_result = ReadComponents[0].split('][')
        # 获取最后符合要求的部分
        CurrentCellID = '[' + ']['.join(split_result[-3:])
        if CurrentCellID not in CellList:
            continue

        if ProcessedCell != CurrentCellID:

            if ProcessedCell != None:
                print(ProcessedCell)
                fp.close()

            ProcessedCell = CurrentCellID
            print(output + '/Result/SCpair/' + str(CurrentCellID) + '.pairs.gz')
            fp = gzip.open(output + '/Result/SCpair/' + str(CurrentCellID) + '.pairs.gz', "w+")

            fp.write(("## pairs format v1.0.0\n#shape: whole matrix\n#genome_assembly: " + GenomeInfo + "\n").encode(
                'utf-8'))
            for line in open(ChromSize):
                fp.write(("#chromsize: " + line.split("\t")[0] + " " + line.split("\t")[1]).encode('utf-8'))

            for line in open(HeadInfo):
                fp.write(("#samheader: " + line).encode('utf-8'))
            fp.write("#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type mapq1 mapq2\n".encode('utf-8'))

        for i in range(1, len(ReadComponents)):
            for j in range(1, i):
                Pair1 = ReadComponents[i].split(":")
                Pair2 = ReadComponents[j].split(":")
                fp.write((ReadComponents[0] + "\t" + \
                          Pair1[0] + "\t" + \
                          Pair1[1] + "\t" + \
                          Pair2[0] + "\t" + \
                          Pair2[1] + "\t" + \
                          Pair1[2] + "\t" + \
                          Pair1[2] + "\t" + \
                          "UU" + "\t" + \
                          Pair1[3] + "\t" + \
                          Pair1[3] + "\n").encode('utf-8'))
    fp.close()


def filter_reads(input_fastq, output_fastq):
    fastq = gzip.open(output_fastq, 'w')
    with gzip.open(input_fastq, "rt") as fp:
        flag = 0
        for line in fp:
            # line = line.decode()
            if '@SRR' in line and '[NOT_FOUND]' not in line:
                flag = 3
                fastq.write(line.encode('utf-8'))
                continue

            while flag != 0:
                flag -= 1
                fastq.write(line.encode('utf-8'))
                break
    fastq.close()


def GenerateInfo(SAMFile):
    HeadInfo = open("HeadInfo.txt", "w")
    ReadsInfo = open("ScSpriteInfo.txt", "w")
    count = 0
    startTime = time.time()
    with open(SAMFile, "r") as pf:
        Info = pf.readlines()
        for line in Info:
            if line.startswith("@"):
                HeadInfo.write(line)
                continue
            ReadComponent = line.split("\t")
            ReadID = ReadComponent[0].split("::")[0]
            ClusterID = ReadComponent[0].split("::")[1]
            CellID = "[" + "][".join(ClusterID[1:-1].split("][")[3:]) + "]"

            ReadsInfo.write(ReadID + "\t" + CellID + "\t" + ClusterID)
            for i in range(1, len(ReadComponent)):
                ReadsInfo.write("\t" + ReadComponent[i])

            count += 1
            if count % 1000000 == 0:
                print("%d reads has processed in %f s" % (count, time.time() - startTime))
                startTime = time.time()
        HeadInfo.close()
        ReadsInfo.close()