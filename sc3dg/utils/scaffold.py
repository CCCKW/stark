import os
import sys
import subprocess
import time
import re
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio import SeqIO
from Levenshtein import hamming
import gzip
import logging
# jar_dir = sys.prefix + '/inHic' 
# jar_path = jar_dir + '/utils/BarcodeIdentification_v1.2.0.jar'
jar_path = '/cluster/home/Kangwen/Hic/seq/4_scSPRITE/Code/BarcodeIdentification_v1.2.0.jar'



def run_cmd_return_time(cmd):
    t = time.time()
    logging.debug(cmd)
    print(cmd)
    subprocess.run(cmd, shell=True, executable="/bin/bash")
    return time.time() - t


def sprite_filter_reads(infastq,outfastq):
    fastq = gzip.open(outfastq,'w+')
    with gzip.open(infastq, "rt") as fp:
        flag = 0
        for line in fp:
            #line = line.decode()
            if '@' in line and '[NOT_FOUND]' not in line and ']' in line[-2]:
                flag = 3
                fastq.write(line.encode('utf-8'))
                continue

            while flag !=0:
                flag-=1
                fastq.write(line.encode('utf-8'))
                break
    fastq.close()

def fastp(fq1,fq2,o1,o2, threads,max=True,adapter:list=None):
    cmd = 'fastp --overrepresentation_analysis'
    if adapter is not None:
        cmd += ' --adapter_sequence ' + adapter[0]
        cmd += ' --adapter_sequence_r2 ' + adapter[1]
    cmd += ' --correction --thread ' + str(threads)
    cmd += ' --html trimmed.fastp.html '
    cmd += ' --json trimmed.fastp.json '
    cmd += ' -i ' + fq1 + ' -I ' + fq2
    cmd += ' --out1 %s --out2 %s' % (o1, o2)
    cmd += ' --failed_out failed.fq.gz'
    if max:
        cmd += ' --max_len1 100'
        cmd += ' --max_len2 100'
    return run_cmd_return_time(cmd)

def samtools_qc(threads,qc,fastq):
    if qc == 0:
         cmd = 'cp %s.bam %s.qc.bam' % (fastq,fastq)
    else:
        cmd = 'samtools view -h -@ %s -b -q %s %s.bam > %s.qc.bam' % (threads, qc, fastq,fastq)
    
    return run_cmd_return_time(cmd)

def bwa(fa, threads,fastq):
    cmd = 'bwa mem ' + fa
    cmd += ' -t ' + str(threads)
    cmd += ' trimmed-pair1.fastq.gz trimmed-pair2.fastq.gz > '  + fastq + '.sam'
    return run_cmd_return_time(cmd)

def bowtie2(fa, threads,fastq):
    cmd = 'bowtie2 -p ' + str(threads)
    cmd += ' -x ' + fa
    cmd += ' -1 trimmed-pair1.fastq.gz -2 trimmed-pair2.fastq.gz -S ' + fastq + '.sam'
    return run_cmd_return_time(cmd)

def sam_to_bam(fastq, threads):
    
    cmd = 'samtools view -@ ' + str(threads) + '-S ' + fastq + '.sam -b -o ' + fastq + '.bam'
    return run_cmd_return_time(cmd)

def pairtools_parse(fastq, genomesize, species, add_columns, threads):

    cmd = 'pairtools parse ' + fastq + '.qc.bam ' + ' -c ' + genomesize
    cmd += ' --drop-sam --drop-seq --output-stats ' + fastq + '.stats '
    cmd += ' --assembly ' + species  + ' --no-flip '
    cmd += ' --add-columns ' + add_columns
    cmd += ' --walks-policy all '
    cmd += '-o ' + fastq + '.pairs.gz'
    cmd += ' --nproc-in ' + str(threads) +  ' --nproc-out ' + str(threads)
    
    return run_cmd_return_time(cmd)

def pairtools_restrict(srr, enzyme):
    cmd = 'pairtools restrict -f %s %s.pairs.gz -o %s.restrict.pairs.gz' % (enzyme, srr, srr)
    return run_cmd_return_time(cmd)

def pairtools_select(srr, select):
    cmd = 'pairtools select ' + '\"' + select + '\"'
    cmd += ' ' + srr + '.restrict.pairs.gz -o ' + srr+ '.select.pairs.gz'
    return run_cmd_return_time(cmd)

def pairtools_sort(fastq):
    cmd = 'pairtools sort ' + fastq + '.select.pairs.gz -o ' + fastq +'.sort.pairs.gz'
    return run_cmd_return_time(cmd)

def pairtools_dedup(fastq, max_mismatch):
    cmd = 'pairtools dedup '
    cmd += '--max-mismatch ' + max_mismatch
    cmd += ' --mark-dups --output ' + '>( pairtools split --output-pairs ' + fastq + \
        '.dedup.pairs.gz --output-sam ' + fastq + \
        '.nodups.bam ) --output-unmapped >( pairtools split --output-pairs  ' + fastq + \
        '.unmapped.pairs.gz --output-sam  ' + fastq + \
        '.unmapped.bam ) --output-dups >( pairtools split --output-pairs  ' + fastq + \
        '.dups.pairs.gz --output-sam ' + fastq + \
        '.dups.bam ) --output-stats  ' + fastq + \
        '.dedup.stats ' + fastq + \
        '.sort.pairs.gz'

    return run_cmd_return_time(cmd)

def cooler_cload_pairs(genomesize, res, fastq, prefix=''):
    cmd = 'cooler cload pairs -c1 2 -p1 3 -c2 4 -p2 5 '
    if prefix is not  None:
        cmd += genomesize + ':' + str(res) + ' ' + fastq + '.%s.pairs.gz '%prefix+ fastq + '_' + str(res) + '.cool'
    else:
        cmd += genomesize + ':' + str(res) + ' ' + fastq + '.pairs.gz '+ fastq + '_' + str(res) + '.cool'
    return run_cmd_return_time(cmd)



def cooelr_zoomify(fastq,res, resz):
    cmd = 'cooler zoomify '  + \
        fastq + '_' +str(res) + \
        '.cool -r ' + resz + ' -o ' + fastq + '.mcool'
    return run_cmd_return_time(cmd)

def find_substring_in_long_string(long_string, substrings):
    long_string = str(long_string)
    pattern = re.compile('|'.join(map(re.escape, substrings)))
    match = pattern.search(long_string)
    if match:
        return match.group()
    else:
        return None

def mv_Result():
    print(os.getcwd())
    for val in os.listdir('./'):
        if val.endswith('.mcool') or val.endswith('.cool'):
            os.system('mv %s Result/' % val)
        if val.endswith('.dedup.pairs.gz') or val.endswith('.restrict.pairs.gz'):
            os.system('mv %s Result/' % val)

def GetFastqIndex(srr ,Fastq1_read, Fastq2_read, Barcode1_read, Barcode2_read):
    t = time.time()
 
    Index_name = 'barcode.index'
    All_reads_count = 0
    Barcoded_Reads_count = 0
    MisBarcode_reads_count = 0
   
    fastq1_barcode = gzip.open(srr + "._1_barcode.fastq.gz", "w")
    fastq2_barcode = gzip.open(srr + "._2_barcode.fastq.gz", "w")
    Index_File = open(Index_name, "w+")

    Fastq1_reads = gzip.open(Fastq1_read, "rt")
    Fastq2_reads = gzip.open(Fastq2_read, "rt") 
    Barcode1_reads =  gzip.open(Barcode1_read, "rt")
    Barcode2_reads = gzip.open(Barcode2_read, "rt") 


    for record1, Barcode1, record2, Barcode2 in zip(SeqIO.parse(Fastq1_reads, "fastq"),
                                                    SeqIO.parse(Barcode1_reads, "fastq"),
                                                    SeqIO.parse(Fastq2_reads, "fastq"),
                                                    SeqIO.parse(Barcode2_reads, "fastq")):
  
        # 在这里，你可以使用 record1、Barcode1、record2、Barcode2 来处理每个文件的记录
        All_reads_count += 1
    
        ####Check Read ID
        if (record1.id != record2.id) or (Barcode1.id != Barcode2.id) or (record1.id != Barcode2.id):
            print("The Fastq Files are not matched")

        ####Check and Select barcode
        if len(str(Barcode1.seq)) != 8 or len(str(Barcode2.seq)) != 8 or "N" in str(Barcode1.seq) or "N" in str(Barcode2.seq):
            MisBarcode_reads_count+=1
            continue
        # print(1)
        Barcoded_Reads_count+=1
        ### Write Read ID
        fastq1_barcode.write(("@"+str(record1.id)+"\t"+str(Barcode1.seq)+"_"+str(Barcode2.seq)+"\n").encode("utf-8"))
        fastq2_barcode.write(("@"+str(record2.id)+"\t"+str(Barcode1.seq)+"_"+str(Barcode2.seq)+"\n").encode("utf-8"))
        Index_File.write((str(record2.id)+"\t"+str(Barcode1.seq)+"_"+str(Barcode2.seq)+"\n"))

        ### Write Seq
        fastq1_barcode.write((str(record1.seq)+"\n").encode("utf-8"))
        fastq2_barcode.write((str(record2.seq)+"\n").encode("utf-8"))

        ### Write append Info
        fastq1_barcode.write(("+" + str(record1.id)+"\n").encode("utf-8"))
        fastq2_barcode.write(("+" + str(record2.id)+"\n").encode("utf-8"))

        ### Write quality
        Read1_Quality = ''.join([chr(phred + 33) for phred in record1.letter_annotations["phred_quality"]])
        Read2_Quality = ''.join([chr(phred + 33) for phred in record2.letter_annotations["phred_quality"]])
        fastq1_barcode.write((str(Read1_Quality)+"\n").encode("utf-8"))
        fastq2_barcode.write((str(Read2_Quality)+"\n").encode("utf-8"))

    fastq1_barcode.close()
    fastq2_barcode.close()
    Index_File.close()
    
    print("All_reads_count: ", All_reads_count)
    print("Barcoded_Reads_count: ", Barcoded_Reads_count)
    print("MisBarcode_reads_count", MisBarcode_reads_count)
    return time.time() - t

def digest(size,fa,  enzyme):
    if ' ' in enzyme:
        enzyme = enzyme.replace(' ', '\ ')
    cmd = 'cooler digest  -o digest.bed %s %s %s' % (size, fa, enzyme)
    print(cmd)
    return run_cmd_return_time(cmd)

def split_cells(fastq):
    print(fastq)
    pair = fastq +   '.restrict.pairs.gz'
    ReadDict, SClist = GenerateSCDict('barcode.index')
    GenerateSCPairs(ReadDict, SClist, pair, './Result/SCpair')

def split_fastqs(r1, r2, r1_o, r2_o, barcodes):
	'''Given paired-end fastq data, split reads based off of an inline 10 (or 11) mer barcode'''
	mismatch = 0
	not_found_R1 = 0
	not_found_R2 = 0
	reads = 0
	#fqr1 = FastqGeneralIterator(r1)
	#fqr2 = FastqGeneralIterator(r2)
	#seqzip = zip(fqr1, fqr2) #Zip up the two iterators for expediency

	with gzip.open(r1, "rt") as handle1:
		with gzip.open(r2, "rt") as handle2:
			for (title1, seq1, qual1), (title2, seq2, qual2) in zip(FastqGeneralIterator(handle1),FastqGeneralIterator(handle2)):
				barcode1 = seq1[:8] #Just look in read 1--barcodes SHOULD be the same on both ends of the molecule)
				barcode2 = seq2[:8] #Check barcode 2 as well; print out how many times they disagree
				test1 = checkHamming(barcodes, barcode1)
				test2 = checkHamming(barcodes, barcode2)
				if test1[0]:
					if test2[0]:
					#if the barcodes match, print out the trimmed / split reads to new files
						if test1[1] == test2[1]:
							print("@%s&%s\n%s\n+\n%s" % (title1,test1[1], seq1[11:], qual1[11:]), file=r1_o)
							print("@%s&%s\n%s\n+\n%s" % (title2,test2[1], seq2[11:], qual2[11:]), file=r2_o)
						else:
							mismatch += 1
					else: #If there isn't a match in R1
						not_found_R2 += 1
				elif test2[0]:
					not_found_R1 += 1
				else:
					not_found_R1 += 1
					not_found_R2 += 1
				reads += 1
	print("mismatch: ",mismatch)
	print("not_found_R1: ",not_found_R1)
	print("not_found_R2: ",not_found_R2)
	print("reads: ",reads)
	handle1.close()
	handle2.close()
	return mismatch, not_found_R1, not_found_R2, reads

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

def awk():
    cmd = 'awk -F\'\\t\' \'{print $1 \"\\t\"  $3 \"_\" $4}\''
    cmd +=  ' '  'read_barcode.index > barcode.index' 
    print(cmd) 
    return run_cmd_return_time(cmd)

def samtools_merge(fastq):
    cmd = 'samtools merge *.bam -o ' + fastq + '.merged.bam -f'
    time1 = run_cmd_return_time(cmd)
    cmd = 'rm *.bam'
    time2 = run_cmd_return_time(cmd)
    cmd = 'mv ' + fastq + '.merged.bam ' + fastq + '.bam'
    time3 = run_cmd_return_time(cmd)
    return time1 + time2 + time3  

def Generate_snHIC_ReadsIndex():
    file1_path = 'trimmed-pair1.fastq.gz'
    file2_path = 'trimmed-pair2.fastq.gz'
    indexName = 'barcode.index'
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
    DirtyFile.close()


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

def addWildcards(string):
    foo= ""
    for i in range(len(string)):
            foo+=string[:i] + '.' + string[i+1:] + "|"
    return foo[:-1]

def rc(seq):

    rcd = ""
    rcs = {"A":"T","C":"G","T":"A","G":"C", "N":"N"}
    for char in seq[::-1]:
        rcd+=rcs[char]
    return rcd

def associateBarcodes(fhi_r1, fhi_r2, fho_r1, fho_r2, internal_bc):
    outfile = open('read_barcode.index', 'w')
    '''Takes a Fastq filehandle in and filehandle out for each read; writes barcode-clipped fastqs to\
    filehandle out and returns a dictionary with barcode pairs associated with given read names'''
    #fhi_r1 =  gzip.open(fhi_r1, "rt")
    #fhi_r2 = gzip.open(fhi_r2, "rt")
    seq_r1 = FastqGeneralIterator(fhi_r1)
    seq_r2 = FastqGeneralIterator(fhi_r2)
    seqzip = zip(seq_r1, seq_r2) #zip up the two fastq iterators to increase performance
    bc_seq = {}
    bc1s = [] #Container to store internal barcode IDs
    bridge = ''
    #Create a list of all possible internal barcodes (necessary
    #for hamming distance calculations)
    for line in internal_bc:
        bc = line.rstrip('\n')
        strip = 'GATC' + bc + 'GAATTC'
        bc1s.append(bc)
    bridge_re = re.compile('GATC[ACGT]{8}GAATTC')
    #We want to find all reads where
    #the internal adaptor is present
    #in either R1 or R2. Accepting that we're
    #throwing out some stuff because of mismatches,
    ##let's just use the re module to do this.
    for i in seqzip:
        title1, seq1, qual1 = i[0]
        title2, seq2, qual2 = i[1]
        new_title = title1.split()[0]
        m1 = bridge_re.search(seq1)
        m2 = bridge_re.search(seq2)
        terminal_bc = title1.split('&')[1]
        if m1:
            if m2:
                pos1 = m1.start()
                pos2 = m2.start()
                trim_seq1 = seq1[:pos1]
                trim_seq2 = seq2[:pos2]
                trim_qual1 = qual1[:pos1]
                trim_qual2 = qual2[:pos2]
                #Check for length of the read
                if len(trim_seq1) < 30 or len(trim_seq2) < 30:
                    continue #and continue if either one is too short
                barcode_seq1 = seq1[pos1 + 4:pos1 + 12]
                barcode_seq2 = seq2[pos2 + 4:pos2 + 12]
                ham_bc1 = checkHamming(bc1s, barcode_seq1)
                ham_bc2 = checkHamming(bc1s, barcode_seq2)
                if ham_bc1[0]:
                    if ham_bc2[0]:
                        print("@%s\n%s\n+\n%s" % (new_title, trim_seq1, trim_qual1), file=fho_r1)
                        print("@%s\n%s\n+\n%s" % (new_title, trim_seq2, trim_qual2), file=fho_r2)
                        bc1 = ham_bc1[1]
                        bc2 = ham_bc2[1]
                        print("%s\t%s\t%s\t%s" % (new_title, bc1, bc2, terminal_bc), file=outfile)
            else:
                pos1 = m1.start()
                pos1_end = m1.end()
                trim_seq1 = seq1[:pos1]
                trim_qual1 = qual1[:pos1]
                #Check for length of the rea
                if len(trim_seq1) < 30:
                    continue #and continue if either one is too short
                if len(seq1[pos1_end:]) < 12:
                    continue
                barcode_seq2 = rc(seq1[pos1_end:pos1_end + 8])
                barcode_seq1 = seq1[pos1 + 4:pos1 + 12]
                ham_bc1 = checkHamming(bc1s, barcode_seq1)
                ham_bc2 = checkHamming(bc1s, barcode_seq2)
                if ham_bc1[0]:
                    if ham_bc2[0]:
                        print("@%s\n%s\n+\n%s" % (new_title, trim_seq1, trim_qual1), file=fho_r1)
                        print("@%s\n%s\n+\n%s" % (new_title, seq2, qual2), file=fho_r2)
                        bc1 = ham_bc1[1]
                        bc2 = ham_bc2[1]
                        print("%s\t%s\t%s\t%s" % (new_title, bc1, bc2, terminal_bc),file=outfile)
        else:
            if m2:
                pos2 = m2.start()
                pos2_end = m2.end()
                trim_seq2 = seq2[:pos2]
                trim_qual2 = qual2[:pos2]
                #Check for length of the read
#				print seq2[pos2_end:]
                if len(trim_seq2) < 30:
                    continue #and continue if either one is too short
                if len(seq2[pos2_end:]) < 12:
                    continue #and continue if either one is too short
                barcode_seq1 = rc(seq2[pos2_end:pos2_end + 8])
                barcode_seq2 = seq2[pos2 + 4:pos2 + 12]
                ham_bc1 = checkHamming(bc1s, barcode_seq1)
                ham_bc2 = checkHamming(bc1s, barcode_seq2)
                if ham_bc1[0]:
                    if ham_bc2[0]:
                        print("@%s\n%s\n+\n%s" % (new_title, seq1, qual1), file=fho_r1)
                        print("@%s\n%s\n+\n%s" % (new_title, trim_seq2, trim_qual2), file=fho_r2)
                        bc1 = ham_bc1[1]
                        bc2 = ham_bc2[1]
                        print("%s\t%s\t%s\t%s" % (new_title, bc1, bc2, terminal_bc),file=outfile)

def analyze_scDHC_V2design(fhi_bc, fhi_r1,fhi_r2,fho_r1,fho_r2):
	fhi_bc = open(fhi_bc )
	# fhi_r1 = gzip.open(fhi_r1,'r')
	# fhi_r2 = gzip.open(fhi_r2,'r' )
	fho_r1 = open(fho_r1 , 'w')
	fho_r2 = open(fho_r2, 'w')
	assoc = associateBarcodes(fhi_r1, fhi_r2, fho_r1, fho_r2, fhi_bc)
	fhi_bc.close()
	# fhi_r1.close()
	# fhi_r2.close()
	fho_r1.close()
	fho_r2.close()
     
def checkHamming(barcodes,barcode):
	'''Given a list of barcodes, check that the given barcode is within edit\
		distance 2 to any of the list of barcodes'''
	for bc in barcodes:
		match = False
		hd = hamming(barcode, bc)
		if hd <= 2:
			match = True
			barcode = bc
			break
	return (match, barcode)



def inline_splitter(r1, r2, barcode_fhi, r1_o, r2_o, log):

	barcode_fhi = open(barcode_fhi) #Barcode filehandle
	#Barcodes
	barcodes = []
	for line in barcode_fhi:
		barcodes.append(line.split()[0])
	r1_o = open(r1_o, 'w') #R1 filehandle out
	r2_o = open(r2_o, 'w') #R2 filehandle out
	mismatch, not_found_R1, not_found_R2, reads = split_fastqs(r1, r2, r1_o, r2_o, barcodes) #Split some fastqs, yo
	#Close the filehandles
	# r1.close()
	# r2.close()
	r1_o.close()
	r2_o.close()
	barcode_fhi.close()

	# 将mismatch, not_found_R1, not_found_R2, reads写入log
	with open(log, 'w') as f:
		f.write('mismatch: ' + str(mismatch) + '\n')
		f.write('not_found_R1: ' + str(not_found_R1) + '\n')
		f.write('not_found_R2: ' + str(not_found_R2) + '\n')
		f.write('reads: ' + str(reads) + '\n')


	f.close()
     
def cload_correct_zoomify(dir_, pair, genomesize,res,resz):
    target = dir_ + pair
    name = pair.split('.')[0]
    cmd = 'cooler cload pairs -c1 2 -p1 3 -c2 4 -p2 5 '
    cmd += genomesize + ':' + str(res) + ' ' + target + ' ' + dir_  +  name + '_' + str(res) + '.cool'
    # print(cmd)
    time1 = run_cmd_return_time(cmd)

    # cmd = 'hicCorrectMatrix correct --matrix  ' + dir_ + name + '_' + str(res) + '.cool'
    # cmd += ' --correctionMethod KR --outFileName '
    # cmd += ' ' + dir_ + name + '_' + str(res) + '.KR.cool'
    # cmd += ' --filterThreshold -1.5 5.0 --chromosomes chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19'
    # # print(cmd)
    # time1 = run_cmd_return_time(cmd)


    cmd = 'cooler zoomify '  + \
         dir_ + name + '_' +str(res) + \
        '.cool -r ' + resz + ' -o ' + dir_ + name + '.mcool'
    time1 = run_cmd_return_time(cmd) 



def snHic_barcode(file1_path, file2_path):
    #####记录barcode对应的reads
    fp = open('barcode.index', "w")
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
    DirtyFile.close()
    
    
    
# droplet





def hictools_combine(fastq_dir, filename, combine_type='atac'):
    """使用hictools合并Hi-C数据"""
    # move _i.fast.gz 转化为 _Ri.fq.gz
 
    cmd = f'hictools combine_hic {combine_type} {fastq_dir}/{filename}'
    tmp = run_cmd_return_time(cmd)
    # 还原
  
    return tmp
def bowtie_align_barcode(fastq_dir, filename, ref_10x, read_type='R1'):
    """使用Bowtie比对barcode"""
    input_file = f'{fastq_dir}/{filename}_{read_type}_combined.fq.gz'
    output_file = f'{fastq_dir}/{filename}_{read_type}_BC.sam'
    cmd = f'zcat {input_file} | bowtie -p 10 -x {ref_10x} - --nofw -m 1 -v 1 -S {output_file}'
    return run_cmd_return_time(cmd)

def modify_sam_format(fastq_dir, filename, read_type='R1'):
    """修改SAM文件格式 - 使用sc3dg.utils模块"""
    from sc3dg.utils.modify_sam import modify_sam_file as modify_sam_main
    
    input_sam = f'{fastq_dir}/{filename}_{read_type}_BC.sam'
    output_sam = f'{fastq_dir}/{filename}_{read_type}_BC_modified.sam'
    
    t = time.time()
    modify_sam_main(input_sam, output_sam)
    return time.time() - t

def merge_and_filter_reads(fastq_dir, filename):
    """合并、排序、筛选有效配对reads"""
    # 合并
    cmd1 = f'samtools merge -@ 20 {fastq_dir}/{filename}_merged.bam {fastq_dir}/{filename}_R1_BC_modified.sam {fastq_dir}/{filename}_R3_BC_modified.sam'
    time1 = run_cmd_return_time(cmd1)
    
    # 排序
    cmd2 = f'samtools sort -@ 32 -n {fastq_dir}/{filename}_merged.bam -o {fastq_dir}/{filename}_sorted_merged_name.bam'
    time2 = run_cmd_return_time(cmd2)
    
    # 筛选配对reads
    cmd3 = f'samtools view -@ 32 -b -F 8 {fastq_dir}/{filename}_sorted_merged_name.bam > {fastq_dir}/{filename}_paired_only.bam'
    time3 = run_cmd_return_time(cmd3)
    
    # 转换为SAM
    cmd4 = f'samtools view -@ 32 -h {fastq_dir}/{filename}_paired_only.bam > {fastq_dir}/{filename}_paired_only.sam'
    time4 = run_cmd_return_time(cmd4)
    
    return time1 + time2 + time3 + time4

def reconstruct_fastq_with_barcode(fastq_dir, filename):
    """重建包含barcode信息的FASTQ文件 - 使用sc3dg.utils模块"""
    from sc3dg.utils.recon_fq import sam_to_paired_fastq as recon_fq_main
    
    input_sam = f'{fastq_dir}/{filename}_paired_only.sam'
    output_r1 = f'{fastq_dir}/{filename}_R1_BC_cov.fq'
    output_r3 = f'{fastq_dir}/{filename}_R3_BC_cov.fq'
    
    t = time.time()
    recon_fq_main(input_sam, output_r1, output_r3)
    return time.time() - t

def cleanup_intermediate_files(fastq_dir, filename):
    """清理中间文件"""
    files_to_remove = [
        f'{fastq_dir}/{filename}_R1_BC.sam',
        f'{fastq_dir}/{filename}_R3_BC.sam', 
        f'{fastq_dir}/{filename}_R1_BC_modified.sam',
        f'{fastq_dir}/{filename}_R3_BC_modified.sam',
        f'{fastq_dir}/{filename}_sorted_merged_name.bam',
        f'{fastq_dir}/{filename}_paired_only.bam',
        f'{fastq_dir}/{filename}_paired_only.sam',
        f'{fastq_dir}/{filename}_R1_BC_cov.fq',
        f'{fastq_dir}/{filename}_R3_BC_cov.fq'
    ]
    
    for file_path in files_to_remove:
        if os.path.exists(file_path):
            os.remove(file_path)

def process_droplet_hic(opt, fastq, log_out):
    """处理droplet Hi-C数据的主函数"""
    
    # 设置路径和参数
    fastq_dir = opt['fastq_dir'][opt['fastq_log'].index(fastq)]
    fastq_dir = os.path.dirname(fastq_dir)
    code_dir = opt.get('code_dir', '/path/to/code')  # 需要在配置中设置
    ref_10x = opt.get('ref_10x', '/path/to/ref')    # 需要在配置中设置
    
    filename = fastq
    threads = opt['thread']
    
    log_out.info(f'开始处理droplet Hi-C数据: {filename}')
    
    try:
        # 1. 使用hictools合并数据
        t = hictools_combine(fastq_dir, filename, 'atac')
        log_out.info(f'hictools_combine: {t}')
        
        # 2. 使用Bowtie比对barcode (并行处理R1和R3)
        import threading
        
        def align_r1():
            return bowtie_align_barcode(fastq_dir, filename, ref_10x, 'R1')
        
        def align_r3():
            return bowtie_align_barcode(fastq_dir, filename, ref_10x, 'R3')
        
        thread1 = threading.Thread(target=align_r1)
        thread2 = threading.Thread(target=align_r3)
        
        thread1.start()
        thread2.start()
        thread1.join()
        thread2.join()
        
        log_out.info('bowtie比对完成')
       
        # 3. 修改SAM文件格式 (并行处理)
        def modify_r1():
            return modify_sam_format(code_dir, fastq_dir, filename, 'R1')
        
        def modify_r3():
            return modify_sam_format(code_dir, fastq_dir, filename, 'R3')
        
        thread1 = threading.Thread(target=modify_r1)
        thread2 = threading.Thread(target=modify_r3)
        
        thread1.start()
        thread2.start()
        thread1.join()
        thread2.join()
        
        log_out.info('SAM格式修改完成')
        
        # 4. 合并、排序、筛选reads
        t = merge_and_filter_reads(fastq_dir, filename)
        log_out.info(f'merge_and_filter_reads: {t}')
        
        # 5. 重建FASTQ文件
        t = reconstruct_fastq_with_barcode(code_dir, fastq_dir, filename)
        log_out.info(f'reconstruct_fastq_with_barcode: {t}')
        
        # 6. 创建临时目录
        trim_dir = f'{fastq_dir}/tmp'
        if not os.path.exists(trim_dir):
            os.makedirs(trim_dir)
        
        # 7. 使用fastp进行质控 (修改输入文件路径)
        input_r1 = f'{fastq_dir}/{filename}_R1_BC_cov.fq'
        input_r3 = f'{fastq_dir}/{filename}_R3_BC_cov.fq'
        output_r1 = f'{trim_dir}/{filename}_R1_BC_cov_val_1.fq.gz'
        output_r3 = f'{trim_dir}/{filename}_R3_BC_cov_val_2.fq.gz'
        
        t = fastp(input_r1, input_r3, output_r1, output_r3, threads, max=False)
        log_out.info(f'fastp: {t}')
        
        # 8. 使用BWA比对到参考基因组
        t = bwa_mem_specific(opt['index'], threads, filename, trim_dir, output_r1, output_r3)
        log_out.info(f'bwa_mem: {t}')
        
        # 后续步骤与原流程相同，从pairtools开始...
        # 这里可以复用原有的处理逻辑
        
        # 9. 清理中间文件
        cleanup_intermediate_files(fastq_dir, filename)
        
        log_out.info(f'droplet Hi-C处理完成: {filename}')
        
    except Exception as e:
        log_out.error(f'处理droplet Hi-C数据时出错: {str(e)}')
        raise

def bwa_mem_specific(bwaindex, threads, filename, trim_dir, input_r1, input_r3):
    """BWA mem比对的特定实现"""
    output_bam = f'{trim_dir}/{filename}.bam'
    cmd = f'bwa mem -SP5M {bwaindex} -t {threads} {input_r1} {input_r3} | samtools view -bhS -@ 5 - > {output_bam}'
    return run_cmd_return_time(cmd)

def process_droplet_barcode_steps(opt, fastq, fastq_dir, threads):
    """处理droplet特有的barcode处理步骤"""
    start_time = time.time()
    
    # 2. Bowtie比对barcode (并行)
    import threading
    
    def align_r1():
        return bowtie_align_barcode(fastq_dir, fastq, opt['ref_10x'], 'R1')
    
    def align_r3():
        return bowtie_align_barcode(fastq_dir, fastq, opt['ref_10x'], 'R3')
    
    thread1 = threading.Thread(target=align_r1)
    thread2 = threading.Thread(target=align_r3)
    thread1.start()
    thread2.start()
    thread1.join()
    thread2.join()
    
    # 3. 修改SAM格式 (并行)
    def modify_r1():
        return modify_sam_format(fastq_dir, fastq, 'R1')
    
    def modify_r3():
        return modify_sam_format(fastq_dir, fastq, 'R3')
    
    thread1 = threading.Thread(target=modify_r1)
    thread2 = threading.Thread(target=modify_r3)
    thread1.start()
    thread2.start()
    thread1.join()
    thread2.join()
    
    # 4. 合并、排序、筛选
    merge_and_filter_reads(fastq_dir, fastq)
    
    # 5. 重建FASTQ
    reconstruct_fastq_with_barcode(fastq_dir, fastq)
    
    return time.time() - start_time

def cleanup_droplet_intermediate_files(fastq_dir, filename):
    """清理droplet中间文件"""
    files_to_remove = [
        f'{fastq_dir}/{filename}_R1_BC.sam',
        f'{fastq_dir}/{filename}_R3_BC.sam', 
        f'{fastq_dir}/{filename}_R1_BC_modified.sam',
        f'{fastq_dir}/{filename}_R3_BC_modified.sam',
        f'{fastq_dir}/{filename}_sorted_merged_name.bam',
        f'{fastq_dir}/{filename}_paired_only.bam',
        f'{fastq_dir}/{filename}_paired_only.sam',
        f'{fastq_dir}/{filename}_R1_BC_cov.fq',
        f'{fastq_dir}/{filename}_R3_BC_cov.fq'
    ]
    
    for file_path in files_to_remove:
        if os.path.exists(file_path):
            os.remove(file_path)

def bwa_mem_droplet(bwaindex, threads, fastq):
    """droplet专用的BWA比对"""
    cmd = f'bwa mem -SP5M {bwaindex} -t {threads} trimmed-pair1.fastq.gz trimmed-pair2.fastq.gz > {fastq}.sam'
    return run_cmd_return_time(cmd)