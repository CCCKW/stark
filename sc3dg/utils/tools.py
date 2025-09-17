import os
import sys
import logging

genome = ['mm10', 'hg38', 'hg19','mm9']
bwa_prefix = ['amb', 'ann']
def check_enzyme(opt):
     import Bio.Restriction as biorst

     if opt['enzyme'] not in biorst.AllEnzymes:
        lists = '['
        for val in list(biorst.AllEnzymes):
       
            lists += str(val )
            lists += ', '
        lists += ']'
        
        sys.exit('enzyme is not exist, you can use: ' + lists)

def get_enzyme_bed(opt):
    enzyme = opt['enzyme']
    index = opt['index']
    index = index.rsplit('/', 1)[0]
    if not os.path.exists(index + '/digest_%s.bed' % enzyme):
        sys.exit('digest_%s.bed is not exist' % enzyme)
    else:
        opt['enzyme_bed'] = index + '/digest_%s.bed' % enzyme

def check_index(opt):
    print(opt['index'])
    if not opt['type'] == 'sn_m3c':
        if opt['aligner']  == 'bowtie2' :
            if not os.path.exists(opt['index'] + '.1.bt2'):
                sys.exit('index is not exist/index is not bowtie2 index')
            
        else:
            if not os.path.exists(opt['index'] + '.amb'):
                sys.exit('index is not exist')
    else:
        if not os.path.exists(opt['index'].rsplit('/', 1)[0] + '/Bisulfite_Genome'):
            print('bismark index is not exist')
  
            if not os.path.exists(opt['index'] + '.1.bt2'):
                sys.exit(opt['index'] + '  is bowtie2 index. fail to create bismark index')
            else:
                print('we will build bismark index')
        else:
            print('bismark index is exist')

def parser_fa_and_chrom_size(opt):
    index = opt['index']
    index = index.rsplit('/', 1)[0]
    flag = 0
    if opt['aligner'] == 'bowtie2':
        for val in os.listdir(index):
            if val.endswith('.fa'):
                opt['fa_index'] = index + '/' + val
                flag += 1
                break
        for val in os.listdir(index):
            if val.endswith('.chrom.sizes'):
                opt['genomesize'] = index + '/' + val
                flag += 1
                break
    else:
        for val in os.listdir(index):
            if val.endswith('.fa'):
                opt['fa_index'] = index + '/' + val
                flag += 1
                break
        for val in os.listdir(index):
            if val.endswith('.chrom.sizes'):
                opt['genomesize'] = index + '/' + val
                flag += 1
                break
    if flag == 2:
        return opt
    raise Exception('fa file is not exist')


def make_workflow(path):
    if os.path.exists(path):
        print('****output folder is already exist*********')
    else:
        os.makedirs(path)

def parse_species(opt):
    for val in genome:
        if val in opt['index'].lower():
            return val

def count_sample(opt:dict):

    fastq_log = []
    for root, dir_, files in os.walk(opt['fastq']):
        for file in files:
            if file.endswith('_1.fastq.gz') or file.endswith('_1.fastq'):
                fastq_log.append(file.split('_1.fastq')[0])
   
    
              

    fastq_dir = []
    for root, dir_, files in os.walk(opt['fastq']):
       for file in files:
            if '_1.fastq' in file:
               fastq_dir.append(root+'/'+file)
      
    
    fastq_dir_final = []

    for val in fastq_log:
        flag, file = return_(val, fastq_dir)
        if flag:
            fastq_dir_final.append(root + '/' + file)
        else:
            sys.exit('fastq_dir is not exist')

 
    opt['fastq_dir'] = fastq_dir

    
    opt['fastq_log'] = fastq_log
            
    opt['run_sample'] = fastq_log
 
    if opt['sample'] is not None:
        with open(opt['sample'], 'r') as f:
            lines = f.readlines()
        lines = [x.strip() for x in lines]
   
        tmp_dir = []
        tmp_log = []
        tmp_run = []
        for i,j,k in zip(opt['fastq_dir'],opt['fastq_log'],opt['run_sample']):
            if k in lines:
                # print(i,j,k, lines.index(k))
                tmp_dir.append(i)
                tmp_log.append(j)
                tmp_run.append(k)
        opt['fastq_dir'] = tmp_dir
        opt['fastq_log'] = tmp_log
        opt['run_sample'] = tmp_run
        
    return opt, opt['fastq_log']

def return_(val, files):
    # print(val,files)
    for f in files:
        if val in f:
            return True, f
    return False, None

def log_(opt):
    log_out = logging.getLogger('log_out')
    log_out.setLevel(logging.DEBUG)
    os.chdir(opt['running_path'])
    fh = logging.FileHandler(filename=opt['logging'],mode='w')

    log_out.addHandler(fh)


    log_out.debug('*********run pipleine*********')
    for k in opt.keys():
        log_out.debug('# ' + k + ': ' + str(opt[k]))
    
    log_out.debug('*********Total: %s*********' % len(opt['fastq_log']))
    return log_out
    
    

def check_index_exists(log_out, opt):
    flag = True
    with gzip.open('trimmed-pair1.fastq.gz','rt') as f:
        for i,l in enumerate(f):
            if i == 100:
                break
            
            if '@' == l[0] and ':' not in l:
                flag = False
    if not flag:
        log_out.debug('fastq index is not exist')
        sys.exit('fastq index is not exist')
