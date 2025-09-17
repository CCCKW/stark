import click
import os
import subprocess
from sc3dg.utils import scSPRITE,sn_m3c,scNano,scMethyl
from sc3dg.utils import help
from sc3dg.utils import tools as tl
import time
from multiprocessing import Process, Queue, Pool, Manager
from sc3dg.utils.scaffold import *
from sc3dg.utils.assembly import *
import logging

Path = os.getcwd()
print(Path)

tech = {
        'scSPRITE':scSPRITE,
        'sn_m3c':sn_m3c,
        'scNano':scNano,
        'scMethyl': scMethyl
       }

def run_exec(type_, opt, fq, log_out):
    if type_ in ['sn_m3c', 'scSPRITE', 'scNano', 'scMethyl']:
        if type_ == 'sn_m3c' and opt['aligner'] == 'bwa':
            log_out.debug('sn_m3c only support bowtie2')
            print('sn_m3c only support bowtie2')
            return
        exec(f"{type_}.pp(opt=opt, fastq=fq, log_out=log_out)")
    else:
        assembly(type_, opt, fq, log_out)

@click.command('count', short_help='the main pipeline of Hi-C')
@click.option('-o', '--output', help='Save path', required=True)
@click.option('-f', '--fastq', required=True, help='Fastq directory, run all if -s is not specified')
@click.option('--logging', help='Logging file path', default='./log.log')
@click.option('-t', '--type', required=True, type=click.Choice(['scHic', 'snHic', 'sciHic', 'scSPRITE', 'dipC', 'sn_m3c', 'HIRES', 'scNano', 'scMethyl', 'LiMAC', 'paired','droplet','GAGE-seq']), help='Type of Hi-C')
@click.option('-e', '--enzyme', help='Enzyme, e.g., mboi', required=True, type=str)
@click.option('-r', '--resolution', help='Resolution', default=10000, type=int)
@click.option('-i', '--index', help='BWA fa file/Bowtie fa file', required=True)
@click.option('-s', '--sample', help='Sample file, txt', default=None)
@click.option('--exist-barcode', is_flag=True, help='Barcode exists')
@click.option('--qc', type=int, default=0, help='Samtools view to QC')
@click.option('--add-columns', default='mapq', help=help.help['add_columns'])
@click.option('--thread', help='Thread count', type=int, default=20)
@click.option('--adaptor-file-bc2', help='Adaptor file for barcode 2, default is None', default=None)
@click.option('--adaptor-file-l2', help='Adaptor file for barcode 2, default is None', default=None)
@click.option('--adaptor-file-bc1', help='Adaptor file for barcode 2, default is None', default=None)
@click.option('--worker', help='Simultaneously run the worker of the pipeline', type=int, default=2)
@click.option('--select', help=help.help['select'], default="mapq1>=30 and mapq2>=30")
@click.option('--max-mismatch', help=help.help['max_mismatch'], type=int, default=3)
@click.option('--aligner', default='bwa', type=click.Choice(['bwa', 'bowtie2']), help='Aligner')
@click.option('--repeat-masker', help='Repeat masker (bed file) for scSPRITE', default=None)
@click.option('--sprite-config', default=None ,help='scSPRITE config')
@click.option('--scnano-barcode', help='scNano barcode for PCR and TN5, stored in a folder and named as TN5.txt and PCR index.txt respectively', default=None)
@click.option('--zoomify-res', help='zoomify',type=str, default='10000,40000,100000,500000,1000000')
@click.option('--ref-10x', help='10x reference index for bowtie alignment (required for droplet)', default=None)
@click.option('--droplet-10x', help='10x Chromium barcode list file', default=None)
@click.option('--paired-10x', help='10x ARC fasta file for bowtie index', default=None)
def count(output, fastq, logging, type_, enzyme, resolution, index, sample, exist_barcode, qc, 
         add_columns, thread, worker, select, max_mismatch, aligner, repeat_masker, sprite_config, 
         scnano_barcode, zoomify_res, ref_10x, droplet_10x, paired_10x, 
         adaptor_file_bc2, adaptor_file_l2, adaptor_file_bc1):
    """Run the Hi-C pipeline with the specified options."""
    opt = {
        'output': output,
        'fastq': fastq,
        'logging': logging,
        'type': type_,
        'enzyme': enzyme,
        'resolution': resolution,
        'index': index,
        'sample': sample,
        'exist_barcode': exist_barcode,
        'qc': qc,
        'add_columns': add_columns,
        'thread': thread,
        'worker': worker,
        'select': select,
        'max_mismatch': max_mismatch,
        'aligner': aligner,
        'repeat_masker': repeat_masker,
        'sprite_config': sprite_config,
        'scNano_barcode': scnano_barcode,
        'zoomify_res': zoomify_res.split(','),
        'ref_10x': ref_10x,
        'droplet_10x': droplet_10x,
        'paired_10x': paired_10x,
        'adaptor-file-bc2': adaptor_file_bc2,
        'adaptor-file-l2': adaptor_file_l2,
        'adaptor-file-bc1': adaptor_file_bc1,
    }
    opt['running_path'] = Path

    # 先创建output文件夹，没有则创建，有则退出;然后移到该工作目录下面
    tl.make_workflow(opt['output'])
    os.chdir(opt['output'])

    # 然后判断是不是到底有多少个样本，如果是一个样本，那么就直接运行，如果是多个样本，那么就要用多进程来运行
    opt, fastq_log = tl.count_sample(opt)

    # 判断index是否存在
    tl.check_index(opt)
    tl.check_enzyme(opt)

    tl.get_enzyme_bed(opt)

    opt = tl.parser_fa_and_chrom_size(opt)

    # 查询是什么物种
    opt['species'] = tl.parse_species(opt)
    
    # droplet类型的特殊检查
    if opt['type'] == 'droplet' or opt['type'] == 'GAGE-seq':
        # check bowtie
        def check_bowtie():
            try:
                subprocess.run(['bowtie', '--version'], 
                            capture_output=True, 
                            check=True)
                return True
            except (subprocess.CalledProcessError, FileNotFoundError):
                return False
        
        def check_hictools():
            try:
                subprocess.run(['hictools', '--help'], 
                            capture_output=True, 
                            check=True)
                return True
            except (subprocess.CalledProcessError, FileNotFoundError):
                return False
        
        if not check_bowtie():
            print("Error: bowtie is not installed or not found in PATH.")
            return
            
        if not check_hictools():
            print("Error: hictools is not installed or not found in PATH.")
            return
            
        # 检查droplet必需的配置
        if opt['type'] == 'droplet':
            required_configs = ['ref_10x']
            for config in required_configs:
                if not opt[config]:
                    print(f"Error: --{config.replace('_', '-')} is required for droplet Hi-C processing")
                    return

    # global logging
    log_out = tl.log_(opt)
    os.chdir(opt['output'])

    if len(fastq_log) <= opt['worker']:
        print('*****num of worker is more than num of sample*****')
        p = Pool(opt['worker'])
        for fq in fastq_log:
           run_exec(opt['type'], opt, fq, log_out)
        p.close()
        p.join()
    else:
        print('*****num of worker is less than num of sample*****')
        divide_list = []
        for i in range(0, len(fastq_log), opt['worker']):
            divide_list.append(fastq_log[i:i + opt['worker']])

        # logging batch 
        log_out.debug('divide_list:{}'.format(len(divide_list)))

        for i, val in enumerate(divide_list):
            p = Pool(len(val))
            for j, fq in enumerate(val):
                log_out.debug('%s || %s.....Dealing fq: %s \n' % (i*opt['worker']+ j, len(fastq_log), fq))
                p.apply_async(run_exec, args=(opt['type'], opt, fq, log_out))
                
            p.close()
            p.join()
            log_out.debug('Batch %s done || %s \n' % (i, len(divide_list)/opt['worker']))
    
    log_out.debug('All done || %s \n' % (len(fastq_log)))