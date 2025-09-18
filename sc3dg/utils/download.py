import os
import sys

genome_url = {
    'mm10':'https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz',
    'mm9':'https://hgdownload.soe.ucsc.edu/goldenPath/mm9/bigZips/mm9.fa.gz',
    'mm39':'https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.fa.gz',
    'hg38':'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz',
    'hg19':'https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz'
}
chromosome_size_url = {
    'mm10':'https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes',
    'mm9':'https://hgdownload.soe.ucsc.edu/goldenPath/mm9/bigZips/mm9.chrom.sizes',
    'mm39':'https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.chrom.sizes',
    'hg38':'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes',
    'hg19':'https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes'
}


def make_index(path, version='hg38', aligner='bwa'):
    if not os.path.exists(path):
        os.makedirs(path)

    os.chdir(path)

    # 下载fa
    cmd = 'wget %s' % genome_url[version]
    os.system(cmd)

    cmd = 'wget %s' % chromosome_size_url[version]
    os.system(cmd)

    if aligner == 'bwa':
        cmd = 'bwa index %s.fa.gz' % version
        os.system(cmd)
    elif aligner == 'bowtie2':
        cmd = 'bowite2 index %s.fa.gz %s' % (version, version)
        os.system(cmd)
    else:
        sys.exit('aligner should be bwa/bowtie2')

    
def make_digest(path, enzyme='MboI',genome='hg38'):
    if not os.path.exists(path):
        sys.exit('path is not exist')
    print(os.path.exists(path + '/%s.fa' % genome))
    if not os.path.exists(path + '/%s.fa.gz' % genome) and not os.path.exists(path + '/%s.fa' % genome):
        sys.exit('you should run index first')
    if not os.path.exists(path + '/%s.fa.chrom.sizes' % genome):
        sys.exit('you should run index first')

    check_enzyme(enzyme)
    
    os.chdir(path)
    size = './%s.fa.chrom.sizes' % genome
    fa = './%s.fa.gz' % genome if os.path.exists('./%s.fa.gz' % genome) else './%s.fa' % genome
    
    cmd = 'cooler digest  -o digest_%s.bed %s %s %s' % (enzyme, size, fa, enzyme)
    print('digesting.....')
    os.system(cmd)
    

def check_enzyme(enzyme):
     import Bio.Restriction as biorst

     if enzyme not in biorst.AllEnzymes:
        lists = '['
        for val in list(biorst.AllEnzymes):
       
            lists += str(val )
            lists += ', '
        lists += ']'
        
        sys.exit('enzyme is not exist, you can use: ' + lists)




