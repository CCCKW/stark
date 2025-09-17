#!/bin/bash
# 指定使用 bash 解释器运行此脚本

set +o posix
# 禁用 POSIX 模式（某些 shell 行为可能更宽松），通常用于避免某些兼容性问题

bwaindex=/cluster/home/jwj/reference/mm10/bwa/mm10.fa
# BWA 比对所用的参考基因组索引文件前缀（小鼠 mm10）

genomesize=/cluster/home/jwj/reference/mm10/bwa/mm10.chrom.sizes
# 基因组染色体大小文件，用于 pairtools 等工具识别染色体长度

CodeDir=/cluster/home/jwj/Project/scHIC_pipeline/3_Pipeline_summary/16_DroplitHiC/Code
# 存放自定义脚本和工具的代码目录

fastq_dir=/cluster/home/jwj/Project/scHIC_pipeline/3_Pipeline_summary/16_DroplitHiC/fastq_dir
# 原始或中间 FASTQ 文件的存储路径

trim_dir=/cluster/home/jwj/Project/scHIC_pipeline/3_Pipeline_summary/16_DroplitHiC/fastq_dir/tmp
# 用于存放质控（trimming）后中间文件的临时目录

FileName=SRR27586280
# 当前处理的样本名称（SRA 编号示例）

Droplet_10X=/cluster/home/jwj/Project/scHIC_pipeline/3_Pipeline_summary/16_DroplitHiC/Code/737K-cratac-v1
# 10x Chromium 的 barcode 列表文件（Droplet 版本）

Paired_10X=/cluster/home/jwj/Project/scHIC_pipeline/3_Pipeline_summary/16_DroplitHiC/Code/737K-arc-v1.fa
# 10x ARC（ATAC-HiC）使用的 fasta 文件（可能用于 bowtie 索引）



#### Droplet Hi-C
${CodeDir}/hictools/hictools combine_hic atac ${fastq_dir}/${FileName}
# 使用自定义工具 hictools 将原始 R1/R2/R3 FASTQ 文件按 10x 结构合并
# 'atac' 模式表示处理的是类似 10x ATAC-HiC 的数据结构
# 输出可能是 _R1_combined.fq.gz 和 _R3_combined.fq.gz（含 barcode 和 read）

#### Paired Hi-C
#${CodeDir}/hictools/hictools combine_hic arc ${fastq_dir}/${FileName}
# 注释掉的命令：用于处理另一种类型的 paired-end 10x ARC Hi-C 数据



#使用 Bowtie 提取并比对 Barcode
zcat ${fastq_dir}/${FileName}_R1_combined.fq.gz | bowtie -p 10 -x ${ref_10X} - --nofw -m 1 -v 1 -S ${fastq_dir}/${FileName}_R1_BC.sam &
zcat ${fastq_dir}/${FileName}_R3_combined.fq.gz | bowtie -p 10 -x ${ref_10X} - --nofw -m 1 -v 1 -S ${fastq_dir}/${FileName}_R3_BC.sam &
wait

#修改 SAM 文件格式（适配下游处理）
python ${CodeDir}/modify_sam.py ${fastq_dir}/${FileName}_R1_BC.sam ${fastq_dir}/${FileName}_R1_BC_modified.sam &
python ${CodeDir}/modify_sam.py ${fastq_dir}/${FileName}_R3_BC.sam ${fastq_dir}/${FileName}_R3_BC_modified.sam &
wait

#合并、排序、筛选有效配对 reads
samtools merge -@ 20 ${fastq_dir}/${FileName}_merged.bam ${fastq_dir}/${FileName}_R1_BC_modified.sam  ${fastq_dir}/${FileName}_R3_BC_modified.sam
samtools sort -@ 32 -n ${fastq_dir}/${FileName}_merged.bam -o ${fastq_dir}/${FileName}_sorted_merged_name.bam
samtools view -@ 32 -b -F 8 ${fastq_dir}/${FileName}_sorted_merged_name.bam > ${fastq_dir}/${FileName}_paired_only.bam
samtools view  -@ 32 -h ${fastq_dir}/${FileName}_paired_only.bam > ${fastq_dir}/${FileName}_paired_only.sam

#重建包含 barcode 信息的 FASTQ 文件
python ${CodeDir}/recon_fq.py ${fastq_dir}/${FileName}_paired_only.sam ${fastq_dir}/${FileName}_R1_BC_cov.fq ${fastq_dir}/${FileName}_R3_BC_cov.fq


#创建临时目录并进行质控（fastp）
mkdir ${trim_dir}

fastp \
  -i ${fastq_dir}/${FileName}_R1_BC_cov.fq \
  -I ${fastq_dir}/${FileName}_R3_BC_cov.fq \
  -o ${trim_dir}/${FileName}_R1_BC_cov_val_1.fq.gz \
  -O ${trim_dir}/${FileName}_R3_BC_cov_val_2.fq.gz \
  --html ${trim_dir}/${FileName}_trim_report.html \
  --json ${trim_dir}/${FileName}_trim_report.json \
  --thread 16

#比对到参考基因组（BWA）
bwa mem -SP5M ${bwaindex}  -t 64  ${trim_dir}/${FileName}_R1_BC_cov_val_1.fq.gz ${trim_dir}/${FileName}_R3_BC_cov_val_2.fq.gz | samtools view -bhS -@ 5 - > ${trim_dir}/${FileName}.bam

####Generate pairs
pairtools parse ${trim_dir}/${FileName}.bam -c ${genomesize} \
  --drop-sam --drop-seq --output-stats ${trim_dir}/${FileName}.stats \
  --assembly hg38 --no-flip \
  --add-columns mapq \
  --walks-policy all \
  -o ${trim_dir}/${FileName}.pairs.gz

##split single cell
python ${CodeDir}/splitCell.py ${trim_dir}/${FileName}.pairs.gz

#QC select
pairtools select "mapq1>=30 and mapq2>=30" ${trim_dir}/${FileName}.pairs.gz -o ${trim_dir}/${FileName}.filted.pairs.gz

#Sort pairs
pairtools sort --nproc 32 ${FileName}.filted.pairs.gz -o ${FileName}.sorted.pairs.gz

#dedup
pairtools dedup \
    --max-mismatch 3 \
    --mark-dups \
    --output >( pairtools split --output-pairs ${FileName}.nodups.pairs.gz --output-sam ${FileName}.nodups.bam ) \
    --output-unmapped >( pairtools split --output-pairs ${FileName}.unmapped.pairs.gz --output-sam ${FileName}.unmapped.bam ) \
    --output-dups >( pairtools split --output-pairs ${FileName}.dups.pairs.gz --output-sam ${FileName}.dups.bam) \
    --output-stats ${FileName}.dedup.stats \
    ${FileName}.sorted.pairs.gz


###
cooler cload pairs -c1 2 -p1 3 -c2 4 -p2 5 ${genomesize}:10000 ${FileName}.nodups.pairs.gz ${FileName}.10000.cool

hicCorrectMatrix correct --matrix  ${FileName}.${resoluation}.cool  --correctionMethod KR --outFileName  ${FileName}.${resoluation}.KR.cool --filterThreshold -1.5 5.0 --chromosomes chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX 

cooler zoomify ${FileName}.${resoluation}.KR.cool -r 10000,50000,100000,500000 -o ${FileName}.mcool



rm ${fastq_dir}/${FileName}_R1_BC.sam ${fastq_dir}/${FileName}_R3_BC.sam ${fastq_dir}/${FileName}_R1_BC_modified.sam ${fastq_dir}/${FileName}_R3_BC_modified.sam 
rm ${fastq_dir}/${FileName}_sorted_merged_name.bam ${fastq_dir}/${FileName}_paired_only.bam ${fastq_dir}/${FileName}_paired_only.sam
rm ${fastq_dir}/${FileName}_R1_BC_cov.fq  ${fastq_dir}/${FileName}_R3_BC_cov.fq
rm ${trim_dir}/${FileName}_R1_BC_cov_val_1.fq.gz ${trim_dir}/${FileName}_R3_BC_cov_val_2.fq.gz














