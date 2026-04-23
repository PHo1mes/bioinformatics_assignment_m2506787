#!/bin/bash
#Before running this pipeline you should have done the following:

#Setup the directory structure
#mkdir data results
#mkdir data/untrimmed_fastq data/trimmedfastq data/aligned_data data/reference data/other
#mkdir results/fastqc_trimmed_reads results/fastqc_untrimmed_reads results/annovar results/snpEff

#Installed all tools needed

#Downloaded the fastq and annotation files
#wget https://s3-eu-west-1.amazonaws.com/workshopdata2017/NGS0001.R1.fastq.qz
#wget https://s3-eu-west-1.amazonaws.com/workshopdata2017/NGS0001.R2.fastq.qz
#mv *fastq.qz ~/assignment_ngs/data/untrimmed_fastq/
#wget https://s3-eu-west-1.amazonaws.com/workshopdata2017/annotation.bed
# mv annotation.bed ~/assignment_ngs/data/

#Downloaded and indexed your reference genome
#wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
#mv hg19.fa.gz ~/assignment_ngs/data/reference
#bwa index ~/assignment_ngs/data/reference/hg19.fa.gz

#Downloaded the dnSNP database for annovar
#cd ~/annovar
#./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar avsnp147 humandb/

#If you have not done any of these steps, please see the instructions for each step. These will not run with the pipeline.

# Pre-Alignment QC

zcat $1 > ~/assignment_ngs/data/untrimmed_fastq/NGS0001.R1.fastq
zcat $2 > ~/assignment_ngs/data/untrimmed_fastq/NGS0001.R2.fastq
rm ~/assignment_ngs/data/untrimmed_fastq/*.fastq.qz

# FastQC (untrimmed)

fastqc -t 4 ~/assignment_ngs/data/untrimmed_fastq/*.fastq
mv ~/assignment_ngs/data/untrimmed_fastq/*fastqc* ~/assignment_ngs/results/fastqc_untrimmed_reads/

#Trimmomattic

trimmomatic PE  \
-threads 4 \
-phred33 \
~/assignment_ngs/data/untrimmed_fastq/NGS0001.R1.fastq ~/assignment_ngs/data/untrimmed_fastq/NGS0001.R2.fastq \
-baseout ~/assignment_ngs/data/trimmed_fastq/NGS0001_trimmed_R \
ILLUMINACLIP:"/home/ubuntu/anaconda3/pkgs/trimmomatic-0.40-hdfd78af_0/share/trimmomatic-0.40-0/adapters/NexteraPE-PE.fa":2:30:10 \
TRAILING:25 MINLEN:50

rm ~/assignment_ngs/data/untrimmed_fastq/*.fastq
rm ~/assignment_ngs/data/trimmed_fastq/*U

#FastQC (trimmed)

fastqc -t 4 ~/assignment_ngs/data/trimmed_fastq/*P
mv ~/assignment_ngs/data/trimmed_fastq/*fastqc* ~/assignment_ngs/results/fastqc_trimmed_reads/

#Alignment

#BWAMem alignment

bwa mem -t 4 -v 1 -I 250,50 ~/assignment_ngs/data/reference/hg19.fa.gz ~/assignment_ngs/data/trimmed_fastq/NGS0001_trimmed_R_1P ~/assignment_ngs/data/trimmed_fastq/NGS0001_trimmed_R_2P > ~/assignment_ngs/data/aligned_data/NGS0001.sam
rm ~/assignment_ngs/data/trimmed_fastq/*
samtools view -h -b ~/assignment_ngs/data/aligned_data/NGS0001.sam > ~/assignment_ngs/data/aligned_data/NGS0001.bam
samtools sort ~/assignment_ngs/data/aligned_data/NGS0001.bam > ~/assignment_ngs/data/aligned_data/NGS0001_sorted.bam
rm ~/assignment_ngs/data/aligned_data/NGS0001.bam

#Duplicate marking usign Picard MarkDuplicates

picard MarkDuplicates I=~/assignment_ngs/data/aligned_data/NGS0001_sorted.bam O=~/assignment_ngs/data/aligned_data/NGS0001_sorted_marked.bam M=marked_dup_metrics.txt
rm ~/assignment_ngs/data/aligned_data/NGS0001_sorted.bam

#Quality filtering using Samtools

samtools view -F 1796  -q 20 -o ~/assignment_ngs/data/aligned_data/NGS0001_sorted_filtered.bam ~/assignment_ngs/data/aligned_data/NGS0001_sorted_marked.bam
rm ~/assignment_ngs/data/aligned_data/NGS0001_sorted_marked.bam
samtools index ~/assignment_ngs/data/aligned_data/NGS0001_sorted_filtered.bam

#Generating alignment statistics

#Flagstats
samtools flagstat ~/assignment_ngs/data/aligned_data/NGS0001_sorted_filtered.bam > ~/assignment_ngs/data/other/NGS0001_flagstat

#idxstats

samtools idxstats ~/assignment_ngs/data/aligned_data/NGS0001_sorted_filtered.bam > ~/assignment_ngs/data/other/NGS0001_idxstats

#Insertsizemetrics

picard CollectInsertSizeMetrics I=~/assignment_ngs/data/aligned_data/NGS0001_sorted_filtered.bam O=~/assignment_ngs/data/other/NGS0001_insert_size_metrics.txt H=~/assignment_ngs/data/other/NGS0001_insert_size_histogram.pdf M=0.5 > ~/assignment_ngs/data/other/NGS0001_CollectInsertSizeMetrics

#Depth of coverage

bedtools coverage -a ~/assignment_ngs/data/annotation.bed -b ~/assignment_ngs/data/aligned_data/NGS0001_sorted_filtered.bam > ~/assignment_ngs/data/other/NGS0001_coverage

#Variant Calling

#Variant calling using Freebayes

zcat ~/assignment_ngs/data/reference/hg19.fa.gz > ~/assignment_ngs/data/reference/hg19.fa
samtools faidx ~/assignment_ngs/data/reference/hg19.fa
freebayes --bam ~/assignment_ngs/data/aligned_data/NGS0001_sorted_filtered.bam --fasta-reference ~/assignment_ngs/data/reference/hg19.fa --vcf ~/assignment_ngs/results/NGS0001.vcf
rm ~/assignment_ngs/data/aligned_data/NGS0001_sorted_filtered.bam
bgzip ~/assignment_ngs/results/NGS0001.vcf

#Quality filtering

tabix -p vcf ~/assignment_ngs/results/NGS0001.vcf.gz
vcffilter -f "QUAL > 1 & QUAL / AO > 10 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1" ~/assignment_ngs/results/NGS0001.vcf.gz > ~/assignment_ngs/results/NGS0001filtered.vcf
bedtools intersect -header -wa -a ~/assignment_ngs/results/NGS0001filtered.vcf -b ~/assignment_ngs/data/annotation.bed > ~/assignment_ngs/results/NGS0001_filtered_intersect.vcf
bgzip ~/assignment_ngs/results/NGS0001_filtered_intersect.vcf
tabix -p vcf ~/assignment_ngs/results/NGS0001_filtered_intersect.vcf.gz

#Variant Annotation and Prioritisation

#Annovar

#Annotation

cd ~/annovar
./convert2annovar.pl -format vcf4 ~/assignment_ngs/results/NGS0001_filtered_intersect.vcf.gz > ~/assignment_ngs/results/NGS0001_filtered_intersect.avinput
./table_annovar.pl ~/assignment_ngs/results/NGS0001_filtered_intersect.avinput humandb/ -buildver hg19 -out ~/assignment_ngs/results/NGS0001_filtered_intersect -remove -protocol refGene,ensGene,clinvar_20180603,exac03,dbnsfp31a_interpro,avsnp147 -operation g,g,f,f,f,f -otherinfo -nastring . -csvout
cd ~

#Filtering for exonic variants not in dbSNP

awk -F',' '$6 == "\"exonic\"" && $30 == "."' ~/assignment_ngs/results/NGS0001_filtered_intersect.hg19_multianno.csv > ~/assignment_ngs/results/NGS0001_annovar_exonic_notindbSNP.csv
cd ~

#snpEff

java -Xmx16g -jar /home/ubuntu/anaconda3/pkgs/snpeff-4.3.1t-0/share/snpeff-4.3.1t-0/snpEff.jar ann hg19 ~/assignment_ngs/results/NGS0001_filtered_intersect.vcf.gz > ~/assignment_ngs/results/NGS0001_snpEff_anno.vcf

#snpSift

java -Xmx4g -jar /home/ubuntu/anaconda3/pkgs/snpsift-4.3.1t-hdfd78af_3/share/snpsift-4.3.1t-3/SnpSift.jar filter "(ANN[*].EFFECT =~ 'exon') & (ID = '.')" ~/assignment_ngs/results/NGS0001_snpEff_anno.vcf > ~/assignment_ngs/results/NGS0001_snpEff_filtered.vcf

