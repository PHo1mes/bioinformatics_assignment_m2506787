#!/bin/bash
##Before running this pipeline you should have done the following:

# Setup the directory structure:
#mkdir data results
#mkdir data/untrimmed_fastq data/trimmedfastq data/aligned_data data/reference data/other
#mkdir results/fastqc_trimmed_reads results/fastqc_untrimmed_reads results/annovar results/snpEff

# Installed all tools needed:
#wget https://repo.anaconda.com/archive/Anaconda3-2022.10-Linux-x86_64.sh
#chmod +x ./Anaconda3-2022.10-Linux-x86_64.sh
#bash ./Anaconda3-2022.10-Linux-x86_64.sh
#source ~/.bashrc
#conda config --add channels defaults
#conda config --add channels bioconda
#conda config --add channels conda-forge
#conda install samtools
#conda install bwa
#conda install picard
#conda install bedtools
#conda install trimmomatic
#conda install fastqc
#conda install -c bioconda snpeff=4.3 #snpEff version 4.3 is used as later versions are not compatible with the version of java installed on the machine
#conda install deepvariant
#sudo apt install libvcflib-tools


# Downloaded the fastq and annotation files, and moved them to the correct locations for this pipeline to run:
#wget https://s3-eu-west-1.amazonaws.com/workshopdata2017/NGS0001.R1.fastq.qz
#wget https://s3-eu-west-1.amazonaws.com/workshopdata2017/NGS0001.R2.fastq.qz
#mv *fastq.qz ~/assignment_ngs/data/untrimmed_fastq/
#wget https://s3-eu-west-1.amazonaws.com/workshopdata2017/annotation.bed
#mv annotation.bed ~/assignment_ngs/data/

# Downloaded and indexed your reference genome:
#wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
#mv hg19.fa.gz ~/assignment_ngs/data/reference
#bwa index ~/assignment_ngs/data/reference/hg19.fa.gz

# Downloaded the dnSNP database for annovar:
#cd ~/annovar
#./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar avsnp147 humandb/

# Downloaded the hg19 genome for snpEff
# snpEff download -v hg19

#If you have not done any of these, please see the instructions for each step. These will not run with the pipeline however are required for it to work.

## 2.2 Pre-Alignment QC

# changes the files out of the fastq.qz format which cannot be read by fastqc or trimmomatic and clears the original files to create space
zcat $1 > ~/assignment_ngs/data/untrimmed_fastq/NGS0001.R1.fastq
zcat $2 > ~/assignment_ngs/data/untrimmed_fastq/NGS0001.R2.fastq
rm ~/assignment_ngs/data/untrimmed_fastq/*.fastq.qz

# FastQC (untrimmed)
# runs fastqc using 4 threads of the cpu and moves the results to the correct directory
fastqc -t 4 ~/assignment_ngs/data/untrimmed_fastq/*.fastq
mv ~/assignment_ngs/data/untrimmed_fastq/*fastqc* ~/assignment_ngs/results/fastqc_untrimmed_reads/

#Trimmomatic
# Trimmomatic is used to remove low quality reads and adaptor sequences from the fastq file
# PE instructs that this is paired end sequencing data, -threads 4 instructs to use 4 CPU threads and -phred33 instructs how the quality score is encoded
# ILLUMINACLIP is used to instruct what the adaptor sequences are.
# TRAILING:25 removes low quality bases from the 3' end by cutting every base up till a base with a quality of 25 is reached.
# MINLEN:50 removes reads that are shorter than 50bp after the trimming is done
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

## 2.3 Alignment

#BWAMem alignment
# Runs BWA mem to align the trimmed sequences to the reference genome
# -t 4 instructs the cpu to use 4 threads for this
# -v 1 sets the verbosity level, meaning how much is written into the command prompt while the command is running
# -R sets the read group information, which has been updated to match the new input files
# -I 250,50 instructs the expected insert size to be 250 bases with a standard deviation of 50
bwa mem -t 4 -v 1 -R '@RG\tID:NGS0001\tSM:NGS0001\tLB:lib1\tPL:ILLUMINA' -I 250,50 ~/assignment_ngs/data/reference/hg19.fa.gz ~/assignment_ngs/data/trimmed_fastq/NGS0001_trimmed_R_1P ~/assignment_ngs/data/trimmed_fastq/NGS0001_trimmed_R_2P > ~/assignment_ngs/data/aligned_data/NGS0001.sam

#converting the file from .sam format to .bam format, which is required for the upcoming steps
samtools view -h -b ~/assignment_ngs/data/aligned_data/NGS0001.sam > ~/assignment_ngs/data/aligned_data/NGS0001.bam
rm ~/assignment_ngs/data/aligned_data/NGS0001.sam
# sorting the .bam file by genomic coordinate, which is required for it to be read in upcoming steps
samtools sort ~/assignment_ngs/data/aligned_data/NGS0001.bam > ~/assignment_ngs/data/aligned_data/NGS0001_sorted.bam
rm ~/assignment_ngs/data/aligned_data/NGS0001.bam

# Duplicate marking usign Picard MarkDuplicates
# performs duplicate marking to identify sequences PCR duplicates
picard MarkDuplicates I=~/assignment_ngs/data/aligned_data/NGS0001_sorted.bam O=~/assignment_ngs/data/aligned_data/NGS0001_sorted_marked.bam M=marked_dup_metrics.txt
rm ~/assignment_ngs/data/aligned_data/NGS0001_sorted.bam

# Quality filtering using Samtools
# removing any low quality reads
# -F 1796 is used as a flag for various reads that should be removed, such as unmapped, secondarily aligned and duplicate reads
# -q 20 sets the minimum quality score to be kept at 20
samtools view -F 1796  -q 20 -o ~/assignment_ngs/data/aligned_data/NGS0001_sorted_filtered.bam ~/assignment_ngs/data/aligned_data/NGS0001_sorted_marked.bam
rm ~/assignment_ngs/data/aligned_data/NGS0001_sorted_marked.bam
samtools index ~/assignment_ngs/data/aligned_data/NGS0001_sorted_filtered.bam
#the bam and index file are kept so that the genome can be viewed, such as by IGV

# Generating alignment statistics

# Flagstats
# generates a summary of what types of reads are found in the file, based on their flags
samtools flagstat ~/assignment_ngs/data/aligned_data/NGS0001_sorted_filtered.bam > ~/assignment_ngs/data/other/NGS0001_flagstat

# idxstats
# generates a summary of which chromosomes the reads have mapped to
samtools idxstats ~/assignment_ngs/data/aligned_data/NGS0001_sorted_filtered.bam > ~/assignment_ngs/data/other/NGS0001_idxstats

# Insertsizemetrics
# generates a summary of the insert sizes that have been kept after filtering
picard CollectInsertSizeMetrics I=~/assignment_ngs/data/aligned_data/NGS0001_sorted_filtered.bam O=~/assignment_ngs/data/other/NGS0001_insert_size_metrics.txt H=~/assignment_ngs/data/other/NGS0001_insert_size_histogram.pdf M=0.5 > ~/assignment_ngs/data/other/NGS0001_CollectInsertSizeMetrics

# Depth of coverage
# generates a summary of the coverage, meaning how much of the genome the reads map to, and the depth meaning how many reads map to a given region
bedtools coverage -a ~/assignment_ngs/data/annotation.bed -b ~/assignment_ngs/data/aligned_data/NGS0001_sorted_filtered.bam > ~/assignment_ngs/data/other/NGS0001_coverage

## 2.4 Variant Calling

# Variant calling using deepvariant
# preparing the reference genome so it can be indexed
zcat ~/assignment_ngs/data/reference/hg19.fa.gz > ~/assignment_ngs/data/reference/hg19.fa
# indexing the reference genome to enable efficient reading by deepvariant
samtools faidx ~/assignment_ngs/data/reference/hg19.fa
# using deepvariant to find the variants between the sequencing data and the reference genome
deepvariant \
  --model_type=WGS \
  --ref=~/assignment_ngs/data/reference/hg19.fa \
  --reads=~/assignment_ngs/data/aligned_data/NGS0001_sorted_filtered.bam \
  --output_vcf=~/assignment_ngs/results/NGS0001.vcf \
  --output_gvcf=~/assignment_ngs/results/NGS0001..g.vcf.gz \
  --num_shards=4
bgzip ~/assignment_ngs/results/NGS0001.vcf

# Quality filtering

tabix -p vcf ~/assignment_ngs/results/NGS0001.vcf.gz
# this filters the vcf file to remove ultra low quality variants, QUAL > 1, variants without suitable supporting reads for their low quality, QUAL / A0 > 10, variants only found the forward or reverse strand, SAF > 0 & SAR > 0, and variants where the reads are not supported on both sides, RPR > 1 & RPL > 1.
vcffilter -f "QUAL > 1 & QUAL / AO > 10 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1" ~/assignment_ngs/results/NGS0001.vcf.gz > ~/assignment_ngs/results/NGS0001filtered.vcf
bedtools intersect -header -wa -a ~/assignment_ngs/results/NGS0001filtered.vcf -b ~/assignment_ngs/data/annotation.bed > ~/assignment_ngs/results/NGS0001_filtered_intersect.vcf
bgzip ~/assignment_ngs/results/NGS0001_filtered_intersect.vcf
tabix -p vcf ~/assignment_ngs/results/NGS0001_filtered_intersect.vcf.gz
rm ~/assignment_ngs/results/NGS0001.vcf
rm ~/assignment_ngs/results/NGS0001.vcf.gz

## 2.5 Variant Annotation and Prioritisation

#Annovar

#Annotation
# changing to the directory containing annovar for easier access to its commands
cd ~/annovar
# converting the vcf file into an annovar input file
./convert2annovar.pl -format vcf4 ~/assignment_ngs/results/NGS0001_filtered_intersect.vcf.gz > ~/assignment_ngs/results/NGS0001_filtered_intersect.avinput
# running annovar using its hg19 database to annotate the converted vcf file based on the hg19 database
# the -protocol and -operation options are used to add extra information based on different databases: refGene for RefSeq and ensGene for ensembl which both includes gene and functional information, clinvar_20180603 for ClinVar which includes clinical information, exac03 for ExAC which contains population frequencies, dbnsfp31a_interpro for dbNSFP and InterPro which contains predictions of clinical significance and avsnp147 for dbSNP which contains rsIDs.  
./table_annovar.pl ~/assignment_ngs/results/NGS0001_filtered_intersect.avinput humandb/ -buildver hg19 -out ~/assignment_ngs/results/NGS0001_filtered_intersect -remove -protocol refGene,ensGene,clinvar_20180603,exac03,dbnsfp31a_interpro,avsnp147 -operation g,g,f,f,f,f -otherinfo -nastring . -csvout
cd ~

#Filtering for exonic variants not in dbSNP
# $6 == "\"exonic\"" filters the 6th column for any exonic variants
# $30 == "." filters the 30th column (dbSNP) for any variants not present in the database
awk -F',' '$6 == "\"exonic\"" && $30 == "."' ~/assignment_ngs/results/NGS0001_filtered_intersect.hg19_multianno.csv > ~/assignment_ngs/results/NGS0001_annovar_exonic_notindbSNP.csv
cd ~

#snpEff
# running snpEff using java, with the -Xmx16g option to give it 16GB of RAM, as conda by default doesn't give it enough RAM to complete
# annotating the vcf file based on the downloaded hg19 genome
java -Xmx16g -jar /home/ubuntu/anaconda3/pkgs/snpeff-4.3.1t-0/share/snpeff-4.3.1t-0/snpEff.jar ann hg19 ~/assignment_ngs/results/NGS0001_filtered_intersect.vcf.gz > ~/assignment_ngs/results/NGS0001_snpEff_anno.vcf

#snpSift
# filtering the snpEff output for any exonic variants, using (ANN[*].EFFECT =~ 'exon') which naturally includes all types of coding variants, such as missense, nonsense, frameshift, etc.
# filtering for any variants not seen in dbSNP using (ID = '.') as a . is used to represent the empty space
java -Xmx4g -jar /home/ubuntu/anaconda3/pkgs/snpsift-4.3.1t-hdfd78af_3/share/snpsift-4.3.1t-3/SnpSift.jar filter "(ANN[*].EFFECT =~ 'exon') & (ID = '.')" ~/assignment_ngs/results/NGS0001_snpEff_anno.vcf > ~/assignment_ngs/results/NGS0001_snpEff_filtered.vcf

#Clearing all unnecessary files
rm ~/assignment_ngs/results/NGS0001_filtered_intersect.vcf
rm ~/assignment_ngs/results/NGS0001_filtered_intersect.vcf.gz
rm ~/assignment_ngs/results/NGS0001_filtered_intersect.avinput
rm ~/assignment_ngs/results/NGS0001_filtered_intersect.hg19_multianno.csv
rm ~/assignment_ngs/results/NGS0001_snpEff_anno.vcf
