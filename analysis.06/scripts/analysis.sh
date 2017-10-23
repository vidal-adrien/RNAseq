#!/bin/sh
#Author: Vidal Adrien
#Dependencies: fastqc, seqtk, bwa, samtools, bedtools

# Study PRJEB10569:
#	- Run ERR990557 (Exp: ERX1071736, Sample: SAMEA3516071)
#	- Run ERR990558 (Exp: ERX1071737, Sample: SAMEA3516072)
#	- Run ERR990559 (Exp: ERX1071738, Sample: SAMEA3516073)
#	- Run ERR990560 (Exp: ERX1071739, Sample: SAMEA3516074)
#
# Tissue: Ovaries, Mutation: Rpp30, Replicates: Yes

#Current dir analysis.xx/scripts
cd ..

######-------------------------------------######
######    Raw data files (fastq/a)
######-------------------------------------######

mkdir sources 

#Downloading (1):
wget -P ./sources ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR990/ERR990557/ERR990557.fastq.gz
wget -P ./sources ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR990/ERR990558/ERR990558.fastq.gz
wget -P ./sources ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR990/ERR990559/ERR990559.fastq.gz
wget -P ./sources ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR990/ERR990560/ERR990560.fastq.gz

#Quality Control:
mkdir fastqc_samp_57 & fastqc -o fastqc_samp_57 ./sources/ERR990557.fastq.gz
mkdir fastqc_samp_58 & fastqc -o fastqc_samp_58 ./sources/ERR990558.fastq.gz
mkdir fastqc_samp_59 & fastqc -o fastqc_samp_59 ./sources/ERR990559.fastq.gz
mkdir fastqc_samp_60 & fastqc -o fastqc_samp_60 ./sources/ERR990560.fastq.gz

#Unzipping (2):
gunzip -c ./sources/ERR990557.fastq.gz > ./sources/ERR990557.fastq
gunzip -c ./sources/ERR990558.fastq.gz > ./sources/ERR990558.fastq
gunzip -c ./sources/ERR990559.fastq.gz > ./sources/ERR990559.fastq
gunzip -c ./sources/ERR990560.fastq.gz > ./sources/ERR990560.fastq

#Sampling (randomly) (3):
mkdir fastq
samp_size=8000000
seqtk sample -s100 ./sources/ERR990557.fastq $samp_size > ./fastq/ERR990557s.fastq
seqtk sample -s100 ./sources/ERR990558.fastq $samp_size > ./fastq/ERR990558s.fastq
seqtk sample -s100 ./sources/ERR990559.fastq $samp_size > ./fastq/ERR990559s.fastq
seqtk sample -s100 ./sources/ERR990560.fastq $samp_size > ./fastq/ERR990560s.fastq
rm ./sources/*.fastq #remove unzipped fastq files (large files).

#Reference sequence:
wget -P ./sources ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.01_FB2017_04/fasta/dmel-all-chromosome-r6.01.fasta.gz
gunzip -c ./sources/dmel-all-chromosome-r6.01.fasta.gz > ./sources/dmel-all-chromosome-r6.01.fasta
rm ./sources/dmel-all-chromosome-r6.01.fasta.gz
ref="./sources/dmel-all-chromosome-r6.01.fasta"

#Genome annotations:
wget -P ./sources ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.01_FB2017_04/gff/dmel-all-r6.01.gff.gz
gunzip -c ./sources/dmel-all-r6.01.gff.gz > ./sources/dmel-all-r6.01.gff
rm dmel-all-r6.01.gff.gz

######-------------------------------------######
######    Indexing - Aligning reads on the genome (4)
######-------------------------------------######

mkdir bwa

#~ Index the fasta file
bwa index $ref

#~ Run bwa:
bwa aln -t 4 $ref ./fastq/ERR990557s.fastq > bwa/ERR990557s.sai
bwa samse $ref bwa/ERR990557s.sai ./fastq/ERR990557s.fastq > bwa/ERR990557s.sam
#~  For ERR990557s

#~ Run bwa:
bwa aln -t 4 $ref ./fastq/ERR990558s.fastq > bwa/ERR990558s.sai
bwa samse $ref bwa/ERR990558s.sai ./fastq/ERR990558s.fastq > bwa/ERR990558s.sam
#~  For ERR990558s

#~ Run bwa:
bwa aln -t 4 $ref ./fastq/ERR990559s.fastq > bwa/ERR990559s.sai
bwa samse $ref bwa/ERR990559s.sai ./fastq/ERR990559s.fastq > bwa/ERR990559s.sam
#~  For ERR990559s

#~ Run bwa:
bwa aln -t 4 $ref ./fastq/ERR990560s.fastq > bwa/ERR990560s.sai
bwa samse $ref bwa/ERR990560s.sai ./fastq/ERR990560s.fastq > bwa/ERR990560s.sam
#~  For ERR990560s

#~  remove sai files:
rm bwa/*.sai

cd ./bwa
#~ Filter uniquely mapped reads:
grep XT:A:U ERR990557s.sam > ERR990557s_filt.sam
grep XT:A:U ERR990558s.sam > ERR990558s_filt.sam
grep XT:A:U ERR990559s.sam > ERR990559s_filt.sam
grep XT:A:U ERR990560s.sam > ERR990560s_filt.sam
#~ The tag XT:A:U means that the read maps to a unique position in the genome.
cd ..

#~ Additional filters + sam to sorted bam conversion:
samtools faidx $ref
samtools view -bS -F 1548 -q 30 -t $ref.fai ./bwa/ERR990557s_filt.sam | samtools sort - ./bwa/ERR990557s_filt_sorted #.bam
samtools view -bS -F 1548 -q 30 -t $ref.fai ./bwa/ERR990558s_filt.sam | samtools sort - ./bwa/ERR990558s_filt_sorted #.bam
samtools view -bS -F 1548 -q 30 -t $ref.fai ./bwa/ERR990559s_filt.sam | samtools sort - ./bwa/ERR990559s_filt_sorted #.bam
samtools view -bS -F 1548 -q 30 -t $ref.fai ./bwa/ERR990560s_filt.sam | samtools sort - ./bwa/ERR990560s_filt_sorted #.bam
#~  -F to filter out PCR duplicates, low quality reads, etc.  http://broadinstitute.github.io/picard/explain-flags.html
#~  -q is to filter out reads mapping with a low quality. http://maq.sourceforge.net/qual.shtml 


######-------------------------------------######
######    Counting reads aligning with genes (5)
######-------------------------------------######

gff2bed --do-not-sort < ./sources/dmel-all-r6.01.gff > ./sources/dmel-all-r6.01.bed #With sorting, this script fills up the system disk (for some reason) when operating and fails because it's limited on my machine.
#~ bedops tool.

samtools index ./bwa/ERR990557s_filt_sorted.bam
samtools index ./bwa/ERR990558s_filt_sorted.bam
samtools index ./bwa/ERR990559s_filt_sorted.bam
samtools index ./bwa/ERR990560s_filt_sorted.bam

bedtools multicov -s -bams ./bwa/ERR990557s_filt_sorted.bam ./bwa/ERR990558s_filt_sorted.bam ./bwa/ERR990559s_filt_sorted.bam ./bwa/ERR990560s_filt_sorted.bam -bed ./sources/dmel-all-r6.01.bed > coverage.txt
