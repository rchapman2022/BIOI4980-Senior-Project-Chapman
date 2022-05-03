#!/bin/bash

# Create a kallisto Index. First, download the latest ensembl Hg38 CDNA assembly.
# Then create hte index using kallisto's built in functionality.
wget http://ftp.ensembl.org/pub/release-105/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
kallisto index -i hg38idx Homo_sapiens.GRCh38.cdna.all.fa


# Subset 1
##########

# Create Directory
mkdir Subset1
cd Subset1

# Download the SRA samples for this subset and extract the fastq files
prefetch --option-file ../TranscriptomicSubset1AccList.txt
fasterq-dump --split-files *

# Create directories for each sample and copy move the corresponding reads into the directory
mkdir 191T
mv SRR5581212/ 191T

mkdir 191N
mv SRR5581213* 191N

mkdir 145T
mv SRR5581016* 145T

mkdir 145N
mv SRR5581015* 145N

mkdir 151T
mv SRR5581012* 151T

mkdir 151N
mv SRR5581011* 151N

mkdir 17T
mv SRR5580997* 17T

mkdir 17N
mv SRR5580996* 17N

mkdir 7T
mv SRR5581157* 7T/

mkdir 7N
mv SRR5581156* 7N/

# Run FASTQC on the samples
for file in *
do
    base=$(basename $file/SRR*_1.fastq _1.fastq)
    fastqc $file/*.fastq
done

# Data is preprocessed, thus no trimming is required.

# Align and quantify the reads using kallisto
for file in * 
do
    base=$(basename $file/*_1.fastq _1.fastq)
    echo $base
    kallisto quant -i ../hg38idx -o $file/ -t 8 $file/$base"_1.fastq" $file/$ base"_2.fastq"
    echo "--------------------"
done

# Rename and move files into an abundance folder
mkdir abundances
for file in *
do
    cp $file/abundance.tsv ./abundances/$file"_abundance.tsv"
done


cd ..

# Subset 2
##########

# Create Directory
mkdir Subset2
cd Subset2

# Download the SRA samples for this subset and extract the fastq files
prefetch --option-file ../TranscriptomicSubset2AccList.txt
fasterq-dump --split-files *

# Create directories for each sample and copy move the corresponding reads into the directory
mkdir 187N
mv SRR5581219* 187N

mkdir 187T
mv SRR5581218* 187T

mkdir 135T
mv SRR5580988* 135T

mkdir 135N
mv SRR5580989* 135N

mkdir 182T
mv SRR5581325* 182T

mkdir 182N
mv SRR5581324* 182N

mkdir 39T
mv SRR5581104* 39T

mkdir 39N
mv SRR5581103* 39N

mkdir 23T
mv SRR5581237* 23T

mkdir 23N
mv SRR5581238* 23N

# Run FASTQC on the samples
for file in *
do
    base=$(basename $file/SRR*_1.fastq _1.fastq)
    fastqc $file/*.fastq
done

# Data is preprocessed, thus no trimming is required.

# Align and quantify the reads using kallisto
for file in * 
do
    base=$(basename $file/*_1.fastq _1.fastq)
    echo $base
    kallisto quant -i ../hg38idx -o $file/ -t 8 $file/$base"_1.fastq" $file/$ base"_2.fastq"
    echo "--------------------"
done

# Rename and move files into an abundance folder
mkdir abundances
for file in *
do
    cp $file/abundance.tsv ./abundances/$file"_abundance.tsv"
done


cd ..

# Subset 3
##########

# Create Directory
mkdir Subset3
cd Subset3

# Download the SRA samples for this subset and extract the fastq files
prefetch --option-file ../TranscriptomicSubset2AccList.txt
fasterq-dump --split-files *

# Create directories for each sample and copy move the corresponding reads into the directory
mkdir 29N
mv SRR5581236* 29N

mkdir 29T
mv SRR5581235* 29T

mkdir 2N
mv SRR5581154* 2N

mkdir 2T
mv SRR5581155* 2T

mkdir 16T
mv SRR5581255* 16T

mkdir 16N
mv SRR5581256* 16N

mkdir 183N
mv SRR5581082* 183N

mkdir 183T
mv SRR5581083* 183T

mkdir 202N
mv SRR5581268* 202N

mkdir 202T
mv SRR5581269* 202T

# Run FASTQC on the samples
for file in *
do
    base=$(basename $file/SRR*_1.fastq _1.fastq)
    fastqc $file/*.fastq
done

# Data is preprocessed, thus no trimming is required.

# Align and quantify the reads using kallisto
for file in * 
do
    base=$(basename $file/*_1.fastq _1.fastq)
    echo $base
    kallisto quant -i ../hg38idx -o $file/ -t 8 $file/$base"_1.fastq" $file/$ base"_2.fastq"
    echo "--------------------"
done

# Rename and move files into an abundance folder
mkdir abundances
for file in *
do
    cp $file/abundance.tsv ./abundances/$file"_abundance.tsv"
done