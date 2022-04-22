#!/bin/bash

# Create Directories for each subset
mkdir s1
mkdir s2
mkdir s3

# Create a kallisto Index. First, download the latest ensembl Hg38 CDNA assembly.
# Then create hte index using kallisto's built in functionality.
wget http://ftp.ensembl.org/pub/release-105/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
kallisto index -i hg38idx Homo_sapiens.GRCh38.cdna.all.fa

# Data is preprocessed, thus no QC or trimming is required.

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
for file in *; do cp $file/abundance.tsv ./abundances/$file"_abundance.tsv"