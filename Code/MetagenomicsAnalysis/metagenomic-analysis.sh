#!/bin/bash

# Create a directory for the raw reads and intermediate data
mkdir rawData
cd rawData

# Download the data from ENA (Project accession: PRJEB34333)
prefetch --option-file ../../PRJEB34333-Accessions.txt
for file in *
do
    echo $file
    fastq-dump --gzip --split-files --outdir $file $file
done

# Extact fastq files for each sample in gzipped format
for file in *
do
    echo $file
    fastq-dump --split-files --gzip --outdir $file $file

# Perform an initial QC Check on the metagenomics samples
for file in *
do 
    fastqc $file/*
done

# Perform Quality Trimming on the Samples using trimmmomatic
# Includes the following options:
#   - Remove Illumina Adapters
#   - Uses a sliding window to trim regions with an average quality 
#     > 20
#   - Removes any reads with < 50 base pairs after the previous operations were performed
for file in *
do
    base=$(basename $file/*_1.fastq.gz _1.fastq.gz)
    echo $base
    trimmomatic PE -threads 6 $file/*_1.fastq.gz $file/*_2.fastq.gz $file/$base"_1.trim.fq" $file/$base"_1.untrim.fq" $file/$base"_2.trim.fq" $file/$base"_2.untrim.fq" ILLUMINACLIP:../../trimmomatic/adapters/TruSeq3-PE-2.fa:2:30:10:8:TRUE SLIDINGWINDOW:4:20 MINLEN:50; fastqc $file/*.trim.fq
    echo "------------"
done


# Download the human reference genome 
mkdir ../../genomicReferences
cd ../../genomicReferences
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.29_GRCh38.p14/GCA_000001405.29_GRCh38.p14_genomic.fna.gz
gunzip GCA_000001405.29_GRCh38.p14_genomic.fna.gz

# Create a bowtie2 index for the genome
bowtie2-build GCA_000001405.29_GRCh38.p14_genomic.fna hg38idx

cd ../metagenomicAnalysis/rawData

# Perform Host Read removal using bowtie2
# The following options are provided:
#   - very-sensitive
#   - write unaligned reads to file
for file in *
do
    base=$(basename $file/*_1.fastq.gz _1.fastq.gz)
    echo $base
    bowtie2 -x ../../genomicReferences/hg38idx -1 $file/$base"_1.trim.fq" -2 $file/$base"_2.trim.fq" --very-sensitive --un-conc $file/$base".unaligned-human.fq" --threads 6 -S $file/$base"-human.sam"
    echo "----------------"
done


# Perform Taxonomic Classification using MetaPHlAn3
for file in *
do
    base=$(basename $file/*_1.fastq.gz _1.fastq.gz)
    echo $base
    metaphlan $file/$base".unaligned-human.1.fq",$file/$base".unaligned-human.2.fq" --bowtie2out $file/$base"-mg-bt2.bt2" --nproc 8 --input_type fastq -o $file/$base"-metagenome.txt" --bowtie2db ~/metaphlanDB -x mpa_v30_CHOCOPhlAn_201901
    echo "-------------------"
done

# Perform Metabolic Pathway analysis using HUMAnN3. The tools requies that 
# paired-end read files be concatenated into a single file, thus this is performed first. 
for file in *
do 
    base=$(basename $file/*_1.fastq.gz _1.fastq.gz)
    echo $base
    cat $file/$base".unaligned-human.1.fq" $file/$base".unaligned-human.2.fq" > $file/$base".unaligned-human.all.fq"
    humann --input $file/$base".unaligned-human.all.fq" --output $file/ --taxonomic-profile $file/$base"-metagenome.txt"
    echo "-------------------------------"
done

cd ..

# Move Taxonmic Classification results to a results directory
mkdir TaxClassResults
cd TaxClassResults
mkdir TumorResults
while read sample
do
    cp ../rawData/$sample/$sample"-metagenome.txt" TumorResults/
done < ../Metagenomics-Tumor-Samples.txt

mkdir HealthyResults
while read sample
do
    cp ../rawData/$sample/$sample"-metagenome.txt" HealthyResults/
done < ../Metagenomics-Healthy-Samples.txt

# Merge Taxonomic Classifications for healthy and tumor samples 
merge_metaphlan_tables.py TumorResults/*-metagenome.txt > TumorResults/Combined-abundances.txt
merge_metaphlan_tables.py HealthyResults/*-metagenome.txt > HealthyResults/Combined-abundances.txt


cd ..

# Move Metabolic Pathway Analysis results to a results directory
mkdir MetabolicPathwayResults
cd MetabolicPathwayResults
mkdir TumorSamples
while read sample
do
    cp ../rawData/$sample/$sample".unaligned-human."*".tsv" TumorSamples/
done < ../Metagenomics-Tumor-Samples.txt

mkdir HealthySamples
while read sample
do
    cp ../rawData/$sample/$sample".unaligned-human."*".tsv" HealthySamples/
done < ../Metagenomics-Healthy-Samples.txt

# Merge Metabolic Pathway Analysis results for healthy and tumor samples
humann_join_tables -i HealthySamples/ -o healthyTissue_genefamilies.tsv --file_name genefamilies
humann_join_tables -i HealthySamples/ -o healthyTissue_pathcoverage.tsv --file_name pathcoverage
humann_join_tables -i HealthySamples/ -o healthyTissue_pathabundance.tsv --file_name pathabundance
humann_join_tables -i TumorSamples/ -o tumorTissue_genefamilies.tsv --file_name genefamilies
humann_join_tables -i TumorSamples/ -o tumorTissue_pathcoverage.tsv --file_name pathcoverage
humann_join_tables -i TumorSamples/ -o tumorTissue_pathabundance.tsv --file_name pathabundance

# Renormalize the genefamilies table to be based on relative abundance
humann_renorm_table -i healthyTissue_genefamilies.tsv -o healthyTissue_genefamilies-relab.tsv --units relab
humann_renorm_table -i tumorTissue_genefamilies.tsv -o tumorTissue_genefamilies-relab.tsv --units relab