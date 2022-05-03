# BIOI4980-Senior-Project-Chapman

<p align="center">
  <img style="background-color: rgb(300, 300, 300);" src="https://github.com/rchapman2022/BIOI4980-Senior-Project-Chapman/blob/main/BIOI4980-Workflow.png">
</p>

The goal of this project is to compare the microbiome composition and gene expression of Colorectal Cancer in hopes of identifying precursory mechanistic steps for pathogenesis. To accomplish this, both a metagenomic and transcriptomic analysis were conducted using data from published sources:

**Metagenomic data:** 
    Debesa-Tur, G., Pérez-Brocal, V., Ruiz-Ruiz, S. et al. Metagenomic analysis of formalin-fixed paraffin-embedded tumor and normal mucosa reveals differences in the microbiome of colorectal cancer patients. Sci Rep 11, 391 (2021). https://doi.org/10.1038/s41598-020-79874-y

**Transcriptomic data:**
    Wu, SM., Tsai, WS., Chiang, SF. et al. Comprehensive transcriptome profiling of Taiwanese colorectal cancer implicates an ethnic basis for pathogenesis. Sci Rep 10, 4526 (2020). https://doi.org/10.1038/s41598-020-61273-y


## Environment Set-up

A separate analysis and visualization eviroment are included because the Graphlan and its associated packages conflicted with other packages in the analysis environment. Thus, YAML files for both conda enviroments can be found in the **Environments** directory.

The following code can be used to create pre-configured conda environments:

```
conda env create -f analysis-env.yml

conda env create -f visualization-env.yml
```

## Metagenomics Analysis

The metagenomics anlaysis consists of 5 distinct steps:
1. Visualization of the QC data using FASTQC
2. Quality Trimming using Trimmomatic
3. Host read removal via Bowtie2
4. Taxonomic Classification using MetaPhlAn3
5. Metabolic Pathway analysis using HUMAnN3

A bash script has been included to automate the process of running this analysis. This script can be found in the **Code** directory. Also included with this script is a text file containing sample SRA accession numbers (PRJEB34333-Accessions.txt) to be downloaded as well as two other sample files (Metagenomics-Healthy-Samples.txt and Metagenomics-Tumor-Samples.txt) designating which samples are tumor and healthy tissue. These file is required for the script to work properly. Thus, place **both the script and sample accession text files in the same directory** when running.

The script can be run in the following manner:
```
conda activate bioi4980-analysis-env

mkdir metagenomics-analysis/
# move the script and all sample files into this folder

source metagenomic-analysis.sh
```

The script will generate the following file structure:
```
.
└── metagenomics-analysis/
    ├── rawData - Contains the raw and intermediate analysis data
    ├── TaxClassResults/ - Contains the results of the Taxonomic Classification Analysis/
    │   ├── HealthyResults/
    │   │   ├── Combined-abundances.txt - contains merged abundances from all healthy samples
    │   │   └── individual sample abundances
    │   └── TumorResults/
    │       ├── Combined-abundances.txt - contains merged abundances from all tumor
    │       └── individual sample abundances
    └── MetabolicPathwayResults/ - Contains the results of the Metabolic Pathway Analysis/
        ├── TumorSamples/
        │   └── individual tumor sample results
        ├── HealthySamples/
        │   └── individual healthy sample results
        ├── healthyTissue_genefamilies.tsv
        ├── healthyTissue_genefamilies-relab.tsv
        ├── healthyTissue_pathcoverage.tsv
        ├── healthyTissue_pathabundance.tsv
        ├── tumorTissue_genefamilies.tsv
        ├── tumorTissue_genefamilies-relab.tsv
        ├── tumorTissue_pathcoverage.tsv
        └── tumorTissue_pathabundance.tsv
```

Results of the Taxonomic Classification analysis can be found in the "Combined-abundances.txt" file for both samples.

Results of the Metabolic Pathway analysis can be found in the .tsv files in the results folder.


## Transcriptomic Analysis

The transcriptomic analysis was separated into to distinct steps:
1. Read alignement and Transcript quantification
2. Differential Expression and Pathway Analysis

Additionally, it is important to note that, due to the size of the dataset and limitations regarding memory/storage, three distinct subsets consisting of 5 paired CRC and healthy tissue samples were analyzed separate from one another. The provided scripts take this into account.

### Transcript Quantification

Transcript Quantification is performed by kallisto, and consists of both a pseudoalignment and quantification step. A bash script has been included to automate the process of running this analysis. Because this analysis consisted of three subsets of data, the bash script will create a file structure to accommodate this. This script can be found in the **Code** directory. Also included are three text files containing SRA accession numbers for the samples in each subset. These files are required for the script to work properly. Thus, place **the script and SRA accession files all in the same directory** when running.

The script can be run in the following manner:
```
conda activate bioi4980-analysis-env

source transcriptomic-analysis.sh
```

Results from this sript will include transcript abundance files in .tsv format for each sample analyzed. Additionally, the names of the files will be replaced with the sample number and tissue type to make this data easier for a human to work with.

### Differential Expression and Pathway Analysis

The Differential Expression and Pathway analysis was performed using R. For each sample, an R script has been provided to automated the process of running these analyses. The scripts can be opened in R studio and walked through step by step to generate the results. As well, a sample list for each subset is also included Each of these files can be found in the respective **Subset** directory of the **Differential Expression Analysis** directory. As well, an R script to identify upregulated and downregulated genes common to the subsets is included in the **Differential Expression Analysis** directory. 

Before using these scripts, the dependencies need to be installed:
```
BiocManager::install("rhdf5")
BiocManager::install("tximport")
BiocManager::install("DESeq2")
BiocManager::install("GenomicFeatures")
BiocManager::install("biomaRt")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("EnhancedVolcano")
BiocManager::install("")

install.packages("pathfindR")
install.packages("RMariaDB")
install.packages("dplyr")
install.packages("tidyverse")
```

Next, tt is important to note that these scripts **requires a specific file organization**. The file organization of the **Differential Expression Analysis** directory mimics the correct layout. 
- Each subset's analysis should be conducted in its own directory named *Subset#*. 
- The script for that subset, the the sample list for that subset, and the abundance.tsv files produced in the transcript quantification step should be placed in the directory. 
- The script to identify the common differentially expressed gene identification should be placed in the parent directory containing the Subset folders.

