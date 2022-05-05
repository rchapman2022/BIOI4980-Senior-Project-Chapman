# BIOI4980-Senior-Project-Chapman

<p align="center">
  <img style="background-color: rgb(300, 300, 300);" src="https://github.com/rchapman2022/BIOI4980-Senior-Project-Chapman/blob/main/BIOI4980-Workflow.png">
</p>

The goal of this project is to compare the microbiome composition and gene expression of Colorectal Cancer in hopes of identifying precursory mechanistic steps for disease pathogenesis. To accomplish this, both a metagenomic and transcriptomic analysis were conducted using data from published sources:

**Metagenomic data:** 
    Debesa-Tur, G., Pérez-Brocal, V., Ruiz-Ruiz, S. et al. Metagenomic analysis of formalin-fixed paraffin-embedded tumor and normal mucosa reveals differences in the microbiome of colorectal cancer patients. Sci Rep 11, 391 (2021). https://doi.org/10.1038/s41598-020-79874-y

**Transcriptomic data:**
    Wu, SM., Tsai, WS., Chiang, SF. et al. Comprehensive transcriptome profiling of Taiwanese colorectal cancer implicates an ethnic basis for pathogenesis. Sci Rep 10, 4526 (2020). https://doi.org/10.1038/s41598-020-61273-y


## Environment Set-up

A separate analysis and visualization environment are included because the Graphlan package and its associated packages (used for visualization) conflicted with other packages in the analysis environment. Thus, YAML files for both conda environments can be found in the **Environments** directory.

The following code can be used to create pre-configured conda environments:

```
conda env create -f analysis-env.yml

conda env create -f visualization-env.yml
```

## Metagenomics Analysis

The metagenomics analysis consists of 5 distinct steps:
1. Visualization of the QC data using FASTQC
2. Quality Trimming using Trimmomatic
3. Host read removal via Bowtie2
4. Taxonomic Classification using MetaPhlAn3
5. Metabolic Pathway analysis using HUMAnN3

A bash script has been included to automate the process of running this analysis. This script can be found in the **Code** directory. Also included with this script is a text file containing sample SRA accession numbers (PRJEB34333-Accessions.txt) to be downloaded as well as two other sample files (Metagenomics-Healthy-Samples.txt and Metagenomics-Tumor-Samples.txt) designating which samples are tumor and healthy tissue. These files are required for the script to work properly. Thus, place **both the script and sample accession text files in the same directory** when running.

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
    ├── rawData/ - Contains the raw and intermediate analysis data
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

Results of the Taxonomic Classification analysis can be found in the "Combined-abundances.txt" file for both tissue types.

Results of the Metabolic Pathway analysis can be found in the .tsv files in the MetabolicPathwayResults folder.

### Visualization
Once the taxonomic classification and metabolic pathway analyses have been complete, the data can be visualized using the Graphlan tool. This creates a circular representation of microbiome composition. A bash script has been included to automate this analysis. **The script is dependent upon the taxonomic classification results (SEE SECTION ABOVE).** 

The script can be run in the following manner (*File paths are based on metagenomics anlaysis output file structure (SEE SECTION ABOVE)*):
```
conda activate bioi4980-visualization-env


# ENSURE THAT YOU ARE IN THE TaxClassResults Directory

source metagenomic-visualization.sh
```

The script will generate .png images inside of the TumorSamples/ and HealthySamples/ Directories. Both a figure and a figure legend will be included:
```
.
└── metagenomics-analysis/
    ├── TaxClassResults/ - Contains the results of the Taxonomic Classification Analysis/
    │   ├── HealthyResults/
    │   │   ├── HealthySamples.png - the circular graph
    │   │   └── HealthySamples_legend.png - an external legend for the figure
    │   └── TumorResults/
    │       ├── TumorSamples.png - the circular graph
    │       └── TumorSamples_legend.png - an exerternal legend for the figure
```

## Transcriptomic Analysis

The transcriptomic analysis was separated into two distinct steps:
1. Read Alignement and Transcript Quantification
2. Differential Expression and Pathway Analysis

Additionally, it is important to note that, due to the size of the dataset and limitations regarding memory/storage, three distinct subsets consisting of 5 paired CRC and healthy tissue samples were analyzed separately from one another. The provided scripts take this into account.

### Read Alignment and Transcript Quantification

Transcript Quantification is performed by kallisto, and consists of both a pseudoalignment and quantification step. A bash script has been included to automate the process of running this analysis. Because this analysis consisted of three subsets of data, the bash script will create a file structure to accommodate this. This script can be found in the **Code** directory. Also included are three text files containing SRA accession numbers for the samples in each subset. These files are required for the script to work properly. Thus, place **the script and SRA accession files all in the same directory** when running.

The script can be run in the following manner:
```
conda activate bioi4980-analysis-env

mkdir transcriptomic-analysis/
# move the script and all sample files into this folder

source transcriptomic-analysis.sh
```

The script will create the following file structure:
```
.
└── transcriptomic-analysis/
    ├── Subset1/
    │   ├── abundances/ - Contain the transcript abundances for the subset
    │   └── Folders containing the intermediate files for each sample
    ├── Subset2/
    │   ├── abundances/ - Contain the transcript abundances for the subset
    │   └── Folders containing the intermediate files for each sample
    ├── Subset3/
    │   ├── abundances/ - Contain the transcript abundances for the subset
    │   └── Folders containing the intermediate files for each sample
    ├── Homo_sapiens.GRCh38.cdna.all.fa.gz - the reference transcriptome used for this analysis
    └── hg38idx - the kallisto index used for alignment
```

Within each subset, the abundances directory contains the primary output of the transcript quantification analysis. These files will be used for the following differential expression analysis (SEE SECTION BELOW).

### Differential Expression and Pathway Analysis

The Differential Expression and Pathway analyses were performed using R. For each sample, an R script has been provided to automate the process of running these analyses. The scripts can be opened in R studio and walked through step by step to generate the results. As well, a sample list for each subset is also included. Each of these files can be found in the respective **Subset** directory of the **Differential Expression Analysis** directory. As well, an R script to identify upregulated and downregulated genes common to the subsets is included in the **Differential Expression Analysis** directory (FindCommonGenes.R). 

#### Dependencies
Before using these scripts, the dependencies need to be installed:
```
BiocManager::install("rhdf5")
BiocManager::install("tximport")
BiocManager::install("DESeq2")
BiocManager::install("GenomicFeatures")
BiocManager::install("biomaRt")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("EnhancedVolcano")

install.packages("pathfindR")
install.packages("RMariaDB")
install.packages("dplyr")
install.packages("tidyverse")
```

Next, it is important to note that these scripts **requires a specific file structure**:
```
.
├── Differential-Expression analysis/
    ├── Subset1/
    │   ├── Sample-abundance.tsv - the transcript abundances generated in the previous step
    │   ├── Subset1Analysis.R - the R script for running the analysis
    │   └── Subset1Samples.txt - a text file containing the sample file names and conditions
    ├── Subset2/
    │   ├── Sample-abundance.tsv - the transcript abundances generated in the previous step
    │   ├── Subset2Analysis.R - the R script for running the analysis
    │   └── Subset2Samples.txt - a text file containing the sample file names and conditions
    ├── Subset3/
    │   ├── Sample-abundance.tsv - the transcript abundances generated in the previous step
    │   ├── Subset3Analysis.R - the R script for running the analysis
    │   └── Subset3Samples.txt - a text file containing the sample file names and conditions
    └── FindCommonGenes.R
```
The file organization of the **Differential Expression Analysis** directory of this repository mimics the correct layout. 
- Each subset's analysis should be conducted in its own directory named *Subset#*. 
- The script for that subset, the sample list for that subset, and the abundance.tsv files produced in the transcript quantification step should be placed in the directory. 
- The script to identify the common differentially expressed gene identification should be placed in the parent directory containing the subset folders.

#### Individual Sample Differential Expression
The individual subset scripts should be run **before** the FindCommonGenes.R script, as the files created serve as the input for this script.

If opened in RStudio, the R scripts can be walked through to conduct the analysis. 

After each analysis is complete, the following files will be added to each subset directory:
- a .csv file containing differentially upregulated genes
- a .csv file containing differentially downregulated genes
- volcano plots for the differentially upregulated and downregulated genes 
- a folder containing the pathfindR results will be added to the subset directory.
- a pathway enrichment plot for the top 15 enriched pathways

```
.
└── Differential-Expression analysis/
    └── Subset#/
        ├── Sample-abundance.tsv - the transcript abundances generated in the previous step
        ├── Subset1Analysis.R - the R script for running the analysis
        ├── Subset1Samples.txt - a text file containing the sample file names and conditions
        ├── AllGenesPathways/
        ├── Subset1_up_regulated_genes.csv
        ├── Subset1_down_regulated_genes.csv
        ├── Subset1UpRegVolcano.png
        ├── Subset1DownRegVolcano.png
        └── Subset1EnrichedPathways.png
```
#### Identify Common Genes
The FindCommonGenes.R script will collect the upregulated and downregulated genes from each subset and identify those in common between the samples

Again, if opened in RStudio, the R script can be walked through to conduct the analysis. 

After this script is complete a list of the common, differentially upregulated and downregulated genes will be added to the parent directory:
```
.
└── Differential-Expression analysis/
    ├── Subset1/
    ├── Subset2/
    ├── Subset3/
    ├── FindCommonGenes.R
    ├── Compared_Down_regulated_genes.csv
    └── Compared_Up_regulated_genes.csv
```
## Dependency References

- [SRA Tools](https://github.com/ncbi/sra-tools)
- [FASTQC](https://github.com/s-andrews/FastQC)
- [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
- [MetaPhlAn](https://github.com/biobakery/MetaPhlAn)
- [HUMAnN3](https://github.com/biobakery/humann)
- [Graphlan](https://github.com/biobakery/graphlan)
- [export2graphlan](https://github.com/SegataLab/export2graphlan)
- [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
- [Kallisto](https://pachterlab.github.io/kallisto/manual)
- [rhdf5](https://github.com/grimbough/rhdf5)
- [tximport](https://github.com/mikelove/tximport)
- [GenomicFeatures](https://bioconductor.org/packages/release/bioc/html/GenomicFeatures.html)
- [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)
- [biomaRt](https://bioconductor.org/packages/release/bioc/html/biomaRt.html)
- [org.Hs.eg.db](https://bioconductor.org/packages/release/data/annotation/html/org.Hs.eg.db.html)
- [EnhancedVolcano](https://github.com/kevinblighe/EnhancedVolcano)
- [pathfindR](https://github.com/egeulgen/pathfindR)
- [RMariaDB](https://github.com/r-dbi/RMariaDB)
- [Tidyverse](https://www.tidyverse.org/)
- [dplyr](https://dplyr.tidyverse.org/)



