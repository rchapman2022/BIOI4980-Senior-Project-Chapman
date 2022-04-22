# Set the working directory to the location of this R script
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Import the required libraries
library("rhdf5")
library("tximport")
library("DESeq2")
library("GenomicFeatures")
library("RMariaDB")
library("biomaRt")
library("org.Hs.eg.db")
library("pathfindR")
library("EnhancedVolcano")
library("dplyr")
library("tidyverse")

# This option allows for an infinite number of overlaps
# in ggplot2 plots. This is necessary for later
# use.
options(ggrepel.max.overlaps=Inf)

# Read in the sample names and abundance TSV files
samples <- read.table(file.path("Subset1Samples.txt"), header=TRUE)
files <- file.path(paste0(samples$sample, "_abundance.tsv"))
names(files) <- paste0(samples$sample, "_abundance.tsv")

# Creates a transcript DB from the ensembl released
# used as the kallisto reference
txdb <- makeTxDbFromEnsembl("Homo sapiens", release=105)
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- AnnotationDbi::select(txdb, k, "GENEID", "TXNAME")

# Imports the transcript quantification data
# from the files using tximport
txi.kallisto <- tximport(files, type = "kallisto", tx2gene = tx2gene, ignoreAfterBar = TRUE, ignoreTxVersion = TRUE)

# Uses the imported tx quantification data to create a 
# DESeq2 object and conducts a differential expression analysis.
ddsTxi <- DESeqDataSetFromTximport(txi.kallisto, colData = samples, design = ~ condition)
dds <- DESeq(ddsTxi)
res <- results(dds)
res.df <- data.frame(res)

# Creates an ensembl gene biomart for the reference assembly
# used in the analysis. This will be used to annotate
# the differential expression results.
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

# Grabs the differentially upregulated genes from the DESeq2
# results. Any gene with a log2FoldChange > 0 was considered
# upregulated and a p-value cutoff of 0.05 was enforced
upRegulatedGenes <- res.df[res.df$log2FoldChange > 0 & res.df$padj < 0.05, ]

# Grabs the gene ids of the upregulated genes and uses them
# to query the ensembl biomart and retrieve the ensembl_gene_id,
# gene name, and entrez_gene_id
upRegulatedGenes$geneid <- rownames(upRegulatedGenes)
upregGeneNames <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "entrezgene_id"), values = upRegulatedGenes$geneid, mart = ensembl)
colnames(upregGeneNames) <- c("geneid", "geneName", "entrez_gene_id")

# Adds the retrieved attributes to the dataframe containing the
# upregulated genes
upReg <- merge(upRegulatedGenes, upregGeneNames)

# Saves the upregulated genes to file
write.csv(upReg, "Subset1_up_regulated_genes.csv")

# Grabs the differentially downregulated genes from the DESeq2
# results. Any gene with a log2FoldChange < 0 was considered
# downregulated and a p-value cutoff of 0.05 was enforced
downRegulatedGenes <- res.df[res.df$log2FoldChange < 0 & res.df$padj < 0.05, ]

# Grabs the gene ids of the downregulated genes and uses them
# to query the ensembl biomart and retrieve the ensembl_gene_id,
# gene name, and entrez_gene_id
downRegulatedGenes$geneid <- rownames(downRegulatedGenes)
downregGeneNames <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "entrezgene_id"), values = downRegulatedGenes$geneid, mart = ensembl)
colnames(downregGeneNames) <- c("geneid", "geneName", "entrez_gene_id")

# Adds the retrieved attributes to the dataframe containing the
# downregulated genes
downReg <- merge(downRegulatedGenes, downregGeneNames)

# Adds the retrieved attributes to the dataframe containing the
# downregulated genes
write.csv(downReg, "Subset1_down_regulated_genes.csvs")



# Generates a Volcano Plot to represent the population of 
# differentially upregulated genes. The top 10 most
# upregulated genes based on log2FoldChange are labeled
upToMark <- head(arrange(upReg,desc(log2FoldChange)), n = 10)$geneName
EnhancedVolcano(upReg, lab=upReg$geneName, selectLab = upToMark, x= "log2FoldChange", y = 'padj',  FCcutoff = 0, pCutoff = 0.05, legendPosition = 'none', pointSize = 2, col=c('black', 'black', 'black', 'cornflowerblue'), labSize = 6, drawConnectors = TRUE, colConnectors = 'black', title = "Differentially Upregulated Genes in CRC Tissue for Subset 1") + coord_flip()


# Generates a Volcano Plot to represent the population of 
# differentially downregulated genes. The top 10 most
# downregulated genes based on log2FoldChange are labeled
downToMark <- head(arrange(downReg,log2FoldChange), n = 10)$geneName
EnhancedVolcano(downReg, lab=downReg$geneName, selectLab = downToMark, x= "log2FoldChange", y = 'padj',  FCcutoff = 0, pCutoff = 0.05, legendPosition = 'none', pointSize = 3, col=c('black', 'black', 'black', 'tomato1'), labSize = 6, drawConnectors = TRUE, colConnectors = 'black', title = "Differentially Downregulated Genes in CRC Tissue for Subset 1") + coord_flip()


# Grabs all of the differentially expressed genes from the
# DESeq2 results and uses the geneid to find the ensembl_gene_id, gene name,
# and entrezgene_id attributes for each.
allgenes <- res.df
allgenes$geneid <- rownames(allgenes)
allgenesGeneNames <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "entrezgene_id"), values = allgenes$geneid, mart = ensembl)
colnames(allgenesGeneNames) <- c("geneid", "geneName", "entrez_gene_id")

# Addes the new queried attributes to the dataframe with the 
# rest of the data.
allGenesFull <- merge(allgenes, allgenesGeneNames)

# Grabs the geneName, log2FoldChange, and adjusted p-value to be
# of the genes and creates a new dataframe. Any data with a NA
# adjusted p-value is dropped
allPRFinput <- allGenesFull[, c("geneName", "log2FoldChange", "padj")]
allPRFinput<- allPRFinput %>% drop_na(padj)

# PathFindR is run using the geneName, log2FoldChagne, and adjusted p-value of 
# all differentially expressed genes. A p-value cutoff of 0.05 is enforced to filter
# for significant differential expression.
allPRF <- run_pathfindR(allPRFinput, output_dir = "AllGenesPathways", p_val_threshold = 0.05)

# Creates an enrichment plot showing the top 15 enriched pathways
enrichment_chart(result_df = allPRF, top_terms = 10)

# Saves the data workspace
save.image("Subset1-image.RData")
