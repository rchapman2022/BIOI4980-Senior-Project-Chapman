# Sets the working directory to the location of the R script 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# import necessary libraries
library(dplyr)

#Read in Up and Down regulated genes for Subset1
S1Up <- read.csv("Subset1/Subset1_up_regulated_genes.csv")
S1Down <- read.csv("Subset1/Subset1_down_regulated_genes.csv")

#Read in Up and Down regulated genes for Subset2
S2Up <- read.csv("Subset2/Subset2_up_regulated_genes.csv")
S2Down <- read.csv("Subset2/Subset2_down_regulated_genes.csv")

#Read in Up and Down regulated genes for Subset3
S3Up <- read.csv("Subset3/Subset3_up_regulated_genes.csv")
S3Down <- read.csv("Subset3/Subset3_down_regulated_genes.csv")

# Join the Subsets on gene name and place in a data frame containing
# the Log2FoldChange for each subset
commonUp <- select(S1Up, geneid, log2FoldChange, geneName) %>% left_join(select(S2Up, geneid, log2FoldChange, geneName), by=c("geneName", "geneid")) %>% left_join(select(S3Up, geneid, log2FoldChange, geneName), by=c("geneName", "geneid"))
colnames(commonUp) <- c("geneid", "S1Log2FC", "geneName", "S2Log2FC", "S3Log2FC")

#Remove any rows with null values (gene isn't found in 1+ sample)
commonUp <- commonUp[complete.cases(commonUp),]

# Add an column containing the average of the Log2Fold Change value
# for the three samples.
commonUp <- mutate(commonUp %>% rowwise(), avgLog2FC = rowMeans(cbind(S1Log2FC, S2Log2FC, S3Log2FC)))

# Create a .csv file output
write.csv(commonUp, file="Compared_Up_regulated_genes.csv")


# Join the Subsets on gene name and place in a data frame containing
# the Log2FoldChange for each subset
commonDown <- select(S1Down, geneid, log2FoldChange, geneName) %>% left_join(select(S2Down, geneid, log2FoldChange, geneName), by=c("geneName", "geneid")) %>% left_join(select(S3Down, geneid, log2FoldChange, geneName), by=c("geneName", "geneid"))
colnames(commonDown) <- c("geneid", "S1Log2FC", "geneName", "S2Log2FC", "S3Log2FC")

#Remove any rows with null values (gene isn't found in 1+ sample)
commonDown <- commonDown[complete.cases(commonDown),]

# Add an column containing the average of the Log2Fold Change value
# for the three samples.
commonDown <- mutate(commonDown %>% rowwise(), avgLog2FC = rowMeans(cbind(S1Log2FC, S2Log2FC, S3Log2FC)))

# Create a .csv file output
write.csv(commonDown, file="Compared_Down_regulated_genes.csv")
