#!/bin/bash

## This script is assuming that it is being run in the TaxClassResults created by 
## the metagenomic-analysis.sh script


# Tumor Tissue Visualization
############################
cd TumorResults

# Reformat the abundance file for graphlan
tail -n +2 Combined-abundances.txt | cut -f1,3- > Combined-abundances-reformat.txt

# Use the export2graphlan package to setup the annotaiton and tree files
export2graphlan.py --skip_rows 1 -i Combined-abundances-reformat.txt --tree Combined-abundances.tree.txt --annotation Combined-abundances.annot.txt --most_abundant 100 --abundance_threshold 1 --least_biomarkers 3 --annotations 5,6 --external_annotations 7 --min_clade_size 1 --title "Tumor Tissue Microbiome"

# Create an xml file for graphlan based on the annotations and tree files
graphlan_annotate.py --annot Combined-abundances.annot.txt Combined-abundances.tree.txt Combined-abundances.xml

# Generate the circular phylogenetic tree visualization
graphlan.py --dpi 300 Combined-abundances.xml TumorSamples.png --external_legends


cd ..

# Healthy Tissue Visualization
##############################
cd HealthyResults

# Reformat the abundance file for graphlan
tail -n +2 Combined-abundances.txt | cut -f1,3- > Combined-abundances-reformat.txt

# Use the export2graphlan package to setup the annotaiton and tree files
export2graphlan.py --skip_rows 1 -i Combined-abundances-reformat.txt --tree Combined-abundances.tree.txt --annotation Combined-abundances.annot.txt --most_abundant 100 --abundance_threshold 1 --least_biomarkers 3 --annotations 5,6 --external_annotations 7 --min_clade_size 1 --title "Healthy Tissue Microbiome"

# Create an xml file for graphlan based on the annotations and tree files
graphlan_annotate.py --annot Combined-abundances.annot.txt Combined-abundances.tree.txt Combined-abundances.xml

# Generate the circular phylogenetic tree visualization
graphlan.py --dpi 3cd 00 Combined-abundances.xml HealthySamples.png --external_legends