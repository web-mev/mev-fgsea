#!/bin/bash

# the dge results file:
DGE_FILE=$1
ORGANISM=$2
GENE_IDS=$3

if [ $ORGANISM = "Human" ]
then
    Rscript /opt/software/run_fgsea.R \
        $DGE_FILE \
        /opt/software/resources/human_genes.tsv \
        /opt/software/resources/human_pathways.gmt \
        $GENE_IDS
elif [ $ORGANISM = "Mouse" ]
then
    Rscript /opt/software/run_fgsea.R \
        $DGE_FILE \
        /opt/software/resources/mouse_genes.tsv \
        /opt/software/resources/mouse_pathways.gmt \
        $GENE_IDS
else
    echo "Not a valid organism choice." >&2
    exit 1;
fi