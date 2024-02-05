#!/bin/bash

DEST_DIR=$1

# Create the gene mapping files:
Rscript /usr/local/bin/create_gene_mappings.R human $DEST_DIR/human_genes.tsv
Rscript /usr/local/bin/create_gene_mappings.R mouse $DEST_DIR/mouse_genes.tsv

# Create the pathway databases:
Rscript /usr/local/bin/create_pathway_databases.R human $DEST_DIR/human_pathways.gmt
Rscript /usr/local/bin/create_pathway_databases.R mouse $DEST_DIR/mouse_pathways.gmt