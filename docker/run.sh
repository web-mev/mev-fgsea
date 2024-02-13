#!/bin/bash

# the dge results file:
DGE_FILE=$1
ORGANISM=$2
GENE_IDS=$3
GENE_SETS_DB=$4

declare -A HUMAN_SETS_MAP
HUMAN_SETS_MAP["MSigDB Hallmark"]="/opt/resources/h.all.v2023.2.Hs.symbols.gmt"
HUMAN_SETS_MAP["Reactome"]="/opt/resources/c2.cp.reactome.v2023.2.Hs.symbols.gmt"
HUMAN_SETS_MAP["WikiPathways"]="/opt/resources/c2.cp.wikipathways.v2023.2.Hs.symbols.gmt"
HUMAN_SETS_MAP["Ontology - Biological process"]="/opt/resources/c5.go.bp.v2023.2.Hs.symbols.gmt"
HUMAN_SETS_MAP["Ontology - Molecular function"]="/opt/resources/c5.go.mf.v2023.2.Hs.symbols.gmt"
HUMAN_SETS_MAP["Ontology - Cellular component"]="/opt/resources/c5.go.cc.v2023.2.Hs.symbols.gmt"
HUMAN_SETS_MAP["(Human only) Oncogenic signatures"]="/opt/resources/c6.all.v2023.2.Hs.symbols.gmt"
HUMAN_SETS_MAP["(Human only) ImmuneSigDB"]="/opt/resources/c7.immunesigdb.v2023.2.Hs.symbols.gmt"


declare -A MOUSE_SETS_MAP
MOUSE_SETS_MAP["MSigDB Hallmark"]="/opt/resources/mh.all.v2023.2.Mm.symbols.gmt"
MOUSE_SETS_MAP["Reactome"]="/opt/resources/m2.cp.reactome.v2023.2.Mm.symbols.gmt"
MOUSE_SETS_MAP["WikiPathways"]="/opt/resources/m2.cp.wikipathways.v2023.2.Mm.symbols.gmt"
MOUSE_SETS_MAP["Ontology - Biological process"]="/opt/resources/m5.go.bp.v2023.2.Mm.symbols.gmt"
MOUSE_SETS_MAP["Ontology - Molecular function"]="/opt/resources/m5.go.mf.v2023.2.Mm.symbols.gmt"
MOUSE_SETS_MAP["Ontology - Cellular component"]="/opt/resources/m5.go.cc.v2023.2.Mm.symbols.gmt"

if [ $ORGANISM = "Human" ]
then
    GENE_SET_FILE=${HUMAN_SETS_MAP[${GENE_SETS_DB}]}
    GENE_MAP_FILE=/opt/resources/human_genes.tsv
elif [ $ORGANISM = "Mouse" ]
then
    GENE_SET_FILE=${MOUSE_SETS_MAP[${GENE_SETS_DB}]}
    GENE_MAP_FILE=/opt/resources/mouse_genes.tsv
else
    echo "Not a valid organism choice." >&2
    exit 1;
fi

if [ ! -z "$GENE_SET_FILE" ]
then
    Rscript /usr/local/bin/run_fgsea.R \
        $DGE_FILE \
        $GENE_MAP_FILE \
        $GENE_SET_FILE \
        $GENE_IDS
else
    echo "The choice of \"$GENE_SETS_DB\" was not a valid choice for your organism ($ORGANISM). Note that some gene sets are organism specific."
    exit 1;
fi