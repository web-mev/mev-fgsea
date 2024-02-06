suppressMessages(suppressWarnings(library("fgsea", character.only=T, warn.conflicts = F, quietly = T)))
suppressMessages(suppressWarnings(library("dplyr", character.only=T, warn.conflicts = F, quietly = T)))
suppressMessages(suppressWarnings(library("rjson", character.only=T, warn.conflicts = F, quietly = T)))

# args from command line:
args<-commandArgs(TRUE)

# the results of a differential gene expression analysis.
# We expect a column giving the log fold change, named "logFC"
# and a p-value column named "pval".
# The first column has the gene identifiers/symbols
# The file SHOULD have a header line
DGE_RESULTS_FILE <- args[1]

# a table which maps gene symbols to entrez IDs
# Includes 'ALIAS', 'SYMBOL', 'ENSEMBL', 'ENTREZID'
GENE_MAPPING_FILE <- args[2]

# a GMT-format file with the pathways:
GMT_PATHWAYS_FILE <- args[3]

# the "type" of gene identifier
GENE_ID_TYPE <- tolower(args[4])

# Gene sets smaller than this will be removed. Defaults to 15
MIN_SETSIZE = 15

# change the working directory to co-locate with the counts file:
working_dir <- dirname(DGE_RESULTS_FILE)
setwd(working_dir)

# load the data:
dge_df = read.table(DGE_RESULTS_FILE, sep = '\t', stringsAsFactors=F, header=T, row.names=1)

# assert that the expected columns are there:
missingCol <- function(df, col){
    return(!(col %in% colnames(df)))
}
LOG_2_FC = 'log2FoldChange'
LOG_FC = 'logFC'
PVAL = 'pval'
PVALUE = 'pvalue'

if (missingCol(dge_df, LOG_2_FC) & missingCol(dge_df, LOG_FC)){
    message('We require a column named "log2FoldChange" or "logFC" in the input matrix.')
    quit(status=1)
}
if (missingCol(dge_df, PVAL) & missingCol(dge_df, PVALUE)){
    message('We require a column named "pval" or "pvalue" in the input matrix.')
    quit(status=1)
}

# figure out which column has the fold change
if (LOG_2_FC %in% colnames(dge_df)){
    fc_col_num = which(colnames(dge_df) == LOG_2_FC)
} else {
    fc_col_num = which(colnames(dge_df) == LOG_FC)
}

# figure out which column has the p-value
if (PVAL %in% colnames(dge_df)){
    pval_col_num = which(colnames(dge_df) == PVAL)
} else {
    pval_col_num = which(colnames(dge_df) == PVALUE)
}


# make a symbol column from the row index:
dge_df$symbol = rownames(dge_df)

# drop any NA's
dge_df = dge_df[complete.cases(dge_df[,c(fc_col_num, pval_col_num)]),]

# calculate the ranking value-- this is the -log10(p-value)*sgn(lfc)
dge_df$rnk = -log10(dge_df[,pval_col_num])*sign(dge_df[,fc_col_num])

# subset the dataframe:
dge_df = dge_df[c('symbol','rnk')]

# check for infinite values originating from the log function. Some DGE programs (DESeq2)
# can assign hard zeros to the p-value. Handle this by finding the largest finite values
# and replacing the +/- Inf with that max PLUS a constant value. The 'rnk' column is 
# simply used to rank the hits, so the value shouldn't matter as long as it preserves the
# approprite order. We could be fancier and assign all the infinite values different ranks,
# but we just assign them all the same here.
finite_idx = is.finite(dge_df$rnk)
min_finite = min(dge_df[finite_idx, 'rnk'])
max_finite = max(dge_df[finite_idx, 'rnk'])
neg_inf_idx = !finite_idx & (dge_df$rnk < 0)
pos_inf_idx = !finite_idx & (dge_df$rnk > 0)
delta = 5 # the offset to add/subtract
dge_df[neg_inf_idx, 'rnk'] = min_finite - delta
dge_df[pos_inf_idx, 'rnk'] = max_finite + delta

# read the dataframe which contains the gene mapping info:
gene_info_df = read.table(GENE_MAPPING_FILE, sep='\t')

# merge to keep only those where we have a mapping:
if (GENE_ID_TYPE != 'symbol') {
dge_df = merge(
    dge_df, gene_info_df[,c(toupper(GENE_ID_TYPE), 'SYMBOL')], by.x=0, by.y=toupper(GENE_ID_TYPE))
} else {
    dge_df['SYMBOL'] = rownames(dge_df)
}

# If no remaining rows, error out
if(dim(dge_df)[1] == 0){
    message('After mapping the gene identifiers, there were no remaining rows. Was the choice of gene identifier correct?')
    quit(status=1)
}

dge_df = dge_df[,c('symbol', 'rnk', 'SYMBOL')]

# drop those without an ID
dge_df <- dge_df[!unlist(lapply(dge_df$SYMBOL, is.null)),]

# remove duplicates:
dge_df <- distinct(dge_df, SYMBOL, .keep_all=T)

# list of the stats named by the entrez ID
# e.g.:
# > head(stats)
#      TSPAN6        TNMD        DPM1       SCYL3       FIRRM         FGR 
# -31.7067988  -0.3059057 -91.4904801  33.8878384 312.2079695   0.2784677 
stats = setNames(dge_df[,'rnk'], dge_df[,'SYMBOL'])

pathways = gmtPathways(GMT_PATHWAYS_FILE)
fgseaRes <- fgsea(pathways = pathways, 
                  stats    = stats,
                  minSize  = MIN_SETSIZE,
                  maxSize  = 500)

if (dim(fgseaRes)[1] == 0){
    message(sprintf('The table of gsea results was empty. Most commonly, this error can also be caused by choosing the incorrect gene identifiers or organism. However, this may also be due to small gene sets not passing the minimum size of %s genes.', MIN_SETSIZE))
    quit(status=1)
}

# for the front-end, we would like to display the GSEA "rug plot" where we show
# the location of a single pathway's genes in the ranked list.
# rnk gives the order, with the negative sign putting the largest positive 
# stats at the top of the list.
# For example:
# > dummyRanks = c("A"=5.5, "B"=-0.2, "C"=3.3, "D"=-3.6, "E"=1.1)
# > dummyRanks
#    A    B    C    D    E 
#  5.5 -0.2  3.3 -3.6  1.1 
# > rnk = rank(-dummyRanks)
# > rnk
# A B C D E 
# 1 4 2 5 3
# Then, given a pathway (e.g. with genes A and D), we get 
rnk <- rank(-stats)

# runs through all the pathways and gets the rank of that pathway's genes
# in the whole ranked list.
# For example,
# > pathway_ranks
# $`186574_Endocrine-committed_Ngn3+_progenitor_cells`
# [1] 3887 4228
#
# $`1368092_Rora_activates_gene_expression`
# [1] 9311 9005 1805 7229 8473
#
# For the first pathway, it has genes with entrezID of 
# [1] "18012" "18088" "18506" "53626"
# Only two of those were found in the ranked list based on the DGE results.
# Those genes (18088, 18506) were found at positions/ranks 3887 and 4228
pathway_ranks <- do.call(rbind, lapply(pathways, function(p) {
    list(ranks=sort(unname(as.vector(na.omit(rnk[p])))))
}))

# merge the fgsea results with the ranks
m = merge(as.data.frame(fgseaRes), pathway_ranks, by.x='pathway', by.y=0)

# the ranks and leadingEdge columns are exported strangely to JSON
# if we use `m` directly. Create a list structure to use for export.

# want to re-map back to the original identifier, as the leadingEdge is 
# given as symbols.
remapping_df <- unique(gene_info_df[,c(toupper(GENE_ID_TYPE), 'SYMBOL')])
remapping_df <- distinct(remapping_df, SYMBOL, .keep_all=T)
rownames(remapping_df) <- remapping_df$SYMBOL

q = apply(
    m,
    1,
    function(r){
        leadingEdge <- remapping_df[r$leadingEdge, toupper(GENE_ID_TYPE)]
        if (length(leadingEdge) == 1){
            leadingEdge = list(leadingEdge)
        }
        list(
            pathway=r$pathway, 
            pval=r$pval, 
            padj=r$padj, 
            log2err=r$log2err,
            ES=r$ES,
            NES=r$NES,
            size=r$size,
            ranks=r$ranks, 
            leadingEdge=leadingEdge
        )
    }
)

# reformat the table as a json string and write to file
results_json_str <- toJSON(q)
output_filename = 'fgsea_results_demo.json'
results_json_file <- paste(working_dir, output_filename, sep='/')
write(results_json_str, results_json_file)

# for WebMEV compatability, need to create an outputs.json file.
json_str = paste0(
       '{"pathway_results":"', results_json_file, '"}'
)
output_json <- paste(working_dir, 'outputs.json', sep='/')
write(json_str, output_json)
