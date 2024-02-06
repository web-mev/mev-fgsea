suppressMessages(suppressWarnings(library("org.Hs.eg.db", character.only=T, warn.conflicts = F, quietly = T)))
suppressMessages(suppressWarnings(library("org.Mm.eg.db", character.only=T, warn.conflicts = F, quietly = T)))

args<-commandArgs(TRUE)
organism <- args[1]
output_filename <- args[2]

if(organism == 'human'){
    db <- org.Hs.eg.db
    keys <- names(as.list(org.Hs.egALIAS2EG))
} else if(organism == 'mouse'){
    db <- org.Mm.eg.db
    keys <- names(as.list(org.Mm.egALIAS2EG))
} else {
    message('Unsupported organism choice.')
    quit(status=1)
}
# Create a mapping table. This will allow us to take the gene identifiers
# and map them to symbols for the fgsea process.
# "ALIAS"	"ENSEMBL"	"SYMBOL"	"REFSEQ"
# "1"	"A1B"	"ENSG00000121410"	"A1BG"	"NM_130786"
# "2"	"A1B"	"ENSG00000121410"	"A1BG"	"NP_570602"
# "3"	"A1B"	"ENSG00000172164"	"SNTB1"	"NM_021021"
# "4"	"A1B"	"ENSG00000172164"	"SNTB1"	"NP_066301"
# "5"	"A1B"	"ENSG00000172164"	"SNTB1"	"XM_011517239"
gene_mappings <- select(
    db, 
    key=keys, 
    columns=c('ALIAS', 'ENSEMBL', 'SYMBOL', 'REFSEQ'), 
    keytype="ALIAS"
)
write.table(gene_mappings, output_filename, sep='\t', quote=T)
