library(reactome.db)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(data.table)

args<-commandArgs(TRUE)
organism <- tolower(args[1])
output_filename <- args[2]

if(organism == 'human'){
    entrez.universe <- keys(org.Hs.eg.db, "ENTREZID")
} else if(organism == 'mouse'){
    entrez.universe <- keys(org.Mm.eg.db, "ENTREZID")
} else {
    message('Unsupported organism choice.')
    quit(status=1)
}

# Selecting reactome gene sets. This gives a table of ENTREZID to PATHID:
# > head(pathways)
#   ENTREZID        PATHID
# 1        1  R-HSA-109582
# 2        1  R-HSA-114608
# 3        1  R-HSA-168249
# 4        1  R-HSA-168256
# 5        1 R-HSA-6798695
pathways <- na.omit(select(reactome.db, keys=entrez.universe, c("PATHID"),
                           keytype = 'ENTREZID'))

# creates a list where the PATHID is the name and it points at an
# array of entrezIDs for that pathway
pathways <- split(pathways$ENTREZID, pathways$PATHID)

# Name of the pathways instead of the cryptic identifiers (e.g. R-HSA-109582)
# > head(pathway2name)
#           PATHID                                      PATHNAME
# 1: R-HSA-1059683         Homo sapiens: Interleukin-6 signaling
# 2:  R-HSA-109581                       Homo sapiens: Apoptosis
# 3:  R-HSA-109582                      Homo sapiens: Hemostasis
pathway2name <- as.data.table(
    na.omit(
        select(reactome.db, 
            names(pathways),
            c("PATHNAME"), 'PATHID')
    )
)

# Remove organism prefix from the pathway name. Takes 
# "Homo sapiens: Apoptosis" ---> "Apoptosis"
pathway2name[, PATHNAME := sub("^[^:]*: ", "", PATHNAME)]
pathway2name <- setNames(pathway2name$PATHNAME, pathway2name$PATHID)
pathway2name <- iconv(pathway2name, "latin1", "ASCII", sub="")

# create the text for a GMT-format file. See:
# https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29
pathway.lines <- sapply(names(pathways), function(p) {
    link <- p
    name <- paste0(p, "_", pathway2name[p])
    name <- gsub("[ ()/]+", "_", name)
    sprintf("%s\t%s\t%s", name, link, paste0(pathways[[p]], collapse="\t"))
})
write(pathway.lines, file=output_filename)

