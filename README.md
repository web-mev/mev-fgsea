## About

This is an WebMEV-compatible process for running fast gene set enrichment analysis (FGSEA) https://www.biorxiv.org/content/10.1101/060012v2.full.pdf

**Generating the required additional files**

To perform gene identifier mapping and associate those genes with biological pathways, we have a "helper" script (`docker/create_gene_mappings.R`) that is distributed with the repository. This script not run in "real time" (i.e. when the fgsea process is invoked by a user), but rather is run during the container build process and is thus packaged with the Docker image. This allows full transparency about the gene mapping and pathways used to run the analysis.

Also note that the pathway database files are pulled from MSigDB on container build and are also packaged with the image. 

Both the gene mappinga and pathway database files are located in `/opt/resources/`

---

### To run external of WebMeV:
 Either:
- build the Docker image using the contents of the `docker/` folder (e.g. `docker build -t myuser/fgsea:v1 .`) 
- pull the docker image from the GitHub container repository (see [https://github.com/web-mev/mev-fgsea/pkgs/container/mev-edger](https://github.com/web-mev/mev-fgsea/pkgs/container/mev-fgsea))

To run, change into the directory containing the differential expression results you wish to run pathway analysis on. Then:
```
docker run -it -v $PWD:/work <IMAGE> Rscript /usr/local/bin/run_fgsea.R \
        /work/<DiffExp_FILE> \
        <GENE ID MAPPING FILE> \
        <GENE SET FILE IN GMT FORMAT> \
        <GENE IDENTIFIER>
```
Note that we mounted your current directory to `/work` inside the Docker container. Hence, your file of differential expression results is relative to `/work`.

This assumes:

- Your differential expression results (in tab-delimited format) have the gene identifiers in the first column (either symbol, Ensembl, or RefSeq) and that the file contains the following columns:
  - A column named "log2FoldChange" or "logFC". This has the estimated fold change for the contrast between sample cohorts
  - A column named "pval" or "pvalue". This has the significance of the differential expression for the contrast. Can be the raw p-value or adjusted p-value.
  
  These columns are used to created a ranked list where the largest/most significant upregulated genes are at the top of the list and the largest/most significant downregulated genes are at the bottom.

- `<GENE ID MAPPING FILE>` is a file that allows us to convert between gene identifiers. It can be:
  - Human: `/opt/resources/human_genes.tsv`
  - Mouse: `/opt/resources/mouse_genes.tsv`
- `<GENE SET FILE IN GMT FORMAT>` is a GMT-format file which provides the pathways (e.g. name of the pathway and the associated genes). See the `Dockerfile` for the pathway files that are built into the image (    in the `/opt/resources/` directory. You can of course also provide your own GMT-format file (https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#Gene_Set_Database_Formats). Note that the main R script assumes the GMT file contains pathways using gene symbols as the identifiers, so your own custom GMT-format file would also need to specify the pathways using gene symbols.
-  `<GENE IDENTIFIER>` tells the script which gene identifiers are in your differential expression results (the identifiers in the first column). Currently this must be one of:
   - `SYMBOL`
   - `ENSEMBL`
   - `REFSEQ`
 
Note that the script needs to map your differential expression results to common gene symbols. Hence, there may be situations where mapping can be ambiguous due to the nature of mapping between multiple gene identifier conventions which are not exactly 1:1. 
