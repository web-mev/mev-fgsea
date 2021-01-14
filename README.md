## About

This is an WebMEV-compatible process for running fast gene set enrichment analysis (FGSEA) https://www.biorxiv.org/content/10.1101/060012v2.full.pdf

**Generating the required additional files**

To perform gene identifier mapping and associate those genes with biological pathways, we have a couple of "helper" scripts that are distributed with the repository. These are not run in "real time" (i.e. when the fgsea process is invoked by a user), but rather should be run prior to committing the repository. This will create the necessary "database-like" files which will then be committed to the repository. This allows full transparency about the gene mapping and pathways used to run the analysis.

To create those files, change into the `docker` directory and run the following:

```
# Build the docker image:
docker build -t <user>/<image> .

# execute the script which creates the necessary files. By mounting
# the current directory (and using the environment variable), those
# files are written into the current docker/ directory on the host machine
docker run -d \
    -v $PWD:/workspace \
    <user>/<image>
    /opt/software/create_required_files.sh /workspace/
```

**Other**

To export the fully-configured environment, we used `conda env export -n fgsea_env > environment.yml`. However, to create that initial environment, we ran:
```
conda install \
    -c conda-forge -c r -c bioconda \
    r-base=4.0.2 r-dplyr bioconductor-fgsea bioconductor-org.hs.eg.db bioconductor-org.mm.eg.db r-rjson bioconductor-reactome.db
```
This latter command is given here so we can build new environments at a later time. Even seemingly small package additions changes can't always be integrated without having conda resolve the full dependency graph. 