FROM mambaorg/micromamba

# for building the conda env
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
# Build the environment based on the yml file
RUN micromamba install -y -n base -f /tmp/environment.yml \
  && micromamba clean --all --yes

USER root
RUN mkdir -p /opt/resources && \
  chmod a+r /opt/resources && \
  chown $MAMBA_USER:$MAMBA_USER /opt/resources

# Allows us to create the required files prior to committing the repo
ADD create_gene_mappings.R /usr/local/bin/

USER $MAMBA_USER
ARG MAMBA_DOCKERFILE_ACTIVATE=1

# Download msigdb database files for human and mouse:

# Human:
RUN curl ftp://ftp.broadinstitute.org/pub/gsea/msigdb/human/gene_sets/h.all.v2023.2.Hs.symbols.gmt > /opt/resources/h.all.v2023.2.Hs.symbols.gmt
RUN curl ftp://ftp.broadinstitute.org/pub/gsea/msigdb/human/gene_sets/c2.cp.reactome.v2023.2.Hs.symbols.gmt > /opt/resources/c2.cp.reactome.v2023.2.Hs.symbols.gmt
RUN curl ftp://ftp.broadinstitute.org/pub/gsea/msigdb/human/gene_sets/c2.cp.wikipathways.v2023.2.Hs.symbols.gmt > /opt/resources/c2.cp.wikipathways.v2023.2.Hs.symbols.gmt
RUN curl ftp://ftp.broadinstitute.org/pub/gsea/msigdb/human/gene_sets/c5.go.bp.v2023.2.Hs.symbols.gmt > /opt/resources/c5.go.bp.v2023.2.Hs.symbols.gmt
RUN curl ftp://ftp.broadinstitute.org/pub/gsea/msigdb/human/gene_sets/c5.go.mf.v2023.2.Hs.symbols.gmt > /opt/resources/c5.go.mf.v2023.2.Hs.symbols.gmt
RUN curl ftp://ftp.broadinstitute.org/pub/gsea/msigdb/human/gene_sets/c5.go.cc.v2023.2.Hs.symbols.gmt > /opt/resources/c5.go.cc.v2023.2.Hs.symbols.gmt
RUN curl ftp://ftp.broadinstitute.org/pub/gsea/msigdb/human/gene_sets/c6.all.v2023.2.Hs.symbols.gmt > /opt/resources/c6.all.v2023.2.Hs.symbols.gmt
RUN curl ftp://ftp.broadinstitute.org/pub/gsea/msigdb/human/gene_sets/c7.immunesigdb.v2023.2.Hs.symbols.gmt > /opt/resources/c7.immunesigdb.v2023.2.Hs.symbols.gmt


# Mouse
RUN curl ftp://ftp.broadinstitute.org/pub/gsea/msigdb/mouse/gene_sets/mh.all.v2023.2.Mm.symbols.gmt > /opt/resources/mh.all.v2023.2.Mm.symbols.gmt
RUN curl ftp://ftp.broadinstitute.org/pub/gsea/msigdb/mouse/gene_sets/m2.cp.reactome.v2023.2.Mm.symbols.gmt > /opt/resources/m2.cp.reactome.v2023.2.Mm.symbols.gmt
RUN curl ftp://ftp.broadinstitute.org/pub/gsea/msigdb/mouse/gene_sets/m2.cp.wikipathways.v2023.2.Mm.symbols.gmt > /opt/resources/m2.cp.wikipathways.v2023.2.Mm.symbols.gmt
RUN curl ftp://ftp.broadinstitute.org/pub/gsea/msigdb/mouse/gene_sets/m5.go.bp.v2023.2.Mm.symbols.gmt > /opt/resources/m5.go.bp.v2023.2.Mm.symbols.gmt
RUN curl ftp://ftp.broadinstitute.org/pub/gsea/msigdb/mouse/gene_sets/m5.go.mf.v2023.2.Mm.symbols.gmt > /opt/resources/m5.go.mf.v2023.2.Mm.symbols.gmt
RUN curl ftp://ftp.broadinstitute.org/pub/gsea/msigdb/mouse/gene_sets/m5.go.cc.v2023.2.Mm.symbols.gmt > /opt/resources/m5.go.cc.v2023.2.Mm.symbols.gmt


# Create the gene mapping files:
RUN Rscript /usr/local/bin/create_gene_mappings.R human /opt/resources/human_genes.tsv
RUN Rscript /usr/local/bin/create_gene_mappings.R mouse /opt/resources/mouse_genes.tsv

# the scripts which actually run the analysis
ADD run.sh /usr/local/bin/
ADD run_fgsea.R /usr/local/bin/