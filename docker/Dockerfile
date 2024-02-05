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
ADD create_required_files.sh /usr/local/bin/
ADD create_gene_mappings.R /usr/local/bin/
ADD create_pathway_databases.R /usr/local/bin/

USER $MAMBA_USER
ARG MAMBA_DOCKERFILE_ACTIVATE=1
RUN /usr/local/bin/create_required_files.sh /opt/resources

# the scripts which actually run the analysis
ADD run.sh /usr/local/bin
ADD run_fgsea.R /usr/local/bin/