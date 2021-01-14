From debian:stretch

SHELL ["/bin/bash", "-c"]

RUN apt-get update && \
  apt-get install -y build-essential \
    wget \
    zlib1g-dev \
    libcurl4-openssl-dev \
    libncurses-dev \
    libbz2-dev \
    liblzma-dev \
    libxml2-dev \
    libopenblas-dev

# Install all software under /opt/software:
RUN mkdir -p /opt/software 
RUN mkdir -p /opt/software/resources 

# Get the miniconda installer script and run it in silent mode:
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /opt/software/miniconda.sh && \
  bash /opt/software/miniconda.sh -b -p /opt/software/miniconda
ENV PATH="/opt/software/miniconda/bin:${PATH}"

# Install the R/bioc dependencies:
ADD environment.yml /opt/software/
RUN conda env create -f /opt/software/environment.yml

# add to the PATH. Otherwise, conda can swallow return codes
ENV PATH="/opt/software/miniconda/envs/fgsea_env/bin:${PATH}"

# Now add the necessary scripts:

# Allows us to create the required files prior to committing the repo
ADD create_required_files.sh /opt/software/
ADD create_gene_mappings.R /opt/software/
ADD create_pathway_databases.R /opt/software/
RUN chmod +x /opt/software/create_required_files.sh && \
  /opt/software/create_required_files.sh /opt/software/resources

# the scripts which actually run the analysis
ADD run.sh /opt/software
RUN chmod +x /opt/software/run.sh
ADD run_fgsea.R /opt/software/

ENTRYPOINT ["/bin/bash"]

