FROM ubuntu:22.04

RUN apt-get -y update && \
    apt-get -y dist-upgrade && \
    DEBIAN_FRONTEND=noninteractive apt-get -y --no-install-recommends --no-install-suggests install \
        build-essential \
        bzip2 \
        cmake \
        git \
        libnss-sss \
        libtbb2 \
        libtbb-dev \
        ncurses-dev \
        nodejs \
        tzdata \
        unzip \
        wget \
        zlib1g \
        zlib1g-dev \
        vim \
        htop && \
    apt-get clean


## R and latest R-studio

RUN apt-get install -y r-base && \
    apt-get clean

RUN apt-get install -y --no-install-recommends --no-install-suggests gdebi-core && \
    apt-get clean

RUN wget https://download2.rstudio.org/server/jammy/amd64/rstudio-server-2022.07.1-554-amd64.deb -O /tmp/rstudio.deb && \
    gdebi -n /tmp/rstudio.deb && \
    rm -rf /tmp/rstudio.deb && \
    apt-get clean

## Installing conda (taken from https://fabiorosado.dev/blog/install-conda-in-docker/)

ENV CONDA_DIR /opt/conda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
   /bin/bash ~/miniconda.sh -b -p /opt/conda

# Put conda in path so we can use conda activate
ENV PATH=$CONDA_DIR/bin:$PATH

RUN conda install -c conda-forge mamba

COPY rna-seq-r.yaml /opt/conda/environment.yaml
RUN mamba env create -f /opt/conda/environment.yaml

## Path to R from conda env
RUN echo "rsession-which-r=/opt/conda/envs/rna-seq-r/bin/R" | \
    tee /etc/rstudio/rserver.conf


## Adding student user

RUN useradd -G rstudio-server -m -d /home/student -s /bin/bash student && \
    echo student:sysbiopass | chpasswd

EXPOSE 8787

USER root

CMD /usr/lib/rstudio-server/bin/rserver --server-daemonize=0 --server-app-armor-enabled=0
