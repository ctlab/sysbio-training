# based on https://github.com/genome/docker-rnaseqget \

FROM ubuntu:16.04

RUN apt-get -y update && \
    apt-get -y dist-upgrade && \
    apt-get -y --no-install-recommends --no-install-suggests install \
		build-essential \
		bzip2 \
		cmake \
		git \
		libnss-sss \
		libtbb2 \
		libtbb-dev \
		ncurses-dev \
		nodejs \
		python-dev \
		python-pip \
		tzdata \
		unzip \
		wget \
		zlib1g \
		zlib1g-dev \
		vim \
        htop && \
    apt-get clean 

############################# R & Rstudio

RUN apt-get install -y --no-install-recommends --no-install-suggests \
        gdebi-core software-properties-common && \
    apt-get clean 

RUN add-apt-repository -y 'deb http://cloud.r-project.org/bin/linux/ubuntu xenial-cran35/' && \
    apt-get -y update

RUN apt-get install -y --allow-unauthenticated r-base r-base-dev && \
    apt-get clean

RUN wget http://download2.rstudio.org/rstudio-server-1.1.463-amd64.deb -O /tmp/rstudio.deb \
      && gdebi -n /tmp/rstudio.deb \
      && rm -rf /tmp/rstudio.deb \
      && apt-get clean

####### R packages

ENV CRAN_URL https://cloud.r-project.org/

RUN R -e 'install.packages("BiocManager")'

RUN apt-get install -y libxml2-dev curl libcurl4-openssl-dev \
      && apt-get clean
RUN R -e 'BiocManager::install(c("RCurl", "XML"))'

RUN apt-get install -y libssl-dev \
      && apt-get clean 

RUN R -e 'BiocManager::install(c("httr", "GenomicFeatures"))'

RUN R -e 'BiocManager::install(c(\
    "limma", "DESeq2", "fgsea", \
    "org.Mm.eg.db", "org.Hs.eg.db"))'


# Rmarkdown in rstudio
RUN R -e 'install.packages(c("rmarkdown", "rprojroot", "caTools"))'

# readr makes tximport faster
RUN R -e 'BiocManager::install(c("tximport", "readr"))'


RUN R -e 'install.packages(c("devtools"))'
RUN R -e 'BiocManager::install("fgsea")'
RUN R -e 'install.packages(c("ggrepel", "msigdbr", "pheatmap"))'
RUN R -e 'BiocManager::install(c("GEOquery", "rhdf5"))'
RUN R -e 'BiocManager::install("pachterlab/sleuth")'



RUN useradd -G rstudio-server -m -d /home/student -s /bin/bash student \
    && echo student:sysbiopass | chpasswd


############## Installing packages with conda


RUN mkdir /opt/conda
RUN chown student:student /opt/conda

RUN \
        echo 'PATH="/opt/conda/bin:${PATH}"' >> /etc/profile && \
        echo 'LD_LIBRARY_PATH="/opt/conda/lib:${LD_LIBRARY_PATH}"' >> /etc/profile

USER student

RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-4.5.11-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -f -p /opt/conda && \
    rm ~/miniconda.sh && \
    /opt/conda/bin/conda clean -tipsy


RUN \
        echo 'PATH="/opt/conda/bin:${PATH}"' >> ~/.profile && \
        echo 'LD_LIBRARY_PATH="/opt/conda/lib:${LD_LIBRARY_PATH}"' >> ~/.profile

ENV PATH="/opt/conda/bin:${PATH}"
ENV LD_LIBRARY_PATH="/opt/conda/lib:${LD_LIBRARY_PATH}"

RUN conda config --add channels conda-forge
RUN conda config --add channels bioconda
RUN conda install numpy scipy
RUN conda install scikit-learn
RUN conda install numba
RUN pip install umap-learn

# hot fixes

USER root
RUN R -e 'install.packages(c("ggplot2", "dplyr", "Seurat"))'
RUN R -e 'BiocManager::install("MAST")'
RUN R -e 'install.packages(c("shinyjs", "RColorBrewer", "scales"))'
RUN R -e 'devtools::install_github("AnalytixWare/ShinySky")'
RUN R -e 'install.packages(c("DT"))'


RUN conda install samtools 
RUN conda install deeptools


EXPOSE 8787

USER root

CMD /usr/lib/rstudio-server/bin/rserver --server-daemonize=0 --server-app-armor-enabled=0
