# Dockerfile for gene expression module

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


# Rmarkdown in rstudio
RUN R -e 'install.packages(c("rmarkdown", "rprojroot", "caTools"))'

RUN R -e 'install.packages(c("devtools"))'


RUN useradd -G rstudio-server -m -d /home/student -s /bin/bash student \
    && echo student:sysbiopass | chpasswd


############## Installing packages with conda


RUN mkdir /opt/conda
RUN chown student:student /opt/conda

RUN \
	echo 'PATH="/opt/conda/bin:${PATH}"' >> /etc/profile && \
	echo 'LD_LIBRARY_PATH="/opt/conda/lib:${LD_LIBRARY_PATH}"' >> /etc/profile

USER student

RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -f -p /opt/conda && \
    rm ~/miniconda.sh && \
    /opt/conda/bin/conda clean -tipsy 


RUN \
	echo 'PATH="/opt/conda/bin:${PATH}"' >> ~/.profile && \
	echo 'LD_LIBRARY_PATH="/opt/conda/lib:${LD_LIBRARY_PATH}"' >> ~/.profile

ENV PATH="/opt/conda/bin:${PATH}"
ENV LD_LIBRARY_PATH="/opt/conda/lib:${LD_LIBRARY_PATH}"

RUN conda config --add channels bioconda
RUN conda install picard samtools bwa

RUN conda install star rsem kallisto
RUN conda install trimmomatic trinity spades transrate-tools cd-hit hisat2 subread fastqc deeptools
RUN conda install -c biocore blast-plus

# hot fixes here
USER root

RUN pip install RSeQC
RUN pip install multiqc

RUN cd /home/student && \
  git clone https://gitlab.com/ezlab/busco.git && \
  cd busco && \
  python setup.py install

# Using R from rstudio, not conda's
RUN conda remove --force r-base 

RUN apt-get -y install less

RUN echo export LC_ALL=C.UTF-8 >> /etc/profile && \
    echo export LANG=C.UTF-8 >> /etc/profile

ENV LC_ALL=C.UTF-8 
ENV LANG=C.UTF-8 

RUN R -e 'install.packages(c("BiocManager", "ggplot2", "reshape2", "MASS", "RColorBrewer", "pheatmap", "ggrepel", "dplyr", "Seurat"))'
RUN R -e 'BiocManager::install(c("affxparser", "affy", "GEOquery", "mouse4302.db", "Biobase", "limma", "sva", "fgsea", "DESeq2", "apeglm", "org.Mm.eg.db", "vsn", "MAST"))'
RUN R -e 'devtools::install_github("https://github.com/ctlab/LinSeed")'

RUN pip install BioPython

# hot fixes end


EXPOSE 8787

USER root

CMD /usr/lib/rstudio-server/bin/rserver --server-daemonize=0 --server-app-armor-enabled=0
