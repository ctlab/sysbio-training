FROM ubuntu:18.04

USER root

# Update all the packages
RUN apt-get update --fix-missing

# Essential development tools and headers
ENV TZ Europe/Moscow
RUN DEBIAN_FRONTEND="noninteractive" apt-get -y --no-install-recommends --no-install-suggests install \
		build-essential \
		bc \
		bzip2 \
		cmake \
		curl \
		git \
		less \
		libcurl4-openssl-dev \
		libnss-sss \
		libssl-dev \
		libxml2-dev \
		mc \
		tree \
		tzdata \
		unzip \
		ucommon-utils \
		util-linux \
		bsdmainutils \
		tree \
		vim \
		nano \
		wget \
		gnupg2 \
		zlib1g \
		zlib1g-dev \
        htop && \
    apt-get autoremove -y && \
    apt-get clean

# R 4.0 & Rstudio
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
RUN DEBIAN_FRONTEND="noninteractive" apt-get install -y --no-install-recommends --no-install-suggests \
    gdebi-core software-properties-common && \
    add-apt-repository -y 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran40/' && \
    apt-get -y update && \
    apt-get install -y --allow-unauthenticated r-base r-base-dev && \
    curl --location https://download2.rstudio.org/server/bionic/amd64/rstudio-server-1.3.1093-amd64.deb \
        --output /tmp/rstudio.deb && \
    gdebi -n /tmp/rstudio.deb && \
    rm -rf /tmp/rstudio.deb && \
    apt-get autoremove -y && \
    apt-get clean

# Install conda
RUN curl --location https://repo.continuum.io/miniconda/Miniconda2-4.6.14-Linux-x86_64.sh --output ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh

ENV PATH /opt/conda/bin:${PATH}
ENV LD_LIBRARY_PATH /opt/conda/lib:${LD_LIBRARY_PATH}

# Conda envs
RUN mkdir -p /root/conda_envs
COPY conda_envs /root/conda_envs

# Bio utils, Python 2.7 with MACS2, SICER
RUN conda env create -q --name bio --file /root/conda_envs/bio.environment.yml
RUN conda env create -q --name snakemake --file /root/conda_envs/snakemake.environment.yml
RUN conda clean -i -l -t -y && rm -rf /usr/local/src/*

# Fix SICER paths
RUN sed -i 's#python#/opt/conda/envs/bio/bin/python#g' /opt/conda/envs/bio/share/sicer-1.1-4/SICER.sh && \
    ln -sf /opt/conda/envs/bio/share/sicer-1.1-4/SICER.sh /usr/local/bin/SICER.sh && \
    sed -i 's#python#/opt/conda/envs/bio/bin/python#g' /opt/conda/envs/bio/share/sicer-1.1-4/SICER-rb.sh && \
    ln -sf /opt/conda/envs/bio/share/sicer-1.1-4/SICER-rb.sh /usr/local/bin/SICER-rb.sh

# Add bio to path
ENV PATH /opt/conda/envs/bio/bin/:${PATH}
ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8

# Fix: R shell uses env variables from ~/profile
RUN echo 'PATH="/opt/conda/bin:${PATH}"' >> /etc/profile && \
	echo 'PATH="/opt/conda/envs/bio/bin:${PATH}"' >> /etc/profile && \
	echo 'LD_LIBRARY_PATH="/opt/conda/lib:${LD_LIBRARY_PATH}"' >> /etc/profile && \
	echo "export LC_ALL=C.UTF-8" >> /etc/profile && \
    echo "export LANG=C.UTF-8" >> /etc/profile

# Add user and configure /opt/conda
RUN useradd -G rstudio-server -m -d /home/student -s /bin/bash student && \
    echo student:sysbiopass | chpasswd && \
    usermod -aG sudo student && \
    chown student:student /opt/conda


USER student

# Fix possible problems with missing LD_LIBRARY_PATH in conda
ENV PATH /opt/conda/bin:${PATH}
ENV PATH /opt/conda/envs/bio/bin/:${PATH}
ENV LD_LIBRARY_PATH /opt/conda/lib:${LD_LIBRARY_PATH}
ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8

# Fix: R shell uses env variables from ~/profile
RUN echo 'PATH="/opt/conda/bin:${PATH}"' >> ~/.profile && \
	echo 'PATH="/opt/conda/envs/bio/bin:${PATH}"' >> ~/.profile && \
	echo 'LD_LIBRARY_PATH="/opt/conda/lib:${LD_LIBRARY_PATH}"' >> ~/.profile && \
	echo "export LC_ALL=C.UTF-8" >> ~/.profile && \
    echo "export LANG=C.UTF-8" >> ~/.profile

# expose default R-Studio port
EXPOSE 8787

# Finish installation
USER root

# TODO: workaround: conda env bio installs R 3.6.1, which overrides R 4.0.3 installed by apt-get
RUN mv /opt/conda/envs/bio/bin/R /opt/conda/envs/bio/bin/R1 && \
    mv /opt/conda/envs/bio/bin/Rscript /opt/conda/envs/bio/bin/Rscript1

## Download all the data
#RUN wget https://artyomovlab.wustl.edu/publications/supp_materials/4Oleg/2020_ITMO_epigenetics_practice/mnt.tar.gz && \
# tar xvf mnt.tar.gz && rm mnt.tar.gz
#
#RUN cd /mnt && \
#    wget https://github.com/broadinstitute/picard/releases/download/2.23.8/picard.jar && \
#    wget https://download.jetbrains.com/biolabs/span/span-0.13.5244.jar && \
#    wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/hg19.chrom.sizes && \
#    wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_30/GRCh37_mapping/gencode.v30lift37.annotation.gtf.gz && \
#    wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/hg19.fa.gz
#

RUN mkdir -p /usr/local/lib/R/site-library && chmod a+rw /usr/local/lib/R/site-library
USER student

# R packages
ENV CRAN_URL https://cloud.r-project.org/
RUN R -e 'install.packages("BiocManager")' && \
    R -e 'BiocManager::install(c("ChIPseeker", "ChIPpeakAnno", "EnsDb.Hsapiens.v75", "ReactomePA", "UpSetR"), version = "3.12")' && \
    R -e 'install.packages("Matrix", repos = "http://cran.r-project.org")' && \
    R -e 'install.packages("boot", repos = "http://cran.r-project.org")' && \
    R -e 'install.packages("MASS", repos = "http://cran.r-project.org")' && \
    R -e 'install.packages("KernSmooth", repos = "http://cran.r-project.org")'

# Install Homer genome package
USER root
RUN chmod -R a+rw /opt/conda/envs/bio/share/homer
USER student
RUN cd ~ && perl /opt/conda/envs/bio/share/homer/configureHomer.pl -keepScript -install hg19

# Install missing package
RUN R -e 'BiocManager::install("org.Hs.eg.db")'

# Launch as root to for correct authentification
USER root
CMD /usr/lib/rstudio-server/bin/rserver --server-daemonize=0 --server-app-armor-enabled=0