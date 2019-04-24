FROM ubuntu:16.04

USER root

# Update all the packages
RUN apt-get update --fix-missing

# Essential development tools and headers
RUN apt-get -y --no-install-recommends --no-install-suggests install \
		build-essential \
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
		wget \
		zlib1g \
		zlib1g-dev \
        htop && \
    apt-get autoremove -y && \
    apt-get clean

# R 3.5 & Rstudio
RUN apt-get install -y --no-install-recommends --no-install-suggests \
    gdebi-core software-properties-common && \
    add-apt-repository -y 'deb http://cloud.r-project.org/bin/linux/ubuntu xenial-cran35/' && \
    apt-get -y update && \
    apt-get install -y --allow-unauthenticated r-base r-base-dev && \
    curl --location http://download2.rstudio.org/rstudio-server-1.1.463-amd64.deb --output /tmp/rstudio.deb && \
    gdebi -n /tmp/rstudio.deb && \
    rm -rf /tmp/rstudio.deb && \
    apt-get autoremove -y && \
    apt-get clean

# R packages
ENV CRAN_URL https://cloud.r-project.org/
RUN R -e 'install.packages("BiocManager")' && \
    R -e 'BiocManager::install("ChIPseeker", version = "3.8")'

# Install conda
RUN curl --location https://repo.continuum.io/miniconda/Miniconda2-4.6.14-Linux-x86_64.sh --output ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh

ENV PATH /opt/conda/bin:${PATH}
ENV LD_LIBRARY_PATH /opt/conda/lib:${LD_LIBRARY_PATH}

# Conda envs
RUN mkdir -p /root/conda_envs
COPY conda_envs /root/conda_envs

# Bio utils, Python 2.7 with MACS2, SICER and snakemake
RUN conda env create -q --name bio --file /root/conda_envs/bio.environment.yml && \
    conda env create -q --name snakemake --file /root/conda_envs/snakemake.environment.yml && \
    conda clean -i -l -t -y && \
    rm -rf /usr/local/src/*

# Fix SICER paths
RUN sed -i 's#python#/opt/conda/envs/bio/bin/python#g' /opt/conda/envs/bio/share/sicer-1.1-3/SICER.sh && \
    ln -sf /opt/conda/envs/bio/share/sicer-1.1-3/SICER.sh /usr/local/bin/SICER.sh && \
    sed -i 's#python#/opt/conda/envs/bio/bin/python#g' /opt/conda/envs/bio/share/sicer-1.1-3/SICER-rb.sh && \
    ln -sf /opt/conda/envs/bio/share/sicer-1.1-3/SICER-rb.sh /usr/local/bin/SICER-rb.sh

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
    chown student:student /opt/conda

USER student

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
USER root
CMD /usr/lib/rstudio-server/bin/rserver --server-daemonize=0 --server-app-armor-enabled=0