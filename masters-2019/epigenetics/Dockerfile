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
		libnss-sss \
		libtbb2 \
		libtbb-dev \
		mc \
		nano \
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

# R & Rstudio
RUN apt-get install -y --no-install-recommends --no-install-suggests \
        gdebi-core software-properties-common && \
    apt-get clean
RUN add-apt-repository -y 'deb http://cloud.r-project.org/bin/linux/ubuntu xenial-cran35/' && \
    apt-get -y update
RUN apt-get install -y --allow-unauthenticated r-base r-base-dev && \
    apt-get clean
RUN curl --location http://download2.rstudio.org/rstudio-server-1.1.463-amd64.deb --output /tmp/rstudio.deb \
      && gdebi -n /tmp/rstudio.deb \
      && rm -rf /tmp/rstudio.deb \
      && apt-get clean

# Install conda
RUN curl --location https://repo.continuum.io/miniconda/Miniconda3-4.5.11-Linux-x86_64.sh --output ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh
RUN ln -snf /bin/bash /bin/sh

# Conda envs
RUN mkdir -p /root/conda_envs
COPY conda_envs /root/conda_envs

# Bio utils, Python 2.7 with MACS2, SICER
RUN conda env create -q --name bio --file /root/conda_envs/bio.environment.yml
RUN sed -i 's#python#/opt/conda/envs/bio/bin/python#g' /opt/conda/envs/bio/share/sicer-1.1-3/SICER.sh
RUN ln -sf /opt/conda/envs/bio/share/sicer-1.1-3/SICER.sh /usr/local/bin/SICER.sh
RUN sed -i 's#python#/opt/conda/envs/bio/bin/python#g' /opt/conda/envs/bio/share/sicer-1.1-3/SICER-rb.sh
RUN ln -sf /opt/conda/envs/bio/share/sicer-1.1-3/SICER-rb.sh /usr/local/bin/SICER-rb.sh
ENV PATH /opt/conda/envs/bio/bin/:${PATH}
ENV LD_LIBRARY_PATH /opt/conda/lib:${LD_LIBRARY_PATH}

# Fix: R shell uses env variables from ~/profile
RUN \
	echo 'PATH="/opt/conda/envs/bio/bin:${PATH}"' >> ~/.profile && \
	echo 'LD_LIBRARY_PATH="/opt/conda/lib:${LD_LIBRARY_PATH}"' >> ~/.profile

ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8

RUN echo "export LC_ALL=C.UTF-8" >> /etc/profile && \
    echo "export LANG=C.UTF-8" >> /etc/profile

# Add user
RUN useradd -G rstudio-server -m -d /home/student -s /bin/bash student \
    && echo student:sysbiopass | chpasswd

RUN chown student:student /opt/conda

# Download Picard tools
RUN curl --location https://github.com/broadinstitute/picard/releases/download/2.10.7/picard.jar --output /opt/picard.jar

# Download SPAN Peak Analyzer, see details: https://research.jetbrains.org/groups/biolabs/tools/span-peak-analyzer
RUN curl --location https://download.jetbrains.com/biolabs/span/span-0.10.0.4787.jar --output /opt/span.jar

# expose default R-Studio port
EXPOSE 8787
USER root
CMD /usr/lib/rstudio-server/bin/rserver --server-daemonize=0 --server-app-armor-enabled=0