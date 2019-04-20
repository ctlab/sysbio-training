FROM ubuntu:latest

USER root

# Update all the packages
RUN apt-get update --fix-missing

# Install conda
RUN apt-get install -y curl bzip2
RUN curl --location https://repo.continuum.io/miniconda/Miniconda3-4.5.11-Linux-x86_64.sh --output ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh
ENV PATH /opt/conda/bin:${PATH}
ENV LD_LIBRARY_PATH /opt/conda/lib:${LD_LIBRARY_PATH}
RUN ln -snf /bin/bash /bin/sh

# Conda envs
RUN mkdir -p /root/conda_envs
COPY conda_envs /root/conda_envs

# Python 2.7 with MACS2, SICER
RUN conda env create -q --name py2.7 --file /root/conda_envs/py2.7.environment.yml
RUN ln -sf /opt/conda/envs/py2.7/bin/macs2 /usr/local/bin/macs2
# SICER is python2 library, force it!
RUN sed -i 's#python#/opt/conda/envs/py2.7/bin/python#g' /opt/conda/envs/py2.7/share/sicer-1.1-3/SICER.sh
RUN ln -sf /opt/conda/envs/py2.7/share/sicer-1.1-3/SICER.sh /usr/local/bin/SICER.sh
RUN ln -sf /opt/conda/envs/py2.7/share/sicer-1.1-3/SICER-rb.sh /usr/local/bin/SICER-rb.sh

# Bio utils available in path
RUN conda env create -q --name bio --file /root/conda_envs/bio.environment.yml
ENV PATH /opt/conda/envs/bio/bin/:$PATH

# R & Rstudio
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

# Add user
RUN useradd -G rstudio-server -m -d /home/student -s /bin/bash student \
    && echo student:sysbiopass | chpasswd
RUN mkdir /opt/conda
RUN chown student:student /opt/conda

USER student

# Download Picard tools
RUN cd ~ && curl --location https://github.com/broadinstitute/picard/releases/download/2.10.7/picard.jar \
    --output picard.jar

# Download SPAN Peak Analyzer, see details: https://research.jetbrains.org/groups/biolabs/tools/span-peak-analyzer
RUN cd ~ && curl --location https://download.jetbrains.com/biolabs/span/span-0.10.0.4787.jar --output span.jar

# hot fixes end
EXPOSE 8787
USER root
CMD /usr/lib/rstudio-server/bin/rserver --server-daemonize=0 --server-app-armor-enabled=0