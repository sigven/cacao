FROM ubuntu:18.04

RUN apt-get update \
  && apt-get install -y python3-pip python3-dev \
  && cd /usr/local/bin \
  && ln -s /usr/bin/python3 python \
  && pip3 install --upgrade pip

RUN apt-get update && apt-get -y install apache2 apt-utils build-essential curl git unzip vim wget sudo libudunits2-dev libgeos-dev libgdal-dev
RUN apt-get update && apt-get install apt-transport-https


RUN sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
RUN sudo apt-get update \
   && sudo apt-get -y install software-properties-common

ENV TZ=Europe/Minsk
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

RUN sudo add-apt-repository 'deb [arch=amd64,i386] https://cloud.r-project.org/bin/linux/ubuntu bionic-cran40/'
RUN apt-get update && apt-get -y install r-base

USER root
WORKDIR /

ENV PACKAGE_BIO="libhts2 bedtools"
ENV PACKAGE_DEV="gfortran gcc-multilib autoconf libmariadb-client-lgpl-dev liblzma-dev libncurses5-dev libblas-dev liblapack-dev libssh2-1-dev libxml2-dev vim libssl-dev libcairo2-dev libbz2-dev libcurl4-openssl-dev"
ENV PYTHON_MODULES="numpy==1.19.2 cython==0.29.21 scipy==1.5.3 pandas==1.1.3 cyvcf2==0.20.9 toml==0.10.1"

RUN apt-get update \
	&& apt-get install -y --no-install-recommends \
		nano ed locales vim-tiny fonts-texgyre \
    $PACKAGE_DEV $PACKAGE_BIO \
	&& rm -rf /var/lib/apt/lists/*

RUN apt-get autoremove

RUN wget https://github.com/samtools/samtools/releases/download/1.10/samtools-1.10.tar.bz2
RUN bunzip2 -dc samtools-1.10.tar.bz2 | tar xvf -
RUN cd samtools-1.10 && ./configure --prefix=/usr/local/bin && make && make install

## ADD pcgrr source package
WORKDIR /
ADD R/ /
WORKDIR /

# RUN apt-get update && apt-get install -y --no-install-recommends \
#   r-base-core \
#   r-recommended \
#  	r-base

## cacao R package dependencies: restore all R package versions and dependencies with renv package
## PINNED versions of all packages can be found in /renv.lock (JSON format)
RUN R -e "install.packages(repos = NULL, 'https://cloud.r-project.org/src/contrib/renv_0.12.3.tar.gz', lib = '/usr/lib/R/library/')"
RUN cp /cacao/renv.lock /
RUN R -e "library(renv); renv::restore(prompt = F, library = '/usr/lib/R/library/')"

# ENV PATH=$PATH:/pcgr
# ENV PYTHONPATH=:/pcgr/lib:${PYTHONPATH}
# WORKDIR /

## PINNED version - cacao
RUN R -e "install.packages(repos = NULL, 'cacao_0.3.1.tar.gz', lib = '/usr/lib/R/library/')"


# RUN R -e "install.packages(c('configr','rmarkdown','dplyr','stringr','tidyr','magrittr','devtools','htmltools','plotly','ggplot2'), dependencies = T, repos = 'http://cran.us.r-project.org')"
# RUN R -e "library(devtools); devtools::install_github('rstudio/DT'); devtools::install_github('mjkallen/rlogging'); devtools::install_github('kent37/summarywidget')"
# RUN R -e "library(devtools); devtools::install_github('rstudio/crosstalk')"
# RUN rm -rf /tmp/R*

## Install tools used for compilation
RUN sudo -H pip install --upgrade pip
RUN sudo -H pip install -U setuptools
RUN sudo -H pip install $PYTHON_MODULES


ENV PATH=/usr/local/bin/bin:$PATH

## Mosdepth installation from Bioconda
WORKDIR /
ENV HOME=/usr/local
RUN wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
RUN bash Miniconda3-latest-Linux-x86_64.sh -b
ENV PATH=$PATH:/usr/local/miniconda3/bin:
RUN conda config --add channels bioconda
RUN conda install mosdepth

## PINNED VERSION: pandoc (for HTML report generation)
RUN wget https://github.com/jgm/pandoc/releases/download/2.11.2/pandoc-2.11.2-1-amd64.deb && \
  dpkg -i pandoc* && \
  rm pandoc* && \
  apt-get clean

## Clean Up
RUN apt-get clean autoclean
RUN rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
RUN rm -rf /var/lib/{dpkg,cache,log}
RUN rm -rf /samtools-1.10.tar.bz2
RUN rm -rf /Miniconda3-latest-Linux-x86_64.sh

WORKDIR /
ADD cacao.R /
ADD cacao.py /

# WORKDIR /
# ADD R /
# RUN R -e "devtools::install('cacao')"

WORKDIR /workdir/output
