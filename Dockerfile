FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y wget r-base snakemake samtools bowtie2 fastqc python3 python3-pandas python-is-python3 libcurl4-openssl-dev pandoc libssl-dev 

RUN Rscript -e "install.packages('BiocManager',repos='http://cran.r-project.org', ask=F); BiocManager::install(c('Biostrings', 'Gviz'), ask=F);" && \
    Rscript -e "install.packages('reshape2', repos='http://cran.r-project.org')" && \
    Rscript -e "install.packages('rmarkdown',repos='http://cran.r-project.org')" && \
    Rscript -e "install.packages('data.table',repos='http://cran.r-project.org')" && \
    Rscript -e "install.packages('ggplot2',repos='http://cran.r-project.org')" && \
    Rscript -e "install.packages('pander',repos='http://cran.r-project.org')" && \
    Rscript -e "install.packages('R.utils',repos='http://cloud.r-project.org')" && \
    Rscript -e "install.packages('gridExtra',repos='http://cran.r-project.org')" && \
    Rscript -e "install.packages('knitr',repos='https://cran.r-project.org')" && \
    Rscript -e "install.packages('tidyr', repos='http://cran.r-project.org')" && \
    Rscript -e "install.packages('cowplot',repos='http://cran.r-project.org')" && \
    Rscript -e "install.packages('bbmle',repos='http://cloud.r-project.org')" && \
    Rscript -e "install.packages('emdbook',repos='http://cran.r-project.org')" && \
    Rscript -e "install.packages('cowplot', repos='http://cran.r-project.org')" && \
    Rscript -e "install.packages('ggbeeswarm', repos='http://cran.r-project.org')" && \
    Rscript -e "install.packages('tinytex',repos='http://cran.r-project.org'); tinytex::install_tinytex()"


WORKDIR /docker

# Install Bismark -- no compilation needed. Requires perl, but this is already installed:
# perl -v  ## This is perl 5, version 30, subversion 0 (v5.30.0)
# Modified from https://hub.docker.com/r/aryeelab/bismark/dockerfile :

RUN wget https://github.com/FelixKrueger/Bismark/archive/0.23.0.tar.gz && \
    tar zxf 0.23.0.tar.gz && \
    cp -p Bismark-0.23.0/bismark* /usr/bin/

# Install BBmap -- no compilation needed. Requires java, but this is already installed:
# java --version  ## openjdk 11.0.19 2023-04-18

RUN wget https://downloads.sourceforge.net/project/bbmap/BBMap_39.01.tar.gz && \
    tar zxf BBMap_39.01.tar.gz  && \
    cp -r -p bbmap /usr/bin/

RUN rm -rf /docker

COPY refs/ /refs/
COPY --chmod=0755 scripts/bismark_odk /usr/bin/
COPY scripts/ /scripts/
COPY demo/ /demo/

## can three of the layers above be combined, using .dockerexclude if needed?
## COPY . ./ 

# random bits to be used by some UNIX functions
RUN openssl enc -aes-256-ctr -pass pass:7 -nosalt </dev/zero 2>/dev/null | head -c 10000000 > /scripts/rand_seed7.txt

## avoid dynamically generating indices if new refs are added?
## RUN /usr/bin/bismark_genome_preparation --bowtie2 --verbose /refs/ref_v5

WORKDIR /bss

## add -n at command line for dry run, should be okay of it's at end, after -s
ENTRYPOINT ["snakemake", "-p", "-c1", "-s", "/scripts/bss_snakefile.py"]





