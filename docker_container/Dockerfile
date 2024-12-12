FROM python:3.7-slim-buster

MAINTAINER Peiming (Peter) Huang, phuang@bcm.edu

ARG SCANDIUM_VER="1.5.1"
ARG HTSLIB_VER="1.19.1"
ARG SAMTOOLS_VER="1.19.2"

#Metadata
LABEL container.base.image="python:3.7-slim-buster"
LABEL software.name="Scandium with HTSlib/Samtools and python wrapper script to handle integration with AWS Batch and S3"
LABEL software.version=${SCANDIUM_VER}
LABEL software.description="An NGS sequencing QC software used internally at the HGSC"
LABEL software.website="NA"
LABEL tags="Genomics"

WORKDIR /scratch

# Install HTSlib/Samtools and necessary packages/dependencies
RUN apt-get -y update && apt-get -y install --no-install-recommends \
  autoconf bzip2 gcc libbz2-dev liblzma-dev libncurses5-dev libncursesw5-dev make wget zlib1g-dev && \
  rm -rf /var/lib/apt/lists/* && \
  wget https://github.com/samtools/htslib/releases/download/${HTSLIB_VER}/htslib-${HTSLIB_VER}.tar.bz2 && \
  tar -xvf htslib-${HTSLIB_VER}.tar.bz2 && \
  cd htslib-${HTSLIB_VER} && autoconf && ./configure && make && make install && \
  mv /scratch/htslib-${HTSLIB_VER} /usr/local/ && \
  cd /scratch && \
  wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VER}/samtools-${SAMTOOLS_VER}.tar.bz2 && \
  tar --bzip2 -xf samtools-${SAMTOOLS_VER}.tar.bz2 && \
  cd samtools-${SAMTOOLS_VER} && \
  autoconf && ./configure && make && make install && \
  rm -rf /scratch/* && apt-get clean

# Install Scandium and necessary packages/dependencies
COPY scandium-${SCANDIUM_VER} /scratch/scandium-${SCANDIUM_VER}

RUN apt-get -y update && apt-get -y install --no-install-recommends \
  default-libmysqlclient-dev libcurl4-openssl-dev libmariadb-dev libssl-dev && \
  cd /scratch/scandium-${SCANDIUM_VER} && make && make install && \
  rm -rf scratch/* && apt-get clean

ENTRYPOINT ["/usr/bin/sandium"]
CMD ["/bin/sh"]
