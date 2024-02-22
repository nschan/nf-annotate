FROM ubuntu:xenial

# install dependencies and clean up apt garbage
RUN apt-get update && apt-get install --no-install-recommends -y \
 libncurses5-dev \
 libbz2-dev \
 liblzma-dev \
 libcurl4-gnutls-dev \
 zlib1g-dev \
 libssl-dev \
 gcc \
 wget \
 make \
 perl \
 bzip2 \
 gnuplot \
 ca-certificates \
 gawk && \
 apt-get autoclean && rm -rf /var/lib/apt/lists/*
 
# for easy upgrade later. ARG variables only persist during image build
ARG SAMTOOLSVER=1.17
RUN wget https://github.com/samtools/samtools/releases/download/${SAMTOOLSVER}/samtools-${SAMTOOLSVER}.tar.bz2 && \
 tar -xjf samtools-${SAMTOOLSVER}.tar.bz2 && \
 rm samtools-${SAMTOOLSVER}.tar.bz2 && \
 cd samtools-${SAMTOOLSVER} && \
 ./configure && \
 make && \
 make install && \
 cd ..


ARG MINIMAP2_VER=2.26
RUN wget https://github.com/lh3/minimap2/releases/download/v${MINIMAP2_VER}/minimap2-${MINIMAP2_VER}_x64-linux.tar.bz2 
RUN tar -xjf minimap2-${MINIMAP2_VER}_x64-linux.tar.bz2 --no-same-owner && \
 rm minimap2-${MINIMAP2_VER}_x64-linux.tar.bz2

ENV PATH="${PATH}:/minimap2-${MINIMAP2_VER}_x64-linux"
