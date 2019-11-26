FROM ubuntu:bionic

################## METADATA ######################
LABEL base.image="ubuntu:bionic-20190204"
LABEL version="0.1"
LABEL software="GenomeChronicler"
LABEL about.tags="GenomeChronicler, PGP-UK"
LABEL authors="chatzipantsiou@gmail.com,a.guerra@ucl.ac.uk" \
      description="Docker image containing base requirements for PGP-UK GenomeChronicler"


################## MAINTAINER ######################
MAINTAINER <chatzipantsiou@gmail.com>

ENV DEBIAN_FRONTEND noninteractive

# sh
RUN apt-get update && \
    apt-get install -y \
    wget \
    tree \
    nano \
    openjdk-8-jdk \
    apt-transport-https \
    software-properties-common && \
    apt-key adv \
    --keyserver keyserver.ubuntu.com \
    --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 && \
    add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/' && \
    add-apt-repository ppa:linuxuprising/java  && \
    apt-get -y install \ 
    r-base-core \
    tk-dev \
    mesa-common-dev \
    libglu1-mesa-dev && \
    apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/*
 
# CPAN perl
RUN wget -O - http://cpanmin.us | perl - --self-upgrade && \
    cpan File::chdir && \
    cpan Config::General && \
    cpan Data::Dumper && \
    cpan Excel::Writer::XLSX && \
    cpan DBI && \
    cpan DBD::SQLite 

# CLI tools
RUN cd /usr/bin && \
    wget https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2 && \
    tar -vxjf htslib-1.9.tar.bz2 && \
    cd htslib-1.9 && \
    make && \
    cd /usr/bin && \
    wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 && \
    tar -vxjf samtools-1.9.tar.bz2 && \
    cd samtools-1.9 && \
    make && \
    cd /usr/bin && \
    wget https://github.com/samtools/bcftools/releases/download/1.9/bcftools-1.9.tar.bz2 && \
    tar -vxjf bcftools-1.9.tar.bz2 && \
    cd bcftools-1.9 && \
    make

# R
RUN cd / &&\
    R --slave -e 'install.packages("RColorBrewer" , dependencies=TRUE, repos = "http://cran.us.r-project.org")' && \
    R --slave -e 'install.packages("TeachingDemos", dependencies=TRUE, repos = "http://cran.us.r-project.org")'

# LaTeX - update needed for repos after R 
RUN apt-get -y update && apt-get -y install \
    texlive-latex-base \
    texlive-fonts-recommended\
    texlive-latex-extra \
    lmodern

# dev utils
RUN apt-get -y install \
    libssl-dev \
    libcurl4-openssl-dev \ 
    git

# git clone GenomeChronicler 
RUN git clone https://github.com/PGP-UK/GenomeChronicler.git && \
    cd GenomeChronicler && \
    bash SetupMeFirst.sh && \
    ln -sf GenomeChronicler_mainDruid.pl genomechronicler && \
    chmod +x GenomeChronicler_mainDruid.pl && \
    chmod +x genomechronicler && \
    find /GenomeChronicler/scripts -type f -iname "*.pl" -exec chmod +x {} \; && \
    find /GenomeChronicler/scripts -type f -iname "*.R" -exec chmod +x {} \;
    
ENV LC_ALL=C
ENV PATH=$PATH:$PWD
ENV PATH="$PATH:/usr/bin/bcftools-1.9"
ENV PATH="$PATH:/usr/bin/samtools-1.9"
ENV PATH="$PATH:/usr/bin/htslib-1.9"
ENV PATH="$PATH:/GenomeChronicler/"
ENV PATH="$PATH:/GenomeChronicler/scripts"
ENV PATH="$PATH:/GenomeChronicler/genomechronicler"

WORKDIR /usr/bin/