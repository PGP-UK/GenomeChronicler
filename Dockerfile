FROM ubuntu:20.04

################## METADATA ######################
LABEL version="0.2"
LABEL software="GenomeChronicler"
LABEL about.tags="GenomeChronicler, PGP-UK"
LABEL authors="ismail.moghul@gmail.com,a.guerra@ucl.ac.uk,chatzipantsiou@gmail.com"
LABEL description="Docker image containing base requirements for PGP-UK GenomeChronicler"

ENV LC_ALL=C
ENV DEBIAN_FRONTEND noninteractive

RUN apt-get update \
    && apt-get install -y \
        software-properties-common \
        wget \
        r-base-core \
        tk-dev \
        mesa-common-dev \
        libhts-dev \
        libglu1-mesa-dev \
    && add-apt-repository ppa:linuxuprising/java \
    && apt-get install -y openjdk-8-jdk

# CPAN perl
RUN wget -O - http://cpanmin.us | perl - --self-upgrade \
    && cpan File::chdir \
    && cpan Config::General \
    && cpan Data::Dumper \
    && cpan Excel::Writer::XLSX \
    && cpan DBI \
    && cpan DBD::SQLite

# R
RUN R --slave -e 'install.packages(c("RColorBrewer", "TeachingDemos"))'

# LaTeX - update needed for repos after R
RUN apt-get install -y \
        texlive-latex-base \
        texlive-fonts-recommended \
        texlive-latex-extra \
        lmodern \
    && apt-get clean \
    && apt-get purge \
    && rm -rf /var/lib/apt/lists/* /tmp/*

RUN wget https://github.com/samtools/htslib/releases/download/1.16/htslib-1.16.tar.bz2 \
    && tar -xjf htslib-1.16.tar.bz2 && rm htslib-1.16.tar.bz2 \
    && cd htslib-1.16 && ./configure && make && make install

RUN wget https://github.com/samtools/samtools/releases/download/1.16.1/samtools-1.16.1.tar.bz2 \
    && tar -xjf samtools-1.16.1.tar.bz2 && rm samtools-1.16.1.tar.bz2 \
    && cd samtools-1.16.1 && ./configure && make && make install

RUN wget https://github.com/samtools/bcftools/releases/download/1.16/bcftools-1.16.tar.bz2 \
    && tar -xjf bcftools-1.16.tar.bz2 && rm bcftools-1.16.tar.bz2 \
    && cd bcftools-1.16 && ./configure && make && make install

WORKDIR /GenomeChronicler
COPY GenomeChronicler_mainDruid.pl /GenomeChronicler/
COPY scripts/ /GenomeChronicler/scripts
COPY templates/ /GenomeChronicler/templates/

RUN curl -Ls https://github.com/PGP-UK/GenomeChronicler/releases/download/0.91/reference.tar.gz | tar -zxvf

RUN curl -Ls https://github.com/PGP-UK/GenomeChronicler/releases/download/0.91/software.tar.gz | \
    tar zxvf - software.linux && mv software.linux software

ENV PATH "$PATH:/GenomeChronicler/scripts"

CMD /GenomeChronicler/scripts/genomechronicler
