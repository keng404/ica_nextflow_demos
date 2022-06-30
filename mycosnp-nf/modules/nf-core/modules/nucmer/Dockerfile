FROM library/ubuntu:20.04
MAINTAINER Ken Eng <keneng87@gmail.com>
LABEL    software="mummer4" \ 
    container="mummer4" \ 
    about.summary="Efficient sequence alignment of full genomes" \ 
    about.home="https://github.com/mummer4" \ 
    software.version="4.0.0rc1" \ 
    version="1" \ 
    extra.identifiers.biotools="mummer" \ 
    about.copyright=" 2018 by Marçais G and The Institute for Genomic Research." \ 
    about.license="Artistic" \ 
    about.license_file="/usr/share/doc/mummer/copyright" \ 
    extra.binaries="/usr/bin/combineMUMs,/usr/bin/delta-filter,/usr/bin/delta2blocks,/usr/bin/delta2maf,/usr/bin/dnadiff,/usr/bin/exact-tandems,/usr/bin/gaps,/usr/bin/mapview,/usr/bin/mgaps,/usr/bin/mummer,/usr/bin/mummer-annotate,/usr/bin/mummerplot,/usr/bin/nucmer,/usr/bin/nucmer2xfig,/usr/bin/promer,/usr/bin/repeat-match,/usr/bin/run-mummer1,/usr/bin/run-mummer3,/usr/bin/show-aligns,/usr/bin/show-coords,/usr/bin/show-diff,/usr/bin/show-snps,/usr/bin/show-tiling" \ 
    extra.identifiers.conda="mummer" \ 
    about.tags="field::biology, field::biology:bioinformatics, implemented-in::c++,:commandline, role::program, scope::utility, use::comparing,:TODO,:plaintext, works-with::TODO" 

ENV DEBIAN_FRONTEND noninteractive
RUN apt-get update -y && apt-get upgrade -y && apt-get install -y wget tar g++ libssl-dev libffi-dev  build-essential yaggo gnuplot && apt-get clean && apt-get purge && rm -rf /var/lib/apt/lists/* /tmp/*
### getting python3
RUN apt-get update -y && apt-get install -y python3-pip python3-dev
### install visation packages
RUN apt-get update -y && apt-get install -y fig2dev xfig
RUN apt-get update -y && apt-get install -y python-dev ruby-full
RUN cd /opt && \ 
    wget --no-check-certificate https://github.com/mummer4/mummer/releases/download/v4.0.0rc1/mummer-4.0.0rc1.tar.gz && \
   tar -xvf mummer-4.0.0rc1.tar.gz && \
    cd /opt//mummer-4.0.0rc1 && \ 
   ./configure  --enable-python-binding --enable-perl-binding --enable-ruby-binding && \
    make -j && \
    make install && \
    ldconfig