FROM ubuntu:xenial
MAINTAINER Jon Palmer <nextgenusfs@gmail.com>

RUN apt-get update \
	&& apt-get install --fix-missing -y build-essential curl file g++ git make ruby2.3 ruby2.3-dev uuid-runtime \
	bioperl sudo wget libboost-all-dev libncurses5-dev zlib1g-dev \
	&& ln -sf ruby2.3 /usr/bin/ruby \
	&& ln -sf gem2.3 /usr/bin/gem

RUN localedef -i en_US -f UTF-8 en_US.UTF-8 \
	&& useradd -m -s /bin/bash linuxbrew \
	&& echo 'linuxbrew ALL=(ALL) NOPASSWD:ALL' >>/etc/sudoers
RUN git clone https://github.com/Linuxbrew/brew.git /home/linuxbrew/.linuxbrew \
    && chown -R linuxbrew: /home/linuxbrew/.linuxbrew \
	&& cd /home/linuxbrew/.linuxbrew \
	&& git remote set-url origin https://github.com/Linuxbrew/brew.git

USER linuxbrew
WORKDIR /home/linuxbrew
ENV PATH=/home/linuxbrew/funannotate:/home/linuxbrew/conda/bin:/home/linuxbrew/.linuxbrew/bin:/home/linuxbrew/.linuxbrew/sbin:/home/linuxbrew/augustus/bin:/home/linuxbrew/RepeatModeler:/home/linuxbrew/RepeatMasker:/home/linuxbrew/RepeatMasker/util:/home/linuxbrew/data/gm_et_linux_64/gmes_petap:/home/linuxbrew/data/signalp-4.1:$PATH \
	SHELL=/bin/bash

RUN brew doctor || true

ENV PERL5LIB=/usr/local/lib/perl5/site_perl:${PERL5LIB}

RUN sudo cpan -i Getopt::Long Pod::Usage File::Basename threads threads::shared \
        Thread::Queue Carp Data::Dumper YAML Hash::Merge Logger::Simple Parallel::ForkManager \
        DBI Text::Soundex Scalar::Util::Numeric
	
RUN brew tap homebrew/science 
RUN brew tap nextgenusfs/tap
RUN brew tap homebrew/dupes
RUN brew update

#install new cmake version for compliation of bamtools
RUN brew install cmake

#gmap-gsnap, bamtools, augustus are failling, try to install it separately
RUN git clone git://github.com/pezmaster31/bamtools.git \
    && cd bamtools && mkdir build && cd build &&\
    cmake .. && make && sudo make install && cd /usr/include &&  sudo ln -f -s ../local/include/bamtools/ &&\
    cd /usr/lib/ &&  sudo ln -f -s /usr/local/lib/bamtools/libbamtools.* .
    
RUN wget http://bioinf.uni-greifswald.de/augustus/binaries/augustus-3.2.3.tar.gz && \
    tar -zxvf augustus-3.2.3.tar.gz && rm augustus-3.2.3.tar.gz && mv augustus-3.2.3 augustus \
    && cd augustus && make clean && make
    
RUN wget http://research-pub.gene.com/gmap/src/gmap-gsnap-2015-12-31.v10.tar.gz \
    && tar -zxvf gmap-gsnap-2015-12-31.v10.tar.gz && rm gmap-gsnap-2015-12-31.v10.tar.gz \
    && mv gmap-2015-12-31// gmap && cd gmap/  && ./configure && make && sudo make install

#conda install
RUN wget --quiet https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    sudo /bin/bash ~/miniconda.sh -b -p /home/linuxbrew/conda && \
    rm ~/miniconda.sh

RUN sudo chown -R linuxbrew: /home/linuxbrew/conda

RUN conda update -y conda
RUN conda install -y numpy pandas scipy seaborn natsort scikit-learn psutil biopython
RUN conda install -y -c etetoolkit ete3 ete3_external_apps
RUN pip install --upgrade goatools fisher

RUN brew install blat kent-tools mummer hmmer exonerate repeatscout trf rmblast recon repeatmasker repeatmodeler trnascan bedtools tbl2asn raxml trimal mafft braker evidencemodeler gag proteinortho diamond

#currently blast v2.6.0 tblastn is giving strange results, default back to v2.2.31
RUN rm /home/linuxbrew/.linuxbrew/bin/tblastn

#grab most recent version of funannotate
RUN git clone git://github.com/nextgenusfs/funannotate.git && \
    cd funannotate && \
    git checkout 5aacd0632c7ed99bde5f3e816d0bdeca751228b9 && \
    cd ..

ENV AUGUSTUS_CONFIG_PATH=/home/linuxbrew/augustus/config \
    EVM_HOME=/home/linuxbrew/.linuxbrew/opt/evidencemodeler \
    GENEMARK_PATH=/home/linuxbrew/data/gm_et_linux_64/gmes_petap \
    BAMTOOLS_PATH=/home/linuxbrew/bamtools/bin

#autosetup funannotate database
RUN funannotate setup -m all -d /home/linuxbrew/DB

RUN mkdir /home/linuxbrew/data

WORKDIR /home/linuxbrew/data


