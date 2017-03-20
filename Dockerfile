FROM ubuntu:xenial
MAINTAINER Jon Palmer <nextgenusfs@gmail.com>

RUN apt-get update \
	&& apt-get install -y build-essential curl file g++ git make ruby2.3 ruby2.3-dev uuid-runtime \
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
ENV PATH=/home/linuxbrew/funannotate:/home/linuxbrew/conda/bin:/home/linuxbrew/.linuxbrew/bin:/home/linuxbrew/.linuxbrew/sbin:/home/linuxbrew/augustus/bin:/home/linuxbrew/RepeatModeler:/home/linuxbrew/RepeatMasker:/home/linuxbrew/RepeatMasker/util:/home/linuxbrew/gmes_petap:$PATH \
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

Run brew install blast blat kent-tools mummer hmmer exonerate repeatscout trf rmblast recon trnascan bedtools tbl2asn raxml trimal mafft braker evidencemodeler gag proteinortho 

#install new version of repeatmasker, don't run configure until after homebrew installs dependencies
RUN wget http://www.repeatmasker.org/RepeatMasker-open-4-0-7.tar.gz && \
    tar -zxvf RepeatMasker-open-4-0-7.tar.gz

#create config file
RUN echo -e "\nenv\n/home/linuxbrew/RepeatMasker\n/home/linuxbrew/.linuxbrew/bin/trf\n2\n/home/linuxbrew/.linuxbrew/bin\nY\n4\n/home/linuxbrew/.linuxbrew/bin\nN\n5" > RepeatMasker/repeatmasker.config

RUN cd RepeatMasker && perl ./configure <repeatmasker.config && cd ..

#install repeatmodeler
RUN wget http://www.repeatmasker.org/RepeatModeler-open-1-0-8.tar.gz &&\
    tar -zxvf RepeatModeler-open-1-0-8.tar.gz

#create config file
RUN echo -e "\n/usr/bin/perl\n/home/linuxbrew/RepeatModeler\n/home/linuxbrew/RepeatMasker\n/home/linuxbrew/.linuxbrew/bin\n/home/linuxbrew/.linuxbrew/opt/repeatscout\n/home/linuxbrew/.linuxbrew/bin\n1\n/home/linuxbrew/.linuxbrew/bin\nY\n3" > RepeatModeler/repeatmodeler.config

RUN cd RepeatModeler && perl ./configure <repeatmodeler.config && cd ..

#grab most recent version of funannotate
RUN git clone git://github.com/nextgenusfs/funannotate.git

ENV AUGUSTUS_CONFIG_PATH=/home/linuxbrew/augustus/config \
    EVM_HOME=/home/linuxbrew/.linuxbrew/opt/evidencemodeler \
    GENEMARK_PATH=/home/linuxbrew/gmes_petap \
    BAMTOOLS_PATH=/home/linuxbrew/bamtools/bin

#autosetup funannotate database
RUN funannotate setup -m all -d /home/linuxbrew/DB

RUN mkdir /home/linuxbrew/data

WORKDIR /home/linuxbrew/data


