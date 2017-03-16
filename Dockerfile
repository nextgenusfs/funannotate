FROM ubuntu:16.04

RUN apt-get clean all && apt-get update && apt-get install -y build-essential apt-utils git wget perl \
    libatlas-base-dev libfreetype6-dev libz-dev libboost-iostreams-dev libpng-dev pkg-config curl \
    python2.7 python-dev python-setuptools python-pip python-numpy python-numpy python-scipy \
    python-pandas python-matplotlib python-biopython python-psutil python-sklearn \
    cmake samtools bedtools zlib1g-dev libc6 libdbd-mysql-perl libdbi-perl libncurses5-dev \
    default-jre nano bioperl cpanminus libcairo2-dev libpango1.0-dev sudo

RUN cpanm Getopt::Long Pod::Usage File::Basename threads threads::shared \
        Thread::Queue Carp Data::Dumper YAML Hash::Merge Logger::Simple Parallel::ForkManager \
        DBI Text::Soundex Scalar::Util::Numeric
        
RUN pip install -U biopython numpy pandas scipy matplotlib seaborn natsort goatools fisher scikit-learn

USER funannotate 
        
WORKDIR /home/funannotate

RUN git clone git://github.com/pezmaster31/bamtools.git && cd bamtools && mkdir build && cd build &&\
    cmake .. && make && sudo make install && cd /usr/include &&  sudo ln -f -s ../local/include/bamtools/ &&\
    cd /usr/lib/ &&  sudo ln -f -s /usr/local/lib/bamtools/libbamtools.* .

RUN git clone https://github.com/nextgenusfs/funannotate.git  

RUN wget http://bioinf.uni-greifswald.de/augustus/binaries/augustus-3.2.3.tar.gz && \
    tar -zxvf augustus-3.2.3.tar.gz && rm augustus-3.2.3.tar.gz && mv augustus-3.2.3 augustus && cd augustus  && make clean && make 
       