FROM nextgenusfs/funannotate

LABEL maintainer="Jon Palmer <nextgenusfs@gmail.com>"

USER linuxbrew

WORKDIR /home/linuxbrew

COPY gm_key_64.gz \
    signalp-4.1f.Linux.tar.gz \
    RepBaseRepeatMaskerEdition-20170127.tar.gz \
    /home/linuxbrew/

RUN zcat gm_key_64.gz > /home/linuxbrew/.gm_key && \
    tar -zxvf signalp-4.1f.Linux.tar.gz && \ 
    sed -i 's,/usr/cbs/bio/src/signalp-4.1,/home/linuxbrew/signalp-4.1,g' signalp-4.1/signalp && \
    sed -i 's,#!/usr/bin/perl,#!/usr/bin/env perl,g' signalp-4.1/signalp
    
RUN tar -zxvf RepBaseRepeatMaskerEdition-20170127.tar.gz -C /home/linuxbrew/repeatmasker && \
    rm -rf RepBaseRepeatMaskerEdition-20170127.tar.gz && \
    cd /home/linuxbrew/repeatmasker && perl ./configure < /home/linuxbrew/repeatmasker.txt && \
    cd /home/linuxbrew/repeatmodeler && perl ./configure < /home/linuxbrew/repeatmodeler.txt && \
    funannotate setup -d /home/linuxbrew/DB && \
    mkdir /home/linuxbrew/data

WORKDIR /home/linuxbrew/data

ENTRYPOINT /bin/bash
