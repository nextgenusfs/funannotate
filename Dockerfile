FROM nextgenusfs/funannotate

LABEL maintainer="Jon Palmer <nextgenusfs@gmail.com>"

WORKDIR /home/linuxbrew

COPY gm_et_linux_64.tar.gz \
    gm_key_64.gz \
    signalp-4.1c.Linux.tar.Z \
    RepBaseRepeatMaskerEdition-20170127.tar.gz \
    /home/linuxbrew/
    
RUN tar -zxvf gm_et_linux_64.tar.gz && \
    gunzip gm_key_64.gz &&  cp gm_key_64 /home/linuxbrew/.gm_key && rm -rf gm_et_linux_64.tar.gz && \
    zcat signalp-4.1c.Linux.tar.Z | tar -xvf - && \
    sed -i 's,/usr/cbs/bio/src/signalp-4.1,/home/linuxbrew/signalp-4.1,g' signalp-4.1/signalp && \
    tar -zxvf RepBaseRepeatMaskerEdition-20170127.tar.gz -C /home/linuxbrew/repeatmasker && \
    cd /home/linuxbrew/repeatmasker && perl ./configure < /home/linuxbrew/repeatmasker.txt && \
    cd /home/linuxbrew/repeatmodeler && perl ./configure < /home/linuxbrew/repeatmodeler.txt && \
    funannotate setup -d /home/linuxbrew/DB --update

WORKDIR /data

ENTRYPOINT "/home/linuxbrew/startup.sh" && /bin/bash