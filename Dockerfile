FROM nextgenusfs/funannotate

WORKDIR /home/linuxbrew

COPY gm_et_linux_64.tar.gz /home/linuxbrew/gm_et_linux_64.tar.gz
COPY gm_key_64.gz /home/linuxbrew/gm_key_64.gz
COPY signalp-4.1c.Linux.tar.Z /home/linuxbrew/signalp-4.1c.Linux.tar.Z

RUN tar zxvf gm_et_linux_64.tar.gz && \
    mv /home/linuxbrew/gm_et_linux_64/gmes_petap $GENEMARK_PATH && \
    gunzip gm_key_64.gz && \
    cp gm_key_64 /home/linuxbrew/.gm_key && \
    rm -rf /home/linuxbrew/gm_et_linux_64 gm_key_64 gm_et_linux_64.tar.gz

RUN zcat signalp-4.1c.Linux.tar.Z | tar -xvf - && \
    sed -i 's,/usr/cbs/bio/src/signalp-4.1,/home/linuxbrew/data/signalp-4.1,g' signalp-4.1/signalp && \
    mv /home/linuxbrew/signalp-4.1 /home/linuxbrew/data/signalp-4.1 

COPY RepBaseRepeatMaskerEdition-20170127.tar.gz /home/linuxbrew/RepBaseRepeatMaskerEdition-20170127.tar.gz

RUN tar zxvf RepBaseRepeatMaskerEdition-20170127.tar.gz -C /home/linuxbrew/.linuxbrew/opt/repeatmasker/libexec
RUN brew postinstall repeatmasker

RUN funannotate setup -m dep -d /home/linuxbrew/DB