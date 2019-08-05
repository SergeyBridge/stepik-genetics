FROM parseq/stepik-variant-calling-tools

RUN apt-get clean | apt-get update
RUN apt-get install -qy snakemake

ENV PIPELINE snakemake

VOLUME /home/snakemake_dir
ARG UID="1000"

RUN useradd -ms /bin/bash $UID
USER $UID

WORKDIR /home/snakemake_dir


#CMD $EDITOR
ENTRYPOINT $PIPELINE