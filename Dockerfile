FROM parseq/stepik-variant-calling-tools

RUN apt-get clean | apt-get update
RUN apt-get install -qy nano

ENV EDITOR nano

#VOLUME /home/stepik
ARG UID="1000"

RUN useradd -ms /bin/bash $UID
USER $UID

WORKDIR /home/stepik


CMD $EDITOR
CMD ["nano", " " ]
