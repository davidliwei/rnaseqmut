FROM ubuntu

RUN apt-get update  && apt-get -y install cmake make g++ zlib1g-dev

COPY . /rnaseqmut

RUN mkdir /rnaseqmut/src/bamtools/build
WORKDIR /rnaseqmut/src/bamtools/build

RUN cmake ..
RUN make

WORKDIR /rnaseqmut/src
RUN make

RUN cp /rnaseqmut/bin/rnaseqmut /usr/bin/

ENTRYPOINT /bin/bash
