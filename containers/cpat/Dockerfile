From python:3.11.14-trixie

RUN apt-get update && apt-get install -y python3-pip procps
RUN cp /usr/bin/python3 /usr/bin/python

RUN pip3 install numpy
RUN pip3 install CPAT

WORKDIR /data/
