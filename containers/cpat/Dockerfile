From r-base:3.6.3

RUN apt-get update && apt-get install -y python3-pip procps
RUN cp /usr/bin/python3 /usr/bin/python

RUN pip3 install numpy
RUN pip3 install CPAT

WORKDIR /data/
ENTRYPOINT ["/bin/bash"]
