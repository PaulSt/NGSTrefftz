FROM ubuntu:latest

WORKDIR /home/app

RUN apt-get update && DEBIAN_FRONTEND="noninteractive" apt-get -y install vim python3 python3-pip python3-distutils python3-tk libpython3-dev libxmu-dev tk-dev tcl-dev cmake git g++ libglu1-mesa-dev libblas-dev liblapack-dev
RUN apt-get -y install libboost-all-dev
## try with clang
#RUN apt-get -y install clang
#ENV CC=/usr/bin/clang
#ENV CXX=/usr/bin/clang
## try with gcc 11
#RUN apt-get -y install software-properties-common
#RUN add-apt-repository -y ppa:ubuntu-toolchain-r/test
#RUN apt-get -y install g++-11 gcc-11 
#ENV CC=/usr/bin/gcc-11 
#ENV CXX=/usr/bin/g++-11
#RUN pip3 install numpy

RUN apt-get update
RUN apt-get update
RUN apt-get install -y software-properties-common
RUN add-apt-repository universe
#RUN add-apt-repository ppa:ngsolve/nightly -y
RUN add-apt-repository ppa:ngsolve/ngsolve -y
RUN apt-get install ngsolve -y
#RUN apt-get install npm nodejs -y

ENV PYTHONPATH=/usr/lib/python3/dist-packages/
ENV NETGENDIR=/usr/bin/

COPY . /home/app/tngs
RUN rm -r /home/app/tngs/make && mkdir /home/app/tngs/make
RUN cmake -B/home/app/tngs/make -S/home/app/tngs/src
RUN make -C/home/app/tngs/make
RUN make -C/home/app/tngs/make install

RUN python3 /home/app/tngs/test_trefftz.py -v
