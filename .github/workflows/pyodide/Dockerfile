#FROM ubuntu:23.04
#ENV DEBIAN_FRONTEND=noninteractive
#ENV PYODIDE_VERSION=0.23.4
#RUN apt-get update && apt-get -y upgrade
#RUN apt-get update && apt-get -y install \
        #nodejs \
        #ccache \
        #python3-pip \
        #python3 \
        #git \
        #make \
        #pkg-config \
        #g++ \
        #lbzip2 \
        #xz-utils \
        #autoconf \
        #libtool \
        #unzip \
        #zip \
        #xxd \
        #wget

#WORKDIR /root
#RUN git clone https://github.com/emscripten-core/emsdk.git
#RUN cd emsdk && ./emsdk install latest && ./emsdk activate latest
#RUN wget https://github.com/pyodide/pyodide/releases/download/${PYODIDE_VERSION}/xbuildenv-${PYODIDE_VERSION}.tar.bz2
#RUN tar xvf xbuildenv-${PYODIDE_VERSION}.tar.bz2
#RUN ls

#RUN mkdir /root/output
#WORKDIR /root/output
#RUN wget https://github.com/pyodide/pyodide/releases/download/${PYODIDE_VERSION}/pyodide-core-${PYODIDE_VERSION}.tar.bz2
#RUN tar xvf pyodide-core-${PYODIDE_VERSION}.tar.bz2
#RUN ln -s /usr/bin/python3 /usr/bin/python
#RUN rm /usr/lib/python3.11/EXTERNALLY-MANAGED
#RUN pip install cmake
#RUN git clone --recursive https://github.com/NGSolve/ngsolve /root/ngsolve
#RUN mv /root/ngsolve/tests/gitlab-ci/pyodide/build_opencascade.sh /root/
#RUN bash /root/build_opencascade.sh
#RUN mkdir /root/ngsolve/build
#WORKDIR /root/ngsolve/build
#RUN bash /root/ngsolve/tests/gitlab-ci/pyodide/build_in_docker.sh

### load above steps 
FROM paulstdocker/ngsolve_pyodide:latest

### build ngstrefftz
#RUN git clone --recursive https://github.com/PaulSt/NGSTrefftz /root/ngstrefftz
#RUN cd /root/ngstrefftz && git checkout ngsolve_update
ADD . /root/ngstrefftz
WORKDIR /root/ngstrefftz/build
RUN bash ../.github/workflows/pyodide/build_in_docker.sh

### output
WORKDIR /root/output/pyodide
RUN rm /opt/netgen/python/*/*.a
WORKDIR /opt/netgen/python
RUN zip -r pyngcore.zip pyngcore
RUN zip -r netgen.zip netgen
RUN zip -r ngsolve.zip ngsolve
RUN zip -r ngstrefftz.zip ngstrefftz
RUN zip -r ngstents.zip ngstents
RUN cp /opt/netgen/python/*.zip /root/output/pyodide
WORKDIR /root/output/pyodide
RUN cp /root/ngsolve/tests/gitlab-ci/pyodide/generate_repodata.js .
RUN node generate_repodata.js
RUN rm *.whl
RUN python3 /root/ngstrefftz/.github/workflows/pyodide/merge.py
RUN cd /root/output/ && tar -cvjSf ngstrefftz_pyodide.tar.bz2 pyodide

