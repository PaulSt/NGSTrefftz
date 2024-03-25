### load above steps 
FROM paulstdocker/ngsolve_pyodide:latest

### build ngstrefftz
#RUN git clone --recursive https://github.com/PaulSt/NGSTrefftz /root/ngstrefftz
#RUN cd /root/ngstrefftz && git checkout ngsolve_update
ADD . /root/ngstrefftz
RUN mkdir -p /root/ngstrefftz/build
WORKDIR /root/ngstrefftz/build
RUN --mount=type=cache,target=/ccache/ bash /root/ngstrefftz/.github/workflows/pyodide/build_in_docker.sh

### output
WORKDIR /root/output/pyodide
RUN rm /opt/netgen/python/*/*.a
RUN apt install -y zip
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
