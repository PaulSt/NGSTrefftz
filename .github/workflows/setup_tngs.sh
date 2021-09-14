#!/bin/bash
source ~/.bashrc
export NETGENDIR=/home/runner/ngsuite/ngsolve-install/bin
export PATH=$NETGENDIR:$PATH
export PYTHONPATH=$NETGENDIR/../lib/python3/dist-packages:$PATH
cmake -B ${GITHUB_WORKSPACE}/build -S ${GITHUB_WORKSPACE}/src/trefftz -DCMAKE_BUILD_TYPE=${BUILD_TYPE} 
cd ${GITHUB_WORKSPACE}/build 
make
