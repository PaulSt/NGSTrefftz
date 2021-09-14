#!/bin/bash
sudo apt-get update && sudo apt-get -y install python3 python3-distutils python3-tk libpython3-dev libxmu-dev tk-dev tcl-dev cmake git g++ libglu1-mesa-dev
sudo apt-get -y install libblas-dev liblapack-dev libboost-all-dev liblapacke-dev 
export BASEDIR=~/ngsuite
mkdir -p $BASEDIR
cd $BASEDIR && git clone https://github.com/PaulSt/ngsolve ngsolve-src
cd $BASEDIR/ngsolve-src && git submodule update --init --recursive
mkdir $BASEDIR/ngsolve-build
mkdir $BASEDIR/ngsolve-install
cd $BASEDIR/ngsolve-build
cmake -DCMAKE_INSTALL_PREFIX=${BASEDIR}/ngsolve-install ${BASEDIR}/ngsolve-src
make
make install
echo "export NETGENDIR=${BASEDIR}/ngsolve-install/bin" >> ~/.bashrc
echo "export PATH=\$NETGENDIR:\$PATH" >> ~/.bashrc
export PYTHONPATH_TMP=`python3 -c "from distutils.sysconfig import get_python_lib; print(get_python_lib(1,0,''))"`
echo "export PYTHONPATH=\$NETGENDIR/../${PYTHONPATH_TMP}:\$PATH" >> ~/.bashrc
echo "{NETGENDIR}={${BASEDIR}/ngsolve-install/bin}" >> $GITHUB_ENV
echo "{PATH}={\$NETGENDIR:\$PATH}" >> $GITHUB_ENV
echo "{PYTHONPATH}={\$NETGENDIR/../${PYTHONPATH_TMP}:\$PATH}" >> $GITHUB_ENV
export NETGENDIR=${BASEDIR}/ngsolve-install/bin
export PATH=$NETGENDIR:$PATH
export PYTHONPATH=\$NETGENDIR/../${PYTHONPATH_TMP}:\$PATH
cmake -B ${GITHUB_WORKSPACE}/build -S ${GITHUB_WORKSPACE}/src/trefftz -DCMAKE_BUILD_TYPE=${BUILD_TYPE} 
cd ${GITHUB_WORKSPACE}/build 
make
