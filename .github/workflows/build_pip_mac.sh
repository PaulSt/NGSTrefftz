#! /bin/bash
set -e
cd ../..
rm -rf _skbuild dist venv_ngs

export PYDIR=$Python3_ROOT_DIR/bin
#export PYDIR=/Library/Frameworks/Python.framework/Versions/$1/bin

#$PYDIR/python3 --version
#$PYDIR/python3 -m venv venv_ngs
#. venv_ngs/bin/activate

export PATH=/Applications/CMake.app/Contents/bin:$PATH
export NETGEN_Dir=$PYDIR/../lib/python$1/site-packages/netgen/cmake
export NGSolve_Dir=$PYDIR/../lib/python$1/site-packages/ngsolve/cmake
export CMAKE_PREFIX_PATH=$CMAKE_PREFIX_PATH:$NGSolve_Dir:$NETGEN_Dir
#export PYDIR=/Library/Frameworks/Python.framework/Versions/$1/bin
#$PYDIR/python3 -m venv ../venv_ngs
#source ../venv_ngs/bin/activate
$PYDIR/pip3 install scikit-build wheel setuptools==69.5.1 setuptools_scm==8.1.0
#$PYDIR/pip3 install -U pytest-check numpy wheel scikit-build mkl==2022.* mkl-devel==2022.* setuptools

export CMAKE_OSX_ARCHITECTURES='arm64;x86_64'
#export CMAKE_OSX_ARCHITECTURES='x86_64'
$PYDIR/pip3 install -r ./.github/workflows/ngsolve_version.txt
#export NETGEN_Dir=$PWD/venv_ngs/lib/make/netgen/
#export NGSolve_Dir=$PWD/venv_ngs/lib/cmake/ngsolve/
#export NETGEN_Dir=$PWD/venv_ngs/lib/python$1/site-packages/netgen/cmake
#export NGSolve_Dir=$PWD/venv_ngs/lib/python$1/site-packages/ngsolve/cmake
#export CMAKE_PREFIX_PATH=$CMAKE_PREFIX_PATH:$NGSolve_Dir:$NETGEN_Dir

#$PYDIR/pip3 wheel .
$PYDIR/python3 setup.py bdist_wheel --plat-name macosx-10.15-universal2 -d wheelhouse

cat src/ngstrefftz.egg-info/SOURCES.txt
