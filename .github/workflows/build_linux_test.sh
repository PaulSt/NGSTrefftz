#! /bin/bash
cd workspace
yum -y update && yum -y install git ninja-build fontconfig-devel tk-devel tcl-devel libXmu-devel mesa-libGLU-devel ccache python-devel lapack blas 

export pyversion=310
export PYDIR="/opt/python/cp${pyversion}-cp${pyversion}/bin"
export PATH=$PATH:$PYDIR

$PYDIR/pip install ngsolve --pre
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/_internal/cpython-3.10.9/lib/

cmake -DCMAKE_CXX_COMPILER=ngscxx -DNETGENDIR=$PYDIR/.. -DCMAKE_PREFIX_PATH=$PYDIR/..  -DPYTHON_EXECUTABLE=$env.PYDIR/python3 -DPYTHON_LIBRARY=$env.PYDIR/../lib -DPYTHON_INCLUDE_DIR=$env.PYDIR/../include -B ./NGSTrefftz/make -S ./NGSTrefftz/src/
make -C ./NGSTrefftz/make
env CTEST_OUTPUT_ON_FAILURE=1 make -C ./NGSTrefftz/make test
make -C ./NGSTrefftz/make install  
