#! /bin/bash
set -e
yum -y update
yum -y install ninja-build fontconfig-devel tk-devel tcl-devel libXmu-devel mesa-libGLU-devel ccache

py=/opt/python/cp39-cp39/bin/python
$py fix_auditwheel_policy.py
cd ../..

git config --global --add safe.directory '*'

export ORIGINAL_PATH="$PATH"

for pyversion in 38 39 310
do
    export PYDIR="/opt/python/cp${pyversion}-cp${pyversion}/bin"
    export PATH="$ORIGINAL_PATH:$PYDIR"
    #echo $PYDIR
    #$PYDIR/pip install -U pytest-check numpy wheel scikit-build mkl==2021.* mkl-devel==2021.*
    #$PYDIR/pip install netgen-mesher
    #NETGENDIR=/opt/_internal/cpython-3.9.13/bin

    #rm -rf /home/app/ngstrefftz/make
    rm -rf _skbuild
    $PYDIR/pip install ngsolve

    $PYDIR/pip wheel -vvv .
    auditwheel repair ngstrefftz*.whl
    rm -rf *.whl
    $PYDIR/pip uninstall -y ngsolve
    $PYDIR/pip uninstall -y netgen-mesher

    # avx2 build:
    rm -rf _skbuild
    $PYDIR/pip install ngsolve-avx2
    NETGEN_ARCH=avx2 $PYDIR/pip wheel -vvv .
    auditwheel repair ngstrefftz*.whl
    rm -rf *.whl
    $PYDIR/pip uninstall -y ngsolve-avx2
    $PYDIR/pip uninstall -y netgen-mesher-avx2
done

$PYDIR/pip install -U twine
#$PYDIR/twine upload wheelhouse/*manylinux*.whl
