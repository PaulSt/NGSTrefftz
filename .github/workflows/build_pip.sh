#! /bin/bash
cd workspace/NGSTrefftz
set -e
yum -y update
yum -y install ninja-build fontconfig-devel tk-devel tcl-devel libXmu-devel mesa-libGLU-devel ccache

py=/opt/python/cp39-cp39/bin/python
cd .github/workflows/ && $py fix_auditwheel_policy.py && cd ../..

rm -rf wheelhouse
mkdir wheelhouse

git config --global --add safe.directory '*'

export ORIGINAL_PATH="$PATH"

for pyversion in 38 39 310 311 312
do
    export PYDIR="/opt/python/cp${pyversion}-cp${pyversion}/bin"
    export PATH="$ORIGINAL_PATH:$PYDIR"
    #echo $PYDIR
    #$PYDIR/pip install -U pytest-check numpy wheel scikit-build mkl==2021.* mkl-devel==2021.*
    #$PYDIR/pip install netgen-mesher
    #NETGENDIR=/opt/_internal/cpython-3.9.13/bin

    #rm -rf /home/app/ngstrefftz/make
    rm -rf _skbuild
    $PYDIR/pip install pytest-check numpy wheel scikit-build mkl==2023.* mkl-devel==2023.* setuptools setuptools_scm
    $PYDIR/pip install -r ./.github/workflows/ngsolve_version.txt

    $PYDIR/pip wheel -vvv .
    #cat src/ngstrefftz.egg-info/SOURCES.txt
    #auditwheel repair ngstrefftz*.whl
    rename linux_ manylinux_2_17_x86_64.manylinux2014_ ngstrefftz*.whl
    mv ngstrefftz*.whl wheelhouse/
    rm -rf *.whl
    $PYDIR/pip uninstall -y ngsolve netgen-mesher setuptools setuptools_scm pytest-check numpy wheel scikit-build mkl mkl-devel 


    # avx2 build:
    #rm -rf _skbuild
    #$PYDIR/pip install ngsolve-avx2 --pre
    #NETGEN_ARCH=avx2 $PYDIR/pip wheel -vvv .
    #auditwheel repair ngstrefftz*.whl
    #rm -rf *.whl
    #$PYDIR/pip uninstall -y ngsolve-avx2
    #$PYDIR/pip uninstall -y netgen-mesher-avx2
done

#$PYDIR/pip install -U twine
#$PYDIR/twine upload wheelhouse/*manylinux*.whl
#$PYDIR/twine upload --repository testpypi wheelhouse/*manylinux*.whl
