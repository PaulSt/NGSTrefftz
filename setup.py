#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from skbuild import setup
from setuptools_scm import get_version

_cmake_args = []
version = get_version().split('+')[0]
with open('.github/workflows/ngsolve_version.txt') as f:
    ngsolve_version = f.readline()
install_requires = [ ngsolve_version, 'ngstents >= 0.0.2.dev39' ]


import os, sys
# try:
    # import netgen
    # from distutils.sysconfig import get_python_lib
    # py_install_dir = get_python_lib(1,0,'').replace('\\','/')
    # root_dir = os.path.abspath(os.path.join(netgen.__file__, '../'*(len(py_install_dir.split('/'))+2)))
# except:
root_dir = sys.prefix
print(f'root_dir: {root_dir}')
if 'darwin' in sys.platform:
    _cmake_args += [
        '-DBUILD_STUB_FILES=ON',
    ]
elif 'linux' in sys.platform:
    _cmake_args += [
        '-DUSE_MKL:BOOL=ON',
        f'-DMKL_ROOT:PATH={root_dir}',
        f'-DMKL_LIBRARY:PATH={root_dir}/lib/libmkl_rt.so.2',
        f'-DMKL_INCLUDE_DIR:PATH={root_dir}/include',
        '-DUSE_CUDA=ON',
        '-DCMAKE_CUDA_ARCHITECTURES=all',
        '-DBUILD_STUB_FILES=ON',
    ]
    install_requires.append('mkl')
    packages = []
elif 'win' in sys.platform:
    _cmake_args += [
        '-DUSE_MKL:BOOL=ON',
        f'-DMKL_ROOT:PATH={root_dir}',
        f'-DMKL_LIBRARY:PATH={root_dir}/Library/lib/mkl_rt.lib',
        f'-DMKL_INCLUDE_DIR:PATH={root_dir}/Library/include',
        f'-DNGSOLVE_INSTALL_DIR_TCL:PATH=Scripts',
        '-DBUILD_STUB_FILES=OFF',
    ]
    install_requires.append('mkl')


setup(
    name='ngstrefftz',
    version=version,
    author='Paul Stocker',
    author_email='p.stocker@math.uni-goettingen.de',
    description='NGSTrefftz is an add-on to NGSolve for Trefftz methods.',
    long_description='NGSTrefftz provides a framework to implement Trefftz finite element spaces for NGSolve, with several Trefftz spaces already implemented. Additionally, Trefftz-DG on tent-pitched meshes for the acoustic wave equation is implemented using meshes provided by ngstents. Furthermore, the package includes an implementation of the embedded Trefftz method.',
    url="https://github.com/PaulSt/ngstrefftz",
    license="LGPL2.1",
    install_requires=install_requires,
    packages=["ngstrefftz"],
    package_dir={"ngstrefftz": "src"},
    # cmake_process_manifest_hook=install_filter,
    cmake_args=_cmake_args,
    cmake_source_dir='src',
)
