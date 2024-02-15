#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import glob
import os
import sys
# import netgen.version
import site

from skbuild import setup
import skbuild.cmaker
from subprocess import check_output
from distutils.sysconfig import get_python_lib
import pkg_resources

from os.path import dirname, isdir, join
import os
import re
import subprocess

def get_version():
    """
    Gets the current version number.
    If in a git repository, it is the current git tag.
    Otherwise it is the one contained in the PKG-INFO file.
    """
    version_re = re.compile('^Version: (.+)$', re.M)
    d = dirname(__file__)

    if isdir(join(d, '.git')):
        # Get the version using "git describe".
        cmd = 'git describe --always --tags --match v[0-9]*'.split()
        try:
            version = subprocess.check_output(cmd).decode().strip()
        except subprocess.CalledProcessError:
            print('Unable to get version number from git tags')
            exit(1)

        # PEP 386 compatibility
        if '-' in version:
            version = '.post'.join(version.split('-')[:2])

        # Don't declare a version "dirty" merely because a time stamp has
        # changed. If it is dirty, append a ".dev1" suffix to indicate a
        # development revision after the release.
        with open(os.devnull, 'w') as fd_devnull:
            subprocess.call(['git', 'status'],
                            stdout=fd_devnull, stderr=fd_devnull)

        cmd = 'git diff-index --name-only HEAD'.split()
        try:
            dirty = subprocess.check_output(cmd).decode().strip()
        except subprocess.CalledProcessError:
            print('Unable to get git index status')
            exit(1)

        # if dirty != '':
            # version += '.dev1'

        # strip the v for pypi
        version = version[1:]

    else:
        # Extract the version from the PKG-INFO file.
        with open(join(d, 'PKG-INFO')) as f:
            version = version_re.search(f.read()).group(1)

    return version


# import netgen.version
# import ngsolve
# import pkg_resources
# netgen_name = netgen.config.NETGEN_PYTHON_PACKAGE_NAME
# # avx2 = netgen_name.replace('netgen-mesher', '') # keep -avx2 suffix
# avx2 = ''
# name = 'ngstrefftz' + avx2
# # ngsolve_version = pkg_resources.get_distribution("ngsolve").version
with open('.github/workflows/ngsolve_version.txt') as f:
    ngsolve_version = f.readline()
install_requires = [ ngsolve_version, 'ngstents >= 0.0.1.post31' ]


# def install_filter(cmake_manifest):
    # return cmake_manifest

# def _patched_parse_manifests(self):
    # paths = \
        # glob.glob(os.path.join(skbuild.cmaker.CMAKE_BUILD_DIR(), "ngstrefftz", "install_manifest*.txt"))
    # try:
        # return [self._parse_manifest(path) for path in paths][0]
    # except IndexError:
        # return []
   
# # we are using the ngsolve superbuild (to download and build some dependencies)
# # patch the parse_manifests function to point to the actual ngsolve cmake project within the superbuild
# skbuild.cmaker.CMaker._parse_manifests = _patched_parse_manifests

py_install_dir = get_python_lib(1,0,'').replace('\\','/')
print("python install dir:")
print(py_install_dir)

# root_dir = os.path.abspath(os.path.join(netgen.__file__, '../'*(len(py_install_dir.split('/'))+2)))

_cmake_args = ['-DCMAKE_CXX_COMPILER=ngscxx']

packages=["ngstrefftz"]

if 'darwin' in sys.platform:
    _cmake_args += ['-DPY_INSTALL_DIR='+py_install_dir]
elif 'linux' in sys.platform:
    install_requires.append('mkl')
    # packages = []
elif 'win' in sys.platform:
    _cmake_args += ['-DPY_INSTALL_DIR='+py_install_dir]
    install_requires.append('mkl')

cmake_prefix_path = ""
if 'PYDIR' in os.environ:
    cmake_prefix_path += os.environ["PYDIR"]
    _cmake_args += [f'-DPYTHON_EXECUTABLE={os.environ["PYDIR"]}/python3']
    _cmake_args += [f'-DPYTHON_LIBRARY={os.environ["PYDIR"]}/../lib']
    _cmake_args += [f'-DPYTHON_INCLUDE_DIR={os.environ["PYDIR"]}/../include']
if 'CMAKE_PREFIX_PATH' in os.environ:
    cmake_prefix_path += os.environ["CMAKE_PREFIX_PATH"]
_cmake_args += ['-DCMAKE_PREFIX_PATH='+cmake_prefix_path]

print("cmake args:")
print(_cmake_args)

setup(
    name='ngstrefftz',
    version=str(get_version()),
    author='Paul Stocker',
    author_email='p.stocker@math.uni-goettingen.de',
    description='NGSTrefftz is an add-on to NGSolve for Trefftz methods.',
    long_description='NGSTrefftz provides a framework to implement Trefftz finite element spaces for NGSolve, with several Trefftz spaces already implemented. Additionally, Trefftz-DG on tent-pitched meshes for the acoustic wave equation is implemented using meshes provided by ngstents. Furthermore, the package includes an implementation of the embedded Trefftz method.',
    url="https://github.com/PaulSt/ngstrefftz",
    license="LGPL2.1",
    install_requires=install_requires,
    packages=packages,
    package_dir={"ngstrefftz": "src"},
    # cmake_process_manifest_hook=install_filter,
    cmake_args=_cmake_args,
    cmake_source_dir='src',
)
