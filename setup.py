#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import re
import sys
import platform
import subprocess
import warnings
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
from distutils.version import LooseVersion

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


class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=''):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)
class CMakeBuild(build_ext):
    def run(self):
        try:
            out = subprocess.check_output(['cmake', '--version'])
        except OSError:
            raise RuntimeError("CMake must be installed to build the following extensions: " +
                               ", ".join(e.name for e in self.extensions))
        if "NETGENDIR" not in os.environ:
            warnings.warn("Could not find NETGENDIR")
        try:
            import ngsolve
        except:
            raise RuntimeError("Could not run NGSolve, is it installed? Did you set pythonpath correctly?")

        for ext in self.extensions:
            self.build_extension(ext)
    def build_extension(self, ext):
        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))
        # required for auto-detection of auxiliary "native" libs
        if not extdir.endswith(os.path.sep):
            extdir += os.path.sep
        cmake_args = ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=' + extdir,
                      '-DPYTHON_EXECUTABLE=' + sys.executable,
                      '-DCMAKE_CXX_COMPILER=ngscxx']
        if 'PYDIR' in os.environ:
            cmake_args += [f'-DCMAKE_PREFIX_PATH={os.environ["PYDIR"]}/..']
            cmake_args += [f'-DPYTHON_EXECUTABLE={os.environ["PYDIR"]}/python3']
            cmake_args += [f'-DPYTHON_LIBRARY={os.environ["PYDIR"]}/../lib']
            cmake_args += [f'-DPYTHON_INCLUDE_DIR={os.environ["PYDIR"]}/../include']
        cfg = 'Debug' if self.debug else 'Release'
        build_args = ['--config', cfg]
        cmake_args += ['-DCMAKE_BUILD_TYPE=' + cfg]
        build_args += ['--', '-j2']
        env = os.environ.copy()
        env['CXXFLAGS'] = '{} -DVERSION_INFO=\\"{}\\"'.format(env.get('CXXFLAGS', ''),
                                                              self.distribution.get_version())
        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)

        subprocess.check_call(['cmake', ext.sourcedir] + cmake_args, cwd=self.build_temp, env=env)
        subprocess.check_call(['cmake', '--build', '.'] + build_args, cwd=self.build_temp)

        # subprocess.check_call(['mv', '_trefftz.so', 'ngstrefftz'], cwd=self.build_lib)
        # subprocess.check_call(['mkdir', 'ngstents'], cwd=self.build_lib)
        # subprocess.check_call(['mv', '_pytents.so', 'ngstents'], cwd=self.build_lib)

import netgen.version
import ngsolve
import pkg_resources
netgen_name = netgen.config.NETGEN_PYTHON_PACKAGE_NAME
avx2 = netgen_name.replace('netgen-mesher', '') # keep -avx2 suffix
name = 'ngstrefftz' + avx2
ngsolve_version = pkg_resources.get_distribution("ngsolve").version
install_requires = [ 'ngsolve'+avx2+'>='+ngsolve_version ]

if sys.argv[1] == "sdist":
    package_data = {"ngstrefftz": ["*"
                                ,"../test/*"
                                ,"../external_dependencies/ngstents/*"\
                                ,"../external_dependencies/ngstents/src/*"\
                                ,"../external_dependencies/ngstents/py/*"\
                                ]}
    name += "-src"
else:
    package_data = {}

setup(
    name=name,
    version=str(get_version()),
    author='Paul Stocker',
    author_email='p.stocker@math.uni-goettingen.de',
    description='NGSTrefftz is an add-on to NGSolve for Trefftz methods.',
    long_description='NGSTrefftz provides a framework to implement Trefftz finite element spaces for NGSolve, with several Trefftz spaces already implemented. Additionally, Trefftz-DG on tent-pitched meshes for the acoustic wave equation is implemented using meshes provided by ngstents. Furthermore, the package includes an implementation of the embedded Trefftz method.',
    url="https://github.com/PaulSt/ngstrefftz",
    install_requires=install_requires,
    ext_modules=[CMakeExtension('ngstrefftz_py','src')],
    cmdclass=dict(build_ext=CMakeBuild),
    packages=["ngstrefftz"],
    package_dir={"ngstrefftz": "src"},
    python_requires='>=3.8',
    package_data=package_data,
)
