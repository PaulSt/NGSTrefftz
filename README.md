# NGSTrefftz
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/PaulSt/NGSTrefftz/HEAD?filepath=doc%2Fnotebooks%2Findex.ipynb)
[![Docker Image Version (latest semver)](https://img.shields.io/docker/v/paulstdocker/ngstrefftz?label=docker&logo=docker&sort=semver)](https://hub.docker.com/r/paulstdocker/ngstrefftz)
[![PyPI](https://img.shields.io/pypi/v/ngstrefftz?color=blue&logo=pypi)](https://pypi.org/project/ngstrefftz/)
[![GitHub Workflow Status](https://img.shields.io/github/workflow/status/PaulSt/NGSTrefftz/build?logo=github)](https://github.com/PaulSt/NGSTrefftz/actions)

NGSTrefftz provides a framework to implement Trefftz finite element spaces for [NGSolve](https://www.ngsolve.com), with several Trefftz spaces already implemented. Additionally, Trefftz-DG on tent-pitched meshes for the acoustic wave equation is implemented using meshes provided by ![ngstents](https://github.com/jayggg/ngstents). Furthermore, the package includes an implementation of the embedded Trefftz method.

## Try it out!
You can try out some jupyter notebooks:
* Launch the Binder here:   
  [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/PaulSt/NGSTrefftz/HEAD?filepath=doc%2Fnotebooks%2Findex.ipynb)
* Or run the docker locally (you need to have docker installed):
```bash
git clone https://github.com/PaulSt/NGSTrefftz
cd NGSTrefftz && docker build -t ngstrefftz_jupyter .
docker run -p 8888:8888 ngstrefftz_jupyter
```

## Installing the package
You can either:
 * install using pip
```bash
pip install ngstrefftz
```
 * or build from source
```bash
git clone --recursive https://github.com/PaulSt/NGSTrefftz
mkdir ./NGSTrefftz/make && cd ./NGSTrefftz/make
cmake ../src && make install
```
### Dependencies
You need to have [NGSolve](https://www.ngsolve.com/) installed. To build the package you need
* cmake  >= 3.1
* gcc >= 9 or clang >= 10
* ngsolve >= 6.2

To access the newest features the nightly version of NGSolve works best and lapack >= 3.9 is required.

## Papers using the code
* Tent pitching and Trefftz-DG method for the acoustic wave equation  
[![arXiv](https://img.shields.io/badge/arXiv-1907.02367-b31b1b.svg)](https://arxiv.org/abs/1907.02367)
* A space-time quasi-Trefftz DG method for the wave equation with piecewise-smooth coefficients  
[![arXiv](https://img.shields.io/badge/arXiv-2011.04617-b31b1b.svg)](https://arxiv.org/abs/2011.04617)
* Embedded Trefftz discontinuous Galerkin methods  
[![arXiv](https://img.shields.io/badge/arXiv-2201.07041-b31b1b.svg)](https://arxiv.org/abs/2201.07041)


![](.github/wave.gif)

