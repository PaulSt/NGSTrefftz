# NGSTrefftz
**an add-on to NGSolve for Trefftz methods**

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/PaulSt/NGSTrefftz/HEAD?filepath=docs%2Fnotebooks%2Findex.ipynb)
[![Docker Image Version (latest semver)](https://img.shields.io/docker/v/paulstdocker/ngstrefftz?label=docker&logo=docker&sort=semver)](https://hub.docker.com/r/paulstdocker/ngstrefftz)
[![PyPI](https://img.shields.io/pypi/v/ngstrefftz?color=blue&logo=pypi)](https://pypi.org/project/ngstrefftz/)
[![GitHub Workflow Status](https://github.com/PaulSt/NGSTrefftz/actions/workflows/build.yml/badge.svg?branch=main)](https://github.com/PaulSt/NGSTrefftz/actions/workflows/build.yml)
[![status](https://joss.theoj.org/papers/c2f4e85b118c22b81aa27d7799265409/status.svg)](https://joss.theoj.org/papers/c2f4e85b118c22b81aa27d7799265409)
[![docs](https://img.shields.io/badge/docs-NGSTrefftz-blue?logo=readthedocs)](https://paulst.github.io/NGSTrefftz/)

NGSTrefftz provides a framework to implement Trefftz finite element spaces for [NGSolve](https://www.ngsolve.org), with several Trefftz spaces already implemented. Additionally, Trefftz-DG on tent-pitched meshes for the acoustic wave equation is implemented using meshes provided by [ngstents](https://github.com/jayggg/ngstents). Furthermore, the package includes an implementation of the embedded Trefftz method.

## Try it out!
You can try out some jupyter notebooks:
* Launch the Binder here:   
  [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/PaulSt/NGSTrefftz/HEAD?filepath=docs%2Fnotebooks%2Findex.ipynb)
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
cmake ../ && make install
```

### Dependencies
 * if you are using pip to install the package you can install the newest ngsolve and other dependencies using
```bash
pip install ngsolve --pre
```
 * if you build the package from source the following software needs to be installed
   * [cmake](https://cmake.org/) >= 3.1
   * [gcc](https://gcc.gnu.org/) >= 9 or [clang](https://clang.llvm.org/) >= 10
   * [lapack](http://www.netlib.org/lapack/) >= 3.9 (required to access the newest features)
   * [ngsolve](https://www.ngsolve.org) >= 6.2  (to access the newest features the nightly version of NGSolve works best)

## News
⚠️ Apr, 2025: `TrefftzEmbedding` has a new interface, please check the documentation for details.

🚀 Jul, 2024: Conforming Trefftz embedding implementation by [@johann-cm](https://github.com/johann-cm) 

⚠️ Oct, 2022: With v0.2.0 the git history has undergone a major cleanup, please make sure to clone the repo anew.

🚀 Oct, 2022: New and improved implementation of the embedded Trefftz method via `EmbeddedTrefftzFES`!

🚀 Aug, 2022: Three different [pip](https://pypi.org/search/?q=ngstrefftz)-installer available, now using wheels! 

🚀 Mar, 2022: NGSTrefftz now has a [website](https://paulst.github.io/NGSTrefftz/)! 

⚠️ Feb, 2022: If you are using NGSolve nightly releases: [NGSolve@eda758d](https://github.com/NGSolve/ngsolve/commit/eda758d99483888851913d8a5c9aff4d0cbc9cc2) breaks a dependency and [NGSolve@3d52ecd](https://github.com/NGSolve/ngsolve/commit/3d52ecd615f2b7c409219eebaba99288ea19c1bc) produces import issue. Make sure to update ngstrefftz submodules and move to newest ngsolve version, at least [NGSolve@5839a09](https://github.com/NGSolve/ngsolve/commit/5839a09810235a938bd85807d8e638d3a0b6c69d).

🚀 Jan, 2022: NGSTrefftz is now available via pip! 

🚀 Nov, 2021: NGSTrefftz now comes in a docker and with binder notebooks! 

## Publications using the code
* *A unified framework for Trefftz-like discretization methods*  
Philip L. Lederer, Christoph Lehrenfeld, Paul Stocker, Igor Voulis  
[![arXiv](https://img.shields.io/badge/arXiv-2412.00806-b31b1b.svg)](https://arxiv.org/abs/2412.00806)
* *Inf-sup stable space-time Local Discontinuous Galerkin method for the heat equation*  
Sergio Gómez, Chiara Perinati, Paul Stocker  
[![arXiv](https://img.shields.io/badge/arXiv-2411.14819-b31b1b.svg)](https://arxiv.org/abs/2411.14819)
* *Polynomial quasi-Trefftz DG for PDEs with smooth coefficients: elliptic problems*  
Lise-Marie Imbert-Gérard, Andrea Moiola, Chiara Perinati, Paul Stocker  
[![arXiv](https://img.shields.io/badge/arXiv-2408.00392-b31b1b.svg)](https://arxiv.org/abs/2408.00392)
* *Trefftz Discontinuous Galerkin discretization for the Stokes problem*  
Philip L. Lederer, Christoph Lehrenfeld, Paul Stocker  
[![arXiv](https://img.shields.io/badge/arXiv-2306.14600-b31b1b.svg)](https://arxiv.org/abs/2306.14600)
* *Unfitted Trefftz discontinuous Galerkin methods for elliptic boundary value problems*  
Fabian Heimann, Christoph Lehrenfeld, Paul Stocker, Henry von Wahl  
[![arXiv](https://img.shields.io/badge/arXiv-2212.12236-b31b1b.svg)](https://arxiv.org/abs/2212.12236)
* *Embedded Trefftz discontinuous Galerkin methods*  
Christoph Lehrenfeld, Paul Stocker   
[![arXiv](https://img.shields.io/badge/arXiv-2201.07041-b31b1b.svg)](https://arxiv.org/abs/2201.07041)
* *A space-time quasi-Trefftz DG method for the wave equation with piecewise-smooth coefficients*   
Lise-Marie Imbert-Gérard, Andrea Moiola, Paul Stocker  
[![arXiv](https://img.shields.io/badge/arXiv-2011.04617-b31b1b.svg)](https://arxiv.org/abs/2011.04617)
* *Tent pitching and Trefftz-DG method for the acoustic wave equation*  
Ilaria Perugia, Joachim Schöberl, Paul Stocker, Christoph Wintersteiger   
[![arXiv](https://img.shields.io/badge/arXiv-1907.02367-b31b1b.svg)](https://arxiv.org/abs/1907.02367)
* *Space-time Trefftz DG methods for parabolic PDEs*  
Constanze Heil, [Master's thesis](https://doi.org/10.25625/ZSA8UU/2L4C1E)
* *Embedded Trefftz Trace DG Methods for PDEs on unfitted Surfaces*  
Erik Schlesinger, [Master's thesis](https://doi.org/10.25625/QTOPWD/93ZYRQ)


If you are using `ngstrefftz` in your academic work, please consider citing 
```
Stocker, P., (2022). NGSTrefftz: Add-on to NGSolve for Trefftz methods. 
Journal of Open Source Software, 7(71), 4135, https://doi.org/10.21105/joss.04135
```
