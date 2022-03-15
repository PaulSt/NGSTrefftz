---
title: '`NGSTrefftz`: Add-on to NGSolve for Trefftz methods'
tags:
  - numerical methods for PDEs
  - finite elements 
  - Trefftz methods
authors:
  - name: Paul Stocker^[corresponding author] 
    orcid: 0000-0001-5073-3366 
    affiliation: "1" # (Multiple affiliations must be quoted)
affiliations:
 - name: Georg-August-Universit&#0228;t, G&#0246;ttingen, Germany
   index: 1
date: 29 January 2021
bibliography: paper.bib

---

# Summary
`NGSTrefftz` is an add-on to [`Netgen/NGSolve`](https://ngsolve.org/), a finite element software for the numerical treatment of partial differential equations (PDEs).
The package implements Trefftz based discontinuous Galerkin (DG) methods in `NGSolve`.
Trefftz methods reduce the number of unknowns in the discretization of PDE problems by injecting knowledge of the PDE into the approximation functions.
Like `NGSolve`, `NGSTrefftz` is written in C++ and integrates seamlessly with the easy-to-use Python interface of `NGSolve`.

Trefftz methods originate from @trefftz1926 and have since been developed for a wide range of problems, for an overview see @cbe1997trefftz, @TrefftzSurvey, @Qin05, @LLHC08, @kk95.
The central principle of Trefftz methods is the construction of a discrete basis of solutions to the differential operator under consideration, making the space of Trefftz functions problem dependent. 
In combination with finite elements the Trefftz basis is constructed locally, on each mesh element, continuity and boundary conditions are then enforced in the variational formulation.

# Statement of need

`NGSTrefftz` provides a framework to implement Trefftz finite element spaces for `NGSolve`, with spaces for Laplace equation, Helmholtz equation, and acoustic wave equation already implemented. 
`NGSolve` provides a flexible framework to implement variational formulations for arbitrary physical models. 
The Trefftz finite element spaces can then be used with the (bi-)linear forms generated in `NGSolve`.
The focus lies on the combination of Trefftz functions with DG methods or least-squares formulations,
this approach is also often referred to as frameless Trefftz elements.

On top of that, the package provides unique features that are:

* A quasi-Trefftz space that provides Trefftz-like properties for the acoustic wave equation with smooth coefficients. The method and results are presented in @qtrefftz.

* A space--time Trefftz method for the acoustic wave equation on tent-pitched meshes. 
Tent-pitching is a space--time meshing strategy that provides mesh elements that conform to the causality constraint of a hyperbolic system. 
The meshes are generated using [`ngstents`](https://github.com/jayggg/ngstents) and can be used with the Trefftz and quasi-Trefftz space for the acoustic wave equation, results are shown in @StockerSchoeberl and @qtrefftz.

* A general framework to produce Trefftz spaces implicitly is provided by an implementation of the embedded Trefftz method, see @embtrefftz.
The approach produces a Galerkin projection of an underlying discontinuous Galerkin method onto a subspace of Trefftz-type. 
It can be applied to very general cases, including inhomogeneous sources and non-constant coefficient differential operators.

The aim of this package is to facilitate research into Trefftz methods and to make them more accessible to a broad audience.
To the best of our knowledge, the only other open source software package that provides Trefftz finite element methods is [`FreeHyTe`](https://www.sites.google.com/site/ionutdmoldovan/freehyte) [@freehyte], implemented in `MATLAB`. 
Examples for the usage of all the features in `NGSTrefftz` are provided in the form of jupyter notebooks.

# Acknowledgements
P. Stocker has been supported by the Austrian Science Fund (FWF) through the projects F 65, W 1245, and by the German Research Foundation (DFG) through grant 432680300 - SFB 1456.
Special thanks to the `NGSolve` team for their continuous support.

# References
