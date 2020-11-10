# trefftzngs
I am working on a DG tent-pitching extension for NGSolve using Trefftz basis functions. The code is under constant development. 

depends on [ngsolve](https://github.com/NGSolve/ngsolve)
and requires to add to [intrule.cpp](https://github.com/NGSolve/ngsolve/blob/master/fem/intrule.cpp)
```
template class SIMD_MappedIntegrationRule<3,4>;
```
