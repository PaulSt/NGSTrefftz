import sys
sys.path.append("..")
from trefftzngs import *
import numpy as np
from ngsolve import *
import netgen.gui
from DGeq import *
import time
from prodmesh import *
from netgen.geom2d import unit_square
from netgen.csg import unit_cube
from ngsolve.TensorProductTools import *
from ngsolve.comp import *


def TestSolution2D(fes,c,timeoffset=0):
    k = 3
    truesol = sin( k*(c*y + x) )
    v0 = c*k*cos(k*(c*y+x))
    sig0 = -k*cos(k*(c*y+x))
    gD = v0
    # truesol = exp(-1.3*k*((x-0.15)-c*y)*((x-0.15)-c*y)) + exp(-k*((x-0.9)+c*y)*((x-0.9)+c*y))
    # truesol = exp(-k*((x-0.5)-c*y)*((x-0.5)-c*y))
    # v0 = grad(U0)[1]
    # sig0 = -grad(U0)[0]
    # gD = 0 #v0
    U0 = GridFunction(fes)
    U0.Set(truesol)
    return [truesol,U0,sig0,v0,gD]


def PostProcess(fes, truesol, sol):
    mesh = fes.mesh
    U = GridFunction(fes)
    U.Set(truesol)
    L2error = sqrt(Integrate((truesol - sol)*(truesol - sol), mesh))
    sH1error = sqrt(Integrate((grad(U) - grad(sol))*(grad(U) - grad(sol)), mesh))
    # print("error=", L2error)
    # print("grad-error=", sH1error)
    Draw(sol,mesh,'sol')
    Draw(grad(sol),mesh,'gradsol')
    return [L2error,sH1error]


def Cartsolve2D(fes,c,fullsys=False,inputsol=None):
    """
    We can solve on a simple rectangle grid
    >>> N = 4
    >>> c = 2
    >>> order = 8
    >>> mesh = CartSquare(N,c*N)

    using Trefftz basis
    >>> fes = FESpace("trefftzfespace", mesh, order = order, wavespeed = c, dgjumps=True, basistype=0)
    >>> Cartsolve2D(fes,c) # doctest:+ELLIPSIS
    [17.0, ..., ...e-09, ...e-08]

    or normal L2 basis, requiring the full system
    >>> fes = L2(mesh, order=order, dgjumps=True)
    >>> Cartsolve2D(fes,c,True) # doctest:+ELLIPSIS
    [81.0, ..., ...e-10, ...e-09]
    """
    if inputsol is None:
        inputsol = TestSolution2D(fes,c)
    [truesol,U0,sig0,v0,gD] = inputsol

    start = time.time()
    [a,f] = DGeqsys(fes,U0,v0,sig0,c,gD,fullsys)
    # print("DGsys: ", str(time.clock()-start))

    start = time.time()
    [gfu,cond] = DGsolve(fes,a,f)
    # print("DGsolve: ", str(time.clock()-start))

    [L2error, sH1error] = PostProcess(fes,truesol,gfu)

    dof=fes.ndof/fes.mesh.ne

    return [dof,cond,L2error,sH1error]


if __name__ == "__main__":
    import doctest
    doctest.testmod()
