from ngstrefftz import *
# from ngstents import TentSlab
from netgen.geom2d import unit_square
from netgen.csg import unit_cube
from ngsolve.TensorProductTools import *
from ngsolve import *
from svdt import *
import time


########################################################################
# Laplace
########################################################################
def testlaptrefftz(order,mesh):
    """
    >>> order = 5
    >>> mesh = Mesh(unit_square.GenerateMesh(maxh=0.3))
    >>> testlaptrefftz(order,mesh) # doctest:+ELLIPSIS
    1...e-08
    >>> mesh = Mesh(unit_cube.GenerateMesh(maxh = 1))
    >>> testlaptrefftz(order,mesh) # doctest:+ELLIPSIS
    2...e-06
    """
    fes = FESpace("trefftzfespace",mesh,order=order,eq="laplace")
    a,f = dglap(fes,exactlap)
    gfu = GridFunction(fes)
    gfu.vec.data = a.mat.Inverse() * f.vec
    return sqrt(Integrate((gfu-exactlap)**2, mesh))

########################################################################
# Helmholtz
########################################################################
def dghelm(fes,fes2,bndc,omega):
    mesh = fes.mesh
    order = fes.globalorder
    n = specialcf.normal(mesh.dim)
    h = specialcf.mesh_size
    alpha = 1/(omega*h)
    beta = omega*h
    delta = omega*h

    u = fes.TrialFunction()
    v = fes.TestFunction()
    if fes2 is not None:
        v = fes2.TestFunction()
    jump_u = (u-u.Other())*n
    jump_v = (v-v.Other())*n
    jump_du = (grad(u)-grad(u.Other()))*n
    jump_dv = (grad(v)-grad(v.Other()))*n
    mean_u = 0.5 * ((u)+(u.Other()))
    mean_du = 0.5 * (grad(u)+grad(u.Other()))
    mean_dv = 0.5 * (grad(v)+grad(v.Other()))

    a = BilinearForm(fes)
    if fes2 is not None:
        a = BilinearForm(fes,fes2)
    a += mean_u*(jump_dv) * dx(skeleton=True)
    a += 1/omega*1j*beta*jump_du*(jump_dv) * dx(skeleton=True)
    a += -mean_du*(jump_v) * dx(skeleton=True)
    a += omega*1j*alpha*jump_u*(jump_v) * dx(skeleton=True)

    a += (1-delta)*u*(grad(v))*n * ds(skeleton=True)
    a += 1/omega*1j*delta*(grad(u)*n)*((grad(v))*n) * ds(skeleton=True)
    a += -delta*grad(u)*n*(v) * ds(skeleton=True)
    a += omega*1j*(1-delta)*u*(v) * ds(skeleton=True)

    f = LinearForm(fes)
    if fes2 is not None:
        f = LinearForm(fes2)
    f += 1/omega*1j*delta*bndc*(grad(v))*n*ds(skeleton=True)
    f += (1-delta)*bndc*(v)*ds(skeleton=True)

    with TaskManager():
        a.Assemble()
        f.Assemble()
    return a,f

def testhelmtrefftz(order,mesh):
    """
    >>> order = 5
    >>> mesh = Mesh(unit_square.GenerateMesh(maxh=0.3))
    >>> testhelmtrefftz(order,mesh) # doctest:+ELLIPSIS
    6...e-10
    """
    omega=1
    exact = exp(1j*sqrt(0.5)*(x+y))
    gradexact = CoefficientFunction((sqrt(0.5)*1j*exact, sqrt(0.5)*1j*exact))
    n = specialcf.normal(mesh.dim)
    bndc = CoefficientFunction((sqrt(0.5)*1j*exact, sqrt(0.5)*1j*exact))*n + 1j*omega*exact

    fes = FESpace("trefftzfespace",mesh,order=order,eq="laplace")
    fes = trefftzfespace(mesh,order=order,eq="helmholtz",complex=True,dgjumps=True)
    fes2 = trefftzfespace(mesh,order=order,eq="helmholtzconj",complex=True,dgjumps=True)
    a,f = dghelm(fes,fes2,bndc,omega)
    gfu = GridFunction(fes)
    gfu.Set(exact)
    with TaskManager():
        gfu.vec.data = a.mat.Inverse() * f.vec
    terror = sqrt(Integrate((gfu-exact)*Conj(gfu-exact), mesh).real)
    return terror


########################################################################
# Waveeq
########################################################################
def TestSolution2D(fes,c,timeoffset=0):
    k = 3
    truesol = sin( k*(c*y + x) )
    v0 = c*k*cos(k*(c*y+x))
    sig0 = -k*cos(k*(c*y+x))
    gD = v0
    U0 = GridFunction(fes)
    U0.Set(truesol)
    return [truesol,U0,sig0,v0,gD]


def PostProcess(fes, truesol, sol):
    mesh = fes.mesh
    U = GridFunction(fes)
    U.Set(truesol)
    L2error = sqrt(Integrate((truesol - sol)*(truesol - sol), mesh))
    sH1error = sqrt(Integrate((grad(U) - grad(sol))*(grad(U) - grad(sol)), mesh))
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
    [81.0, ..., ...e-12, ...e-10]
    """
    if inputsol is None:
        inputsol = TestSolution2D(fes,c)
    [truesol,U0,sig0,v0,gD] = inputsol

    start = time.time()
    [a,f] = DGwaveeqsys(fes,U0,v0,sig0,c,gD,fullsys,False,0.5,0.5,1)
    # print("DGsys: ", str(time.clock()-start))

    start = time.time()
    gfu = DGsolve(fes,a,f)
    cond = 0
    # print("DGsolve: ", str(time.clock()-start))

    [L2error, sH1error] = PostProcess(fes,truesol,gfu)

    dof=fes.ndof/fes.mesh.ne

    return [dof,cond,L2error,sH1error]



def TestQTrefftz(order, mesh, t_step,qtrefftz=1):
    """
    Solve with quasi-Trefftz basis functions
    >>> order = 4
    >>> SetNumThreads(1)
    >>> t_step = 1
    >>> for h in [4,8,16,32]:
    ...        mesh = CartSquare(h,h)
    ...        TestQTrefftz(order,mesh,t_step) # doctest:+ELLIPSIS
    0.001...
    0.0001...
    ...e-05
    ...e-06
    """
    ca=2.5
    bdd = CoefficientFunction((
            (x+1)**ca * exp(-sqrt(ca*(ca-1))*y),
            ca*(x+1)**(ca-1) * exp(-sqrt(ca*(ca-1))*y),
            -sqrt(ca*(ca-1)) * (x+1)**ca * exp(-sqrt(ca*(ca-1))*y)
        ))
    wavespeed=CoefficientFunction((x+1))

    U0=bdd[0]
    gD=bdd[2]
    v0=bdd[2]
    sig0=-bdd[1]

    fes = trefftzfespace(mesh, order=order, dgjumps=True, basistype=0, useshift=True, eq="qtwave")
    fes.SetWavespeed(wavespeed)
    [a,f] = DGwaveeqsys(fes,U0,v0,sig0,wavespeed,gD,True,False,alpha=0.5,beta=0.5,gamma=1,mu=0.5)
    gfu = DGsolve(fes,a,f)
    dgerror = DGnormerror(fes,gfu,bdd[1:3],wavespeed,alpha=0.5,beta=0.5)

    return dgerror


def TestBessel(order, mesh, t_step):
    """
    Solve using quasi-Trefftz basis functions
    >>> order = 4
    >>> SetNumThreads(1)
    >>> t_step = 1
    >>> for h in [4,8,16,32]:
    ...        mesh = CartSquare(h,h,xshift=3)
    ...        TestBessel(order,mesh,t_step) # doctest:+ELLIPSIS
    9...e-06
    ...e-07
    ...e-08
    ...e-09
    """

    D = mesh.dim
    t = CoordCF(D)
    t_start = 0

    c=0
    bdd = CoefficientFunction((
        ((x+c)**(-2)*sin(x+c)-(x+c)**(-1)*cos(x+c)) * cos(y),
        (-2*(x+c)**(-3)*sin(x+c) + 2*(x+c)**(-2)*cos(x+c) + (x+c)**(-1)*sin(x+c)) * cos(y),
        -((x+c)**(-2)*sin(x+c)-(x+c)**(-1)*cos(x+c)) * sin(y)
        ))
    wavespeed=1/sqrt((x+c)**2-2)
    BB = (x+c)**2

    U0=bdd[0]
    gD=bdd[2]
    v0=bdd[2]
    sig0=-bdd[1]

    fes = trefftzfespace(mesh, order=order, dgjumps=True, basistype=0, useshift=True,usescale=False, eq="qtwave")
    fes.SetWavespeed(wavespeed,BB)
    [a,f] = DGwaveeqsys(fes,U0,v0,sig0,wavespeed,gD,True,False,alpha=0,beta=0,gamma=1,mu=0,BB=BB)
    gfu = DGsolve(fes,a,f)
    dgerror = DGnormerror(fes,gfu,bdd[1:3],wavespeed,alpha=0,beta=0,BB=BB)

    return dgerror


if __name__ == "__main__":
    import doctest
    doctest.testmod()
