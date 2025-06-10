# import sys, os
# sys.path.append(os.path.join(os.path.dirname(sys.path[0])))
from ngstrefftz import *
# from ngstents import TentSlab
from netgen.geom2d import unit_square
from netgen.csg import unit_cube
from ngsolve.TensorProductTools import *
from ngsolve import *
from dg import *
import time
ngsglobals.msg_level=0
SetNumThreads(4)


import netgen.meshing as ngm
def CartSquare(N,t_steps,xshift=0,bndc="dirichlet"):
	ngmesh = ngm.Mesh()
	ngmesh.SetGeometry(unit_square)
	ngmesh.dim = 2
	pnums = []
	for j in range(t_steps + 1):
		for i in range(N + 1):
			pnums.append(ngmesh.Add(ngm.MeshPoint(ngm.Pnt(i / N + xshift, j / t_steps, 0))))

	foo = ngm.FaceDescriptor(surfnr=1,domin=1,bc=1)
	ngmesh.Add (foo)
	ngmesh.SetMaterial(1, "mat")
	for j in range(t_steps):
		for i in range(N):
			ngmesh.Add(ngm.Element2D(1, [pnums[i + j * (N + 1)],
										pnums[i + 1 + j * (N + 1)],
										pnums[i + 1 + (j + 1) * (N + 1)],
										pnums[i + (j + 1) * (N + 1)]]))

	fde = ngm.FaceDescriptor(surfnr=1,domin=1,bc=1)
	fde.bcname = "inflow"
	fdid = ngmesh.Add(fde)
	for i in range(N):
		ngmesh.Add(ngm.Element1D([pnums[i], pnums[i + 1]], index=1))

	fde = ngm.FaceDescriptor(surfnr=2,domin=1,bc=2)
	fde.bcname = "outflow"
	fdid = ngmesh.Add(fde)
	for i in range(N):
		ngmesh.Add(ngm.Element1D([pnums[i + t_steps * (N + 1)], pnums[i + 1 + t_steps * (N + 1)]], index=2))

	fde = ngm.FaceDescriptor(surfnr=3,domin=1,bc=3)
	fde.bcname = bndc
	fdid = ngmesh.Add(fde)
	for i in range(t_steps):
		ngmesh.Add(ngm.Element1D([pnums[N + i * (N + 1)], pnums[N + (i + 1) * (N + 1)]], index=3))
		ngmesh.Add(ngm.Element1D([pnums[0 + i * (N + 1)], pnums[0 + (i + 1) * (N + 1)]], index=3))


	ngmesh.SetBCName(0,"inflow")
	ngmesh.SetBCName(1,"outflow")
	ngmesh.SetBCName(2,bndc)

	mesh = Mesh(ngmesh)
	# print("boundaries" + str(mesh.GetBoundaries()))
	return mesh


########################################################################
# Laplace
########################################################################
def testlaptrefftz(order,mesh):
    """
    >>> order = 5
    >>> mesh = Mesh(unit_square.GenerateMesh(maxh=0.3))
    >>> testlaptrefftz(order,mesh) # doctest:+ELLIPSIS
    8...e-09
    >>> mesh = Mesh(unit_cube.GenerateMesh(maxh = 1))
    >>> testlaptrefftz(order,mesh) # doctest:+ELLIPSIS
    2...e-06
    """
    fes = trefftzfespace(mesh,order=order,eq="laplace")
    a,f = dgell(fes,exactlap)
    gfu = GridFunction(fes)
    gfu.vec.data = a.mat.Inverse() * f.vec
    return sqrt(Integrate((gfu-exactlap)**2, mesh))

########################################################################
# Helmholtz
########################################################################

def testhelmtrefftz(order,mesh):
    """
    >>> order = 5
    >>> mesh = Mesh(unit_square.GenerateMesh(maxh=0.4))
    >>> [testhelmtrefftz(order,mesh)] # doctest:+ELLIPSIS
    [...e-09]
    """
    omega=1
    exact = exp(1j*sqrt(0.5)*(x+y))
    gradexact = CoefficientFunction((sqrt(0.5)*1j*exact, sqrt(0.5)*1j*exact))
    n = specialcf.normal(mesh.dim)
    bndc = CoefficientFunction((sqrt(0.5)*1j*exact, sqrt(0.5)*1j*exact))*n + 1j*omega*exact

    fes = trefftzfespace(mesh,order=order,eq="helmholtz",complex=True,dgjumps=True)
    fes2 = trefftzfespace(mesh,order=order,eq="helmholtzconj",complex=True,dgjumps=True)
    a,f = dghelm(fes,fes2,bndc,omega)
    gfu = GridFunction(fes)
    with TaskManager():
        gfu.vec.data = a.mat.Inverse() * f.vec
    terror = sqrt(Integrate((gfu-exact)*Conj(gfu-exact), mesh).real)
    return terror


########################################################################
# Waveeq
########################################################################

def Cartsolve2D(fes,c,fullsys=False):
    """
    We can solve on a simple rectangle grid
    >>> N = 4
    >>> c = 2
    >>> mesh = CartSquare(N,c*N)

    using Trefftz basis
    >>> fes = trefftzfespace(mesh, order = 8, dgjumps=True, eq="wave")
    >>> fes.SetCoeff(c)
    >>> Cartsolve2D(fes,c) # doctest:+ELLIPSIS
    [17.0, ..., ...e-09, ...e-07]

    or normal L2 basis, requiring the full system
    >>> fes = monomialfespace(mesh, order=8, dgjumps=True)
    >>> Cartsolve2D(fes,c,True) # doctest:+ELLIPSIS
    [45.0, ..., ...e-09, ...e-07]
    """
    k = 3
    sol = sin( k*(c*y + x) )
    v0 = c*k*cos(k*(c*y+x))
    sig0 = -k*cos(k*(c*y+x))
    gD = v0
    U0=sol

    start = time.time()
    [a,f] = DGwaveeqsys(fes,U0,v0,sig0,c,gD,fullsys,False,0.5,0.5,1)
    # print("DGsys: ", str(time.clock()-start))

    start = time.time()
    gfu = GridFunction(fes, name="uDG")
    gfu.vec.data = a.mat.Inverse()*f.vec
    cond = 0
    # print("DGsolve: ", str(time.clock()-start))

    mesh = fes.mesh
    gts = CF((sol.Diff(x),sol.Diff(y)))
    L2error = sqrt(Integrate((sol - gfu)*(sol - gfu), mesh))
    sH1error = sqrt(Integrate((gts - grad(gfu))*(gts - grad(gfu)), mesh))

    dof=fes.ndof/fes.mesh.ne

    return [dof,cond,L2error,sH1error]



def TestQTrefftz(order, mesh, t_step,qtrefftz=1):
    """
    Solve with quasi-Trefftz basis functions
    >>> order = 4
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

    fes = trefftzfespace(mesh, order=order, dgjumps=True, eq="qtwave")
    fes.SetCoeff(wavespeed)
    [a,f] = DGwaveeqsys(fes,U0,v0,sig0,wavespeed,gD,True,False,alpha=0.5,beta=0.5,gamma=1,mu=0.5)
    gfu = GridFunction(fes, name="uDG")
    gfu.vec.data = a.mat.Inverse()*f.vec
    dgerror = DGnormerror(fes,gfu,bdd[1:3],wavespeed,alpha=0.5,beta=0.5)

    return dgerror


def TestBessel(order, mesh, t_step):
    """
    Solve using quasi-Trefftz basis functions
    >>> order = 4
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

    fes = trefftzfespace(mesh, order=order, dgjumps=True, eq="qtwave")
    fes.SetCoeff(wavespeed,BB)
    [a,f] = DGwaveeqsys(fes,U0,v0,sig0,wavespeed,gD,True,False,alpha=0,beta=0,gamma=1,mu=0,BB=BB)
    gfu = GridFunction(fes, name="uDG")
    gfu.vec.data = a.mat.Inverse()*f.vec
    dgerror = DGnormerror(fes,gfu,bdd[1:3],wavespeed,alpha=0,beta=0,BB=BB)

    return dgerror



########################################################################
# Heateq
########################################################################
import ngsolve.meshes as ngsm
def TestHeat():
    """
    Solve heat equation using Trefftz fcts
    >>> TestHeat() # doctest:+ELLIPSIS
    0.0003...
    """
    mesh = ngsm.MakeStructured2DMesh(nx=32, ny=32, periodic_x=False)
    diffusion = 10

    order = 4
    #fes = L2(mesh, order=order, dgjumps=True)
    fes = trefftzfespace(mesh, order=2*order+1, dgjumps=True, eq="heat")
    fes.SetCoeff(diffusion)
    ubnd = exp(y*diffusion)*exp(x)
    gfu = dgheat(fes, diffusion, ubnd)

    return sqrt(Integrate((gfu-ubnd)**2, mesh))


########################################################################
# QTelliptic
########################################################################

def qtell(mesh,order):

    """
    Solve qtell equation using Trefftz fcts
    >>> mesh = Mesh(unit_square.GenerateMesh(maxh=0.3))
    >>> for i in range(3):
    ...     qtell(mesh,4) # doctest:+ELLIPSIS
    ...     mesh.Refine()
    3...e-05
    ...e-06
    ...e-08
    """
    exact = (1+x+y)**5
    A = CF((1+x+y,0,0,1+x+y),dims=(2,2))
    B = CF((1,0))
    C = 3/(1+x+y)
    rhs = -sum( (A*CF((exact.Diff(x),exact.Diff(y))))[i].Diff(var) for i,var in enumerate([x,y])) + B*CF((exact.Diff(x),exact.Diff(y))) + C*exact
    rhs.Compile()
    Dbndc = exact
    Dbnd=".*"
    Nbnd=""
    Nbndc = 0

    fes = trefftzfespace(mesh,order=order,eq="qtelliptic")
    with TaskManager():
        fes.SetCoeff(A,B,C)
        uf = fes.GetParticularSolution(rhs)

    a,f = dgell(fes,Dbndc,rhs=rhs,uf=uf,A=A,B=B,C=C)
    with TaskManager():
        gfu = GridFunction(fes)
        gfu.vec.data = a.mat.Inverse() * f.vec
    gfu += uf

    error = sqrt(Integrate((gfu-exact)**2,mesh))
    return error


if __name__ == "__main__":
    import doctest
    doctest.testmod()
