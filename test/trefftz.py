# import sys, os
# sys.path.append(os.path.join(os.path.dirname(sys.path[0])))
from ngstrefftz import *
from embt import dglap
# from ngstents import TentSlab
from netgen.geom2d import unit_square
from netgen.csg import unit_cube
from ngsolve.TensorProductTools import *
from ngsolve import *
from embt import *
import time
ngsglobals.msg_level=0


def DGwaveeqsys(fes,U0,v0,sig0,c,gD,fullsys=False, applyrhs = False,alpha=0.5,beta=0.5,gamma=1,mu=0,BB=1):
    D = fes.mesh.dim - 1
    # if not isinstance(BB,(int,float)): BB.spacedim = D
    U = fes.TrialFunction()
    V = fes.TestFunction()
    gU = grad(U)
    gV = grad(V)

    v = gU[D]
    sig = CoefficientFunction(tuple([-gU[i] for i in  range(D)]))*BB
    w = gV[D]
    tau = CoefficientFunction(tuple([-gV[i] for i in  range(D)]))*BB

    vo = gU.Other()[D]
    sigo = CoefficientFunction(tuple([-gU.Other()[i] for i in  range(D)]))*BB
    wo = gV.Other()[D]
    tauo = CoefficientFunction(tuple([-gV.Other()[i] for i in  range(D)]))*BB

    h = specialcf.mesh_size
    n = specialcf.normal(D+1)
    n_t = n[D]/Norm(n)
    n_x = CoefficientFunction( tuple([n[i]/Norm(n) for i in  range(D)]) )

    mean_v = 0.5*(v+vo)
    mean_w = 0.5*(w+wo)
    mean_sig = 0.5*(sig+sigo)
    mean_tau = 0.5*(tau+tauo)

    jump_vx = ( v - vo ) * n_x
    jump_wx = ( w - wo ) * n_x
    jump_sigx = (( sig - sigo ) * n_x)
    jump_taux = (( tau - tauo ) * n_x)

    jump_vt = ( v - vo ) * n_t
    jump_wt = ( w - wo ) * n_t
    jump_sigt = ( sig - sigo ) * n_t
    jump_taut = ( tau - tauo ) * n_t

    jump_Ut = (U - U.Other()) * n_t
    jump_Vt = (V - V.Other()) * n_t

    timelike = n_x*n_x # n_t=0
    spacelike = n_t**2 # n_x=0

    delta=0.5
    theta=1
    gR=0

    a = BilinearForm(fes)
    if(fullsys==True):
        HV = V.Operator("hesse")
        a += SymbolicBFI( - v * (BB*sum([-HV[i*(D+2)] for i in range(D)]) + pow(c,-2)*HV[(D+1)*(D+1)-1]) )
        if not isinstance(BB,(int,float)):
            a += SymbolicBFI( - v * (BB.Diff(x)*(-gV[0]) + BB.Diff(y)*(-gV[1])) )
        # a += SymbolicBFI(  sig*(CoefficientFunction([-HV[i,D] for i in range(D)]) + CoefficientFunction([HV[D,i] for i in range(D)])) )
        HU = U.Operator("hesse")
        a += SymbolicBFI(  mu * pow(c,2)
                              * (-sum([HU[i*(D+2)] for i in range(D)]) + pow(c,-2)*HU[(D+1)*(D+1)-1])
                              * (-sum([HV[i*(D+2)] for i in range(D)]) + pow(c,-2)*HV[(D+1)*(D+1)-1])
                        )
    #space like faces, w/o x jump ASSUME TENSOR MESH
    if(applyrhs == False):
        a += SymbolicBFI( spacelike * ( pow(c,-2)*IfPos(n_t,v,vo)*jump_wt + IfPos(n_t,sig,sigo)*jump_taut/BB ), VOL, skeleton=True )
    # a += SymbolicBFI( spacelike * ( IfPos(n_t,v,vo)*(pow(c,-2)*jump_wt+jump_taux) + IfPos(n_t,sig,sigo)*(jump_wx+jump_taut) ), VOL, skeleton=True )
    #time like faces
    a += SymbolicBFI( timelike * ( mean_v*jump_taux + mean_sig*jump_wx + alpha*jump_vx*jump_wx + beta*jump_sigx*jump_taux ), VOL, skeleton=True )        #t=T (or *x)
    a += SymbolicBFI( ( pow(c,-2)*v*w + sig*tau/BB ), BND, definedon=fes.mesh.Boundaries("outflow"), skeleton=True)
    #dirichlet boundary 'timelike'
    a += SymbolicBFI( ( sig*n_x*w + alpha*v*w ), BND, definedon=fes.mesh.Boundaries("dirichlet"), skeleton=True)
    #impedence boundary
    a += SymbolicBFI( (1-delta)*theta/c*v*w+(1-delta)*v*(tau*n_x)+delta*(sig*n_x)*w+delta*c/theta*(sig*n_x)*(tau*n_x) , BND, definedon=fes.mesh.Boundaries("robin"), skeleton=True)
    #correction term to recover sol of second order system
    if(applyrhs == False):
        a += SymbolicBFI( spacelike * ( gamma * jump_Ut*jump_Vt ), VOL, skeleton=True )
    # a += SymbolicBFI( spacelike * ( gamma * (n_t*jump_Ut)*IfPos(n_t,V.Other(),V) ), VOL, skeleton=True )
    # a += SymbolicBFI( spacelike * ( gamma * (-n_t*jump_Ut)*IfPos(n_t,V.Other(),V) ), VOL, skeleton=True )
    # a += SymbolicBFI( spacelike * ( gamma * (-jump_Ut)*IfPos(n_t,V.Other(),V) ), VOL, skeleton=True )
    a += SymbolicBFI( ( gamma * U*V ), BND, definedon=fes.mesh.Boundaries("inflow"), skeleton=True )
    a.Assemble()

    if(applyrhs == False):

        f = LinearForm(fes)
        f += SymbolicLFI( ( pow(c,-2)*v0*w + sig0*tau ), BND, definedon=fes.mesh.Boundaries("inflow"), skeleton=True) #t=0 (or *(1-x))
        f += SymbolicLFI( ( gD * (alpha*w - tau*n_x) ), BND, definedon=fes.mesh.Boundaries("dirichlet"), skeleton=True) #dirichlet boundary 'timelike'
        f += SymbolicLFI( gamma * ( (U0)*V ), BND, definedon=fes.mesh.Boundaries("inflow"),  skeleton=True ) #rhs correction term to recover sol of second order system
        f += SymbolicLFI( ( gR * ((1-delta)*w-delta*c/theta*tau*n_x) ), BND, definedon=fes.mesh.Boundaries("robin"), skeleton=True) #robin boundary 'timelike'
        f.Assemble()
    else:
        # f = GridFunction(fes)
        # rhs = BilinearForm(fes)
        # rhs += SymbolicBFI( ( pow(c,-2)*vo*w + sigo*tau ), BND, definedon=fes.mesh.Boundaries("inflow"), skeleton=True) #t=0 (or *(1-x))
        # rhs += SymbolicBFI( ( gD * (alpha*w - tau*n_x) ), BND, definedon=fes.mesh.Boundaries("dirichlet"), skeleton=True) #dirichlet boundary 'timelike'
        # rhs += SymbolicBFI( gamma * ( U.Other()*V ), BND, definedon=fes.mesh.Boundaries("inflow"),  skeleton=True ) #rhs correction term to recover sol of second order system
        # rhs.Apply(U0.vec,f.vec)

        # a = BilinearForm(fes)
        # if(fullsys==True):
            # HV = V.Operator("hesse")
            # a += SymbolicBFI(  -v*(-sum([HV[i*(D+2)] for i in range(D)]) + pow(c,-2)*HV[(D+1)*(D+1)-1]) )
        # a += SymbolicBFI( timelike * ( mean_v*jump_taux + mean_sig*jump_wx + alpha*jump_vx*jump_wx + beta*jump_sigx*jump_taux ), VOL, skeleton=True ) #time like faces
        # a += SymbolicBFI( IfPos(n_t,1,0) * spacelike * ( pow(c,-2)*v*w + sig*tau ), element_boundary=True) #t=T (or *x)
        # a += SymbolicBFI( IfPos(n_t,0,1) * spacelike * ( gamma * U*V ), element_boundary=True ) #BND correction term to recover sol of second order system
        # a += SymbolicBFI( ( sig*n_x*w + alpha*v*w ), BND, definedon=fes.mesh.Boundaries("dirichlet"), skeleton=True) #dirichlet boundary 'timelike'
        # a += SymbolicBFI( (1-delta)*theta/c*v*w+(1-delta)*v*(tau*n_x)+delta*(sig*n_x)*w+delta*c/theta*(sig*n_x)*(tau*n_x) , BND, definedon=fes.mesh.Boundaries("robin"), skeleton=True)
        # a.Assemble()
        f = GridFunction(fes)
        rhs = BilinearForm(fes)
        rhs += SymbolicBFI( IfPos(n_t,0,1) * spacelike * ( pow(c,-2)*vo*w + sigo*tau ), element_boundary=True ) #space like faces, w/o x jump
        rhs += SymbolicBFI( IfPos(n_t,0,1) * spacelike * ( U.Other()*V ), element_boundary=True )
        rhs += SymbolicBFI( ( gD * (alpha*w - tau*n_x) ), BND, definedon=fes.mesh.Boundaries("dirichlet"), skeleton=True) #dirichlet boundary 'timelike'
        rhs += SymbolicBFI( ( gR * ((1-delta)*w-delta*c/theta*tau*n_x) ), BND, definedon=fes.mesh.Boundaries("robin"), skeleton=True) #robin boundary 'timelike'
        rhs.Apply(U0.vec,f.vec)

    return [a,f]

def DGnormerror(fes,uh,gradtruesol,c,alpha,beta,BB=1):
    D = fes.mesh.dim - 1
    U = fes.TrialFunction()
    V = fes.TestFunction()
    gU = grad(U)
    gV = grad(V)

    v = gU[D]
    sig = CoefficientFunction(tuple([-gU[i] for i in  range(D)]))
    w = gV[D]
    tau = CoefficientFunction(tuple([-gV[i] for i in  range(D)]))

    vo = gU.Other()[D]
    sigo = CoefficientFunction(tuple([-gU.Other()[i] for i in  range(D)]))
    wo = gV.Other()[D]
    tauo = CoefficientFunction(tuple([-gV.Other()[i] for i in  range(D)]))

    h = specialcf.mesh_size
    n = specialcf.normal(D+1)
    n_t = n[D]/Norm(n)
    n_x = CoefficientFunction( tuple([n[i]/Norm(n) for i in  range(D)]) )

    mean_v = 0.5*(v+vo)
    mean_w = 0.5*(w+wo)
    mean_sig = 0.5*(sig+sigo)
    mean_tau = 0.5*(tau+tauo)

    jump_vx = ( v - vo ) * n_x
    jump_wx = ( w - wo ) * n_x
    jump_sigx = (( sig - sigo ) * n_x)
    jump_taux = (( tau - tauo ) * n_x)

    jump_vt = ( v - vo ) * n_t
    jump_wt = ( w - wo ) * n_t
    jump_sigt = ( sig - sigo ) * n_t
    jump_taut = ( tau - tauo ) * n_t

    jump_Ut = (U - U.Other()) * n_t
    jump_Vt = (V - V.Other()) * n_t

    timelike = n_x*n_x # n_t=0
    spacelike = n_t**2 # n_x=0

    # make use of the fact that on a Cart mesh we have A(w,tau;w,tau)=||(w,tau)||^2_DG
    a = BilinearForm(fes)
    # if(fullsys==True):
        # HV = V.Operator("hesse")
        # # a += SymbolicBFI(  -v*(-HV[0]+pow(c,-2)*HV[3]) ) #- sig*(-HV[1]+HV[2])  )
        # a += SymbolicBFI(  -v*(-sum([HV[i*(D+2)] for i in range(D)]) + pow(c,-2)*HV[(D+1)*(D+1)-1]) )
        # HU = U.Operator("hesse")
        # a += SymbolicBFI(  mu * pow(c,2)
                              # * (-sum([HU[i*(D+2)] for i in range(D)]) + pow(c,-2)*HU[(D+1)*(D+1)-1])
                              # * (-sum([HV[i*(D+2)] for i in range(D)]) + pow(c,-2)*HV[(D+1)*(D+1)-1])
                        # )
    # space like faces, w/o x jump
    a += SymbolicBFI( spacelike * 0.5 * pow(c,-2)*jump_wt*jump_vt , VOL, skeleton=True)
    a += SymbolicBFI( spacelike * 0.5 * 1/BB*jump_taut * jump_sigt , VOL, skeleton=True )
    # time like faces
    a += SymbolicBFI( timelike * alpha * jump_vx * jump_wx , VOL, skeleton=True )
    a += SymbolicBFI( timelike * beta * jump_sigx * jump_taux , VOL, skeleton=True )
    jumppart = uh.vec.CreateVector()
    a.Apply(uh.vec,jumppart)
    # smooth solution vanishes on jumppart
    norm = 0
    norm += uh.vec.InnerProduct(jumppart)
    norm += 0.5 * Integrate((BoundaryFromVolumeCF(grad(uh)[D]) - gradtruesol[D])**2 / c, fes.mesh, definedon=fes.mesh.Boundaries("outflow|inflow"))
    norm += 0.5 * Integrate(sum((BoundaryFromVolumeCF(grad(uh)[i]) - gradtruesol[i])**2 / sqrt(BB) for i in range(D)), fes.mesh, definedon=fes.mesh.Boundaries("outflow|inflow"))
    norm += Integrate(alpha * (BoundaryFromVolumeCF(grad(uh)[D]) - gradtruesol[D])**2 , fes.mesh, definedon=fes.mesh.Boundaries("dirichlet"))

    return sqrt(norm)



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
    >>> fes = trefftzfespace(mesh, order = order, dgjumps=True)
    >>> fes.SetCoeff(c)
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
    gfu = GridFunction(fes, name="uDG")
    gfu.vec.data = a.mat.Inverse()*f.vec
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
    # import pdb; pdb.set_trace()
    u, v = fes.TnT()

    eps = diffusion
    b = CoefficientFunction((0, 1))
    ubnd = exp(y*diffusion)*exp(x)
    lambd = 10
    h = specialcf.mesh_size
    n = specialcf.normal(mesh.dim)

    space_jump_u = (u-u.Other())
    space_jump_v = (v-v.Other())
    space_mean_dudn = 0.5* (grad(u)[0]+grad(u.Other())[0])*n[0]
    space_mean_dvdn = 0.5* (grad(v)[0]+grad(v.Other())[0])*n[0]

    space_diffusion = grad(u)[0]*grad(v)[0] * dx \
        + lambd*(order*n[0])**2/h*space_jump_u*space_jump_v*dx(skeleton=True) \
        + (-space_mean_dudn*space_jump_v - space_mean_dvdn*space_jump_u)*dx(skeleton=True) \
        + (-grad(u)[0]*n[0] * v -grad(v)[0]*n[0] * u)*ds(definedon=mesh.Boundaries("left|right"),skeleton=True) \
        + lambd*(order*n[0])**2/h*u*v*ds(definedon=mesh.Boundaries("left|right"),skeleton=True)

    space_time_convection = -b * u * grad(v)*dx \
        + b*n*IfPos(b*n, u, u.Other()) * (v-v.Other()) * dx(skeleton=True) \
        + b*n*u*v * ds(definedon=mesh.Boundaries("top"), skeleton=True)

    a = BilinearForm(fes, symmetric=False)
    a += eps * space_diffusion + space_time_convection
    a.Assemble()

    #f_coef = (ubnd.Diff(y)-eps*ubnd.Diff(x).Diff(x) )

    f = LinearForm(fes)
    #f += f_coef * v * dx
    f += -b*n* ubnd * v * ds(definedon=mesh.Boundaries("bottom"),skeleton=True)
    f += -eps*grad(v)[0]*n[0] * ubnd * ds(definedon=mesh.Boundaries("left|right"),skeleton=True)
    f += eps*lambd*(order*n[0])**2/h*ubnd*v*ds(definedon=mesh.Boundaries("left|right"),skeleton=True)
    f.Assemble()

    gfu = GridFunction(fes)
    gfu.vec.data = a.mat.Inverse() * f.vec

    return sqrt(Integrate((gfu-ubnd)**2, mesh))



if __name__ == "__main__":
    import doctest
    doctest.testmod()
