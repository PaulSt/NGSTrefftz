from ngsolve import *
from ngstrefftz import *
from netgen.geom2d import unit_square
from netgen.csg import unit_cube
import time

Lap = lambda u : sum(Trace(u.Operator('hesse')))

exactlap = exp(x)*sin(y)
exactpoi = sin(x)*sin(y)

eps = 10**-8
mesh2d = Mesh(unit_square.GenerateMesh(maxh=0.3))
mesh3d = Mesh(unit_cube.GenerateMesh(maxh = 1))


########################################################################
# Laplace
########################################################################
def dgell(fes,Dbndc,rhs=0,uf=0,A=1,B=None,C=0,Dbnd=".*",Nbnd="",Nbndc=0):
    """
    >>> fes = L2(mesh2d, order=5,  dgjumps=True)#,all_dofs_together=True)
    >>> a,f = dgell(fes,exactlap)
    >>> gfu = GridFunction(fes)
    >>> gfu.vec.data = a.mat.Inverse() * f.vec
    >>> sqrt(Integrate((gfu-exactlap)**2, mesh2d)) # doctest:+ELLIPSIS
    3...e-09
    """
    if not B: B = CF(tuple([0]*fes.mesh.dim))
    mesh = fes.mesh
    order = fes.globalorder
    n = specialcf.normal(mesh.dim)
    h = specialcf.mesh_size
    alpha = 4*order**2/h
    u = fes.TrialFunction()
    v = fes.TestFunction()

    u = fes.TrialFunction()
    v = fes.TestFunction()
    jump = lambda u: (u-u.Other())*n
    mean_d = lambda u: 0.5 * A * (grad(u)+grad(u).Other())
    mean_B = lambda u: 0.5 * B * (u+u.Other())

    a = BilinearForm(fes)
    a += A*grad(u)*grad(v) * dx \
        +alpha*jump(u)*jump(v) * dx(skeleton=True) \
        +(-mean_d(u)*jump(v)-mean_d(v)*jump(u)) * dx(skeleton=True) \
        +alpha*u*v * ds(skeleton=True,definedon=mesh.Boundaries(Dbnd)) \
        +(-A*grad(u)*n*v-A*grad(v)*n*u)* ds(skeleton=True,definedon=mesh.Boundaries(Dbnd))
    a += (-B*u*grad(v) + C*u*v) * dx \
        + (mean_B(u) * jump(v) + 0.5*sqrt((B*n)**2)*jump(u)*jump(v)) * dx(skeleton=True) \
        + B*u*n*v * ds(skeleton=True,definedon=mesh.Boundaries(Nbnd))

    f = LinearForm(fes)
    f += Dbndc * (-A*grad(v)*n + alpha*v - B*n*v) * ds(skeleton=True,definedon=mesh.Boundaries(Dbnd)) \
         - Nbndc * v * ds(skeleton=True,definedon=mesh.Boundaries(Nbnd)) \
         + rhs*v*dx
    if uf:
        f += -A*grad(uf)*grad(v) * dx \
            -alpha*jump(uf)*jump(v) * dx(skeleton=True) \
            -(-mean_d(uf)*jump(v)-mean_d(v)*jump(uf)) * dx(skeleton=True) \
            -alpha*uf*v * ds(skeleton=True,definedon=mesh.Boundaries(Dbnd)) \
            -(-A*grad(uf)*n*v-A*grad(v)*n*uf)* ds(skeleton=True,definedon=mesh.Boundaries(Dbnd))
        f += (B*uf*grad(v) - C*uf*v) * dx \
            - (mean_B(uf) * jump(v) + 0.5*sqrt((B*n)**2)*jump(uf)*jump(v)) * dx(skeleton=True) \
            - B*uf*n*v * ds(skeleton=True,definedon=mesh.Boundaries(Nbnd))
        # f.vec.data = f.vec - auf.vec

    with TaskManager():
        a.Assemble()
        f.Assemble()

    return a,f


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


########################################################################
# Wave
########################################################################
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

########################################################################
# Heat
########################################################################

def dgheat(fes, diffusion, ubnd):
    mesh = fes.mesh
    order = fes.globalorder
    u, v = fes.TnT()
    eps = diffusion
    b = CoefficientFunction((0, 1))
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

    f = LinearForm(fes)
    #f += f_coef * v * dx
    f += -b*n* ubnd * v * ds(definedon=mesh.Boundaries("bottom"),skeleton=True)
    f += -eps*grad(v)[0]*n[0] * ubnd * ds(definedon=mesh.Boundaries("left|right"),skeleton=True)
    f += eps*lambd*(order*n[0])**2/h*ubnd*v*ds(definedon=mesh.Boundaries("left|right"),skeleton=True)
    f.Assemble()

    gfu = GridFunction(fes)
    gfu.vec.data = a.mat.Inverse() * f.vec
    return gfu


########################################################################
# LDG
########################################################################
def SolveLapLDG(mesh,order=1,bndc=0,rhs=0,trefftz=False,condense=True):
    """
    >>> mesh2d = Mesh(unit_square.GenerateMesh(maxh=0.3))
    >>> SolveLapLDG(mesh2d,order=5,bndc=exactlap,condense=True) # doctest:+ELLIPSIS
    2...e-09
    """
    # >>> SolveLapLDG(mesh2d,order=5,bndc=exactlap,condense=False) # doctest:+ELLIPSIS
    # 2...e-09
    SetNumThreads(1)
    f1 = L2(mesh,order=order,dgjumps=True)
    # f1 = trefftzfespace(mesh,order=order,eq="laplace")
    f2 = VectorL2(mesh,order=order-1,dgjumps=True)
    fes = f1 * f2
    #print(fes.FreeDofs(coupling=True))
    #print(VDofs)

    n = specialcf.normal(mesh.dim)
    h = specialcf.mesh_size
    u,q = fes.TrialFunction()
    v,r = fes.TestFunction()
 
    C11 = 1/h
    C12 = CF((1,1))
    
    jump_u = (u-u.Other())*n
    jump_v = (v-v.Other())*n
    jump_r = (r-r.Other())*n
    jump_q = (q-q.Other())*n
    mean_u = 0.5 * (u+u.Other())
    mean_v = 0.5 * (v+v.Other())

    a = BilinearForm(fes)
    a += q*r * dx 
    # b
    a += u*div(r)*dx - (mean_u + C12*jump_u)*jump_r*dx(skeleton=True)
    # -b
    a += -1*( v*div(q)*dx - (mean_v + C12*jump_v)*jump_q*dx(skeleton=True) )
    # c
    a += C11*jump_u*jump_v*dx(skeleton=True) + C11*u*v * ds(skeleton=True)

    f = LinearForm(fes)
    f += bndc*r*n * ds(skeleton=True)
    f += C11*bndc*v * ds(skeleton=True)
    f += rhs*v * dx

    with TaskManager():
        a.Assemble()
        f.Assemble()
    # aspy(a.mat)

    # if trefftz:  
        # eps=10**-9
        # Lap = lambda u : sum(Trace(u.Operator('hesse')))
        # u = fes.TrialFunction()
        # fes2 = L2(mesh, order=order-2,  dgjumps=True) * VectorL2(mesh,order=order-1,dgjumps=True) 
        # (v,r) = fes2.TestFunction()
        # op = Lap(u)*v*dx
        # with TaskManager():
            # PP = TrefftzEmbedding(op,fes,test_fes=fes2)
        # PPT = PP.CreateTranspose()
        # with TaskManager():
            # TA = PPT@a.mat@PP
            # TU = TA.Inverse(inverse='sparsecholesky')*(PPT*f.vec)
            # tgfu = GridFunction(fes)
            # tgfu.vec.data = PP*TU
    # else:
    timer = time.time()
    if condense:
        with TaskManager():
            for i in range(f1.ndof,fes.ndof):
               fes.SetCouplingType(i,COUPLING_TYPE.HIDDEN_DOF)
            # VDofs = BitArray(fes.ndof)
            # VDofs[:] = 1
            # for i in range(fes.ndof):
               # VDofs[i] = 0 if fes.CouplingType(i)==COUPLING_TYPE.HIDDEN_DOF else 1
            amat = CondenseDG(a.mat,f.vec,fes) #.FreeDofs(coupling=True))
            gfu = GridFunction(fes)    
            gfu.components[0].vec.data = amat.Inverse() * f.vec[0:f1.ndof] #freedofs=VDofs
    else:
        with TaskManager():
            gfu = GridFunction(fes)    
            gfu.vec.data = a.mat.Inverse() * f.vec
    # print("Time:",time.time()-timer)

    error = sqrt(Integrate((gfu.components[0]-exactlap)**2, mesh))
    return error

########################################################################
# Stokes
########################################################################

def StokesDG(fes, nu, rhs, ubnd=None, bndname="inflow"):
    mesh = fes.components[0].mesh
    k = fes.components[0].globalorder
    u, v = fes.TrialFunction()[0], fes.TestFunction()[0]
    p, q = fes.TrialFunction()[1], fes.TestFunction()[1]
    lam, mu = fes.TrialFunction()[-1], fes.TestFunction()[-1]

    alpha = 20  # interior penalty param
    stab = 1e-9

    n = specialcf.normal(mesh.dim)
    h = specialcf.mesh_size

    jump_u = u - u.Other()
    jump_v = v - v.Other()
    mean_dudn = 0.5 * (grad(u) + grad(u.Other())) * n
    mean_dvdn = 0.5 * (grad(v) + grad(v.Other())) * n
    mean_q = 0.5 * n * (q + q.Other())
    mean_p = 0.5 * n * (p + p.Other())

    a = nu * InnerProduct(grad(u), grad(v)) * dx
    a += nu * alpha * k**2 / h * jump_u * jump_v * dx(skeleton=True)
    a += nu * (-mean_dudn * jump_v - mean_dvdn * jump_u) * dx(skeleton=True)
    a += nu * alpha * k**2 / h * u * v * ds(skeleton=True)
    a += nu * (-grad(u) * n * v - grad(v) * n * u) * ds(skeleton=True)
    a += (mean_p * jump_v + mean_q * jump_u) * dx(skeleton=True)
    a += (p * v * n + q * u * n) * ds(skeleton=True)
    a += (-div(u) * q - div(v) * p) * dx
    a += (p * mu + q * lam) * dx

    c = a
    c += -stab * p * q * dx
    c += stab * lam * mu * dx

    f = rhs * v * dx(bonus_intorder=5)
    if ubnd:
        f += nu * alpha * k**2 / h * ubnd * v * ds(skeleton=True, definedon=mesh.Boundaries(bndname))
        f += nu * (- grad(v) * n * ubnd) * ds(skeleton=True, definedon=mesh.Boundaries(bndname))

    return a,c,f

def SolveStokesDG(mesh, k):
    """ 
    >>> SolveStokesDG(mesh2d, 5) # doctest:+ELLIPSIS
    [7...e-06, 0.0004...]
    """
    nu = 1.0
    zeta = cos(pi*x*(1-x)*y*(1-y))
    pexact = sin(pi*(x+y))
    uexact = CF((zeta.Diff(y), - zeta.Diff(x)))
    graduexact = CF((uexact.Diff(x),uexact.Diff(y)),dims=(2,2)).trans
    f1 = - nu*uexact[0].Diff(x).Diff(x) - nu*uexact[0].Diff(y).Diff(y) + pexact.Diff(x)
    f2 = - nu*uexact[1].Diff(x).Diff(x) - nu*uexact[1].Diff(y).Diff(y) + pexact.Diff(y)
    rhs = CF((f1,f2))

    V = VectorL2(mesh, order=k, dgjumps=True)
    Q = L2(mesh, order=k - 1, dgjumps=True)
    Z = NumberSpace(mesh)
    fes = V * Q * Z

    ah,ch,fh = StokesDG(fes, nu, rhs)

    a = BilinearForm(ah)
    a.Assemble()
    c = BilinearForm(ch)
    c.Assemble()
    f = LinearForm(fh)
    f.Assemble()

    gfu = GridFunction(fes)
    gfu.vec.data = CGSolver(a.mat, c.mat.Inverse()) * f.vec
    ndof = fes.ndof

    uh, ph = gfu.components[0:2]
    return [sqrt(Integrate(InnerProduct(uexact-uh,uexact-uh),mesh)) , sqrt(Integrate(InnerProduct(pexact-ph,pexact-ph),mesh)) ]


if __name__ == "__main__":
    import doctest
    doctest.testmod()
