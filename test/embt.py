# import sys, os
# sys.path.append(os.path.join(os.path.dirname(sys.path[0])))
from ngsolve import *
from ngstrefftz import *
import time
from netgen.geom2d import unit_square
from netgen.csg import unit_cube

order = 5
exactlap = exp(x)*sin(y)
exactpoi = sin(x)*sin(y)
eps = 10**-8
mesh2d = Mesh(unit_square.GenerateMesh(maxh=0.3))
mesh3d = Mesh(unit_cube.GenerateMesh(maxh = 1))
SetNumThreads(3)

Lap = lambda u : sum(Trace(u.Operator('hesse')))

########################################################################
# L2
########################################################################
def dglap(fes,bndc,rhs=0,uf=0):
    """
    >>> fes = L2(mesh2d, order=order,  dgjumps=True)#,all_dofs_together=True)
    >>> a,f = dglap(fes,exactlap)
    >>> gfu = GridFunction(fes)
    >>> gfu.vec.data = a.mat.Inverse() * f.vec
    >>> sqrt(Integrate((gfu-exactlap)**2, mesh2d)) # doctest:+ELLIPSIS
    3...e-09
    """
    mesh = fes.mesh
    order = fes.globalorder
    alpha = 4
    n = specialcf.normal(mesh.dim)
    h = specialcf.mesh_size
    u = fes.TrialFunction()
    v = fes.TestFunction()

    jump_u = u-u.Other()
    jump_v = v-v.Other()
    mean_dudn = 0.5*n * (grad(u)+grad(u.Other()))
    mean_dvdn = 0.5*n * (grad(v)+grad(v.Other()))

    a = BilinearForm(fes)
    a += grad(u)*grad(v) * dx \
        +alpha*order**2/h*jump_u*jump_v * dx(skeleton=True) \
        +(-mean_dudn*jump_v-mean_dvdn*jump_u) * dx(skeleton=True) \
        +alpha*order**2/h*u*v * ds(skeleton=True) \
        +(-n*grad(u)*v-n*grad(v)*u)* ds(skeleton=True)
    a.Assemble()

    f = LinearForm(fes)
    f += alpha*order**2/h*bndc*v * ds(skeleton=True) \
         +(-n*grad(v)*bndc)* ds(skeleton=True) \
         +rhs*v*dx
    f.Assemble()

    if uf:
        fes2 = L2(fes.mesh, order=fes.globalorder, dgjumps=True)
        u = fes2.TrialFunction()

        jump_u = (u-u.Other())*n
        jump_v = (v-v.Other())*n
        mean_dudn = 0.5 * (grad(u)+grad(u.Other()))
        mean_dvdn = 0.5 * (grad(v)+grad(v.Other()))

        a2 = BilinearForm(fes2,fes,nonassemble=True)
        a2 += grad(u)*grad(v) * dx \
            +alpha*order**2/h*jump_u*jump_v * dx(skeleton=True) \
            +(-mean_dudn*jump_v-mean_dvdn*jump_u) * dx(skeleton=True) \
            +alpha*order**2/h*u*v * ds(skeleton=True) \
            +(-n*grad(u)*v-n*grad(v)*u)* ds(skeleton=True)
        auf = f.vec.CreateVector()
        a2.Apply(uf.vec,auf)
        f.vec.data = f.vec - auf
    return a,f

########################################################################
# PySVDTrefftz
########################################################################
try:
    import scipy as sp
    import numpy as np
except ImportError as e:
    pass

# import matplotlib.pylab as plt
# def trunc(values, decs=0):
    # return np.round(values*10**decs)/(10**decs)
# def spspy(sparsemat):
    # rows,cols,vals = sparsemat.COO()
    # A = sp.sparse.csr_matrix((vals,(rows,cols)))
    # # A = A.todense()
    # plt.spy(A); plt.plot(); plt.show()


def LocalP(A,eps=10**-9,printP=0):
    U,s,V = np.linalg.svd(A)
    nz = (s<eps).sum()
    nnz = (s>=eps).sum()
    S=np.diag(s)
    P = U[:,nnz:]
    # if printP==True:
        # print(trunc(A,decs=0))
        # print(trunc(U,decs=3))
        # print(trunc(P,decs=3))
    return P

def PySVDTrefftz(op,fes,eps):
    """
    # >>> fes = L2(mesh2d, order=order,  dgjumps=True)#,all_dofs_together=True)
    # >>> u,v = fes.TnT()
    # >>> uh = u.Operator("hesse")
    # >>> vh = v.Operator("hesse")
    # >>> op = (uh[0,0]+uh[1,1])*(vh[0,0]+vh[1,1])*dx
    # >>> P = PySVDTrefftz(op,fes,eps)
    # >>> a,f = dglap(fes,exactlap)
    # >>> rows,cols,vals = a.mat.COO()
    # >>> A = sp.sparse.csr_matrix((vals,(rows,cols)))
    # >>> A = A.todense()
    # >>> TA = P.transpose()@A@P

    # >>> tsol = np.linalg.solve(TA,P.transpose()@f.vec.FV())
    # >>> tpgfu = GridFunction(fes)
    # >>> tpgfu.vec.data = P@tsol
    # >>> sqrt(Integrate((tpgfu-exactlap)**2, mesh2d)) # doctest:+ELLIPSIS
    # 1...e-08
    """
    mesh = fes.mesh
    opbf = BilinearForm(fes)
    opbf += op
    opbf.Assemble()
    rows,cols,vals = opbf.mat.COO()
    A = sp.sparse.csr_matrix((vals,(rows,cols)))
    A = A.todense()
    localndof = int(fes.ndof/mesh.ne)
    trefftzndof = 2*order+1
    P = np.zeros([A.shape[1],(trefftzndof)*mesh.ne])
    for i in range(mesh.ne):
        # U,s,V = sp.linalg.svd(A[i*localndof:(i+1)*localndof,i*localndof:(i+1)*localndof])
        elmat = A[i*localndof:(i+1)*localndof,i*localndof:(i+1)*localndof]
        LP = LocalP(elmat,eps) #,i+1==mesh.ne)
        P[i*localndof:(i+1)*localndof,i*trefftzndof:(i+1)*trefftzndof] = LP
    # plt.spy(P); plt.plot(); plt.show()
    return P





########################################################################
# EmbTrefftz
########################################################################
def testembtrefftz(fes):
    """
    >>> fes = L2(mesh2d, order=order,  dgjumps=True)#,all_dofs_together=True)
    >>> testembtrefftz(fes) # doctest:+ELLIPSIS
    8...e-09
    """
    start = time.time()
    mesh = fes.mesh
    u,v = fes.TnT()
    uh = u.Operator("hesse")
    vh = v.Operator("hesse")
    op = (uh[0,0]+uh[1,1])*(vh[0,0]+vh[1,1])*dx
    # op = grad(u)*grad(v) * dx
    # op = InnerProduct(u.Operator("hesse"),v.Operator("hesse"))*dx
    with TaskManager():
        PP = TrefftzEmbedding(op,fes,eps)
    # spspy(PP)
    PPT = PP.CreateTranspose()
    a,f = dglap(fes,exactlap)
    TA = PPT@a.mat@PP
    TU = TA.Inverse()*(PPT*f.vec)
    tpgfu = GridFunction(fes)
    tpgfu.vec.data = PP*TU
    # print("trefftz time: ", time.time()-start)
    # for t in Timers():
        # if 'svdt' in t['name']:
            # print(t)
    return sqrt(Integrate((tpgfu-exactlap)**2, mesh))


def testembtrefftz_mixed(fes):
    """
    >>> fes = L2(mesh2d, order=order,  dgjumps=True)#,all_dofs_together=True)
    >>> testembtrefftz_mixed(fes) # doctest:+ELLIPSIS
    8...e-09
    """
    mesh = fes.mesh
    test_fes = L2(mesh, order=fes.globalorder-2,  dgjumps=True)#,all_dofs_together=True)
    u=fes.TrialFunction()
    v=test_fes.TestFunction()
    op = Lap(u)*(v)*dx
    startsvd = time.time()
    with TaskManager():
        PP = TrefftzEmbedding(op,fes,test_fes=test_fes)
    PPT = PP.CreateTranspose()
    a,f = dglap(fes,exactlap)
    TA = PPT@a.mat@PP
    TU = TA.Inverse()*(PPT*f.vec)
    tpgfu = GridFunction(fes)
    tpgfu.vec.data = PP*TU
    return sqrt(Integrate((tpgfu-exactlap)**2, mesh))


def testembtrefftznonsym(fes):
    """
    >>> fes = L2(mesh2d, order=order,  dgjumps=True)#,all_dofs_together=True)
    >>> testembtrefftznonsym(fes) # doctest:+ELLIPSIS
    8...e-09
    """
    start = time.time()
    mesh = fes.mesh
    u,v = fes.TnT()
    uh = u.Operator("hesse")
    vh = v.Operator("hesse")
    op = (uh[0,0]+uh[1,1])*v*dx
    with TaskManager():
        PP = TrefftzEmbedding(op,fes,eps)
    PPT = PP.CreateTranspose()
    a,f = dglap(fes,exactlap)
    TA = PPT@a.mat@PP
    TU = TA.Inverse()*(PPT*f.vec)
    tpgfu = GridFunction(fes)
    tpgfu.vec.data = PP*TU
    return sqrt(Integrate((tpgfu-exactlap)**2, mesh))

def testembtrefftzpoi(fes):
    """
    >>> fes = L2(mesh2d, order=order,  dgjumps=True)#,all_dofs_together=True)
    >>> testembtrefftzpoi(fes) # doctest:+ELLIPSIS
    3...e-09
    """
    mesh = fes.mesh
    start = time.time()
    u,v = fes.TnT()
    uh = u.Operator("hesse")
    vh = v.Operator("hesse")
    op = (uh[0,0]+uh[1,1])*(vh[0,0]+vh[1,1])*dx
    rhs = -exactpoi.Diff(x).Diff(x)-exactpoi.Diff(y).Diff(y)
    lop = -rhs*(vh[0,0]+vh[1,1])*dx
    with TaskManager():
        PP,ufv = TrefftzEmbedding(op,fes,lop,eps)
    PPT = PP.CreateTranspose()
    a,f = dglap(fes,exactpoi,rhs)
    TA = PPT@a.mat@PP
    TU = TA.Inverse()*(PPT*(f.vec-a.mat*ufv))
    tpgfu = GridFunction(fes)
    tpgfu.vec.data = PP*TU+ufv
    return sqrt(Integrate((tpgfu-exactpoi)**2, mesh))


def testembtrefftzpoi_mixed(fes):
    """
    >>> fes = L2(mesh2d, order=order,  dgjumps=True)#,all_dofs_together=True)
    >>> testembtrefftzpoi_mixed(fes) # doctest:+ELLIPSIS
    3...e-09
    """
    mesh = fes.mesh
    test_fes = L2(mesh, order=fes.globalorder-2,  dgjumps=True)#,all_dofs_together=True)
    u=fes.TrialFunction()
    v=test_fes.TestFunction()
    op = Lap(u)*(v)*dx
    rhs = -exactpoi.Diff(x).Diff(x)-exactpoi.Diff(y).Diff(y)
    lop = -rhs*v*dx
    with TaskManager():
        PP,ufv = TrefftzEmbedding(op,fes,lop,test_fes=test_fes)
    PPT = PP.CreateTranspose()
    a,f = dglap(fes,exactpoi,rhs)
    TA = PPT@a.mat@PP
    TU = TA.Inverse()*(PPT*(f.vec-a.mat*ufv))
    tpgfu = GridFunction(fes)
    tpgfu.vec.data = PP*TU+ufv
    return sqrt(Integrate((tpgfu-exactpoi)**2, mesh))


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

def testembtrefftzhelm(fes):
    """
    >>> fes = L2(mesh2d, order=order,  dgjumps=True, complex=True)
    >>> testembtrefftzhelm(fes) # doctest:+ELLIPSIS
    8...e-10
    """
    mesh = fes.mesh
    omega=1
    exact = exp(1j*sqrt(0.5)*(x+y))
    n = specialcf.normal(mesh.dim)
    bndc = CoefficientFunction((sqrt(0.5)*1j*exact, sqrt(0.5)*1j*exact))*n + 1j*omega*exact
    # test_fes = L2(mesh2d, order=order-2,  dgjumps=True, complex=True)
    u=fes.TrialFunction()
    v=fes.TestFunction()
    # v=test_fes.TestFunction()
    op = (-Lap(u))*(Lap(v))*dx#+1j*(-Lap(u)-u)*(Lap(v))*dx
    op += (-u)*(Lap(v))*dx#+1j*(-Lap(u)-u)*(Lap(v))*dx
    # op = (-Lap(u)-u)*((v))*dx+1j*(-Lap(u)-u)*((v))*dx
    with TaskManager():
        PP = TrefftzEmbedding(op,fes,eps=10**-8)
        # PP = TrefftzEmbedding(op,fes,test_fes=test_fes)
    PPT = PP.CreateTranspose()
    a,f = dghelm(fes,None,bndc,omega)
    TA = PPT@a.mat@PP
    TU = TA.Inverse()*(PPT*f.vec)
    tpgfu = GridFunction(fes)
    tpgfu.vec.data = PP*TU
    return sqrt(Integrate(InnerProduct(tpgfu-exact,tpgfu-exact).real, mesh))


########################################################################
# EmbTrefftzFESpace
########################################################################


def testembtrefftzfes(mesh,order):
    """
    >>> testembtrefftzfes(mesh2d,order) # doctest:+ELLIPSIS
    3...e-09
    """
    # exact = sin(10*x)*sin(10*y)
    # rhs = 200*sin(10*x)*sin(10*y)
    rhs = -exactpoi.Diff(x).Diff(x)-exactpoi.Diff(y).Diff(y)
    eps = 10**-7

    fes = L2(mesh, order=order, dgjumps=True)
    etfes = EmbeddedTrefftzFES(fes)

    u,v = fes.TnT()
    op = Lap(u)*Lap(v)*dx
    lop = -rhs*Lap(v)*dx

    uf = GridFunction(fes)
    uf.vec.data = etfes.SetOp(op,lf=lop,eps=eps)

    # etfes = Compress(etfes, etfes.FreeDofs())
    a,f = dglap(etfes,exactpoi,rhs,uf)
    # import pdb; pdb.set_trace()

    inv = a.mat.Inverse(inverse="sparsecholesky")
    gfu = GridFunction(etfes)
    gfu.vec.data = inv * f.vec
    nsol = gfu + uf

    return sqrt(Integrate((nsol-exactpoi)**2, mesh))



if __name__ == "__main__":
    import doctest
    doctest.testmod()
