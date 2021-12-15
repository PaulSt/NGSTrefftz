from ngsolve import *
from tngs import *
import time
from netgen.geom2d import unit_square
import scipy as sp
import numpy as np

# import matplotlib.pylab as plt
# def trunc(values, decs=0):
    # return np.round(values*10**decs)/(10**decs)
# def spspy(sparsemat):
    # rows,cols,vals = sparsemat.COO()
    # A = sp.sparse.csr_matrix((vals,(rows,cols)))
    # # A = A.todense()
    # plt.spy(A); plt.plot(); plt.show()

order = 5
exact = exp(x)*sin(y)
eps = 10**-8
mesh = Mesh(unit_square.GenerateMesh(maxh=0.3))
SetNumThreads(3)

########################################################################
# L2
########################################################################
def dglap(fes,bndc):
    """
    >>> fes = L2(mesh, order=order,  dgjumps=True)#,all_dofs_together=True)
    >>> a,f = dglap(fes,exact)
    >>> gfu = GridFunction(fes)
    >>> gfu.vec.data = a.mat.Inverse() * f.vec
    >>> sqrt(Integrate((gfu-exact)**2, mesh)) # doctest:+ELLIPSIS
    4...e-09
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
         +(-n*grad(v)*bndc)* ds(skeleton=True)
    f.Assemble()
    return a,f


########################################################################
# PySVDTrefftz
########################################################################

def LocalP(A,eps=10**-9,printP=0):
    U,s,V = np.linalg.svd(A)
    nz = (s<eps).sum()
    nnz = (s>=eps).sum()
    S=np.diag(s)
    P = U[:,nnz:]
    if printP==True:
        print(trunc(A,decs=0))
        print(trunc(U,decs=3))
        print(trunc(P,decs=3))
    return P

def PySVDTrefftz(op,fes,eps):
    """
    >>> fes = L2(mesh, order=order,  dgjumps=True)#,all_dofs_together=True)
    >>> u,v = fes.TnT()
    >>> uh = u.Operator("hesse")
    >>> vh = v.Operator("hesse")
    >>> op = (uh[0,0]+uh[1,1])*(vh[0,0]+vh[1,1])*dx
    >>> P = PySVDTrefftz(op,fes,eps)
    >>> a,f = dglap(fes,exact)
    >>> rows,cols,vals = a.mat.COO()
    >>> A = sp.sparse.csr_matrix((vals,(rows,cols)))
    >>> A = A.todense()
    >>> TA = P.transpose()@A@P

    >>> tsol = np.linalg.solve(TA,P.transpose()@f.vec.FV())
    >>> tpgfu = GridFunction(fes)
    >>> tpgfu.vec.data = P@tsol
    >>> sqrt(Integrate((tpgfu-exact)**2, mesh)) # doctest:+ELLIPSIS
    1...e-08
    """

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
# SVDTrefftz
########################################################################
def testsvdtrefftz():
    """
    >>> testsvdtrefftz() # doctest:+ELLIPSIS
    1...e-08
    """
    start = time.time()
    fes = L2(mesh, order=order,  dgjumps=True)#,all_dofs_together=True)
    u,v = fes.TnT()
    uh = u.Operator("hesse")
    vh = v.Operator("hesse")
    op = (uh[0,0]+uh[1,1])*(vh[0,0]+vh[1,1])*dx
    # op = grad(u)*grad(v) * dx
    # op = InnerProduct(u.Operator("hesse"),v.Operator("hesse"))*dx
    with TaskManager():
        PP = SVDTrefftz(op,fes,eps)
    # spspy(PP)
    PPT = PP.CreateTranspose()
    a,f = dglap(fes,exact)
    TA = PPT@a.mat@PP
    TU = TA.Inverse()*(PPT*f.vec)
    tpgfu = GridFunction(fes)
    tpgfu.vec.data = PP*TU
    # print("trefftz time: ", time.time()-start)
    # for t in Timers():
        # if 'svdt' in t['name']:
            # print(t)
    return sqrt(Integrate((tpgfu-exact)**2, mesh))




########################################################################
# Trefftz
########################################################################
def testlaptrefftz():
    """
    >>> testlaptrefftz() # doctest:+ELLIPSIS
    1...e-08
    """
    fes = FESpace("trefftzfespace",mesh,order=order,eq="laplace")
    a,f = dglap(fes,exact)
    gfu = GridFunction(fes)
    gfu.vec.data = a.mat.Inverse() * f.vec
    return sqrt(Integrate((gfu-exact)**2, mesh))

if __name__ == "__main__":
    import doctest
    doctest.testmod()
