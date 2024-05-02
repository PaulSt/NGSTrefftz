# import sys, os
# sys.path.append(os.path.join(os.path.dirname(sys.path[0])))
from ngsolve import *
from ngstrefftz import *
from dg import *
import time
from netgen.geom2d import unit_square
from netgen.csg import unit_cube
SetNumThreads(3)

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

def PySVDConstraintTrefftz(op,fes,cop_lhs,cop_rhs,fes_constraint,eps):
    mesh = fes.mesh
    order = fes.globalorder
    opbf = BilinearForm(fes)
    opbf += op
    opbf.Assemble()
    rows,cols,vals = opbf.mat.COO()
    A = sp.sparse.csr_matrix((vals,(rows,cols)))
    A = A.todense()

    copbf_lhs = BilinearForm(trialspace=fes, testspace=fes_constraint)
    copbf_lhs += cop_lhs # fes.TrialFunction()*fes_constraint.TestFunction()*dx(element_boundary=True) # cop_lhs
    copbf_lhs.Assemble()
    copbf_rhs = BilinearForm(fes_constraint)
    copbf_rhs += cop_rhs
    copbf_rhs.Assemble()


    rows,cols,vals = copbf_lhs.mat.COO()
    B1 = sp.sparse.csr_matrix((vals,(rows,cols)))
    B1 = B1.todense()
    print(B1.shape)

    rows,cols,vals = copbf_rhs.mat.COO()
    B2 = sp.sparse.csr_matrix((vals,(rows,cols)))
    B2 = B2.todense()

    #localndof = int(fes.ndof/mesh.ne)
    trefftzndof = 2*order+1-3

    P = np.zeros([A.shape[1],(trefftzndof)*mesh.ne+fes_constraint.ndof])


    for el, el_c in zip(fes.Elements(),fes_constraint.Elements()):
        nr = el.nr
        dofs = el.dofs
        dofs_c = el_c.dofs
        print("dofs:", dofs)
        print("dofs_c:", dofs_c)
        elmat_l = A[dofs,:][:,dofs]
        elmat_b1 = B1[dofs_c,:][:,dofs]
        elmat_b2 = B2[dofs_c,:][:,dofs_c]

        print("elmat_b1",elmat_b1.shape)
        print("elmat_l",elmat_l.shape)
        elmat_lhs = np.vstack([elmat_b1,elmat_l])

        elmat_rhs = np.vstack([elmat_b2,np.zeros([len(dofs),len(dofs_c)])])

        print("elmat_lhs",elmat_lhs)
        print("elmat_rhs",elmat_rhs)

        U,s,V = np.linalg.svd(elmat_lhs)
        T1 = V.T * np.hstack([np.diag(s),np.zeros((V.shape[0],U.shape[0]-V.shape[0]))]) * U.T * elmat_rhs
        print("T1",T1)

        P[dofs,:][:,dofs_c] = T1[:,0:len(dofs_c)] 

        input("")
    return True





if __name__ == "__main__":

    #eps = 10**-8
    mesh2d = Mesh(unit_square.GenerateMesh(maxh=0.4))

    fes = L2(mesh2d, order=1,  dgjumps=True)#,all_dofs_together=True)
    u,v = fes.TnT()
    uh = u.Operator("hesse")
    vh = v.Operator("hesse")
    op = (uh[0,0]+uh[1,1])*(vh[0,0]+vh[1,1])*dx

    fes_constraint = FacetFESpace(mesh2d, order=0)#,all_dofs_together=True)
    uF, vF = fes_constraint.TnT()
    cop_lhs = u*vF*dx(element_boundary=True)
    cop_rhs = uF*vF*dx(element_boundary=True)

    P = PySVDConstraintTrefftz(op,fes,cop_lhs,cop_rhs,fes_constraint,eps)

    a,f = dgell(fes,exactlap)
    rows,cols,vals = a.mat.COO()
    A = sp.sparse.csr_matrix((vals,(rows,cols)))
    A = A.todense()
    
    TA = P.transpose()@A@P  
    tsol = np.linalg.solve(TA,P.transpose()@f.vec.FV())
    tpgfu = GridFunction(fes)
    tpgfu.vec.data = P@tsol
    print("error : ", sqrt(Integrate((tpgfu-exactlap)**2, mesh2d)) )