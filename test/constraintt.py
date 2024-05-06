# import sys, os
# sys.path.append(os.path.join(os.path.dirname(sys.path[0])))
from ngsolve import *
from ngstrefftz import *

# from dg import *
import dg
import time
from netgen.geom2d import unit_square
from netgen.csg import unit_cube
import scipy as sp
import numpy as np
import matplotlib.pyplot as plt


SetNumThreads(3)

########################################################################
# PySVDTrefftz
########################################################################
# import matplotlib.pylab as plt
# def trunc(values, decs=0):
# return np.round(values*10**decs)/(10**decs)
# def spspy(sparsemat):
# rows,cols,vals = sparsemat.COO()
# A = sp.sparse.csr_matrix((vals,(rows,cols)))
# # A = A.todense()
# plt.spy(A); plt.plot(); plt.show()


def LocalP(A, eps=10**-9, printP=0):
    U, s, V = np.linalg.svd(A)
    nz = (s < eps).sum()
    nnz = (s >= eps).sum()
    S = np.diag(s)
    P = U[:, nnz:]
    # if printP==True:
    # print(trunc(A,decs=0))
    # print(trunc(U,decs=3))
    # print(trunc(P,decs=3))
    return P


def PySVDConstraintTrefftz(
    op: comp.SumOfIntegrals,
    fes: FESpace,
    cop_lhs: comp.SumOfIntegrals,
    cop_rhs: comp.SumOfIntegrals,
    fes_constraint: FESpace,
    eps: float,
    debug: bool = False,
) -> np.ndarray:
    """
    produces an embedding matrix P

    `op`: the differential operation

    `fes`: the finite element space of `op`

    `cop_lhs`: left hand side of the constraint operation

    `cop_rhs`: right hand side of the constraint operation

    `fes_constraint`: finite element space of the constraint operation

    `eps`: eigenvalues less than eps will be treated as 0

    `debug`: if `True`, print debug messages. Default: `False`

    returns: P

    raises: LinAlgError: if the imput is not sound, a non-invertible matrix may
        be tried to be inverted
    """
    mesh = fes.mesh
    order = fes.globalorder

    # let L be the matrix corrensponding to
    # the differential operator op
    opbf = BilinearForm(fes)
    opbf += op
    opbf.Assemble()
    rows, cols, vals = opbf.mat.COO()
    L = sp.sparse.csr_matrix((vals, (rows, cols)))
    L = L.todense()

    # let B1 be the matrix corrensponding to
    # the left hand side constraint operator cop_lhs
    copbf_lhs = BilinearForm(trialspace=fes, testspace=fes_constraint)
    copbf_lhs += cop_lhs
    copbf_lhs.Assemble()
    rows, cols, vals = copbf_lhs.mat.COO()
    B1 = sp.sparse.csr_matrix((vals, (rows, cols)))
    B1 = B1.todense()
    if debug:
        print("B1.shape", B1.shape)

    # let B2 be the matrix corrensponding to
    # the right hand side constraint operator cop_rhs
    copbf_rhs = BilinearForm(fes_constraint)
    copbf_rhs += cop_rhs
    copbf_rhs.Assemble()
    rows, cols, vals = copbf_rhs.mat.COO()
    B2 = sp.sparse.csr_matrix((vals, (rows, cols)))
    B2 = B2.todense()

    # localndof = int(fes.ndof/mesh.ne)
    trefftzndof = 2 * order + 1 - 3

    P = np.zeros([L.shape[1], (trefftzndof) * mesh.ne + fes_constraint.ndof])

    # solve the following linear system in an element-wise fashion:
    # L @ T1 = B for the unknown matrix T1,
    # with the given matrices:
    #     /   \    /   \
    #  A= |B_1| B= |B_2|
    #     | L |    | 0 |
    #     \   /    \   /
    for el, el_c in zip(fes.Elements(), fes_constraint.Elements()):
        nr = el.nr
        dofs = el.dofs
        dofs_c = el_c.dofs

        if debug:
            print("dofs:", dofs)
            print("dofs_c:", dofs_c)

        # produce local sub-matrices from the global matrices L, B1, B2
        elmat_l = L[dofs, :][:, dofs]
        elmat_b1 = B1[dofs_c, :][:, dofs]
        elmat_b2 = B2[dofs_c, :][:, dofs_c]

        if debug:
            print("elmat_b1", elmat_b1.shape)
            print("elmat_l", elmat_l.shape)

        #     /   \    /   \
        #  A= |B_1| B= |B_2|
        #     | L |    | 0 |
        #     \   /    \   /
        elmat_a = np.vstack([elmat_b1, elmat_l])
        elmat_b = np.vstack([elmat_b2, np.zeros([len(dofs), len(dofs_c)])])

        if debug:
            print("elmat_a", elmat_a)
            print("elmat_b", elmat_b)

        # A = U @ s @ V, singular value decomposition
        U, s, V = np.linalg.svd(elmat_a)

        # pseudo inverse of s
        s_inv = np.hstack(
            [np.diag(1.0 / s), np.zeros((V.shape[0], U.shape[0] - V.shape[0]))]
        )

        if debug:
            print("U", U)
            print("s_inv", s_inv)
            print("V", V)

        # solve A @ T1 = B
        # i.e. T1 = A^{-1} @ B
        # for the unknown T1
        T1 = V.T @ s_inv @ U.T @ elmat_b

        if debug:
            print("T1", T1)

        # place the local solution T1
        # into the global solution P
        P[dofs, :][:, dofs_c] = T1[:, 0 : len(dofs_c)]

        # gfu = GridFunction(fes)
        # print("T1 slice:\n", T1[:, 0].flatten())
        # gfu.vec.data[dofs, :][:, dofs_c] = T1[:, 0 : len(dofs_c)]
        # Draw(gfu)
        if debug:
            plt.plot(s, marker="x")
            plt.show()
            # input("")
    return P


if __name__ == "__main__":

    # eps = 10**-8
    mesh2d = Mesh(unit_square.GenerateMesh(maxh=0.4))

    fes = L2(mesh2d, order=1, dgjumps=True)  # ,all_dofs_together=True)
    u, v = fes.TnT()
    uh = u.Operator("hesse")
    vh = v.Operator("hesse")
    op = (uh[0, 0] + uh[1, 1]) * (vh[0, 0] + vh[1, 1]) * dx

    fes_constraint = FacetFESpace(mesh2d, order=0)  # ,all_dofs_together=True)
    uF, vF = fes_constraint.TnT()
    cop_lhs = u * vF * dx(element_boundary=True)
    cop_rhs = uF * vF * dx(element_boundary=True)

    P = PySVDConstraintTrefftz(
        op, fes, cop_lhs, cop_rhs, fes_constraint, dg.eps, debug=True
    )

    a, f = dg.dgell(fes, dg.exactlap)
    rows, cols, vals = a.mat.COO()
    A = sp.sparse.csr_matrix((vals, (rows, cols)))
    A = A.todense()

    TA = P.transpose() @ A @ P
    tsol = np.linalg.solve(TA, P.transpose() @ f.vec.FV())
    tpgfu = GridFunction(fes)
    tpgfu.vec.data = P @ tsol
    print("error : ", sqrt(Integrate((tpgfu - dg.exactlap) ** 2, mesh2d)))
