from ngsolve import *
from trefftzngs import *
import numpy as np
import scipy.linalg as spla


def DGeqsys(fes,U0,v0,sig0,c,gD,fullsys=False):
	D = sig0.dim
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
	#gfu.vec.data = a.mat.Inverse() * f.vec

	timelike = n_x*n_x#n_x**2 #IfPos(n_t,0,IfPos(-n_t,0,1)) # n_t=0
	spacelike = n_t**2 #IfPos(n_x,0,IfPos(-n_x,0,1)) # n_x=0

	alpha = 0 #pow(10,5)
	beta = 0 #pow(10,5)
	gamma = 1

	a = BilinearForm(fes)
	if(fullsys==True):
	    HV = V.Operator("hesse")
	    a += SymbolicBFI(  -v*(-HV[0]+pow(c,-2)*HV[3]) - sig*(-HV[1]+HV[2])  )
	# a += SymbolicBFI( spacelike * ( IfPos(n_t,v,vo)*(pow(c,-2)*jump_wt+jump_taux) + IfPos(n_t,sig,sigo)*(jump_wx+jump_taut) ) ,VOL,  skeleton=True ) #space like faces
	a += SymbolicBFI( spacelike * ( IfPos(n_t,v,vo)*(pow(c,-2)*jump_wt) + IfPos(n_t,sig,sigo)*(jump_taut) ) ,VOL,  skeleton=True ) #space like faces, no jump in x since horizontal
	a += SymbolicBFI( timelike 	* ( mean_v*jump_taux + mean_sig*jump_wx + alpha*jump_vx*jump_wx + beta*jump_sigx*jump_taux ) ,VOL, skeleton=True ) #time like faces
	a += SymbolicBFI( spacelike * IfPos(n_t,1,0) * ( pow(c,-2)*v*w + sig*tau ), BND, skeleton=True) #t=T (or *x)
	a += SymbolicBFI( timelike 	* ( sig*n_x*w + alpha*v*w ), BND, skeleton=True) #dirichlet boundary 'timelike'
	#a += SymbolicBFI( spacelike * ( gamma * (-n_t*jump_Ut)*IfPos(n_t,V.Other(),V) ) ,VOL,  skeleton=True ) #correction term to recover sol of second order system
	a += SymbolicBFI( spacelike * ( gamma * (-jump_Ut)*IfPos(n_t,V.Other(),V) ) ,VOL,  skeleton=True ) #correction term to recover sol of second order system
	a += SymbolicBFI( spacelike * ( gamma * IfPos(-n_t,1,0) * U*V ) ,BND,  skeleton=True ) #BND correction term to recover sol of second order system

	a.Assemble()

	f = LinearForm(fes)
	f += SymbolicLFI( spacelike * IfPos(-n_t,1,0) *  ( pow(c,-2)*v0*w + sig0*tau ), BND, skeleton=True) #t=0 (or *(1-x))
	f += SymbolicLFI( timelike 	* ( gD * (alpha*w - tau*n_x) ), BND, skeleton=True) #dirichlet boundary 'timelike'
	f += SymbolicLFI( spacelike * gamma * IfPos(-n_t,1,0) *  ( (U0)*V ) ,BND,  skeleton=True ) #rhs correction term to recover sol of second order system
	f.Assemble()

	return [a,f]




def DGsolve(fes,a,f):
	gfu = GridFunction(fes, name="uDG")

	nmat = np.zeros((a.mat.height,a.mat.width))
	nvec = np.zeros(a.mat.width)

	for i in range(a.mat.width):#gfu.vec.data = a.mat.Inverse() * f.vec
		nvec[i] = f.vec[i]/sqrt(a.mat[i,i])
		for j in range(a.mat.height):
			nmat[j,i] = a.mat[j,i]/sqrt(a.mat[i,i]*a.mat[j,j])
	#nvec = f.vec.FV().NumPy()
	sol = spla.solve(nmat,nvec)

	for i in range(a.mat.height):
		gfu.vec[i] = sol[i]/sqrt(a.mat[i,i])

    #print(min(np.linalg.eigvalsh(0.5*(nmat + nmat.transpose()))))
	#nmatinv = np.linalg.inv(nmat)
    #print(min(np.linalg.eigvalsh(0.5*(nmatinv + nmatinv.transpose()))))
    #nvec = f.vec.FV().NumPy() #nvec
    # print("cond nmat: ", np.linalg.cond(nmat))
	cond = np.linalg.cond(nmat)
	return [gfu,cond]
