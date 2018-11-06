# -*- mode: python-mode; python-indent-offset: 4 -*-
import netgen.meshing as ngm
from ngsolve import *
from netgen.geom2d import unit_square
from ngsolve import TensorProductTools
import netgen.geom2d as ngeom2d

def QadSegMesh(n,x0,x1):
    ngmesh = ngm.Mesh(dim=1)
    pids = []
    for i in range(n+1):
        pids.append (ngmesh.Add (ngm.MeshPoint(ngm.Pnt(x0+(x1-x0)*i/n*i/n, 0, 0))))
    TensorProductTools.AddEdgeEls(0,1,1,n,pids,ngmesh)
    ngmesh.Add (ngm.Element0D( pids[0], index=1))
    ngmesh.Add (ngm.Element0D( pids[n], index=2))
    ngmesh.SetBCName(0,"left")
    ngmesh.SetBCName(1,"right")
    return ngmesh


def PeriodicProdMesh(ngmeshbase,t_steps):
    ngmesh = ngm.Mesh()
    ngmesh.dim=3
    pnums = []

    for i,p in enumerate(ngmeshbase.Points()):
        x,y,z = p.p
        pnums.append( ngmesh.Add(ngm.MeshPoint(ngm.Pnt(x,y,0))) )
        pnums.append( ngmesh.Add(ngm.MeshPoint(ngm.Pnt(x,y,t_steps))) )
        # pnums.append( ngmesh.Add(ngm.MeshPoint(ngm.Pnt(x,y,2*t_steps))) )
        ngmesh.AddPointIdentification(pnums[-2],pnums[-1], identnr=1,type=2) #master, slave
    El1d = ngmeshbase.Elements1D()
    El2d = ngmeshbase.Elements2D()

    ngmesh.SetMaterial(1, "mat")
    for el in El2d:
        ngmesh.Add(ngm.Element3D(1, [pnums[2*(el.points[0].nr-1)]
                                    ,pnums[2*(el.points[1].nr-1)]
                                    ,pnums[2*(el.points[2].nr-1)]
                                    ,pnums[2*(el.points[0].nr-1)+1]
                                    ,pnums[2*(el.points[1].nr-1)+1]
                                    ,pnums[2*(el.points[2].nr-1)+1] ]))
        # ngmesh.Add(ngm.Element3D(1, [pnums[3*(el.points[0].nr-1)+1]
        #                             ,pnums[3*(el.points[1].nr-1)+1]
        #                             ,pnums[3*(el.points[2].nr-1)+1]
        #                             ,pnums[3*(el.points[0].nr-1)+2]
        #                             ,pnums[3*(el.points[1].nr-1)+2]
        #                             ,pnums[3*(el.points[2].nr-1)+2] ]))
    fde = ngm.FaceDescriptor(surfnr=1,domin=1,bc=1)
    fde.bcname = "inflow"
    fdid = ngmesh.Add(fde)
    for el in El2d:
        ngmesh.Add(ngm.Element2D(fdid, [pnums[2*(el.points[2].nr-1)]
                                       ,pnums[2*(el.points[1].nr-1)]
                                       ,pnums[2*(el.points[0].nr-1)] ]))
    fde = ngm.FaceDescriptor(surfnr=2,domin=1,bc=2)
    fde.bcname = "outflow"
    fdid = ngmesh.Add(fde)
    for el in El2d:
         ngmesh.Add(ngm.Element2D(fdid, [pnums[2*(el.points[0].nr-1)+1]
                                        ,pnums[2*(el.points[1].nr-1)+1]
                                        ,pnums[2*(el.points[2].nr-1)+1] ]))
    fde = ngm.FaceDescriptor(surfnr=3,domin=1,bc=3)
    fde.bcname = "dirichlet"
    fdid = ngmesh.Add(fde)
    for el in El1d:
        ngmesh.Add(ngm.Element2D(fdid, [pnums[2*(el.points[0].nr-1)]
                                       ,pnums[2*(el.points[1].nr-1)]
                                       ,pnums[2*(el.points[1].nr-1)+1]
                                       ,pnums[2*(el.points[0].nr-1)+1]]))
        # ngmesh.Add(ngm.Element2D(fdid, [pnums[3*(el.points[0].nr-1)+1]
        #                                 ,pnums[3*(el.points[1].nr-1)+1]
        #                                 ,pnums[3*(el.points[1].nr-1)+2]
        #                                 ,pnums[3*(el.points[0].nr-1)+2]]))


    ngmesh.SetBCName(0,"inflow")
    ngmesh.SetBCName(1,"outflow")
    ngmesh.SetBCName(2,"dirichlet")
    mesh = Mesh(ngmesh)
    return mesh



def ProdMesh(ngmeshbase,t_steps):
	ngmesh = ngm.Mesh()
	ngmesh.dim=3
	pnums = []

	for p in ngmeshbase.Points():
		x,y,z = p.p
		ngmesh.Add(ngm.MeshPoint(ngm.Pnt(x,y,0)))
	for i,p in enumerate(ngmeshbase.Points()):
		x,y,z = p.p
		pnums.append(ngmesh.Add(ngm.MeshPoint(ngm.Pnt(x,y,t_steps))))

	El1d = ngmeshbase.Elements1D()
	El2d = ngmeshbase.Elements2D()

	ngmesh.SetMaterial(1, "mat")
	for el in El2d:
		ngmesh.Add(ngm.Element3D(1, [el.points[0].nr
									,el.points[1].nr
									,el.points[2].nr
									,pnums[el.points[0].nr-1]
									,pnums[el.points[1].nr-1]
									,pnums[el.points[2].nr-1]]))
	fde = ngm.FaceDescriptor(surfnr=1,domin=1,bc=1)
	fde.bcname = "inflow"
	fdid = ngmesh.Add(fde)
	for el in El2d:
		ngmesh.Add(ngm.Element2D(fdid, [el.points[2].nr
									,el.points[1].nr
									,el.points[0].nr]))
	fde = ngm.FaceDescriptor(surfnr=2,domin=1,bc=2)
	fde.bcname = "outflow"
	fdid = ngmesh.Add(fde)
	for el in El2d:
		ngmesh.Add(ngm.Element2D(fdid, [pnums[el.points[0].nr-1]
									,pnums[el.points[1].nr-1]
									,pnums[el.points[2].nr-1]]))
	fde = ngm.FaceDescriptor(surfnr=3,domin=1,bc=3)
	fde.bcname = "dirichlet"
	fdid = ngmesh.Add(fde)
	for el in El1d:
		ngmesh.Add(ngm.Element2D(fdid, [el.points[0].nr
									,el.points[1].nr
									,pnums[el.points[1].nr-1]
									,pnums[el.points[0].nr-1]]))

	ngmesh.SetBCName(0,"inflow")
	ngmesh.SetBCName(1,"outflow")
	ngmesh.SetBCName(2,"dirichlet")
	mesh = Mesh(ngmesh)
	return mesh


def CartSquare(N,t_steps):
	ngmesh = ngm.Mesh()
	ngmesh.SetGeometry(unit_square)
	ngmesh.dim = 2
	pnums = []
	for j in range(t_steps + 1):
		for i in range(N + 1):
			pnums.append(ngmesh.Add(ngm.MeshPoint(ngm.Pnt(i / N, j / t_steps, 0))))

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
	fde.bcname = "dirichlet"
	fdid = ngmesh.Add(fde)
	for i in range(t_steps):
		ngmesh.Add(ngm.Element1D([pnums[N + i * (N + 1)], pnums[N + (i + 1) * (N + 1)]], index=3))
		ngmesh.Add(ngm.Element1D([pnums[0 + i * (N + 1)], pnums[0 + (i + 1) * (N + 1)]], index=3))


	ngmesh.SetBCName(0,"inflow")
	ngmesh.SetBCName(1,"outflow")
	ngmesh.SetBCName(2,"dirichlet")

	mesh = Mesh(ngmesh)
	# print("boundaries" + str(mesh.GetBoundaries()))
	return mesh

def LshapeMesh(maxh = 0.5):
    geo = ngeom2d.SplineGeometry()
    p1 = geo.AppendPoint (0,0)
    p2 = geo.AppendPoint (0.5,0)
    p3 = geo.AppendPoint (0.5,0.5)
    p4 = geo.AppendPoint (1,0.5)
    p5 = geo.AppendPoint (1,1)
    p6 = geo.AppendPoint (0,1)
    geo.Append (["line", p1, p2])
    geo.Append (["line", p2, p3])
    geo.Append (["line", p3, p4])
    geo.Append (["line", p4, p5])
    geo.Append (["line", p5, p6])
    geo.Append (["line", p6, p1])
    mesh = geo.GenerateMesh (maxh=maxh)
    return mesh

def CircleMesh(maxh=0.5):
    geo = ngeom2d.SplineGeometry()
    geo.AddCircle((0.5,0.5),1,bc="circle")
    ngmesh = geo.GenerateMesh(maxh=maxh)
    return ngmesh

def TunnelMesh(maxh = 0.5):
    geo = ngeom2d.SplineGeometry()
    p1 = geo.AppendPoint (-0.5,-1.25)
    p2 = geo.AppendPoint (0.5,-1.25)
    p3 = geo.AppendPoint (0.5,-0.3)
    p4 = geo.AppendPoint (0.002,-0.3)
    p5 = geo.AppendPoint (0.002,-0.2)
    p6 = geo.AppendPoint (0.5,-0.2)
    p7 = geo.AppendPoint (0.5,0.75)
    p8 = geo.AppendPoint (-0.5,0.75)
    p9 = geo.AppendPoint (-0.5,-0.2)
    p10 = geo.AppendPoint (-0.002,-0.2)
    p11 = geo.AppendPoint (-0.002,-0.3)
    p12 = geo.AppendPoint (-0.5,-0.3)
    geo.Append (["line", p1, p2])
    geo.Append (["line", p2, p3])
    geo.Append (["line", p3, p4])
    geo.Append (["line", p4, p5])
    geo.Append (["line", p5, p6])
    geo.Append (["line", p6, p7])
    geo.Append (["line", p7, p8])
    geo.Append (["line", p8, p9])
    geo.Append (["line", p9, p10])
    geo.Append (["line", p10, p11])
    geo.Append (["line", p11, p12])
    geo.Append (["line", p12, p1])
    mesh = geo.GenerateMesh (maxh=maxh)
    return mesh

def RefineAround(center,r,mesh):
    for el in mesh.Elements():
        for v in el.vertices:
            dist = 0
            p = mesh.ngmesh.Points()[v.nr+1]
            for n,pn in enumerate(p):
                dist += (pn-center[n])*(pn-center[n])
            dist = sqrt(dist)
            if(dist<r):
                mesh.SetRefinementFlag(el,1)
                break
            else:
                mesh.SetRefinementFlag(el,0)
    mesh.Refine()
    return mesh
