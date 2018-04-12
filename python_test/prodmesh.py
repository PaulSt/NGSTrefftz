import netgen.meshing as ngm
from ngsolve import *
from netgen.geom2d import unit_square

def ProdMesh(ngmeshbase,t_steps):
	ngmesh = ngm.Mesh()
	ngmesh.dim=3
	pnums = []

	for p in ngmeshbase.Points():
		x,y,z = p.p
		ngmesh.Add(ngm.MeshPoint(ngm.Pnt(x,y,0)))
	for p in ngmeshbase.Points():
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

	ngmesh.Add(ngm.FaceDescriptor(surfnr=1,domin=1,bc=1))
	for el in El2d:
		ngmesh.Add(ngm.Element2D(1, [el.points[2].nr
									,el.points[1].nr
									,el.points[0].nr]))
		ngmesh.Add(ngm.Element2D(1, [pnums[el.points[0].nr-1]
									,pnums[el.points[1].nr-1]
									,pnums[el.points[2].nr-1]]))
	for el in El1d:
		ngmesh.Add(ngm.Element2D(1, [el.points[0].nr
									,el.points[1].nr
									,pnums[el.points[1].nr-1]
									,pnums[el.points[0].nr-1]]))

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
	for i in range(t_steps):
		ngmesh.Add(ngm.Element1D([pnums[N + i * (N + 1)], pnums[N + (i + 1) * (N + 1)]], index=1))
		ngmesh.Add(ngm.Element1D([pnums[0 + i * (N + 1)], pnums[0 + (i + 1) * (N + 1)]], index=1))
	for i in range(N):
		ngmesh.Add(ngm.Element1D([pnums[i], pnums[i + 1]], index=2))
		ngmesh.Add(ngm.Element1D([pnums[i + t_steps * (N + 1)], pnums[i + 1 + t_steps * (N + 1)]], index=2))

	mesh = Mesh(ngmesh)
	# print("boundaries" + str(mesh.GetBoundaries()))
	return mesh
