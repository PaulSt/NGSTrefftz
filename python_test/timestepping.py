import netgen.gui
import tkinter
from math import pi
from ngsolve import *
from netgen.geom2d import SplineGeometry
import numpy as np

from prodmesh import *
from netgen.geom2d import unit_square
from netgen.meshing import MeshingParameters

mesh = Mesh(unit_square.GenerateMesh(maxh=0.08))
maxh=0.01
minh=0.01
mp = MeshingParameters (maxh = maxh)
refpoints = 500
for i in range(0, refpoints+1):
    for j in range(0, refpoints+1):
        xk = i/refpoints
        yk = j/refpoints
        mp.RestrictH (x=xk, y=yk, z=0, h=max(minh,sqrt(0.005*((xk-0.5)*(xk-0.5)+(yk-0.5)*(yk-0.5)))))

mesh = Mesh( LshapeMesh(maxh,mp) )

fes = H1(mesh, order=2) #, dirichlet="bottom|right|left|top")

u,v = fes.TnT() # TnT : Trial and Test function

dt = 0.001

a = BilinearForm(fes)
a += SymbolicBFI (grad(u)*grad(v))
a.Assemble()

m = BilinearForm(fes)
m += SymbolicBFI (u*v)
m.Assemble()
# invm = m.mat.Inverse(freedofs=fes.FreeDofs())
invm = m.mat.Inverse()

u0 = GridFunction(fes)
u1 = GridFunction(fes)
u2 = GridFunction(fes)
uD = GridFunction(fes)

delta = 800
u0.Set( exp(-delta*((x-0.25)*(x-0.25))) )
u1.Set( exp(-delta*((x-0.25)*(x-0.25))) )
# u0.Set( exp(-delta*(((x-0.5)*(x-0.5))+((y-0.5)*(y-0.5))) ))
# u1.Set( exp(-delta*(((x-0.5)*(x-0.5))+((y-0.5)*(y-0.5))) ))
# sq = sqrt(0.5);
# u0.Set(sin( 0+sq*(x+y) ))
# u1.Set(sin( dt+sq*(x+y) ))
Draw(u1,mesh,"u")
input()

# time =  0.0
time = dt
res = u0.vec.CreateVector()
while time < 1:
    time += dt
    # uD.Set(sin(time+sq*(x+y)),BND)
    res.data = dt*dt * (- a.mat * u1.vec  )
    u2.vec.data = 2*u1.vec - u0.vec + invm * res
    # for i,t in enumerate(fes.FreeDofs()):
        # if t==False:
            # u2.vec[i] = uD.vec[i]
    u0.vec.data = u1.vec.data
    u1.vec.data = u2.vec.data

    if(int(1000*time)%100==0):
        Redraw(blocking=True)
    print("\r",time,"   ",time%0.1,end="")
