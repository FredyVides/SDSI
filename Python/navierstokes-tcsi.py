from ngsolve import *

# Program based on program navierstokes.py
# Code by Joachim Schoeberl Netgen/NGSolve

# Example: netgen navierstokes-tcsi.py

# default values nu = 0.001, Maxh=0.02,tau=0.0001

# viscosity
nu = 0.001

# timestepping parameters
tau = 0.0001
tend = 10
Maxh=0.02

from netgen.geom2d import SplineGeometry
geo = SplineGeometry()
geo.AddRectangle( (0, 0), (2, 0.41), bcs = ("wall", "outlet", "wall", "inlet"))
geo.AddRectangle ( (0.16, 0.16), (0.24,0.24), leftdomain=0, rightdomain=1, bcs=("cyl","cyl","cyl","cyl"), maxh=Maxh/2)
mesh = Mesh( geo.GenerateMesh(maxh=Maxh))

mesh.Curve(3)

V = VectorH1(mesh,order=3, dirichlet="wall|cyl|inlet")
Q = H1(mesh,order=2)

X = FESpace([V,Q])

u,p = X.TrialFunction()
v,q = X.TestFunction()

stokes = nu*InnerProduct(grad(u), grad(v))+div(u)*q+div(v)*p - 1e-10*p*q
a = BilinearForm(X)
a += stokes*dx
a.Assemble()

# nothing here ...
f = LinearForm(X)   
f.Assemble()

# gridfunction for the solution
gfu = GridFunction(X)

# parabolic inflow at inlet:
uin = CoefficientFunction( (1.5*4*y*(0.41-y)/(0.41*0.41), 0) )
gfu.components[0].Set(uin, definedon=mesh.Boundaries("inlet"))

# solve Stokes problem for initial conditions:
inv_stokes = a.mat.Inverse(X.FreeDofs())

res = f.vec.CreateVector()
res.data = f.vec - a.mat*gfu.vec
gfu.vec.data += inv_stokes * res


# matrix for implicit Euler 
mstar = BilinearForm(X)
mstar += SymbolicBFI(u*v + tau*stokes)
mstar.Assemble()
inv = mstar.mat.Inverse(X.FreeDofs(), inverse="sparsecholesky")

# the non-linear term 
conv = BilinearForm(X, nonassemble = True)
conv += (grad(u) * u) * v * dx

# for visualization
Draw (Norm(gfu.components[0]), mesh, "velocity", sd=3)

# implicit Euler/explicit Euler splitting method:
t = 0

import numpy as np

GFUdata=np.array([gfu.vec.data])

counter = 1

with TaskManager():
    while t < tend:
        print ("t=", t, end="\r")

        conv.Apply (gfu.vec, res)
        res.data += a.mat*gfu.vec
        gfu.vec.data -= tau * inv * res 
        
        
        if counter%200==1: GFUdata=np.append(GFUdata,[gfu.vec.data],axis=0)

        t = t + tau
        counter=counter+1
        
        Redraw()

np.savetxt("GFUdata.csv", GFUdata.T, delimiter=",")
