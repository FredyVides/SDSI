from pandas import read_csv
from numpy import linspace
import sys

GFU=read_csv('../DataSets/GFUdata.csv',header=None,index_col=0)
GFU=GFU.values
N=GFU.shape[1]

t=linspace(0, 10,N)
Maxh=0.02

from ngsolve import *
from netgen import gui

from netgen.geom2d import SplineGeometry
geo = SplineGeometry()
geo.AddRectangle( (0, 0), (2, 0.41), bcs = ("wall", "outlet", "wall", "inlet"))
geo.AddRectangle ( (0.16, 0.16), (0.24,0.24), leftdomain=0, rightdomain=1, bcs=("cyl","cyl","cyl","cyl"), maxh=Maxh/2)
mesh = Mesh( geo.GenerateMesh(maxh=Maxh))

mesh.Curve(3)

V = VectorH1(mesh,order=3, dirichlet="wall|cyl|inlet")
Q = H1(mesh,order=2)

X = FESpace([V,Q])

# gridfunction for the visualization
gfu = GridFunction(X)

# for visualization
gfu.vec.data = GFU[:,0]
Draw (Norm(gfu.components[0]), mesh, "velocity", sd=4)

counter = 0

while counter <= (N-2):
        counter=counter+1
        print ("t=", t[counter], end="\r")        
        gfu.vec.data = GFU[:,counter]
        Redraw(blocking=True)
