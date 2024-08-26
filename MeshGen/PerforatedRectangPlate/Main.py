
import numpy as np
import sys
sys.path.append("../include")
import pygmsh
import meshio
from WriteMesh import *


##
msh = meshio.read("./Mesh/QuadMesh1.msh")
print(msh)

Node = msh.points
Element = msh.cells[0].data

WriteNodeAndElementFromSandiaMesh("./QuadMesh1", Node, Element)

##
msh = meshio.read("./Mesh/QuadMesh2.msh")
print(msh)

Node = msh.points
Element = msh.cells[0].data

WriteNodeAndElementFromSandiaMesh("./QuadMesh2", Node, Element)

##
msh = meshio.read("./Mesh/QuadMesh3.msh")
print(msh)

Node = msh.points
Element = msh.cells[0].data

WriteNodeAndElementFromSandiaMesh("./QuadMesh3", Node, Element)

##
msh = meshio.read("./Mesh/QuadMesh4.msh")
print(msh)

Node = msh.points
Element = msh.cells[0].data

WriteNodeAndElementFromSandiaMesh("./QuadMesh4", Node, Element)

##
msh = meshio.read("./Mesh/TriMesh1.msh")
print(msh)

Node = msh.points
Element = msh.cells[0].data

WriteNodeAndElementFromSandiaMesh("./TriMesh1", Node, Element)

##
msh = meshio.read("./Mesh/TriMesh2.msh")
print(msh)

Node = msh.points
Element = msh.cells[0].data

WriteNodeAndElementFromSandiaMesh("./TriMesh2", Node, Element)

##
msh = meshio.read("./Mesh/TriMesh3.msh")
print(msh)

Node = msh.points
Element = msh.cells[0].data

WriteNodeAndElementFromSandiaMesh("./TriMesh3", Node, Element)

