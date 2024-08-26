 
import numpy as np
import sys
sys.path.append("../include")
import pygmsh
import meshio
from WriteMesh import *

#https://github.com/nschloe/meshio

## Example2
## Circle
with pygmsh.geo.Geometry() as geom:
    geom.add_polygon(
        [
            [0.0, 0.0],
            [1.0, 0.0],
            [1.0, 1.0],
            [0.0, 1.0],
        ],
        mesh_size=0.3,
    )
    mesh = geom.generate_mesh()
    pygmsh.write("./Patch2D.msh")
mesh.write("./Patch2D.vtk")
print(mesh)
msh = meshio.read("./Patch2D.msh")
Node = msh.points

#for ii in msh.cells:
    #print(ii,"\t",ii.data)
#exit(1)
Element = msh.cells[-1].data
#print(msh.cells)
#print(Node.shape)
#print(Element.shape)

WriteNodeAndElementFromSandiaMesh("./Patch2D", Node, Element)
