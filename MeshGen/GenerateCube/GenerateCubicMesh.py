 
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
    #add_box(0, 1, 0, 1, 0, 1, 0.05)
    geom.add_box(0, 5, 0, 10, 0, 30, mesh_size=1)
    mesh = geom.generate_mesh()
    pygmsh.write("./cubic.msh")
print(mesh)
exit(1)
mesh.write("./cubic.vtk")
