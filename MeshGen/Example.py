 
import numpy as np
import sys
#sys.path.append("./include")
import pygmsh
import meshio
from include.WriteMesh import *
#import include.WriteMesh as Write

#https://github.com/nschloe/meshio

## Example1
## simple geometry
with pygmsh.geo.Geometry() as geom:
    geom.add_polygon(
        [
            [0.0, 0.0],
            [1.0, -0.2],
            [1.1, 1.2],
            [0.1, 0.7],
        ],
        mesh_size=0.1,
    )
    mesh = geom.generate_mesh()
    pygmsh.write("./Result/Example1.msh")
mesh.write("./Result/Example1.vtk")
# *.msh file is for gmsh
# *.vtk file is for paraview


## Example2
## Circle
with pygmsh.geo.Geometry() as geom:
    geom.add_circle([0.0, 0.0], 1.0, mesh_size=0.2)
    mesh = geom.generate_mesh()
    pygmsh.write("./Result/Example2.msh")
mesh.write("./Result/Example2.vtk")


## Example3
## input mesh in the meshio
# two triangles and one quad
points = [
    [0.0, 0.0],
    [1.0, 0.0],
    [0.0, 1.0],
    [1.0, 1.0],
    [2.0, 0.0],
    [2.0, 1.0],
]
cells = [
    ("triangle", [[0, 1, 2], [1, 3, 2]]),
    ("quad", [[1, 4, 5, 3]]),
]

mesh = meshio.Mesh(
    points,
    cells,
    # Optionally provide extra data on points, cells, etc.
    point_data={"T": [0.3, -1.2, 0.5, 0.7, 0.0, -3.0]},
    # Each item in cell data must match the cells array
    cell_data={"a": [[0.1, 0.2], [0.4]]},
)
mesh.write(
    "./Result/Example3.vtk",  # str, os.PathLike, or buffer/open file
    # file_format="vtk",  # optional if first argument is a path; inferred from extension
)

# Alternative with the same options
#meshio.write_points_cells("./Result/Example3.vtk", points, cells)

#exit(1)
## Example4
## Write node and element from Meshio
directory = "./Result/Example3"
MeshType  = "q4"
WriteNodeAndElementFromMeshio(mesh, directory, MeshType)
exit(1)
