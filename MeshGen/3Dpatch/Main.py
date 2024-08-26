
import numpy as np
import sys
sys.path.append("../include")
import pygmsh
import meshio
from WriteMesh import *


##
msh = meshio.read("./Mesh/1111111.msh")
print(msh)

Node = msh.points
Element = msh.cells[0].data

WriteNodeAndElementFromSandiaMesh("./1111111", Node, Element,3)

