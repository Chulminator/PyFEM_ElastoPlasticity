
import numpy as np
import sys
sys.path.append("../include")
import pygmsh
import meshio
from WriteMesh import *


##
msh = meshio.read("./Mesh/HSMesh3.msh")
print(msh)

Node = msh.points
Element = msh.cells[0].data

WriteNodeAndElementFromSandiaMesh("./HSMesh3", Node, Element)

