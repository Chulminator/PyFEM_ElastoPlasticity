import numpy as np
import math
from FemBulk import *
from MeshHandle import *

def ApplyBC(Model):
    Dim      = Fem.Dimension
    for ind, x in enumerate(Node.Coord):
        Ndof = NodeDof(Node, [Node.Id[ind]])
        if math.fabs(x[0] - 0.0 ) < 1e-05:  # Left boundary
            Node.BC_E[Ndof[0]*Dim] = 0.0
        if math.fabs(x[1] - 0.0 ) < 1e-05:  # Bottom boundary
            Node.BC_E[Ndof[0]*Dim+1] = 0.0

    #Facet = GenerateFacet(Fem, Element)
    Facet = Fem.Facet
    if Fem.ElmentType == 'q4' or Fem.ElmentType == 't3':
        GP1d = [[-np.sqrt(1./3), 1],
                [np.sqrt(1./3), 1]]
    else:
        assert False, "Boundary condition is not ready for the "+Fem.ElmentType+"Mesh"


    Node.BC_N      = np.zeros([Node.NNode*Fem.Dimension])
    ApplyingTraction = 100.
    for FacetNode in Facet.AdjacNode:
        Ndof = NodeDof(Node, FacetNode)
        x1 = Node.Coord[Ndof[0]][0]
        y1 = Node.Coord[Ndof[0]][1]
        x2 = Node.Coord[Ndof[1]][0]
        y2 = Node.Coord[Ndof[1]][1]

        if math.fabs(x1 - 10.0 ) < 1e-05 and math.fabs(x2 - 10.0 ) < 1e-05:
            length = np.sqrt((x1-x2)**2+(y1-y2)**2)
            for ii, gp in enumerate(GP1d):
                s      = gp[0]
                weight = gp[1]
                N = np.zeros(2)
                N[0] = (1.0 - s)*0.5
                N[1] = (1.0 + s)*0.5
                Node.BC_N[Ndof[0]*Dim] += N[0] * ApplyingTraction * weight * Fem.width * (length*0.5)/Fem.totalstep
                Node.BC_N[Ndof[1]*Dim] += N[1] * ApplyingTraction * weight * Fem.width * (length*0.5)/Fem.totalstep
                # (length*0.5) is jacobian
    return

