import numpy as np
import math
from FemBulk import *
from MeshHandle import *

def ApplyBC(Fem, Node, Element):
    Dim    = Fem.Dimension
    Facet = GenerateFacet(Fem, Element)


    ######### gauss point should be different for the different elemeny type #########
    # ex
    # t3, q4 -> 1d 2point
    # t6 -> 1d 3point
    # etc
    #
    # now only t3, q4 are prepared
    ######### gauss point should be different for the different elemeny type #########

    GP1d = [[-np.sqrt(1./3), 1],
            [np.sqrt(1./3), 1]]

    ApplyingTraction = 100
    for FacetNode in Facet.AdjacNode:
        Ndof = NodeDof(Node, FacetNode)
        x1 = Node.Coord[Ndof[0]][0]
        y1 = Node.Coord[Ndof[0]][1]
        x2 = Node.Coord[Ndof[1]][0]
        y2 = Node.Coord[Ndof[1]][1]
        if math.fabs(y1 - 1.0 ) < 1e-05 and math.fabs(y2 - 1.0 ) < 1e-05:
            length = np.sqrt((x1-x2)**2+(y1-y2)**2)
            for ii, gp in enumerate(GP1d):
                s      = gp[0]
                weight = gp[1]
                N = np.zeros(2)
                N[0] = (1.0 - s)*0.5
                N[1] = (1.0 + s)*0.5
                Node.BC_N[Ndof[0]*Dim+1] += N[0] * ApplyingTraction * weight * Fem.width * (length*0.5) /Fem.totalstep
                Node.BC_N[Ndof[1]*Dim+1] += N[1] * ApplyingTraction * weight * Fem.width * (length*0.5) /Fem.totalstep



    for ind, x in enumerate(Node.Coord):
        Ndof = NodeDof(Node, [Node.Id[ind]])
        ## uniaxial tension ###########################
        if math.fabs(x[1] - 0.0 ) < 1e-05:
            Node.BC_E[Ndof[0]*Dim+1] = 0.0
            if math.fabs(x[0] - 0.0 ) < 1e-05:
                Node.BC_E[Ndof[0]*Dim] = 0.0
        #if math.fabs(x[1] - 1.0 ) < 1e-05:
            #Node.BC_N[Ndof[0]*Dim+1] = 100.0/Fem.totalstep # force control
        ###############################################




        #if math.fabs(x[1] - 0.0 ) < 1e-05:
            #Node.BC_E[Ndof[0]*2+1] = 0.0
        #if math.fabs(x[0] - 0.0 ) < 1e-05:
            #Node.BC_E[Ndof[0]*2] = 0.0
        #if math.fabs(x[1] - 1.0 ) < 1e-05:
            #Node.BC_E[Ndof[0]*2] = 0.03/Fem.totalstep
    return
