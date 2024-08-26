import numpy as np
import math
from FemBulk import *

#def ApplyBC_PlateWithHole(Fem, Node, Element):
def ApplyBC(Fem, Node, Element):
    #print(Node.NNode)1
    #exit(1)
    Dim    = Fem.Dimension
    for ind, x in enumerate(Node.Coord):
        Ndof = NodeDof(Node, [Node.Id[ind]])
        ## uniaxial tension ###########################
        if math.fabs(x[2] - 0.0 ) < 1e-05:          # fixed displacement at the bottom
            Node.BC_E[Ndof[0]*Dim+2] = 0.0
            #print("\t\t\t\t",Node.Id[ind])
            if math.fabs(x[0] - 0.0 ) < 1e-05:      # fixed displacement at the bottom
                Node.BC_E[Ndof[0]*Dim] = 0.0
                #print("\t\t",Node.Id[ind])
            if math.fabs(x[1] - 0.0 ) < 1e-05:      # fixed displacement at the bottom
                Node.BC_E[Ndof[0]*Dim+1] = 0.0
                #print("\t\t\t",Node.Id[ind])

        if math.fabs(x[2] - 1.0 ) < 1e-05:
            Node.BC_N[Ndof[0]*Dim+2] = 100.0/Fem.totalstep # force control
            #print("\t",Node.Id[ind])
        ##############################################
    return
