import numpy as np
import math
from FemBulk import *
from MeshHandler import *

def ApplyBC(Model):
    Dim    = Model.Dim
    for ind, node in enumerate(Model.Node):
        if math.fabs(node.Coord[0] - 0.0 ) < 1e-05:
            node.BC_E[0] = 0.0

        if math.fabs(node.Coord[1] - 0.0 ) < 1e-05:
            node.BC_E[1] = 0.0

    #Facet = GenerateFacet(Fem, Element)
    if Model.Facet == None:
        Model.Facet = GenerateFacet(Model)
        Model.NFacet  = len(Model.Facet)

    if Model.Edge == None:
        Model.Edge = GenerateEdge(Model)
        Model.NEdge   = len(Model.Edge)

    GP1d = [[-np.sqrt(1./3), 1],
            [np.sqrt(1./3), 1]]

    ApplyingTraction = 600.
    for facet in Model.Facet:
        FacetNode = facet.AdjNode
        #Ndof = Model.GetNodeDof(FacetNode)
        x1 = FacetNode[0].Coord[0]
        y1 = FacetNode[0].Coord[1]
        x2 = FacetNode[1].Coord[0]
        y2 = FacetNode[1].Coord[1]

        if math.fabs(y1 - 1.0 ) < 1e-05 and math.fabs(y2 - 1.0 ) < 1e-05:
            length = np.sqrt((x1-x2)**2+(y1-y2)**2)
            for ii, gp in enumerate(GP1d):
                s      = gp[0]
                weight = gp[1]
                N = np.zeros(2)
                N[0] = (1.0 - s)*0.5
                N[1] = (1.0 + s)*0.5
                FacetNode[0].BC_N[1] += N[0] * ApplyingTraction * weight * Model.width * (length*0.5) / Model.totalstep
                FacetNode[1].BC_N[1] += N[1] * ApplyingTraction * weight * Model.width * (length*0.5) / Model.totalstep
    #print(str(ind)+"="*40)
    return
