# ---------------------------------------------------------------- 
# Written by: Hyoung Suk Suh (h.suh@columbia.edu)     
# sub by: Chulmin Kweon      (Chulmin.Kweon@columbia.edu)
# ----------------------------------------------------------------       

import numpy as np
import sys
sys.path.append("./include")
import time as ct

import Read
import WriteResult
import LinearElasticity as LinEla
import FemBulk as Bulk

tic = ct.time()
# ---------------------------------------------------------------- 
# Input parameters
# ----------------------------------------------------------------
Fem = Read.Model()
Node, Element = Read.ReadInput("./input/PlateWithHole/PlateWithHole_113_93_elastic.inp",Fem)
Fem.G_Edof = Fem.Dimension* len(Element.Connectivity[0])

LinEla.ConstructStiffness(Fem, Node, Element)
#LinEla.ApplyBC(Fem, Node, Element)
Bulk.ApplyBC_PlateWithHole(Fem, Node, Element)
Global_K = Bulk.Assembleage(Fem, Node, Element)
LinEla.LinearSolve(Node, Global_K)
#Node.u = np.squeeze(Node.u)
#print(Node.u)
LinEla.CalculateStrainAndStressAtGP(Fem, Node, Element)
Bulk.CalculateStressAtNode(Fem, Node, Element)
Bulk.CalculateVonMisesStressAtNode(Fem, Node)
Bulk.CalculatePrincipalStressAtNode(Node)

Attribute = np.zeros((Node.NNode,3))
Attribute[:,0] = Node.stress[:,0]
Attribute[:,1] = Node.sigma1
Attribute[:,2] = Node.sigmaVM
WriteResult.WriteCustomizedAttribute(Fem, Node, Element, Attribute)
WriteResult.WriteDisp(Fem, Node, Element)
WriteResult.WriteStress(Fem, Node, Element)
toc = ct.time() - tic
print(' ')
print('Elapsed CPU time: ', toc, '[sec]')

