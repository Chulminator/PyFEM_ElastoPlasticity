# ---------------------------------------------------------------- 
# Written by: CK in promechanics.org at Columbia University
# ----------------------------------------------------------------       
import numpy as np
import sys
import time as ct
import os.path
import math
from numpy import ix_
from FemBulk import *

def ConstructStiffness(Fem, Node, Element):
    # Construct stiffness
    G_Edof = Fem.G_Edof
    Dim    = Fem.Dimension
    if Dim == 3:
        IxI          = np.zeros([6,6])
        IxI[0:3,0:3] = np.ones([3,3])
        II = np.diag([1., 1., 1., 0.5, 0.5, 0.5])
        D = Fem.lamda*IxI + 2.*Fem.mu*II
    elif Dim == 2:
        IxI          = np.zeros([3,3])
        IxI[0:2,0:2] = np.ones([2,2])
        II = np.diag([1., 1., 0.5])
        D = Fem.lamda*IxI + 2.*Fem.mu*II

    for ind1, Edof in enumerate(Element.Connectivity):
        ENdof = ElementalNodeDof(Fem, Node, Edof)
        E_K = np.zeros([G_Edof, G_Edof])
        B    = Element.B_matrix[ind1]
        Jacc = Element.Jacc[ind1]
        for ind, gp in enumerate(Fem.GP):
            B_GP    = B[ind]
            Jacc_GP = Jacc[ind]
            if Fem.Dimension == 2:
                GP_K = B_GP.T @ D @ B_GP * Fem.width * Jacc_GP * gp[-1]
            elif Fem.Dimension == 3:
                GP_K = B_GP.T @ D @ B_GP * Jacc_GP * gp[-1]
            E_K += GP_K
        Element.Stiffness.append(E_K)
    return

def LinearSolve(Node, Global_K):
    IndexBCN  = list(np.where(np.isnan(Node.BC_E))[0])
    IndexBCE  = list(np.where(np.invert(np.isnan(Node.BC_E)))[0])
    ValueBCE  = Node.BC_E[IndexBCE]
    
    Sliced_K12 = Global_K[ix_(IndexBCN, IndexBCE)]
    Sliced_K11 = Global_K[ix_(IndexBCN, IndexBCN)]
    
    
    BCE = ValueBCE
    BCN = Node.BC_N[IndexBCN]
    #F = BCN - np.expand_dims(Sliced_K12 @ BCE, axis=1)?
    F = BCN - Sliced_K12 @ BCE
    #print(Sliced_K12 @ BCE)
    
    Node.u[IndexBCE] = ValueBCE
    Node.u[IndexBCN] = np.linalg.solve(Sliced_K11, F)
    return
