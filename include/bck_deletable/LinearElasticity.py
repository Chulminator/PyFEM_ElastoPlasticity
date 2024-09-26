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

def CalculateStrainAndStressAtGP(Fem, Node, Element):
    # Calculate strain and stress at gauss point
    D = ConstitutiveLaw(Fem)
    Element.GPstrain = []
    Element.GPstress = []
    for ind1, Edof in enumerate(Element.Connectivity):
        ENdof = ElementalNodeDof(Node, Edof)
        ElemStrain = []
        ElemStress = []
        B = Element.B_matrix[ind1]
        for ind, gp in enumerate(Fem.GP):
            B_GP = B[ind]
            strain = B_GP @ Node.u[ENdof]
            ElemStrain.append(strain)
            ElemStress.append(D @ strain)
        Element.GPstrain.append(ElemStrain)
        Element.GPstress.append(ElemStress)

    return

def ConstructStiffness(Fem, Node, Element):
    # Construct stiffness
    G_Edof = Fem.G_Edof
    Dim    = Fem.Dimension
    D = ConstitutiveLaw(Fem)
    Element.B_matrix =[]
    Element.Jacc  =[]
    Element.Stiffness  =[]
    Element.Area  =[]
    for Edof in Element.Connectivity:
        #ENdof = ElementalNodeDof(Edof)
        NodeTarget = []
        for dof in Edof:
            NodeTarget.append(Node.Coord[Node.Id.index(dof)])
        E_K = np.zeros([G_Edof, G_Edof])
        Area = 0
        B_matrix = []
        Jacc_Elem = []
        for gp in Fem.GP:
            Jacc, B = Strain(Fem, NodeTarget, gp[0], gp[1])
            GP_K = B.T @ D @ B * Fem.width * Jacc * gp[-1]
            E_K += GP_K
            Area = Area + Jacc * gp[-1]
            B_matrix.append(B)
            Jacc_Elem.append(Jacc)
        Element.B_matrix.append(B_matrix)
        Element.Jacc.append(Jacc_Elem)
        Element.Stiffness.append(E_K)
        Element.Area.append(Area)
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

def ConstitutiveLaw(Fem):
    #D_matrix = np.zeros([3,3])
    tmp1 = np.zeros([3,3])
    tmp1[0,0] = 1.
    tmp1[0,1] = 1.
    tmp1[1,0] = 1.
    tmp1[1,1] = 1.
    
    tmp2 = np.zeros([3,3])
    tmp2[0,0] = 2.
    tmp2[1,1] = 2.
    tmp2[2,2] = 1.
    
    return tmp1 *Fem.lamda + tmp2 * Fem.mu
