# ---------------------------------------------------------------- 
# Written by: CK in promechanics.org at Columbia University
# ----------------------------------------------------------------       
import numpy as np
import sys
import time as ct
import os.path
import math
from numpy import ix_
from LinearElasticity import *

def CalculateVonMisesStressAtNode(Fem, Node):
    # Calculate Von Mises stress at node
    # the stress at nodes should be calculated first
    # https://www.comsol.com/blogs/what-is-the-difference-between-plane-stress-and-plane-strain
    if Fem.twoD == 'planestrain':
        x  = Node.stress[:,0]
        y  = Node.stress[:,1]
        z  = Fem.nu*(x+y)
        xy = Node.stress[:,2]
        yz = np.zeros((Node.NNode))
        zx = np.zeros((Node.NNode))
    elif Fem.twoD == 'planestress':
        x  = Node.stress[:,0]
        y  = Node.stress[:,1]
        z  = np.zeros((Node.NNode))
        xy = Node.stress[:,2]
        yz = np.zeros((Node.NNode))
        zx = np.zeros((Node.NNode))
    else:
        #print("Plane Stress is not ready")
        assert(0, "Plane Stress or plane strain")
    Node.sigmaVM = np.sqrt(0.5*((x-y)**2 + (y-z)**2 + (z-x)**2 + 6.0*(xy**2 + yz**2 + zx**2)))
    return


def CalculatePrincipalStressAtNode(Node):
    # Calculate Von Mises stress at node
    # the stress at nodes should be calculated first
    Node.sigma1 = 0.5*(Node.stress[:,0] + Node.stress[:,1]) + np.sqrt(0.25*(Node.stress[:,0]-Node.stress[:,1])**2 + Node.stress[:,2]**2)
    return



def CalculateStressAtNode(Fem, Node, Element):
    # Calculate stress at nodes
    # the stress at gauss points should be calculated first
    # this should be done in plasticity "ConstructStiffness"
    G_Edof = Fem.G_Edof
    Dim    = Fem.Dimension
    Node.stress = np.zeros([Node.NNode,3])
    LocalNElem  = np.zeros([Node.NNode])

    GP_SR = []
    for gp in Fem.GP:
        N, _ = ShapeFunction(1./gp[0], 1./gp[1])
        GP_SR.append(N)
    GP_SR = np.array(GP_SR)
    #print(GP_SR)
    for ind1, Edof in enumerate(Element.Connectivity):

        ENdof = ElementalNodeDof(Node, Edof)
        Ndof  = NodeDof(Node, Edof)

        B = Element.B_matrix[ind1]
        ElementGPstress = np.array(Element.GPstress[ind1])
        #print( ElementGPstress )
        #print(GP_SR @ ElementGPstress)
        #print(Node.stress[Ndof,:])
        Node.stress[Ndof,:] += GP_SR @ ElementGPstress
        LocalNElem[Ndof] += 1
    Node.stress[:,0] = Node.stress[:,0]/LocalNElem
    Node.stress[:,1] = Node.stress[:,2]/LocalNElem
    Node.stress[:,2] = Node.stress[:,2]/LocalNElem

def Assembleage(Fem, Node, Element):
    # Stiffness matrix is assembled in Global scale
    Global_K = np.zeros([Node.NNode * Fem.Dimension, Node.NNode * Fem.Dimension])
    for EID, Edof in enumerate(Element.Connectivity):
        ENdof = ElementalNodeDof(Node, Edof)
        for Li, Gi in enumerate(ENdof):
            for Lj, Gj in enumerate(ENdof):
                Global_K[Gi,Gj] += Element.Stiffness[EID][Li,Lj]
        
    return Global_K
    
def ApplyBC(Fem, Node, Element):
    # Apply boundary condition
    ################### have to be edited
    G_Edof = Fem.G_Edof
    Dim    = Fem.Dimension

    EBC =[]
    for ind,x in enumerate(Node.Coord):
        if math.fabs(x[0] - 0.0 ) < 1e-05:
            Node.BC_E[Node.Id[ind]*2-2] = 0
            Node.BC_E[Node.Id[ind]*2-1]   = 0
        if math.fabs(x[0] - 2.0 ) < 1e-05 and math.fabs(x[1] - 0.05 ) < 1e-05 :
            #print(Node.Id[ind])
            Node.BC_N[Node.Id[ind]*2-1] = (3.*Fem.E*0.1**4./12.)/8.
    return

def ApplyBC_PlateWithHole(Fem, Node, Element):
    EBC =[]
    Dim    = Fem.Dimension
    for ind,x in enumerate(Node.Coord):
        Ndof = NodeDof(Node, [Node.Id[ind]])
        if math.fabs(x[0] - 0.0 ) < 1e-05:
            Node.BC_E[Ndof[0]*2] = 0
        if math.fabs(x[1] - 0.0 ) < 1e-05:
            Node.BC_E[Ndof[0]*2+1] = 0.0
        if math.fabs(x[0] - 10.0 ) < 1e-05:
            Node.BC_E[Ndof[0]*2] = 0.02/Fem.totalstep
    return

def ElementalNodeDof(Node, Connectivity): ################################## Generalization check
    tmp = NodeDof(Node, Connectivity)
    tmp1 = np.array(tmp)*2
    tmp2 = np.array(tmp)*2 + 1

    c = np.empty((tmp1.size + tmp2.size), dtype =tmp1.dtype)
    c[0::2] = tmp1
    c[1::2] = tmp2

    return c

def NodeDof(Node, Connectivity):

    Ndof1 =[]
    for ii in Connectivity:
        Ndof1.append(Node.Id.index(ii))
    return Ndof1

def ShapeFunction(s,t):
    N_matrix = np.array([(s-1.)*(t-1.)/4.,\
                        (-s-1.)*(t-1.)/4.,\
                        (-s-1.)*(-t-1.)/4.,\
                        (s-1.)*(-t-1.)/4.])
    
    dN   = np.array([[t-1., 1.-t, 1.+t, -1.-t],
                    [s-1., -1.-s, 1.+s, 1.-s]])
    dN *= 0.25
    return N_matrix, dN 

def Strain(Fem, NodeTarget,s,t):
    G_Edof = Fem.G_Edof
    Dim    = Fem.Dimension
    _, dN = ShapeFunction(s,t)

    Jacobian = np.matmul(dN, NodeTarget)

    B1= np.array([[1., 0., 0., 0.],
                  [0., 0., 0., 1],
                  [0., 1., 1., 0]])

    B2 = np.zeros([4,4])
    
    B2[2:4,2:4] = np.linalg.inv(Jacobian)
    B2[0:2,0:2] = np.linalg.inv(Jacobian)
    
    B3 = np.zeros([4,G_Edof])
    
    for ind in range(int(G_Edof/Dim)):
        B3[0:2,2*ind]   = dN[:,ind]
        B3[2:4,2*ind+1] = dN[:,ind]
        
    B_matrix= B1 @ B2 @ B3
    Jacc = np.linalg.det(Jacobian)
    return Jacc, B_matrix
