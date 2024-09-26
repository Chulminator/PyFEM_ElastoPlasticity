# ----------------------------------------------------------------
# Written by: CK in promechanics.org at Columbia University
# ----------------------------------------------------------------       
import numpy as np
import sys
import time as ct
import os.path
import math
from numpy import ix_
#from DataStructure import *


def InitialHSP(Fem, Node, Element):
    G_Edof = Fem.G_Edof
    Dim    = Fem.Dimension
    #IndexBCE  = list(np.where(np.invert(np.isnan(Node.BC_E)))[0])
    Node.F_HSP = np.zeros_like(Node.F_int)
    for ind1, Edof in enumerate(Element.Connectivity):
        ENdof = ElementalNodeDof(Fem, Node, Edof)
        B    = Element.B_matrix[ind1]
        Jacc = Element.Jacc[ind1]
        GPF_HSP = np.zeros((len(ENdof)))
        for ind, gp in enumerate(Fem.GP):
            B_GP    = B[ind]
            Jacc_GP = Jacc[ind]
            if Fem.Dimension == 2:
                HSP_stress = np.array([Fem.HSP, Fem.HSP, 0.])
                GPF_HSP += B_GP.T @ HSP_stress * Fem.width * Jacc_GP * gp[-1]
            elif Fem.Dimension == 3:
                HSP_stress = np.array([-Fem.HSP, -Fem.HSP, -Fem.HSP, 0., 0., 0.])
                GPF_HSP += B_GP.T @ HSP_stress * Jacc_GP * gp[-1]
            else:
                assert False, "Dimension error"
        Node.F_HSP[ENdof] += GPF_HSP
    return






#from LinearElasticity import *
def GetGPSR(Fem):
    GP_SR = []
    Dim = Fem.Dim
    if Dim == 2:
        if Fem.ElmentType == "q4":
            for gp in Fem.GP:
                N, _ = ShapeFunction(Fem,1./gp[0], 1./gp[1])
                GP_SR.append(N)
            GP_SR = np.array(GP_SR)
        elif Fem.ElmentType == "t3":
            N = [[1.],
                 [1.],
                 [1.]]
            GP_SR = np.array(N)
            #print(GP_SR)
            #assert 0, "T3 has to be checked one more time "

        elif Fem.ElmentType == "t6":
            s = [-1/3,  5/3, -1/3,  2/3,  2/3, -1/3]
            t = [-1/3, -1/3,  5/3, -1/3,  2/3,  2/3]
            for ii in range(6):
                N= np.array([s[ii], t[ii], 1.-s[ii]-t[ii]])
                GP_SR.append(N)
            GP_SR = np.array(GP_SR)
        else:
            print("Element type: ", Fem.ElmentType)
            assert False, "Element type is not ready"

    elif Dim ==3:
        if Fem.ElmentType == "hex8":
            for gp in Fem.GP:
                N, _ = ShapeFunction(Fem,1./gp[0], 1./gp[1], 1./gp[2])
                GP_SR.append(N)
            GP_SR = np.array(GP_SR)

        elif Fem.ElmentType == "tet4":
            N = [[1.],
                 [1.],
                 [1.],
                 [1.]]
            GP_SR = np.array(N)

        else:
            assert False, "Element type is not ready"

    else:
        assert False, "Dimension Error"

    return GP_SR


def CalculatePlasticMultiplierAtNode(Model):
    Dim               = Model.Dim
    LocalNElem        = np.zeros([Model.NNode])
    PlasticMultiplier = np.zeros((Model.NNode))

    for ii, element in enumerate(Model.Element):
        GP_SR  = GetGPSR(element)
        Ndof  = Model.GetNodeDof(element.Connectivity)

        ElementGPlamda1 = np.array(element.GPlamda)

        PlasticMultiplier[Ndof] += GP_SR @ ElementGPlamda1

        LocalNElem[Ndof] += 1

    PlasticMultiplier /= LocalNElem
    Model.NodalPlasticMultiplier = PlasticMultiplier



def CalculateNodalDisp(Model):
    Dim = Model.Dim
    Nodalu = np.zeros((Model.NNode*Dim))
    for ii, node in enumerate(Model.Node):
        for jj in range(Dim):
            Nodalu[ii*Dim+jj] = node.u[jj]
    Model.Nodalu = Nodalu


def CalculateVonMisesStressAtNode(Model):
    # Calculate Von Mises stress at node
    # the stress at nodes should be calculated first
    # https://www.comsol.com/blogs/what-is-the-difference-between-plane-stress-and-plane-strain
    x  = np.zeros((Model.NNode))
    y  = np.zeros((Model.NNode))
    z  = np.zeros((Model.NNode))+Model.HSP
    xy = np.zeros((Model.NNode))
    yz = np.zeros((Model.NNode))
    zx = np.zeros((Model.NNode))
    #for ii, node in enumerate(Model.Node):
    if Model.Dim == 2:
        if Model.twoD == 'planestrain':
            # This is not appropriate for the heterogenous material
            lam = Model.lamda
            mu  = Model.mu
            x   = Model.Nodalstress[:,0]
            y   = Model.Nodalstress[:,1]
            z   = lam*(Model.Nodalstrain_e[:,0] + Model.Nodalstrain_e[:,1]) + (lam+2.*mu)*Model.Nodalstrain_e[:,3]
            xy  = Model.Nodalstress[:,2]

        elif Model.twoD == 'planestress':
            x   = Model.Nodalstress[:,0]
            y   = Model.Nodalstress[:,1]
            xy  = Model.Nodalstress[:,2]

        else:
            #print("Plane Stress is not ready")
            assert False, "Plane Stress or plane strain"

    elif Model.Dim == 3:
        x  = Model.Nodalstress[:,0]
        y  = Model.Nodalstress[:,1]
        z  = Model.Nodalstress[:,2]
        xy = Model.Nodalstress[:,3]
        yz = Model.Nodalstress[:,4]
        zx = Model.Nodalstress[:,5]

    else:
        assert False, "Dimension Error"

    Model.NodalsigmaVM = np.sqrt(0.5*((x-y)**2 + (y-z)**2 + (z-x)**2 + 6.0*(xy**2 + yz**2 + zx**2)))


#def CalculatePrincipalStressAtNode(Node):
    ## Calculate Von Mises stress at node
    ## the stress at nodes should be calculated first
    #Node.sigma1 = 0.5*(Node.stress[:,0] + Node.stress[:,1]) + np.sqrt(0.25*(Node.stress[:,0]-Node.stress[:,1])**2 + Node.stress[:,2]**2)
    #return



def CalculateStrainAtNode(Model):
    # Calculate strain at nodes
    # the strain at gauss points should be calculated first
    # this should be done in plasticity "ConstructStiffness"
    Dim    = Model.Dim
    LocalNElem  = np.zeros([Model.NNode])
    if Dim == 2:
        if Model.twoD == "planestress":
            strain_e = np.zeros((Model.NNode,3))
            strain_p = np.zeros((Model.NNode,3))
            Dim2 = 3
        elif Model.twoD == "planestrain":
            strain_e = np.zeros((Model.NNode,4))
            strain_p = np.zeros((Model.NNode,4))
            Dim2 = 4
    elif Dim ==3:
        strain_e = np.zeros((Model.NNode,6))
        strain_p = np.zeros((Model.NNode,6))
        Dim2 = 6
    else:
        assert False, "Dimension Error"

    for ii, element in enumerate(Model.Element):
        GP_SR  = GetGPSR(element)
        Ndof  = Model.GetNodeDof(element.Connectivity)

        ElementGPstrain_e = np.array(element.GPstrain_e)
        ElementGPstrain_p = np.array(element.GPstrain_p)
        strain_e[Ndof,:] += GP_SR @ ElementGPstrain_e
        strain_p[Ndof,:] += GP_SR @ ElementGPstrain_p
        LocalNElem[Ndof] += 1
    for ii in range(Dim2):
        strain_e[:,ii] = strain_e[:,ii]/LocalNElem
        strain_p[:,ii] = strain_p[:,ii]/LocalNElem

    Model.Nodalstrain_e = strain_e
    Model.Nodalstrain_p = strain_p


def CalculateStressAtNode(Model):
    # Calculate stress at nodes
    # the stress at gauss points should be calculated first in plasticity "ConstructStiffness"
    # the stress at nodes is recovered at node using stress recovery

    Dim    = Model.Dim

    if Dim == 2:
        if Model.twoD == "planestress":
            stress = np.zeros([Model.NNode,3])
            Dim2 = 3
        elif Model.twoD == "planestrain":
            stress = np.zeros([Model.NNode,3])
            Dim2 = 3
    elif Dim == 3:
        stress = np.zeros([Model.NNode,6])
        Dim2 = 6
    else:
        assert 0, "Dimension check"
    LocalNElem  = np.zeros([Model.NNode])


    for ii, element in enumerate(Model.Element):
        GP_SR = GetGPSR(element)
        Ndof  = Model.GetNodeDof(element.Connectivity)

        ElementGPstress = np.array(element.GPstress)
        stress[Ndof,:] += GP_SR @ ElementGPstress
        LocalNElem[Ndof] += 1

    for ii in range(Dim2):
        stress[:,ii] = stress[:,ii]/LocalNElem

    Model.Nodalstress = stress


def Assembleage(Model):
    # Stiffness matrix is assembled in Global scale
    Global_K = np.zeros([Model.NNode * Model.Dim, Model.NNode * Model.Dim])

    for EID, element in enumerate(Model.Element):
        ENdof = Model.ElementalNodeDof(element.Connectivity)
        Global_K[ix_(ENdof,ENdof)] += element.Stiffness
    return Global_K
    

def BC_SetUpAtStep(Model):

    Dim       = Model.Dim

    for ind, Node in enumerate(Model.Node):
        tmp = np.isnan(Node.BC_E)
        Node.du = np.zeros(Dim)
        for jj in range(Dim):
            if tmp[jj]:
                Node.F_ext[jj]  = Node.BC_N_init[jj]
                Node.F_ext[jj] += Node.BC_N[jj] * Model.step
            else:
                if Model.step != 0:                            ########################## check
                    Node.du[jj]    = Node.BC_E[jj]
                Node.u[jj]        += Node.du[jj]
    return


#def ElementalNodeDof(Fem, Node, Connectivity): ################################## Generalization check
    #tmp = NodeDof(Node, Connectivity)
    #if Fem.Dimension == 2:
        #tmp1 = np.array(tmp)*2
        #tmp2 = np.array(tmp)*2 + 1

        #c = np.empty((tmp1.size + tmp2.size), dtype =tmp1.dtype)
        #c[0::2] = tmp1
        #c[1::2] = tmp2
    #elif Fem.Dimension == 3:
        #tmp1 = np.array(tmp)*3
        #tmp2 = np.array(tmp)*3 + 1
        #tmp3 = np.array(tmp)*3 + 2

        #c = np.empty((tmp1.size + tmp2.size + tmp3.size), dtype =tmp1.dtype)
        #c[0::3] = tmp1
        #c[1::3] = tmp2
        #c[2::3] = tmp3

    #else:
        #assert False,"Dimension Generalization"
    #return c


def ShapeFunction(Fem, s,t,u=0.):

    if Fem.ElmentType == "q4":
        N_matrix = np.array([(s-1.)*(t-1.)/4.,\
                            (-s-1.)*(t-1.)/4.,\
                            (-s-1.)*(-t-1.)/4.,\
                            (s-1.)*(-t-1.)/4.])

        dN   = np.array([[t-1., 1.-t, 1.+t, -1.-t],
                        [s-1., -1.-s, 1.+s, 1.-s]])
        dN *= 0.25

    elif Fem.ElmentType == 't3':
    # T3
        N_matrix = np.array([s,\
                            t,\
                            1. - s - t])

        dN   = np.array([[1., 0., -1.],
                        [0., 1., -1.]])


    elif Fem.ElmentType == 't6':
    #T6
        N_matrix = np.zeros((6))
        u = 1.0 - s - t
        N_matrix[0] = 2.0 * s * s - s
        N_matrix[1] = 2.0 * t * t - t
        N_matrix[2] = 2.0 * u * u - u
        N_matrix[3] = 4.0 * s * t
        N_matrix[4] = 4.0 * t * u
        N_matrix[5] = 4.0 * u * s

        dN = np.zeros((2,6))
        dN[0,0] = 4.0 * s - 1.0
        dN[0,1] = 0.0
        dN[0,2] = 4.0 * (s + t) - 3.0
        dN[0,3] = 4.0 * t
        dN[0,4] = -4.0 * t
        dN[0,5] = 4.0 - 4.0 * t - 8.0 * s

        dN[1,0] = 0.0
        dN[1,1] = 4.0 * t - 1.0
        dN[1,2] = 4.0 * (s + t) - 3.0
        dN[1,3] = 4.0 * s
        dN[1,4] = 4.0 - 4.0 * s - 8.0 * t
        dN[1,5] = -4.0 * s

    elif Fem.ElmentType == 'hex8':
        N_matrix = np.zeros((8))
        N_matrix[0] = (1.-s)*(1.-t)*(1.+u)
        N_matrix[1] = (1.-s)*(1.-t)*(1.-u)
        N_matrix[2] = (1.-s)*(1.+t)*(1.-u)
        N_matrix[3] = (1.-s)*(1.+t)*(1.+u)
        N_matrix[4] = (1.+s)*(1.-t)*(1.+u)
        N_matrix[5] = (1.+s)*(1.-t)*(1.-u)
        N_matrix[6] = (1.+s)*(1.+t)*(1.-u)
        N_matrix[7] = (1.+s)*(1.+t)*(1.+u)

        N_matrix *= 0.125

        #dN = np.zeros((3,8))
        dN = np.array([[-(1.-t)*(1.+u), -(1.-t)*(1.-u), -(1.+t)*(1.-u), -(1.+t)*(1.+u),  (1.-t)*(1.+u),  (1.-t)*(1.-u),  (1.+t)*(1.-u),  (1.+t)*(1.+u)],
                       [-(1.-s)*(1.+u), -(1.-s)*(1.-u),  (1.-s)*(1.-u),  (1.-s)*(1.+u), -(1.+s)*(1.+u), -(1.+s)*(1.-u),  (1.+s)*(1.-u),  (1.+s)*(1.+u)],
                       [ (1.-s)*(1.-t), -(1.-s)*(1.-t), -(1.-s)*(1.+t),  (1.-s)*(1.+t),  (1.+s)*(1.-t), -(1.+s)*(1.-t), -(1.+s)*(1.+t),  (1.+s)*(1.+t)]])
        dN *= 0.125

    elif Fem.ElmentType == 'tet4':
        N_matrix = np.array([s,\
                            t,\
                            u,\
                            1.-s-t-u])

        dN   = np.array([[1., 0., 0.,-1.],
                         [0., 1., 0.,-1.],
                         [0., 0., 1.,-1.]])

    else:
        assert 0,"check the Fem.ElementType"

    return N_matrix, dN 

def Strain(Fem, NodeTarget,s,t,u=0.):
    # Cook, R.D., 2007. Concepts and applications of finite element analysis. John wiley & sons.
    G_Edof = Fem.G_Edof
    Dim    = Fem.Dim
    _, dN = ShapeFunction(Fem,s,t,u)
    Jacobian = np.matmul(dN, NodeTarget)

    if Dim == 2:
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
    else:
        B1= np.array([[1., 0., 0., 0., 0., 0., 0., 0., 0.],
                      [0., 0., 0., 0., 1., 0., 0., 0., 0.],
                      [0., 0., 0., 0., 0., 0., 0., 0., 1.],
                      [0., 1., 0., 1., 0., 0., 0., 0., 0.],
                      [0., 0., 0., 0., 0., 1., 0., 1., 0.],
                      [0., 0., 1., 0., 0., 0., 1., 0., 0.]])

        B2 = np.zeros([9,9])

        B2[0:3,0:3] = np.linalg.inv(Jacobian)
        B2[3:6,3:6] = np.linalg.inv(Jacobian)
        B2[6:9,6:9] = np.linalg.inv(Jacobian)

        B3 = np.zeros([9,G_Edof])

        for ind in range(int(G_Edof/Dim)):
            B3[0    :Dim*1,Dim*ind]   = dN[:,ind]
            B3[Dim*1:Dim*2,Dim*ind+1] = dN[:,ind]
            B3[Dim*2:Dim*3,Dim*ind+2] = dN[:,ind]

        
    B_matrix= B1 @ B2 @ B3
    Jacc = np.linalg.det(Jacobian)

    return Jacc, B_matrix
