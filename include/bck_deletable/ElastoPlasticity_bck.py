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
from util.tensor_operations import *
from util.coordinate_transforms import *
import importlib


class GPstate:
    eps_e: np.array
    eps_p: np.array
    deps: np.array
    lamda: float
    stress: np.array

def ElementSetUp(Fem, Node, Element):
    G_Edof = Fem.G_Edof
    Dim    = Fem.Dimension

    Element.B_matrix =[]
    Element.Jacc  =[]
    Element.Area  =[]
    Element.GPstrain0_e = []
    Element.GPstrain0_p = []
    Element.GPstress0 = []
    Element.GPlamda0 = []
    Element.GPstrain1_e = []
    Element.GPstrain1_p = []
    Element.GPstress1 = []
    Element.GPlamda1 = []
    for Edof in Element.Connectivity:
        #ENdof = ElementalNodeDof(Edof)
        NodeTarget = []
        for dof in Edof:
            NodeTarget.append(Node.Coord[Node.Id.index(dof)])

        Init_tensor = []
        Init_tensor1 = []

        B_matrix = []
        Jacc_Elem = []

        Area = 0.
        for gp in Fem.GP:
            Jacc, B = Strain(Fem, NodeTarget, gp[0], gp[1])
            Area = Area + Jacc * gp[-1]
            B_matrix.append(B)
            Jacc_Elem.append(Jacc)
            Init_tensor.append(np.zeros((3,1)))
            Init_tensor1.append(0.0)
        Element.B_matrix.append(B_matrix)
        Element.Jacc.append(Jacc_Elem)
        Element.Area.append(Area)

        Element.GPstrain0_e.append(Init_tensor)
        Element.GPstrain0_p.append(Init_tensor)
        Element.GPstress0.append(Init_tensor)
        Element.GPlamda0.append(Init_tensor1)
        Element.GPstrain1_e.append(Init_tensor)
        Element.GPstrain1_p.append(Init_tensor)
        Element.GPstress1.append(Init_tensor)
        Element.GPlamda1.append(Init_tensor1)
    return


def ConstructStiffness(Fem, Node, Element):
    G_Edof = Fem.G_Edof
    Dim    = Fem.Dimension
    IndexBCE  = list(np.where(np.invert(np.isnan(Node.BC_E)))[0])
    #print( Node.du[IndexBCE])

    ConstitutiveModel = ConstitutiveLaw(Fem)
    Element.Stiffness  =[]
    Element.GPstrain0_e = Element.GPstrain1_e
    Element.GPstrain0_p = Element.GPstrain1_p
    Element.GPlamda0    = Element.GPlamda1
    Element.GPstrain1_e = []
    Element.GPstrain1_p = []
    Element.GPstress    = []
    Element.GPlamda1    = []
    Node.F_int = np.zeros_like(Node.F_int)
    for ind1, Edof in enumerate(Element.Connectivity):
        ENdof = ElementalNodeDof(Node, Edof)
        E_K = np.zeros([G_Edof, G_Edof])
        ElemStrain_e = []
        ElemStrain_p = []
        ElemStress   = []
        Elemlamda    = []
        B    = Element.B_matrix[ind1]
        Jacc = Element.Jacc[ind1]
        GPeps_e = Element.GPstrain0_e[ind1]
        GPeps_p = Element.GPstrain0_p[ind1]
        GPlamda = Element.GPlamda0[ind1]
        GPF_int = np.zeros((len(ENdof)))

        #print(ind1)
        for ind, gp in enumerate(Fem.GP):
            B_GP    = B[ind]
            Jacc_GP = Jacc[ind]
            GPdata  =  GPstate()
            GPdata.eps_e  = GPeps_e[ind]
            GPdata.eps_p  = GPeps_p[ind]
            GPdata.lamda  = GPlamda[ind]
            GPdata.deps = B_GP @ Node.du[ENdof]
            ConstitutiveModel.ReturnMapping(GPdata)
            GP_K = B_GP.T @ GPdata.D @ B_GP * Fem.width * Jacc_GP * gp[-1]
            E_K += GP_K
            GPF_int += B_GP.T @ GPdata.stress * Fem.width * Jacc_GP * gp[-1]
            ElemStrain_e.append(GPdata.eps_e)
            ElemStrain_p.append(GPdata.eps_p)
            Elemlamda.append(GPdata.lamda)
            ElemStress.append(GPdata.stress)
        Element.GPstrain1_e.append(ElemStrain_e)
        Element.GPstrain1_p.append(ElemStrain_p)
        Element.GPstress.append(ElemStress)
        Element.GPlamda1.append(Elemlamda)
        Element.Stiffness.append(E_K)
        Node.F_int[ENdof] += GPF_int
    return


def Solve(Node, Global_K):
    # https://mae.ufl.edu/nkim/egm6352/Chap2.pdf
    # Kim, N.H., 2014. Introduction to nonlinear finite element analysis. Springer Science & Business Media.

    IndexBCN  = list(np.where(np.isnan(Node.BC_E))[0])
    IndexBCE  = list(np.where(np.invert(np.isnan(Node.BC_E)))[0])

    Node.u1 = np.copy(Node.u)

    Sliced_K12 = Global_K[ix_(IndexBCN, IndexBCE)]
    Sliced_K11 = Global_K[ix_(IndexBCN, IndexBCN)]

    BCE = Node.u[IndexBCE]
    BCN = Node.BC_N[IndexBCN]
    F_out = np.zeros_like(Node.BC_E)
    F_out[IndexBCN] = BCN
    F_total = F_out[IndexBCN] - Node.F_int[IndexBCN]
    #print("\n\n\n",F_total)
    #print(Node.F_int[IndexBCN])
    #print(F_out[IndexBCN])

    #print(IndexBCN)
    #print(IndexBCE)
    #print(Node.du,"\n")
    Node.du[IndexBCN] = np.linalg.solve(Sliced_K11, F_total)
    Node.u1[IndexBCN] = Node.u[IndexBCN] + Node.du[IndexBCN]
    Node.du[IndexBCE] = Node.u1[IndexBCE] - Node.u[IndexBCE]
    #print(Node.u,"\n")
    #print(Node.u1,"\n")
    #print(Node.du,"\n")
    #print((Global_K @ Node.du)[IndexBCN])
    #input("")
    Node.u  = np.copy(Node.u1)
    return np.linalg.norm(F_total)**2/(1.+np.linalg.norm(F_out[IndexBCN])**2)

class ConstitutiveLaw():
    def __init__(self, Fem):
        self.twoD = Fem.twoD
        self.K  = Fem.K
        self.mu = Fem.mu
        self.lam = Fem.lamda
        self.E  = Fem.E
        self.nu = Fem.nu
        self.ReturnMappingMode = Fem.ReturnMappingMode
        tmp = importlib.import_module(Fem.PlasticModel)
        self.PlasticityModel = tmp.MyPlasticity(Fem)
        #print(self.PlasticityModel.E)

    def Voigt2Tensor2DStrain(self, voigt):
        # https://www.comsol.com/blogs/what-is-the-difference-between-plane-stress-and-plane-strain
        eps_e = np.zeros((3,3))
        if self.twoD == 'planestrain':
            eps_e[0,0] = voigt[0]
            eps_e[0,1] = voigt[2]
            eps_e[1,0] = voigt[2]
            eps_e[1,1] = voigt[1]

        elif self.twoD == 'planestress':
            nu = self.nu
            eps_e[0,0] = voigt[0]
            eps_e[0,1] = voigt[2]
            eps_e[1,0] = voigt[2]
            eps_e[1,1] = voigt[1]
            eps_e[2,2] = -nu/(1.0-nu)*(eps_e[0,0] + eps_e[1,1])
        else:
            assert 0,"Check Fem.twoD"
        return eps_e

    def Tensor2Voigt2D(self, Tensor, flag='strain'):
        # https://www.comsol.com/blogs/what-is-the-difference-between-plane-stress-and-plane-strain
        voigt = np.zeros((3))
        if flag == 'strain':
            voigt[0] = Tensor[0,0]
            voigt[1] = Tensor[1,1]
            voigt[2] = Tensor[0,1]
        elif flag == 'stress':
            voigt[0] = Tensor[0,0]
            voigt[1] = Tensor[1,1]
            voigt[2] = Tensor[0,1]*0.5
        else:
            assert(0,"check the flag in Tensor2Voigt2D")
        return voigt

    def ReturnMapping(self, GPstate):
        if (self.ReturnMappingMode == 'eigenspace'):
            self.ReturnMappingEigenSpace(GPstate)
        elif(self.ReturnMappingMode == 'tensorspace' ):
            self.ReturnMappingTensorSpace(GPstate)
            #print("not prepared yet")
        elif(self.ReturnMappingMode == 'pqspace' ):
            print("not prepared yet")
            exit(1)
        #print(GPstate.sigma)
        #print(GPstate.D)
        #exit(1)
        return

    def ReturnMappingTensorSpace(self, GPstate):
        K  = self.K
        lam= self.lam
        mu = self.mu

        I   = np.eye(3)
        II  = identity_4(3)
        IxI = tensor_oMult(I,I)

        D = K*IxI+2.0*mu*(II-1.0/3.0*IxI)
        #D = lam*IxI + 2.*mu*II

        eps_e = self.Voigt2Tensor2DStrain( GPstate.eps_e )
        eps_p = self.Voigt2Tensor2DStrain( GPstate.eps_p )
        deps  = self.Voigt2Tensor2DStrain( GPstate.deps  )
        eps   = eps_e + eps_p   # strain
        sigma = np.zeros((3,3))
        lamda = GPstate.lamda

        # [1] Compute trial strain
        eps_e_tr = eps_e + deps

        eps_e_11_tr = eps_e_tr[0,0]
        eps_e_22_tr = eps_e_tr[1,1]
        eps_e_33_tr = eps_e_tr[2,2]
        eps_e_12_tr = eps_e_tr[0,1]
        eps_e_23_tr = eps_e_tr[1,2]
        eps_e_31_tr = eps_e_tr[2,0]

        # [2] Compute trial stress
        sigma_tr = sigma + np.tensordot(D, deps)

        sig11_tr = sigma_tr[0,0]
        sig22_tr = sigma_tr[1,1]
        sig33_tr = sigma_tr[2,2]
        sig12_tr = sigma_tr[0,1]
        sig23_tr = sigma_tr[1,2]
        sig31_tr = sigma_tr[2,0]

        YF = self.PlasticityModel.f(sig11_tr, sig22_tr, sig33_tr, sig12_tr, sig23_tr, sig31_tr, lamda) # yield surface
        #print(self.E)
        #print(self.nu)
        #print(self.lam)
        #print(self.mu)
        #print(eps_e_tr)
        #print(sigma_tr)
        #print(YF)
        #exit(1)
        #print(YF)
        #if YF <= 0.:
            #print(">> Elastic!")

            # Update stress & strain
        sigma = sigma_tr

        eps_e = eps_e_tr
        eps   = eps_e + eps_p

        ddsdde = FthTensor2Voigt(D)
        idx = [0,1,3]
        D = ddsdde[ix_(idx, idx)]
        #else:
            ##print(">> Plastic!")

            ## Initialize variables
            #eps_e_11 = eps_e[0,0]
            #eps_e_22 = eps_e[1,1]
            #eps_e_33 = eps_e[2,2]
            #eps_e_12 = eps_e[0,1]
            #eps_e_23 = eps_e[1,2]
            #eps_e_31 = eps_e[2,0]
            #dlamda = 0

            #x = np.zeros(7) # target variables
            #x[0] = eps_e_11
            #x[1] = eps_e_22
            #x[2] = eps_e_33
            #x[3] = eps_e_12
            #x[4] = eps_e_23
            #x[5] = eps_e_31
            #x[6] = dlamda

            ## Newton-Raphson iteration (return mapping)
            #for ii in range(20):
                ## Initialize residual and jacobian
                #res = np.zeros(7)
                #jac = np.zeros((7,7))

                ## Current strain
                #eps_e_11_current = x[0]
                #eps_e_22_current = x[1]
                #eps_e_33_current = x[2]
                #eps_e_12_current = x[3]
                #eps_e_23_current = x[4]
                #eps_e_31_current = x[5]


                ## Current stress
                #sig11_current = (lam+2*mu)*eps_e_11_current + lam*eps_e_22_current + lam*eps_e_33_current
                #sig22_current = lam*eps_e_11_current + (lam+2*mu)*eps_e_22_current + lam*eps_e_33_current
                #sig33_current = lam*eps_e_11_current + lam*eps_e_22_current + (lam+2*mu)*eps_e_33_current
                #sig12_current = 2*mu*eps_e_12_current
                #sig23_current = 2*mu*eps_e_23_current
                #sig31_current = 2*mu*eps_e_31_current

                ## Current lamda
                #lamda_current    = lamda + x[6]

                ## Update derivatives
                ## >> First order derivatives
                #dfdsig11, dfdsig22, dfdsig33, dfdsig12, dfdsig23, dfdsig31 = self.PlasticityModel.df(sig11_current, sig22_current, sig33_current, sig12_current, sig23_current, sig31_current)

                ## >> Second order derivatives
                #d2fdsig11sig11, d2fdsig11sig22, d2fdsig11sig33, d2fdsig11sig12, d2fdsig11sig23, d2fdsig11sig31,  \
                                #d2fdsig22sig22, d2fdsig22sig33, d2fdsig22sig12, d2fdsig22sig23, d2fdsig22sig31,  \
                                                #d2fdsig33sig33, d2fdsig33sig12, d2fdsig33sig23, d2fdsig33sig31,  \
                                                                #d2fdsig12sig12, d2fdsig12sig23, d2fdsig12sig31,  \
                                                                                #d2fdsig23sig23, d2fdsig23sig31,  \
                                                                                                #d2fdsig31sig31 = \
                #self.PlasticityModel.df2(sig11_current, sig22_current, sig33_current, sig12_current, sig23_current, sig31_current)

                ## Update residual
                #res[0] = x[0] - eps_e_11_tr + x[6]*dfdsig11
                #res[1] = x[1] - eps_e_22_tr + x[6]*dfdsig22
                #res[2] = x[2] - eps_e_33_tr + x[6]*dfdsig33
                #res[3] = x[3] - eps_e_12_tr + x[6]*dfdsig12
                #res[4] = x[4] - eps_e_23_tr + x[6]*dfdsig23
                #res[5] = x[5] - eps_e_31_tr + x[6]*dfdsig31
                #res[6] = self.PlasticityModel.f(sig11_current, sig22_current, sig33_current, sig12_current, sig23_current, sig31_current, lamda_current)

                ## Update jacobian ***
                #jac[0,0] = 1. + x[6]*(d2fdsig11sig11*dsig11deps11 + d2fdsig11sig22*dsig22deps11 + d2fdsig11sig33*dsig33deps11)
                #jac[0,1] =      x[6]*(d2fdsig11sig11*dsig11deps22 + d2fdsig11sig22*dsig22deps22 + d2fdsig11sig33*dsig33deps22)
                #jac[0,2] =      x[6]*(d2fdsig11sig11*dsig11deps33 + d2fdsig11sig22*dsig22deps33 + d2fdsig11sig33*dsig33deps33)
                #jac[0,3] =      x[6]*d2fdsig11sig12*dsig12deps12
                #jac[0,4] =      x[6]*d2fdsig11sig23*dsig23deps23
                #jac[0,5] =      x[6]*d2fdsig11sig31*dsig31deps31
                #jac[0,6] = dfdsig11

                #jac[1,0] =      x[6]*(d2fdsig11sig22*dsig11deps11 + d2fdsig22sig22*dsig22deps11 + d2fdsig22sig33*dsig33deps11)
                #jac[1,1] = 1. + x[6]*(d2fdsig11sig22*dsig11deps22 + d2fdsig22sig22*dsig22deps22 + d2fdsig22sig33*dsig33deps22)
                #jac[1,2] =      x[6]*(d2fdsig11sig22*dsig11deps33 + d2fdsig22sig22*dsig22deps33 + d2fdsig22sig33*dsig33deps33)
                #jac[1,3] =      x[6]*d2fdsig22sig12*dsig12deps12
                #jac[1,4] =      x[6]*d2fdsig22sig23*dsig23deps23
                #jac[1,5] =      x[6]*d2fdsig22sig31*dsig31deps31
                #jac[1,6] = dfdsig22

                #jac[2,0] =      x[6]*(d2fdsig11sig33*dsig11deps11 + d2fdsig22sig33*dsig22deps11 + d2fdsig33sig33*dsig33deps11)
                #jac[2,1] =      x[6]*(d2fdsig11sig33*dsig11deps22 + d2fdsig22sig33*dsig22deps22 + d2fdsig33sig33*dsig33deps22)
                #jac[2,2] = 1. + x[6]*(d2fdsig11sig33*dsig11deps33 + d2fdsig22sig33*dsig22deps33 + d2fdsig33sig33*dsig33deps33)
                #jac[2,3] =      x[6]*d2fdsig33sig12*dsig12deps12
                #jac[2,4] =      x[6]*d2fdsig33sig23*dsig23deps23
                #jac[2,5] =      x[6]*d2fdsig33sig31*dsig31deps31
                #jac[2,6] = dfdsig33

                #jac[3,0] =      x[6]*(d2fdsig11sig12*dsig11deps11 + d2fdsig22sig12*dsig22deps11 + d2fdsig33sig12*dsig33deps11)
                #jac[3,1] =      x[6]*(d2fdsig11sig12*dsig11deps22 + d2fdsig22sig12*dsig22deps22 + d2fdsig33sig12*dsig33deps22)
                #jac[3,2] =      x[6]*(d2fdsig11sig12*dsig11deps33 + d2fdsig22sig12*dsig22deps33 + d2fdsig33sig12*dsig33deps33)
                #jac[3,3] = 1. + x[6]*d2fdsig12sig12*dsig12deps12
                #jac[3,4] =      x[6]*d2fdsig12sig23*dsig23deps23
                #jac[3,5] =      x[6]*d2fdsig12sig31*dsig31deps31
                #jac[3,6] = dfdsig12

                #jac[4,0] =      x[6]*(d2fdsig11sig23*dsig11deps11 + d2fdsig22sig23*dsig22deps11 + d2fdsig33sig23*dsig33deps11)
                #jac[4,1] =      x[6]*(d2fdsig11sig23*dsig11deps22 + d2fdsig22sig23*dsig22deps22 + d2fdsig33sig23*dsig33deps22)
                #jac[4,2] =      x[6]*(d2fdsig11sig23*dsig11deps33 + d2fdsig22sig23*dsig22deps33 + d2fdsig33sig23*dsig33deps33)
                #jac[4,3] =      x[6]*d2fdsig12sig23*dsig12deps12
                #jac[4,4] = 1. + x[6]*d2fdsig23sig23*dsig23deps23
                #jac[4,5] =      x[6]*d2fdsig23sig31*dsig31deps31
                #jac[4,6] = dfdsig23

                #jac[5,0] =      x[6]*(d2fdsig11sig31*dsig11deps11 + d2fdsig22sig31*dsig22deps11 + d2fdsig33sig31*dsig33deps11)
                #jac[5,1] =      x[6]*(d2fdsig11sig31*dsig11deps22 + d2fdsig22sig31*dsig22deps22 + d2fdsig33sig31*dsig33deps33)
                #jac[5,2] =      x[6]*(d2fdsig11sig31*dsig11deps33 + d2fdsig22sig31*dsig22deps33 + d2fdsig33sig31*dsig33deps33)
                #jac[5,3] =      x[6]*d2fdsig12sig31*dsig12deps12
                #jac[5,4] =      x[6]*d2fdsig23sig31*dsig23deps23
                #jac[5,5] = 1. + x[6]*d2fdsig31sig31*dsig31deps31
                #jac[5,6] = dfdsig31

                #jac[6,0] = dfdsig11*dsig11deps11 + dfdsig22*dsig22deps11 + dfdsig33*dsig33deps11
                #jac[6,1] = dfdsig11*dsig11deps22 + dfdsig22*dsig22deps22 + dfdsig33*dsig33deps22
                #jac[6,2] = dfdsig11*dsig11deps33 + dfdsig22*dsig22deps33 + dfdsig33*dsig33deps33
                #jac[6,3] = dfdsig12*dsig12deps12
                #jac[6,4] = dfdsig23*dsig23deps23
                #jac[6,5] = dfdsig31*dsig31deps31
                #jac[6,6] = 0.

                ## Solve system of equations
                #dx = np.linalg.solve(jac, -res) # increment of target variables

                ## Update x
                #x = x + dx

                ## Compute error
                #err = np.linalg.norm(dx)

                #print(" Newton iter.",ii, ": err =", err)

                #if err < 1e-7:
                    #input("stop==")
                    #break

            ## Update strain
            #eps   += deps
            #eps_e = np.array([[x[0], x[3], x[5]],
                              #[x[3], x[1], x[4]],
                              #[x[5], x[4], x[2]]])
            #eps_p = eps - eps_e
            #lamda = lamda + x[6]

            ## Update stress
            #sigma = np.tensordot(D, eps_e)

            ## D tangent operator has to be implemented ######################################

            #ddsdde = FthTensor2Voigt(D)
            #idx = [0,1,3]
            #D = ddsdde[ix_(idx, idx)]

        GPstate.eps_e  = self.Tensor2Voigt2D(eps_e)
        GPstate.eps_p  = self.Tensor2Voigt2D(eps_p)
        GPstate.lamda  = lamda
        GPstate.stress = self.Tensor2Voigt2D(sigma,'stress')
        GPstate.D      = D
        return





    def ReturnMappingEigenSpace(self, GPstate):

        K  = self.K
        mu = self.mu
        nu = self.nu
        lam= self.lam

        I   = np.eye(3)
        II  = identity_4(3)
        IxI = tensor_oMult(I,I)

        #D = K*IxI+2.0*mu*(II-1.0/3.0*IxI)
        D = lam*IxI + 2.*mu*II

        a = K + (4/3)*mu
        b = K - (2/3)*mu

        Ce_principal = np.array([[a, b, b],
                            [b, a, b],
                            [b, b, a]])
        dsig1depse1 = a
        dsig1depse2 = b
        dsig1depse3 = b
        dsig2depse1 = b
        dsig2depse2 = a
        dsig2depse3 = b
        dsig3depse1 = b
        dsig3depse2 = b
        dsig3depse3 = a

        eps_e = self.Voigt2Tensor2DStrain( GPstate.eps_e )
        eps_p = self.Voigt2Tensor2DStrain( GPstate.eps_p )
        deps  = self.Voigt2Tensor2DStrain( GPstate.deps  )
        eps   = eps_e + eps_p   # strain
        sigma = np.zeros((3,3))
        lamda = GPstate.lamda

        eps_e_tr = eps_e + deps

        eps_e_tr_principal_mag, eps_e_tr_principal_vec = np.linalg.eig(eps_e_tr)

        eps_e_tr1 = eps_e_tr_principal_mag[0]
        eps_e_tr2 = eps_e_tr_principal_mag[1]
        eps_e_tr3 = eps_e_tr_principal_mag[2]

        n1 = eps_e_tr_principal_vec[:,0]
        n2 = eps_e_tr_principal_vec[:,1]
        n3 = eps_e_tr_principal_vec[:,2]

        # [2] Compute trial stress
        sigma_tr_principal_mag = np.inner(Ce_principal, eps_e_tr_principal_mag)

        sigma_tr1 = sigma_tr_principal_mag[0]
        sigma_tr2 = sigma_tr_principal_mag[1]
        sigma_tr3 = sigma_tr_principal_mag[2]

        sigma_tr = sigma_tr1*np.tensordot(n1,n1,axes=0) + sigma_tr2*np.tensordot(n2,n2,axes=0) + sigma_tr3*np.tensordot(n3,n3,axes=0)

        YF = self.PlasticityModel.f(sigma_tr1, sigma_tr2, sigma_tr3, lamda) # yield surface

        #print(YF)
        if YF <= 0.:
            #print(">> Elastic!")

            # Update stress & strain
            sigma = sigma_tr

            eps_e = eps_e_tr
            eps   = eps_e + eps_p

            ddsdde = FthTensor2Voigt(D)
            idx = [0,1,3]
            D = ddsdde[ix_(idx, idx)]
            #print(D)
        else:
            #print(">> Plastic!")
            # Initialize variables
            eps_e_principal_mag, eps_e_principal_vec = np.linalg.eig(eps_e)

            eps_e1 = eps_e_principal_mag[0]
            eps_e2 = eps_e_principal_mag[1]
            eps_e3 = eps_e_principal_mag[2]
            dlamda  = 0.

            x = np.zeros(4) # target variables
            x[0] = eps_e1
            x[1] = eps_e2
            x[2] = eps_e3
            x[3] = dlamda

            # Newton-Raphson iteration (return mapping)
            for ii in range(20):
                # Initialize residual and jacobian
                res = np.zeros(4)
                jac = np.zeros((4,4))

                # Current strain
                eps_e1_current = x[0]
                eps_e2_current = x[1]
                eps_e3_current = x[2]

                # Current stress
                sigma1_current = a*eps_e1_current + b*eps_e2_current + b*eps_e3_current
                sigma2_current = b*eps_e1_current + a*eps_e2_current + b*eps_e3_current
                sigma3_current = b*eps_e1_current + b*eps_e2_current + a*eps_e3_current

                sigma1_current = sigma1_current
                sigma2_current = sigma2_current
                sigma3_current = sigma3_current

                # Current lamda
                lamda_current = lamda + x[3]

                # Update derivatives
                # >> First order derivatives
                dfdsig1, dfdsig2, dfdsig3 = self.PlasticityModel.df(sigma1_current, sigma2_current, sigma3_current)

                # >> Second order derivatives
                d2fdsig1dsig1, d2fdsig2dsig2, d2fdsig3dsig3, d2fdsig1dsig2, d2fdsig2dsig3, d2fdsig3dsig1 \
                    = self.PlasticityModel.df2(sigma1_current, sigma2_current, sigma3_current)

                # Update residual
                res[0] = x[0] - eps_e_tr1 + x[3]*dfdsig1
                res[1] = x[1] - eps_e_tr2 + x[3]*dfdsig2
                res[2] = x[2] - eps_e_tr3 + x[3]*dfdsig3
                res[3] = self.PlasticityModel.f(sigma1_current, sigma2_current, sigma3_current, lamda_current)

                # Update Jacobian ***
                jac[0,0] = 1 + x[3]*(d2fdsig1dsig1*dsig1depse1 + d2fdsig1dsig2*dsig2depse1 + d2fdsig3dsig1*dsig3depse1)
                jac[0,1] =     x[3]*(d2fdsig1dsig1*dsig1depse2 + d2fdsig1dsig2*dsig2depse2 + d2fdsig3dsig1*dsig3depse2)
                jac[0,2] =     x[3]*(d2fdsig1dsig1*dsig1depse3 + d2fdsig1dsig2*dsig2depse3 + d2fdsig3dsig1*dsig3depse3)
                jac[0,3] = dfdsig1

                jac[1,0] =     x[3]*(d2fdsig1dsig2*dsig1depse1 + d2fdsig2dsig2*dsig2depse1 + d2fdsig2dsig3*dsig3depse1)
                jac[1,1] = 1 + x[3]*(d2fdsig1dsig2*dsig1depse2 + d2fdsig2dsig2*dsig2depse2 + d2fdsig2dsig3*dsig3depse2)
                jac[1,2] =     x[3]*(d2fdsig1dsig2*dsig1depse3 + d2fdsig2dsig2*dsig2depse3 + d2fdsig2dsig3*dsig3depse3)
                jac[1,3] = dfdsig2

                jac[2,0] =     x[3]*(d2fdsig3dsig1*dsig1depse1 + d2fdsig2dsig3*dsig2depse1 + d2fdsig3dsig3*dsig3depse1)
                jac[2,1] =     x[3]*(d2fdsig3dsig1*dsig1depse2 + d2fdsig2dsig3*dsig2depse2 + d2fdsig3dsig3*dsig3depse2)
                jac[2,2] = 1 + x[3]*(d2fdsig3dsig1*dsig1depse3 + d2fdsig2dsig3*dsig2depse3 + d2fdsig3dsig3*dsig3depse3)
                jac[2,3] = dfdsig3

                jac[3,0] = dfdsig1*dsig1depse1 + dfdsig2*dsig2depse1 + dfdsig3*dsig3depse1
                jac[3,1] = dfdsig1*dsig1depse2 + dfdsig2*dsig2depse2 + dfdsig3*dsig3depse2
                jac[3,2] = dfdsig1*dsig1depse3 + dfdsig2*dsig2depse3 + dfdsig3*dsig3depse3
                jac[3,3] = 0

                # Solve system of equations
                dx = np.linalg.solve(jac, -res) # increment of target variables

                # Update x
                x = x + dx

                # Compute error
                err = np.linalg.norm(dx)

                #print(" Newton iter.",ii, ": err =", err)

                if err < 1e-7:
                    break

            # Update strain
            eps   = eps + deps
            eps_e = x[0]*np.tensordot(n1,n1,axes=0) + x[1]*np.tensordot(n2,n2,axes=0) + x[2]*np.tensordot(n3,n3,axes=0)
            eps_p = eps - eps_e
            lamda = lamda + x[3]

            # Update stress
            sigma1 = a*x[0] + b*x[1] + b*x[2]
            sigma2 = b*x[0] + a*x[1] + b*x[2]
            sigma3 = b*x[0] + b*x[1] + a*x[2]
            sigma  = sigma1*np.tensordot(n1,n1,axes=0) + sigma2*np.tensordot(n2,n2,axes=0) + sigma3*np.tensordot(n3,n3,axes=0)

            aAB = np.zeros((3,3))
            mab = np.zeros((3,3,3,3))
            dsigdepse = np.zeros((3,3))
            dsigdepse[0,0] = dsig1depse1
            dsigdepse[0,1] = dsig1depse2
            dsigdepse[0,2] = dsig1depse3
            dsigdepse[1,0] = dsig2depse1
            dsigdepse[1,1] = dsig2depse2
            dsigdepse[1,2] = dsig2depse3
            dsigdepse[2,0] = dsig3depse1
            dsigdepse[2,1] = dsig3depse2
            dsigdepse[2,2] = dsig3depse3

            #aAB = np.linalg.solve(jac[0:3,0:3],dsigdepse)
            aAB = dsigdepse
            bAB = np.zeros((3,3))
            eps_e_tr = np.zeros((3))
            sig_e_cur = np.zeros((3))
            eps_e_tr[0] = eps_e_tr1
            eps_e_tr[1] = eps_e_tr2
            eps_e_tr[2] = eps_e_tr3
            sig_e_cur[0] = sigma[0,0]
            sig_e_cur[1] = sigma[1,1]
            sig_e_cur[2] = sigma[2,2]

            for ii in range(3):
                for jj in range(3):
                    if(np.abs(eps_e_tr[jj] - eps_e_tr[ii])) < 1e-10:
                        bAB[ii,jj] = 0.0
                    else:
                        bAB[ii,jj] = (sig_e_cur[jj] - sig_e_cur[ii])/(eps_e_tr[jj] - eps_e_tr[ii])

            for ii in range(3):
                for jj in range(3):
                    na = eps_e_principal_vec[:,ii].reshape(-1,1)
                    nb = eps_e_principal_vec[:,jj].reshape(-1,1)
                    mab[ii,jj,:,:] = na @ nb.T

            Calgo = np.zeros((3,3,3,3,))

            for ii in range(3):
                for jj in range(3):
                    for mm in range(3):
                        for nn in range(3):
                            Calgo[ii,jj,mm,nn] += aAB[0,0]*mab[0,0,ii,jj]*mab[0,0,mm,nn]\
                                                + aAB[0,1]*mab[0,0,ii,jj]*mab[1,1,mm,nn]\
                                                + aAB[0,2]*mab[0,0,ii,jj]*mab[2,2,mm,nn]\
                                                + aAB[1,0]*mab[1,1,ii,jj]*mab[0,0,mm,nn]\
                                                + aAB[1,1]*mab[1,1,ii,jj]*mab[1,1,mm,nn]\
                                                + aAB[1,2]*mab[1,1,ii,jj]*mab[2,2,mm,nn]\
                                                + aAB[2,1]*mab[2,2,ii,jj]*mab[0,0,mm,nn]\
                                                + aAB[2,2]*mab[2,2,ii,jj]*mab[1,1,mm,nn]\
                                                + aAB[2,2]*mab[2,2,ii,jj]*mab[2,2,mm,nn]\
                                            +0.5* bAB[0,1]*(mab[0,1,ii,jj]*mab[0,1,mm,nn]+mab[0,1,ii,jj]*mab[1,0,mm,nn])\
                                            +0.5* bAB[0,2]*(mab[0,2,ii,jj]*mab[0,2,mm,nn]+mab[0,1,ii,jj]*mab[2,0,mm,nn])\
                                            +0.5* bAB[1,0]*(mab[1,0,ii,jj]*mab[1,0,mm,nn]+mab[1,0,ii,jj]*mab[0,1,mm,nn])\
                                            +0.5* bAB[1,2]*(mab[1,2,ii,jj]*mab[1,2,mm,nn]+mab[1,2,ii,jj]*mab[2,1,mm,nn])\
                                            +0.5* bAB[2,0]*(mab[2,0,ii,jj]*mab[2,0,mm,nn]+mab[2,0,ii,jj]*mab[0,2,mm,nn])\
                                            +0.5* bAB[2,1]*(mab[2,1,ii,jj]*mab[2,1,mm,nn]+mab[2,1,ii,jj]*mab[1,2,mm,nn])

            ddsdde = FthTensor2Voigt(Calgo)
            idx = [0,1,3]
            D = ddsdde[ix_(idx, idx)]

        GPstate.eps_e  = self.Tensor2Voigt2D(eps_e)
        GPstate.eps_p  = self.Tensor2Voigt2D(eps_p)
        GPstate.lamda  = lamda
        GPstate.stress = self.Tensor2Voigt2D(sigma,'stress')
        GPstate.D      = D
        return

