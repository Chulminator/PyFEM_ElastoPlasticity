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
    deps_e: np.array
    deps_p: np.array
    dstress: np.array
    lamda: float
    stress: np.array

def ElementSetUp(Fem, Node, Element):
    G_Edof = Fem.G_Edof
    Dim    = Fem.Dimension

    Element.B_matrix =[]
    Element.Jacc  =[]
    Element.Area  =[]
    Element.GPstrain_e = []
    Element.GPstrain_p = []
    Element.GPlamda0 = []
    Element.GPlamda1 = []
    Element.GPlamda1 = []
    Element.GPstress = []
    #i = 0
    for Edof in Element.Connectivity:
        #i +=1
        #print(i)
        #ENdof = ElementalNodeDof(Edof)
        NodeTarget = []
        for dof in Edof:
            NodeTarget.append(Node.Coord[Node.Id.index(dof)])
        #print(*NodeTarget, sep='\n')
        Init_tensor = []
        Init_tensor1 = []
        Init_tensor2 = []
        B_matrix = []
        Jacc_Elem = []
        Area = 0.

        for gp in Fem.GP:
            if Fem.Dimension == 2:
                Jacc, B = Strain(Fem, NodeTarget, gp[0], gp[1])
            elif Fem.Dimension ==3:
                Jacc, B = Strain(Fem, NodeTarget, gp[0], gp[1], gp[2])
            else:
                assert(0,"Check Fem.Dimension")

            Area = Area + Jacc * gp[-1]
            B_matrix.append(B)
            Jacc_Elem.append(Jacc)
            if Fem.Dimension ==2:
                tmp_tensor = np.zeros((3))
                tmp_tensor1 = np.zeros((3))
                tmp_tensor1[0] = Fem.HSP
                tmp_tensor1[1] = Fem.HSP
            elif Fem.Dimension ==3:
                tmp_tensor = np.zeros((6))
                tmp_tensor1 = np.zeros((6))
                tmp_tensor1[0] = Fem.HSP
                tmp_tensor1[1] = Fem.HSP
                tmp_tensor1[2] = Fem.HSP
            else:
                assert(0,"Dimension error")
            Init_tensor.append(tmp_tensor)
            Init_tensor1.append(0.0)
            Init_tensor2.append(tmp_tensor1)
        #print(*B_matrix,sep="\n")
        #print(*Init_tensor,sep="\n")
        #exit(1)
        Element.B_matrix.append(B_matrix)
        Element.Jacc.append(Jacc_Elem)
        Element.Area.append(Area)

        Element.GPstrain_e.append(Init_tensor)
        Element.GPstrain_p.append(Init_tensor)
        Element.GPlamda0.append(Init_tensor1)
        Element.GPlamda1.append(Init_tensor1)
        Element.GPstress.append(Init_tensor2)

    #print(Element.GPstress)
    #exit(1)
    return


def ConstructStiffness(Fem, Node, Element, ConstitutiveModel):
    G_Edof = Fem.G_Edof
    Dim    = Fem.Dimension
    #IndexBCE  = list(np.where(np.invert(np.isnan(Node.BC_E)))[0])

    #ConstitutiveModel = ConstitutiveLaw(Fem)
    Element.Stiffness  =[]
    Element.GPlamda0    = Element.GPlamda1
    Element.GPstrain1_e = []
    Element.GPstrain1_p = []
    Element.GPstress1   = []
    Element.GPlamda1    = []
    Node.F_int = np.zeros_like(Node.F_int)
    for ind1, Edof in enumerate(Element.Connectivity):
        ENdof = ElementalNodeDof(Fem, Node, Edof)
        E_K = np.zeros([G_Edof, G_Edof])
        Elemlamda    = []
        B    = Element.B_matrix[ind1]
        Jacc = Element.Jacc[ind1]
        GPstress= Element.GPstress[ind1]
        GPeps_e = Element.GPstrain_e[ind1]
        GPeps_p = Element.GPstrain_p[ind1]
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
            GPdata.deps   = B_GP @ Node.du[ENdof]
            GPdata.stress = GPstress[ind]
            
            ConstitutiveModel.ReturnMapping(GPdata)

            #GPstress[ind] += GPdata.dstress
            GPeps_e[ind]  += GPdata.deps_e
            GPeps_p[ind]  += GPdata.deps_p
            GPstress[ind] =  GPdata.stress
            #GPeps_e[ind]  = GPdata.eps_e
            #GPeps_p[ind]  = GPdata.eps_p

            if Fem.Dimension == 2:
                GP_K = B_GP.T @ GPdata.D @ B_GP * Fem.width * Jacc_GP * gp[-1]
                GPF_int += B_GP.T @ GPstress[ind] * Fem.width * Jacc_GP * gp[-1]
            elif Fem.Dimension == 3:
                GP_K = B_GP.T @ GPdata.D @ B_GP * Jacc_GP * gp[-1]
                GPF_int += B_GP.T @ GPstress[ind] * Jacc_GP * gp[-1]
            else:
                assert 0, "Dimension error"

            E_K += GP_K
            Elemlamda.append(GPdata.lamda)
        Element.GPlamda1.append(Elemlamda)
        Element.Stiffness.append(E_K)
        Node.F_int[ENdof] += GPF_int
    return


def Solve(Node, Global_K):
    # https://mae.ufl.edu/nkim/egm6352/Chap2.pdf
    # Kim, N.H., 2014. Introduction to nonlinear finite element analysis. Springer Science & Business Media.
    # [[K11 K12],  [[d1],     [[Fout1],
    #  [K21 K22]]   [d2]]   =  [Fout2]]
    #
    # K12 * d2 + K11 * d1 = F1
    # d1 = K11 \ ( Fout1 - K12 * d2 )

    IndexBCN  = list(np.where(np.isnan(Node.BC_E))[0])
    IndexBCE  = list(np.where(np.invert(np.isnan(Node.BC_E)))[0])

    Node.u1 = np.copy(Node.u)

    Sliced_K12 = Global_K[ix_(IndexBCN, IndexBCE)]
    Sliced_K11 = Global_K[ix_(IndexBCN, IndexBCN)]

    F_total = Node.F_ext[IndexBCN] - Node.F_int[IndexBCN]

    Node.du[IndexBCN] = np.linalg.solve(Sliced_K11, F_total)
    Node.u1[IndexBCN] = Node.u[IndexBCN] + Node.du[IndexBCN]
    Node.du[IndexBCE] = Node.u1[IndexBCE] - Node.u[IndexBCE]

    Node.u  = np.copy(Node.u1)
    return np.linalg.norm(F_total)**2/(1.+np.linalg.norm(Node.F_ext[IndexBCN])**2)

class ConstitutiveLaw():
    def __init__(self, Fem):
        self.twoD = Fem.twoD
        self.Dim  = Fem.Dimension
        self.K  = Fem.K
        self.mu = Fem.mu
        self.lam = Fem.lamda
        self.E  = Fem.E
        self.nu = Fem.nu
        self.ReturnMappingMode = Fem.ReturnMappingMode
        tmp = importlib.import_module(Fem.PlasticModel)
        self.PlasticityModel = tmp.MyPlasticity(Fem)
        self.MatProp = Fem.MatProp
        self.HSP = Fem.HSP

    def Voigt2Tensor(self, voigt, flag='strain'):
        # https://www.comsol.com/blogs/what-is-the-difference-between-plane-stress-and-plane-strain
        Tensor = np.zeros((3,3))
        if flag.lower() == 'strain':
            if self.twoD == 'planestrain':
                Tensor[0,0] = voigt[0]
                Tensor[0,1] = voigt[2]
                Tensor[1,0] = voigt[2]
                Tensor[1,1] = voigt[1]

            elif self.twoD == 'planestress':
                nu = self.nu
                Tensor[0,0] = voigt[0]
                Tensor[0,1] = voigt[2]
                Tensor[1,0] = voigt[2]
                Tensor[1,1] = voigt[1]
                Tensor[2,2] = -nu/(1.0-nu)*(Tensor[0,0] + Tensor[1,1])
            else:
                if self.Dim == 3:
                    Tensor[0,0] = voigt[0]
                    Tensor[1,1] = voigt[1]
                    Tensor[2,2] = voigt[2]

                    Tensor[0,1] = voigt[3] # xy
                    Tensor[1,0] = voigt[3]

                    Tensor[1,2] = voigt[4] # yz
                    Tensor[2,1] = voigt[4]

                    Tensor[0,2] = voigt[5] # zx
                    Tensor[2,0] = voigt[5]
                else:
                    assert 0,"Check Fem.twoD or Fem.Dimension"
        elif flag.lower() == 'stress':
            if self.twoD == 'planestrain':
                nu = self.nu
                Tensor[0,0] = voigt[0]
                Tensor[0,1] = voigt[2]
                Tensor[1,0] = voigt[2]
                Tensor[1,1] = voigt[1]
                Tensor[2,2] = nu*(Tensor[0,0] + Tensor[1,1])

            elif self.twoD == 'planestress':
                nu = self.nu
                Tensor[0,0] = voigt[0]
                Tensor[0,1] = voigt[2]
                Tensor[1,0] = voigt[2]
                Tensor[1,1] = voigt[1]
                
            else:
                if self.Dim == 3:
                    Tensor[0,0] = voigt[0]
                    Tensor[1,1] = voigt[1]
                    Tensor[2,2] = voigt[2]

                    Tensor[0,1] = voigt[3] # xy
                    Tensor[1,0] = voigt[3]

                    Tensor[1,2] = voigt[4] # yz
                    Tensor[2,1] = voigt[4]

                    Tensor[0,2] = voigt[5] # zx
                    Tensor[2,0] = voigt[5]
                else:
                    assert 0,"Check Fem.twoD or Fem.Dimension"

        return Tensor

    def Tensor2Voigt(self, Tensor, flag='strain'):        
        if self.Dim ==2:
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
                assert(0,"check the flag in Tensor2Voigt")
        elif self.Dim ==3:
            voigt = np.zeros((6))
            if flag == 'strain':
                voigt[0] = Tensor[0,0]
                voigt[1] = Tensor[1,1]
                voigt[2] = Tensor[2,2]
                voigt[3] = Tensor[0,1]
                voigt[4] = Tensor[1,2]
                voigt[5] = Tensor[2,0]
            elif flag == 'stress':
                voigt[0] = Tensor[0,0]
                voigt[1] = Tensor[1,1]
                voigt[2] = Tensor[2,2]
                voigt[3] = Tensor[0,1]*0.5
                voigt[4] = Tensor[1,2]*0.5
                voigt[5] = Tensor[2,0]*0.5
            else:
                assert(0,"check the flag in Tensor2Voigt")
        else:
            assert 0,"Check the dimension"
        return voigt

    def InitialaPressure(self):
        if self.Dim == 2:
            tmp =np.array([self.HSP, self.HSP, 0.])
        elif self.Dim == 3:
            tmp =np.array([self.HSP, self.HSP, self.HSP, 0., 0., 0.])
        else:
            assert False, "Dimesion Check!"

        return self.Voigt2Tensor(tmp, flag='stress')

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
        assert Falser, "Not implemented"
        K  = self.K
        lam= self.lam
        mu = self.mu

        dsig11deps11 = lam + 2.*mu
        dsig11deps22 = lam
        dsig11deps33 = lam

        dsig22deps11 = lam
        dsig22deps22 = lam + 2.*mu
        dsig22deps33 = lam

        dsig33deps11 = lam
        dsig33deps22 = lam
        dsig33deps33 = lam + 2.*mu

        dsig12deps12 = 2.*mu
        dsig23deps23 = 2.*mu
        dsig31deps31 = 2.*mu

        I   = np.eye(3)
        II  = identity_4(3)
        IxI = tensor_oMult(I,I)

        #D = K*IxI+2.0*mu*(II-1.0/3.0*IxI)
        D = lam*IxI + 2.*mu*II

        eps_e = self.Voigt2Tensor( GPstate.eps_e )
        eps_p = self.Voigt2Tensor( GPstate.eps_p )
        deps  = self.Voigt2Tensor( GPstate.deps  )
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
        sigma_tr = np.tensordot(D, eps_e_tr)

        sig11_tr = sigma_tr[0,0]
        sig22_tr = sigma_tr[1,1]
        sig33_tr = sigma_tr[2,2]
        sig12_tr = sigma_tr[0,1]
        sig23_tr = sigma_tr[1,2]
        sig31_tr = sigma_tr[2,0]

        YF = self.PlasticityModel.f(sig11_tr, sig22_tr, sig33_tr, sig12_tr, sig23_tr, sig31_tr, lamda) # yield surface

        if YF <= 0.:
            #print(">> Elastic!")
            # Update stress & strain
            sigma = sigma_tr

            eps_e = eps_e_tr
            eps   = eps_e + eps_p

            ddsdde = FthTensor2Voigt(D)
            if self.Dim ==2:
                idx = [0,1,3]
                D = ddsdde[ix_(idx, idx)]
            elif self.Dim ==3:
                D = ddsdde
                #print(D)
                #print(D.shape)
                #print("under construction here!!!")
                #exit(1)
        else:
            #print(">> Plastic!")
            # Initialize variables
            eps_e_11 = eps_e[0,0]
            eps_e_22 = eps_e[1,1]
            eps_e_33 = eps_e[2,2]
            eps_e_12 = eps_e[0,1]
            eps_e_23 = eps_e[1,2]
            eps_e_31 = eps_e[2,0]
            dlamda = 0.

            x = np.zeros(7) # target variables
            x[0] = eps_e_11
            x[1] = eps_e_22
            x[2] = eps_e_33
            x[3] = eps_e_12
            x[4] = eps_e_23
            x[5] = eps_e_31
            x[6] = dlamda

            # Newton-Raphson iteration (return mapping)
            for ii in range(20):
                # Initialize residual and jacobian
                res = np.zeros(7)
                jac = np.zeros((7,7))

                # Current strain
                eps_e_11_current = x[0]
                eps_e_22_current = x[1]
                eps_e_33_current = x[2]
                eps_e_12_current = x[3]
                eps_e_23_current = x[4]
                eps_e_31_current = x[5]


                # Current stress
                sig11_current = (lam+2.*mu)*eps_e_11_current + lam*eps_e_22_current + lam*eps_e_33_current
                sig22_current = lam*eps_e_11_current + (lam+2.*mu)*eps_e_22_current + lam*eps_e_33_current
                sig33_current = lam*eps_e_11_current + lam*eps_e_22_current + (lam+2.*mu)*eps_e_33_current
                sig12_current = 2.*mu*eps_e_12_current
                sig23_current = 2.*mu*eps_e_23_current
                sig31_current = 2.*mu*eps_e_31_current

                # Current lamda
                lamda_current    = lamda + x[6]
                if lamda_current <0:
                    print("lamda_current\t",lamda_current)
                    input("")
                    lamda_current = 0.
                    x[6] = -lamda

                # Update derivatives
                # >> First order derivatives
                dfdsig11, dfdsig22, dfdsig33, dfdsig12, dfdsig23, dfdsig31, dfdlamda = self.PlasticityModel.df(sig11_current, sig22_current, sig33_current, sig12_current, sig23_current, sig31_current, lamda_current)

                # >> Second order derivatives
                d2fdsig11sig11, d2fdsig11sig22, d2fdsig11sig33, d2fdsig11sig12, d2fdsig11sig23, d2fdsig11sig31,  \
                                d2fdsig22sig22, d2fdsig22sig33, d2fdsig22sig12, d2fdsig22sig23, d2fdsig22sig31,  \
                                                d2fdsig33sig33, d2fdsig33sig12, d2fdsig33sig23, d2fdsig33sig31,  \
                                                                d2fdsig12sig12, d2fdsig12sig23, d2fdsig12sig31,  \
                                                                                d2fdsig23sig23, d2fdsig23sig31,  \
                                                                                                d2fdsig31sig31 = \
                self.PlasticityModel.df2(sig11_current, sig22_current, sig33_current, sig12_current, sig23_current, sig31_current, lamda_current)

                # Update residual
                res[0] = x[0] - eps_e_11_tr + x[6]*dfdsig11
                res[1] = x[1] - eps_e_22_tr + x[6]*dfdsig22
                res[2] = x[2] - eps_e_33_tr + x[6]*dfdsig33
                res[3] = x[3] - eps_e_12_tr + x[6]*dfdsig12
                res[4] = x[4] - eps_e_23_tr + x[6]*dfdsig23
                res[5] = x[5] - eps_e_31_tr + x[6]*dfdsig31
                res[6] = self.PlasticityModel.f(sig11_current, sig22_current, sig33_current, sig12_current, sig23_current, sig31_current, lamda_current)

                # Update jacobian ***
                jac[0,0] = 1. + x[6]*(d2fdsig11sig11*dsig11deps11 + d2fdsig11sig22*dsig22deps11 + d2fdsig11sig33*dsig33deps11)
                jac[0,1] =      x[6]*(d2fdsig11sig11*dsig11deps22 + d2fdsig11sig22*dsig22deps22 + d2fdsig11sig33*dsig33deps22)
                jac[0,2] =      x[6]*(d2fdsig11sig11*dsig11deps33 + d2fdsig11sig22*dsig22deps33 + d2fdsig11sig33*dsig33deps33)
                jac[0,3] =      x[6]*d2fdsig11sig12*dsig12deps12
                jac[0,4] =      x[6]*d2fdsig11sig23*dsig23deps23
                jac[0,5] =      x[6]*d2fdsig11sig31*dsig31deps31
                jac[0,6] = dfdsig11

                jac[1,0] =      x[6]*(d2fdsig11sig22*dsig11deps11 + d2fdsig22sig22*dsig22deps11 + d2fdsig22sig33*dsig33deps11)
                jac[1,1] = 1. + x[6]*(d2fdsig11sig22*dsig11deps22 + d2fdsig22sig22*dsig22deps22 + d2fdsig22sig33*dsig33deps22)
                jac[1,2] =      x[6]*(d2fdsig11sig22*dsig11deps33 + d2fdsig22sig22*dsig22deps33 + d2fdsig22sig33*dsig33deps33)
                jac[1,3] =      x[6]*d2fdsig22sig12*dsig12deps12
                jac[1,4] =      x[6]*d2fdsig22sig23*dsig23deps23
                jac[1,5] =      x[6]*d2fdsig22sig31*dsig31deps31
                jac[1,6] = dfdsig22

                jac[2,0] =      x[6]*(d2fdsig11sig33*dsig11deps11 + d2fdsig22sig33*dsig22deps11 + d2fdsig33sig33*dsig33deps11)
                jac[2,1] =      x[6]*(d2fdsig11sig33*dsig11deps22 + d2fdsig22sig33*dsig22deps22 + d2fdsig33sig33*dsig33deps22)
                jac[2,2] = 1. + x[6]*(d2fdsig11sig33*dsig11deps33 + d2fdsig22sig33*dsig22deps33 + d2fdsig33sig33*dsig33deps33)
                jac[2,3] =      x[6]*d2fdsig33sig12*dsig12deps12
                jac[2,4] =      x[6]*d2fdsig33sig23*dsig23deps23
                jac[2,5] =      x[6]*d2fdsig33sig31*dsig31deps31
                jac[2,6] = dfdsig33

                jac[3,0] =      x[6]*(d2fdsig11sig12*dsig11deps11 + d2fdsig22sig12*dsig22deps11 + d2fdsig33sig12*dsig33deps11)
                jac[3,1] =      x[6]*(d2fdsig11sig12*dsig11deps22 + d2fdsig22sig12*dsig22deps22 + d2fdsig33sig12*dsig33deps22)
                jac[3,2] =      x[6]*(d2fdsig11sig12*dsig11deps33 + d2fdsig22sig12*dsig22deps33 + d2fdsig33sig12*dsig33deps33)
                jac[3,3] = 1. + x[6]*d2fdsig12sig12*dsig12deps12
                jac[3,4] =      x[6]*d2fdsig12sig23*dsig23deps23
                jac[3,5] =      x[6]*d2fdsig12sig31*dsig31deps31
                jac[3,6] = dfdsig12

                jac[4,0] =      x[6]*(d2fdsig11sig23*dsig11deps11 + d2fdsig22sig23*dsig22deps11 + d2fdsig33sig23*dsig33deps11)
                jac[4,1] =      x[6]*(d2fdsig11sig23*dsig11deps22 + d2fdsig22sig23*dsig22deps22 + d2fdsig33sig23*dsig33deps22)
                jac[4,2] =      x[6]*(d2fdsig11sig23*dsig11deps33 + d2fdsig22sig23*dsig22deps33 + d2fdsig33sig23*dsig33deps33)
                jac[4,3] =      x[6]*d2fdsig12sig23*dsig12deps12
                jac[4,4] = 1. + x[6]*d2fdsig23sig23*dsig23deps23
                jac[4,5] =      x[6]*d2fdsig23sig31*dsig31deps31
                jac[4,6] = dfdsig23

                jac[5,0] =      x[6]*(d2fdsig11sig31*dsig11deps11 + d2fdsig22sig31*dsig22deps11 + d2fdsig33sig31*dsig33deps11)
                jac[5,1] =      x[6]*(d2fdsig11sig31*dsig11deps22 + d2fdsig22sig31*dsig22deps22 + d2fdsig33sig31*dsig33deps33)
                jac[5,2] =      x[6]*(d2fdsig11sig31*dsig11deps33 + d2fdsig22sig31*dsig22deps33 + d2fdsig33sig31*dsig33deps33)
                jac[5,3] =      x[6]*d2fdsig12sig31*dsig12deps12
                jac[5,4] =      x[6]*d2fdsig23sig31*dsig23deps23
                jac[5,5] = 1. + x[6]*d2fdsig31sig31*dsig31deps31
                jac[5,6] = dfdsig31

                jac[6,0] = dfdsig11*dsig11deps11 + dfdsig22*dsig22deps11 + dfdsig33*dsig33deps11
                jac[6,1] = dfdsig11*dsig11deps22 + dfdsig22*dsig22deps22 + dfdsig33*dsig33deps22
                jac[6,2] = dfdsig11*dsig11deps33 + dfdsig22*dsig22deps33 + dfdsig33*dsig33deps33
                jac[6,3] = dfdsig12*dsig12deps12
                jac[6,4] = dfdsig23*dsig23deps23
                jac[6,5] = dfdsig31*dsig31deps31
                jac[6,6] = dfdlamda

                # Solve system of equations
                dx = np.linalg.solve(jac, -res) # increment of target variables

                # Update x
                x += dx

                # Compute error
                err = np.linalg.norm(dx)

                #print("\tNewton iter.",ii, ": err =", err)

                if err < 1e-7:
                    #input("stop==")
                    break

            # Update strain
            #print("\n")
            eps   += deps
            eps_e = np.array([[x[0], x[3], x[5]],
                              [x[3], x[1], x[4]],
                              [x[5], x[4], x[2]]])
            eps_p = eps - eps_e
            lamda = lamda + x[6]

            # Update stress
            sigma = np.tensordot(D, eps_e)

            # D tangent operator has to be implemented ######################################

            ddsdde = FthTensor2Voigt(D)
            if self.Dim == 2:
                idx = [0,1,3]
                D = ddsdde[ix_(idx, idx)]
            elif self.Dim ==3:
                D = ddsdde
            else:
                assert 0, "Dimension Check"

        GPstate.eps_e  = self.Tensor2Voigt(eps_e)
        GPstate.eps_p  = self.Tensor2Voigt(eps_p)
        GPstate.lamda  = lamda
        GPstate.stress = self.Tensor2Voigt(sigma,'stress')
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

        eps_e_init  = self.Voigt2Tensor( GPstate.eps_e )
        eps_p_init  = self.Voigt2Tensor( GPstate.eps_p )
        deps        = self.Voigt2Tensor( GPstate.deps  )
        sigma_init  = self.Voigt2Tensor( GPstate.stress ,'stress' )
        lamda = GPstate.lamda   # plastic multiplier

        HSP = self.InitialaPressure()

        eps_e_tr = eps_e_init + deps # trial strain
        eps_e_tr_principal_mag, eps_e_tr_principal_vec = np.linalg.eigh(eps_e_tr)

        eps_e_tr1 = eps_e_tr_principal_mag[0]
        eps_e_tr2 = eps_e_tr_principal_mag[1]
        eps_e_tr3 = eps_e_tr_principal_mag[2]

        n1 = eps_e_tr_principal_vec[:,0]
        n2 = eps_e_tr_principal_vec[:,1]
        n3 = eps_e_tr_principal_vec[:,2]

        # [2] Compute trial stress
        sigma_tr_principal_mag = np.inner(Ce_principal, eps_e_tr_principal_mag)

        sigma_tr1 = sigma_tr_principal_mag[0] + HSP[0,0]
        sigma_tr2 = sigma_tr_principal_mag[1] + HSP[1,1]
        sigma_tr3 = sigma_tr_principal_mag[2] + HSP[2,2]

        sigma_tr = sigma_tr1*np.tensordot(n1,n1,axes=0) + sigma_tr2*np.tensordot(n2,n2,axes=0) + sigma_tr3*np.tensordot(n3,n3,axes=0)

        YF = self.PlasticityModel.f(sigma_tr1, sigma_tr2, sigma_tr3, lamda) # yield surface

        if YF <= 0.:
            #print(">> Elastic!")
            # Update stress & strain
            sigma  = sigma_tr
            dsigma = sigma_tr - sigma_init

            eps_e = eps_e_init + deps
            eps_p = eps_p_init
            eps   = eps_e + eps_p

            deps_e = eps_e - eps_e_init
            deps_p = eps_p - eps_p_init

            ddsdde = FthTensor2Voigt(D)

            if self.Dim == 2:
                idx = [0,1,3]
                D = ddsdde[ix_(idx, idx)]
            elif self.Dim ==3:
                D = ddsdde
            else:
                assert 0, "Dimension Check"
            #print(D)
        else:
            #print(">> Plastic!")
            # Initialize variables
            eps_e_principal_mag, eps_e_principal_vec = np.linalg.eigh(eps_e_init)

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

                sigma1_current = sigma1_current + HSP[0,0]
                sigma2_current = sigma2_current + HSP[1,1]
                sigma3_current = sigma3_current + HSP[2,2]

                # Current lamda
                lamda_current = lamda + x[3]

                # Update derivatives
                # >> First order derivatives
                dfdsig1, dfdsig2, dfdsig3, dfdlamda = self.PlasticityModel.df(sigma1_current, sigma2_current, sigma3_current, lamda_current)

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
                jac[3,3] = dfdlamda

                # Solve system of equations
                dx = np.linalg.solve(jac, -res) # increment of target variables

                # Update x
                x = x + dx
                if x[3] < 0. :
                    #print("residual values")
                    #print(x)
                    x[3] = 0.
                    #input("")


                # Compute error
                err = np.linalg.norm(dx)

                #print(" Newton iter.",ii, ": err =", err)

                if err < 1e-7:
                    break

            #input("==============")
            # Update strain
            #deps_e = (eps_e_tr1-x[0])*np.tensordot(n1,n1,axes=0)\
                   #+ (eps_e_tr2-x[1])*np.tensordot(n2,n2,axes=0)\
                   #+ (eps_e_tr3-x[2])*np.tensordot(n3,n3,axes=0)
            #deps_p = deps - deps_e

            ##eps   = eps_e_init + eps_p_init + deps

            ##eps_e = x[0]*np.tensordot(n1,n1,axes=0)\
                  ##+ x[1]*np.tensordot(n2,n2,axes=0)\
                  ##+ x[2]*np.tensordot(n3,n3,axes=0)

            ##deps_e = eps_e - eps_e_init

            ##eps_p  = eps - eps_e

            ##deps_p = eps_p - eps_p_init

            ##lamda = lamda + x[3]

            ### Update stress
            ##dsigma1 = a*deps_e[0] + b*deps_e[1] + b*deps_e[2]
            ##dsigma2 = b*deps_e[0] + a*deps_e[1] + b*deps_e[2]
            ##dsigma3 = b*deps_e[0] + b*deps_e[1] + a*deps_e[2]

            ##dsigma  = dsigma1*np.tensordot(n1,n1,axes=0)\
                    ##+ dsigma2*np.tensordot(n2,n2,axes=0)\
                    ##+ dsigma3*np.tensordot(n3,n3,axes=0)

            ###############
            eps   = eps_e_init + eps_p_init + deps
            eps_e = x[0]*np.tensordot(n1,n1,axes=0) + x[1]*np.tensordot(n2,n2,axes=0) + x[2]*np.tensordot(n3,n3,axes=0)
            eps_p = eps - eps_e
            lamda = lamda + x[3]

            # Update stress
            sigma1 = a*x[0] + b*x[1] + b*x[2]
            sigma2 = b*x[0] + a*x[1] + b*x[2]
            sigma3 = b*x[0] + b*x[1] + a*x[2]
            sigma  = sigma1*np.tensordot(n1,n1,axes=0)\
                   + sigma2*np.tensordot(n2,n2,axes=0)\
                   + sigma3*np.tensordot(n3,n3,axes=0)
            sigma += HSP


            deps_e = eps_e - eps_e_init
            deps_p = eps_p - eps_p_init
            dsigma = sigma - sigma_init
            ###############


            #sigma  = sigma_tr
            #dsigma = sigma_tr - sigma_init

            #eps_e = eps_e_init + deps
            #eps_p = eps_p_init
            #eps   = eps_e + eps_p

            #deps_e = eps_e - eps_e_init
            #deps_p = eps_p - eps_p_init



            #eps   = eps_e_init + eps_p_init + deps
            #eps_e = x[0]*np.tensordot(n1,n1,axes=0)\
                  #+ x[1]*np.tensordot(n2,n2,axes=0)\
                  #+ x[2]*np.tensordot(n3,n3,axes=0)
            #eps_p = eps - eps_e
            #lamda = lamda + x[3]

            ## Update stress
            #sigma1 = a*x[0] + b*x[1] + b*x[2]
            #sigma2 = b*x[0] + a*x[1] + b*x[2]
            #sigma3 = b*x[0] + b*x[1] + a*x[2]
            #sigma  = sigma1*np.tensordot(n1,n1,axes=0) + sigma2*np.tensordot(n2,n2,axes=0) + sigma3*np.tensordot(n3,n3,axes=0)

            #deps_e = eps_e - eps_e_init
            #deps_p = eps_p - eps_p_init
            #dsigma = sigma - sigma_init


            #aAB = np.zeros((3,3))
            #mab = np.zeros((3,3,3,3))
            #dsigdepse = np.zeros((3,3))
            #dsigdepse[0,0] = dsig1depse1
            #dsigdepse[0,1] = dsig1depse2
            #dsigdepse[0,2] = dsig1depse3
            #dsigdepse[1,0] = dsig2depse1
            #dsigdepse[1,1] = dsig2depse2
            #dsigdepse[1,2] = dsig2depse3
            #dsigdepse[2,0] = dsig3depse1
            #dsigdepse[2,1] = dsig3depse2
            #dsigdepse[2,2] = dsig3depse3

            ##aAB = np.linalg.solve(jac[0:3,0:3],dsigdepse)
            #aAB = dsigdepse
            #bAB = np.zeros((3,3))
            #eps_e_tr = np.zeros((3))
            #sig_e_cur = np.zeros((3))
            #eps_e_tr[0] = eps_e_tr1
            #eps_e_tr[1] = eps_e_tr2
            #eps_e_tr[2] = eps_e_tr3
            #sig_e_cur[0] = sigma[0,0]
            #sig_e_cur[1] = sigma[1,1]
            #sig_e_cur[2] = sigma[2,2]

            #for ii in range(3):
                #for jj in range(3):
                    #if(np.abs(eps_e_tr[jj] - eps_e_tr[ii])) < 1e-10:
                        #bAB[ii,jj] = 0.0
                    #else:
                        #bAB[ii,jj] = (sig_e_cur[jj] - sig_e_cur[ii])/(eps_e_tr[jj] - eps_e_tr[ii])

            #for ii in range(3):
                #for jj in range(3):
                    #na = eps_e_principal_vec[:,ii].reshape(-1,1)
                    #nb = eps_e_principal_vec[:,jj].reshape(-1,1)
                    #mab[ii,jj,:,:] = na @ nb.T

            #Calgo = np.zeros((3,3,3,3,))

            #for ii in range(3):
                #for jj in range(3):
                    #for mm in range(3):
                        #for nn in range(3):
                            #Calgo[ii,jj,mm,nn] += aAB[0,0]*mab[0,0,ii,jj]*mab[0,0,mm,nn]\
                                                #+ aAB[0,1]*mab[0,0,ii,jj]*mab[1,1,mm,nn]\
                                                #+ aAB[0,2]*mab[0,0,ii,jj]*mab[2,2,mm,nn]\
                                                #+ aAB[1,0]*mab[1,1,ii,jj]*mab[0,0,mm,nn]\
                                                #+ aAB[1,1]*mab[1,1,ii,jj]*mab[1,1,mm,nn]\
                                                #+ aAB[1,2]*mab[1,1,ii,jj]*mab[2,2,mm,nn]\
                                                #+ aAB[2,1]*mab[2,2,ii,jj]*mab[0,0,mm,nn]\
                                                #+ aAB[2,2]*mab[2,2,ii,jj]*mab[1,1,mm,nn]\
                                                #+ aAB[2,2]*mab[2,2,ii,jj]*mab[2,2,mm,nn]\
                                            #+0.5* bAB[0,1]*(mab[0,1,ii,jj]*mab[0,1,mm,nn]+mab[0,1,ii,jj]*mab[1,0,mm,nn])\
                                            #+0.5* bAB[0,2]*(mab[0,2,ii,jj]*mab[0,2,mm,nn]+mab[0,1,ii,jj]*mab[2,0,mm,nn])\
                                            #+0.5* bAB[1,0]*(mab[1,0,ii,jj]*mab[1,0,mm,nn]+mab[1,0,ii,jj]*mab[0,1,mm,nn])\
                                            #+0.5* bAB[1,2]*(mab[1,2,ii,jj]*mab[1,2,mm,nn]+mab[1,2,ii,jj]*mab[2,1,mm,nn])\
                                            #+0.5* bAB[2,0]*(mab[2,0,ii,jj]*mab[2,0,mm,nn]+mab[2,0,ii,jj]*mab[0,2,mm,nn])\
                                            #+0.5* bAB[2,1]*(mab[2,1,ii,jj]*mab[2,1,mm,nn]+mab[2,1,ii,jj]*mab[1,2,mm,nn])

            ddsdde = FthTensor2Voigt(D)
            #ddsdde = FthTensor2Voigt(Calgo)
            if self.Dim == 2:
                idx = [0,1,3]
                D = ddsdde[ix_(idx, idx)]
            elif self.Dim ==3:
                D = ddsdde
            else:
                assert 0, "Dimension Check"

        GPstate.eps_e   = self.Tensor2Voigt(eps_e)
        GPstate.eps_p   = self.Tensor2Voigt(eps_p)
        GPstate.lamda   = lamda
        GPstate.dstress = self.Tensor2Voigt(dsigma,'stress')
        GPstate.stress  = self.Tensor2Voigt(sigma,'stress')
        GPstate.deps_e  = self.Tensor2Voigt(deps_e)
        GPstate.deps_p  = self.Tensor2Voigt(deps_p)

        #print(GPstate.stress)
        GPstate.D      = D
        return

