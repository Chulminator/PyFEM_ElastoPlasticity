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
from util.EigCalc import *
import importlib
from scipy import linalg
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import spsolve


class GPstate:
    eps_e: np.array
    eps_p: np.array
    deps: np.array
    deps_e: np.array
    deps_p: np.array
    dstress: np.array
    lamda: float
    stress: np.array


def ElementSetUp(Program, Model):

    #ModelElement =
    Dim          = Model.Dim

    for Element in Model.Element:
        Edof = Element.Connectivity

        NodeTarget = []
        for Node in Edof:
            NodeTarget.append(Node.Coord)

        Init_tensor = []
        Init_tensor1 = []
        Init_tensor2 = []
        Init_tensor3 = []
        B_matrix  = []
        Jacc_Elem = []
        Area = 0.

        for gp in Element.GP:
            if Model.Dim   == 2:
                Jacc, B = Strain(Element, NodeTarget, gp[0], gp[1])

            elif Model.Dim ==3:
                Jacc, B = Strain(Element, NodeTarget, gp[0], gp[1], gp[2])

            else:
                assert False, "Check Model.Dim"

            Area = Area + Jacc * gp[-1]
            B_matrix.append(B)
            Jacc_Elem.append(Jacc)
            if Element.Dim ==2:
                if Model.twoD == "planestrain":
                    tmp_tensor = np.zeros((4))
                elif Model.twoD == "planestress":
                    tmp_tensor = np.zeros((3))
                else:
                    assert False, "Check Model.twoD"
                tmp_tensor1 = np.zeros((3))
                tmp_tensor1[0] = Model.HSP
                tmp_tensor1[1] = Model.HSP
            elif Element.Dim ==3:
                tmp_tensor = np.zeros((6))
                tmp_tensor1 = np.zeros((6))
                tmp_tensor1[0] = Model.HSP
                tmp_tensor1[1] = Model.HSP
                tmp_tensor1[2] = Model.HSP
            else:
                assert False, "Check Fem.Dimension"
            Init_tensor.append(np.copy(tmp_tensor))
            Init_tensor1.append(0.0)
            Init_tensor2.append(tmp_tensor1)
            Init_tensor3.append(np.copy(tmp_tensor))

        Element.B_matrix = B_matrix
        Element.Jacc     = Jacc_Elem
        Element.Area     = Area

        Element.GPstrain_e = Init_tensor
        Element.GPlamda    = Init_tensor1
        Element.GPstress   = Init_tensor2
        Element.GPstrain_p = Init_tensor3
    return


def ConstructStiffness(Model):
    Dim    = Model.Dim

    for node in Model.Node:
        node.F_int = np.zeros(Dim)

    for ind1, Element in enumerate(Model.Element):
        G_Edof = Element.G_Edof
        E_K = np.zeros([G_Edof, G_Edof])
        B    = Element.B_matrix
        Jacc = Element.Jacc
        GPstress = Element.GPstress
        GPeps_e  = Element.GPstrain_e
        GPeps_p  = Element.GPstrain_p
        GPlamda  = Element.GPlamda
        GPF_int  = np.zeros(G_Edof)
        du       = np.zeros(G_Edof)

        for ind, Node in enumerate(Element.Connectivity):
            du[ind*Dim:(ind+1)*Dim] = Node.du
            #Node.F_int = np.zeros(Dim)

        for ind, gp in enumerate(Element.GP):
            B_GP    = B[ind]
            Jacc_GP = Jacc[ind]
            GPdata  =  GPstate()
            GPdata.eps_e  = GPeps_e[ind]
            GPdata.eps_p  = GPeps_p[ind]
            GPdata.lamda  = GPlamda[ind]
            GPdata.deps   = B_GP @ du
            GPdata.stress = GPstress[ind]
            
            Element.ConstitutiveModel.ReturnMapping(GPdata, Model.GlobalNRStep)

            GPeps_e[ind]  = GPdata.eps_e
            GPeps_p[ind]  = GPdata.eps_p
            GPstress[ind] =  GPdata.stress
            GPlamda[ind]  =  GPdata.lamda

            if Model.Dim == 2:
                GP_K = B_GP.T @ GPdata.D @ B_GP * Model.width * Jacc_GP * gp[-1]
                GPF_int += B_GP.T @ GPstress[ind] * Model.width * Jacc_GP * gp[-1]
            elif Model.Dim == 3:
                GP_K = B_GP.T @ GPdata.D @ B_GP * Jacc_GP * gp[-1]
                GPF_int += B_GP.T @ GPstress[ind] * Jacc_GP * gp[-1]
            else:
                assert False, "Dimension error"

            E_K += GP_K

        for ind, Node in enumerate(Element.Connectivity):
            Node.F_int += GPF_int[ind*Dim:(ind+1)*Dim]

        Element.Stiffness = E_K
    return


def Solve(Model, Global_K):
    # https://mae.ufl.edu/nkim/egm6352/Chap2.pdf
    # Kim, N.H., 2014. Introduction to nonlinear finite element analysis. Springer Science & Business Media.
    # [[K11 K12],  [[d1],     [[Fout1],
    #  [K21 K22]]   [d2]]   =  [Fout2]]
    #
    # K12 * d2 + K11 * d1 = F1
    # d1 = K11 \ ( Fout1 - K12 * d2 )

    Dim = Model.Dim
    IndexBCN =[]
    IndexBCE =[]

    F_ext = np.zeros(Model.NNode * Dim)
    F_int = np.zeros(Model.NNode * Dim)
    u     = np.zeros(Model.NNode * Dim)
    du    = np.zeros(Model.NNode * Dim)
    for ii, node in enumerate(Model.Node):
        for jj in range(Dim):
            u[ii*Dim+jj]     = node.u[jj]
            F_ext[ii*Dim+jj] = node.F_ext[jj]
            F_int[ii*Dim+jj] = node.F_int[jj]

            if np.isnan( node.BC_E[jj] ):
                IndexBCN.append(ii*Dim+jj)
            else:
                IndexBCE.append(ii*Dim+jj)

    u1 = np.copy(u)

    Sliced_K11 = Global_K[ix_(IndexBCN, IndexBCN)]

    F_total = F_ext[IndexBCN] - F_int[IndexBCN]

    #du[IndexBCN] = linalg.solve(Sliced_K11, F_total)
    Sliced_K11_csc = csc_matrix(Sliced_K11)
    du[IndexBCN] = spsolve(Sliced_K11_csc, F_total)


    u1[IndexBCN] = u[IndexBCN]  + du[IndexBCN]
    du[IndexBCE] = u1[IndexBCE] -  u[IndexBCE]

    u  = np.copy(u1)

    for ii, node in enumerate( Model.Node):
        for jj in range(Dim):
            node.u[jj]  = u[ii*Dim + jj]
            node.du[jj] = du[ii*Dim + jj]

    return np.linalg.norm(F_total)**2/(1.+np.linalg.norm(F_ext[IndexBCN])**2)



class ConstitutiveLaw():
    def __init__(self, Model, Element):
        self.twoD = Model.twoD
        self.ReturnMappingMode = Model.ReturnMappingMode
        tmp = importlib.import_module(Model.PlasticModel)
        self.PlasticityModel = tmp.MyPlasticity(Model)
        self.GlobalNRStep = Model.GlobalNRStep
        self.HSP = Model.HSP
        self.Dim  = Model.Dim

        self.K       = Element.K
        self.mu      = Element.mu
        self.lam     = Element.lamda
        self.E       = Element.E
        self.nu      = Element.nu
        self.MatProp = Element.MatProp

    def Voigt2Tensor(self, voigt, flag='strain'):
        # https://www.comsol.com/blogs/what-is-the-difference-between-plane-stress-and-plane-strain
        Tensor = np.zeros((3,3))
        if flag.lower() == 'strain':
            if self.twoD == 'planestrain':
                Tensor[0,0] = voigt[0]
                Tensor[0,1] = voigt[2]*0.5
                Tensor[1,0] = voigt[2]*0.5
                Tensor[1,1] = voigt[1]
                assert False, "this function is not available"

            elif self.twoD == 'planestress':
                nu = self.nu
                Tensor[0,0] = voigt[0]
                Tensor[0,1] = voigt[2]*0.5
                Tensor[1,0] = voigt[2]*0.5
                Tensor[1,1] = voigt[1]
                Tensor[2,2] = -nu/(1.0-nu)*(Tensor[0,0] + Tensor[1,1])
                assert False, "this function is not available"
            else:
                if self.Dim == 3:
                    Tensor[0,0] = voigt[0]
                    Tensor[1,1] = voigt[1]
                    Tensor[2,2] = voigt[2]

                    Tensor[0,1] = voigt[3]*0.5 # xy
                    Tensor[1,0] = voigt[3]*0.5

                    Tensor[1,2] = voigt[4]*0.5 # yz
                    Tensor[2,1] = voigt[4]*0.5

                    Tensor[0,2] = voigt[5]*0.5 # zx
                    Tensor[2,0] = voigt[5]*0.5
                else:
                    print("self.twoD", self.twoD)
                    print("self.Dim",  self.Dim)
                    assert False,"Check Fem.twoD or Fem.Dimension"

        elif flag.lower() == 'stress':
            if self.twoD == 'planestrain':
                nu = self.nu
                Tensor[0,0] = voigt[0]
                Tensor[0,1] = voigt[2]
                Tensor[1,0] = voigt[2]
                Tensor[1,1] = voigt[1]
                Tensor[2,2] = nu*(Tensor[0,0] + Tensor[1,1])
                assert False, "this function is not available"

            elif self.twoD == 'planestress':
                nu = self.nu
                Tensor[0,0] = voigt[0]
                Tensor[0,1] = voigt[2]
                Tensor[1,0] = voigt[2]
                Tensor[1,1] = voigt[1]
                assert False, "this function is not available"
                
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
                    assert False,"Check Fem.twoD or Fem.Dimension"
        return Tensor


    def Tensor2Voigt(self, Tensor, flag='strain'):        
        if self.Dim ==2:
            voigt = np.zeros((3))
            if flag == 'strain':
                voigt[0] = Tensor[0,0]
                voigt[1] = Tensor[1,1]
                voigt[2] = Tensor[0,1] + Tensor[1,0]
            elif flag == 'stress':
                voigt[0] = Tensor[0,0]
                voigt[1] = Tensor[1,1]
                voigt[2] = Tensor[0,1]
            else:
                assert False ,"check the flag in Tensor2Voigt"
        elif self.Dim ==3:
            voigt = np.zeros((6))
            if flag == 'strain':
                voigt[0] = Tensor[0,0]
                voigt[1] = Tensor[1,1]
                voigt[2] = Tensor[2,2]
                voigt[3] = Tensor[0,1]+Tensor[1,0]
                voigt[4] = Tensor[1,2]+Tensor[2,1]
                voigt[5] = Tensor[2,0]+Tensor[0,2]

            elif flag == 'stress':
                voigt[0] = Tensor[0,0]
                voigt[1] = Tensor[1,1]
                voigt[2] = Tensor[2,2]
                voigt[3] = Tensor[0,1]
                voigt[4] = Tensor[1,2]
                voigt[5] = Tensor[2,0]
            else:
                assert False ,"check the flag in Tensor2Voigt"
        else:
            assert False,"Check the dimension"
        return voigt


    def InitialPressure(self):
        return np.diag([self.HSP, self.HSP, self.HSP])


    def ReturnMapping(self, GPstate, GlobalNRStep):
        self.GlobalNRStep = GlobalNRStep
        if (self.ReturnMappingMode.lower() == 'eigenspace'):
            if self.twoD == 'planestrain':
                self.ReturnMappingEigenSpacePE(GPstate)

            elif self.twoD == 'planestress':
                self.ReturnMappingEigenSpacePS(GPstate)

            elif self.Dim == 3:
                self.ReturnMappingEigenSpace(GPstate)

            else:
                assert False, "Check dimension or Fem.twoD"


        elif(self.ReturnMappingMode == 'tensorspace' ):
            if self.Dim == 3:
                self.ReturnMappingTensorSpace(GPstate)

            elif self.twoD == 'planestrain':
                self.ReturnMappingTensorSpacePE(GPstate)

            else:
                assert False, "PlaneStress is not implemented for tensor space - ReturnMappingMode"

        elif(self.ReturnMappingMode == 'pqspace' ):
            self.ReturnMappingPQSpace(GPstate)
        else:
            assert False, "ReturnMappingMode"
        return



    def ReturnMappingTensorSpacePE(self, GPstate):
        #https://getfem.org/userdoc/model_plasticity_small_strain.html
        #assert False, "Not implemented"
        K   = self.K
        lam = self.lam
        mu  = self.mu

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

        D = lam*IxI + 2.*mu*II

        eps_e_init  = GPstate.eps_e
        eps_p_init  = GPstate.eps_p
        sigma_init  = GPstate.stress
        deps        = GPstate.deps
        lamda       = GPstate.lamda   # plastic equivalent strain

        eps_e_init_tensor      = np.zeros((3,3))
        eps_e_init_tensor[0,0] = eps_e_init[0]
        eps_e_init_tensor[1,1] = eps_e_init[1]
        eps_e_init_tensor[2,2] = eps_e_init[3]

        eps_e_init_tensor[0,1] = eps_e_init[2]*0.5
        eps_e_init_tensor[1,0] = eps_e_init[2]*0.5

        deps_tensor      = np.zeros((3,3))
        deps_tensor[0,0] = deps[0]
        deps_tensor[1,1] = deps[1]
        deps_tensor[0,1] = deps[2]*0.5
        deps_tensor[1,0] = deps[2]*0.5  # deps_[2,2] shoule be zero <- plane strain

        HSP_tensor       = self.InitialPressure()

        #[1] compute trial elastic strain
        eps_e_tr_tensor  = eps_e_init_tensor + deps_tensor # trial strain

        eps_e_11_tr = eps_e_tr_tensor[0,0]
        eps_e_22_tr = eps_e_tr_tensor[1,1]
        eps_e_33_tr = eps_e_tr_tensor[2,2]
        eps_e_12_tr = eps_e_tr_tensor[0,1]

        # [2] Compute trial stress
        sigma_tr_tensor = np.tensordot(D, eps_e_tr_tensor) + HSP_tensor

        sig11_tr = sigma_tr_tensor[0,0]
        sig22_tr = sigma_tr_tensor[1,1]
        sig33_tr = sigma_tr_tensor[2,2]
        sig12_tr = sigma_tr_tensor[0,1]

        YF = self.PlasticityModel.f(sig11_tr, sig22_tr, sig33_tr, sig12_tr, 0., 0., lamda) # yield surface

        if YF <= 0. or self.GlobalNRStep == 0:
            sigma = np.zeros(3)
            sigma[0] = sigma_tr_tensor[0,0]
            sigma[1] = sigma_tr_tensor[1,1]
            sigma[2] = sigma_tr_tensor[0,1]

            dsigma = sigma - sigma_init

            eps_e       = eps_e_init
            eps_e[0:3] += deps
            eps_p       = eps_p_init
            eps         = eps_e + eps_p

            deps_e = eps_e - eps_e_init
            deps_p = eps_p - eps_p_init

            ddsdde = FthTensor2Voigt(D)
            if self.Dim == 2:
                idx = [0,1,3]
                D = ddsdde[ix_(idx, idx)]
            elif self.Dim ==3:
                D = ddsdde
            else:
                assert False, "Dimension Check"

        else:
            #print(">> Plastic!")
            # Initialize variables
            eps_e_11 = eps_e_init_tensor[0,0]
            eps_e_22 = eps_e_init_tensor[1,1]
            eps_e_33 = eps_e_init_tensor[2,2]
            eps_e_12 = eps_e_init_tensor[0,1]

            dlamda = 0.

            x = np.zeros(5) # target variables
            x[0] = eps_e_11
            x[1] = eps_e_22
            x[2] = eps_e_33
            x[3] = eps_e_12
            x[4] = dlamda

            # Newton-Raphson iteration (return mapping)
            for ii in range(20):
                # Initialize residual and jacobian
                res = np.zeros(5)
                jac = np.zeros((5,5))

                # Current strain
                eps_e_11_current = x[0]
                eps_e_22_current = x[1]
                eps_e_33_current = x[2]
                eps_e_12_current = x[3]

                # Current stress
                sig11_current = (lam+2.*mu)*eps_e_11_current + lam*eps_e_22_current + lam*eps_e_33_current
                sig22_current = lam*eps_e_11_current + (lam+2.*mu)*eps_e_22_current + lam*eps_e_33_current
                sig33_current = lam*eps_e_11_current + lam*eps_e_22_current + (lam+2.*mu)*eps_e_33_current
                sig12_current = 2.*mu*eps_e_12_current

                # Current lamda
                lamda_current = lamda + x[4]

                # Update derivatives
                # >> First order derivatives
                dfdsig11, dfdsig22, dfdsig33, dfdsig12, dfdsig23, dfdsig31, dfdlamda = self.PlasticityModel.df(sig11_current, sig22_current, sig33_current, sig12_current, 0., 0., lamda_current)

                # >> Second order derivatives
                d2fdsig11sig11, d2fdsig11sig22, d2fdsig11sig33, d2fdsig11sig12, d2fdsig11sig23, d2fdsig11sig31,  \
                                d2fdsig22sig22, d2fdsig22sig33, d2fdsig22sig12, d2fdsig22sig23, d2fdsig22sig31,  \
                                                d2fdsig33sig33, d2fdsig33sig12, d2fdsig33sig23, d2fdsig33sig31,  \
                                                                d2fdsig12sig12, d2fdsig12sig23, d2fdsig12sig31,  \
                                                                                d2fdsig23sig23, d2fdsig23sig31,  \
                                                                                                d2fdsig31sig31 = \
                self.PlasticityModel.df2(sig11_current, sig22_current, sig33_current, sig12_current, 0., 0.)

                # Update residual
                res[0] = x[0] - eps_e_11_tr + x[4]*dfdsig11
                res[1] = x[1] - eps_e_22_tr + x[4]*dfdsig22
                res[2] = x[2] - eps_e_33_tr + x[4]*dfdsig33
                res[3] = x[3] - eps_e_12_tr + x[4]*dfdsig12
                res[4] = self.PlasticityModel.f(sig11_current, sig22_current, sig33_current, sig12_current, 0., 0., lamda_current)

                # Update jacobian ***
                jac[0,0] = 1. + x[4]*(d2fdsig11sig11*dsig11deps11 + d2fdsig11sig22*dsig22deps11 + d2fdsig11sig33*dsig33deps11)
                jac[0,1] =      x[4]*(d2fdsig11sig11*dsig11deps22 + d2fdsig11sig22*dsig22deps22 + d2fdsig11sig33*dsig33deps22)
                jac[0,2] =      x[4]*(d2fdsig11sig11*dsig11deps33 + d2fdsig11sig22*dsig22deps33 + d2fdsig11sig33*dsig33deps33)
                jac[0,3] =      x[4]*d2fdsig11sig12*dsig12deps12
                jac[0,4] = dfdsig11

                jac[1,0] =      x[4]*(d2fdsig11sig22*dsig11deps11 + d2fdsig22sig22*dsig22deps11 + d2fdsig22sig33*dsig33deps11)
                jac[1,1] = 1. + x[4]*(d2fdsig11sig22*dsig11deps22 + d2fdsig22sig22*dsig22deps22 + d2fdsig22sig33*dsig33deps22)
                jac[1,2] =      x[4]*(d2fdsig11sig22*dsig11deps33 + d2fdsig22sig22*dsig22deps33 + d2fdsig22sig33*dsig33deps33)
                jac[1,3] =      x[4]*d2fdsig22sig12*dsig12deps12
                jac[1,4] = dfdsig22

                jac[2,0] =      x[4]*(d2fdsig11sig33*dsig11deps11 + d2fdsig22sig33*dsig22deps11 + d2fdsig33sig33*dsig33deps11)
                jac[2,1] =      x[4]*(d2fdsig11sig33*dsig11deps22 + d2fdsig22sig33*dsig22deps22 + d2fdsig33sig33*dsig33deps22)
                jac[2,2] = 1. + x[4]*(d2fdsig11sig33*dsig11deps33 + d2fdsig22sig33*dsig22deps33 + d2fdsig33sig33*dsig33deps33)
                jac[2,3] =      x[4]*d2fdsig33sig12*dsig12deps12
                jac[2,4] = dfdsig33

                jac[3,0] =      x[4]*(d2fdsig11sig12*dsig11deps11 + d2fdsig22sig12*dsig22deps11 + d2fdsig33sig12*dsig33deps11)
                jac[3,1] =      x[4]*(d2fdsig11sig12*dsig11deps22 + d2fdsig22sig12*dsig22deps22 + d2fdsig33sig12*dsig33deps22)
                jac[3,2] =      x[4]*(d2fdsig11sig12*dsig11deps33 + d2fdsig22sig12*dsig22deps33 + d2fdsig33sig12*dsig33deps33)
                jac[3,3] = 1. + x[4]*d2fdsig12sig12*dsig12deps12
                jac[3,4] = dfdsig12

                jac[4,0] = dfdsig11*dsig11deps11 + dfdsig22*dsig22deps11 + dfdsig33*dsig33deps11
                jac[4,1] = dfdsig11*dsig11deps22 + dfdsig22*dsig22deps22 + dfdsig33*dsig33deps22
                jac[4,2] = dfdsig11*dsig11deps33 + dfdsig22*dsig22deps33 + dfdsig33*dsig33deps33
                jac[4,3] = dfdsig12*dsig12deps12
                jac[4,4] = dfdlamda

                # Solve system of equations
                dx = linalg.solve(jac, -res) # increment of target variables

                # Update x
                x += dx

                # Compute error
                err = np.linalg.norm(dx)

                #print("\tNewton iter.",ii, ": err =", err)

                if err < 1e-7:
                    #input("stop==")
                    break


            #eps_tensor   = eps_e_init_tensor + deps_tensor + np.array([[eps_p_init[0], eps_p_init[2], 0.]
                                                                #[eps_p_init[2], eps_p_init[1], 0.]
                                                                #[ 0.,           0., eps_p_init[3]]])

            eps   = eps_e_init[0:3] + eps_p_init[0:3] + deps

            eps_e_tensor = np.array([[x[0], x[3], 0.],
                                     [x[3], x[1], 0.],
                                     [0., 0., x[2]]])

            eps_e = np.array([x[0], x[1], x[3]*2., x[2] ])

            eps_p      =  np.zeros(4)
            eps_p[0:3] =  eps - eps_e[0:3]
            eps_p[3]   = -eps_e[3]
            lamda      = lamda + x[4]

            # Update stress
            sigma_tensor  = np.tensordot(D,eps_e_tensor)
            sigma_tensor += HSP_tensor

            sigma    = np.zeros(3)
            sigma[0] = sigma_tensor[0,0]
            sigma[1] = sigma_tensor[1,1]
            sigma[2] = sigma_tensor[0,1]

            deps_e = eps_e - eps_e_init
            deps_p = eps_p - eps_p_init
            dsigma = sigma - sigma_init


            ############## Cep #############
            #deveps = eps_e_tensor - 1.0/3.0 *np.trace(eps_e_tensor)*I
            #n      = deveps / np.linalg.norm(deveps,'fro')
            #nxn  = tensor_oMult(n,n)
            #D = K*IxI+2.0*mu*(II-1.0/3.0*IxI - nxn)
            ############## Cep #############
            # D tangent operator has to be implemented ######################################
            ddsdde = FthTensor2Voigt(D)
            if self.Dim == 2:
                idx = [0,1,3]
                D = ddsdde[ix_(idx, idx)]
            elif self.Dim ==3:
                D = ddsdde
            else:
                assert 0, "Dimension Check"

        GPstate.eps_e   = eps_e
        GPstate.eps_p   = eps_p
        GPstate.lamda   = lamda
        GPstate.dstress = dsigma
        GPstate.stress  = sigma
        GPstate.deps_e  = deps_e
        GPstate.deps_p  = deps_p

        GPstate.D      = D
        return




    def ReturnMappingEigenSpacePE(self, GPstate):
        # plane strain in eigen space
        K  = self.K
        mu = self.mu
        nu = self.nu
        lam= self.lam

        I   = np.eye(3)
        II  = identity_4(3)
        IxI = tensor_oMult(I,I)

        D = lam*IxI + 2.*mu*II

        aa = K + (4/3)*mu
        bb = K - (2/3)*mu

        dsig1depse1 = aa
        dsig1depse2 = bb
        dsig1depse3 = bb
        dsig2depse1 = bb
        dsig2depse2 = aa
        dsig2depse3 = bb
        dsig3depse1 = bb
        dsig3depse2 = bb
        dsig3depse3 = aa

        eps_e_init  = GPstate.eps_e
        eps_p_init  = GPstate.eps_p
        sigma_init  = GPstate.stress
        deps        = GPstate.deps
        lamda       = GPstate.lamda   # plastic equivalent strain

        eps_e_init_tensor      = np.zeros((3,3))
        eps_e_init_tensor[0,0] = eps_e_init[0]
        eps_e_init_tensor[1,1] = eps_e_init[1]
        eps_e_init_tensor[2,2] = eps_e_init[3]

        eps_e_init_tensor[0,1] = eps_e_init[2]*0.5
        eps_e_init_tensor[1,0] = eps_e_init[2]*0.5

        deps_tensor      = np.zeros((3,3))
        deps_tensor[0,0] = deps[0]
        deps_tensor[1,1] = deps[1]
        deps_tensor[0,1] = deps[2]*0.5
        deps_tensor[1,0] = deps[2]*0.5  # deps_[2,2] shoule be zero <- plane strain

        HSP_tensor       = self.InitialPressure()

        #[1] compute trial elastic strain
        eps_e_tr_tensor  = eps_e_init_tensor + deps_tensor # trial strain

        eps_e_tr_mag, eps_e_tr_vec = np.linalg.eigh(eps_e_tr_tensor[0:2,0:2])
        eps_e_tr_mag, eps_e_tr_vec = SortEig2D(eps_e_tr_mag, eps_e_tr_vec)

        eps_e_tr1 = eps_e_tr_mag[0]
        eps_e_tr2 = eps_e_tr_mag[1]
        eps_e_tr3 = eps_e_tr_tensor[2,2]

        k1 = eps_e_tr_vec[:,0]
        k2 = eps_e_tr_vec[:,1]

        sigma1_tr = aa*eps_e_tr1 + bb*eps_e_tr2 + bb*eps_e_tr3
        sigma2_tr = bb*eps_e_tr1 + aa*eps_e_tr2 + bb*eps_e_tr3
        sigma3_tr = bb*eps_e_tr1 + bb*eps_e_tr2 + aa*eps_e_tr3

        YF = self.PlasticityModel.f(sigma1_tr, sigma2_tr, sigma3_tr, lamda) # yield surface

        if YF <= 0. or self.GlobalNRStep == 0:
            sigma_tr_tensor = sigma1_tr*np.tensordot(k1,k1,axes=0)\
                            + sigma2_tr*np.tensordot(k2,k2,axes=0)

            sigma = np.zeros(3)
            sigma[0] = sigma_tr_tensor[0,0]
            sigma[1] = sigma_tr_tensor[1,1]
            sigma[2] = sigma_tr_tensor[0,1]

            dsigma = sigma - sigma_init

            eps_e       = eps_e_init
            eps_e[0:3] += deps
            eps_p       = eps_p_init
            eps         = eps_e + eps_p

            deps_e = eps_e - eps_e_init
            deps_p = eps_p - eps_p_init

            ddsdde = FthTensor2Voigt(D)
            if self.Dim == 2:
                idx = [0,1,3]
                D = ddsdde[ix_(idx, idx)]
            elif self.Dim ==3:
                D = ddsdde
            else:
                assert False, "Dimension Check"

        else:
            #print(">> Plastic!")
            # Initialize variables
            eps_e_init_mag, eps_e_init_vec = np.linalg.eigh(eps_e_init_tensor[0:2,0:2])
            eps_e_init_mag, eps_e_init_vec = SortEig2D(eps_e_init_mag, eps_e_init_vec)

            k1 = eps_e_init_vec[:,0]
            k2 = eps_e_init_vec[:,1]

            eps_e1 = eps_e_init_mag[0]
            eps_e2 = eps_e_init_mag[1]
            eps_e3 = eps_e_init_tensor[2,2]

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
                # eps_e3_current = 0.
                sigma1_current = aa*eps_e1_current + bb*eps_e2_current + bb*eps_e3_current
                sigma2_current = bb*eps_e1_current + aa*eps_e2_current + bb*eps_e3_current
                sigma3_current = bb*eps_e1_current + bb*eps_e2_current + aa*eps_e3_current

                sigma1_current = sigma1_current + HSP_tensor[0,0]
                sigma2_current = sigma2_current + HSP_tensor[1,1]
                sigma3_current = sigma3_current + HSP_tensor[2,2]

                # Current lamda
                lamda_current = lamda + x[3]

                # Update derivatives
                # >> First order derivatives
                dfdsig1, dfdsig2, dfdsig3, dfdlamda = self.PlasticityModel.df(sigma1_current, sigma2_current, sigma3_current, lamda_current)

                # >> Second order derivatives
                d2fdsig1dsig1, d2fdsig2dsig2, d2fdsig3dsig3, d2fdsig1dsig2, d2fdsig2dsig3, d2fdsig3dsig1\
                    = self.PlasticityModel.df2(sigma1_current, sigma2_current, sigma3_current)

                # Update residual
                res[0] = x[0] - eps_e_tr1 + x[3]*dfdsig1
                res[1] = x[1] - eps_e_tr2 + x[3]*dfdsig2
                res[2] = x[2] - eps_e_tr3 + x[3]*dfdsig3
                res[3] = self.PlasticityModel.f(sigma1_current, sigma2_current, sigma3_current, lamda_current)

                # Update Jacobian ***
                jac[0,0] = 1. + x[3]*(d2fdsig1dsig1*dsig1depse1 + d2fdsig1dsig2*dsig2depse1 + d2fdsig3dsig1*dsig3depse1)
                jac[0,1] =      x[3]*(d2fdsig1dsig1*dsig1depse2 + d2fdsig1dsig2*dsig2depse2 + d2fdsig3dsig1*dsig3depse2)
                jac[0,2] =      x[3]*(d2fdsig1dsig1*dsig1depse3 + d2fdsig1dsig2*dsig2depse3 + d2fdsig3dsig1*dsig3depse3)
                jac[0,3] = dfdsig1

                jac[1,0] =      x[3]*(d2fdsig1dsig2*dsig1depse1 + d2fdsig2dsig2*dsig2depse1 + d2fdsig2dsig3*dsig3depse1)
                jac[1,1] = 1. + x[3]*(d2fdsig1dsig2*dsig1depse2 + d2fdsig2dsig2*dsig2depse2 + d2fdsig2dsig3*dsig3depse2)
                jac[1,2] =      x[3]*(d2fdsig1dsig2*dsig1depse3 + d2fdsig2dsig2*dsig2depse3 + d2fdsig2dsig3*dsig3depse3)
                jac[1,3] = dfdsig2

                jac[2,0] =      x[3]*(d2fdsig3dsig1*dsig1depse1 + d2fdsig2dsig3*dsig2depse1 + d2fdsig3dsig3*dsig3depse1)
                jac[2,1] =      x[3]*(d2fdsig3dsig1*dsig1depse2 + d2fdsig2dsig3*dsig2depse2 + d2fdsig3dsig3*dsig3depse2)
                jac[2,2] = 1. + x[3]*(d2fdsig3dsig1*dsig1depse3 + d2fdsig2dsig3*dsig2depse3 + d2fdsig3dsig3*dsig3depse3)
                jac[2,3] = dfdsig3

                jac[3,0] = dfdsig1*dsig1depse1 + dfdsig2*dsig2depse1 + dfdsig3*dsig3depse1
                jac[3,1] = dfdsig1*dsig1depse2 + dfdsig2*dsig2depse2 + dfdsig3*dsig3depse2
                jac[3,2] = dfdsig1*dsig1depse3 + dfdsig2*dsig2depse3 + dfdsig3*dsig3depse3
                jac[3,3] = dfdlamda

                # Solve system of equations
                dx = linalg.solve(jac, -res) # increment of target variables

                # Update x
                x = x + dx

                # Compute error
                err = np.linalg.norm(dx)

                if err < 1e-12:
                    break

            #if x[3] < 0 or dx[3] < 0:
            if x[3] < 0 < 0:
                print("\t","x")
                print("\t",x)
                print(" Newton iter.",ii, ": err =", err)
                print("WARNING: Maybe you should check the algorithmic tangent IN ./include/ElastoPlasticity.py")
            eps   = eps_e_init[0:3] + eps_p_init[0:3] + deps

            eps_e_tensor = np.zeros((3,3))
            eps_e_tensor[0:2,0:2] = x[0]*np.tensordot(k1,k1,axes=0)\
                                  + x[1]*np.tensordot(k2,k2,axes=0)
            eps_e_tensor[2,2]     = x[2]

            eps_e    = np.zeros(4)
            eps_e[0] = eps_e_tensor[0,0]
            eps_e[1] = eps_e_tensor[1,1]
            eps_e[2] = eps_e_tensor[0,1] * eps_e_tensor[1,0]
            eps_e[3] = eps_e_tensor[2,2]

            eps_p      =  np.zeros(4)
            eps_p[0:3] =  eps - eps_e[0:3]
            eps_p[3]   = -eps_e[3]
            lamda      = lamda + x[3]

            # Update stress
            sigma_tensor  = np.tensordot(D,eps_e_tensor)
            sigma_tensor += HSP_tensor
            sigma    = np.zeros(3)
            sigma[0] = sigma_tensor[0,0]
            sigma[1] = sigma_tensor[1,1]
            sigma[2] = sigma_tensor[0,1]

            ###########################
            ###########################################
            #This does not work for the plasticity force control. which makes sense
            ############# Cep #############
            #deveps = eps_e_tensor - 1.0/3.0 *np.trace(eps_e_tensor)*I
            #n      = deveps / np.linalg.norm(deveps,'fro')
            #nxn    = tensor_oMult(n,n)
            #D      = K*IxI+2.0*mu*(II-1.0/3.0*IxI - nxn)
            ############# Cep #############
            ######################################################################

            deps_e = eps_e - eps_e_init
            deps_p = eps_p - eps_p_init
            dsigma = sigma - sigma_init
            ###############


            ddsdde = FthTensor2Voigt(D)
            if self.Dim == 2:
                idx = [0,1,3]
                D = ddsdde[ix_(idx, idx)]
            elif self.Dim ==3:
                D = ddsdde
            else:
                assert 0, "Dimension Check"

        #print(GPstate.stress)
        GPstate.D      = D
        GPstate.eps_e   = eps_e
        GPstate.eps_p   = eps_p
        GPstate.lamda   = lamda
        GPstate.dstress = dsigma
        GPstate.stress  = sigma
        GPstate.deps_e  = deps_e
        GPstate.deps_p  = deps_p

        return



    def ReturnMappingEigenSpacePS(self, GPstate):
        K   = self.K
        mu  = self.mu
        E   = self.E
        nu  = self.nu
        lam = self.lam

        I   = np.eye(3)
        II  = identity_4(3)
        IxI = tensor_oMult(I,I)

        dsig1depse1 = lam + 2. *mu
        dsig1depse2 = lam
        dsig2depse1 = lam
        dsig2depse2 = lam + 2. *mu
        dsig3depse3 = mu

        D = np.array([[1., nu, 0.],
                      [nu, 1., 0.],
                      [0., 0., (1.-nu)*0.5]])*E/(1.-nu**2)

        eps_e_init   = GPstate.eps_e
        eps_p_init   = GPstate.eps_p
        sigma_init   = GPstate.stress
        deps         = GPstate.deps
        lamda        = GPstate.lamda   # plastic multiplier


        eps_e_tr     = eps_e_init + deps # trial strain

        # [2] Compute trial stress
        HSP          = np.array([self.HSP, self.HSP, 0.])

        sigma_tr     = D @ eps_e_tr + HSP
        tmp      = np.zeros((2,2))
        tmp[0,0] = sigma_tr[0]
        tmp[1,1] = sigma_tr[1]
        tmp[0,1] = sigma_tr[2]
        tmp[1,0] = sigma_tr[2]
        sigma_tr_mag, sigma_tr_vec = np.linalg.eigh(tmp)
        sigma_tr_mag, sigma_tr_vec = SortEig2D(sigma_tr_mag, sigma_tr_vec)


        YF = self.PlasticityModel.f(sigma_tr_mag[0], sigma_tr_mag[1], self.HSP,  lamda) # yield surface
        if YF <= 0. or self.GlobalNRStep == 0:
            #print(">> Elastic!")
            # Update stress & strain
            sigma  = sigma_tr
            dsigma = sigma_tr - sigma_init

            eps_e = eps_e_init + deps
            eps_p = eps_p_init
            eps   = eps_e + eps_p

            deps_e = eps_e - eps_e_init
            deps_p = eps_p - eps_p_init

        else:
            #print(">> Plastic!")
            # Initialize variables
            tmp = np.zeros((2,2))
            tmp[0,0] = eps_e_init[0]
            tmp[1,1] = eps_e_init[1]
            tmp[0,1] = eps_e_init[2]*0.5
            tmp[1,0] = eps_e_init[2]*0.5
            eps_e_init_mag, eps_e_init_vec = np.linalg.eigh(tmp)
            eps_e_init_mag, eps_e_init_vec = SortEig2D(eps_e_init_mag, eps_e_init_vec)

            k1 = eps_e_init_vec[:,0]
            k2 = eps_e_init_vec[:,1]

            eps_e1 = eps_e_init_mag[0]
            eps_e2 = eps_e_init_mag[1]

            dlamda  = 0.

            tmp = np.zeros((2,2))
            tmp[0,0] = eps_e_tr[0]
            tmp[1,1] = eps_e_tr[1]
            tmp[0,1] = eps_e_tr[2]*0.5
            tmp[1,0] = eps_e_tr[2]*0.5

            eps_e_tr_mag, eps_e_tr_vec = np.linalg.eigh(tmp)
            eps_e_tr_mag, eps_e_tr_vec = SortEig2D(eps_e_tr_mag, eps_e_tr_vec)

            eps_e_tr1 = eps_e_tr_mag[0]
            eps_e_tr2 = eps_e_tr_mag[1]

            x = np.zeros(3) # target variables
            x[0] = eps_e1
            x[1] = eps_e2
            x[2] = dlamda

            # Newton-Raphson iteration (return mapping)
            for ii in range(20):
                # Initialize residual and jacobian
                res = np.zeros(3)
                jac = np.zeros((3,3))

                # Current strain
                eps_e1_current = x[0]
                eps_e2_current = x[1]

                # Current stress
                eps_e_current = np.zeros((2,2))
                eps_e_current = x[0]*np.tensordot(k1,k1,axes=0) + x[1]*np.tensordot(k2,k2,axes=0)

                sigma_current = D @ np.array([eps_e_current[0,0], eps_e_current[1,1], eps_e_current[0,1]*2. ])

                sigma1_current = sigma_current[0] + HSP[0]
                sigma2_current = sigma_current[1] + HSP[1]

                # Current lamda
                lamda_current = lamda + x[2]

                # Update derivatives
                # >> First order derivatives
                dfdsig1, dfdsig2, _, dfdlamda = self.PlasticityModel.df(sigma1_current, sigma2_current, self.HSP, lamda_current)

                # >> Second order derivatives
                d2fdsig1dsig1, d2fdsig2dsig2, _, d2fdsig1dsig2, _, _\
                    = self.PlasticityModel.df2(sigma1_current, sigma2_current, self.HSP)

                # Update residual
                res[0] = x[0] - eps_e_tr1 + x[2]*dfdsig1
                res[1] = x[1] - eps_e_tr2 + x[2]*dfdsig2
                res[2] = self.PlasticityModel.f(sigma1_current, sigma2_current, self.HSP, lamda_current)

                # Update Jacobian ***
                jac[0,0] = 1. + x[2]*(d2fdsig1dsig1*dsig1depse1 + d2fdsig1dsig2*dsig2depse1)
                jac[0,1] =      x[2]*(d2fdsig1dsig1*dsig1depse2 + d2fdsig1dsig2*dsig2depse2)
                jac[0,2] = dfdsig1

                jac[1,0] =      x[2]*(d2fdsig1dsig2*dsig1depse1 + d2fdsig2dsig2*dsig2depse1)
                jac[1,1] = 1. + x[2]*(d2fdsig1dsig2*dsig1depse2 + d2fdsig2dsig2*dsig2depse2)
                jac[1,2] = dfdsig2

                jac[2,0] = dfdsig1*dsig1depse1 + dfdsig2*dsig2depse1
                jac[2,1] = dfdsig1*dsig1depse2 + dfdsig2*dsig2depse2
                jac[2,2] = dfdlamda

                # Solve system of equations
                dx = linalg.solve(jac, -res) # increment of target variables

                # Update x
                x = x + dx

                # Compute error
                err = np.linalg.norm(dx)

                #print(" Newton iter.",ii, ": err =", err)

                if err < 1e-7:
                    break


            eps   = eps_e_init + eps_p_init + deps
            tmp   = x[0]*np.tensordot(k1,k1,axes=0) + x[1]*np.tensordot(k2,k2,axes=0)
            eps_e = np.zeros(3)
            eps_e[0] = tmp[0,0]
            eps_e[1] = tmp[1,1]
            eps_e[2] = tmp[0,1] + tmp[1,0]
            eps_p = eps - eps_e
            lamda = lamda + x[2]

            # Update stress
            sigma  = D @ eps_e
            sigma += HSP

            deps_e = eps_e - eps_e_init
            deps_p = eps_p - eps_p_init
            dsigma = sigma - sigma_init

        GPstate.eps_e   = eps_e
        GPstate.eps_p   = eps_p
        GPstate.lamda   = lamda
        GPstate.dstress = dsigma
        GPstate.stress  = sigma
        GPstate.deps_e  = deps_e
        GPstate.deps_p  = deps_p

        #print(GPstate.stress)
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

        HSP = self.InitialPressure()

        eps_e_tr = eps_e_init + deps # trial strain
        eps_e_tr_principal_mag, eps_e_tr_principal_vec = np.linalg.eigh(eps_e_tr)
        eps_e_tr_principal_mag, eps_e_tr_principal_vec = SortEig(eps_e_tr_principal_mag, eps_e_tr_principal_vec)

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

        sigma_tr = sigma_tr1*np.tensordot(n1,n1,axes=0)\
                 + sigma_tr2*np.tensordot(n2,n2,axes=0)\
                 + sigma_tr3*np.tensordot(n3,n3,axes=0)

        YF = self.PlasticityModel.f(sigma_tr1, sigma_tr2, sigma_tr3, lamda) # yield surface

        if YF <= 0. or self.GlobalNRStep == 0:
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
            eps_e_init_mag, eps_e_init_vec = np.linalg.eigh(eps_e_init)
            eps_e_init_mag, eps_e_init_vec = SortEig(eps_e_init_mag, eps_e_init_vec)

            eps_e1 = eps_e_init_mag[0]
            eps_e2 = eps_e_init_mag[1]
            eps_e3 = eps_e_init_mag[2]
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
                jac[0,0] = 1. + x[3]*(d2fdsig1dsig1*dsig1depse1 + d2fdsig1dsig2*dsig2depse1 + d2fdsig3dsig1*dsig3depse1)
                jac[0,1] =      x[3]*(d2fdsig1dsig1*dsig1depse2 + d2fdsig1dsig2*dsig2depse2 + d2fdsig3dsig1*dsig3depse2)
                jac[0,2] =      x[3]*(d2fdsig1dsig1*dsig1depse3 + d2fdsig1dsig2*dsig2depse3 + d2fdsig3dsig1*dsig3depse3)
                jac[0,3] = dfdsig1

                jac[1,0] =      x[3]*(d2fdsig1dsig2*dsig1depse1 + d2fdsig2dsig2*dsig2depse1 + d2fdsig2dsig3*dsig3depse1)
                jac[1,1] = 1. + x[3]*(d2fdsig1dsig2*dsig1depse2 + d2fdsig2dsig2*dsig2depse2 + d2fdsig2dsig3*dsig3depse2)
                jac[1,2] =      x[3]*(d2fdsig1dsig2*dsig1depse3 + d2fdsig2dsig2*dsig2depse3 + d2fdsig2dsig3*dsig3depse3)
                jac[1,3] = dfdsig2

                jac[2,0] =      x[3]*(d2fdsig3dsig1*dsig1depse1 + d2fdsig2dsig3*dsig2depse1 + d2fdsig3dsig3*dsig3depse1)
                jac[2,1] =      x[3]*(d2fdsig3dsig1*dsig1depse2 + d2fdsig2dsig3*dsig2depse2 + d2fdsig3dsig3*dsig3depse2)
                jac[2,2] = 1. + x[3]*(d2fdsig3dsig1*dsig1depse3 + d2fdsig2dsig3*dsig2depse3 + d2fdsig3dsig3*dsig3depse3)
                jac[2,3] = dfdsig3

                jac[3,0] = dfdsig1*dsig1depse1 + dfdsig2*dsig2depse1 + dfdsig3*dsig3depse1
                jac[3,1] = dfdsig1*dsig1depse2 + dfdsig2*dsig2depse2 + dfdsig3*dsig3depse2
                jac[3,2] = dfdsig1*dsig1depse3 + dfdsig2*dsig2depse3 + dfdsig3*dsig3depse3
                jac[3,3] = dfdlamda


                # Solve system of equations
                dx = linalg.solve(jac, -res) # increment of target variables

                # Update x
                x = x + dx

                # Compute error
                err = np.linalg.norm(dx)

                #print(" Newton iter.",ii, ": err =", err)

                if err < 1e-7:
                    break

            ###############
            eps   = eps_e_init + eps_p_init + deps
            eps_e = x[0]*np.tensordot(n1,n1,axes=0)\
                  + x[1]*np.tensordot(n2,n2,axes=0)\
                  + x[2]*np.tensordot(n3,n3,axes=0)
            eps_p = eps - eps_e

            # Update stress
            sigma1 = a*x[0] + b*x[1] + b*x[2]
            sigma2 = b*x[0] + a*x[1] + b*x[2]
            sigma3 = b*x[0] + b*x[1] + a*x[2]
            sigma  = sigma1*np.tensordot(n1,n1,axes=0)\
                   + sigma2*np.tensordot(n2,n2,axes=0)\
                   + sigma3*np.tensordot(n3,n3,axes=0)
            sigma += HSP
            lamda = lamda + x[3]


            deps_e = eps_e - eps_e_init
            deps_p = eps_p - eps_p_init
            dsigma = sigma - sigma_init
            ###############

            ######################################################################
            ############## Cep #############
            deveps = eps_e - 1.0/3.0 *np.trace(eps_e)*I
            n      = deveps / np.linalg.norm(deveps,'fro')
            nxn  = tensor_oMult(n,n)
            D = K*IxI+2.0*mu*(II-1.0/3.0*IxI - nxn)
            ############## Cep #############
            ######################################################################

            ddsdde = FthTensor2Voigt(D)
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
        GPstate.D      = D
        return


    def ReturnMappingPQSpace(self, GPstate):
        K  = self.K
        lam= self.lam
        mu = self.mu

        I   = np.eye(3)
        II  = identity_4(3)
        IxI = tensor_oMult(I,I)

        D = lam*IxI + 2.*mu*II

        eps_e_init  = self.Voigt2Tensor( GPstate.eps_e )
        eps_p_init  = self.Voigt2Tensor( GPstate.eps_p )
        sigma_init  = self.Voigt2Tensor( GPstate.stress ,'stress' )
        deps        = self.Voigt2Tensor( GPstate.deps  )
        lamda = GPstate.lamda   # plastic multiplier

        HSP = self.InitialPressure()

        # [1] Compute trial strain
        eps_e_tr   = eps_e_init + deps
        eps_e_v_tr = np.trace(eps_e_tr)
        eps_e_e_tr = eps_e_tr - 1./3.*eps_e_v_tr*I
        eps_e_s_tr = np.sqrt(2./3.) * np.linalg.norm(eps_e_e_tr)

        if eps_e_s_tr < 1e-15:
            n_hat = np.zeros((3,3))
        else:
            n_hat = eps_e_e_tr / np.linalg.norm(eps_e_e_tr)

        # [2] Compute trial stress
        sigma_tr = sigma_init + np.tensordot(D, deps) + HSP
        p_tr     = (1/3)*np.trace(sigma_tr)
        s_tr     = sigma_tr - p_tr*I
        q_tr     = np.sqrt(3./2.)*np.linalg.norm(s_tr)

        YF = self.PlasticityModel.f(p_tr, q_tr, lamda) # yield surface


        if YF <= 0. or self.GlobalNRStep == 0:
            sigma  = sigma_tr
            dsigma = sigma_tr - sigma_init

            eps_e = eps_e_tr
            eps_p = eps_p_init
            eps   = eps_e + eps_p

            deps_e = eps_e - eps_e_init
            deps_p = eps_p - eps_p_init

            ddsdde = FthTensor2Voigt(D)
            if self.Dim ==2:
                idx = [0,1,3]
                D = ddsdde[ix_(idx, idx)]
            elif self.Dim ==3:
                D = ddsdde

        else:
            #print("plastic range")
            eps_e_v = np.trace(eps_e_init)
            eps_e_e = eps_e_init - (1/3)*eps_e_v*I
            eps_e_s = np.sqrt(2/3)*np.linalg.norm(eps_e_e)
            dlamda  = 0

            x = np.zeros(3) # target variables
            x[0] = eps_e_v
            x[1] = eps_e_s
            x[2] = dlamda
            for ii in range(20):

                # Initialize residual and jacobian
                res = np.zeros(3)
                jac = np.zeros((3,3))

                # Current strain
                eps_e_v = x[0]
                eps_e_s = x[1]
                eps_e = (1/3)*eps_e_v*I + np.sqrt(3/2)*eps_e_s*n_hat

                sigma = np.tensordot(D, eps_e) + HSP
                p     = (1./3.)*np.trace(sigma)
                s     = sigma - p*I
                q     = np.sqrt(3./2.)*np.linalg.norm(s)

                # Current lamda
                lamda_current = lamda + x[2]

                # >> First order derivatives
                dfdp, dfdq, dfdlamda = self.PlasticityModel.df(p, q, lamda_current)

                dpdepsev = K
                dpdepses = 0

                dqdepsev = 0
                dqdepses = 3*mu

                res[0] = x[0] - eps_e_v_tr + x[2]*dfdp
                res[1] = x[1] - eps_e_s_tr + x[2]*dfdq
                res[2] = self.PlasticityModel.f(p, q, lamda) # yield surface

                jac[0,0] = 1.
                jac[0,1] = 0.
                jac[0,2] = dfdp
                jac[1,0] = 0.
                jac[1,1] = 1.
                jac[1,2] = dfdq
                jac[2,0] = dfdp*dpdepsev + dfdq*dqdepsev
                jac[2,1] = dfdp*dpdepses + dfdq*dqdepses
                jac[2,2] = dfdlamda

                dx = linalg.solve(jac, -res) # increment of target variables

                # Update x
                x = x + dx

                # Compute error
                err = np.linalg.norm(dx)

                if err < 1e-10:
                    break

            ###############
            eps   = eps_e_init + eps_p_init + deps
            eps_e_v = x[0]
            eps_e_s = x[1]
            eps_e   = (1/3)*eps_e_v*I + np.sqrt(3/2)*eps_e_s*n_hat
            eps_p = eps - eps_e
            lamda = lamda + x[2]

            # Update stress
            sigma = np.tensordot(D, eps_e)
            sigma += HSP

            deps_e = eps_e - eps_e_init
            deps_p = eps_p - eps_p_init
            dsigma = sigma - sigma_init
            ###############
            ############### Cep #############
            nxn   = tensor_oMult(n_hat,n_hat)
            #gamma =  2 * mu * x[2] / np.linalg.norm(s_tr)
            D = K*IxI+2.0*mu*(II-1.0/3.0*IxI - nxn)
            #D = K*IxI+2.0*mu*gamma*(II-1.0/3.0*IxI - nxn)
            #D = K*IxI+2.0*mu*(II-1.0/3.0*IxI - nxn) - 2.0*mu*gamma*(II-1.0/3.0*IxI - nxn)
            ############### Cep #############


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
        GPstate.D      = D
        return

    def ReturnMappingTensorSpace(self, GPstate):
        #assert Falser, "Not implemented"
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

        D = lam*IxI + 2.*mu*II

        eps_e_init  = self.Voigt2Tensor( GPstate.eps_e )
        eps_p_init  = self.Voigt2Tensor( GPstate.eps_p )
        sigma_init  = self.Voigt2Tensor( GPstate.stress ,'stress' )
        deps        = self.Voigt2Tensor( GPstate.deps  )
        lamda = GPstate.lamda   # plastic multiplier

        HSP = self.InitialPressure()

        # [1] Compute trial strain
        eps_e_tr = eps_e_init + deps

        eps_e_11_tr = eps_e_tr[0,0]
        eps_e_22_tr = eps_e_tr[1,1]
        eps_e_33_tr = eps_e_tr[2,2]
        eps_e_12_tr = eps_e_tr[0,1]
        eps_e_23_tr = eps_e_tr[1,2]
        eps_e_31_tr = eps_e_tr[2,0]

        # [2] Compute trial stress
        sigma_tr = np.tensordot(D, eps_e_tr) + HSP

        sig11_tr = sigma_tr[0,0]
        sig22_tr = sigma_tr[1,1]
        sig33_tr = sigma_tr[2,2]
        sig12_tr = sigma_tr[0,1]
        sig23_tr = sigma_tr[1,2]
        sig31_tr = sigma_tr[2,0]

        YF = self.PlasticityModel.f(sig11_tr, sig22_tr, sig33_tr, sig12_tr, sig23_tr, sig31_tr, lamda) # yield surface

        if YF <= 0. or self.GlobalNRStep == 0:
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

            if self.Dim ==2:
                idx = [0,1,3]
                D = ddsdde[ix_(idx, idx)]
            elif self.Dim ==3:
                D = ddsdde
        else:
            #print(">> Plastic!")
            # Initialize variables
            eps_e_11 = eps_e_init[0,0]
            eps_e_22 = eps_e_init[1,1]
            eps_e_33 = eps_e_init[2,2]
            eps_e_12 = eps_e_init[0,1]
            eps_e_23 = eps_e_init[1,2]
            eps_e_31 = eps_e_init[2,0]
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


                # Update derivatives
                # >> First order derivatives
                dfdsig11, dfdsig22, dfdsig33, dfdsig12, dfdsig23, dfdsig31, dfdlamda = self.PlasticityModel.df(sig11_current, sig22_current, sig33_current, sig12_current, sig23_current, sig31_current, lamda_current)

                ## >> Second order derivatives
                #d2fdsig11sig11, d2fdsig11sig22, d2fdsig11sig33, d2fdsig11sig12, d2fdsig11sig23, d2fdsig11sig31,  \
                                #d2fdsig22sig22, d2fdsig22sig33, d2fdsig22sig12, d2fdsig22sig23, d2fdsig22sig31,  \
                                                #d2fdsig33sig33, d2fdsig33sig12, d2fdsig33sig23, d2fdsig33sig31,  \
                                                                #d2fdsig12sig12, d2fdsig12sig23, d2fdsig12sig31,  \
                                                                                #d2fdsig23sig23, d2fdsig23sig31,  \
                                                                                                #d2fdsig31sig31 = \
                #self.PlasticityModel.df2(sig11_current, sig22_current, sig33_current, sig12_current, sig23_current, sig31_current, lamda_current)

                d2fdsig11sig11, d2fdsig11sig22, d2fdsig11sig33, d2fdsig11sig12, d2fdsig11sig23, d2fdsig11sig31,  \
                                d2fdsig22sig22, d2fdsig22sig33, d2fdsig22sig12, d2fdsig22sig23, d2fdsig22sig31,  \
                                                d2fdsig33sig33, d2fdsig33sig12, d2fdsig33sig23, d2fdsig33sig31,  \
                                                                d2fdsig12sig12, d2fdsig12sig23, d2fdsig12sig31,  \
                                                                                d2fdsig23sig23, d2fdsig23sig31,  \
                                                                                                d2fdsig31sig31 = \
                self.PlasticityModel.df2(sig11_current, sig22_current, sig33_current, sig12_current, sig23_current, sig31_current)

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
                dx = linalg.solve(jac, -res) # increment of target variables

                # Update x
                x += dx

                # Compute error
                err = np.linalg.norm(dx)
                #print("\tNewton iter.",ii, ": err =", err)

                if err < 1e-7:
                    #input("stop==")
                    break

            eps   = eps_e_init + eps_p_init + deps
            eps_e = np.array([[x[0], x[3], x[5]],
                              [x[3], x[1], x[4]],
                              [x[5], x[4], x[2]]])
            eps_p = eps - eps_e
            lamda = lamda + x[6]

            # Update stress
            sigma  = np.tensordot(D, eps_e)
            sigma += HSP

            deps_e = eps_e - eps_e_init
            deps_p = eps_p - eps_p_init
            dsigma = sigma - sigma_init


            ############## Cep #############
            deveps = eps_e - 1.0/3.0 *np.trace(eps_e)*I
            n      = deveps / np.linalg.norm(deveps,'fro')
            nxn  = tensor_oMult(n,n)
            D = K*IxI+2.0*mu*(II-1.0/3.0*IxI - nxn)
            ############## Cep #############
            # D tangent operator has to be implemented ######################################
            ddsdde = FthTensor2Voigt(D)
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

        GPstate.D      = D
        return

