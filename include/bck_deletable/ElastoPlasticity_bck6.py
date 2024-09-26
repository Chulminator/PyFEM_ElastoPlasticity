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
    Element.GPlamda = []
    #Element.GPlamda1 = []
    Element.GPstress = []
    #i = 0
    for Edof in Element.Connectivity:
        NodeTarget = []
        for dof in Edof:
            NodeTarget.append(Node.Coord[Node.Id.index(dof)])
        Init_tensor = []
        Init_tensor1 = []
        Init_tensor2 = []
        Init_tensor3 = []
        B_matrix = []
        Jacc_Elem = []
        Area = 0.

        for gp in Fem.GP:
            if Fem.Dimension == 2:
                Jacc, B = Strain(Fem, NodeTarget, gp[0], gp[1])
            elif Fem.Dimension ==3:
                Jacc, B = Strain(Fem, NodeTarget, gp[0], gp[1], gp[2])
            else:
                assert False, "Check Fem.Dimension"

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
                assert False, "Check Fem.Dimension"
            Init_tensor.append(np.copy(tmp_tensor))
            Init_tensor3.append(np.copy(tmp_tensor))
            Init_tensor1.append(0.0)
            Init_tensor2.append(tmp_tensor1)
        Element.B_matrix.append(B_matrix)
        Element.Jacc.append(Jacc_Elem)
        Element.Area.append(Area)

        Element.GPstrain_e.append(Init_tensor)
        Element.GPlamda.append(Init_tensor1)
        Element.GPstress.append(Init_tensor2)
        Element.GPstrain_p.append(Init_tensor3)

    return


def ConstructStiffness(Fem, Node, Element, ConstitutiveModel):
    G_Edof = Fem.G_Edof
    Dim    = Fem.Dimension

    Element.Stiffness  =[]
    Node.F_int = np.zeros_like(Node.F_int)
    for ind1, Edof in enumerate(Element.Connectivity):
        #print("-"*50)
        #print("Element id: ",ind1)
        ENdof = ElementalNodeDof(Fem, Node, Edof)
        E_K = np.zeros([G_Edof, G_Edof])
        #Elemlamda    = []
        B    = Element.B_matrix[ind1]
        Jacc = Element.Jacc[ind1]
        GPstress= Element.GPstress[ind1]
        GPeps_e = Element.GPstrain_e[ind1]
        GPeps_p = Element.GPstrain_p[ind1]
        GPlamda = Element.GPlamda[ind1]
        GPF_int = np.zeros((len(ENdof)))

        #print(ind1)
        for ind, gp in enumerate(Fem.GP):
            #print("\tGauss point: ",ind)
            B_GP    = B[ind]
            Jacc_GP = Jacc[ind]
            GPdata  =  GPstate()
            GPdata.eps_e  = GPeps_e[ind]
            GPdata.eps_p  = GPeps_p[ind]
            GPdata.lamda  = GPlamda[ind]
            GPdata.deps   = B_GP @ Node.du[ENdof]
            GPdata.stress = GPstress[ind]
            
            ConstitutiveModel.ReturnMapping(GPdata)

            #GPstress[ind] += GPdata.dstress # -> not working
            GPeps_e[ind]  = GPdata.eps_e
            GPeps_p[ind]  = GPdata.eps_p

            #GPeps_e[ind]  += GPdata.deps_e
            #GPeps_p[ind]  += GPdata.deps_p
            GPstress[ind] =  GPdata.stress
            GPlamda[ind]  =  GPdata.lamda

            if Fem.Dimension == 2:
                GP_K = B_GP.T @ GPdata.D @ B_GP * Fem.width * Jacc_GP * gp[-1]
                GPF_int += B_GP.T @ GPstress[ind] * Fem.width * Jacc_GP * gp[-1]
            elif Fem.Dimension == 3:
                GP_K = B_GP.T @ GPdata.D @ B_GP * Jacc_GP * gp[-1]
                GPF_int += B_GP.T @ GPstress[ind] * Jacc_GP * gp[-1]
            else:
                assert 0, "Dimension error"

            E_K += GP_K
            #Elemlamda.append(GPdata.lamda)
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

    #print("\nNode.F_ext\n",Node.F_ext[IndexBCN])
    #print("\nNode.F_int\n",Node.F_int[IndexBCN])
    #print("\nF_total\n",F_total)

    #print("\nNode.du\n",Node.du[IndexBCN])
    #input("")
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
                Tensor[0,1] = voigt[2]*0.5
                Tensor[1,0] = voigt[2]*0.5
                Tensor[1,1] = voigt[1]

            elif self.twoD == 'planestress':
                nu = self.nu
                Tensor[0,0] = voigt[0]
                Tensor[0,1] = voigt[2]*0.5
                Tensor[1,0] = voigt[2]*0.5
                Tensor[1,1] = voigt[1]
                #Tensor[2,2] = -nu*(Tensor[0,0] + Tensor[1,1])
                Tensor[2,2] = -nu/(1.0-nu)*(Tensor[0,0] + Tensor[1,1])
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
                voigt[3] = Tensor[0,1]*0.5
                voigt[4] = Tensor[1,2]*0.5
                voigt[5] = Tensor[2,0]*0.5
            else:
                assert False ,"check the flag in Tensor2Voigt"
        else:
            assert False,"Check the dimension"
        return voigt

    def InitialaPressure(self):
        return np.diag([self.HSP, self.HSP, self.HSP])

    def ReturnMapping(self, GPstate):
        #print(self.ReturnMappingMode)
        #print(self.twoD)
        #exit(1)
        if (self.ReturnMappingMode.lower() == 'eigenspace'):
            if self.twoD == 'planestrain':
                self.ReturnMappingEigenSpacePE(GPstate)
            elif self.twoD == 'planestress':
                self.ReturnMappingEigenSpacePS(GPstate)
            else:
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


    def ReturnMappingEigenSpacePE(self, GPstate):

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
        sigma_init  = self.Voigt2Tensor( GPstate.stress ,'stress' )
        deps        = self.Voigt2Tensor( GPstate.deps  )
        lamda = GPstate.lamda   # plastic multiplier

        HSP = self.InitialaPressure()

        eps_e_tr = eps_e_init + deps # trial strain
        eps_e_tr_principal_mag, eps_e_tr_principal_vec = np.linalg.eigh(eps_e_tr)
        eps_e_tr_principal_mag, eps_e_tr_principal_vec = SortEig(eps_e_tr_principal_mag, eps_e_tr_principal_vec)

        eps_e_tr1 = eps_e_tr_principal_mag[0]
        eps_e_tr2 = eps_e_tr_principal_mag[1]
        eps_e_tr3 = eps_e_tr_principal_mag[2]

        #n1 = eps_e_tr_principal_vec[:,0]
        #n2 = eps_e_tr_principal_vec[:,1]
        #n3 = eps_e_tr_principal_vec[:,2]

        # [2] Compute trial stress

        #sigma_tr_principal_mag = np.inner(Ce_principal, eps_e_tr_principal_mag)

        #sigma_tr1 = sigma_tr_principal_mag[0] + HSP[0,0]
        #sigma_tr2 = sigma_tr_principal_mag[1] + HSP[1,1]
        #sigma_tr3 = sigma_tr_principal_mag[2] + HSP[2,2]


        #sigma_tr = sigma_tr1*np.tensordot(n1,n1,axes=0) + sigma_tr2*np.tensordot(n2,n2,axes=0) + sigma_tr3*np.tensordot(n3,n3,axes=0)
        sigma_tr = np.tensordot(D,eps_e_tr) + HSP

        sigma_tr_principal_mag, sigma_tr_principal_vec = np.linalg.eigh(sigma_tr)
        sigma_tr1 = sigma_tr_principal_mag[0]
        sigma_tr2 = sigma_tr_principal_mag[1]
        sigma_tr3 = sigma_tr_principal_mag[2]

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
            eps_e_principal_mag, eps_e_principal_vec = np.linalg.eigh(eps_e_init[0:2,0:2])
            eps_e_principal_mag, eps_e_principal_vec = SortEig2D(eps_e_principal_mag, eps_e_principal_vec)
            #print(eps_e_principal_mag, eps_e_principal_vec)
            k1 = eps_e_principal_vec[:,0]
            k2 = eps_e_principal_vec[:,1]

            eps_e1 = eps_e_principal_mag[0]
            eps_e2 = eps_e_principal_mag[1]
            dlamda  = 0.

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
                # eps_e3_current = 0.
                sigma1_current = a*eps_e1_current + b*eps_e2_current
                sigma2_current = b*eps_e1_current + a*eps_e2_current
                sigma3_current = (sigma1_current+ sigma2_current)*nu

                sigma1_current = sigma1_current + HSP[0,0]
                sigma2_current = sigma2_current + HSP[1,1]
                sigma3_current = sigma3_current + HSP[2,2]

                # Current lamda
                lamda_current = lamda + x[2]

                # Update derivatives
                # >> First order derivatives
                dfdsig1, dfdsig2, dfdsig3, dfdlamda = self.PlasticityModel.df(sigma1_current, sigma2_current, sigma3_current, lamda_current)

                # >> Second order derivatives
                d2fdsig1dsig1, d2fdsig2dsig2, d2fdsig3dsig3, d2fdsig1dsig2, d2fdsig2dsig3, d2fdsig3dsig1\
                    = self.PlasticityModel.df2(sigma1_current, sigma2_current, sigma3_current)

      #jac[0,0] = 1 + x[3]*(d2fdsig1dsig1*dsig1depse1 + d2fdsig1dsig2*dsig2depse1 + d2fdsig3dsig1*dsig3depse1)
      #jac[0,1] =     x[3]*(d2fdsig1dsig1*dsig1depse2 + d2fdsig1dsig2*dsig2depse2 + d2fdsig3dsig1*dsig3depse2)
      #jac[0,2] =     x[3]*(d2fdsig1dsig1*dsig1depse3 + d2fdsig1dsig2*dsig2depse3 + d2fdsig3dsig1*dsig3depse3)
      #jac[0,3] = dfdsig1

      #jac[1,0] =     x[3]*(d2fdsig1dsig2*dsig1depse1 + d2fdsig2dsig2*dsig2depse1 + d2fdsig2dsig3*dsig3depse1)
      #jac[1,1] = 1 + x[3]*(d2fdsig1dsig2*dsig1depse2 + d2fdsig2dsig2*dsig2depse2 + d2fdsig2dsig3*dsig3depse2)
      #jac[1,2] =     x[3]*(d2fdsig1dsig2*dsig1depse3 + d2fdsig2dsig2*dsig2depse3 + d2fdsig2dsig3*dsig3depse3)
      #jac[1,3] = dfdsig2

      #jac[2,0] =     x[3]*(d2fdsig3dsig1*dsig1depse1 + d2fdsig2dsig3*dsig2depse1 + d2fdsig3dsig3*dsig3depse1)
      #jac[2,1] =     x[3]*(d2fdsig3dsig1*dsig1depse2 + d2fdsig2dsig3*dsig2depse2 + d2fdsig3dsig3*dsig3depse2)
      #jac[2,2] = 1 + x[3]*(d2fdsig3dsig1*dsig1depse3 + d2fdsig2dsig3*dsig2depse3 + d2fdsig3dsig3*dsig3depse3)
      #jac[2,3] = dfdsig3

      #jac[3,0] = dfdsig1*dsig1depse1 + dfdsig2*dsig2depse1 + dfdsig3*dsig3depse1
      #jac[3,1] = dfdsig1*dsig1depse2 + dfdsig2*dsig2depse2 + dfdsig3*dsig3depse2
      #jac[3,2] = dfdsig1*dsig1depse3 + dfdsig2*dsig2depse3 + dfdsig3*dsig3depse3

                # Update residual
                res[0] = x[0] - eps_e_tr1 + x[2]*dfdsig1
                res[1] = x[1] - eps_e_tr2 + x[2]*dfdsig2
                res[2] = self.PlasticityModel.f(sigma1_current, sigma2_current, sigma3_current, lamda_current)

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
                dx = np.linalg.solve(jac, -res) # increment of target variables

                # Update x
                x = x + dx

                # Compute error
                err = np.linalg.norm(dx)

                #print(" Newton iter.",ii, ": err =", err)

                if err < 1e-7:
                    break


            eps   = eps_e_init + eps_p_init + deps
            eps_e = np.zeros((3,3))
            eps_e[0:2,0:2] = x[0]*np.tensordot(k1,k1,axes=0) + x[1]*np.tensordot(k2,k2,axes=0)
            eps_p = eps - eps_e
            lamda = lamda + x[2]

            # Update stress
            sigma  = np.tensordot(D,eps_e)
            sigma += HSP

            #tmp    = np.tensordot(D,eps_e)
            #tmp1   = np.tensordot(D,eps_e_init)
            ######################################################################
            ############# Cep #############
            deveps = eps_e - 1.0/3.0 *np.trace(eps_e)*I
            n      = deveps / np.linalg.norm(deveps,'fro')
            nxn    = tensor_oMult(n,n)
            D      = K*IxI+2.0*mu*(II-1.0/3.0*IxI - nxn)
            ############# Cep #############
            ######################################################################


            deps_e = eps_e - eps_e_init
            deps_p = eps_p - eps_p_init
            dsigma = sigma - sigma_init
            ###############

            ######################## DELETE ########################
            #print("*="*25)

            print("\nitr \n",ii)
            print("\nerr \n",err)
            print("\nHPS \n",HSP)
            print("\nx[0], x[1] \n",x[0], x[1])
            print("\nk1, k2 \n",k1, k2)
            print("\neps_e_init\n",eps_e_init)
            print("\neps_e\n",eps_e)
            print("\nsigma_init\n",sigma_init)
            print("\nsigma\n",sigma)
            print("\ndsigma\n",dsigma)
            print("\nVM\n",np.sqrt(0.5*((sigma[0,0]-sigma[1,1])**2 + (sigma[1,1]-sigma[2,2])**2 + (sigma[2,2]-sigma[0,0])**2)+ 6.0*(sigma[0,1]**2)    ) )
            input("")
            ######################## DELETE ########################

            ddsdde = FthTensor2Voigt(D)
            #ddsdde = FthTensor2Voigt(Calgo)
            #ddsdde = FthTensor2Voigt(Cep)
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



    def ReturnMappingEigenSpacePS(self, GPstate):
        #print("I am here")
        #exit(1)
        K  = self.K
        mu = self.mu
        E = self.E
        nu = self.nu
        lam= self.lam
        I   = np.eye(3)
        II  = identity_4(3)
        IxI = tensor_oMult(I,I)

        dsig1depse1 = lam + 2. *mu
        dsig1depse2 = lam
        dsig2depse1 = lam
        dsig2depse2 = lam + 2. *mu

        D = np.array([[1., nu, 0.],
                      [nu, 1., 0.],
                      [0., 0., (1.-nu)*0.5]])*E/(1.-nu**2)

        #eps_e_init   = self.Voigt2Tensor( GPstate.eps_e )
        #eps_p_init   = self.Voigt2Tensor( GPstate.eps_p )
        #sigma_init   = self.Voigt2Tensor( GPstate.stress, 'stress')
        #deps         = self.Voigt2Tensor( GPstate.deps )

        eps_e_init   = GPstate.eps_e
        eps_p_init   = GPstate.eps_p
        sigma_init   = GPstate.stress
        deps         = GPstate.deps

        #print( D @ GPstate.eps_e )
        #print( D @ GPstate.eps_e )

        lamda        = GPstate.lamda   # plastic multiplier

        HSP          = np.array([self.HSP, self.HSP, 0.])

        eps_e_tr     = eps_e_init + deps # trial strain

        # [2] Compute trial stress
        sigma_tr     = D @ eps_e_tr + HSP

        YF = self.PlasticityModel.f(sigma_tr[0], sigma_tr[1], 0., lamda) # yield surface

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
            #k1 = np.array([1, 0])
            #k2 = np.array([0, 1])

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
            #print(eps_e_tr_mag)
            print(eps_e_init_vec)
            print(eps_e_tr_vec)
            input("="*20)
            #print("="*40)

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
                # eps_e3_current = 0.
                eps_e_current = np.zeros((2,2))
                eps_e_current[0:2,0:2] = x[0]*np.tensordot(k1,k1,axes=0) + x[1]*np.tensordot(k2,k2,axes=0)
                #eps_e_current[2,2] = -nu/(1.0-nu)*(eps_e_current[0,0] + eps_e_current[1,1])

                sigma_current = D @ np.array([eps_e_current[0,0], eps_e_current[1,1], eps_e_current[0,1]*2. ])

                tmp = np.zeros((2,2))
                tmp[0,0] = sigma_current[0]
                tmp[1,1] = sigma_current[1]
                tmp[0,1] = sigma_current[2]
                tmp[1,0] = sigma_current[2]

                sigma_current_principal_mag, sigma_current_principal_vec = np.linalg.eigh(tmp)
                sigma_current_principal_mag, sigma_current_principal_vec = SortEig2D(sigma_current_principal_mag, sigma_current_principal_vec)

                #sigma1_current = (sigma_current[0] + sigma_current[1])*0.5 + np.sqrt(((sigma_current[0] - sigma_current[1])*0.5)**2 + sigma_current[2]**2 )
                #sigma2_current = (sigma_current[0] + sigma_current[1])*0.5 - np.sqrt(((sigma_current[0] - sigma_current[1])*0.5)**2 + sigma_current[2]**2 )

                sigma1_current = sigma_current_principal_mag[0] +  HSP[0]
                sigma2_current = sigma_current_principal_mag[1] +  HSP[1]

                # Current lamda
                lamda_current = lamda + x[2]

                # Update derivatives
                # >> First order derivatives
                #dfdsig1, dfdsig2, _, dfdlamda = self.PlasticityModel.df(sigma1_current, sigma2_current, 0., lamda_current)

                # >> Second order derivatives
                d2fdsig1dsig1, d2fdsig2dsig2, _, d2fdsig1dsig2, _, _\
                    = self.PlasticityModel.df2(sigma1_current, sigma2_current, 0.)

                # Update residual
                res[0] = x[0] - eps_e_tr1 + x[2]*dfdsig1
                res[1] = x[1] - eps_e_tr2 + x[2]*dfdsig2
                res[2] = self.PlasticityModel.f(sigma1_current, sigma2_current, 0., lamda_current)

                # Update Jacobian ***
                jac[0,0] = 1. + x[2]
                jac[0,1] =      x[2]
                jac[0,2] = dfdsig1

                jac[1,0] =      x[2]
                jac[1,1] = 1. + x[2]
                jac[1,2] = dfdsig2

                jac[2,0] = dfdsig1*dsig1depse1 + dfdsig2*dsig2depse1
                jac[2,1] = dfdsig1*dsig1depse2 + dfdsig2*dsig2depse2
                jac[2,2] = dfdlamda

                # Solve system of equations
                dx = np.linalg.solve(jac, -res) # increment of target variables

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
            #print(eps_e)
            #print(D)
            sigma  = D @ eps_e
            sigma += HSP

            #tmp    = np.tensordot(D,eps_e)
            #tmp1   = np.tensordot(D,eps_e_init)
            ######################################################################
            ############# Cep #############
            #deveps = eps_e - 1.0/3.0 *np.trace(eps_e)*I
            #n      = deveps / np.linalg.norm(deveps,'fro')
            #nxn    = tensor_oMult(n,n)
            #D      = K*IxI+2.0*mu*(II-1.0/3.0*IxI - nxn)
            ############# Cep #############
            ######################################################################


            deps_e = eps_e - eps_e_init
            deps_p = eps_p - eps_p_init
            dsigma = sigma - sigma_init
            ###############

            ######################## DELETE ########################
            #print("*="*25)
            #print("EigenspacePS")
            #print("\nitr \n",ii)
            #print("\nerr \n",err)
            #print("\nHPS \n",HSP)
            #print("\nx[0], x[1] \n",x[0], x[1])
            #print("\nk1, k2 \n",k1, k2)
            #print("\neps_e_init\n",eps_e_init)
            #print("\neps_e\n",eps_e)
            #print("\nsigma_init\n",sigma_init)
            #print("\nsigma\n",sigma)
            #print("\ndsigma\n",dsigma)
            ##print("\nVM\n",np.sqrt(0.5*((sigma[0,0]-sigma[1,1])**2 + (sigma[1,1]-sigma[2,2])**2 + (sigma[2,2]-sigma[0,0])**2)+ 6.0*(sigma[0,1]**2)    ) )
            #print("\nVM\n",np.sqrt(0.5*((sigma[0]-sigma[1])**2 + (sigma[1])**2 + (sigma[0])**2)+ 6.0*(sigma[2]**2)    ) )
            #input("")
            ######################## DELETE ########################

            #ddsdde = FthTensor2Voigt(D)
            #ddsdde = FthTensor2Voigt(Calgo)
            #ddsdde = FthTensor2Voigt(Cep)
            #if self.Dim == 2:
                #idx = [0,1,3]
                #D = ddsdde[ix_(idx, idx)]
            #elif self.Dim ==3:
                #D = ddsdde
            #else:
                #assert 0, "Dimension Check"

        #GPstate.eps_e   = self.Tensor2Voigt(eps_e)
        #GPstate.eps_p   = self.Tensor2Voigt(eps_p)
        #GPstate.lamda   = lamda
        #GPstate.dstress = self.Tensor2Voigt(dsigma,'stress')
        #GPstate.stress  = self.Tensor2Voigt(sigma,'stress')
        #GPstate.deps_e  = self.Tensor2Voigt(deps_e)
        #GPstate.deps_p  = self.Tensor2Voigt(deps_p)

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
        #print(HSP)
        #exit(1)

        eps_e_tr = eps_e_init + deps # trial strain
        eps_e_tr_principal_mag, eps_e_tr_principal_vec = np.linalg.eigh(eps_e_tr)
        eps_e_tr_principal_mag, eps_e_tr_principal_vec = SortEig(eps_e_tr_principal_mag, eps_e_tr_principal_vec)
        #print("\neps_e_tr_principal_vec\n",eps_e_tr_principal_vec)

        eps_e_tr1 = eps_e_tr_principal_mag[0]
        eps_e_tr2 = eps_e_tr_principal_mag[1]
        eps_e_tr3 = eps_e_tr_principal_mag[2]

        n1 = eps_e_tr_principal_vec[:,0]
        n2 = eps_e_tr_principal_vec[:,1]
        n3 = eps_e_tr_principal_vec[:,2]

        # [2] Compute trial stress
        sigma_tr_principal_mag = np.inner(Ce_principal, eps_e_tr_principal_mag)
        #print(sigma_tr_principal_mag)
        #print("*"*50)

        #print(sigma_tr_principal_mag)
        #print(eps_e_tr_principal_mag)
        #print(eps_e_tr_principal_vec)
        ##print("-"*60)
        #input("-"*60)
        sigma_tr1 = sigma_tr_principal_mag[0] + HSP[0,0]
        sigma_tr2 = sigma_tr_principal_mag[1] + HSP[1,1]
        sigma_tr3 = sigma_tr_principal_mag[2] + HSP[2,2]



        sigma_tr = sigma_tr1*np.tensordot(n1,n1,axes=0) + sigma_tr2*np.tensordot(n2,n2,axes=0) + sigma_tr3*np.tensordot(n3,n3,axes=0)

        #print(sigma_tr1)
        #print(sigma_tr2)
        #print(sigma_tr3)
        #print("*-"*30)
        YF = self.PlasticityModel.f(sigma_tr1, sigma_tr2, sigma_tr3, lamda) # yield surface
        #print(Ce_principal)
        #print(eps_e_tr_principal_mag)
        #print(sigma_tr1, sigma_tr2, sigma_tr3, lamda)
        #print(YF)
        #print(sigma_tr1+100., sigma_tr2+100., sigma_tr3+100., lamda)
        #print(self.PlasticityModel.f(sigma_tr1+100., sigma_tr2+100., sigma_tr3+100., lamda))
        #print("*="*20)


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
            eps_e_principal_mag, eps_e_principal_vec = SortEig(eps_e_principal_mag, eps_e_principal_vec)

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

                #print(eps_e_tr_principal_mag)
                #print(sigma_tr1, sigma_tr2, sigma_tr3)
                #print(self.PlasticityModel.f(sigma_tr1, sigma_tr2, sigma_tr3, lamda))
                #print(sigma1_current, sigma2_current, sigma3_current)
                #print(self.PlasticityModel.f(sigma1_current, sigma2_current, sigma3_current, lamda_current))
                #input("=="*25)


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
            #print("\nsigma\n",sigma)
            ##print(sigma+HSP)
            #print("\nYield surface\t", self.PlasticityModel.f(sigma1, sigma2, sigma3, 0.) )
            ##print(self.PlasticityModel.f(sigma1+100., sigma2+100., sigma3+100., 0.))
            #sigma1, sigma2, sigma3
            #print("\nVonMises\t", np.sqrt(0.5*((sigma1-sigma2)**2 + (sigma2-sigma3)**2 + (sigma3-sigma1)**2))    )
            #print("*="*25)
            sigma += HSP


            deps_e = eps_e - eps_e_init
            deps_p = eps_p - eps_p_init
            dsigma = sigma - sigma_init
            ###############



            ######################################################################
            ############## Cep #############
            #deveps = eps_e - 1.0/3.0 *np.trace(eps_e)*I
            #n      = deveps / np.linalg.norm(deveps,'fro')
            #nxn  = tensor_oMult(n,n)
            #D = K*IxI+2.0*mu*(II-1.0/3.0*IxI - nxn)
            ############## Cep #############
            ######################################################################

            ######################## DELETE ########################
            #print("*="*25)

            #print("\nitr \n",ii)
            #print("\nerr \n",err)
            #print("\nHPS \n",HSP)
            #print("\nx[0], x[1], x[2] \n",x[0], x[1], x[2])
            ##print("\nk1, k2 \n",k1, k2)
            #print("\neps_e_init\n",eps_e_init)
            #print("\neps_e\n",eps_e)
            #print("\nsigma_init\n",sigma_init)
            #print("\nsigma\n",sigma)
            #print("\ndsigma\n",dsigma)
            #print("\nVM\n",np.sqrt(0.5*((sigma[0,0]-sigma[1,1])**2 + (sigma[1,1]-sigma[2,2])**2 + (sigma[2,2]-sigma[0,0])**2)+ 6.0*(sigma[0,1]**2+ sigma[1,2]**2 + sigma[2,0]**2)    ) )
            #print("\nVM_init\n",np.sqrt(0.5*((sigma_init[0,0]-sigma_init[1,1])**2 + (sigma_init[1,1]-sigma_init[2,2])**2 + (sigma_init[2,2]-sigma_init[0,0])**2)+ 6.0*(sigma_init[0,1]**2+ sigma_init[1,2]**2 + sigma_init[2,0]**2)    ) )
            #print("\nVM_tr\n",np.sqrt(0.5*((sigma_tr[0,0]-sigma_tr[1,1])**2 + (sigma_tr[1,1]-sigma_tr[2,2])**2 + (sigma_tr[2,2]-sigma_tr[0,0])**2)+ 6.0*(sigma_tr[0,1]**2+ sigma_tr[1,2]**2 + sigma_tr[2,0]**2)    ) )
            #input("")
            ######################## DELETE ########################
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

        #D = K*IxI+2.0*mu*(II-1.0/3.0*IxI)
        D = lam*IxI + 2.*mu*II


        eps_e_init  = self.Voigt2Tensor( GPstate.eps_e )
        eps_p_init  = self.Voigt2Tensor( GPstate.eps_p )
        sigma_init  = self.Voigt2Tensor( GPstate.stress ,'stress' )
        deps        = self.Voigt2Tensor( GPstate.deps  )
        lamda = GPstate.lamda   # plastic multiplier

        HSP = self.InitialaPressure()

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

                # >> Second order derivatives
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
                dx = np.linalg.solve(jac, -res) # increment of target variables

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

            ######################## DELETE ########################
            #print("*="*25)

            print("\nitr \n",ii)
            print("\nerr \n",err)
            print("\nHPS \n",HSP)
            #print("\nx[0], x[1] \n",x[0], x[1])
            #print("\nk1, k2 \n",k1, k2)
            print("\neps_e_init\n",eps_e_init)
            print("\neps_e\n",eps_e)
            print("\nsigma_init\n",sigma_init)
            print("\nsigma\n",sigma)
            print("\ndsigma\n",dsigma)
            print("\nVM\n",np.sqrt(0.5*((sigma[0,0]-sigma[1,1])**2 + (sigma[1,1]-sigma[2,2])**2 + (sigma[2,2]-sigma[0,0])**2) + 6.0*(sigma[0,1]**2) + 6.0*(sigma[1,2]**2) + 6.0*(sigma[2,0]**2)    ) )
            input("")
            ######################## DELETE ########################

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

