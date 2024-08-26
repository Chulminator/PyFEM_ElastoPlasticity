
import numpy as np
#import autograd.numpy as np
#from autograd import elementwise_grad as egrad

from util.coordinate_transforms import *
from util.tensor_operations import *
from util.EigCalc import *


class MyPlasticity():
    def __init__(self, Fem):
        # Material properties ------------------------------------------
        # Parameters of Yld2004-18p yield function
        self.K  = Fem.K
        self.mu = Fem.mu
        self.E  = Fem.E
        self.nu = Fem.nu
        self.lam = Fem.lamda
        self.sigma_y = 387.5    # [MPa]
        self.gamma_y = 3074.6e3 # [MPa]
        self.n_y     = 0.075
        self.a     = Fem.MatProp[2]
        self.cp_12 = Fem.MatProp[3]
        self.cp_13 = Fem.MatProp[4]
        self.cp_21 = Fem.MatProp[5]
        self.cp_23 = Fem.MatProp[6]
        self.cp_31 = Fem.MatProp[7]
        self.cp_32 = Fem.MatProp[8]
        self.cp_44 = Fem.MatProp[9]
        self.cp_55 = Fem.MatProp[10]
        self.cp_66 = Fem.MatProp[11]
        self.cpp_12 = Fem.MatProp[12]
        self.cpp_13 = Fem.MatProp[13]
        self.cpp_21 = Fem.MatProp[14]
        self.cpp_23 = Fem.MatProp[15]
        self.cpp_31 = Fem.MatProp[16]
        self.cpp_32 = Fem.MatProp[17]
        self.cpp_44 = Fem.MatProp[18]
        self.cpp_55 = Fem.MatProp[19]
        self.cpp_66 = Fem.MatProp[20]

        #self.dfdsig11_g = egrad(self.f, 0)
        #self.dfdsig22_g = egrad(self.f, 1)
        #self.dfdsig33_g = egrad(self.f, 2)
        #self.dfdsig12_g = egrad(self.f, 3)
        #self.dfdsig23_g = egrad(self.f, 4)
        #self.dfdsig31_g = egrad(self.f, 5)

        #exit(1)

    def get_transform_matrix(self):
        cp_12 = self.cp_12
        cp_13 = self.cp_13
        cp_21 = self.cp_21
        cp_23 = self.cp_23
        cp_31 = self.cp_31
        cp_32 = self.cp_32
        cp_44 = self.cp_44
        cp_55 = self.cp_55
        cp_66 = self.cp_66
        cpp_12 = self.cpp_12
        cpp_13 = self.cpp_13
        cpp_21 = self.cpp_21
        cpp_23 = self.cpp_23
        cpp_31 = self.cpp_31
        cpp_32 = self.cpp_32
        cpp_44 = self.cpp_44
        cpp_55 = self.cpp_55
        cpp_66 = self.cpp_66

        Lp = 1. / 3. * np.array([
        [cp_12 + cp_13, -2. * cp_12 + cp_13, cp_12 - 2. * cp_13, 0., 0., 0.],
        [-2. * cp_21 + cp_23, cp_21 + cp_23, cp_21 - 2. * cp_23, 0., 0., 0.],
        [-2. * cp_31 + cp_32, cp_31 - 2. * cp_32, cp_31 + cp_32, 0., 0., 0.],
        [0., 0., 0., 3. * cp_44, 0., 0.],
        [0., 0., 0., 0., 3. * cp_55, 0.],
        [0., 0., 0., 0., 0., 3. * cp_66]])

        Lpp = 1. / 3. * np.array([
        [cpp_12 + cpp_13, -2. * cpp_12 + cpp_13, cpp_12 - 2. * cpp_13, 0., 0., 0.],
        [-2. * cpp_21 + cpp_23, cpp_21 + cpp_23, cpp_21 - 2. * cpp_23, 0., 0., 0.],
        [-2. * cpp_31 + cpp_32, cpp_31 - 2. * cpp_32, cpp_31 + cpp_32, 0., 0., 0.],
        [0., 0., 0., 3. * cpp_44, 0., 0.],
        [0., 0., 0., 0., 3. * cpp_55, 0.],
        [0., 0., 0., 0., 0., 3. * cpp_66]])

        #self.Lp = Lp
        #self.Lpp = Lpp

        return Lp, Lpp

    def get_rotation_tensor(self):

        alpha = 0.0
        beta  = 'None'
        gamma = 'None'

        if alpha != 'None':
          omega = np.deg2rad(0.)
          eta = np.deg2rad(90. - alpha)
          psi = np.deg2rad(0.)
        elif beta != 'None':
          omega = np.deg2rad(-90.)
          eta = np.deg2rad(90. - beta)
          psi = np.deg2rad(0.)
        elif gamma != 'None':
          omega = np.deg2rad(-gamma)
          eta = np.deg2rad(0.)
          psi = np.deg2rad(0.)

        D = np.array([[np.cos(omega), np.sin(omega), 0.],
                      [-np.sin(omega), np.cos(omega), 0.],
                      [0., 0., 1.]])
        C = np.array([[1., 0., 0.],
                      [0., np.cos(eta), np.sin(eta)],
                      [0., -np.sin(eta), np.cos(eta)]])
        B = np.array([[np.cos(psi), np.sin(psi), 0.],
                      [-np.sin(psi), np.cos(psi), 0.],
                      [0., 0., 1.]])
        A = np.dot(B, np.dot(C, D))

        return A

    def f(self, sig11, sig22, sig33, sig12, sig23, sig31, lamda):

        a       = self.a
        Rot     = self.get_rotation_tensor()
        Lp, Lpp = self.get_transform_matrix()

        sigma = np.array([[sig11, sig12, sig31],
                          [sig12, sig22, sig23],
                          [sig31, sig23, sig33]])

        sig_hat     = np.dot(Rot.T, np.dot(sigma, Rot))
        sig_hat_vec = np.array([sig_hat[0,0], sig_hat[1,1], sig_hat[2,2], sig_hat[0,1], sig_hat[1,2], sig_hat[2,0]])
        #sig_hat_vec = np.array([sigma[0,0], sigma[1,1], sigma[2,2], sigma[0,1], sigma[1,2], sigma[2,0]])

        sp =  np.dot(Lp,  sig_hat_vec)
        spp = np.dot(Lpp, sig_hat_vec)

        # Calculate principal values of sp and spp
        sp = np.array([[sp[0], sp[3], sp[5]],
                        [sp[3], sp[1], sp[4]],
                        [sp[5], sp[4], sp[2]]])
        spp = np.array([[spp[0], spp[3], spp[5]],
                        [spp[3], spp[1], spp[4]],
                        [spp[5], spp[4], spp[2]]])

        Hp_1 = (sp[0, 0] + sp[1, 1] + sp[2, 2]) / 3.
        Hp_2 = (sp[1, 2]**2. + sp[2, 0]**2. + sp[0, 1]**2. \
            - sp[1, 1] * sp[2, 2] - sp[2, 2] * sp[0, 0]   \
            - sp[0, 0] * sp[1, 1]) / 3.
        Hp_3 = (2. * sp[1, 2] * sp[2, 0] * sp[0, 1]                     \
            + sp[0, 0] * sp[1, 1] * sp[2, 2] - sp[0, 0] * sp[1, 2]**2. \
            - sp[1, 1] * sp[2, 0]**2. - sp[2, 2] * sp[0, 1]**2.) / 2.
        Hpp_1 = (spp[0, 0] + spp[1, 1] + spp[2, 2]) / 3.
        Hpp_2 = (spp[1, 2]**2. + spp[2, 0]**2. + spp[0, 1]**2. \
              - spp[1, 1] * spp[2, 2] - spp[2, 2] * spp[0, 0]  \
              - spp[0, 0] * spp[1, 1]) / 3.
        Hpp_3 = (2. * spp[1, 2] * spp[2, 0] * spp[0, 1]                       \
              + spp[0, 0] * spp[1, 1] * spp[2, 2] - spp[0, 0] * spp[1, 2]**2. \
              - spp[1, 1] * spp[2, 0]**2. - spp[2, 2] * spp[0, 1]**2.) / 2.

        p_p = Hp_1**2. + Hp_2
        q_p = (2. * Hp_1**3. + 3. * Hp_1 * Hp_2 + 2. * Hp_3) / 2.

        if np.fabs(q_p) > np.fabs(p_p**(3. / 2.)) or np.fabs(p_p) < 1e-6:
          q_p = np.sign(q_p)
          p_p = 1.

        xi_p = np.arccos(q_p / p_p**(3. / 2.))


        p_pp = Hpp_1**2. + Hpp_2
        q_pp = (2. * Hpp_1**3. + 3. * Hpp_1 * Hpp_2 + 2. * Hpp_3) / 2.

        if np.fabs(q_pp) > np.fabs(p_pp**(3. / 2.)) or np.fabs(p_pp) < 1e-6:
          q_pp = np.sign(q_pp)
          p_pp = 1.
        xi_pp = np.arccos(q_pp / p_pp**(3. / 2.))

        Sp_1 = 2. * np.sqrt(Hp_1**2. + Hp_2) * np.cos(xi_p / 3.) + Hp_1
        Sp_2 = 2. * np.sqrt(Hp_1**2. + Hp_2) * np.cos((xi_p + 4. * np.pi) / 3.) + Hp_1
        Sp_3 = 2. * np.sqrt(Hp_1**2. + Hp_2) * np.cos((xi_p + 2. * np.pi) / 3.) + Hp_1

        Spp_1 = 2. * np.sqrt(Hpp_1**2. + Hpp_2) * np.cos(xi_pp / 3.) + Hpp_1
        Spp_2 = 2. * np.sqrt(Hpp_1**2. + Hpp_2) * np.cos((xi_pp + 4. * np.pi) / 3.) + Hpp_1
        Spp_3 = 2. * np.sqrt(Hpp_1**2. + Hpp_2) * np.cos((xi_pp + 2. * np.pi) / 3.) + Hpp_1

        # define equivalent stress function
        phi = (1./4. * (np.abs(Sp_1 - Spp_1)**a + np.abs(Sp_1 - Spp_2)**a \
                      + np.abs(Sp_1 - Spp_3)**a + np.abs(Sp_2 - Spp_1)**a \
                      + np.abs(Sp_2 - Spp_2)**a + np.abs(Sp_2 - Spp_3)**a \
                      + np.abs(Sp_3 - Spp_1)**a + np.abs(Sp_3 - Spp_2)**a \
                      + np.abs(Sp_3 - Spp_3)**a) )**(1/a)

        kappa = self.sigma_y * (1. + self.gamma_y*lamda/self.sigma_y)**self.n_y
        #print(self.sigma_y)
        #print(self.gamma_y)
        #print(lamda)
        #print(sigma)
        #print("*"*10)

        return phi - kappa



    def df(self, sigma1, sigma2, sigma3, lamda):

        pp, rho, theta = convert_123_to_prt(sigma1, sigma2, sigma3)


        dfdprt[1] = 1.0
        conv_mat = convert_dfdprt_to_dfd123(p, rho, theta)
        dfdsigi = np.linalg.inv(conv_mat) @ dfdprt

        dfdsig1 = dfdsigi[0]
        dfdsig2 = dfdsigi[1]
        dfdsig3 = dfdsigi[2]
        dfdlamda = -self.n_y*self.gamma_y*(1.0 + self.gamma_y*(lamda)/self.sigma_y)**(self.n_y-1)

        return dfdsig1, dfdsig2, dfdsig3, dfdlamda

    def df2(self, sigma1, sigma2, sigma3):
        x, y, z= sigma1, sigma2, sigma3
        denom = np.sqrt(2)*((x-y)**2+(y-z)**2+(z-x)**2)**(3/2)
        if denom < 1e-6:
            denom = 1e-6
        d2fdsig1dsig1 = 3*(y - z)**2 / denom
        d2fdsig2dsig2 = 3*(z - x)**2 / denom
        d2fdsig3dsig3 = 3*(x - y)**2 / denom
        d2fdsig1dsig2 = 3*(z-x)*(y-z) / denom
        d2fdsig2dsig3 = 3*(x-y)*(z-x) / denom
        d2fdsig3dsig1 = 3*(y-x)*(z-y) / denom

        return d2fdsig1dsig1, d2fdsig2dsig2, d2fdsig3dsig3, d2fdsig1dsig2, d2fdsig2dsig3, d2fdsig3dsig1



    def df2(self, sigma1, sigma2, sigma3, lamda = 0.):

        dist = 1e-5

        dfdsig1, dfdsig2, dfdsig3, _ \
          = self.df(sigma1, sigma2, sigma3, lamda)

        dfdsig1_dist1, dfdsig2_dist1, dfdsig3_dist1, _ \
          = self.df(sigma1+dist, sigma2, sigma3, lamda)

        dfdsig1_dist2, dfdsig2_dist2, dfdsig3_dist2, _ \
          = self.df(sigma1, sigma2+dist, sigma3, lamda)

        dfdsig1_dist3, dfdsig2_dist3, dfdsig3_dist3, _ \
          = self.df(sigma1, sigma2, sigma3+dist, lamda)

        d2fdsig1dsig1 = (dfdsig1_dist1 - dfdsig1) / dist
        d2fdsig2dsig2 = (dfdsig2_dist2 - dfdsig2) / dist
        d2fdsig3dsig3 = (dfdsig3_dist3 - dfdsig3) / dist
        d2fdsig1dsig2 = (dfdsig1_dist2 - dfdsig1) / dist
        d2fdsig2dsig3 = (dfdsig2_dist3 - dfdsig2) / dist
        d2fdsig3dsig1 = (dfdsig3_dist1 - dfdsig3) / dist

        return d2fdsig1dsig1, d2fdsig2dsig2, d2fdsig3dsig3, d2fdsig1dsig2, d2fdsig2dsig3, d2fdsig3dsig1
