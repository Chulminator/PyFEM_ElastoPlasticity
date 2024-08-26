
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


    def df(self, sig11, sig22, sig33, sig12, sig23, sig31, lamda):
        dist = 0.01
        ################################################################################
        f           = self.f(sig11,  sig22,  sig33,  sig12,  sig23,  sig31, lamda)
        f_sig11dist = self.f(sig11+dist, sig22, sig33, sig12, sig23, sig31, lamda)
        f_sig22dist = self.f(sig11, sig22+dist, sig33, sig12, sig23, sig31, lamda)
        f_sig33dist = self.f(sig11, sig22, sig33+dist, sig12, sig23, sig31, lamda)

        f_sig12dist = self.f(sig11, sig22, sig33, sig12+dist, sig23, sig31, lamda)
        f_sig23dist = self.f(sig11, sig22, sig33, sig12, sig23+dist, sig31, lamda)
        f_sig31dist = self.f(sig11, sig22, sig33, sig12, sig23, sig31+dist, lamda)

        dfdsig11 = (f_sig11dist-f)/dist
        dfdsig22 = (f_sig22dist-f)/dist
        dfdsig33 = (f_sig33dist-f)/dist
        dfdsig12 = (f_sig12dist-f)/dist
        dfdsig23 = (f_sig23dist-f)/dist
        dfdsig31 = (f_sig31dist-f)/dist
        ################################################################################

        #dfdsig11_tmp = self.dfdsig11_g(sig11, sig22, sig33, sig12, sig23, sig31, lamda)
        #dfdsig22_tmp = self.dfdsig22_g(sig11, sig22, sig33, sig12, sig23, sig31, lamda)
        #dfdsig33_tmp = self.dfdsig33_g(sig11, sig22, sig33, sig12, sig23, sig31, lamda)
        #dfdsig12_tmp = self.dfdsig12_g(sig11, sig22, sig33, sig12, sig23, sig31, lamda)
        #dfdsig23_tmp = self.dfdsig23_g(sig11, sig22, sig33, sig12, sig23, sig31, lamda)
        #dfdsig31_tmp = self.dfdsig31_g(sig11, sig22, sig33, sig12, sig23, sig31, lamda)

        dfdlamda = -self.n_y*self.gamma_y*(1.0 + self.gamma_y*(lamda)/self.sigma_y)**(self.n_y-1)

        return dfdsig11, dfdsig22, dfdsig33, dfdsig12, dfdsig23, dfdsig31, dfdlamda




    #def df(self, sig11, sig22, sig33, sig12, sig23, sig31, lamda):

        #dfdprt = np.zeros(3)

        #sigma = np.array([[sig11, sig12, sig31],
                          #[sig12, sig22, sig23],
                          #[sig31, sig23, sig33]])

        ##D, V = CalcEig(sigma)
        #D, V = np.linalg.eig(sigma)
        ##D = np.diag(D)
        #V = ReconMat(V)
        #p, rho, theta = convert_123_to_prt(D[0], D[1], D[2])

        #dfdprt[1] = 1.0
        #conv_mat = convert_dfdprt_to_dfd123(p, rho, theta)
        #dfdsigi = np.linalg.inv(conv_mat) @ dfdprt

        #tmp1 = np.array([[1., 0., 0.],
                        #[0., 0., 0.],
                        #[0., 0., 0.]])

        #tmp2 = np.array([[0., 0., 0.],
                        #[0., 1., 0.],
                        #[0., 0., 0.]])

        #tmp3 = np.array([[0., 0., 0.],
                        #[0., 0., 0.],
                        #[0., 0., 1.]])

        #tmp4 = np.array([[0., 1., 0.],
                        #[1., 0., 0.],
                        #[0., 0., 0.]])

        #tmp5 = np.array([[0., 0., 0.],
                        #[0., 0., 1.],
                        #[0., 1., 0.]])

        #tmp6 = np.array([[0., 0., 1.],
                        #[0., 0., 0.],
                        #[1., 0., 0.]])

        #dsigdsig11 = np.dot(V.T, np.dot(tmp1, V))
        #dsigdsig22 = np.dot(V.T, np.dot(tmp2, V))
        #dsigdsig33 = np.dot(V.T, np.dot(tmp3, V))
        #dsigdsig12 = np.dot(V.T, np.dot(tmp4, V))
        #dsigdsig23 = np.dot(V.T, np.dot(tmp5, V))
        #dsigdsig31 = np.dot(V.T, np.dot(tmp6, V))

        #df_dsig_11 = dfdsigi[0] * dsigdsig11[0,0] + dfdsigi[1] * dsigdsig11[1,1] + dfdsigi[2] * dsigdsig11[2,2]
        #df_dsig_22 = dfdsigi[0] * dsigdsig22[0,0] + dfdsigi[1] * dsigdsig22[1,1] + dfdsigi[2] * dsigdsig22[2,2]
        #df_dsig_33 = dfdsigi[0] * dsigdsig33[0,0] + dfdsigi[1] * dsigdsig33[1,1] + dfdsigi[2] * dsigdsig33[2,2]
        #df_dsig_12 = dfdsigi[0] * dsigdsig12[0,0] + dfdsigi[1] * dsigdsig12[1,1] + dfdsigi[2] * dsigdsig12[2,2]
        #df_dsig_23 = dfdsigi[0] * dsigdsig23[0,0] + dfdsigi[1] * dsigdsig23[1,1] + dfdsigi[2] * dsigdsig23[2,2]
        #df_dsig_31 = dfdsigi[0] * dsigdsig31[0,0] + dfdsigi[1] * dsigdsig31[1,1] + dfdsigi[2] * dsigdsig31[2,2]


        #dfdlamda = -self.n_y*self.gamma_y*(1.0 + self.gamma_y*(lamda)/self.sigma_y)**(self.n_y-1)
        #return df_dsig_11, df_dsig_22, df_dsig_33, df_dsig_12, df_dsig_23, df_dsig_31, dfdlamda




    def df2(self, sig11, sig22, sig33, sig12, sig23, sig31, lamda = 0.):

        dist = 1e-6

        dfdsig11, dfdsig22, dfdsig33, dfdsig12, dfdsig23, dfdsig31, _ \
          = self.df(sig11, sig22, sig33, sig12, sig23, sig31,lamda)

        dfdsig11_sig11dist, dfdsig22_sig11dist, dfdsig33_sig11dist, dfdsig12_sig11dist, dfdsig23_sig11dist, dfdsig31_sig11dist, _  \
          = self.df(sig11+dist, sig22, sig33, sig12, sig23, sig31,lamda)

        dfdsig11_sig22dist, dfdsig22_sig22dist, dfdsig33_sig22dist, dfdsig12_sig22dist, dfdsig23_sig22dist, dfdsig31_sig22dist, _ \
          = self.df(sig11, sig22+dist, sig33, sig12, sig23, sig31,lamda)

        dfdsig11_sig33dist, dfdsig22_sig33dist, dfdsig33_sig33dist, dfdsig12_sig33dist, dfdsig23_sig33dist, dfdsig31_sig33dist, _ \
          = self.df(sig11, sig22, sig33+dist, sig12, sig23, sig31,lamda)

        dfdsig11_sig12dist, dfdsig22_sig12dist, dfdsig33_sig12dist, dfdsig12_sig12dist, dfdsig23_sig12dist, dfdsig31_sig12dist, _ \
          = self.df(sig11, sig22, sig33, sig12+dist, sig23, sig31,lamda)

        dfdsig11_sig23dist, dfdsig22_sig23dist, dfdsig33_sig23dist, dfdsig12_sig23dist, dfdsig23_sig23dist, dfdsig31_sig23dist, _ \
          = self.df(sig11, sig22, sig33, sig12, sig23+dist, sig31,lamda)

        dfdsig11_sig31dist, dfdsig22_sig31dist, dfdsig33_sig31dist, dfdsig12_sig31dist, dfdsig23_sig31dist, dfdsig31_sig31dist, _ \
          = self.df(sig11, sig22, sig33, sig12, sig23, sig31+dist,lamda)

        d2fdsig11sig11 = (dfdsig11 - dfdsig11_sig11dist) / dist
        d2fdsig11sig22 = (dfdsig11 - dfdsig11_sig22dist) / dist
        d2fdsig11sig33 = (dfdsig11 - dfdsig11_sig33dist) / dist
        d2fdsig11sig12 = (dfdsig11 - dfdsig11_sig12dist) / dist
        d2fdsig11sig23 = (dfdsig11 - dfdsig11_sig23dist) / dist
        d2fdsig11sig31 = (dfdsig11 - dfdsig11_sig31dist) / dist

        d2fdsig22sig22 = (dfdsig22 - dfdsig22_sig22dist) / dist
        d2fdsig22sig33 = (dfdsig22 - dfdsig22_sig33dist) / dist
        d2fdsig22sig12 = (dfdsig22 - dfdsig22_sig12dist) / dist
        d2fdsig22sig23 = (dfdsig22 - dfdsig22_sig23dist) / dist
        d2fdsig22sig31 = (dfdsig22 - dfdsig22_sig31dist) / dist

        d2fdsig33sig33 = (dfdsig33 - dfdsig33_sig33dist) / dist
        d2fdsig33sig12 = (dfdsig33 - dfdsig33_sig12dist) / dist
        d2fdsig33sig23 = (dfdsig33 - dfdsig33_sig23dist) / dist
        d2fdsig33sig31 = (dfdsig33 - dfdsig33_sig31dist) / dist

        d2fdsig12sig12 = (dfdsig12 - dfdsig12_sig12dist) / dist
        d2fdsig12sig23 = (dfdsig12 - dfdsig12_sig23dist) / dist
        d2fdsig12sig31 = (dfdsig12 - dfdsig12_sig31dist) / dist

        d2fdsig23sig23 = (dfdsig23 - dfdsig23_sig23dist) / dist
        d2fdsig23sig31 = (dfdsig23 - dfdsig23_sig31dist) / dist

        d2fdsig31sig31 = (dfdsig31 - dfdsig31_sig31dist) / dist

        return d2fdsig11sig11, d2fdsig11sig22, d2fdsig11sig33, d2fdsig11sig12, d2fdsig11sig23, d2fdsig11sig31, \
                              d2fdsig22sig22, d2fdsig22sig33, d2fdsig22sig12, d2fdsig22sig23, d2fdsig22sig31, \
                                              d2fdsig33sig33, d2fdsig33sig12, d2fdsig33sig23, d2fdsig33sig31, \
                                                              d2fdsig12sig12, d2fdsig12sig23, d2fdsig12sig31, \
                                                                              d2fdsig23sig23, d2fdsig23sig31, \
                                                                                              d2fdsig31sig31


