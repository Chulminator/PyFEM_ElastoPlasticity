
from PlasticModel.Models_lv.NN_separate import NeuralNetwork_f
import joblib
import torch

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

        self.f_INPUT_scaler  = joblib.load("/root/Program/PyFEM_ElastoPlasticity/include/PlasticModel/Models_lv/Yld2004-18p_f_INPUT_scaler.pkl")
        self.f_OUTPUT_scaler = joblib.load("/root/Program/PyFEM_ElastoPlasticity/include/PlasticModel/Models_lv/Yld2004-18p_f_OUTPUT_scaler.pkl")
        self.model_f = torch.load("/root/Program/PyFEM_ElastoPlasticity/include/PlasticModel/Models_lv/Yld2004-18p_f_NN.pth", map_location=torch.device("cpu")).to("cpu")

        #self.dfdsig_INPUT_scaler  = joblib.load("/root/Program/PyFEM_ElastoPlasticity/include/PlasticModel/Models_lv/Yld2004-18p_dfdsig_INPUT_scaler.pkl")
        #self.dfdsig_OUTPUT_scaler = joblib.load("/root/Program/PyFEM_ElastoPlasticity/include/PlasticModel/Models_lv/Yld2004-18p_dfdsig_OUTPUT_scaler.pkl")
        #self.model_dfdsig = torch.load("/root/Program/PyFEM_ElastoPlasticity/include/PlasticModel/Models_lv/Yld2004-18p_dfdsig_NN.pth", map_location=torch.device("cpu")).to("cpu")

    def f(self, sig11, sig22, sig33, sig12, sig23, sig31, lamda):

        sigma = np.array([[sig11, sig12, sig31],
                          [sig12, sig22, sig23],
                          [sig31, sig23, sig33]])

        D, V = CalcEig(sigma)
        #D, V = np.linalg.eigh(sigma)
        D = np.diag(D)
        V = ReconMat(V)

        p, rho, theta = convert_123_to_prt(D[0], D[1], D[2])
        alp, bet, gam = rotation_angles(V)

        if (alp > np.pi or alp < -np.pi) or (bet < 0. or bet > np.pi /2. ) or (bet < 0. or bet > np.pi ):
            print(alp, bet ,gam)
            exit(1)

        RT = np.array([rho, theta, alp, bet, gam, lamda]).reshape(1,6)
        RT = self.f_INPUT_scaler.transform(RT)
        RT = torch.tensor(RT, dtype=torch.float, requires_grad = False)
        f = self.model_f(RT)
        f = self.f_OUTPUT_scaler.inverse_transform(f.detach().numpy())
        f = f[0]

        return f

    #def df(self, sig11, sig22, sig33, sig12, sig23, sig31, lamda):
        #b1 = 22.7324112957445
        #b2 = 23.5365796733886
        #b3 = 1.27708914659306e-05
        #b4 = 61.7878384488134
        #b5 = 336.393606057586
        #b6 = 0.0586474286739119
        #Y0 = 0.00689476 * 54000.0


        #dist = 0.01

        #f =           self.f(sig11, sig22, sig33, sig12, sig23, sig31, lamda)
        #f_sig11dist = self.f(sig11+dist, sig22, sig33, sig12, sig23, sig31, lamda)
        #f_sig22dist = self.f(sig11, sig22+dist, sig33, sig12, sig23, sig31, lamda)
        #f_sig33dist = self.f(sig11, sig22, sig33+dist, sig12, sig23, sig31, lamda)

        #f_sig12dist = self.f(sig11, sig22, sig33, sig12+dist, sig23, sig31, lamda)
        #f_sig23dist = self.f(sig11, sig22, sig33, sig12, sig23+dist, sig31, lamda)
        #f_sig31dist = self.f(sig11, sig22, sig33, sig12, sig23, sig31+dist, lamda)


        #dfdsig11 = (f_sig11dist-f)/dist
        #dfdsig22 = (f_sig22dist-f)/dist
        #dfdsig33 = (f_sig33dist-f)/dist
        #dfdsig12 = (f_sig12dist-f)/dist
        #dfdsig23 = (f_sig23dist-f)/dist
        #dfdsig31 = (f_sig31dist-f)/dist

        #dfdlamda = - (b1*b2/(b2*lamda+b3) + b4*b5/(b5*lamda+b6))

        #return dfdsig11, dfdsig22, dfdsig33, dfdsig12, dfdsig23, dfdsig31, dfdlamda



    def df(self, sig11, sig22, sig33, sig12, sig23, sig31, lamda):

        dfdprt = np.zeros(3)

        sigma = np.array([[sig11, sig12, sig31],
                          [sig12, sig22, sig23],
                          [sig31, sig23, sig33]])

        D, V = CalcEig(sigma)
        D = np.diag(D)
        #D, V = np.linalg.eigh(sigma)
        #print("="*40)
        #print(sigma)
        #print(D)
        #print(V)
        V = ReconMat(V)
        p, rho, theta = convert_123_to_prt(D[0], D[1], D[2])
        alp, bet ,gam = rotation_angles(V)

        if (alp > np.pi or alp < -np.pi) or (bet < 0. or bet > np.pi /2. ) or (bet < 0. or bet > np.pi ):
            print(alp, bet ,gam)
            exit(1)

        #RT = np.array([rho, theta, alp, bet, gam]).reshape(1,5)
        #RT = self.dfdsig_INPUT_scaler.transform(RT)
        #RT = torch.tensor(RT, dtype=torch.float)

        #df = self.model_dfdsig(RT)
        #df = self.dfdsig_OUTPUT_scaler.inverse_transform(df.detach().numpy())
        #df = df[0]

        #print( df)
        #exit(1)
        #dfdprt[1] = df[0]
        dfdprt[1] = 1.
        #dfdprt[2] = df[1]
        conv_mat = convert_dfdprt_to_dfd123(p, rho, theta)
        dfdsigi = np.linalg.inv(conv_mat) @ dfdprt

        tmp1 = np.array([[1., 0., 0.],
                        [0., 0., 0.],
                        [0., 0., 0.]])

        tmp2 = np.array([[0., 0., 0.],
                        [0., 1., 0.],
                        [0., 0., 0.]])

        tmp3 = np.array([[0., 0., 0.],
                        [0., 0., 0.],
                        [0., 0., 1.]])

        tmp4 = np.array([[0., 1., 0.],
                        [1., 0., 0.],
                        [0., 0., 0.]])

        tmp5 = np.array([[0., 0., 0.],
                        [0., 0., 1.],
                        [0., 1., 0.]])

        tmp6 = np.array([[0., 0., 1.],
                        [0., 0., 0.],
                        [1., 0., 0.]])

        dsigdsig11 = np.dot(V.T, np.dot(tmp1, V))
        dsigdsig22 = np.dot(V.T, np.dot(tmp2, V))
        dsigdsig33 = np.dot(V.T, np.dot(tmp3, V))
        dsigdsig12 = np.dot(V.T, np.dot(tmp4, V))
        dsigdsig23 = np.dot(V.T, np.dot(tmp5, V))
        dsigdsig31 = np.dot(V.T, np.dot(tmp6, V))

        df_dsig_11 = dfdsigi[0] * dsigdsig11[0,0] + dfdsigi[1] * dsigdsig11[1,1] + dfdsigi[2] * dsigdsig11[2,2]
        df_dsig_22 = dfdsigi[0] * dsigdsig22[0,0] + dfdsigi[1] * dsigdsig22[1,1] + dfdsigi[2] * dsigdsig22[2,2]
        df_dsig_33 = dfdsigi[0] * dsigdsig33[0,0] + dfdsigi[1] * dsigdsig33[1,1] + dfdsigi[2] * dsigdsig33[2,2]
        df_dsig_12 = dfdsigi[0] * dsigdsig12[0,0] + dfdsigi[1] * dsigdsig12[1,1] + dfdsigi[2] * dsigdsig12[2,2]
        df_dsig_23 = dfdsigi[0] * dsigdsig23[0,0] + dfdsigi[1] * dsigdsig23[1,1] + dfdsigi[2] * dsigdsig23[2,2]
        df_dsig_31 = dfdsigi[0] * dsigdsig31[0,0] + dfdsigi[1] * dsigdsig31[1,1] + dfdsigi[2] * dsigdsig31[2,2]

        b1 = 22.7324112957445
        b2 = 23.5365796733886
        b3 = 1.27708914659306e-05
        b4 = 61.7878384488134
        b5 = 336.393606057586
        b6 = 0.0586474286739119
        Y0 = 0.00689476 * 54000.0

        dfdlamda = - (b1*b2/(b2*lamda+b3) + b4*b5/(b5*lamda+b6))

        return df_dsig_11, df_dsig_22, df_dsig_33, df_dsig_12, df_dsig_23, df_dsig_31, dfdlamda


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


