
import numpy as np

class MyPlasticity():
    def __init__(self, Fem):
        self.K  = Fem.K
        self.mu = Fem.mu
        self.E  = Fem.E
        self.nu = Fem.nu
        self.twoD = Fem.twoD

    def f(self, sigma1, sigma2, sigma3, lamda):
        return -1.

    def df(self, sigma1, sigma2, sigma3, lamda):
        dfdsig1  = 0.
        dfdsig2  = 0.
        dfdsig3  = 0.

        dfdlamda = 0.

        return dfdsig1, dfdsig2, dfdsig3, dfdlamda

    def df2(self, sigma1, sigma2, sigma3):
        d2fdsig1dsig1 = 0.
        d2fdsig2dsig2 = 0.
        d2fdsig3dsig3 = 0.
        d2fdsig1dsig2 = 0.
        d2fdsig2dsig3 = 0.
        d2fdsig3dsig1 = 0.
        return d2fdsig1dsig1, d2fdsig2dsig2, d2fdsig3dsig3, d2fdsig1dsig2, d2fdsig2dsig3, d2fdsig3dsig1
