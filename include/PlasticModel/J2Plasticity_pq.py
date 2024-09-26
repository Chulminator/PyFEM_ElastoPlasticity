
import numpy as np

class MyPlasticity():
    def __init__(self, Fem):
        self.K  = Fem.K
        self.mu = Fem.mu
        self.E  = Fem.E
        self.nu = Fem.nu
        self.twoD = Fem.twoD
        self.sigma_y = Fem.MatProp[2]
        self.n_hard  = Fem.MatProp[3]

    def f(self, p, q, lamda):
        N       = 2.0
        n_hard  = self.n_hard
        sigma_y = self.sigma_y
        E       = self.E

        phi = q
        kappa = sigma_y*(1 + E*lamda/sigma_y)**n_hard

        return phi - kappa

    def df(self, p, q, lamda):
        N       = 2.0
        n_hard  = self.n_hard
        sigma_y = self.sigma_y
        E       = self.E

        dfdp = 0.
        dfdq = 1.

        dfdlamda = -n_hard*E*(1.0 + E*(lamda)/sigma_y)**(n_hard-1)

        return dfdp, dfdq, dfdlamda

    def df2(self, p, q):
        d2fdpdp = 0.
        d2fdpdq = 0.
        d2fdpdq = 0.

        return d2fdpdp, d2fdpdq, d2fdqdq
