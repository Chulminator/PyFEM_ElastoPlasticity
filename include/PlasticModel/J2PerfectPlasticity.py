
import numpy as np

class MyPlasticity():
    def __init__(self, Fem):
        self.K  = Fem.K
        self.mu = Fem.mu
        self.E  = Fem.E
        self.nu = Fem.nu
        self.twoD = Fem.twoD
        self.sigma_y = Fem.MatProp[2]

    def f(self, sigma1, sigma2, sigma3, lamda):
        N       = 2.0
        sigma_y = self.sigma_y
        E       = self.E

        phi = (1.0/2.0)*(np.abs(sigma1-sigma2)**N + \
                         np.abs(sigma2-sigma3)**N + \
                         np.abs(sigma3-sigma1)**N)
        kappa = sigma_y

        return phi**(1/N) - kappa

    def df(self, sigma1, sigma2, sigma3, lamda):
        N       = 2.0
        sigma_y = self.sigma_y
        E       = self.E
        x, y, z= sigma1, sigma2, sigma3
        denom = np.sqrt(2)*((x-y)**2+(y-z)**2+(z-x)**2)**(1/2)
        #print(denom)
        #input("*="*30)
        if denom < 1e-6:
            denom = 1e-6

        dfdsig1  = (2*x-y-z)/ denom
        dfdsig2  = (2*y-z-x)/ denom
        dfdsig3  = (2*z-x-y)/ denom

        #print("\nsigma_y\n",sigma_y)
        #print("\nn_hard\n",n_hard)
        #print("\nE\n",E)
        #print("\nlamda\n",lamda)

        #dfdlamda = -n_hard*E*(1.0 + E*(lamda)/sigma_y)**(n_hard-1)
        dfdlamda = 0.

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
