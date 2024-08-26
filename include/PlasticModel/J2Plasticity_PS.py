
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

    def f(self, sigmax, sigmay, sigmaxy, lamda):
        N       = 2.0
        n_hard  = self.n_hard
        sigma_y = self.sigma_y
        E       = self.E

        phi = (1.0/2.0)*(np.abs(sigmax-sigmay)**N + \
                         np.abs(sigmay)**N + \
                         np.abs(sigmax)**N + 6*sigmaxy**2 )
        kappa = sigma_y*(1 + E*lamda/sigma_y)**n_hard

        return phi**(1/N) - kappa

    def df(self, sigmax, sigmay, sigmaxy, lamda):
        N       = 2.0
        n_hard  = self.n_hard
        sigma_y = self.sigma_y
        E       = self.E
        x, y = sigmax, sigmay
        a    = sigmaxy

        denom = np.sqrt(2)  *  ((x-y)**2+(y)**2+(x)**2 + 6.*a**2 )**(1/2)
        #print(denom)
        #input("*="*30)
        if denom < 1e-6:
            denom = 1e-6

        dfdsigx  = (2*x-y)/ denom
        dfdsigy  = (2*y-x)/ denom
        dfdsigxy  = 6.*a/ denom

        #print("\nsigma_y\n",sigma_y)
        #print("\nn_hard\n",n_hard)
        #print("\nE\n",E)
        #print("\nlamda\n",lamda)

        dfdlamda = -n_hard*E*(1.0 + E*(lamda)/sigma_y)**(n_hard-1)
        #dfdlamda = 0.

        return dfdsigx, dfdsigy, dfdsigxy, dfdlamda

    def df2(self, sigmax, sigmay, sigmaxy):
        x, y, a= sigmax, sigmay, sigmaxy

        d2fdsig1dsig1 = 2**(1/2)/((x - y)**2 + 6*a**2 + x**2 + y**2)**(1/2) - (2**(1/2)*(4*x - 2*y)**2)/(8*((x - y)**2 + 6*a**2 + x**2 + y**2)**(3/2))
        d2fdsig1dsig2 = (2**(1/2)*(2*x - 4*y)*(4*x - 2*y))/(8*((x - y)**2 + 6*a**2 + x**2 + y**2)**(3/2)) - 2**(1/2)/(2*((x - y)**2 + 6*a**2 + x**2 + y**2)**(1/2))
        d2fdsig3dsig1 =-(3*2**(1/2)*a*(4*x - 2*y))/(2*((x - y)**2 + 6*a**2 + x**2 + y**2)**(3/2))

        d2fdsig2dsig2 = 2**(1/2)/((x - y)**2 + 6*a**2 + x**2 + y**2)**(1/2) - (2**(1/2)*(2*x - 4*y)**2)/(8*((x - y)**2 + 6*a**2 + x**2 + y**2)**(3/2))
        d2fdsig2dsig3 = (3*2**(1/2)*a*(2*x - 4*y))/(2*((x - y)**2 + 6*a**2 + x**2 + y**2)**(3/2))

        d2fdsig3dsig3 = (3*2**(1/2))/((x - y)**2 + 6*a**2 + x**2 + y**2)**(1/2) - (18*2**(1/2)*a**2)/((x - y)**2 + 6*a**2 + x**2 + y**2)**(3/2)

        return d2fdsig1dsig1, d2fdsig2dsig2, d2fdsig3dsig3, d2fdsig1dsig2, d2fdsig2dsig3, d2fdsig3dsig1
