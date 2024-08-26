
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

    def f(self, sigma11, sigma22, sigma33,  sigma12,  sigma23,  sigma31, lamda):
        N       = 2.0
        n_hard  = self.n_hard
        sigma_y = self.sigma_y
        E       = self.E

        phi = (1.0/2.0)*(np.abs(sigma11-sigma22)**N + \
                         np.abs(sigma22-sigma33)**N + \
                         np.abs(sigma33-sigma11)**N + \
                         6.0*(sigma12**2 + sigma23**2 + sigma31**2 ))
        kappa = sigma_y*(1 + E*lamda/sigma_y)**n_hard

        return phi**(1/N) - kappa

    def df(self, sigma11, sigma22, sigma33,  sigma12,  sigma23,  sigma31, lamda):
        N       = 2.0
        n_hard  = self.n_hard
        sigma_y = self.sigma_y
        E       = self.E
        x, y, z= sigma11, sigma22, sigma33
        a, b, c= sigma12, sigma23, sigma31

        denom = np.sqrt(2)*((x-y)**2 + (y-z)**2 + (z-x)**2 + 6.0*(a**2 + b**2 + c**2))**(1/2)

        #print(denom)
        #input("*="*30)
        if denom < 1e-6:
            denom = 1e-6

        dfdsig11  = (2.*x-y-z) / denom
        dfdsig22  = (2.*y-z-x) / denom
        dfdsig33  = (2.*z-x-y) / denom
        dfdsig12  = 3.*a*np.sqrt(2.0) / denom
        dfdsig23  = 3.*b*np.sqrt(2.0) / denom
        dfdsig31  = 3.*c*np.sqrt(2.0) / denom

        dfdlamda = -n_hard*E*(1.0 + E*(lamda)/sigma_y)**(n_hard-1)

        return dfdsig11, dfdsig22, dfdsig33, dfdsig12, dfdsig23, dfdsig31, dfdlamda

    def df2(self, sigma11, sigma22, sigma33,  sigma12,  sigma23,  sigma31):

        x, y, z= sigma11, sigma22, sigma33
        a, b, c= sigma12, sigma23, sigma31

        d2fd1111 = 1/((x - y)**2/2 + (x - z)**2/2 + (y - z)**2/2 + 3*a**2 + 3*b**2 + 3*c**2)**(1/2) - (y - 2*x + z)**2/(4*((x - y)**2/2 + (x - z)**2/2 + (y - z)**2/2 + 3*a**2 + 3*b**2 + 3*c**2)**(3/2))
        d2fd1122 =- 1/(2*((x - y)**2/2 + (x - z)**2/2 + (y - z)**2/2 + 3*a**2 + 3*b**2 + 3*c**2)**(1/2)) - ((x - 2*y + z)*(y - 2*x + z))/(4*((x - y)**2/2 + (x - z)**2/2 + (y - z)**2/2 + 3*a**2 + 3*b**2 + 3*c**2)**(3/2))
        d2fd1133 = - 1/(2*((x - y)**2/2 + (x - z)**2/2 + (y - z)**2/2 + 3*a**2 + 3*b**2 + 3*c**2)**(1/2)) - ((x + y - 2*z)*(y - 2*x + z))/(4*((x - y)**2/2 + (x - z)**2/2 + (y - z)**2/2 + 3*a**2 + 3*b**2 + 3*c**2)**(3/2))

        d2fd1112 = (3*a*(y - 2*x + z))/(2*((x - y)**2/2 + (x - z)**2/2 + (y - z)**2/2 + 3*a**2 + 3*b**2 + 3*c**2)**(3/2))
        d2fd1123 = (3*b*(y - 2*x + z))/(2*((x - y)**2/2 + (x - z)**2/2 + (y - z)**2/2 + 3*a**2 + 3*b**2 + 3*c**2)**(3/2))
        d2fd1131 = (3*c*(y - 2*x + z))/(2*((x - y)**2/2 + (x - z)**2/2 + (y - z)**2/2 + 3*a**2 + 3*b**2 + 3*c**2)**(3/2))


        d2fd2222 = 1/((x - y)**2/2 + (x - z)**2/2 + (y - z)**2/2 + 3*a**2 + 3*b**2 + 3*c**2)**(1/2) - (x - 2*y + z)**2/(4*((x - y)**2/2 + (x - z)**2/2 + (y - z)**2/2 + 3*a**2 + 3*b**2 + 3*c**2)**(3/2))
        d2fd2233 =- 1/(2*((x - y)**2/2 + (x - z)**2/2 + (y - z)**2/2 + 3*a**2 + 3*b**2 + 3*c**2)**(1/2)) - ((x + y - 2*z)*(x - 2*y + z))/(4*((x - y)**2/2 + (x - z)**2/2 + (y - z)**2/2 + 3*a**2 + 3*b**2 + 3*c**2)**(3/2))
        d2fd2212 = (3*a*(x - 2*y + z))/(2*((x - y)**2/2 + (x - z)**2/2 + (y - z)**2/2 + 3*a**2 + 3*b**2 + 3*c**2)**(3/2))
        d2fd2223 = (3*b*(x - 2*y + z))/(2*((x - y)**2/2 + (x - z)**2/2 + (y - z)**2/2 + 3*a**2 + 3*b**2 + 3*c**2)**(3/2))
        d2fd2231 = (3*c*(x - 2*y + z))/(2*((x - y)**2/2 + (x - z)**2/2 + (y - z)**2/2 + 3*a**2 + 3*b**2 + 3*c**2)**(3/2))

        d2fd3333 = 1/((x - y)**2/2 + (x - z)**2/2 + (y - z)**2/2 + 3*a**2 + 3*b**2 + 3*c**2)**(1/2) - (x + y - 2*z)**2/(4*((x - y)**2/2 + (x - z)**2/2 + (y - z)**2/2 + 3*a**2 + 3*b**2 + 3*c**2)**(3/2))
        d2fd3312 = (3*a*(x + y - 2*z))/(2*((x - y)**2/2 + (x - z)**2/2 + (y - z)**2/2 + 3*a**2 + 3*b**2 + 3*c**2)**(3/2))
        d2fd3323 = (3*b*(x + y - 2*z))/(2*((x - y)**2/2 + (x - z)**2/2 + (y - z)**2/2 + 3*a**2 + 3*b**2 + 3*c**2)**(3/2))
        d2fd3331 = (3*c*(x + y - 2*z))/(2*((x - y)**2/2 + (x - z)**2/2 + (y - z)**2/2 + 3*a**2 + 3*b**2 + 3*c**2)**(3/2))


        d2fd1212 = 3/((x - y)**2/2 + (x - z)**2/2 + (y - z)**2/2 + 3*a**2 + 3*b**2 + 3*c**2)**(1/2) - (9*a**2)/((x - y)**2/2 + (x - z)**2/2 + (y - z)**2/2 + 3*a**2 + 3*b**2 + 3*c**2)**(3/2)
        d2fd1223 = -(9*a*b)/((x - y)**2/2 + (x - z)**2/2 + (y - z)**2/2 + 3*a**2 + 3*b**2 + 3*c**2)**(3/2)
        d2fd1231 = -(9*a*c)/((x - y)**2/2 + (x - z)**2/2 + (y - z)**2/2 + 3*a**2 + 3*b**2 + 3*c**2)**(3/2)

        d2fd2323 = 3/((x - y)**2/2 + (x - z)**2/2 + (y - z)**2/2 + 3*a**2 + 3*b**2 + 3*c**2)**(1/2) - (9*b**2)/((x - y)**2/2 + (x - z)**2/2 + (y - z)**2/2 + 3*a**2 + 3*b**2 + 3*c**2)**(3/2)
        d2fd2331 = -(9*b*c)/((x - y)**2/2 + (x - z)**2/2 + (y - z)**2/2 + 3*a**2 + 3*b**2 + 3*c**2)**(3/2)

        d2fd3131 = 3/((x - y)**2/2 + (x - z)**2/2 + (y - z)**2/2 + 3*a**2 + 3*b**2 + 3*c**2)**(1/2) - (9*c**2)/((x - y)**2/2 + (x - z)**2/2 + (y - z)**2/2 + 3*a**2 + 3*b**2 + 3*c**2)**(3/2)

        return d2fd1111, d2fd1122, d2fd1133, d2fd1112, d2fd1123, d2fd1131, \
                         d2fd2222, d2fd2233, d2fd2212, d2fd2223, d2fd2231, \
                                   d2fd3333, d2fd3312, d2fd3323, d2fd3331, \
                                             d2fd1212, d2fd1223, d2fd1231, \
                                                       d2fd2323, d2fd2331, \
                                                                 d2fd3131
