
class Myplasticity():
    def __init__(self, Fem):
        self.K  = Fem.K
        self.mu = Fem.mu
        self.E  = Fem.E
        self.nu = Fem.nu
        self.twoD = Fem.twoD
        self.sigma_y = Fem.MatProp[2]
        self.n_hard  = Fem.MatProp[3]

    def f(self, sigma1, sigma2, sigma3, lamda):
        N       = 2.0
        n_hard  = self.n_hard
        sigma_y = self.sigma_y
        E       = self.E

        phi = (1.0/2.0)*(np.abs(sigma1-sigma2)**N + \
                    np.abs(sigma2-sigma3)**N + \
                    np.abs(sigma3-sigma1)**N)
        kappa = sigma_y*(1 + E*lamda/sigma_y)**n_hard

        return phi**(1/N) - kappa

    def df(self, sigma1, sigma2, sigma3):
        x, y, z= sigma1, sigma2, sigma3
        denom = np.sqrt(2)*((x-y)**2+(y-z)**2+(z-x)**2)**(1/2)
        dfdsig1 = (2*x-y-z)/ denom
        dfdsig2 = (2*y-z-x)/ denom
        dfdsig3 = (2*z-x-y)/ denom
        return dfdsig1, dfdsig2, dfdsig3

    def df2(self, sigma1, sigma2, sigma3):
        x, y, z= sigma1, sigma2, sigma3
        denom = np.sqrt(2)*((x-y)**2+(y-z)**2+(z-x)**2)**(3/2)
        d2fdsig1dsig1 = 3*(y - z)**2 / denom
        d2fdsig2dsig2 = 3*(z - x)**2 / denom
        d2fdsig3dsig3 = 3*(x - y)**2 / denom
        d2fdsig1dsig2 = 3*(z-x)*(y-z) / denom
        d2fdsig2dsig3 = 3*(x-y)*(z-x) / denom
        d2fdsig3dsig1 = 3*(y-x)*(z-y) / denom
        return d2fdsig1dsig1, d2fdsig2dsig2, d2fdsig3dsig3, d2fdsig1dsig2, d2fdsig2dsig3, d2fdsig3dsig1

    def Voigt2Tensor2DStrain(self, voigt):
        # https://www.comsol.com/blogs/what-is-the-difference-between-plane-stress-and-plane-strain
        eps_e = np.zeros((3,3))
        if self.twoD == 'planestrain':
            eps_e[0,0] = voigt[0]
            eps_e[0,1] = voigt[2]
            eps_e[1,0] = voigt[2]
            eps_e[1,1] = voigt[1]

        elif self.twoD == 'planestress':
            nu = self.nu
            eps_e[0,0] = voigt[0]
            eps_e[0,1] = voigt[2]
            eps_e[1,0] = voigt[2]
            eps_e[1,1] = voigt[1]
            eps_e[2,2] = -nu/(1.0-nu)*(eps_e[0,0] + eps_e[1,1])
        else:
            assert 0,"Check Fem.twoD"
        return eps_e

    def Tensor2Voigt2D(self, Tensor, flag='strain'):
        # https://www.comsol.com/blogs/what-is-the-difference-between-plane-stress-and-plane-strain
        voigt = np.zeros((3))
        if flag == 'strain':
            voigt[0] = Tensor[0,0]
            voigt[1] = Tensor[1,1]
            voigt[2] = Tensor[0,1]
        elif flag == 'stress':
            voigt[0] = Tensor[0,0]
            voigt[1] = Tensor[1,1]
            voigt[2] = Tensor[0,1]*0.5
        else:
            assert(0,"check the flag in Tensor2Voigt2D")
        return voigt
