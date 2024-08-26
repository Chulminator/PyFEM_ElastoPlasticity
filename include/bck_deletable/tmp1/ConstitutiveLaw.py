
    # strain
def epsilon(u):

    return sym(grad(u))

def epsilon_v(u):

    strain = epsilon(u)

    return tr(strain) 
    

    # effective stress
def sigma_e(u):

    return lamda*tr(epsilon(u))*Identity(len(u)) + 2.0*mu*epsilon(u)

    # total stress
def sigma(u, p):

    stress_e = sigma_e(u)
    stress   = stress_e - p*Identity(2)

    return stress

    # Darcy's velocity
def w_Darcy(p):

    w_vel = (-k_w/mu_w)*grad(p)

    return w_vel
