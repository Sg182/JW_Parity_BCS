def bcs_overlap(epsilon,eta):   # computes overlap between two BCS states
    overlap = 1
    for k in range(Nx*Ny):
        overlap *= ( 1+ (epsilon[k].conjugate())*eta[k])
    return overlap


def S_x (epsilon,eta,i):
    prefactor = (0.5)*((epsilon[i].conjugate() + eta[i])/(1 + epsilon[i].conjuagte()*eta[i]))
    prefactor_Sx = prefactor*bcs_overlap(epsilon,eta)
    return prefactor_Sx


def S_zS_z(epsilon,eta,i):
    prefactor = (1/4)*((-1 + epsilon[i-1].conjugate()*eta[i-1])/(1 + epsilon[i-1].conjugate()*eta[i-1]))*\
        ((-1 + epsilon[i+1].conjugate()*eta[i+1])/(1 + epsilon[i+1].conjugate()*eta[i+1]))
    prefactor_SzS_z = prefactor*bcs_overlap(epsilon,eta)
    return prefactor_SzS_z

def SzSxSz(epsilon,eta,i):
    prefactor = (1/4)*((-1 + epsilon[i-1].conjugate()*eta[i-1])/(1 + epsilon[i-1].conjugate()*eta[i-1]))*\
    ((epsilon[i].conjugate() + eta[i])/(1 + epsilon[i].conjuagte()*eta[i]))*\
    ((-1 + epsilon[i+1].conjugate()*eta[i+1])/(1 + epsilon[i+1].conjugate()*eta[i+1]))
    prefactor_SzSxSz = prefactor*bcs_overlap
    return prefactor_SzSxSz
    
    
    