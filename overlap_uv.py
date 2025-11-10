import numpy as np
def bcs_overlap(theta,Nsites):   # computes overlap between two BCS states
    overlap = 1
    
    for k in range(Nsites):
        overlap *= ((np.cos(theta[k]))**(2) + (np.sin(theta[k]))**(2))
    return overlap


def S_x(theta,phi,N,p):
    cos_part_p = (np.cos(theta[p]))
    sin_part_p = (np.sin(theta[p]))
    prefactor =  (cos_part_p*sin_part_p)*np.cos(phi[p])
    prefactor_Sx = prefactor#*bcs_overlap(theta,N)
    return 0.5*(np.sin(2*theta[p]))*np.cos(phi[p])

def Sz(theta,N,j):
    """
    Safe Sz factor at site j in spin normalization.
    If j < 0, return 1 (identity).
    """
    if j < 0:
        return 1.0
    return 0.5 * -np.cos(2*theta[j])

def S_xS_x(theta,phi,p,q):
    return 0.25*(np.sin(2*theta[p]))*(np.sin(2*theta[q]))*(np.cos(phi[p]))*(np.cos(phi[q]))


def S_zS_yS_yS_z(theta,phi,p,q,r,s):
    SySy = (np.sin(2*theta[q]))*(np.sin(2*theta[r]))*(np.sin(phi[q]))*(np.sin(phi[r]))
    SzSz = (np.cos(2*theta[p]))*(np.cos(2*theta[s]))
    return 0.0625*SySy*SzSz

def S_yS_yS_z(theta,phi,p,q,r):
    SySy = 0.25*(np.sin(2*theta[p]))*(np.sin(2*theta[q]))*(np.sin(phi[p]))*(np.sin(phi[q]))
    Sz = -0.5*(np.cos(2*theta[r]))
    return SySy*Sz

def S_zS_zS_zS_z(theta,phi,p,q,r,s):
    return 0.0625*(np.cos(2*theta[p]))*(np.cos(2*theta[q]))*(np.cos(2*theta[r]))*(np.cos(2*theta[s]))

def S_zS_zS_z(theta,phi,p,q,r):
    return 0.125*(np.cos(2*theta[p]))*(np.cos(2*theta[q]))*(np.cos(2*theta[r]))


#def Sz(theta,Nsites,i):     # calculates Sz for site i
 #   prefactor_Sz = (1/2)*((-1*(np.cos(theta[i]))**(2) + (np.sin(theta[i]))**(2)))
         
    return bcs_overlap(theta,Nsites)*prefactor_Sz

def S_xS_z(theta,phi,p,q):
    return -0.25*(np.sin(2*theta[p]))*(np.cos(2*theta[q]))*(np.cos(phi[p]))

def S_zS_z(theta,p,q):
    cos_part_p = (np.cos(theta[p]))**2
    sin_part_p = (np.sin(theta[p]))**2
    cos_part_q =  (np.cos(theta[q]))**2
    sin_part_q = (np.sin(theta[q]))**2  
    numerator = ((-1*cos_part_p + sin_part_p)*(-1*cos_part_q + sin_part_q))
    denominator = ((cos_part_p + sin_part_p)*(cos_part_q+sin_part_q))
    #prefactor = (0.25)*((numerator)*1)
    prefactor_SzS_z = 0.25*(np.cos(2*theta[p]))*(np.cos(2*theta[q]))
    return prefactor_SzS_z
    
'''def S_zS_z(theta, p, q):
    return 0.25 * np.cos(2*theta[p]) * np.cos(2*theta[q])'''
def SzSxSz(theta,phi,N,p,i,r):
    cos_part_p = (np.cos(theta[p]))**2
    sin_part_p = (np.sin(theta[p]))**2
    cos_part_i = (np.cos(theta[i]))
    sin_part_i = (np.sin(theta[i])) 
    cos_part_r =  (np.cos(theta[r]))**2
    sin_part_r = (np.sin(theta[r]))**2
    numerator = (-cos_part_p + sin_part_p)*(2*cos_part_i*sin_part_i)*(-cos_part_r+sin_part_r)
    denominator = (cos_part_p + sin_part_p)*(cos_part_i**2+sin_part_i**2)*(cos_part_r+sin_part_r)    #Three body term
    prefactor =  (0.125)*(np.cos(2*theta[p]))*(np.sin(2*theta[i]))*(np.cos(2*theta[r]))*(np.cos(phi[i]))
    prefactor_SzSxSz = prefactor*bcs_overlap(theta,N)
    return prefactor_SzSxSz



'''def SzSxSz(theta, phi,N, p, q, r):
    return 0.25 * np.cos(2*theta[p]) * np.cos(theta[q]) * np.sin(theta[q]) * np.cos(phi[q]) * np.cos(2*theta[r])'''

##This only valid for OBC for now, so the if and else statement doesn't matter here 



def XXZ_1D_overlap(theta, phi, N, Delta, periodic=False):
    E = 0.0
    # onsite term at all sites
    for i in range(N-1):
        E += 0.5*S_x(theta, phi,N, i)

    if not periodic:
        for i in range(N-1):
            p = (i - 1)
            r = (i + 1)
            if p == -1:
                E += -S_xS_z(theta, phi, i, r)
                E +=  0.5*Delta * Sz(theta,N, r)
            else:
                E += -2.0 * SzSxSz(theta, phi,N, p, i, r)
                E +=  Delta * S_zS_z(theta, p, r)

    else:
        for i in range(N-1):
            p = (i - 1)
            r = (i + 1)
            if p == -1:
                E += -S_xS_z(theta, phi, i, r)
                E +=  0.5*Delta * Sz(theta,N, r)
            else:
                E += -2.0 * SzSxSz(theta, phi,N, p, i, r)
                E +=  Delta * S_zS_z(theta, p, r)

        #adding the extra PBC terms
        X_1_X_M = (2**(-(N-1)))*np.prod(np.sin(2.0 * theta[:-1]) * np.cos(phi[:-1])) 
        Z_1_Z_M = 2**(-3)*np.cos(2*theta[0])* np.cos(2*theta[-2]) * np.cos(2*theta[-1])

        prod_X = np.prod(np.sin(2.0*theta[1:-2]) * np.cos(phi[1:-2])) if N > 3 else 1.0
        Y1   = np.sin(2.0*theta[0]) * np.sin(phi[0])
        YM1  = np.sin(2.0*theta[-2])* np.sin(phi[-2])  # site M-1
        ZM   = np.cos(2.0*theta[-1]) 
        Y_1_Y_M = 2**(-N)*prod_X*Y1*ZM*YM1

        E += 2**(N-3)*X_1_X_M + 2*Z_1_Z_M + 2**(N-2)*Y_1_Y_M

        

    return E

def XXZ_2D_overlap(theta,phi,N,Delta,periodic=False):
    E = 0.0
    #if not periodic:
        

    
def J1J2_1D_overlap(theta,phi,N,J2):
    E_J_1 = 0.0
    E_J_2 = 0.0

    '''for i in range(N-1):
        E_J_1 += 0.5*S_x(theta, phi,N, i)

    if True:
        for i in range(N-1):
            p = (i - 1)
            r = (i + 1)
            if p == -1:
                E_J_1 += -S_xS_z(theta, phi, i, r)
                E_J_1 +=  0.5* Sz(theta, r)
            else:
                E_J_1 += -2.0 * SzSxSz(theta, phi,N, p, i, r)
                E_J_1 +=  S_zS_z(theta, p, r)'''
    E_J_1 = XXZ_1D_overlap(theta,phi,N,Delta=1,periodic=False)

    # For J2 part

    for i in range(N-2):
        j = i+1
        E_J_2 += S_xS_x(theta,phi,i,j)
    
    if True:
        for i in range(N-2):
            p = i -1
            r = i+1
            s = i+2

            if p == -1:
                E_J_2 += 4*0.5*S_yS_yS_z(theta,phi,i,r,s)
                E_J_2 += 4*0.5*S_zS_zS_z(theta,phi,i,r,s)
            else:
                E_J_2 += 4*S_zS_yS_yS_z(theta,phi,p,i,r,s)
                E_J_2 += 4*S_zS_zS_zS_z(theta,phi,p,i,r,s)

    return E_J_1 + E_J_2*J2


    

'''def XXZ_1D_overlap(theta, phi, N, Delta, periodic=True):
    E = 0.0
    # onsite term at all sites
    for i in range(N-1):
        E += 0.5*S_x(theta, phi,N, i)

    if periodic:
        for i in range(N-1):
            p = (i - 1) % N
            r = (i + 1) % N
            E += -2.0 * SzSxSz(theta, phi,N, p, i, r)
            E +=  Delta * S_zS_z(theta, p, r)
    else:
        for i in range(1, N-1):
            p, r = i - 1, i + 1
            E += -2.0 * SzSxSz(theta, phi,N, p, i, r)
            E +=  Delta * S_zS_z(theta, p, r)
    return E'''

 