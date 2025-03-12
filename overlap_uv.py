import numpy as np
def bcs_overlap(theta,Nsites):   # computes overlap between two BCS states
    overlap = 1
    for k in range(Nsites):
        overlap *= (np.cos(theta[k]).conjugate()*np.cos(theta[k]) + np.sin(theta[k]).conjugate()*np.sin(theta[k]))
    return overlap


def S_x(theta,N,i):
    prefactor = (1/2)*((np.cos(theta[i]).conjugate()*np.sin(theta[i]))+(np.cos(theta[i])*np.sin(theta[i]).conjugate()))/\
        ((np.cos(theta[i]).conjugate()*np.cos(theta[i])+np.sin(theta[i]).conjugate()*np.sin(theta[i])))
    prefactor_Sx = prefactor*bcs_overlap(theta,N)
    return prefactor_Sx

def Sz(theta,Nsites,i):     # calculates Sz for site i
    prefactor_Sz = (1/2)*((-1*np.cos(theta[i-1]).conjugate()*np.cos(theta[i-1]) + np.cos(theta[i-1]).conjugate()*np.cos(theta[i-1])))/\
        ((np.cos(theta[i-1]).conjugate()*np.cos(theta[i-1])+np.sin(theta[i-1]).conjugate()*np.sin(theta[i-1])))
    
    return bcs_overlap(theta,Nsites)*prefactor_Sz


def S_zS_z(theta,N,p,q):
    
    prefactor = (1/4)*(((-1*np.cos(theta[p]).conjugate()*np.cos(theta[p]) + np.cos(theta[p]).conjugate()*np.cos(theta[p])))*\
                       ((-1*np.cos(theta[q]).conjugate()*np.cos(theta[q]) + np.cos(theta[q]).conjugate()*np.cos(theta[q]))))/\
                       (((np.cos(theta[p]).conjugate()*np.cos(theta[p])+np.sin(theta[p]).conjugate()*np.sin(theta[p])))*\
                        ((np.cos(theta[q]).conjugate()*np.cos(theta[q])+np.sin(theta[q]).conjugate()*np.sin(theta[q]))))
    prefactor_SzS_z = prefactor*bcs_overlap(theta,N)
    return prefactor_SzS_z

def SzSxSz(theta,N,p,i,q):    #Three body term
    prefactor = (1/8)*(((-1*np.cos(theta[p]).conjugate()*np.cos(theta[p]) + np.cos(theta[p]).conjugate()*np.cos(theta[p])))*\
                       ((np.cos(theta[i]).conjugate()*np.sin(theta[i]))+(np.cos(theta[i])*np.sin(theta[i]).conjugate()))*\
                        ((-1*np.cos(theta[q]).conjugate()*np.cos(theta[q]) + np.cos(theta[q]).conjugate()*np.cos(theta[q]))))/\
                       (((np.cos(theta[p]).conjugate()*np.cos(theta[p])+np.sin(theta[p]).conjugate()*np.sin(theta[p])))*\
                        ((np.cos(theta[i]).conjugate()*np.cos(theta[i])+np.sin(theta[i]).conjugate()*np.sin(theta[i])))*\
                        ((np.cos(theta[q]).conjugate()*np.cos(theta[q])+np.sin(theta[q]).conjugate()*np.sin(theta[q]))))
    prefactor_SzSxSz = prefactor*bcs_overlap(theta,N)
    return prefactor_SzSxSz

def XXZ_1D_overlap(theta,Nsites,Delta):  # function to calculate XXZ_Energy in 1D
    sum_energy = 0
    for i in range(Nsites):
        p = (i-1)%(Nsites)   # setting the PBC
        q = (i+1)%(Nsites)   # setting the PBC
        

        sum_Sx = (0.5)*S_x(theta,Nsites,i)
        sum_SzSxSz = -2*SzSxSz(theta,Nsites,p,i,q)
        sum_szSz = Delta*S_zS_z(theta,Nsites,p,q)
        sum_energy = sum_energy + sum_Sx + sum_szSz + sum_SzSxSz

    return sum_energy
    
    
    