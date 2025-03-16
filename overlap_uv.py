import numpy as np
def bcs_overlap(theta,Nsites):   # computes overlap between two BCS states
    overlap = 1
    for k in range(Nsites):
        overlap *= ((np.cos(theta[k]))**(2) + (np.sin(theta[k]))**(2))
    return overlap


def S_x(theta,N,i):
    prefactor = (1/2)*((np.cos(theta[i])*np.sin(theta[i]))+(np.cos(theta[i])*np.sin(theta[i])))/\
        ((np.cos(theta[i])*np.cos(theta[i])+np.sin(theta[i])*np.sin(theta[i])))
    prefactor_Sx = prefactor*bcs_overlap(theta,N)
    return prefactor_Sx

def Sz(theta,Nsites,i):     # calculates Sz for site i
    prefactor_Sz = (1/2)*((-1*(np.cos(theta[i]))**(2) + (np.cos(theta[i]))**(2)))/\
        ((np.cos(theta[i]))**(2) + (np.cos(theta[i]))**(2))
    return bcs_overlap(theta,Nsites)*prefactor_Sz


def S_zS_z(theta,N,p,q):
    
    prefactor = (1/4)*(((-1*(np.cos(theta[p]))**(2) + (np.sin(theta[p]))**(2)))*\
                       ((-1*(np.cos(theta[q]))**(2) + np.sin(theta[q])**(2))))/\
                       ((((np.cos(theta[p]))**2+(np.sin(theta[p]))**(2)))*\
                        (((np.cos(theta[q]))**(2)+(np.sin(theta[q]))**(2))))
    prefactor_SzS_z = prefactor*bcs_overlap(theta,N)
    return prefactor_SzS_z

def SzSxSz(theta,N,p,i,q):    #Three body term
    prefactor = (1/8)*(((-1*(np.cos(theta[p]))**(2) + (np.cos(theta[p]))**(2)))*\
                       ((np.cos(theta[i])*np.sin(theta[i]))+(np.cos(theta[i])*np.sin(theta[i])))*\
                        ((-1*(np.cos(theta[q]))**(2) + (np.cos(theta[q]))**(2))))/\
                       ((((np.cos(theta[p]))**(2) + (np.sin(theta[p]))**(2)))*\
                        (((np.cos(theta[i]))**(2) + (np.sin(theta[i]))**(2)))*\
                        ((((np.cos(theta[i]))**(2) + (np.sin(theta[i]))**(2)))))
    prefactor_SzSxSz = prefactor*bcs_overlap(theta,N)
    return prefactor_SzSxSz

def XXZ_1D_overlap(theta,Nsites,Delta):  # function to calculate XXZ_Energy in 1D
    sum_energy = 0
    for i in range(Nsites):
        p = (i-1)%(Nsites)   # setting the PBC
        q = (i+1)%(Nsites)   # setting the PBC
        

        sum_Sx = (1/2)*S_x(theta,Nsites,i)
        sum_SzSxSz = -2*SzSxSz(theta,Nsites,p,i,q)
        sum_SzSz = Delta*S_zS_z(theta,Nsites,p,q)
        sum_energy = sum_energy + sum_Sx + sum_SzSz + sum_SzSxSz

    return sum_energy
    
