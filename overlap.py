import numpy as np
def bcs_overlap(epsilon,eta,Nsites):   # computes overlap between two BCS states
    overlap = 1
    for k in range(Nsites):
        overlap *= ( 1+ (epsilon[k].conjugate())*eta[k])
    return overlap


def S_x(epsilon,eta,N,i):
    prefactor = (0.5)*((epsilon[i].conjugate() + eta[i])/(1 + epsilon[i].conjugate()*eta[i]))
    prefactor_Sx = prefactor*bcs_overlap(epsilon,eta,N)
    return prefactor_Sx

def Sz(epsilon,eta,Nsites,i):     # calculates Sz for site i
    prefactor_Sz = (1/2)*(-1+epsilon[i-1].conjugate()*eta[i-1])/(1 + epsilon[i-1].conjugate()*eta[i-1]) 
    return bcs_overlap(epsilon,eta,Nsites)*prefactor_Sz


def S_zS_z(epsilon,eta,N,p,q):
    
    prefactor = (1/4)*((-1 + epsilon[p].conjugate()*eta[p])/(1 + epsilon[p].conjugate()*eta[p]))*\
        ((-1 + epsilon[q].conjugate()*eta[q])/(1 + epsilon[q].conjugate()*eta[q]))
    prefactor_SzS_z = prefactor*bcs_overlap(epsilon,eta,N)
    return prefactor_SzS_z

def SzSxSz(epsilon,eta,N,p,i,q):
    prefactor = (1/4)*((-1 + epsilon[p].conjugate()*eta[p])/(1 + epsilon[p].conjugate()*eta[p]))*\
    ((epsilon[i].conjugate() + eta[i])/(1 + epsilon[i].conjugate()*eta[i]))*\
    ((-1 + epsilon[q].conjugate()*eta[q])/(1 + epsilon[q].conjugate()*eta[q]))
    prefactor_SzSxSz = prefactor*bcs_overlap(epsilon,eta,N)
    return prefactor_SzSxSz

def XXZ_1D_overlap(epsilon,eta,Nsites,Delta):  # function to calculate XXZ_Energy in 1D
    sum_energy = 0
    for i in range(Nsites):
        p = (i-1)%(Nsites)   # setting the PBC
        q = (i+1)%(Nsites)   # setting the PBC
        

        sum_Sx = (0.5)*S_x(epsilon,eta,Nsites,i)
        sum_SzSxSz = SzSxSz(epsilon,eta,Nsites,p,i,q)
        sum_szSz = Delta*S_zS_z(epsilon,eta,Nsites,p,q)
        sum_energy = sum_energy + sum_Sx + sum_szSz + sum_SzSxSz

    return sum_energy
    
    
    