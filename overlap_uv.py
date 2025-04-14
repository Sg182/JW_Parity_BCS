import numpy as np
def bcs_overlap(theta,Nsites):   # computes overlap between two BCS states
    overlap = 1
    
    for k in range(Nsites):
        overlap *= ((np.cos(theta[k]))**(2) + (np.sin(theta[k]))**(2))
    return overlap


def S_x(theta,N,p):
    cos_part_p = (np.cos(theta[p]))
    sin_part_p = (np.sin(theta[p]))
    prefactor =  (1/2)*(2*cos_part_p*sin_part_p)
    prefactor_Sx = prefactor*bcs_overlap(theta,N)
    return prefactor_Sx
def Sz(theta,Nsites,i):     # calculates Sz for site i
    prefactor_Sz = (1/2)*((-1*(np.cos(theta[i]))**(2) + (np.sin(theta[i]))**(2)))/\
        ((np.cos(theta[i]))**(2) + (np.sin(theta[i]))**(2))
    return bcs_overlap(theta,Nsites)*prefactor_Sz


def S_zS_z(theta,p,q):
    cos_part_p = (np.cos(theta[p]))**2
    sin_part_p = (np.sin(theta[p]))**2
    cos_part_q =  (np.cos(theta[q]))**2
    sin_part_q = (np.sin(theta[q]))**2  
    numerator = ((-1*cos_part_p + sin_part_p)*(-1*cos_part_q + sin_part_q))
    denominator = ((cos_part_p + sin_part_p)*(cos_part_q+sin_part_q))
    prefactor = (1/4)*((numerator)/(denominator))
    prefactor_SzS_z = prefactor
    return prefactor_SzS_z


def SzSxSz(theta,N,p,i,r):
    cos_part_p = (np.cos(theta[p]))**2
    sin_part_p = (np.sin(theta[p]))**2
    cos_part_i = (np.cos(theta[i]))
    sin_part_i = (np.sin(theta[i])) 
    cos_part_r =  (np.cos(theta[r]))**2
    sin_part_r = (np.sin(theta[r]))**2
    numerator = (-cos_part_p + sin_part_p)*(2*cos_part_i*sin_part_i)*(-cos_part_r+sin_part_r)
    denominator = (cos_part_p + sin_part_p)*(cos_part_i**2+sin_part_i**2)*(cos_part_r+sin_part_r)    #Three body term
    prefactor =  (1/8)*(numerator/denominator)
    prefactor_SzSxSz = prefactor*bcs_overlap(theta,N)
    return prefactor_SzSxSz

def XXZ_1D_overlap(theta,Nsites,Delta):  # function to calculate XXZ_Energy in 1D
    sum_energy = 0
    for i in range(Nsites):
        p = (i-1)%(Nsites)   # setting the PBC
        q = (i+1)%(Nsites)   # setting the PBC
        

        sum_Sx = (1/2)*S_x(theta,Nsites,i)
        sum_SzSxSz = -2*SzSxSz(theta,Nsites,p,i,q)
        sum_SzSz = Delta*S_zS_z(theta,p,q)
        sum_energy = sum_energy + sum_Sx + sum_SzSz + sum_SzSxSz
        #print("The each component are:", sum_Sx,sum_SzSxSz,sum_SzSz)
    return sum_energy
    
