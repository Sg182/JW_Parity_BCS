import numpy as np
def bcs_overlap(theta,Nsites):   # computes overlap between two BCS states
    overlap = 1
    
    for k in range(Nsites):
        overlap *= ((np.cos(theta[k]))**(2) + (np.sin(theta[k]))**(2))
    return overlap


def S_x(theta,N,p):
    cos_part_p = (np.cos(theta[p]))
    sin_part_p = (np.sin(theta[p]))
    prefactor =  (cos_part_p*sin_part_p)*np.cos(phi[p])
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
    prefactor =  (1/8)*(np.cos(2*theta[p]))*(np.sin(2*theta[i]))*(np.cos(2*theta[r]))*(np.cos(phi[i]))

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
    
'''from neighbor import *
import numpy as np'''

'''class BCSHamiltonian:
    def __init__(self, theta,phi, Nx, Ny):
        self.theta = theta
        self.phi = phi
        self.Nx = Nx
        self.Ny = Ny
        self.Nsites = Nx * Ny

    def bcs_overlap(self):
        overlap = 1
        for k in range(self.Nsites):
            overlap *= (np.cos(self.theta[k])**2 + np.sin(self.theta[k])**2)
        return overlap

    def Sz(self, i):
     
        cos_sq = np.cos(self.theta[i])**2
        sin_sq = np.sin(self.theta[i])**2
        numerator = -cos_sq + sin_sq
        denominator = cos_sq + sin_sq
        prefactor_Sz = 0.5 * numerator / denominator
        return self.bcs_overlap() * prefactor_Sz

    def S_x(self,i):
    cos_part_p = (np.cos(self.theta[i]))
    sin_part_p = (np.sin(self.theta[i]))
    prefactor =  (cos_part_p*sin_part_p)*np.cos(self.phi[i])
    prefactor_Sx = prefactor*bcs_overlap()
    return prefactor_Sx

    def Splus_Sminus(self, p, q):
        return np.cos(self.theta[p]) * np.sin(self.theta[p]) * \
               np.cos(self.theta[q]) * np.sin(self.theta[q])*(np.exp(1j*(self.phi[q]-self.phi[p]) ))

    def S_zS_z(self, p, q):
        cp, sp = np.cos(self.theta[p])**2, np.sin(self.theta[p])**2
        cq, sq = np.cos(self.theta[q])**2, np.sin(self.theta[q])**2
        numerator = (-cp + sp) * (-cq + sq)
        denominator = (cp + sp) * (cq + sq)
        return 0.25 * numerator / denominator

    def SzSxSz(theta,N,p,i,r):
    cos_part_p = (np.cos(theta[p]))**2
    sin_part_p = (np.sin(theta[p]))**2
    cos_part_i = (np.cos(theta[i]))
    sin_part_i = (np.sin(theta[i])) 
    cos_part_r =  (np.cos(theta[r]))**2
    sin_part_r = (np.sin(theta[r]))**2
    numerator = (-cos_part_p + sin_part_p)*(2*cos_part_i*sin_part_i)*(-cos_part_r+sin_part_r)
    denominator = (cos_part_p + sin_part_p)*(cos_part_i**2+sin_part_i**2)*(cos_part_r+sin_part_r)    #Three body term
    prefactor =  (1/8)*(np.cos(2*theta[p]))*(np.sin(2*theta[i]))*(np.cos(2*theta[r]))*(np.cos(phi[i]))

    prefactor_SzSxSz = prefactor*bcs_overlap(theta,N)
    return prefactor_SzSxSz

    def XXZ_overlap(self, Delta):   # add another argument 'periodic=True' and divide the body of the function with periodic/non periodic 
        sum_energy = 0
        for i in range(1, self.Nsites + 1):
            x, y = inverse_mapping(i, self.Ny)
            neighbors = neighbor_square(x, y, self.Nx, self.Ny)
            for j in neighbors:
                if i >= j:
                    continue
                p, q = i - 1, j - 1
                sum_energy += self.Splus_Sminus(p, q) + Delta * self.S_zS_z(p, q)
        return sum_energy

    def J1J2_2D_overlap(self, J):
        total_overlap_J1 = 0
        total_overlap_J2 = 0

        # First nearest neighbors
        for i in range(1, self.Nsites + 1):
            x, y = inverse_mapping(i, self.Ny)
            neighbors = neighbor_square(x, y, self.Nx, self.Ny)
            for j in neighbors:
                if i >= j:
                    continue
                p, q = i - 1, j - 1
                total_overlap_J1 += self.Splus_Sminus(p, q) + self.S_zS_z(p, q)

        # Second nearest neighbors
        for i in range(1, self.Nsites + 1):
            x, y = inverse_mapping(i, self.Ny)
            neighbors = second_neighbor(x, y, self.Nx, self.Ny)
            for j in neighbors:
                if i >= j:
                    continue
                p, q = i - 1, j - 1
                total_overlap_J2 += J * (self.Splus_Sminus(p, q) + self.S_zS_z(p, q))

        return total_overlap_J1 + total_overlap_J2'''