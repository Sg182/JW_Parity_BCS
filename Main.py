from overlap_uv import *
from parameter import *
import numpy as np
from Gradient import Total_gradient_XXZ 
from scipy.optimize import minimize
import os
import sys




def XXZ_1D_energy(theta,Nsites,Delta):  # function to calculate XXZ_Energy in 1D
    return (XXZ_1D_overlap(theta,Nsites,Delta)/bcs_overlap(theta,Nsites))

    
def Sz_sum(theta,Nsites):
    Sz_global= 0
    for i in range(Nsites):  # This calculates Sz for each site 
    #print(Sz(eta_optimized,eta_optimized,i)/bcs_overlap(eta_optimized,eta_optimized))  
        Sz_global += Sz(theta,Nsites,i)
    return Sz_global/bcs_overlap(theta,Nsites)


#---------------------------------------------------------------------------------------------------------#
'''This section loads the etas of the previous calculation as an initial guess for the current calculation.
If the file doesn't exist, it creates a folder'''


theta_file = 'theta_opt.txt'

if os.path.exists(theta_file):
    theta_initial = np.loadtxt(theta_file)

else:
     theta_initial = np.random.uniform(-np.pi,np.pi,Nsites)

best_obj = np.inf
best_theta = None
theta0 = theta_initial
for i in range(1):
    #theta0 = np.random.uniform(-np.pi,np.pi, Nsites)
    #theta0 = np.array([np.pi/4,-np.pi/4,np.pi/4,-np.pi/4,np.pi/4,-np.pi/4,
                      # np.pi/4,-np.pi/4,np.pi/4,-np.pi/4,np.pi/4,-np.pi/4])
    
# ----------------------------------------------------------------------------------------------------------------#
# Prepare the result data to write to file
    '''log_file = 'log.txt'
    sys.stdout = open(log_file,'w')'''
# ----------------------------------------------------------------------------------------------------------------#



#-----------------------------------------------------------------------------------------------------------------#
 
# Define constraint dictionary
    constraint = {'type': 'eq', 'fun': lambda theta :Sz_sum(theta,Nsites)}  #Using Lambda function since Sz_sum has two arguments

# Perform minimization with Sz = 0 constraint
    #result = minimize(XXZ_1D_energy, theta0, args=(Nsites,Delta,), method='trust-constr',jac=Total_gradient_XXZ, constraints=constraint,\
                  #options={ 'maxiter':2000})
    result = minimize(XXZ_1D_energy, theta0, args=(Nsites,Delta,), constraints=constraint,\
                  options={'xtol':1e-12,'maxiter':2000})
    #theta_optimized = result.x
    final_energy = result.fun
    print(f"run {i+1}: Energy : {final_energy}")

    if final_energy < best_obj:
        best_obj = final_energy
        best_theta = result.x

#--------------------------------------------------------------------------------------------------------------------#

 
with open(theta_file,'w') as file:     #overwriting the etas in the .txt file
    for theta in best_theta:
        file.write(f"{theta}\n")

#sys.stdout.close()
'''for theta in theta_optimized:
    print(f"{np.cos(theta)}\t{np.sin(theta)}\n")'''

with open('energy_XXZ_1D_12_JW_parity.txt',"a") as file:      # writing the energy to a text file
    file.write(f"{Delta}   {best_obj:.12f}\n")

sys.stdout = sys.__stdout__

if result.success:
    print('optimization successful!')
else:
    print("Warning! optimization failed: ", result.message)
 
#print(bcs_overlap(theta_optimized,Nsites))
print(f"The Energy for {Delta}:  {best_obj:.12f}")
print(f"The global Sz is: {Sz_sum(best_theta,Nsites)}")



'''from neighbor import *
import numpy as np
# Only for Real u,v'''



 
 

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

    def Splus_Sminus(self, p, q):
        return np.cos(self.theta[p]) * np.sin(self.theta[p]) * \
               np.cos(self.theta[q]) * np.sin(self.theta[q])*(np.exp(1j*(self.phi[q]-self.phi[p]) ))

    def S_zS_z(self, p, q):
        cp, sp = np.cos(self.theta[p])**2, np.sin(self.theta[p])**2
        cq, sq = np.cos(self.theta[q])**2, np.sin(self.theta[q])**2
        numerator = (-cp + sp) * (-cq + sq)
        denominator = (cp + sp) * (cq + sq)
        return 0.25 * numerator / denominator

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