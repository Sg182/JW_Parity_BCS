from overlap import *
from parameter import *
import numpy as np
from scipy.optimize import minimize
import os


def XXZ_1D_energy(eta,Nsites,Delta):  # function to calculate XXZ_Energy in 1D
    return (XXZ_1D_overlap(eta,eta,Nsites,Delta)/bcs_overlap(eta,eta,Nsites))

     
def Sz_sum(eta,Nsites):
    Sz_global= 0
    for i in range(1,Nsites+1):  # This calculates Sz for each site 
    #print(Sz(eta_optimized,eta_optimized,i)/bcs_overlap(eta_optimized,eta_optimized))  
        Sz_global += Sz(eta,eta,Nsites,i)/bcs_overlap(eta,eta,Nsites)
    return Sz_global


#--------------------------------------------------------------------------------------#
'''This section loads the etas of the previous calculation as an initial guess for the current calculation.
If the file doesn't exist, it creates a folder'''


eta_file = 'eta_opt.txt'

if os.path.exists(eta_file):
    eta_initial = np.loadtxt(eta_file)

else:
    eta_initial = np.random.uniform(-1,1,Nsites)

#---------------------------------------------------------------------------------------#
#print(eta_initial)

# Define constraint dictionary
constraint = {'type': 'eq', 'fun': lambda eta :Sz_sum(eta,Nsites)}  #Using Lambda function since Sz_sum has two arguments

# Perform minimization with Sz = 0 constraint
result = minimize(XXZ_1D_energy, eta_initial, args=(Nsites,Delta,), method='SLSQP', constraints=constraint)

eta_optimized = result.x
final_energy = result.fun

print(f"The Energy for {Delta}:  {final_energy}")
#print(f"The global Sz is: {Sz_sum(eta_optimized,Nsites)}")

with open(eta_file,'w') as file:     #overwriting the etas in the .txt file
    for eta in eta_optimized:
        file.write(f"{eta}\n")

'''with open('energy_XXZ_1D_12_JW_parity.txt',"a") as file:      # writing the energy to a text file
    file.write(f"{Delta}   {final_energy:.12f}\n")'''




    







 