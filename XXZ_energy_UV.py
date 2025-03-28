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


# ----------------------------------------------------------------------------------------------------------------#
# Prepare the result data to write to file
log_file = 'log.txt'
sys.stdout = open(log_file,'w')
# ----------------------------------------------------------------------------------------------------------------#



#-----------------------------------------------------------------------------------------------------------------#
 
# Define constraint dictionary
constraint = {'type': 'eq', 'fun': lambda theta :Sz_sum(theta,Nsites)}  #Using Lambda function since Sz_sum has two arguments

# Perform minimization with Sz = 0 constraint
result = minimize(XXZ_1D_energy, theta_initial, args=(Nsites,Delta,), method='trust-constr',jac=Total_gradient_XXZ, constraints=constraint,\
                  options={'verbose':2})
theta_optimized = result.x
final_energy = result.fun

#--------------------------------------------------------------------------------------------------------------------#

 
with open(theta_file,'w') as file:     #overwriting the etas in the .txt file
    for theta in theta_optimized:
        file.write(f"{theta}\n")

#sys.stdout.close()
'''for theta in theta_optimized:
    print(f"{np.cos(theta)}\t{np.sin(theta)}\n")'''

with open('energy_XXZ_1D_12_JW_parity.txt',"a") as file:      # writing the energy to a text file
    file.write(f"{Delta}   {final_energy:.12f}\n")

sys.stdout = sys.__stdout__

if result.success:
    print('optimization successful!')
else:
    print("Warning! optimization failed: ", result.message)
 
#print(bcs_overlap(theta_optimized,Nsites))
print(f"The Energy for {Delta}:  {final_energy:.12f}")
print(f"The global Sz is: {Sz_sum(theta_optimized,Nsites)}")
