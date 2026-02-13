import pickle
from overlap_uv import *
from parameter import *
import numpy as np
#from Hub2Pair import hub2pair
#from NBCSCC_X3 import nbcscc
from Gradient import Total_gradient_XXZ 
from scipy.optimize import minimize
import os
import sys




def XXZ_1D_energy(var,Nsites,Delta):
    theta = var[:Nsites]
    phi = var[Nsites:]
    numerator = XXZ_1D_overlap(theta, phi, Nsites, Delta, periodic=False)
    denominator = bcs_overlap(theta, Nsites)
    return (numerator / denominator)  # function to calculate XXZ_Energy in 1D
    #return (XXZ_1D_overlap(theta,phi,Nsites,Delta,periodic=False)/bcs_overlap(theta,Nsites))

    
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
#theta0 = np.random.uniform(-np.pi,np.pi,Nsites)
for i in range(100):
    theta0 = np.random.uniform(-np.pi,np.pi, Nsites)
    phi0 = np.random.uniform(-np.pi, np.pi, Nsites)
    var0 = np.concatenate([theta0, phi0])
    #theta0 = np.array([np.pi/4,-np.pi/4,np.pi/4,-np.pi/4,np.pi/4,-np.pi/4,
                      # np.pi/4,-np.pi/4,np.pi/4,-np.pi/4,np.pi/4,-np.pi/4])
    
# ----------------------------------------------------------------------------------------------------------------#
# Prepare the result data to write to file
    '''log_file = 'log.txt'
    sys.stdout = open(log_file,'w')'''
# ----------------------------------------------------------------------------------------------------------------#

    E_before = XXZ_1D_energy(var0,Nsites,Delta)


#-----------------------------------------------------------------------------------------------------------------#
    '''Sz_target = Nsites/2   # USE THIS CONSTRAINT AT FERROMAGNETIC REGIME
    constraint = {
    'type': 'eq',
    'fun': lambda var: Sz_sum(var[:Nsites], Nsites) - Sz_target
}'''

# Define constraint dictionary
    constraint = {'type': 'eq', 'fun': lambda theta :Sz_sum(theta,Nsites)}  #Using Lambda function since Sz_sum has two arguments

# Perform minimization with Sz = 0 constraint
     
    '''result = minimize(
    XXZ_1D_energy, var0, args=(Nsites, Delta),
    method='SLSQP',
    constraints=[constraint],
    options={'ftol': 1e-12, 'maxiter': 2000})'''


    result = minimize(XXZ_1D_energy, var0, args=(Nsites,Delta,), constraints=constraint,\
                  options={'xtol':1e-12,'maxiter':2000})
    #theta_optimized = result.x
    final_energy = result.fun
    print(f"run {i+1}: Energy : {final_energy}")

    if final_energy < best_obj:
        best_obj = final_energy
        best_theta = result.x[:Nsites]
        best_phi = result.x[Nsites:]
        best_var = np.concatenate([best_theta,best_phi])
        U = np.cos(best_theta)
        V = np.sin(best_theta)
#--------------------------------------------------------------------------------------------------------------------#

with open("UV.p","wb") as f:
    pickle.dump({"U":U,"V":V},f)
 
with open(theta_file,'w') as file:     #overwriting the etas in the .txt file
    for theta in best_theta:
        file.write(f"{theta}\n")

#sys.stdout.close()
'''for theta in theta_optimized:
    print(f"{np.cos(theta)}\t{np.sin(theta)}\n")'''

with open('energy_XXZ_1D_12_JW_parity_OBC.txt',"a") as file:      # writing the energy to a text file
    file.write(f"{Delta}   {best_obj:.12f}\n")

sys.stdout = sys.__stdout__

if result.success:
    print('optimization successful!')
else:
    print("Warning! optimization failed: ", result.message)
 
#print(bcs_overlap(theta_optimized,Nsites))
print(f"The Energy for {Delta}:  {best_obj:.12f}")
print(f"The global Sz is: {Sz_sum(best_theta,Nsites)}")



 
 
 

 