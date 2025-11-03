import numpy as np
from scipy.optimize import minimize
from overlap_pair import *    # <-- import here
from parameter import *
from SFS_Ham_2 import *
#from overlap_uv import Sz
import sys


import warnings
warnings.filterwarnings("ignore")
   # for example

best_obj = np.inf
best_theta = None
ham = SFS_1D_XXZ(Nsites,Delta)
for i in range(100):
    theta0 = np.random.uniform(-np.pi, np.pi, Nsites)

    E_before = energy_theta(
        theta0,
        ham.h001, ham.h100, ham.h010, ham.h011,
        ham.h110, ham.h020, ham.h021, ham.h120
    )

    result = minimize(
        energy_theta,
        theta0,
        args=(ham.h001, ham.h100, ham.h010, ham.h011, ham.h110, ham.h020, ham.h021, ham.h120),
        method='BFGS',
        options={'gtol': 1e-12, 'maxiter': 2000}
    )

    E_final = result.fun

    #print(f"run {i+1}: Energy_prev = {E_before:.12f} , Energy_final = {E_final:.12f}")
    print(f"run {i+1}: Energy_final = {E_final:.12f}")
    

    if E_final < best_obj:
        best_obj = E_final
        best_theta = result.x.copy()

sys.stdout = sys.__stdout__
print("Enuc =", ham.h000 )
Total_energy = best_obj + ham.h000
print(f"\nThe Energy for {Delta}: {(best_obj + ham.h000).real:.14f}")

#print(f"Global Sz: {sum(Sz(best_theta, Nsites, i) for i in range(Nsites))}")
