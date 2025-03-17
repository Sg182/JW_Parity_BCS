from overlap_uv import *
from parameter import *
import numpy as np
from Gradient import Total_gradient_XXZ, Total_Gradient_Sz
import os
import sys
import cyipopt

def XXZ_1D_energy(theta, Nsites, Delta):  
    return (XXZ_1D_overlap(theta, Nsites, Delta) / bcs_overlap(theta, Nsites))

def Sz_sum(theta, Nsites):
    Sz_global = 0
    for i in range(Nsites):
        Sz_global += Sz(theta, Nsites, i)
    return Sz_global / bcs_overlap(theta, Nsites)

# Load initial guess or generate random values for theta
theta_file = 'theta_opt.txt'
if os.path.exists(theta_file):
    theta_initial = np.loadtxt(theta_file)
else:
    theta_initial = np.random.uniform(-np.pi, np.pi, Nsites)

# Prepare the result data to write to file
log_file = 'log.txt'
sys.stdout = open(log_file, 'w')

class XXZ_1D_NLP(cyipopt.Problem):
    def __init__(self, Nsites, Delta):
        self.Nsites = Nsites
        self.Delta = Delta

        lb = np.full(Nsites, -np.pi)  # Lower bound for theta
        ub = np.full(Nsites, np.pi)   # Upper bound for theta

        super(XXZ_1D_NLP, self).__init__(
            n=Nsites,                  
            m=1,                       
            lb=lb,                     
            ub=ub,                     
            cl=np.zeros(1),            
            cu=np.zeros(1)             
        )

    def objective(self, theta):
        return XXZ_1D_energy(theta, self.Nsites, self.Delta)

    def gradient(self, theta):
        return Total_gradient_XXZ(theta, self.Nsites, self.Delta)

    def constraints(self, theta):
        return np.array([Sz_sum(theta, self.Nsites)])

    def jacobian(self, theta):
        return Total_Gradient_Sz(theta, self.Nsites)

    def intermediate(self, alg_mod, iter_count, obj_value, inf_pr, inf_du, mu,
                     d_norm, regularization_size, alpha_du, alpha_pr, ls_trials):
        print(f"Iteration {iter_count}: Energy = {obj_value}, Constraint violation = {inf_pr}")
        return True

# Create the problem instance
nlp = XXZ_1D_NLP(Nsites, Delta)

# Set solver options if needed
nlp.add_option('max_iter', 1000)  # Example: maximum number of iterations
nlp.add_option('tol', 1e-8)       # Example: tolerance for stopping criteria
nlp.add_option('mu_strategy', 'adaptive')  # Adaptive mu strategy
nlp.add_option('print_level', 5)   # Print level for logging

# Solve the problem
theta_optimized = (nlp.solve(theta_initial))[0]

final_energy = XXZ_1D_energy(theta_optimized, Nsites, Delta)

# Save the optimized theta values
with open(theta_file, 'w') as file:
    for theta in theta_optimized:
        file.write(f"{theta}\n")

# Save the final energy
with open('energy_XXZ_1D_12_JW_parity.txt', "a") as file:
    file.write(f"{Delta}   {final_energy:.11f}\n")

sys.stdout = sys.__stdout__

# Print results
print(f"The Energy for {Delta}:  {final_energy:.12f}")
print(f"The global Sz is: {Sz_sum(theta_optimized, Nsites)}")
