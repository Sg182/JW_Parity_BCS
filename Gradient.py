import numpy as np
from parameter import Nsites

def Grad_Sx(theta,i,p):
    if i== p:
        gradient_Sx = np.cos(2*theta[i])
    else:
        gradient_Sx = 0
    return gradient_Sx


def Grad_Sz(theta,i,p):
    if i== p:
        gradient_Sz = np.sin(2*theta[i])
    else:
        gradient_Sz = 0
    return gradient_Sz

def Grad_SzSz(theta,i,p,q):

    if i == p:
        gradient_SzSz = -(1/2)*(np.sin(2*theta[i]))*(np.cos(2*theta[q]))
    elif i ==q:
        gradient_SzSz = -(1/2)*(np.cos(2*theta[p]))*(np.sin(2*theta[i]))
    else:
        gradient_SzSz = 0

    return gradient_SzSz

def Grad_SzSxSz(theta,i,p,q,r):
    
    if i == p:
        gradient_SzSxSz = -(1/4)*(np.sin(2*theta[i]))*(np.sin(2*theta[q]))*(np.cos(2*theta[r]))
    elif i == q:
        gradient_SzSxSz = (1/4)*(np.cos(2*theta[p]))*(np.cos(2*theta[i]))*(np.cos(2*theta[r]))
    elif i==r:
        gradient_SzSxSz = -(1/4)*(np.cos(2*theta[p]))*(np.sin(2*theta[q]))*(np.sin(2*theta[i]))
    else:
        gradient_SzSxSz = 0
    return gradient_SzSxSz


def Gradient_energy_XXZ_theta(theta,i,Nsites,Delta):  # calculates the gradient wrt theta_i

    grad_SzSz = 0
    grad_SzSxSz = 0
    grad_Sx = 0
#---------------------------------------------------------------------------------------------------------------------#
    '''According to the functions defined above, the code in this section
      calculates the gradient wrt theta_i, and keeps the terms 
      that depends only on i  by comparing the p,q with i'''
    
    for j in range(Nsites):  
        p = (j-1)%Nsites
        q= (j+1)%Nsites
        grad_Sx += Grad_Sx(theta,i,j)*(0.5)
        grad_SzSz += Grad_SzSz(theta,i,p,q)*Delta
        grad_SzSxSz += -2*Grad_SzSxSz(theta,i,p,j,q)

    return grad_SzSxSz+grad_SzSz+grad_Sx

#---------------------------------------------------------------------------------------------------------------------#        

def Total_Grad_Sz(theta,i,Nsites):
    total_sz_grad = 0
    
    for j in range(Nsites):
         total_sz_grad += Grad_Sz(theta,i,j)

def Total_Gradient_Sz(theta,Nsites):
    total_sz = []
    for i in range(Nsites):
        total_sz.append((Total_Grad_Sz(theta,i,Nsites)))
    return np.array(total_sz)


def Total_gradient_XXZ(theta,Nsites,Delta):
    Total_gradient = []

    for i in range(Nsites):
        Total_gradient.append(Gradient_energy_XXZ_theta(theta,i,Nsites,Delta))
    
    return np.array(Total_gradient)



    