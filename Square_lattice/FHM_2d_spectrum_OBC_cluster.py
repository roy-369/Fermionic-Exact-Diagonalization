#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  8 18:54:51 2021

@author: roy.369
"""

from quspin.operators import hamiltonian # Hamiltonians and operators
from quspin.basis import spinful_fermion_basis_general # Hilbert space spinful fermion basis
#from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np # generic math functions
from numpy import linalg as LA
#import matplotlib.pyplot as plt
import math
import time
import os
#from matplotlib import rc
from scipy import integrate
import pickle
import sys
import primme

    
def hamil_cons(t,t1,i_x,f_x,i_y,f_y,U,basis,N):
    
    hop_right=[[t,i_x[i],f_x[i]] for i in range(len(i_x))] + [[t,i_y[j],f_y[j]] for j in range(len(i_y))]
    hop_left =[[-t,i_x[i],f_x[i]] for i in range(len(i_x))] + [[-t,i_y[j],f_y[j]] for j in range(len(i_y))] 
    
    hop_right_x_polarized_1 =[[t1,i_x[i],f_x[i]] for i in range(len(i_x))]
    hop_right_x_polarized_2 =[[-t1,f_x[i],i_x[i]] for i in range(len(i_x))]
    hop_left_x_polarized_1 =[[-t1,i_x[i],f_x[i]] for i in range(len(i_x))]
    hop_left_x_polarized_2 =[[t1,f_x[i],i_x[i]] for i in range(len(i_x))]
    
    
    hop_right_y_polarized_1 =[[1j*t1,i_y[j],f_y[j]] for j in range(len(i_y))]
    hop_right_y_polarized_2 =[[1j*t1,f_y[j],i_y[j]] for j in range(len(i_y))]
    hop_left_y_polarized_1 =[[1j*t1,i_y[j],f_y[j]] for j in range(len(i_y))]
    hop_left_y_polarized_2 =[[1j*t1,f_y[j],i_y[j]] for j in range(len(i_y))]
    

     
    #hop_right = [[tf_x,i,(i+1)] for i in range(N-1)] # hopping to the right OBC
    #hop_left = [[-tr_x,i,(i+1)] for i in range(N-1)] # hopping to the left OBC
    int_list = [[U,i,i] for i in range(N)] # onsite interaction
    
    # static and dynamic lists
    static= [	
    		["+-|", hop_left], # up hop left
    		["-+|", hop_right], # up hop right
    		["|+-", hop_left], # down hop left
    		["|-+", hop_right], # down hop right
            ["-|+", hop_right_x_polarized_1],
            ["+|-", hop_right_x_polarized_2],
            ["+|-",hop_left_x_polarized_1],
            ["-|+", hop_left_x_polarized_2],
            
            ["-|+",hop_right_y_polarized_1],
            ["+|-",hop_right_y_polarized_2],
            ["+|-",hop_left_y_polarized_1],
            ["-|+",hop_left_y_polarized_2],

    		["n|n", int_list], # onsite interaction
    		]
    dynamic=[]
    H=hamiltonian(static,dynamic,dtype=np.complex64,basis=basis,check_symm=False,check_pcon = False)
    print("Hamiltonian constructed")
    return H#energies,vectors

    
def spectrum_con(t,t1,i_x,f_x,i_y,f_y,U,basis,Nx,Ny,N,Np,Raw_nos,state_no):
    
       
    H = hamil_cons(t,t1,i_x,f_x,i_y,f_y,U,basis,N)
    hh = H.tocsr()
    
    if(basis.Ns>state_no):
        en,ev = primme.eigsh(hh, state_no, tol=1e-7,which = 'SA')
    else:
        en,ev = primme.eigsh(hh, basis.Ns, tol=1e-7,which = 'SA')    
    

    filename_energy = '%s/Energies_%d_x_sites_%d_y_sites_%d_particles_%0.2f_U_%0.2f_t1.pkl' %(Raw_nos,Nx,Ny,Np,U,t1)
    data_energy = en
    with open(filename_energy, 'wb') as outfile:
       pickle.dump(data_energy, outfile, pickle.HIGHEST_PROTOCOL)
       
    filename_states = '%s/Eigen_states_%d_x_sites_%d_y_sites_%d_particles_%0.2f_U_%0.2f_t1.pkl' %(Raw_nos,Nx,Ny,Np,U,t1)
    data_eigenfunctions = ev
    with open(filename_states, 'wb') as outfile:
       pickle.dump(data_eigenfunctions, outfile, pickle.HIGHEST_PROTOCOL)

# =============================================================================
#     band_edges = []
#     for i in range(1,len(e)):
#         if(np.abs(e[i]-e[i-1])>0.2*U):
#            band_edges.append(i-1) 
#     #print(e)
# =============================================================================
    
    #print(" --- %s seconds for hamiltonian diagnolization "  %(time.time()-start_time_1))
    #return en,ev

    

def main(total,cmdargs):
    if(total!=6):
        raise ValueError('missing args')
    
        
    
    
    
    Nx = int(cmdargs[1])
    Ny = int(cmdargs[2])
    u = float(cmdargs[3])
    
    t1 = float(cmdargs[4])
    state_no = int(cmdargs[5])

    N = Nx*Ny
    t = 1
    
    s = np.arange(N) # sites [0,1,2,...,N_2d-1] in simple notation
    x = s%Nx # x positions for sites
    y = s//Nx # y positions for sites
    T_x = (x+1)%Nx + Nx*y # translation along x-direction
    T_y = x + Nx*((y+1)%Ny)
    i_x = []
    f_x = []
    i_y = []
    f_y = []

       
            
    for i in range(N):
       
       if(i%Nx!= Nx-1):
         i_x.append(i)
         f_x.append(i+1)
    
       if(int(i/Nx) != Ny-1):
         i_y.append(i)   
         f_y.append(i+Nx)  
         
  
                
# =============================================================================
    print(i_x,"i_x")
    print(f_x,"f_x")
    print(i_y,"i_y")
    print(f_y,"f_y")

#     
# =============================================================================
    filling_fraction = 1
    method_velocity = 'Commutator'
    Mag_sec = 'None'
    mz = 0
    Np = N*filling_fraction
    
        
    Nf_list = [(Np-i,i) for i in range(Np+1)]
    
    basis = spinful_fermion_basis_general(N,Nf=Nf_list)
    
    
    print(basis.Ns, "No of states")

        
    Raw_nos = "Text_files"  
    
    if not os.path.exists(Raw_nos):
        os.makedirs(Raw_nos)
  
    spectrum_con(t,t1,i_x,f_x,i_y,f_y,u,basis,Nx,Ny,N,Np,Raw_nos,state_no)
        

if __name__ == '__main__':
    sys.argv
    total = len(sys.argv)  
    print("No of sys arguments",total)
    cmdargs = sys.argv
    main(total,cmdargs)   
