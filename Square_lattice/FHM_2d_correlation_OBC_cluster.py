#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  8 18:54:51 2021

@author: roy.369
"""

from quspin.operators import hamiltonian # Hamiltonians and operators
from quspin.basis import spinful_fermion_basis_general # Hilbert space spinful fermion basis
import numpy as np # generic math functions
from numpy import linalg as LA
import math
import time
import os
from scipy import integrate
import pickle
import sys


    

def sz_sz_cons(site_index_1,site_index_2,basis,N):
    sz_sz_list1 = [[0.25,site_index_1,site_index_2]]
    sz_sz_list2a = [[-0.25,site_index_1,site_index_2]]
    sz_sz_list2b = [[-0.25,site_index_2,site_index_1]]
    static = [["nn|",sz_sz_list1],
              ["|nn",sz_sz_list1],
              ["n|n",sz_sz_list2a],
              ["n|n",sz_sz_list2b]]
    dynamic = []
    Sz_corr = hamiltonian(static,dynamic,basis = basis, dtype = np.complex128,check_symm = False,check_pcon=False)
    
    return Sz_corr
    
def corr_func_cons(en,ev,Site_1,Site_2,basis,N):
   
    Sz_correlation = sz_sz_cons(Site_1,Site_2,basis,N)
    sz_correlation = Sz_correlation.toarray()
    Sz_Sz_ij = 0
    for n in range(len(en)):
        p_n = ev[:,n]
        if(n ==0):
            f_n = 1
        else:
            f_n  =0
            Sz_Sz_ij = Sz_Sz_ij+f_n*np.inner(np.conj(p_n),np.dot(sz_correlation,p_n))
    
    return Sz_Sz_ij


def correlation_grid_cons(data_dir,Nx,Ny,Np,U,t1,basis,N):
    
    filename_energy = "%s/Energies_%d_x_sites_%d_y_sites_%d_particles_%0.2f_U_%0.2f_t1.pkl" %(data_dir,Nx,Ny,Np,U,t1)
    with open(filename_energy, 'rb') as infile:
       en = pickle.load(infile)
    
    filename_states = '%s/Eigen_states_%d_x_sites_%d_y_sites_%d_particles_%0.2f_U_%0.2f_t1.pkl' %(data_dir,Nx,Ny,Np,U,t1)
    with open(filename_states, 'rb') as infile:
       ev = pickle.load(infile)
    
    Sz_Sz_correlation_grid = np.zeros((N,N),dtype = np.complex128)
    for site_1  in range(N):
        for site_2 in range(N):
            if(site_1>=site_2):
               Sz_Sz_correlation_grid[site_1][site_2] = corr_func_cons(en,ev,site_1,site_2,basis,N)

    filename_mz_mz = '%s/Sz_Sz_correlation_y_%d_x_sites_%d_y_sites_%d_particles_%0.2f_U_%0.2f_t1.pkl' %(data_dir,Nx,Ny,Np,U,t1)
    data_mz_mz = Sz_Sz_correlation_grid
    with open(filename_mz_mz, 'wb') as outfile:
       pickle.dump(data_mz_mz, outfile, pickle.HIGHEST_PROTOCOL)



 
def main(total,cmdargs):
    if(total!=6):
        raise ValueError('missing args')
    
      
    Nx_dir = cmdargs[1]
    Ny_dir = cmdargs[2]
    u_dir = cmdargs[3]
    t1_dir = cmdargs[4]
    
    Nx = int(Nx_dir)
    Ny = int(Ny_dir)
    u= float(u_dir)
    t1 = float(t1_dir)
    state_no = int(cmdargs[5])

    N = Nx*Ny
    t = 1
    s = np.arange(N) # sites [0,1,2,...,N_2d-1] in simple notation
    x = s%Nx # x positions for sites
    y = s//Nx # y positions for sites
    T_x = (x+1)%Nx + Nx*y # translation along x-direction
    T_y = x + Nx*((y+1)%Ny)

    rx = np.zeros(N)
    ry = np.zeros(N)
# =============================================================================
    
    a  = 1
    for i in range(N):
         rx[i] = a*x[i]
         ry[i] = a*y[i]

    filling_fraction = 1
    Mag_sec = 'None'
    mz = 0
    Np = N*filling_fraction
    
        
    Nf_list = [(Np-i,i) for i in range(Np+1)]
    
    basis = spinful_fermion_basis_general(N,Nf=Nf_list)
    
    
    print(basis.Ns, "No of states")

        
    data_dir = "Text_files" 
    correlation_grid_cons(data_dir,Nx,Ny,Np,u,t1,basis,N)


if __name__ == '__main__':
    sys.argv
    total = len(sys.argv)  
    print("No of sys arguments",total)
    cmdargs = sys.argv
    main(total,cmdargs)   
