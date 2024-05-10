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


    

def dipole_x_cons(rx,a,basis,N):
    dipole = [[a*rx[i],i,i] for i in range(N)]
    static = [["+-|",dipole],
              ["|+-",dipole]]
    dynamic = []
    dipole_x = hamiltonian(static,dynamic,basis = basis, dtype = np.complex128,check_symm = False,check_pcon=False)
    return dipole_x

def dipole_y_cons(ry,a,basis,N):
    dipole = [[a*ry[i],i,i] for i in range(N)]
    static = [["+-|",dipole],
              ["|+-",dipole]]
    dynamic = []
    dipole_y = hamiltonian(static,dynamic,basis = basis, dtype = np.complex128,check_symm = False,check_pcon=False)
    return dipole_y

def dipole_op_con(t,t1,a,rx,ry,U,basis,Nx,Ny,N,Np,Raw_nos,data_dir):    
    Dx = dipole_x_cons(rx,a,basis,N)
    Dy = dipole_x_cons(ry,a,basis,N)
    dx = Dx.tocsr()
    dy = Dy.tocsr()

    filename_energy = "%s/Energies_%d_x_sites_%d_y_sites_%d_particles_%0.2f_U_%0.2f_t1.pkl" %(data_dir,Nx,Ny,Np,U,t1)
    with open(filename_energy, 'rb') as infile:
       en = pickle.load(infile)

    filename_states = '%s/Eigen_states_%d_x_sites_%d_y_sites_%d_particles_%0.2f_U_%0.2f_t1.pkl' %(data_dir,Nx,Ny,Np,U,t1)
    with open(filename_states, 'rb') as infile:
       ev = pickle.load(infile)

    dipole_x_mat = np.zeros((len(en),len(en)),dtype = np.complex128)
    dipole_y_mat = np.zeros((len(en),len(en)),dtype = np.complex128)
    lambda_x_mat = np.zeros((len(en),len(en)),dtype = np.complex128)
    lambda_y_mat = np.zeros((len(en),len(en)),dtype = np.complex128)
    for n in range(len(en)):
        p_n = ev[:,n]
        bra = p_n.conj().T
        for m in range(len(en)):
            p_m = ev[:,m]
            ket = p_m.copy()
            dipole_x_mat[n][m] = bra.dot(dx.dot(ket))
            dipole_y_mat[n][m] = bra.dot(dy.dot(ket))
            lambda_x_mat[n][m] = -1j*(en[n]-en[m])*dipole_x_mat[n][m]
            lambda_y_mat[n][m] = -1j*(en[n]-en[m])*dipole_y_mat[n][m]

    filename_dipole_x = '%s/Dipole_x_%d_x_sites_%d_y_sites_%d_particles_%0.2f_U_%0.2f_t1.pkl' %(Raw_nos,Nx,Ny,Np,U,t1)
    data_dipole_x = dipole_x_mat
    with open(filename_dipole_x, 'wb') as outfile:
       pickle.dump(data_dipole_x, outfile, pickle.HIGHEST_PROTOCOL)

    filename_dipole_y = '%s/Dipole_y_%d_x_sites_%d_y_sites_%d_particles_%0.2f_U_%0.2f_t1.pkl' %(Raw_nos,Nx,Ny,Np,U,t1)
    data_dipole_y = dipole_y_mat
    with open(filename_dipole_y, 'wb') as outfile:
       pickle.dump(data_dipole_y, outfile, pickle.HIGHEST_PROTOCOL)

    filename_current_x = '%s/Lambda_x_%d_x_sites_%d_y_sites_%d_particles_%0.2f_U_%0.2f_t1.pkl' %(Raw_nos,Nx,Ny,Np,U,t1)
    data_current_x = lambda_x_mat
    with open(filename_current_x, 'wb') as outfile:
       pickle.dump(data_current_x, outfile, pickle.HIGHEST_PROTOCOL)

    filename_current_y = '%s/Lambda_y_%d_x_sites_%d_y_sites_%d_particles_%0.2f_U_%0.2f_t1.pkl' %(Raw_nos,Nx,Ny,Np,U,t1)
    data_current_y = lambda_y_mat
    with open(filename_current_y, 'wb') as outfile:
       pickle.dump(data_current_y, outfile, pickle.HIGHEST_PROTOCOL)

    lambda_xx_mat = np.zeros((len(en),len(en)),dtype = np.complex128)
    lambda_xy_mat = np.zeros((len(en),len(en)),dtype = np.complex128)
    lambda_yx_mat = np.zeros((len(en),len(en)),dtype = np.complex128)
    lambda_yy_mat = np.zeros((len(en),len(en)),dtype = np.complex128)

    for n in range(en):
        for m in range(en):
            for l in range(en):
                lambda_xx_mat[n][m] = lambda_xx_mat[n][m]-(((en[n]-en[l])*dipole_x_mat[n]*dipole_x_mat[l])-((en[l]-en[m])*dipole_x_mat[n][l]*dipole_x_mat[l][m]))
                lambda_xy_mat[n][m] = lambda_xy_mat[n][m]-(((en[n]-en[l])*dipole_x_mat[n]*dipole_y_mat[l])-((en[l]-en[m])*dipole_y_mat[n][l]*dipole_x_mat[l][m]))
                lambda_yx_mat[n][m] = lambda_yx_mat[n][m]-(((en[n]-en[l])*dipole_y_mat[n]*dipole_x_mat[l])-((en[l]-en[m])*dipole_x_mat[n][l]*dipole_y_mat[l][m]))
                lambda_yy_mat[n][m] = lambda_yy_mat[n][m]-(((en[n]-en[l])*dipole_y_mat[n]*dipole_y_mat[l])-((en[l]-en[m])*dipole_y_mat[n][l]*dipole_y_mat[l][m]))



    filename_current_xx = '%s/Lambda_xx_%d_x_sites_%d_y_sites_%d_particles_%0.2f_U_%0.2f_t1.pkl' %(Raw_nos,Nx,Ny,Np,U,t1)
    data_current_xx = lambda_xx_mat
    with open(filename_current_xx, 'wb') as outfile:
       pickle.dump(data_current_xx, outfile, pickle.HIGHEST_PROTOCOL)

    filename_current_xy = '%s/Lambda_xy_%d_x_sites_%d_y_sites_%d_particles_%0.2f_U_%0.2f_t1.pkl' %(Raw_nos,Nx,Ny,Np,U,t1)
    data_current_xy = lambda_xy_mat
    with open(filename_current_xy, 'wb') as outfile:
       pickle.dump(data_current_xy, outfile, pickle.HIGHEST_PROTOCOL)

    filename_current_yx = '%s/Lambda_yx_%d_x_sites_%d_y_sites_%d_particles_%0.2f_U_%0.2f_t1.pkl' %(Raw_nos,Nx,Ny,Np,U,t1)
    data_current_yx = lambda_yx_mat
    with open(filename_current_yx, 'wb') as outfile:
       pickle.dump(data_current_yx, outfile, pickle.HIGHEST_PROTOCOL)

    filename_current_yy = '%s/Lambda_yy_%d_x_sites_%d_y_sites_%d_particles_%0.2f_U_%0.2f_t1.pkl' %(Raw_nos,Nx,Ny,Np,U,t1)
    data_current_yy = lambda_yy_mat
    with open(filename_current_yy, 'wb') as outfile:
       pickle.dump(data_current_yy, outfile, pickle.HIGHEST_PROTOCOL)

 
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
# =============================================================================
    rx = np.zeros(N)
    ry = np.zeros(N)
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
    Raw_nos = "Text_files"
    dipole_op_con(t,t1,a,rx,ry,u,basis,Nx,Ny,N,Np,Raw_nos,data_dir)    

if __name__ == '__main__':
    sys.argv
    total = len(sys.argv)  
    print("No of sys arguments",total)
    cmdargs = sys.argv
    main(total,cmdargs)   
