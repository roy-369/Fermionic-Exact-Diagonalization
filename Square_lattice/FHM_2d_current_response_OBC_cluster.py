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


    

def Delta(x,zeta):
    d = (zeta**2)/(np.pi*(x**2+zeta**2))
    return d

def Principal(x,zeta):
    d = x/((x**2+zeta**2))
    return d

def paramagnetic_current(energies,omega,eta,zeta,Lambda_1x,Lambda_1y):
     

     sigma_xx = 0
     sigma_xy = 0
     sigma_yx = 0
     sigma_yy = 0
     
     for n in range(len(energies)):
              if(n==0):
                  f_n = 1
              else:
                  f_n = 0
              for m in range(len(energies)):

                  if(m==0):
                      f_m = 1
                  else:
                      f_m = 0
                  occ = f_n-f_m
                  sigma_xx = sigma_xx + (Lambda_1x[n][m]*Lambda_1x[m][n]*occ)/((energies[n]-energies[m]+(omega+1j*eta))*(1j*(omega+1j*eta)))
                  sigma_xy = sigma_xy + (Lambda_1x[n][m]*Lambda_1y[m][n]*occ)/((energies[n]-energies[m]+(omega+1j*eta))*(1j*(omega+1j*eta)))
                  sigma_yx = sigma_yx + (Lambda_1y[n][m]*Lambda_1x[m][n]*occ)/((energies[n]-energies[m]+(omega+1j*eta))*(1j*(omega+1j*eta)))
                  sigma_yy = sigma_yy + (Lambda_1y[n][m]*Lambda_1y[m][n]*occ)/((energies[n]-energies[m]+(omega+1j*eta))*(1j*(omega+1j*eta)))
                  
     return sigma_xx,sigma_xy,sigma_yx,sigma_yy


def diamagnetic_current_1(energies,omega,eta,zeta,Lambda_1x,Lambda_1y):
     

     sigma_xx_dia = 0
     sigma_xy_dia = 0
     sigma_yx_dia = 0
     sigma_yy_dia = 0
     for n in range(len(energies)):
              if(n==0):
                  f_n = 1
              else:
                  f_n = 0
              for m in range(len(energies)):
                  if(m==0):
                      f_m = 1
                  else:
                      f_m = 0
                  occ = f_n-f_m
                  e_diff = energies[n]-energies[m]
                  num   = Principal(e_diff,zeta)
                  sigma_xx_dia = sigma_xx_dia + (Lambda_1x[n][m]*Lambda_1x[m][n]*occ*num)/((1j*(omega+1j*eta)))
                  sigma_xy_dia = sigma_xy_dia + (Lambda_1x[n][m]*Lambda_1y[m][n]*occ*num)/((1j*(omega+1j*eta)))
                  sigma_yx_dia = sigma_yx_dia + (Lambda_1y[n][m]*Lambda_1x[m][n]*occ*num)/((1j*(omega+1j*eta)))
                  sigma_yy_dia = sigma_yy_dia + (Lambda_1y[n][m]*Lambda_1y[m][n]*occ*num)/((1j*(omega+1j*eta)))
                  
     return sigma_xx_dia, sigma_xy_dia, sigma_yx_dia, sigma_yy_dia

    
def linear_currents(Nx,Ny,Np,U,t1,N,Raw_nos,Lambda_matrix_1x,Lambda_matrix_1y,energies, omega,eta,zeta):

    
    sigma_xx_para = np.zeros(len(omega),dtype = np.complex128)
    sigma_xy_para = np.zeros(len(omega),dtype = np.complex128)
    sigma_yx_para = np.zeros(len(omega),dtype = np.complex128)
    sigma_yy_para = np.zeros(len(omega),dtype = np.complex128)
    sigma_xx_dia = np.zeros(len(omega),dtype = np.complex128)
    sigma_xy_dia = np.zeros(len(omega),dtype = np.complex128)
    sigma_yx_dia = np.zeros(len(omega),dtype = np.complex128)
    sigma_yy_dia = np.zeros(len(omega),dtype = np.complex128)
    
    for i in range(len(omega)):
        sigma_xx_para[i],sigma_xy_para[i],sigma_yx_para[i],sigma_yy_para[i] = paramagnetic_current(energies,omega[i],eta,zeta,Lambda_matrix_1x,Lambda_matrix_1y)
        sigma_xx_dia[i],sigma_xy_dia[i],sigma_yx_dia[i],sigma_yy_dia[i]= diamagnetic_current_1(energies,omega[i],eta,zeta,Lambda_matrix_1x,Lambda_matrix_1y)
    
    
    filename_sigma_xx_para = '%s/Sigma_xx_para_%d_x_sites_%d_y_sites_%d_particles_%0.2f_U_%0.2f_t1.pkl' %(Raw_nos,Nx,Ny,Np,U,t1)
    data_sigma_xx_para = sigma_xx_para
    with open(filename_sigma_xx_para, 'wb') as outfile:
       pickle.dump(data_sigma_xx_para, outfile, pickle.HIGHEST_PROTOCOL)
     
    filename_sigma_xy_para = '%s/Sigma_xy_para_%d_x_sites_%d_y_sites_%d_particles_%0.2f_U_%0.2f_t1.pkl' %(Raw_nos,Nx,Ny,Np,U,t1)
    data_sigma_xy_para = sigma_xy_para
    with open(filename_sigma_xy_para, 'wb') as outfile:
       pickle.dump(data_sigma_xy_para, outfile, pickle.HIGHEST_PROTOCOL)
       
       
    filename_sigma_yx_para = '%s/Sigma_yx_para_%d_x_sites_%d_y_sites_%d_particles_%0.2f_U_%0.2f_t1.pkl' %(Raw_nos,Nx,Ny,Np,U,t1)
    data_sigma_yx_para = sigma_yx_para
    with open(filename_sigma_yx_para, 'wb') as outfile:
       pickle.dump(data_sigma_yx_para, outfile, pickle.HIGHEST_PROTOCOL)
       
    filename_sigma_yy_para = '%s/Sigma_yy_para_%d_x_sites_%d_y_sites_%d_particles_%0.2f_U_%0.2f_t1.pkl' %(Raw_nos,Nx,Ny,Np,U,t1)
    data_sigma_yy_para = sigma_yy_para
    with open(filename_sigma_yy_para, 'wb') as outfile:
       pickle.dump(data_sigma_yy_para, outfile, pickle.HIGHEST_PROTOCOL)
    

 
    filename_sigma_xx_dia = '%s/Sigma_xx_dia_%d_x_sites_%d_y_sites_%d_particles_%0.2f_U_%0.2f_t1.pkl' %(Raw_nos,Nx,Ny,Np,U,t1)
    data_sigma_xx_dia = sigma_xx_dia
    with open(filename_sigma_xx_dia, 'wb') as outfile:
       pickle.dump(data_sigma_xx_dia, outfile, pickle.HIGHEST_PROTOCOL)    
     
    filename_sigma_xy_dia = '%s/Sigma_xy_dia_%d_x_sites_%d_y_sites_%d_particles_%0.2f_U_%0.2f_t1.pkl' %(Raw_nos,Nx,Ny,Np,U,t1)
    data_sigma_xy_dia = sigma_xy_dia
    with open(filename_sigma_xy_dia, 'wb') as outfile:
       pickle.dump(data_sigma_xy_dia, outfile, pickle.HIGHEST_PROTOCOL)
       
    filename_sigma_yx_dia = '%s/Sigma_yx_dia_%d_x_sites_%d_y_sites_%d_particles_%0.2f_U_%0.2f_t1.pkl' %(Raw_nos,Nx,Ny,Np,U,t1)
    data_sigma_yx_dia = sigma_yx_dia
    with open(filename_sigma_yx_dia, 'wb') as outfile:
       pickle.dump(data_sigma_yx_dia, outfile, pickle.HIGHEST_PROTOCOL)
     
    filename_sigma_yy_dia = '%s/Sigma_yy_dia_%d_x_sites_%d_y_sites_%d_particles_%0.2f_U_%0.2f_t1.pkl' %(Raw_nos,Nx,Ny,Np,U,t1)
    data_sigma_yy_dia = sigma_yy_dia
    with open(filename_sigma_yy_dia, 'wb') as outfile:
       pickle.dump(data_sigma_yy_dia, outfile, pickle.HIGHEST_PROTOCOL)

    
def nonlinear_current_lp(Np,N,U,Raw_nos,Lambda_matrix_1x,Lambda_matrix_1y,dipole_moment_x,dipole_moment_y,energies,omega,eta):

    chi_x_xx_lp = np.zeros(len(omega),dtype = np.complex128)
    chi_x_xy_lp = np.zeros(len(omega),dtype = np.complex128)
    chi_x_yx_lp = np.zeros(len(omega),dtype = np.complex128)
    chi_x_yy_lp = np.zeros(len(omega),dtype = np.complex128)
    chi_y_xx_lp = np.zeros(len(omega),dtype = np.complex128)
    chi_y_xy_lp = np.zeros(len(omega),dtype = np.complex128)
    chi_y_yx_lp = np.zeros(len(omega),dtype = np.complex128)
    chi_y_yy_lp = np.zeros(len(omega),dtype = np.complex128)
    for i in range(len(omega)):
        for n in range(len(energies)):
            if(n == 0):
                  f_n = 1
            else:
                  f_n = 0
                 
            for m in range(len(energies)):
                 if(m == 0):
                      f_m = 1
                 else:
                      f_m = 0
                 occ = f_n-f_m
                 chi_x_xx_lp[i] = chi_x_xx_lp[i]+(-1j*occ*Lambda_matrix_1x[n][m]*Lambda_matrix_1x[m][n]*(dipole_moment_x[n][n]-dipole_moment_x[m][m]))/(energies[n]-energies[m]+(omega[i]+1j*eta))
                 chi_x_xy_lp[i] = chi_x_xy_lp[i]+(-1j*occ*Lambda_matrix_1x[n][m]*Lambda_matrix_1y[m][n]*(dipole_moment_x[n][n]-dipole_moment_x[m][m]))/(energies[n]-energies[m]+(omega[i]+1j*eta))
                 chi_x_yx_lp[i] = chi_x_yx_lp[i]+(-1j*occ*Lambda_matrix_1y[n][m]*Lambda_matrix_1x[m][n]*(dipole_moment_x[n][n]-dipole_moment_x[m][m]))/(energies[n]-energies[m]+(omega[i]+1j*eta))
                 chi_x_yy_lp[i] = chi_x_yy_lp[i]+(-1j*occ*Lambda_matrix_1y[n][m]*Lambda_matrix_1y[m][n]*(dipole_moment_x[n][n]-dipole_moment_x[m][m]))/(energies[n]-energies[m]+(omega[i]+1j*eta))
                 chi_y_xx_lp[i] = chi_y_xx_lp[i]+(-1j*occ*Lambda_matrix_1x[n][m]*Lambda_matrix_1x[m][n]*(dipole_moment_y[n][n]-dipole_moment_y[m][m]))/(energies[n]-energies[m]+(omega[i]+1j*eta))
                 chi_y_xy_lp[i] = chi_y_xy_lp[i]+(-1j*occ*Lambda_matrix_1x[n][m]*Lambda_matrix_1y[m][n]*(dipole_moment_y[n][n]-dipole_moment_y[m][m]))/(energies[n]-energies[m]+(omega[i]+1j*eta))
                 chi_y_yx_lp[i] = chi_y_yx_lp[i]+(-1j*occ*Lambda_matrix_1y[n][m]*Lambda_matrix_1x[m][n]*(dipole_moment_y[n][n]-dipole_moment_y[m][m]))/(energies[n]-energies[m]+(omega[i]+1j*eta))
                 chi_y_yy_lp[i] = chi_y_yy_lp[i]+(-1j*occ*Lambda_matrix_1y[n][m]*Lambda_matrix_1y[m][n]*(dipole_moment_y[n][n]-dipole_moment_y[m][m]))/(energies[n]-energies[m]+(omega[i]+1j*eta))

   
    
    filename_chi_x_xx_lp = '%s/Chi_x_xx_shift_linear_polarized_%d_x_sites_%d_y_sites_%d_particles_%0.2f_U_%0.2f_t1.pkl' %(Raw_nos,Nx,Ny,Np,U,t1)
    data_chi_x_xx_lp = chi_x_xx_lp
    with open(filename_chi_x_xx_lp, 'wb') as outfile:
       pickle.dump(data_chi_x_xx_lp, outfile, pickle.HIGHEST_PROTOCOL) 
    
    filename_chi_x_xy_lp = '%s/Chi_x_xy_shift_linear_polarized_%d_x_sites_%d_y_sites_%d_particles_%0.2f_U_%0.2f_t1.pkl' %(Raw_nos,Nx,Ny,Np,U,t1)
    data_chi_x_xy_lp = chi_x_xy_lp
    with open(filename_chi_x_xy_lp, 'wb') as outfile:
       pickle.dump(data_chi_x_xy_lp, outfile, pickle.HIGHEST_PROTOCOL) 
       
    filename_chi_x_yx_lp = '%s/Chi_x_yx_shift_linear_polarized_%d_x_sites_%d_y_sites_%d_particles_%0.2f_U_%0.2f_t1.pkl' %(Raw_nos,Nx,Ny,Np,U,t1)
    data_chi_x_yx_lp = chi_x_yx_lp
    with open(filename_chi_x_yx_lp, 'wb') as outfile:
       pickle.dump(data_chi_x_yx_lp, outfile, pickle.HIGHEST_PROTOCOL) 
       
    filename_chi_x_yy_lp = '%s/Chi_x_yy_shift_linear_polarized_%d_x_sites_%d_y_sites_%d_particles_%0.2f_U_%0.2f_t1.pkl' %(Raw_nos,Nx,Ny,Np,U,t1)
    data_chi_x_yy_lp = chi_x_yy_lp
    with open(filename_chi_x_yy_lp, 'wb') as outfile:
       pickle.dump(data_chi_x_yy_lp, outfile, pickle.HIGHEST_PROTOCOL) 
       
    filename_chi_y_xx_lp = '%s/Chi_y_xx_shift_linear_polarized_%d_x_sites_%d_y_sites_%d_particles_%0.2f_U_%0.2f_t1.pkl' %(Raw_nos,Nx,Ny,Np,U,t1)
    data_chi_y_xx_lp = chi_y_xx_lp
    with open(filename_chi_y_xx_lp, 'wb') as outfile:
       pickle.dump(data_chi_y_xx_lp, outfile, pickle.HIGHEST_PROTOCOL) 
       
    filename_chi_y_xy_lp = '%s/Chi_y_xy_shift_linear_polarized_%d_x_sites_%d_y_sites_%d_particles_%0.2f_U_%0.2f_t1.pkl' %(Raw_nos,Nx,Ny,Np,U,t1)
    data_chi_y_xy_lp = chi_y_xy_lp
    with open(filename_chi_y_xy_lp, 'wb') as outfile:
       pickle.dump(data_chi_y_xy_lp, outfile, pickle.HIGHEST_PROTOCOL) 
       
    filename_chi_y_yx_lp = '%s/Chi_y_yx_shift_linear_polarized_%d_x_sites_%d_y_sites_%d_particles_%0.2f_U_%0.2f_t1.pkl' %(Raw_nos,Nx,Ny,Np,U,t1)
    data_chi_y_yx_lp = chi_y_yx_lp
    with open(filename_chi_y_yx_lp, 'wb') as outfile:
       pickle.dump(data_chi_y_yx_lp, outfile, pickle.HIGHEST_PROTOCOL) 
       
    filename_chi_y_yy_lp = '%s/Chi_y_yy_shift_linear_polarized_%d_x_sites_%d_y_sites_%d_particles_%0.2f_U_%0.2f_t1.pkl' %(Raw_nos,Nx,Ny,Np,U,t1)
    data_chi_y_yy_lp = chi_y_yy_lp
    with open(filename_chi_y_yy_lp, 'wb') as outfile:
       pickle.dump(data_chi_y_yy_lp, outfile, pickle.HIGHEST_PROTOCOL) 
       
       
    #print(chi_x_xx_lp)
    #print(np.amax(np.real(chi_x_xx_lp)))


def nl_currents_SHG(Np,N,U,Raw_nos,energies,Lambda_matrix_1x,Lambda_matrix_1y,dipole_moment_x,dipole_moment_y,omega,eta):
    
    chi_x_xx_SHG = np.zeros(len(omega),dtype = np.complex128)
    chi_x_xy_SHG = np.zeros(len(omega),dtype = np.complex128)
    chi_x_yx_SHG = np.zeros(len(omega),dtype = np.complex128)
    chi_x_yy_SHG = np.zeros(len(omega),dtype = np.complex128)
    chi_y_xx_SHG = np.zeros(len(omega),dtype = np.complex128)
    chi_y_xy_SHG = np.zeros(len(omega),dtype = np.complex128)
    chi_y_yx_SHG = np.zeros(len(omega),dtype = np.complex128)
    chi_y_yy_SHG = np.zeros(len(omega),dtype = np.complex128)

    
    for i in range(len(omega)):
        for n in range(len(energies)):
            if(energies[n] == np.min(energies)):
                  f_n = 1
            else:
                  f_n = 0
                 
            for m in range(len(energies)):
                 if(energies[m] == np.min(energies)):
                      f_m = 1
                 else:
                      f_m = 0
                 occ = f_n-f_m
                 A_a = -2*occ*(omega[i]+1j*eta)
                 A_b = 1/(energies[n]-energies[m]+(omega[i]+1j*eta))
                 for l in range(len(energies)):
                     
                     
                      a_x_xx_1 = ((Lambda_matrix_1x[n][l]*dipole_moment_x[l][m]-dipole_moment_x[n][l]*Lambda_matrix_1x[l][m])*dipole_moment_x[m][n])/(2*(energies[n]-energies[m]+2*(omega[i]+1j*eta)))          
                      a_x_xx_2 = (((energies[m]-energies[l])*dipole_moment_x[m][l]*dipole_moment_x[l][n])*Lambda_matrix_1x[n][m])/(energies[n]-energies[l]+2*(omega[i]+1j*eta))
                      a_x_xx_3 = (((energies[l]-energies[n])*dipole_moment_x[m][l]*dipole_moment_x[l][n])*Lambda_matrix_1x[n][m])/(energies[l]-energies[m]+2*(omega[i]+1j*eta))
                      
                      chi_x_xx_SHG[i] = chi_x_xx_SHG[i]+A_a*(a_x_xx_1+A_b*(a_x_xx_2-a_x_xx_3))
                      
                      a_x_xy_1 = ((Lambda_matrix_1x[n][l]*dipole_moment_y[l][m]-dipole_moment_y[n][l]*Lambda_matrix_1x[l][m])*dipole_moment_x[m][n])/(2*(energies[n]-energies[m]+2*(omega[i]+1j*eta)))          
                      a_x_xy_2 = (((energies[m]-energies[l])*dipole_moment_y[m][l]*dipole_moment_x[l][n])*Lambda_matrix_1x[n][m])/(energies[n]-energies[l]+2*(omega[i]+1j*eta))
                      a_x_xy_3 = (((energies[l]-energies[n])*dipole_moment_x[m][l]*dipole_moment_y[l][n])*Lambda_matrix_1x[n][m])/(energies[l]-energies[m]+2*(omega[i]+1j*eta))
                      
                      chi_x_xy_SHG[i] = chi_x_xy_SHG[i]+A_a*(a_x_xy_1+A_b*(a_x_xy_2-a_x_xy_3))
                      
                      a_x_yx_1 = ((Lambda_matrix_1y[n][l]*dipole_moment_x[l][m]-dipole_moment_x[n][l]*Lambda_matrix_1y[l][m])*dipole_moment_x[m][n])/(2*(energies[n]-energies[m]+2*(omega[i]+1j*eta)))          
                      a_x_yx_2 = (((energies[m]-energies[l])*dipole_moment_x[m][l]*dipole_moment_x[l][n])*Lambda_matrix_1y[n][m])/(energies[n]-energies[l]+2*(omega[i]+1j*eta))
                      a_x_yx_3 = (((energies[l]-energies[n])*dipole_moment_x[m][l]*dipole_moment_x[l][n])*Lambda_matrix_1y[n][m])/(energies[l]-energies[m]+2*(omega[i]+1j*eta))
                      
                      chi_x_yx_SHG[i] = chi_x_yx_SHG[i]+A_a*(a_x_yx_1+A_b*(a_x_yx_2-a_x_yx_3))
                      
                      a_x_yy_1 = ((Lambda_matrix_1y[n][l]*dipole_moment_y[l][m]-dipole_moment_y[n][l]*Lambda_matrix_1y[l][m])*dipole_moment_x[m][n])/(2*(energies[n]-energies[m]+2*(omega[i]+1j*eta)))          
                      a_x_yy_2 = (((energies[m]-energies[l])*dipole_moment_y[m][l]*dipole_moment_x[l][n])*Lambda_matrix_1y[n][m])/(energies[n]-energies[l]+2*(omega[i]+1j*eta))
                      a_x_yy_3 = (((energies[l]-energies[n])*dipole_moment_y[m][l]*dipole_moment_x[l][n])*Lambda_matrix_1y[n][m])/(energies[l]-energies[m]+2*(omega[i]+1j*eta))
                      
                      chi_x_yy_SHG[i] = chi_x_yy_SHG[i]+A_a*(a_x_yy_1+A_b*(a_x_yy_2-a_x_yy_3))
                      
                      
                      a_y_xx_1 = ((Lambda_matrix_1x[n][l]*dipole_moment_x[l][m]-dipole_moment_x[n][l]*Lambda_matrix_1x[l][m])*dipole_moment_y[m][n])/(2*(energies[n]-energies[m]+2*(omega[i]+1j*eta)))          
                      a_y_xx_2 = (((energies[m]-energies[l])*dipole_moment_x[m][l]*dipole_moment_y[l][n])*Lambda_matrix_1x[n][m])/(energies[n]-energies[l]+2*(omega[i]+1j*eta))
                      a_y_xx_3 = (((energies[l]-energies[n])*dipole_moment_y[m][l]*dipole_moment_x[l][n])*Lambda_matrix_1x[n][m])/(energies[l]-energies[m]+2*(omega[i]+1j*eta))
                      
                      chi_y_xx_SHG[i] = chi_y_xx_SHG[i]+A_a*(a_y_xx_1+A_b*(a_y_xx_2-a_y_xx_3))
                      
                      a_y_xy_1 = ((Lambda_matrix_1x[n][l]*dipole_moment_y[l][m]-dipole_moment_y[n][l]*Lambda_matrix_1x[l][m])*dipole_moment_y[m][n])/(2*(energies[n]-energies[m]+2*(omega[i]+1j*eta)))          
                      a_y_xy_2 = (((energies[m]-energies[l])*dipole_moment_y[m][l]*dipole_moment_y[l][n])*Lambda_matrix_1x[n][m])/(energies[n]-energies[l]+2*(omega[i]+1j*eta))
                      a_y_xy_3 = (((energies[l]-energies[n])*dipole_moment_y[m][l]*dipole_moment_y[l][n])*Lambda_matrix_1x[n][m])/(energies[l]-energies[m]+2*(omega[i]+1j*eta))
                      
                      chi_x_xy_SHG[i] = chi_y_xy_SHG[i]+A_a*(a_y_xy_1+A_b*(a_y_xy_2-a_y_xy_3))
                      
                      a_y_yx_1 = ((Lambda_matrix_1y[n][l]*dipole_moment_x[l][m]-dipole_moment_x[n][l]*Lambda_matrix_1y[l][m])*dipole_moment_y[m][n])/(2*(energies[n]-energies[m]+2*(omega[i]+1j*eta)))          
                      a_y_yx_2 = (((energies[m]-energies[l])*dipole_moment_x[m][l]*dipole_moment_y[l][n])*Lambda_matrix_1y[n][m])/(energies[n]-energies[l]+2*(omega[i]+1j*eta))
                      a_y_yx_3 = (((energies[l]-energies[n])*dipole_moment_y[m][l]*dipole_moment_x[l][n])*Lambda_matrix_1y[n][m])/(energies[l]-energies[m]+2*(omega[i]+1j*eta))
                      
                      chi_y_yx_SHG[i] = chi_y_yx_SHG[i]+A_a*(a_y_yx_1+A_b*(a_y_yx_2-a_y_yx_3))
                      
                      a_y_yy_1 = ((Lambda_matrix_1y[n][l]*dipole_moment_y[l][m]-dipole_moment_y[n][l]*Lambda_matrix_1y[l][m])*dipole_moment_y[m][n])/(2*(energies[n]-energies[m]+2*(omega[i]+1j*eta)))          
                      a_y_yy_2 = (((energies[m]-energies[l])*dipole_moment_y[m][l]*dipole_moment_y[l][n])*Lambda_matrix_1y[n][m])/(energies[n]-energies[l]+2*(omega[i]+1j*eta))
                      a_y_yy_3 = (((energies[l]-energies[n])*dipole_moment_y[m][l]*dipole_moment_y[l][n])*Lambda_matrix_1y[n][m])/(energies[l]-energies[m]+2*(omega[i]+1j*eta))
                      
                      chi_y_yy_SHG[i] = chi_y_yy_SHG[i]+A_a*(a_y_yy_1+A_b*(a_y_yy_2-a_y_yy_3))
        
    
    filename_chi_x_xx_SHG = '%s/Chi_x_xx_SHG_%d_x_sites_%d_y_sites_%d_particles_%0.2f_U_%0.2f_t1.pkl' %(Raw_nos,Nx,Ny,Np,U,t1)
    data_chi_x_xx_SHG = chi_x_xx_SHG
    with open(filename_chi_x_xx_SHG, 'wb') as outfile:
       pickle.dump(data_chi_x_xx_SHG, outfile, pickle.HIGHEST_PROTOCOL)
    
    filename_chi_x_xy_SHG = '%s/Chi_x_xy_SHG_%d_x_sites_%d_y_sites_%d_particles_%0.2f_U_%0.2f_t1.pkl' %(Raw_nos,Nx,Ny,Np,U,t1)
    data_chi_x_xy_SHG = chi_x_xy_SHG
    with open(filename_chi_x_xy_SHG, 'wb') as outfile:
       pickle.dump(data_chi_x_xy_SHG, outfile, pickle.HIGHEST_PROTOCOL)
       
    filename_chi_x_yx_SHG = '%s/Chi_x_yx_SHG_%d_x_sites_%d_y_sites_%d_particles_%0.2f_U_%0.2f_t1.pkl' %(Raw_nos,Nx,Ny,Np,U,t1)
    data_chi_x_yx_SHG = chi_x_yx_SHG
    with open(filename_chi_x_yx_SHG, 'wb') as outfile:
       pickle.dump(data_chi_x_yx_SHG, outfile, pickle.HIGHEST_PROTOCOL)
       
    filename_chi_x_yy_SHG = '%s/Chi_x_yy_SHG_%d_x_sites_%d_y_sites_%d_particles_%0.2f_U_%0.2f_t1.pkl' %(Raw_nos,Nx,Ny,Np,U,t1)
    data_chi_x_yy_SHG = chi_x_yy_SHG
    with open(filename_chi_x_yy_SHG, 'wb') as outfile:
       pickle.dump(data_chi_x_yy_SHG, outfile, pickle.HIGHEST_PROTOCOL)
       
       
    filename_chi_y_xx_SHG = '%s/Chi_y_xx_SHG_%d_x_sites_%d_y_sites_%d_particles_%0.2f_U_%0.2f_t1.pkl' %(Raw_nos,Nx,Ny,Np,U,t1)
    data_chi_y_xx_SHG = chi_y_xx_SHG
    with open(filename_chi_y_xx_SHG, 'wb') as outfile:
       pickle.dump(data_chi_y_xx_SHG, outfile, pickle.HIGHEST_PROTOCOL)
    
    filename_chi_y_xy_SHG = '%s/Chi_y_xy_SHG_%d_x_sites_%d_y_sites_%d_particles_%0.2f_U_%0.2f_t1.pkl' %(Raw_nos,Nx,Ny,Np,U,t1)
    data_chi_y_xy_SHG = chi_y_xy_SHG
    with open(filename_chi_y_xy_SHG, 'wb') as outfile:
       pickle.dump(data_chi_y_xy_SHG, outfile, pickle.HIGHEST_PROTOCOL)
       
    filename_chi_y_yx_SHG = '%s/Chi_y_yx_SHG_%d_x_sites_%d_y_sites_%d_particles_%0.2f_U_%0.2f_t1.pkl' %(Raw_nos,Nx,Ny,Np,U,t1)
    data_chi_y_yx_SHG = chi_y_yx_SHG
    with open(filename_chi_y_yx_SHG, 'wb') as outfile:
       pickle.dump(data_chi_y_yx_SHG, outfile, pickle.HIGHEST_PROTOCOL)
       
    filename_chi_y_yy_SHG = '%s/Chi_y_yy_SHG_%d_x_sites_%d_y_sites_%d_particles_%0.2f_U_%0.2f_t1.pkl' %(Raw_nos,Nx,Ny,Np,U,t1)
    data_chi_y_yy_SHG = chi_y_yy_SHG
    with open(filename_chi_y_yy_SHG, 'wb') as outfile:
       pickle.dump(data_chi_y_yy_SHG, outfile, pickle.HIGHEST_PROTOCOL)



 
def main(total,cmdargs):
    if(total!=6):
        raise ValueError('missing args')
    
        
    
    
    
    Nx_dir = cmdargs[1]
    Ny_dir = cmdargs[2]
    u_dir = cmdargs[3]
    t1_dir = cmdargs[4]
    
    Nx = int(Nx_dir)
    Ny = int(Ny_dir)
    U= float(u_dir)
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
    Raw_nos = "Text_files" 
    
    data_dir = "Text_files"
   
    

    filename_energy = '%s/Energies_%d_x_sites_%d_y_sites_%d_particles_%0.2f_U_%0.2f_t1.pkl' %(Raw_nos,Nx,Ny,Np,U,t1)
    with open(filename_energy, 'rb') as infile:
       en = pickle.load(infile)

    filename_dipole_x = '%s/Dipole_x_%d_x_sites_%d_y_sites_%d_particles_%0.2f_U_%0.2f_t1.pkl' %(Raw_nos,Nx,Ny,Np,U,t1)
    with open(filename_dipole_x, 'rb') as infile:
       dipole_x_mat = pickle.load(infile)

    filename_dipole_y = '%s/Dipole_y_%d_x_sites_%d_y_sites_%d_particles_%0.2f_U_%0.2f_t1.pkl' %(Raw_nos,Nx,Ny,Np,U,t1)
    with open(filename_dipole_y, 'rb') as infile:
       dipole_y_mat = pickle.load(infile)

    filename_current_x = '%s/Lambda_x_%d_x_sites_%d_y_sites_%d_particles_%0.2f_U_%0.2f_t1.pkl' %(Raw_nos,Nx,Ny,Np,U,t1)
    with open(filename_current_x, 'rb') as infile:
       lambda_x_mat = pickle.load(infile)

    filename_current_y = '%s/Lambda_y_%d_x_sites_%d_y_sites_%d_particles_%0.2f_U_%0.2f_t1.pkl' %(Raw_nos,Nx,Ny,Np,U,t1)
    with open(filename_current_y, 'rb') as infile:
       lambda_y_mat = pickle.load(infile)
    
    en_min = np.nanmin(np.real(en))
    en_max = np.nanmax(np.real(en))
    en_width = np.abs(en_max-en_min)
    d_omega = 0.005
    eta = 0.001
    zeta = 4*d_omega
    omega = np.arange(-2*en_width,2*en_width,d_omega)
    
    filename_omega = '%s/Frequency_space_%d_x_sites_%d_y_sites_%d_particles_%0.2f_U_%0.2f_t1.pkl' %(Raw_nos,Nx,Ny,Np,U,t1)
    data_omega = omega
    with open(filename_omega, 'wb') as outfile:
       pickle.dump(data_omega, outfile, pickle.HIGHEST_PROTOCOL)

    linear_currents(Nx,Ny,Np,U,t1,N,Raw_nos,lambda_x_mat,lambda_y_mat,en,omega,eta,zeta)
    
if __name__ == '__main__':
    sys.argv
    total = len(sys.argv)  
    print("No of sys arguments",total)
    cmdargs = sys.argv
    main(total,cmdargs)   
