#!/usr/bin/python
#
# This script provides the calculation of the susceptibility in the different channels from the ED vertex
#
# In the summation over the internal indices, the kernel 1 extracted in the center of the vertex function (scanning over the bosonic frequency at (nu, nu')=(0,0)) has been used for the vertex asymptotic
# The different suceptibilities will be will be stored in the HDF5 FILE which already contains all the info about the ED calculation of the vertex
#
########################################################################################
#
# WARNING: modify the path and the filename to load !!
#
#########################################################################################

#--------------------------------------IMPORTS ------------------------------------------

import h5py
import matplotlib.pyplot as pl
import numpy as np
import os
import sys
import subprocess
#import matplotlib.cm as cm
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from os.path import isfile, join
import math
import numpy 
from agneselib.mymath import *
from agneselib.translate_notation import *
from agneselib.inf_Msum import *
from _functools import partial

#-----------------------------------Read parameters----------------------------------------------
# ------U---------
U=1.0
##-------BETA-------
beta=20.0
#----------------------------------------Read HDF5 files-----------------------------------------

#fRG 2loop full treatement for the vertex asymptotics data 

if ('../dat/'):
    f = h5py.File('../dat/fRG/full/2loop/dat_U1_Beta20_PFCB4096_EDB_OMFL_KAT_2LOOP_SU2.h5', 'r+')   # Read (and write) the hdf5 file in the directory "dat" if existing
else:
    sys.exit("No data file")

#---------  PARAMETERS  ---------------

pi = math.pi

#---------------------- Read objects involved in the calculation from HDF5 file -------------------------------

#--------------------------------- READING ----------------------------------------------------------

#GF

re_g_iw = f["Giw/RE"][:,0,0,0]
im_g_iw = f["Giw/IM"][:,0,0,0]           
print im_g_iw.shape[0] 
fgrid_gf = im_g_iw.shape[0]
N_fermi_gf = fgrid_gf/2


#VERTEX

#Full Vertex

re_f_updo_pf = np.array(f["Vert/RE"][:,:,:,0,0,0,0,0,0,0])

print re_f_updo_pf.shape

im_f_updo_pf = np.array(f["Vert/IM"][:,:,:,0,0,0,0,0,0,0])

# We assume all the channels to have the same B/F grids

fgrid = f["Vert/RE"][:].shape[0] 

N_fermi_pf = fgrid/2

N_bose = fgrid/3  # to create a bosonic frequency grid from -N_bose to N_bose (double the number of fermionic frequency)
N_fermi= (fgrid)/3   # to create a fermionic frequency grid from -N_fermi to N_fermi

fgrid_arr = np.array([(2*n +1)*pi/beta for n in range(-N_fermi, N_fermi)])
bgrid_arr = np.array([(2*n)*pi/beta for n in range(-N_bose, N_bose+1)])

print N_bose, N_fermi

#-------------------------------Reading asymptotics: in the karrash approx there is only K1_effective------------------------
#Kernel 1 effective

rek_updo_ph =  np.array(f["/chi_func/RE_PH"][:,0,0,0,0,0])
imk_updo_ph = np.array(f["/chi_func/IM_PH"][:,0,0,0,0,0])
rek_updo_pp =  np.array(f["/chi_func/RE_PP"][:,0,0,0,0,0])
imk_updo_pp = np.array(f["/chi_func/IM_PP"][:,0,0,0,0,0])
rek_updo_xph =  np.array(f["/chi_func/RE_XPH"][:,0,0,0,0,0])
imk_updo_xph = np.array(f["/chi_func/IM_XPH"][:,0,0,0,0,0])

rek_upup_ph =  np.array(f["/chi_func/RE_PH"][:,0,0,0,0,0]-f["/chi_func/RE_XPH"][:,0,0,0,0,0])
imk_upup_ph =  np.array(f["/chi_func/IM_PH"][:,0,0,0,0,0]-f["/chi_func/IM_XPH"][:,0,0,0,0,0])
rek_upup_xph = -rek_upup_ph
imk_upup_xph = -imk_upup_ph

N_bose_k = (rek_updo_ph.shape[0]-1)/2 #We assume all channels to have the same bosonic range

#Kernel 2
 
rep_updo_ph =  np.array(f["/P_func/RE_PH"][:,:,0,0,0,0,0,0])
imp_updo_ph = np.array(f["/P_func/IM_PH"][:,:,0,0,0,0,0,0])
rep_updo_pp =  np.array(f["/P_func/RE_PP"][:,:,0,0,0,0,0,0])
imp_updo_pp = np.array(f["/P_func/IM_PP"][:,:,0,0,0,0,0,0])
rep_updo_xph =  np.array(f["/P_func/RE_XPH"][:,:,0,0,0,0,0,0])
imp_updo_xph = np.array(f["/P_func/IM_XPH"][:,:,0,0,0,0,0,0])


N_bose_p = (rep_updo_ph.shape[0]-1)/2 #We assume all channels to have the same bosonic range
N_fermi_p = (rep_updo_ph.shape[1])/2 #We assume all channels to have the same fermionic range

#Printing of the different ranges

print 'GF grid:', N_fermi_gf
print 'N_bose_k:', N_bose_k
print 'N_bose_p, N_fermi_p:', N_bose_p, N_fermi_p

#----------------------- FREQUENCY RANGES OF THE ASYMPTOTIC STRUCTURES-----------------------
# THE FFREQUENCY RANGE AND THE SUMMATION RANGE HAS TO BE >> THAN THE BFREQUENCY RANGE IN ORDER TO RESOLVE THE STRUCTURE

#Susceptibility function box:
N_bose_suscept = 100
bgrid_suscept= 2*N_bose_suscept + 1

#Internal double summation fermionic range
N_fermi_sum = 2*N_bose_suscept
N_fermi_sum1 = 500
print N_fermi_sum,N_fermi_sum1,N_bose_suscept

wNum = int(N_fermi_sum/100.0*10.0) # 10% of the the range/2 of my finite Matsubara summation
iMin = N_fermi_sum-wNum

wNum1 = int(N_fermi_sum1/100.0*10.0) # 10% of the the range/2 of my finite Matsubara summation
iMin1 = N_fermi_sum1-wNum1

Nl = 5  #Order of my fitting function (see library inf_Msum)

print wNum, wNum1

sum_fit = generate_sum_func(iMin1,wNum1,Nl) #It returns the Matsubara summation function {-infty,infty}
dblsum_fit = generate_dblsum_func(iMin,wNum,Nl) #It returns double Matsubara summation function {-infty,Infty}

#--------------------------------------- Function definitions------------------------

#Kernel 1

def K_upup_ph(wb):
    if (abs(wb) <= N_bose_k):
        return rek_updo_ph[wb+N_bose_k]+1j*imk_updo_ph[wb+N_bose_k]-rek_updo_xph[wb+N_bose_k]-1j*imk_updo_xph[wb+N_bose_k]
    else:
        return 0.0
def K_updo_ph(wb):
    if (abs(wb) <= N_bose_k):
        return rek_updo_ph[wb+N_bose_k]+1j*imk_updo_ph[wb+N_bose_k]
    else:
        return 0.0

def K_upup_pp(wb):
    return 0.0

def K_updo_pp(wb):
    if (abs(wb) <= N_bose_k):
        return rek_updo_pp[wb+N_bose_k]+1j*imk_updo_pp[wb+N_bose_k]
    else:
        return 0.0
def K_upup_xph(wb):
    if (abs(wb) <= N_bose_k):
        return -K_upup_ph(wb)
    else:
        return 0.0
def K_updo_xph(wb):
    if (abs(wb) <= N_bose_k):
        return rek_updo_xph[wb+N_bose_k]+1j*imk_updo_xph[wb+N_bose_k]
    else:
        return 0.0

#Kernel 2 funcion

def P_upup_ph(wb,wf):
    if (abs(wb) <= N_bose_p and wf >= -N_fermi_p and wf < N_fermi_p):
        return rep_updo_ph[wb+N_bose_p, wf+N_fermi_p]+1j*imp_updo_ph[wb+N_bose_p, wf+N_fermi_p]-rep_updo_xph[wb+N_bose_p, wf+N_fermi_p]-1j*imp_updo_xph[wb+N_bose_p,wf+N_fermi_p]
    else:
        return 0.0

def P_updo_ph(wb,wf):
    if (abs(wb) <= N_bose_p and wf >= -N_fermi_p and wf < N_fermi_p):
        return rep_updo_ph[wb+N_bose_p, wf+N_fermi_p]+1j*imp_updo_ph[wb+N_bose_p, wf+N_fermi_p]
    else:
        return 0.0

def P_upup_pp(wb,wf):
    return 0.0

def P_updo_pp(wb,wf):
    if (abs(wb) <= N_bose_p and wf >= -N_fermi_p and wf < N_fermi_p):
        return rep_updo_pp[wb+N_bose_p, wf+N_fermi_p]+1j*imp_updo_pp[wb+N_bose_p, wf+N_fermi_p]
    else:
        return 0.0

def P_upup_xph(wb,wf):
    return - P_upup_ph(wb,wf)

def P_updo_xph(wb,wf):
    if (abs(wb) <= N_bose_p and wf >= -N_fermi_p and wf < N_fermi_p):
        return rep_updo_xph[wb+N_bose_p, wf+N_fermi_p]+1j*imp_updo_xph[wb+N_bose_p, wf+N_fermi_p]
    else:
        return 0.0

#Full vertex
#purely fermionic function

def f_updo_fun_pf(i,j,k):
    return re_f_updo_pf[i+N_fermi_pf, j+N_fermi_pf, k+N_fermi_pf]+1j*im_f_updo_pf[i+N_fermi_pf, j+N_fermi_pf, k+N_fermi_pf]

def f_upup_fun_pf(i,j,k):
    return re_f_updo_pf[i+N_fermi_pf, j+N_fermi_pf, k+N_fermi_pf]+1j*im_f_updo_pf[i+N_fermi_pf, j+N_fermi_pf, k+N_fermi_pf]-re_f_updo_pf[i+N_fermi_pf, k+N_fermi_pf, j+N_fermi_pf]-1j*im_f_updo_pf[i+N_fermi_pf, k+N_fermi_pf, j+N_fermi_pf]

def isInside(i,j,k):
    return abs(i) <= N_bose and j >= -N_fermi and j < N_fermi and k >= -N_fermi and k < N_fermi

#Full VERTEX FUNCTIONS in the different channels

def f_upup_fun_ph(i,j,k):
    if isInside(i,j,k):    
        return f_upup_fun_pf(PHtoPF((i,j,k))[0], PHtoPF((i,j,k))[1], PHtoPF((i,j,k))[2])
    else:
        return K_upup_ph(i) + P_upup_ph(i,j)+ P_upup_ph(i,k)+ K_upup_xph(PHtoXPH((i,j,k))[0]) + K_upup_pp(PHtoPP((i,j,k))[0])+P_upup_xph(PHtoXPH((i,j,k))[0],PHtoXPH((i,j,k))[1])+P_upup_xph(PHtoXPH((i,j,k))[0],PHtoXPH((i,j,k))[2])+P_upup_pp(PHtoPP((i,j,k))[0],PHtoPP((i,j,k))[1])+P_upup_pp(PHtoPP((i,j,k))[0],PHtoPP((i,j,k))[2]) 

def f_upup_fun_pp(i,j,k):
    if isInside(i,j,k):
        return f_upup_fun_pf(PPtoPF((i,j,k))[0], PPtoPF((i,j,k))[1], PPtoPF((i,j,k))[2])
    else:
        return  K_upup_pp(i)+ P_upup_pp(i,j)+ P_upup_pp(i,k)+ K_upup_ph(PPtoPH((i,j,k))[0]) + K_upup_xph(PPtoXPH((i,j,k))[0])+P_upup_xph(PPtoXPH((i,j,k))[0],PPtoXPH((i,j,k))[1])+P_upup_xph(PPtoXPH((i,j,k))[0],PPtoXPH((i,j,k))[2])+P_upup_ph(PPtoPH((i,j,k))[0],PPtoPH((i,j,k))[1])+P_upup_ph(PPtoPH((i,j,k))[0],PPtoPH((i,j,k))[2]) 

def f_upup_fun_xph(i,j,k):
    if isInside(i,j,k):
        return f_upup_fun_pf(XPHtoPF((i,j,k))[0], XPHtoPF((i,j,k))[1], XPHtoPF((i,j,k))[2])
    
    else:
        return  K_upup_xph(i) +  P_upup_xph(i,j)+ P_upup_xph(i,k)+K_upup_ph(XPHtoPH((i,j,k))[0]) + K_upup_pp(XPHtoPP((i,j,k))[0])+P_upup_pp(XPHtoPP((i,j,k))[0],XPHtoPP((i,j,k))[1])+P_upup_pp(XPHtoPP((i,j,k))[0],XPHtoPP((i,j,k))[2])+P_upup_ph(XPHtoPH((i,j,k))[0],XPHtoPH((i,j,k))[1])+P_upup_ph(XPHtoPH((i,j,k))[0],XPHtoPH((i,j,k))[2]) 

def f_updo_fun_ph(i,j,k):
    if isInside(i,j,k):    
        return f_updo_fun_pf(PHtoPF((i,j,k))[0], PHtoPF((i,j,k))[1], PHtoPF((i,j,k))[2])
    else:
        return -U+ K_updo_ph(i) + P_updo_ph(i,j)+ P_updo_ph(i,k)+ K_updo_xph(PHtoXPH((i,j,k))[0]) + K_updo_pp(PHtoPP((i,j,k))[0])+P_updo_xph(PHtoXPH((i,j,k))[0],PHtoXPH((i,j,k))[1])+P_updo_xph(PHtoXPH((i,j,k))[0],PHtoXPH((i,j,k))[2])+P_updo_pp(PHtoPP((i,j,k))[0],PHtoPP((i,j,k))[1])+P_updo_pp(PHtoPP((i,j,k))[0],PHtoPP((i,j,k))[2]) 

def f_updo_fun_pp(i,j,k):
    if isInside(i,j,k):
        return f_updo_fun_pf(PPtoPF((i,j,k))[0], PPtoPF((i,j,k))[1], PPtoPF((i,j,k))[2])
    else:
        return  -U + K_updo_pp(i)+ P_updo_pp(i,j)+ P_updo_pp(i,k)+ K_updo_ph(PPtoPH((i,j,k))[0]) + K_updo_xph(PPtoXPH((i,j,k))[0])+P_updo_xph(PPtoXPH((i,j,k))[0],PPtoXPH((i,j,k))[1])+P_updo_xph(PPtoXPH((i,j,k))[0],PPtoXPH((i,j,k))[2])+P_updo_ph(PPtoPH((i,j,k))[0],PPtoPH((i,j,k))[1])+P_updo_ph(PPtoPH((i,j,k))[0],PPtoPH((i,j,k))[2]) 


def f_updo_fun_xph(i,j,k):
    if isInside(i,j,k):
        return f_updo_fun_pf(XPHtoPF((i,j,k))[0], XPHtoPF((i,j,k))[1], XPHtoPF((i,j,k))[2])
    
    else:
        return  -U + K_updo_xph(i) +  P_updo_xph(i,j)+ P_updo_xph(i,k)+K_updo_ph(XPHtoPH((i,j,k))[0]) + K_updo_pp(XPHtoPP((i,j,k))[0])+P_updo_pp(XPHtoPP((i,j,k))[0],XPHtoPP((i,j,k))[1])+P_updo_pp(XPHtoPP((i,j,k))[0],XPHtoPP((i,j,k))[2])+P_updo_ph(XPHtoPH((i,j,k))[0],XPHtoPH((i,j,k))[1])+P_updo_ph(XPHtoPH((i,j,k))[0],XPHtoPH((i,j,k))[2]) 


#GF

def G(wf):
    if (wf >= -N_fermi_gf and wf < N_fermi_gf):
        return re_g_iw[wf+N_fermi_gf]+1j*im_g_iw[wf+N_fermi_gf]
    else:
        wMat = 1j*(2*wf +1)*pi/beta     #GF asymtotic behavior
        return 1.0/wMat


# -------------------------------- Diagrammatic calculation of the susceptibility--------------------------------------------------
#                                                                                    
#------------------------------------------------------------------------------------------------------------------------------
#
#                                               PH
#                                                                                    
#
#------------------------------------------------------------------------------------------------------------------------------
print "PH channel diagrams"
#-------------------------------Susceptibilities in the ph notation----------------------------------------
    
def int_chi_0_ph(i,k):
    return -G(k-myfloor_div2(i))*G(k+myceil_div2(i))*1.0/beta   # (-) for internal loop

def int_chi_upup_ph(i,j,k):
    return G(j-myfloor_div2(i))*G(j+myceil_div2(i))*G(k-myfloor_div2(i))*G(k+myceil_div2(i))*f_upup_fun_ph(i,j,k)*1.0/beta*1.0/beta #(-)^2 two internal loops
    
def int_chi_updo_ph(i,j,k):
    return G(j-myfloor_div2(i))*G(j+myceil_div2(i))*G(k-myfloor_div2(i))*G(k+myceil_div2(i))*f_updo_fun_ph(i,j,k)*1.0/beta*1.0/beta 
    
# Arrays creation 

print "chi_updo_ph:"
chi_updo_ph_arr = np.array([dblsum_fit(partial(int_chi_updo_ph,i)) for i in range(-N_bose_suscept,N_bose_suscept+1)])
    
#-------------------------- Susceptibility function arrays PH-------------------------------------------
print "Susceptibiliy pp channel:"

Suscept_updo_ph = chi_updo_ph_arr
    
#------------------------------------------------------------------------------------------------------------------------------
#
#                                               PP
#                                                                                    
#
#------------------------------------------------------------------------------------------------------------------------------
print "PP channel diagrams"
#-------------------------------Susceptibilities in the pp notation----------------------------------------
    
def int_chi_0_pp(i,k):
    return G(myfloor_div2(i)-k-1)*G(k+myceil_div2(i))*1.0/beta  # 2 for the spin, 1/2 for the equivalent lines
    
def int_chi_updo_pp(i,j,k):
    return G(myfloor_div2(i)-j-1)*G(j+myceil_div2(i))*G(myfloor_div2(i)-k-1)*G(k+myceil_div2(i))*f_updo_fun_pp(i,j,k)*1.0/beta*1.0/beta # 2 for the spin, 1/4 for 2 pairs of equivalent lines, 2 from the diagram where we have f_xupdo_pp
    
    
# Arrays creation 

print "chi_0_pp:"
chi_0_pp_arr = np.array([sum_fit(partial(int_chi_0_pp,i)) for i in range(-N_bose_suscept,N_bose_suscept+1)])

print "chi_updo_pp:"
chi_updo_pp_arr = np.array([dblsum_fit(partial(int_chi_updo_pp,i)) for i in range(-N_bose_suscept,N_bose_suscept+1)])

#-------------------------Susceptibility in the pp channel; --------------------------

print "Susceptibiliy pp channel:"

Suscept_updo_pp = (chi_0_pp_arr+chi_updo_pp_arr)

#------------------------------------------------------------------------------------------------------------------------------
#
#                                               XPH
#                                                                                    
#
#------------------------------------------------------------------------------------------------------------------------------
print "XPH channel diagrams"
    
#-------------------------------Susceptibilities in the xph notation----------------------------------------
    
def int_chi_0_xph(i,k):
    return G(k-myfloor_div2(i))*G(k+myceil_div2(i))*1.0/beta #(-) from loop * (-) from swapping-> UPUP ; no loops no swapping-> UPDO 

def int_chi_updo_xph(i,j,k):
    return G(j-myfloor_div2(i))*G(j+myceil_div2(i))*G(k-myfloor_div2(i))*G(k+myceil_div2(i))*f_updo_fun_xph(i,k,j)*1.0/(beta*beta) # no loops no swapping 
    
# Arrays creation 

print "chi_0_xph:"
chi_0_xph_arr = np.array([sum_fit(partial(int_chi_0_xph,i)) for i in range(-N_bose_suscept,N_bose_suscept+1)])

print "chi_updo_xph:"
chi_updo_xph_arr = np.array([dblsum_fit(partial(int_chi_updo_xph,i)) for i in range(-N_bose_suscept,N_bose_suscept+1)])
    
#-------------------------Definition Susceptibility in XPH--------------------------

print "Susceptibiliy xph channel:"

Suscept_updo_xph = (chi_0_xph_arr + chi_updo_xph_arr)

#--------------------------STORING IN THE hdf5 FILE----------------------------------

# HDF5 NEW GROUP CREATION

if ('suscept_fun' in f):
    del f['suscept_fun']

#e HDF5 NEW SUBGROUP CREATION

suscept_grp = f.require_group("suscept_fun")

# HDF5 NEW DATASET CREATION

suscept_grp.create_dataset("RE_PH", data= Suscept_updo_ph.real, dtype='float64', compression="gzip", compression_opts=4)
suscept_grp.create_dataset("IM_PH", data= Suscept_updo_ph.imag, dtype='float64', compression="gzip", compression_opts=4)

suscept_grp.create_dataset("RE_PP", data= Suscept_updo_pp.real, dtype='float64', compression="gzip", compression_opts=4)
suscept_grp.create_dataset("IM_PP", data= Suscept_updo_pp.imag, dtype='float64', compression="gzip", compression_opts=4)
suscept_grp.create_dataset("RE_XPH", data= Suscept_updo_xph.real, dtype='float64', compression="gzip", compression_opts=4)
suscept_grp.create_dataset("IM_XPH", data= Suscept_updo_xph.imag, dtype='float64', compression="gzip", compression_opts=4)

f.close()
