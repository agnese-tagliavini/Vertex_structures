#!/usr/bin/python

#------------------------------------------------------------------------------------
#
# CONTENTS:
# This python script provides all the post processing operations from the (ED-> from  Vienna's group program) two-particle GF (2PGF) to other vertex objects:
#
# - Generalized Susceptibility -> Chi
# - Full vertex -> F
# 
# The input 2PGF are given in (PH, PP) mixed (bosonic, fermionic) shifted notation 
# in order to get a centralized vertex structure 
# 
# NOTE: One can easly translate one notation to another (either purely fermionic or mixed)
#       just using the functions defined in the library: agneselib (~/usr/include/agneselib)
#
# All the quantities are stored in a HDF5 file 
# 
# 
# The fermionic-bosonic coupling (The so-called kernel 2 asymptotic function K2) is also exactly 
# calculated from the ED program and here converted in the HDF5 file
#
# The kernel 1 function (connected with the physical susceptibility) is also provided from the 
# ED program and stored in the HDF5 file
#
# OUTLOOK:
#
# The HDF5 file serves as an input for the invesion of the Bethe-Salpeter equations and the
# calculation of the 2P irreducible vertices as well as the fully irreducible vertex. The program to perform this calculation is a C++ program (the main file is in src/ode.cpp) of the repository.
#
#------------------------------------------------------------------------------------
#------------------------------------IMPORTS------------------------------------


import numpy as np
import matplotlib.pyplot as pl
import os
import sys
import subprocess
import cmath
import math
import h5py
from os.path import isfile, join
import matplotlib.cm as cm
from numpy.linalg import inv
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from agneselibrary.mymath import *                          # mylibrary check in ~/usr/include/agneselib
from _functools import partial

#---------------------------------------------------------------------------------
#
#           PARAMETER DEFINITION
#
#---------------------------------------------------------------------------------

U=1.0           # U
MU= 0.0         # half-filling
BETA=50.0       # BETA

FFREQ_ED = 105      #positive fermionic frequencies in the mixed notation form the ED program
BFREQ_ED = 90       #positive bosonic frequency transfer in the mixed notation from the ED program,

FFREQ_SHIFT = 60   # positive fermionic frequencies in the mixed notation SHIFTED  
BFREQ_SHIFT = 90   # positive bosonic frequency transfer in the mixed notation SHIFTED

PATCH_COUNT = 1    # (DMFT) LOCAL OBJECTS

QN_COUNT = 1       # SINGLE-ORBITAL (SU2 SYMMETRIC CASE -> ONLY UPDO CONFIGURATION NEEDED )

pi = math.pi       

#---------------------------------------------------------------------------------

def run(command):
        output = subprocess.check_output(command, shell=True)
        return output

#----------------------------------------Create HDF5 files-----------------------------------------

if ('/home/agnese/Coding/Vertex_structures/dat/H5FILES/BETA50/4SITES/U1/PREPROC/U'+ str(U)+'_BETA_'+ str(BETA)+'_FFREQ_'+ str(FFREQ_SHIFT)+'_BFREQ_' + str(BFREQ_SHIFT)+'.h5'):
    os.system('rm -r /home/agnese/Coding/Vertex_structures/dat/H5FILES/BETA50/4SITES/U1/PREPROC/U_'+ str(U)+'_BETA_'+ str(BETA)+'_FFREQ_'+ str(FFREQ_SHIFT)+'_BFREQ_' + str(BFREQ_SHIFT)+'.h5')

f = h5py.File('/home/agnese/Coding/Vertex_structures/dat/H5FILES/BETA50/4SITES/U1/PREPROC/U_'+ str(U)+'_BETA_'+ str(BETA)+'_FFREQ_'+ str(FFREQ_SHIFT)+'_BFREQ_' + str(BFREQ_SHIFT)+'.h5', 'w')   # Create the datafile.h5

#---------------------------------------------------------------------------------
#
#                       GREEN'S FUNCTION AND SELF-ENERGY
#
#---------------------------------------------------------------------------------

# GF

g_iw  = np.loadtxt('/home/agnese/Coding/Vertex_structures/dat/ED/4SITES/U_'+str(U)+'_BETA_'+str(BETA)+'_FFREQ_'+ str(FFREQ_ED)+'_BFREQ_'+ str(BFREQ_ED)+'/gm_wim')


N_fermi_gf = g_iw.shape[0]
print ("Number of fermionic frequencies for the GF:")
print N_fermi_gf

def G(w):                           # imaginary part of the GF
 if (w >= 0):
     return g_iw[w,2]+1j*g_iw[w,3]
 else:
     return g_iw[-w-1,2]-1j*g_iw[-w-1,3]

# SELF-ENERGY

def Sig(w):
    return 1j*(2.0*w+1.0)*pi/BETA-MU-1./G(w) # Impurity SE 

g_arr_re = np.array([[[[G(n).real for i in range (0,QN_COUNT)] for j in range (0,QN_COUNT)] for k in range (0,PATCH_COUNT)] for n in range (-N_fermi_gf,N_fermi_gf)])
g_arr_im = np.array([[[[G(n).imag for i in range (0,QN_COUNT)] for j in range (0,QN_COUNT)] for k in range (0,PATCH_COUNT)] for n in range (-N_fermi_gf,N_fermi_gf)])

se_arr_re = np.array([[[[Sig(n).real for i in range (0,QN_COUNT)] for j in range (0,QN_COUNT)] for k in range (0,PATCH_COUNT)]for n in range (-N_fermi_gf,N_fermi_gf)])
se_arr_im = np.array([[[[Sig(n).imag for i in range (0,QN_COUNT)] for j in range (0,QN_COUNT)] for k in range (0,PATCH_COUNT)] for n in range (-N_fermi_gf,N_fermi_gf)])
fgrid_arr = np.array([(2*m +1)*pi/BETA for m in range (-N_fermi_gf,N_fermi_gf)])

#----------------------------------------Create HDF5 files-----------------------------------------

f.create_dataset('Giw/RE', data=g_arr_re, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('Giw/IM', data=g_arr_im, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('Giw/fgrid', data=fgrid_arr, dtype='float64', compression="gzip", compression_opts=4)

f.create_dataset('Sig/RE', data=se_arr_re, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('Sig/IM', data=se_arr_im, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('Sig/fgrid', data=fgrid_arr, dtype='float64', compression="gzip", compression_opts=4)


#--------------------------------------------------------------------------------------
#
#                               TWO-PARTICLE OBJECTS
#
#---------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------
#
#                               THREE-FREQUENCY OBJECTS
#
#--------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
#
#                                       PP
#
#---------------------------------------------------------------------------------------

#-----------------------------GENERALIZED SUSCEPTIBILITY PP-----------------------------

genchi_pp = np.loadtxt('/home/agnese/Coding/Vertex_structures/dat/ED/4SITES/U_'+str(U)+'_BETA_'+str(BETA)+'_FFREQ_'+ str(FFREQ_ED)+'_BFREQ_'+ str(BFREQ_ED)+'/vert_chi_pp')

ffreq_pp = FFREQ_ED
print ffreq_pp
bfreq_pp = BFREQ_ED
print bfreq_pp


def re_genchi_updo_pp(wb,wf,wf1):
	return 	genchi_pp[(2*ffreq_pp)*(2*ffreq_pp)*(wb+bfreq_pp)+(2*ffreq_pp)*(wf+ffreq_pp)+ (wf1+ffreq_pp), 5]

def re_genchi_upup_pp(wb,wf,wf1):
	return 	genchi_pp[(2*ffreq_pp)*(2*ffreq_pp)*(wb+bfreq_pp)+(2*ffreq_pp)*(wf+ffreq_pp)+ (wf1+ffreq_pp), 3]

def im_genchi_updo_pp(wb,wf,wf1):
	return 	genchi_pp[(2*ffreq_pp)*(2*ffreq_pp)*(wb+bfreq_pp)+(2*ffreq_pp)*(wf+ffreq_pp)+ (wf1+ffreq_pp), 6]

def im_genchi_upup_pp(wb,wf,wf1):
	return 	genchi_pp[(2*ffreq_pp)*(2*ffreq_pp)*(wb+bfreq_pp)+(2*ffreq_pp)*(wf+ffreq_pp)+ (wf1+ffreq_pp), 4]

# ----------define arrays to store in hdf5 file -> FROM HERE OBJECT IN THE SHIFTED NOTATION!!

re_genchi_arr_upup_pp = np.array([[[[[[[[[[re_genchi_upup_pp(m,n+myceil_div2(m),p+myceil_div2(m)) for i in range (0,QN_COUNT)] for j in range (0,QN_COUNT)]for ip in range(0,QN_COUNT)] for jp in range (0,QN_COUNT)] for q in range (0,PATCH_COUNT)] for k in range (0,PATCH_COUNT)] for kp in range (0,PATCH_COUNT)] for p in range (-FFREQ_SHIFT,FFREQ_SHIFT )] for n in range (-FFREQ_SHIFT,FFREQ_SHIFT)] for m in range (-BFREQ_SHIFT,BFREQ_SHIFT+1)] )

im_genchi_arr_upup_pp = np.array([[[[[[[[[[im_genchi_upup_pp(m,n+myceil_div2(m),p+myceil_div2(m)) for i in range (0,QN_COUNT)] for j in range (0,QN_COUNT)]for ip in range (0,QN_COUNT)] for jp in range (0,QN_COUNT)] for q in range (0,PATCH_COUNT)] for k in range (0,PATCH_COUNT)] for kp in range (0,PATCH_COUNT)] for p in range (-FFREQ_SHIFT,FFREQ_SHIFT )] for n in range (-FFREQ_SHIFT,FFREQ_SHIFT)] for m in range (-BFREQ_SHIFT,BFREQ_SHIFT+1)] )

re_genchi_arr_updo_pp = np.array([[[[[[[[[[re_genchi_updo_pp(m,n+myceil_div2(m),p+myceil_div2(m)) for i in range (0,QN_COUNT)] for j in range (0,QN_COUNT)] for ip in range (0,QN_COUNT)] for jp in range (0,QN_COUNT)] for q in range (0,PATCH_COUNT)] for k in range (0,PATCH_COUNT)] for kp in range (0,PATCH_COUNT)] for p in range (-FFREQ_SHIFT,FFREQ_SHIFT )] for n in range (-FFREQ_SHIFT,FFREQ_SHIFT)] for m in range (-BFREQ_SHIFT,BFREQ_SHIFT+1)] )

im_genchi_arr_updo_pp = np.array([[[[[[[[[[im_genchi_updo_pp(m,n+myceil_div2(m),p+myceil_div2(m)) for i in range (0,QN_COUNT)] for j in range (0,QN_COUNT)]for ip in range (0,QN_COUNT)] for jp in range (0,QN_COUNT)] for q in range (0,PATCH_COUNT)] for k in range (0,PATCH_COUNT)] for kp in range (0,PATCH_COUNT)] for p in range (-FFREQ_SHIFT,FFREQ_SHIFT )] for n in range (-FFREQ_SHIFT,FFREQ_SHIFT)] for m in range (-BFREQ_SHIFT,BFREQ_SHIFT+1)] )

fgrid_arr_pp = np.array([(2*m +1)*pi/BETA for m in range (-FFREQ_SHIFT,FFREQ_SHIFT)])
bgrid_arr_pp = np.array([(2*n)*pi/BETA for n in range (-BFREQ_SHIFT,BFREQ_SHIFT+1)])

#----------------------------------------Create HDF5 files-----------------------------------------

f.create_dataset('GENCHI/PP/RE_GENCHI_UPUP', data=re_genchi_arr_upup_pp, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('GENCHI/PP/IM_GENCHI_UPUP', data=im_genchi_arr_upup_pp, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('GENCHI/PP/RE_GENCHI_UPDO', data=re_genchi_arr_updo_pp, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('GENCHI/PP/IM_GENCHI_UPDO', data=im_genchi_arr_updo_pp, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('GENCHI/PP/fgrid', data=fgrid_arr_pp, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('GENCHI/PP/bgrid', data=bgrid_arr_pp, dtype='float64', compression="gzip", compression_opts=4)


#----------------------------------- FULL VERTEX PP -----------------------------------------------------------

def genchi_updo_shift_pp(wb,wf,wf1):
	return 	re_genchi_updo_pp(wb,wf+myceil_div2(wb),wf1+myceil_div2(wb))+1.j*im_genchi_updo_pp(wb,wf+myceil_div2(wb),wf1+myceil_div2(wb))

def genchi_upup_shift_pp(wb,wf,wf1):
	return 	re_genchi_upup_pp(wb,wf+myceil_div2(wb),wf1+myceil_div2(wb))+1.j*im_genchi_upup_pp(wb,wf+myceil_div2(wb),wf1+myceil_div2(wb))


def chi_l_pp(wb,wf,wf1):  # Legs on one side of the diagram
    return BETA*G(wf+myceil_div2(wb))*G(myfloor_div2(wb)-wf1-1)

def chi_0_pp(wb,wf,wf1):
    if (wf == wf1):
        return BETA*G(wf+myceil_div2(wb))*G(myfloor_div2(wb)-wf-1)
    else:
        return 0.0

def f_upup_pp(wb,wf,wf1):
    return BETA*BETA*(1.0/(chi_l_pp(wb,wf,wf)))*(genchi_upup_shift_pp(wb,wf,wf1)+chi_0_pp(wb,wf,wf1))*(1.0/(chi_l_pp(wb,wf1,wf1)))

def f_updo_pp(wb,wf,wf1):
    return BETA*BETA*(1.0/(chi_l_pp(wb,wf,wf)))*(genchi_updo_shift_pp(wb,wf,wf1))*(1.0/(chi_l_pp(wb,wf1,wf1)))

# ----------define arrays to store in hdf5 file

re_f_arr_upup_pp = np.array([[[[[[[[[[f_upup_pp(m,n,p).real for i in range (0,QN_COUNT)] for j in range (0,QN_COUNT)]for ip in range(0,QN_COUNT)] for jp in range (0,QN_COUNT)] for q in range (0,PATCH_COUNT)] for k in range (0,PATCH_COUNT)] for kp in range (0,PATCH_COUNT)] for p in range (-FFREQ_SHIFT,FFREQ_SHIFT )] for n in range (-FFREQ_SHIFT,FFREQ_SHIFT)] for m in range (-BFREQ_SHIFT,BFREQ_SHIFT+1)] )

im_f_arr_upup_pp = np.array([[[[[[[[[[f_upup_pp(m,n,p).imag for i in range (0,QN_COUNT)] for j in range (0,QN_COUNT)]for ip in range(0,QN_COUNT)] for jp in range (0,QN_COUNT)] for q in range (0,PATCH_COUNT)] for k in range (0,PATCH_COUNT)] for kp in range (0,PATCH_COUNT)] for p in range (-FFREQ_SHIFT,FFREQ_SHIFT )] for n in range (-FFREQ_SHIFT,FFREQ_SHIFT)] for m in range (-BFREQ_SHIFT,BFREQ_SHIFT+1)] )

re_f_arr_updo_pp = np.array([[[[[[[[[[f_updo_pp(m,n,p).real for i in range (0,QN_COUNT)] for j in range (0,QN_COUNT)]for ip in range(0,QN_COUNT)] for jp in range (0,QN_COUNT)] for q in range (0,PATCH_COUNT)] for k in range (0,PATCH_COUNT)] for kp in range (0,PATCH_COUNT)] for p in range (-FFREQ_SHIFT,FFREQ_SHIFT )] for n in range (-FFREQ_SHIFT,FFREQ_SHIFT)] for m in range (-BFREQ_SHIFT,BFREQ_SHIFT+1)] )

im_f_arr_updo_pp = np.array([[[[[[[[[[f_updo_pp(m,n,p).imag for i in range (0,QN_COUNT)] for j in range (0,QN_COUNT)]for ip in range(0,QN_COUNT)] for jp in range (0,QN_COUNT)] for q in range (0,PATCH_COUNT)] for k in range (0,PATCH_COUNT)] for kp in range (0,PATCH_COUNT)] for p in range (-FFREQ_SHIFT,FFREQ_SHIFT )] for n in range (-FFREQ_SHIFT,FFREQ_SHIFT)] for m in range (-BFREQ_SHIFT,BFREQ_SHIFT+1)] )

#----------------------------------------Create HDF5 files-----------------------------------------

f.create_dataset('VERT/PP/RE_F_UPUP', data=re_f_arr_upup_pp, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('VERT/PP/IM_F_UPUP', data=im_f_arr_upup_pp, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('VERT/PP/RE_F_UPDO', data=re_f_arr_updo_pp, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('VERT/PP/IM_F_UPDO', data=im_f_arr_updo_pp, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('VERT/PP/fgrid', data=fgrid_arr_pp, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('VERT/PP/bgrid', data=bgrid_arr_pp, dtype='float64', compression="gzip", compression_opts=4)

#---------------------------------------------------------------------------------------
#
#                                       PH
#
#---------------------------------------------------------------------------------------

#-----------------------------GENERALIZED SUSCEPTIBILITY PH-----------------------------

genchi_ph = np.loadtxt('/home/agnese/Coding/Vertex_structures/dat/ED/4SITES/U_'+str(U)+'_BETA_'+str(BETA)+'_FFREQ_'+ str(FFREQ_ED)+'_BFREQ_'+ str(BFREQ_ED)+'/vert_chi')

ffreq_ph = FFREQ_ED
print ffreq_pp
bfreq_ph = BFREQ_ED
print bfreq_pp


def re_genchi_updo_ph(wb,wf,wf1):
	return 	genchi_ph[(2*ffreq_ph)*(2*ffreq_ph)*(wb+bfreq_ph)+(2*ffreq_ph)*(wf+ffreq_ph)+ (wf1+ffreq_ph), 5]

def re_genchi_upup_ph(wb,wf,wf1):
	return 	genchi_ph[(2*ffreq_ph)*(2*ffreq_ph)*(wb+bfreq_ph)+(2*ffreq_ph)*(wf+ffreq_ph)+ (wf1+ffreq_ph), 3]

def im_genchi_updo_ph(wb,wf,wf1):
	return 	genchi_ph[(2*ffreq_ph)*(2*ffreq_ph)*(wb+bfreq_ph)+(2*ffreq_ph)*(wf+ffreq_ph)+ (wf1+ffreq_ph), 6]

def im_genchi_upup_ph(wb,wf,wf1):
	return 	genchi_ph[(2*ffreq_ph)*(2*ffreq_ph)*(wb+bfreq_ph)+(2*ffreq_ph)*(wf+ffreq_ph)+ (wf1+ffreq_ph), 4]

# ----------define arrays to store in hdf5 file -> FROM HERE OBJECT IN THE SHIFTED NOTATION!!

re_genchi_arr_upup_ph = np.array([[[[[[[[[[re_genchi_upup_ph(m,n-myfloor_div2(m),p-myfloor_div2(m)) for i in range (0,QN_COUNT)] for j in range (0,QN_COUNT)]for ip in range(0,QN_COUNT)] for jp in range (0,QN_COUNT)] for q in range (0,PATCH_COUNT)] for k in range (0,PATCH_COUNT)] for kp in range (0,PATCH_COUNT)] for p in range (-FFREQ_SHIFT,FFREQ_SHIFT )] for n in range (-FFREQ_SHIFT,FFREQ_SHIFT)] for m in range (-BFREQ_SHIFT,BFREQ_SHIFT+1)] )

im_genchi_arr_upup_ph = np.array([[[[[[[[[[im_genchi_upup_ph(m,n-myfloor_div2(m),p-myfloor_div2(m)) for i in range (0,QN_COUNT)] for j in range (0,QN_COUNT)]for ip in range (0,QN_COUNT)] for jp in range (0,QN_COUNT)] for q in range (0,PATCH_COUNT)] for k in range (0,PATCH_COUNT)] for kp in range (0,PATCH_COUNT)] for p in range (-FFREQ_SHIFT,FFREQ_SHIFT )] for n in range (-FFREQ_SHIFT,FFREQ_SHIFT)] for m in range (-BFREQ_SHIFT,BFREQ_SHIFT+1)] )

re_genchi_arr_updo_ph = np.array([[[[[[[[[[re_genchi_updo_ph(m,n-myfloor_div2(m),p-myfloor_div2(m)) for i in range (0,QN_COUNT)] for j in range (0,QN_COUNT)] for ip in range (0,QN_COUNT)] for jp in range (0,QN_COUNT)] for q in range (0,PATCH_COUNT)] for k in range (0,PATCH_COUNT)] for kp in range (0,PATCH_COUNT)] for p in range (-FFREQ_SHIFT,FFREQ_SHIFT )] for n in range (-FFREQ_SHIFT,FFREQ_SHIFT)] for m in range (-BFREQ_SHIFT,BFREQ_SHIFT+1)] )

im_genchi_arr_updo_ph = np.array([[[[[[[[[[im_genchi_updo_ph(m,n-myfloor_div2(m),p-myfloor_div2(m)) for i in range (0,QN_COUNT)] for j in range (0,QN_COUNT)]for ip in range (0,QN_COUNT)] for jp in range (0,QN_COUNT)] for q in range (0,PATCH_COUNT)] for k in range (0,PATCH_COUNT)] for kp in range (0,PATCH_COUNT)] for p in range (-FFREQ_SHIFT,FFREQ_SHIFT )] for n in range (-FFREQ_SHIFT,FFREQ_SHIFT)] for m in range (-BFREQ_SHIFT,BFREQ_SHIFT+1)] )

fgrid_arr_ph = np.array([(2*m +1)*pi/BETA for m in range (-FFREQ_SHIFT,FFREQ_SHIFT)])
bgrid_arr_ph = np.array([(2*n)*pi/BETA for n in range (-BFREQ_SHIFT,BFREQ_SHIFT+1)])

#----------------------------------------Create HDF5 files-----------------------------------------

f.create_dataset('GENCHI/PH/RE_GENCHI_UPUP', data=re_genchi_arr_upup_ph, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('GENCHI/PH/IM_GENCHI_UPUP', data=im_genchi_arr_upup_ph, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('GENCHI/PH/RE_GENCHI_UPDO', data=re_genchi_arr_updo_ph, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('GENCHI/PH/IM_GENCHI_UPDO', data=im_genchi_arr_updo_ph, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('GENCHI/PH/fgrid', data=fgrid_arr_ph, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('GENCHI/PH/bgrid', data=bgrid_arr_ph, dtype='float64', compression="gzip", compression_opts=4)

f.create_dataset('GENCHI/XPH/RE_GENCHI_UPUP', data=-re_genchi_arr_upup_xph, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('GENCHI/XPH/IM_GENCHI_UPUP', data=-im_genchi_arr_upup_xph, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('GENCHI/XPH/RE_GENCHI_UPDO', data=re_genchi_arr_updo_xph-re_genchi_arr_upup_xph, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('GENCHI/XPH/IM_GENCHI_UPDO', data=im_genchi_arr_updo_xph-im_genchi_arr_upup_xph, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('GENCHI/XPH/fgrid', data=fgrid_arr_ph, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('GENCHI/XPH/bgrid', data=bgrid_arr_ph, dtype='float64', compression="gzip", compression_opts=4)

#----------------------------------- FULL VERTEX PH -----------------------------------------------------------

def genchi_updo_shift_ph(wb,wf,wf1):
	return 	re_genchi_updo_ph(wb,wf-myfloor_div2(wb),wf1-myfloor_div2(wb))+1.j*im_genchi_updo_ph(wb,wf-myfloor_div2(wb),wf1-myfloor_div2(wb))

def genchi_upup_shift_ph(wb,wf,wf1):
	return 	re_genchi_upup_ph(wb,wf-myfloor_div2(wb),wf1-myfloor_div2(wb))+1.j*im_genchi_upup_ph(wb,wf-myfloor_div2(wb),wf1-myfloor_div2(wb))


def chi_lup_ph(wb,wf,wf1):
    return BETA*G(wf-myfloor_div2(wb))*G(wf1+myceil_div2(wb))

def chi_0_ph(wb,wf,wf1):
    if (wf == wf1):
        return BETA*G(wf+myceil_div2(wb))*G(wf-myfloor_div2(wb))
    else:
        return 0.0


def f_upup_ph(wb,wf,wf1):
    return BETA*BETA*(1.0/chi_lup_ph(wb,wf,wf))*(genchi_upup_shift_ph(wb,wf,wf1)+chi_0_ph(wb,wf,wf1))*(1.0/chi_lup_ph(wb,wf1,wf1))

def f_updo_ph(wb,wf,wf1):
    return BETA*BETA*(1.0/chi_lup_ph(wb,wf,wf))*(genchi_updo_shift_ph(wb,wf,wf1))*(1.0/chi_lup_ph(wb,wf1,wf1))

# ----------define arrays to store in hdf5 file

re_f_arr_upup_ph = np.array([[[[[[[[[[f_upup_ph(m,n,p).real for i in range (0,QN_COUNT)] for j in range (0,QN_COUNT)]for ip in range(0,QN_COUNT)] for jp in range (0,QN_COUNT)] for q in range (0,PATCH_COUNT)] for k in range (0,PATCH_COUNT)] for kp in range (0,PATCH_COUNT)] for p in range (-FFREQ_SHIFT,FFREQ_SHIFT )] for n in range (-FFREQ_SHIFT,FFREQ_SHIFT)] for m in range (-BFREQ_SHIFT,BFREQ_SHIFT+1)] )

im_f_arr_upup_ph = np.array([[[[[[[[[[f_upup_ph(m,n,p).imag for i in range (0,QN_COUNT)] for j in range (0,QN_COUNT)]for ip in range(0,QN_COUNT)] for jp in range (0,QN_COUNT)] for q in range (0,PATCH_COUNT)] for k in range (0,PATCH_COUNT)] for kp in range (0,PATCH_COUNT)] for p in range (-FFREQ_SHIFT,FFREQ_SHIFT )] for n in range (-FFREQ_SHIFT,FFREQ_SHIFT)] for m in range (-BFREQ_SHIFT,BFREQ_SHIFT+1)] )

re_f_arr_updo_ph = np.array([[[[[[[[[[f_updo_ph(m,n,p).real for i in range (0,QN_COUNT)] for j in range (0,QN_COUNT)]for ip in range(0,QN_COUNT)] for jp in range (0,QN_COUNT)] for q in range (0,PATCH_COUNT)] for k in range (0,PATCH_COUNT)] for kp in range (0,PATCH_COUNT)] for p in range (-FFREQ_SHIFT,FFREQ_SHIFT )] for n in range (-FFREQ_SHIFT,FFREQ_SHIFT)] for m in range (-BFREQ_SHIFT,BFREQ_SHIFT+1)] )

im_f_arr_updo_ph = np.array([[[[[[[[[[f_updo_ph(m,n,p).imag for i in range (0,QN_COUNT)] for j in range (0,QN_COUNT)]for ip in range(0,QN_COUNT)] for jp in range (0,QN_COUNT)] for q in range (0,PATCH_COUNT)] for k in range (0,PATCH_COUNT)] for kp in range (0,PATCH_COUNT)] for p in range (-FFREQ_SHIFT,FFREQ_SHIFT )] for n in range (-FFREQ_SHIFT,FFREQ_SHIFT)] for m in range (-BFREQ_SHIFT,BFREQ_SHIFT+1)] )

#----------------------------------------Create HDF5 files-----------------------------------------

f.create_dataset('VERT/PH/RE_F_UPUP', data=re_f_arr_upup_ph, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('VERT/PH/IM_F_UPUP', data=im_f_arr_upup_ph, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('VERT/PH/RE_F_UPDO', data=re_f_arr_updo_ph, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('VERT/PH/IM_F_UPDO', data=im_f_arr_updo_ph, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('VERT/PH/fgrid', data=fgrid_arr_ph, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('VERT/PH/bgrid', data=bgrid_arr_ph, dtype='float64', compression="gzip", compression_opts=4)

f.create_dataset('VERT/XPH/RE_F_UPUP', data=-re_f_arr_upup_ph, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('VERT/XPH/IM_F_UPUP', data=-im_f_arr_upup_ph, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('VERT/XPH/RE_F_UPDO', data=re_f_arr_updo_ph-re_f_arr_upup_ph, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('VERT/XPH/IM_F_UPDO', data=im_f_arr_updo_ph-im_f_arr_upup_ph, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('VERT/XPH/fgrid', data=fgrid_arr_ph, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('VERT/XPH/bgrid', data=bgrid_arr_ph, dtype='float64', compression="gzip", compression_opts=4)

#-------------------------------------------------------------------------------------
#
#                               ONE-FREQUENCY OBJECTS
#
#--------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
#
#                                       PP/PH/XPH
#
#---------------------------------------------------------------------------------------

#----------------------------- CHI PP PH XPH -----------------------------

chi = np.loadtxt('/home/agnese/Coding/Vertex_structures/dat/ED/4SITES/U_'+str(U)+'_BETA_'+str(BETA)+'_FFREQ_'+ str(FFREQ_ED)+'_BFREQ_'+ str(BFREQ_ED)+'/suscept')

bfreq_pp = CHI_BFREQ_ED
print bfreq_pp


def re_chi_updo_pp(wb):
    if(w >= 0):
        return chi[wb, 5]
    else:
        return -chi[wb, 5] 

def im_chi_updo_pp(wb):
    if(w >= 0):
        return chi[wb, 6]
    else:
        return -chi[wb, 6]

def re_chi_updo_ph(wb):
    if(w >= 0):
        return chi[wb, 3]
    else:
        return -chi[wb, 3] 

def im_chi_updo_ph(wb):
    if(w >= 0):
        return chi[wb, 4]
    else:
        return -chi[wb, 4]

def re_chi_upup_ph(wb):
    if(w >= 0):
        return chi[wb, 1]
    else:
        return -chi[wb, 1] 

def im_chi_upup_ph(wb):
    if(w >= 0):
        return chi[wb, 2]
    else:
        return -chi[wb, 2]

# ----------define arrays to store in hdf5 file -> FROM HERE OBJECT IN THE SHIFTED NOTATION!!

re_chi_arr_upup_pp = np.array([[[[[[0 for i in range (0,QN_COUNT)] for j in range (0,QN_COUNT)]for ip in range(0,QN_COUNT)] for jp in range (0,QN_COUNT)] for q in range (0,PATCH_COUNT)] for m in range (-TRI_BFREQ_SHIFT,TRI_BFREQ_SHIFT+1)] )

im_chi_arr_upup_pp = np.array([[[[[[0 for i in range (0,QN_COUNT)] for j in range (0,QN_COUNT)]for ip in range (0,QN_COUNT)] for jp in range (0,QN_COUNT)] for q in range (0,PATCH_COUNT)]for m in range (-TRI_BFREQ_SHIFT,TRI_BFREQ_SHIFT+1)] )

re_chi_arr_updo_pp = np.array([[[[[[re_chi_updo_pp(m) for i in range (0,QN_COUNT)] for j in range (0,QN_COUNT)] for ip in range (0,QN_COUNT)] for jp in range (0,QN_COUNT)] for q in range (0,PATCH_COUNT)] for m in range (-TRI_BFREQ_SHIFT,TRI_BFREQ_SHIFT+1)] )

im_chi_arr_updo_pp = np.array([[[[[[im_chi_updo_pp(m) for i in range (0,QN_COUNT)] for j in range (0,QN_COUNT)]for ip in range (0,QN_COUNT)] for jp in range (0,QN_COUNT)] for q in range (0,PATCH_COUNT)] for m in range (-TRI_BFREQ_SHIFT,TRI_BFREQ_SHIFT+1)] )

bgrid_arr_pp = np.array([(2*n)*pi/BETA for n in range (-CHI_BFREQ_SHIFT,CHI_BFREQ_SHIFT+1)])

# ----------define arrays to store in hdf5 file -> FROM HERE OBJECT IN THE SHIFTED NOTATION!!

re_chi_arr_upup_ph = np.array([[[[[[re_chi_upup_ph(m) for i in range (0,QN_COUNT)] for j in range (0,QN_COUNT)]for ip in range(0,QN_COUNT)] for jp in range (0,QN_COUNT)] for q in range (0,PATCH_COUNT)] for m in range (-TRI_BFREQ_SHIFT,TRI_BFREQ_SHIFT+1)] )

im_chi_arr_upup_ph = np.array([[[[[[im_chi_upup_ph(m) for i in range (0,QN_COUNT)] for j in range (0,QN_COUNT)]for ip in range (0,QN_COUNT)] for jp in range (0,QN_COUNT)] for q in range (0,PATCH_COUNT)]for m in range (-TRI_BFREQ_SHIFT,TRI_BFREQ_SHIFT+1)] )

re_chi_arr_updo_ph = np.array([[[[[[re_chi_updo_ph(m) for i in range (0,QN_COUNT)] for j in range (0,QN_COUNT)] for ip in range (0,QN_COUNT)] for jp in range (0,QN_COUNT)] for q in range (0,PATCH_COUNT)] for m in range (-TRI_BFREQ_SHIFT,TRI_BFREQ_SHIFT+1)] )

im_chi_arr_updo_ph = np.array([[[[[[im_chi_updo_ph(m) for i in range (0,QN_COUNT)] for j in range (0,QN_COUNT)]for ip in range (0,QN_COUNT)] for jp in range (0,QN_COUNT)] for q in range (0,PATCH_COUNT)] for m in range (-TRI_BFREQ_SHIFT,TRI_BFREQ_SHIFT+1)] )

bgrid_arr_ph = np.array([(2*n)*pi/BETA for n in range (-CHI_BFREQ_SHIFT,CHI_BFREQ_SHIFT+1)])

#----------------------------------------Create HDF5 files-----------------------------------------

f.create_dataset('CHI/PP/RE_CHI_UPUP', data=re_chi_arr_upup_pp, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('CHI/PP/IM_CHI_UPUP', data=im_chi_arr_upup_pp, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('CHI/PP/RE_CHI_UPDO', data=re_chi_arr_updo_pp, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('CHI/PP/IM_CHI_UPDO', data=im_chi_arr_updo_pp, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('CHI/PP/bgrid', data=bgrid_arr_pp, dtype='float64', compression="gzip", compression_opts=4)

f.create_dataset('CHI/PH/RE_CHI_UPUP', data=re_chi_arr_upup_ph, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('CHI/PH/IM_CHI_UPUP', data=im_chi_arr_upup_ph, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('CHI/PH/RE_CHI_UPDO', data=re_chi_arr_updo_ph, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('CHI/PH/IM_CHI_UPDO', data=im_chi_arr_updo_ph, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('CHI/PH/bgrid', data=bgrid_arr_ph, dtype='float64', compression="gzip", compression_opts=4)

f.create_dataset('CHI/XPH/RE_CHI_UPUP', data=-re_chi_arr_upup_ph, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('CHI/XPH/IM_CHI_UPUP', data=-im_chi_arr_upup_ph, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('CHI/XPH/RE_CHI_UPDO', data=re_chi_arr_updo_ph-re_chi_arr_upup_ph, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('CHI/XPH/IM_CHI_UPDO', data=im_chi_arr_updo_ph-im_chi_arr_upup_ph, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('CHI/XPH/bgrid', data=bgrid_arr_ph, dtype='float64', compression="gzip", compression_opts=4)

#-------------------------------------------------------------------------------------
#
#                               TWO-FREQUENCY OBJECTS
#
#--------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
#
#                                       PP
#
#---------------------------------------------------------------------------------------

#----------------------------- TRILEG PP -----------------------------

trileg_pp = np.loadtxt('/home/agnese/Coding/Vertex_structures/dat/ED/4SITES/U_'+str(U)+'_BETA_'+str(BETA)+'_FFREQ_'+ str(FFREQ_ED)+'_BFREQ_'+ str(BFREQ_ED)+'/trileg_pp')

ffreq_pp = TRI_FFREQ_ED
print ffreq_pp
bfreq_pp = TRI_BFREQ_ED
print bfreq_pp


def re_trileg_updo_pp(wb,wf):
	return 	U*(trileg_pp[(2*ffreq_pp)*(wb+bfreq_pp) + (wf+ffreq_pp), 2] - 1) + U**2*re_chi_updo_pp(wb)

def re_trileg_upup_pp(wb,wf):
	return 	0.0

def im_trileg_updo_pp(wb,wf):
	return 	U*(trileg_pp[(2*ffreq_pp)*(wb+bfreq_pp)+ (wf+ffreq_pp), 3] - 1) + U**2*im_chi_updo _pp(wb) 

def im_trileg_upup_pp(wb,wf):
	return 	0.0

# ----------define arrays to store in hdf5 file -> FROM HERE OBJECT IN THE SHIFTED NOTATION!!

re_trileg_arr_upup_pp = np.array([[[[[[[[re_trileg_upup_pp(m,n+myceil_div2(m)) for i in range (0,QN_COUNT)] for j in range (0,QN_COUNT)]for ip in range(0,QN_COUNT)] for jp in range (0,QN_COUNT)] for q in range (0,PATCH_COUNT)] for k in range (0,PATCH_COUNT)] for n in range (-TRI_FFREQ_SHIFT,TRI_FFREQ_SHIFT)] for m in range (-TRI_BFREQ_SHIFT,TRI_BFREQ_SHIFT+1)] )

im_trileg_arr_upup_pp = np.array([[[[[[[[im_trileg_upup_pp(m,n+myceil_div2(m)) for i in range (0,QN_COUNT)] for j in range (0,QN_COUNT)]for ip in range (0,QN_COUNT)] for jp in range (0,QN_COUNT)] for q in range (0,PATCH_COUNT)] for k in range (0,PATCH_COUNT)] for n in range (-TRI_FFREQ_SHIFT,TRI_FFREQ_SHIFT)] for m in range (-TRI_BFREQ_SHIFT,TRI_BFREQ_SHIFT+1)] )

re_trileg_arr_updo_pp = np.array([[[[[[[[[[re_trileg_updo_pp(m,n+myceil_div2(m)) for i in range (0,QN_COUNT)] for j in range (0,QN_COUNT)] for ip in range (0,QN_COUNT)] for jp in range (0,QN_COUNT)] for q in range (0,PATCH_COUNT)] for k in range (0,PATCH_COUNT)] for n in range (-TRI_FFREQ_SHIFT,TRI_FFREQ_SHIFT)] for m in range (-TRI_BFREQ_SHIFT,TRI_BFREQ_SHIFT+1)] )

im_trileg_arr_updo_pp = np.array([[[[[[[[im_trileg_updo_pp(m,n+myceil_div2(m)) for i in range (0,QN_COUNT)] for j in range (0,QN_COUNT)]for ip in range (0,QN_COUNT)] for jp in range (0,QN_COUNT)] for q in range (0,PATCH_COUNT)] for k in range (0,PATCH_COUNT)] for n in range (-TRI_FFREQ_SHIFT,TRI_FFREQ_SHIFT)] for m in range (-TRI_BFREQ_SHIFT,TRI_BFREQ_SHIFT+1)] )

fgrid_arr_pp = np.array([(2*m +1)*pi/BETA for m in range (-TRI_FFREQ_SHIFT,TRI_FFREQ_SHIFT)])
bgrid_arr_pp = np.array([(2*n)*pi/BETA for n in range (-TRI_BFREQ_SHIFT,TRI_BFREQ_SHIFT+1)])

#----------------------------------------Create HDF5 files-----------------------------------------

f.create_dataset('TRILEG/PP/RE_TRILEG_UPUP', data=re_trileg_arr_upup_pp, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('TRILEG/PP/IM_TRILEG_UPUP', data=im_trileg_arr_upup_pp, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('TRILEG/PP/RE_TRILEG_UPDO', data=re_trileg_arr_updo_pp, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('TRILEG/PP/IM_TRILEG_UPDO', data=im_trileg_arr_updo_pp, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('TRILEG/PP/fgrid', data=fgrid_arr_pp, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('TRILEG/PP/bgrid', data=bgrid_arr_pp, dtype='float64', compression="gzip", compression_opts=4)


#---------------------------------------------------------------------------------------
#
#                                       PH
#
#---------------------------------------------------------------------------------------

trileg_ph = np.loadtxt('/home/agnese/Coding/Vertex_structures/dat/ED/4SITES/U_'+str(U)+'_BETA_'+str(BETA)+'_FFREQ_'+ str(FFREQ_ED)+'_BFREQ_'+ str(BFREQ_ED)+'/trileg')

ffreq_ph = TRI_FFREQ_ED
print ffreq_ph
bfreq_ph = TRI_BFREQ_ED
print bfreq_ph


def re_trileg_updo_ph(wb,wf):
	return 	U*(trileg_ph[(2*ffreq_ph)*(wb+bfreq_ph) + (wf+ffreq_ph), 2]-trileg_ph[(2*ffreq_ph)*(wb+bfreq_ph) + (wf+ffreq_ph), 2] + 1) - U**2*re_chi_updo_ph(wb) 

def re_trileg_upup_ph(wb,wf):
	return 	U*(trileg_ph[(2*ffreq_ph)*(wb+bfreq_ph) + (wf+ffreq_ph), 2]+trileg_ph[(2*ffreq_ph)*(wb+bfreq_ph) + (wf+ffreq_ph), 4]) - U**2*re_chi_upup_ph(wb) 

def im_trileg_updo_ph(wb,wf):
	return 	U*(trileg_ph[(2*ffreq_ph)*(wb+bfreq_ph)+ (wf+ffreq_ph), 3]+trileg_ph[(2*ffreq_ph)*(wb+bfreq_ph)+ (wf+ffreq_ph), 3]+1) - U**2*im_chi_updo(wb)

def im_trileg_upup_ph(wb,wf):
	return U*trileg_ph[(2*ffreq_ph)*(wb+bfreq_ph)+ (wf+ffreq_ph), 3]-trileg_ph[(2*ffreq_ph)*(wb+bfreq_ph)+ (wf+ffreq_ph), 5]-U**2*im_chi_upup(wb)

# ----------define arrays to store in hdf5 file -> FROM HERE OBJECT IN THE SHIFTED NOTATION!!

re_trileg_arr_upup_ph = np.array([[[[[[[[re_trileg_upup_ph(m,n-myfloor_div2(m)) for i in range (0,QN_COUNT)] for j in range (0,QN_COUNT)]for ip in range(0,QN_COUNT)] for jp in range (0,QN_COUNT)] for q in range (0,PATCH_COUNT)] for k in range (0,PATCH_COUNT)] for n in range (-TRI_FFREQ_SHIFT,TRI_FFREQ_SHIFT)] for m in range (-TRI_BFREQ_SHIFT,TRI_BFREQ_SHIFT+1)] )

im_trileg_arr_upup_ph = np.array([[[[[[[[im_trileg_upup_ph(m,n-myfloor_div2(m)) for i in range (0,QN_COUNT)] for j in range (0,QN_COUNT)]for ip in range (0,QN_COUNT)] for jp in range (0,QN_COUNT)] for q in range (0,PATCH_COUNT)] for k in range (0,PATCH_COUNT)] for n in range (-TRI_FFREQ_SHIFT,TRI_FFREQ_SHIFT)] for m in range (-TRI_BFREQ_SHIFT,TRI_BFREQ_SHIFT+1)] )

re_trileg_arr_updo_ph = np.array([[[[[[[[[[re_trileg_updo_ph(m,n-myfloor_div2(m)) for i in range (0,QN_COUNT)] for j in range (0,QN_COUNT)] for ip in range (0,QN_COUNT)] for jp in range (0,QN_COUNT)] for q in range (0,PATCH_COUNT)] for k in range (0,PATCH_COUNT)] for n in range (-TRI_FFREQ_SHIFT,TRI_FFREQ_SHIFT)] for m in range (-TRI_BFREQ_SHIFT,TRI_BFREQ_SHIFT+1)] )

im_trileg_arr_updo_ph = np.array([[[[[[[[im_trileg_updo_ph(m,n-myfloor_div2(m)) for i in range (0,QN_COUNT)] for j in range (0,QN_COUNT)]for ip in range (0,QN_COUNT)] for jp in range (0,QN_COUNT)] for q in range (0,PATCH_COUNT)] for k in range (0,PATCH_COUNT)] for n in range (-TRI_FFREQ_SHIFT,TRI_FFREQ_SHIFT)] for m in range (-TRI_BFREQ_SHIFT,TRI_BFREQ_SHIFT+1)] )

fgrid_arr_ph = np.array([(2*m +1)*pi/BETA for m in range (-TRI_FFREQ_SHIFT,TRI_FFREQ_SHIFT)])
bgrid_arr_ph = np.array([(2*n)*pi/BETA for n in range (-TRI_BFREQ_SHIFT,TRI_BFREQ_SHIFT+1)])

#----------------------------------------Create HDF5 files-----------------------------------------

f.create_dataset('TRILEG/PH/RE_TRILEG_UPUP', data=re_trileg_arr_upup_ph, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('TRILEG/PH/IM_TRILEG_UPUP', data=im_trileg_arr_upup_ph, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('TRILEG/PH/RE_TRILEG_UPDO', data=re_trileg_arr_updo_ph, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('TRILEG/PH/IM_TRILEG_UPDO', data=im_trileg_arr_updo_ph, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('TRILEG/PH/fgrid', data=fgrid_arr_pp, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('TRILEG/PH/bgrid', data=bgrid_arr_pp, dtype='float64', compression="gzip", compression_opts=4)

f.create_dataset('TRILEG/XPH/RE_TRILEG_UPUP', data=-re_trileg_arr_upup_ph, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('TRILEG/XPH/IM_TRILEG_UPUP', data=-im_trileg_arr_upup_ph, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('TRILEG/XPH/RE_TRILEG_UPDO', data=re_trileg_arr_updo_ph-re_trileg_arr_upup_ph, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('TRILEG/XPH/IM_TRILEG_UPDO', data=im_trileg_arr_updo_ph-im_trileg_arr_upup_ph, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('TRILEG/XPH/fgrid', data=fgrid_arr_pp, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('TRILEG/XPH/bgrid', data=bgrid_arr_pp, dtype='float64', compression="gzip", compression_opts=4)


f.close()

