#!/usr/bin/python

#------------------------------------------------------------------------------------
#
# CONTENTS:
# This python script provides all the post processing operations from the (ED-> pomerol program) two-particle GF (2PGF) to other vertex objects:
#
# - Generalized Susceptibility -> Chi
# - Full vertex -> F
# - Fully Irredicible Vertex -> Lambda, Gamma, and Phi  (This will be done once we have the vertex asymptotics)
# 
# The input 2PGF are given in (PH, XPH, PP) mixed (bosonic, fermionic) shifted notation 
# in order to get a centralized vertex structure 
# 
# NOTE: One can easly translate one notation to another (either purely fermionic or mixed)
#       just using the functions defined in the library: agneselib (~/usr/include/agneselib)
#
# All the quantities are stored in a HDF5 file 
# 
# The final part of the program allows for the evaluation of the vertex asymptotics 
# Karrash and Plus functions by means of the script "vertex_asympt.py"
#
# To plot all the quantities use "plot.py"
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

U=1.0
MU= 0.0 # half-filling
BETA=50.0

FFREQ = 120 #fermionic frequencies in the mixed notation
BFREQ = 20   #bosonic frequency transfer in the mixed notation

CHI_BFREQ_ED = 1000
TRI_FFREQ_ED = 210 
TRI_FFREQ_SHIFT = 120
TRI_BFREQ_ED = 180
TRI_BFREQ_SHIFT = 180


PATCH_COUNT = 1

QN_COUNT = 1

pi = math.pi

#---------------------------------------------------------------------------------

def run(command):
        output = subprocess.check_output(command, shell=True)
        return output

#----------------------------------------Create HDF5 files-----------------------------------------

if ('/home/agnese/Coding/Vertex_structures/dat/H5FILES/BETA50/4SITES/reference_U1_W20/PREPROC/U'+ str(U)+'_BETA_'+ str(BETA)+'_FFREQ_'+ str(FFREQ)+'_BFREQ_' + str(BFREQ)+'.h5'):
    os.system('rm -r /home/agnese/Coding/Vertex_structures/dat/H5FILES/BETA50/4SITES/reference_U1_W20/PREPROC/U_'+ str(U)+'_BETA_'+ str(BETA)+'_FFREQ_'+ str(FFREQ)+'_BFREQ_' + str(BFREQ)+'.h5')

f = h5py.File('/home/agnese/Coding/Vertex_structures/dat/H5FILES/BETA50/4SITES/reference_U1_W20/PREPROC/U_'+ str(U)+'_BETA_'+ str(BETA)+'_FFREQ_'+ str(FFREQ)+'_BFREQ_' + str(BFREQ)+'.h5', 'w')   # Create the datafile.h5

#---------------------------------------------------------------------------------
#GF
g_iw  = np.loadtxt('/home/agnese/Coding/Vertex_structures/dat/pomerol/4SITES/U_'+str(U)+'_BETA_'+str(BETA)+'_FFREQ_'+ str(FFREQ)+'_BFREQ_'+ str(BFREQ)+'/gw_imfreq_00.dat')
N_fermi_gf = g_iw.shape[0]-FFREQ
print ("Number of fermionic frequencies for the GF:")
print N_fermi_gf

def G(w):                           # imaginary part of the GF
 if (w >= 0):
     return g_iw[w+FFREQ,2]+1j*g_iw[w+FFREQ,3]
 else:
     return g_iw[-w+FFREQ-1,2]-1j*g_iw[-w+FFREQ-1,3]

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

#---------------------------------------------------------------------------------------------------------
#
#                               PP
#
#---------------------------------------------------------------------------------------------------------
#--------------------------------------- 2PGF PP------------------------------------------------------------

twopgf_updo = np.loadtxt('/home/agnese/Coding/Vertex_structures/dat/pomerol/4SITES/U_'+str(U)+'_BETA_'+str(BETA)+'_FFREQ_'+ str(FFREQ)+'_BFREQ_'+ str(BFREQ)+'/2pgf_updo_pp_shift_W20.dat')
print twopgf_updo.shape
ffreq_pp = FFREQ
#int(0.5*(np.transpose(twopgf_updo)[1,:].max()*BETA/pi))+1
print ffreq_pp
bfreq_pp = BFREQ
#int(0.5*np.transpose(twopgf_updo)[0,:].max()*BETA/pi) 
print bfreq_pp


def re_2pgf_updo_pp(wb,wf,wf1):
    if(wb == 20):
        return 	twopgf_updo[(2*ffreq_pp)*(wf+ffreq_pp)+ (wf1+ffreq_pp), 3]

def re_2pgf_upup_pp(wb,wf,wf1):
    if(wb == 20):
        return 	re_2pgf_updo_pp(wb,wf,wf1)-re_2pgf_updo_pp(wb,wf,-wf1-1)

def im_2pgf_updo_pp(wb,wf,wf1):
    if(wb ==20):
        return 	twopgf_updo[(2*ffreq_pp)*(wf+ffreq_pp)+ (wf1+ffreq_pp), 4]

def im_2pgf_upup_pp(wb,wf,wf1):
    if(wb == 20):
        return 	im_2pgf_updo_pp(wb,wf,wf1)-im_2pgf_updo_pp(wb,wf,-wf1-1)

# ----------define arrays to store in hdf5 file

re_2pgf_arr_upup_pp = np.array([[[[[[[[[[re_2pgf_upup_pp(m,n,p) for i in range (0,QN_COUNT)] for j in range (0,QN_COUNT)]for ip in range(0,QN_COUNT)] for jp in range (0,QN_COUNT)] for q in range (0,PATCH_COUNT)] for k in range (0,PATCH_COUNT)] for kp in range (0,PATCH_COUNT)] for p in range (-ffreq_pp,ffreq_pp )] for n in range (-ffreq_pp,ffreq_pp)] for m in range (bfreq_pp,bfreq_pp+1)] )

im_2pgf_arr_upup_pp = np.array([[[[[[[[[[im_2pgf_upup_pp(m,n,p) for i in range (0,QN_COUNT)] for j in range (0,QN_COUNT)]for ip in range (0,QN_COUNT)] for jp in range (0,QN_COUNT)] for q in range (0,PATCH_COUNT)] for k in range (0,PATCH_COUNT)] for kp in range (0,PATCH_COUNT)] for p in range (-ffreq_pp,ffreq_pp )] for n in range (-ffreq_pp,ffreq_pp)] for m in range (bfreq_pp,bfreq_pp+1)] )

re_2pgf_arr_updo_pp = np.array([[[[[[[[[[re_2pgf_updo_pp(m,n,p) for i in range (0,QN_COUNT)] for j in range (0,QN_COUNT)] for ip in range (0,QN_COUNT)] for jp in range (0,QN_COUNT)] for q in range (0,PATCH_COUNT)] for k in range (0,PATCH_COUNT)] for kp in range (0,PATCH_COUNT)] for p in range (-ffreq_pp,ffreq_pp )] for n in range (-ffreq_pp,ffreq_pp)] for m in range (bfreq_pp,bfreq_pp+1)] )

im_2pgf_arr_updo_pp = np.array([[[[[[[[[[im_2pgf_updo_pp(m,n,p) for i in range (0,QN_COUNT)] for j in range (0,QN_COUNT)]for ip in range (0,QN_COUNT)] for jp in range (0,QN_COUNT)] for q in range (0,PATCH_COUNT)] for k in range (0,PATCH_COUNT)] for kp in range (0,PATCH_COUNT)] for p in range (-ffreq_pp,ffreq_pp )] for n in range (-ffreq_pp,ffreq_pp)] for m in range (bfreq_pp,bfreq_pp+1)] )

fgrid_arr_pp = np.array([(2*m +1)*pi/BETA for m in range (-ffreq_pp,ffreq_pp)])
bgrid_arr_pp = np.array([(2*n)*pi/BETA for n in range (bfreq_pp,bfreq_pp+1)])

#----------------------------------------Create HDF5 files-----------------------------------------

f.create_dataset('2PGF/PP/RE_2PGF_UPUP', data=re_2pgf_arr_upup_pp, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('2PGF/PP/IM_2PGF_UPUP', data=im_2pgf_arr_upup_pp, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('2PGF/PP/RE_2PGF_UPDO', data=re_2pgf_arr_updo_pp, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('2PGF/PP/IM_2PGF_UPDO', data=im_2pgf_arr_updo_pp, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('2PGF/PP/fgrid', data=fgrid_arr_pp, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('2PGF/PP/bgrid', data=bgrid_arr_pp, dtype='float64', compression="gzip", compression_opts=4)

#------------------------------------- GENERALIZED SUSCEPTIBILITY CHI PP --------------------------------------------

# Cutting one disconnected diagram -> generalized susceptibility


def chi_upup_pp(wb,wf,wf1):
    if (wf == -wf1-1-mymod_abs(wb)):
        return re_2pgf_upup_pp(wb,wf,wf1)+1j*im_2pgf_upup_pp(wb,wf,wf1)-BETA*(G(wf+myceil_div2(wb))*G(myfloor_div2(wb)-wf-1))
    else:
        return re_2pgf_upup_pp(wb,wf,wf1)+1j*im_2pgf_upup_pp(wb,wf,wf1)

def chi_updo_pp(wb,wf,wf1):
    if (wf == -wf1-1-mymod_abs(wb)):
        return re_2pgf_updo_pp(wb,wf,wf1)+1j*im_2pgf_updo_pp(wb,wf,wf1)-BETA*(G(wf+myceil_div2(wb))*G(myfloor_div2(wb)-wf-1))
    else:
        return re_2pgf_updo_pp(wb,wf,wf1)+1j*im_2pgf_updo_pp(wb,wf,wf1)

def chi_s_pp(wb,wf,wf1):
    return -chi_upup_pp(wb,wf,wf1) + 2*chi_updo_pp(wb,wf,wf1)

def chi_t_pp(wb,wf,wf1):
    return chi_upup_pp(wb,wf,wf1)


# ----------define arrays to store in hdf5 file

re_chi_arr_upup_pp = np.array([[[[[[[[[[chi_upup_pp(m,n,p).real for i in range (0,QN_COUNT)] for j in range (0,QN_COUNT)]for ip in range(0,QN_COUNT)] for jp in range (0,QN_COUNT)] for q in range (0,PATCH_COUNT)] for k in range (0,PATCH_COUNT)] for kp in range (0,PATCH_COUNT)] for p in range (-ffreq_pp,ffreq_pp )] for n in range (-ffreq_pp,ffreq_pp)] for m in range (bfreq_pp,bfreq_pp+1)] )

im_chi_arr_upup_pp = np.array([[[[[[[[[[chi_upup_pp(m,n,p).imag for i in range (0,QN_COUNT)] for j in range (0,QN_COUNT)]for ip in range(0,QN_COUNT)] for jp in range (0,QN_COUNT)] for q in range (0,PATCH_COUNT)] for k in range (0,PATCH_COUNT)] for kp in range (0,PATCH_COUNT)] for p in range (-ffreq_pp,ffreq_pp )] for n in range (-ffreq_pp,ffreq_pp)] for m in range (bfreq_pp,bfreq_pp+1)] )

re_chi_arr_updo_pp = np.array([[[[[[[[[[chi_updo_pp(m,n,p).real for i in range (0,QN_COUNT)] for j in range (0,QN_COUNT)]for ip in range(0,QN_COUNT)] for jp in range (0,QN_COUNT)] for q in range (0,PATCH_COUNT)] for k in range (0,PATCH_COUNT)] for kp in range (0,PATCH_COUNT)] for p in range (-ffreq_pp,ffreq_pp )] for n in range (-ffreq_pp,ffreq_pp)] for m in range (bfreq_pp,bfreq_pp+1)] )

im_chi_arr_updo_pp = np.array([[[[[[[[[[chi_updo_pp(m,n,p).imag for i in range (0,QN_COUNT)] for j in range (0,QN_COUNT)]for ip in range(0,QN_COUNT)] for jp in range (0,QN_COUNT)] for q in range (0,PATCH_COUNT)] for k in range (0,PATCH_COUNT)] for kp in range (0,PATCH_COUNT)] for p in range (-ffreq_pp,ffreq_pp )] for n in range (-ffreq_pp,ffreq_pp)] for m in range (bfreq_pp,bfreq_pp+1)] )

#----------------------------------------Create HDF5 files-----------------------------------------

f.create_dataset('GENCHI/PP/RE_GENCHI_UPUP', data=re_chi_arr_upup_pp, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('GENCHI/PP/IM_GENCHI_UPUP', data=im_chi_arr_upup_pp, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('GENCHI/PP/RE_GENCHI_UPDO', data=re_chi_arr_updo_pp, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('GENCHI/PP/IM_GENCHI_UPDO', data=im_chi_arr_updo_pp, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('GENCHI/PP/fgrid', data=fgrid_arr_pp, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('GENCHI/PP/bgrid', data=bgrid_arr_pp, dtype='float64', compression="gzip", compression_opts=4)


#----------------------------------- FULL VERTEX PP -----------------------------------------------------------

def chi_l_pp(wb,wf,wf1):  # Legs on one side of the diagram
    return BETA*G(wf+myceil_div2(wb))*G(myfloor_div2(wb)-wf1-1)

def chi_0_pp(wb,wf,wf1):
    if (wf == wf1):
        return BETA*G(wf+myceil_div2(wb))*G(myfloor_div2(wb)-wf-1)
    else:
        return 0.0

def f_upup_pp(wb,wf,wf1):
    return BETA*BETA*(1.0/(chi_l_pp(wb,wf,wf)))*(chi_upup_pp(wb,wf,wf1)+chi_0_pp(wb,wf,wf1))*(1.0/(chi_l_pp(wb,wf1,wf1)))

def f_updo_pp(wb,wf,wf1):
    return BETA*BETA*(1.0/(chi_l_pp(wb,wf,wf)))*(chi_updo_pp(wb,wf,wf1))*(1.0/(chi_l_pp(wb,wf1,wf1)))

# ----------define arrays to store in hdf5 file

re_f_arr_upup_pp = np.array([[[[[[[[[[f_upup_pp(m,n,p).real for i in range (0,QN_COUNT)] for j in range (0,QN_COUNT)]for ip in range(0,QN_COUNT)] for jp in range (0,QN_COUNT)] for q in range (0,PATCH_COUNT)] for k in range (0,PATCH_COUNT)] for kp in range (0,PATCH_COUNT)] for p in range (-ffreq_pp,ffreq_pp )] for n in range (-ffreq_pp,ffreq_pp)] for m in range (bfreq_pp,bfreq_pp+1)] )

im_f_arr_upup_pp = np.array([[[[[[[[[[f_upup_pp(m,n,p).imag for i in range (0,QN_COUNT)] for j in range (0,QN_COUNT)]for ip in range(0,QN_COUNT)] for jp in range (0,QN_COUNT)] for q in range (0,PATCH_COUNT)] for k in range (0,PATCH_COUNT)] for kp in range (0,PATCH_COUNT)] for p in range (-ffreq_pp,ffreq_pp )] for n in range (-ffreq_pp,ffreq_pp)] for m in range (bfreq_pp,bfreq_pp+1)] )

re_f_arr_updo_pp = np.array([[[[[[[[[[f_updo_pp(m,n,p).real for i in range (0,QN_COUNT)] for j in range (0,QN_COUNT)]for ip in range(0,QN_COUNT)] for jp in range (0,QN_COUNT)] for q in range (0,PATCH_COUNT)] for k in range (0,PATCH_COUNT)] for kp in range (0,PATCH_COUNT)] for p in range (-ffreq_pp,ffreq_pp )] for n in range (-ffreq_pp,ffreq_pp)] for m in range (bfreq_pp,bfreq_pp+1)] )

im_f_arr_updo_pp = np.array([[[[[[[[[[f_updo_pp(m,n,p).imag for i in range (0,QN_COUNT)] for j in range (0,QN_COUNT)]for ip in range(0,QN_COUNT)] for jp in range (0,QN_COUNT)] for q in range (0,PATCH_COUNT)] for k in range (0,PATCH_COUNT)] for kp in range (0,PATCH_COUNT)] for p in range (-ffreq_pp,ffreq_pp )] for n in range (-ffreq_pp,ffreq_pp)] for m in range (bfreq_pp,bfreq_pp+1)] )

#----------------------------------------Create HDF5 files-----------------------------------------

f.create_dataset('VERT/PP/RE_F_UPUP', data=re_f_arr_upup_pp, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('VERT/PP/IM_F_UPUP', data=im_f_arr_upup_pp, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('VERT/PP/RE_F_UPDO', data=re_f_arr_updo_pp, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('VERT/PP/IM_F_UPDO', data=im_f_arr_updo_pp, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('VERT/PP/fgrid', data=fgrid_arr_pp, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('VERT/PP/bgrid', data=bgrid_arr_pp, dtype='float64', compression="gzip", compression_opts=4)

#---------------------------------------------------------------------------------------------------------
#
#                               PH
#
#---------------------------------------------------------------------------------------------------------

#--------------------------------------- 2PGF PH------------------------------------------------------------

twopgf_updo_ph_l = np.loadtxt('/home/agnese/Coding/Vertex_structures/dat/pomerol/4SITES/U_'+str(U)+'_BETA_'+str(BETA)+'_FFREQ_'+ str(FFREQ)+'_BFREQ_'+ str(BFREQ)+'/2pgf_updo_ph_shift_W20.dat')
twopgf_updo_xph_l = np.loadtxt('/home/agnese/Coding/Vertex_structures/dat/pomerol/4SITES/U_'+str(U)+'_BETA_'+str(BETA)+'_FFREQ_'+ str(FFREQ)+'_BFREQ_'+ str(BFREQ)+'/2pgf_updo_xph_shift_W20.dat')
ffreq_ph = FFREQ 
#int(0.5*(np.transpose(twopgf_updo_ph_l)[1,:].max()*BETA/pi-1))+1
print ffreq_ph
bfreq_ph = BFREQ
#int(0.5*np.transpose(twopgf_updo_ph_l)[0,:].max()*BETA/pi)+1
print bfreq_ph

def re_2pgf_upup_ph(wb,wf,wf1):
    if(wb==20):
        return 	twopgf_updo_ph_l[(2*ffreq_ph)*(wf+ffreq_ph) + (wf1+ffreq_ph), 3]-twopgf_updo_xph_l[(2*ffreq_ph)*(wf+ffreq_ph) + (wf1+ffreq_ph), 3]

def re_2pgf_updo_ph(wb,wf,wf1):
    if(wb == 20):
        return 	twopgf_updo_ph_l[(2*ffreq_ph)*(wf+ffreq_ph)+ (wf1+ffreq_ph), 3]

def im_2pgf_upup_ph(wb,wf,wf1):
    if(wb == 20):
        return 	twopgf_updo_ph_l[(2*ffreq_ph)*(wf+ffreq_ph) + (wf1+ffreq_ph), 4]-twopgf_updo_xph_l[(2*ffreq_ph)*(wf+ffreq_ph) + (wf1+ffreq_ph), 4]

def im_2pgf_updo_ph(wb,wf,wf1):
    if(wb == 20):
        return 	twopgf_updo_ph_l[(2*ffreq_ph)*(wf+ffreq_ph)+ (wf1+ffreq_ph), 4]

# ----------define arrays to store in hdf5 file

re_2pgf_arr_upup_ph = np.array([[[[[[[[[[re_2pgf_upup_ph(m,n,p) for i in range (0,QN_COUNT)] for j in range (0,QN_COUNT)]for ip in range(0,QN_COUNT)] for jp in range (0,QN_COUNT)] for q in range (0,PATCH_COUNT)] for k in range (0,PATCH_COUNT)] for kp in range (0,PATCH_COUNT)] for p in range (-ffreq_ph,ffreq_ph )] for n in range (-ffreq_ph,ffreq_ph)] for m in range (bfreq_ph,bfreq_ph+1)] )

im_2pgf_arr_upup_ph = np.array([[[[[[[[[[im_2pgf_upup_ph(m,n,p) for i in range (0,QN_COUNT)] for j in range (0,QN_COUNT)]for ip in range(0,QN_COUNT)] for jp in range (0,QN_COUNT)] for q in range (0,PATCH_COUNT)] for k in range (0,PATCH_COUNT)] for kp in range (0,PATCH_COUNT)] for p in range (-ffreq_ph,ffreq_ph )] for n in range (-ffreq_ph,ffreq_ph)] for m in range (bfreq_ph,bfreq_ph+1)] )

re_2pgf_arr_updo_ph = np.array([[[[[[[[[[re_2pgf_updo_ph(m,n,p) for i in range (0,QN_COUNT)] for j in range (0,QN_COUNT)]for ip in range(0,QN_COUNT)] for jp in range (0,QN_COUNT)] for q in range (0,PATCH_COUNT)] for k in range (0,PATCH_COUNT)] for kp in range (0,PATCH_COUNT)] for p in range (-ffreq_ph,ffreq_ph )] for n in range (-ffreq_ph,ffreq_ph)] for m in range (bfreq_ph,bfreq_ph+1)] )

im_2pgf_arr_updo_ph = np.array([[[[[[[[[[im_2pgf_updo_ph(m,n,p) for i in range (0,QN_COUNT)] for j in range (0,QN_COUNT)]for ip in range(0,QN_COUNT)] for jp in range (0,QN_COUNT)] for q in range (0,PATCH_COUNT)] for k in range (0,PATCH_COUNT)] for kp in range (0,PATCH_COUNT)] for p in range (-ffreq_ph,ffreq_ph )] for n in range (-ffreq_ph,ffreq_ph)] for m in range (bfreq_ph,bfreq_ph+1)] )

fgrid_arr_ph = np.array([(2*m +1)*pi/BETA for m in range (-ffreq_ph,ffreq_ph)])
bgrid_arr_ph = np.array([(2*n)*pi/BETA for n in range (bfreq_ph,bfreq_ph+1)])

#----------------------------------------Create HDF5 files-----------------------------------------

f.create_dataset('2PGF/PH/RE_2PGF_UPUP', data=re_2pgf_arr_upup_ph, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('2PGF/PH/IM_2PGF_UPUP', data=im_2pgf_arr_upup_ph, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('2PGF/PH/RE_2PGF_UPDO', data=re_2pgf_arr_updo_ph, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('2PGF/PH/IM_2PGF_UPDO', data=im_2pgf_arr_updo_ph, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('2PGF/PH/fgrid', data=fgrid_arr_ph, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('2PGF/PH/bgrid', data=bgrid_arr_ph, dtype='float64', compression="gzip", compression_opts=4)

#------------------------------------- GENERALIZED SUSCEPTIBILITY CHI PH --------------------------------------------

# Cutting one disconnected diagram -> generalized susceptibility

def chi_upup_ph(wb,wf,wf1):
    if (wb== 0):
        return re_2pgf_upup_ph(wb,wf,wf1)+1j*im_2pgf_upup_ph(wb,wf,wf1)-BETA*(G(wf)*G(wf1))
    else:
        return re_2pgf_upup_ph(wb,wf,wf1)+1j*im_2pgf_upup_ph(wb,wf,wf1)

def chi_updo_ph(wb,wf,wf1):
    if (wb== 0):
        return re_2pgf_updo_ph(wb,wf,wf1)+1j*im_2pgf_updo_ph(wb,wf,wf1)-BETA*(G(wf)*G(wf1))
    else:
        return re_2pgf_updo_ph(wb,wf,wf1)+1j*im_2pgf_updo_ph(wb,wf,wf1)

def tpgf_d_ph(wb,wf,wf1):
    return re_2pgf_upup_ph(wb,wf,wf1)+1j*im_2pgf_upup_ph(wb,wf,wf1)+re_2pgf_updo_ph(wb,wf,wf1)+1j*im_2pgf_updo_ph(wb,wf,wf1)

def tpgf_m_ph(wb,wf,wf1):
    return re_2pgf_upup_ph(wb,wf,wf1)+1j*im_2pgf_upup_ph(wb,wf,wf1) -re_2pgf_updo_ph(wb,wf,wf1)-1j*im_2pgf_updo_ph(wb,wf,wf1)

# ----------define arrays to store in hdf5 file

re_chi_arr_upup_ph = np.array([[[[[[[[[[chi_upup_ph(m,n,p).real for i in range (0,QN_COUNT)] for j in range (0,QN_COUNT)]for ip in range(0,QN_COUNT)] for jp in range (0,QN_COUNT)] for q in range (0,PATCH_COUNT)] for k in range (0,PATCH_COUNT)] for kp in range (0,PATCH_COUNT)] for p in range (-ffreq_ph,ffreq_ph )] for n in range (-ffreq_ph,ffreq_ph)] for m in range (bfreq_ph,bfreq_ph+1)] )

im_chi_arr_upup_ph = np.array([[[[[[[[[[chi_upup_ph(m,n,p).imag for i in range (0,QN_COUNT)] for j in range (0,QN_COUNT)]for ip in range(0,QN_COUNT)] for jp in range (0,QN_COUNT)] for q in range (0,PATCH_COUNT)] for k in range (0,PATCH_COUNT)] for kp in range (0,PATCH_COUNT)] for p in range (-ffreq_ph,ffreq_ph )] for n in range (-ffreq_ph,ffreq_ph)] for m in range (bfreq_ph,bfreq_ph+1)] )

re_chi_arr_updo_ph = np.array([[[[[[[[[[chi_updo_ph(m,n,p).real for i in range (0,QN_COUNT)] for j in range (0,QN_COUNT)]for ip in range(0,QN_COUNT)] for jp in range (0,QN_COUNT)] for q in range (0,PATCH_COUNT)] for k in range (0,PATCH_COUNT)] for kp in range (0,PATCH_COUNT)] for p in range (-ffreq_ph,ffreq_ph )] for n in range (-ffreq_ph,ffreq_ph)] for m in range (bfreq_ph,bfreq_ph+1)] )

im_chi_arr_updo_ph = np.array([[[[[[[[[[chi_updo_ph(m,n,p).imag for i in range (0,QN_COUNT)] for j in range (0,QN_COUNT)]for ip in range(0,QN_COUNT)] for jp in range (0,QN_COUNT)] for q in range (0,PATCH_COUNT)] for k in range (0,PATCH_COUNT)] for kp in range (0,PATCH_COUNT)] for p in range (-ffreq_ph,ffreq_ph )] for n in range (-ffreq_ph,ffreq_ph)] for m in range (bfreq_ph,bfreq_ph+1)] )

#----------------------------------------Create HDF5 files-----------------------------------------

f.create_dataset('GENCHI/PH/RE_GENCHI_UPUP', data=re_chi_arr_upup_ph, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('GENCHI/PH/IM_GENCHI_UPUP', data=im_chi_arr_upup_ph, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('GENCHI/PH/RE_GENCHI_UPDO', data=re_chi_arr_updo_ph, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('GENCHI/PH/IM_GENCHI_UPDO', data=im_chi_arr_updo_ph, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('GENCHI/PH/fgrid', data=fgrid_arr_ph, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('GENCHI/PH/bgrid', data=bgrid_arr_ph, dtype='float64', compression="gzip", compression_opts=4)

#----------------------------------- FULL VERTEX PH -----------------------------------------------------------

def chi_lup_ph(wb,wf,wf1):
    return BETA*G(wf-myfloor_div2(wb))*G(wf1+myceil_div2(wb))

def chi_0_ph(wb,wf,wf1):
    if (wf == wf1):
        return BETA*G(wf+myceil_div2(wb))*G(wf-myfloor_div2(wb))
    else:
        return 0.0


def f_upup_ph(wb,wf,wf1):
    return BETA*BETA*(1.0/chi_lup_ph(wb,wf,wf))*(chi_upup_ph(wb,wf,wf1)+chi_0_ph(wb,wf,wf1))*(1.0/chi_lup_ph(wb,wf1,wf1))

def f_updo_ph(wb,wf,wf1):
    return BETA*BETA*(1.0/chi_lup_ph(wb,wf,wf))*(chi_updo_ph(wb,wf,wf1))*(1.0/chi_lup_ph(wb,wf1,wf1))


# ----------define arrays to store in hdf5 file

re_f_arr_upup_ph = np.array([[[[[[[[[[f_upup_ph(m,n,p).real for i in range (0,QN_COUNT)] for j in range (0,QN_COUNT)]for ip in range(0,QN_COUNT)] for jp in range (0,QN_COUNT)] for q in range (0,PATCH_COUNT)] for k in range (0,PATCH_COUNT)] for kp in range (0,PATCH_COUNT)] for p in range (-ffreq_pp,ffreq_pp )] for n in range (-ffreq_ph,ffreq_ph)] for m in range (bfreq_ph,bfreq_ph+1)] )

im_f_arr_upup_ph = np.array([[[[[[[[[[f_upup_ph(m,n,p).imag for i in range (0,QN_COUNT)] for j in range (0,QN_COUNT)]for ip in range(0,QN_COUNT)] for jp in range (0,QN_COUNT)] for q in range (0,PATCH_COUNT)] for k in range (0,PATCH_COUNT)] for kp in range (0,PATCH_COUNT)] for p in range (-ffreq_pp,ffreq_pp )] for n in range (-ffreq_ph,ffreq_ph)] for m in range (bfreq_ph,bfreq_ph+1)] )

re_f_arr_updo_ph = np.array([[[[[[[[[[f_updo_ph(m,n,p).real for i in range (0,QN_COUNT)] for j in range (0,QN_COUNT)]for ip in range(0,QN_COUNT)] for jp in range (0,QN_COUNT)] for q in range (0,PATCH_COUNT)] for k in range (0,PATCH_COUNT)] for kp in range (0,PATCH_COUNT)] for p in range (-ffreq_pp,ffreq_pp )] for n in range (-ffreq_ph,ffreq_ph)] for m in range (bfreq_ph,bfreq_ph+1)] )

im_f_arr_updo_ph = np.array([[[[[[[[[[f_updo_ph(m,n,p).imag for i in range (0,QN_COUNT)] for j in range (0,QN_COUNT)]for ip in range(0,QN_COUNT)] for jp in range (0,QN_COUNT)] for q in range (0,PATCH_COUNT)] for k in range (0,PATCH_COUNT)] for kp in range (0,PATCH_COUNT)] for p in range (-ffreq_pp,ffreq_pp )] for n in range (-ffreq_ph,ffreq_ph)] for m in range (bfreq_ph,bfreq_ph+1)] )

#----------------------------------------Create HDF5 files-----------------------------------------

f.create_dataset('VERT/PH/RE_F_UPUP', data=re_f_arr_upup_ph, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('VERT/PH/IM_F_UPUP', data=im_f_arr_upup_ph, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('VERT/PH/RE_F_UPDO', data=re_f_arr_updo_ph, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('VERT/PH/IM_F_UPDO', data=im_f_arr_updo_ph, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('VERT/PH/fgrid', data=fgrid_arr_ph, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('VERT/PH/bgrid', data=bgrid_arr_ph, dtype='float64', compression="gzip", compression_opts=4)

#---------------------------------------------------------------------------------------------------------
#
#                               XPH
#
#---------------------------------------------------------------------------------------------------------
#--------------------------------------- 2PGF XPH------------------------------------------------------------

ffreq_xph = FFREQ
#int(0.5*(np.transpose(twopgf_updo_xph_l)[1,:].max()*BETA/pi-1))+1
print ffreq_xph
bfreq_xph = BFREQ
#int(0.5*np.transpose(twopgf_updo_xph_l)[0,:].max()*BETA/pi)+1
print bfreq_xph

def re_2pgf_upup_xph(wb,wf,wf1):
    if(wb == 20):
        return 	twopgf_updo_xph_l[(2*ffreq_xph)*(wf+ffreq_xph) + (wf1+ffreq_xph), 3]-twopgf_updo_ph_l[(2*ffreq_xph)*(wf+ffreq_xph) + (wf1+ffreq_xph), 3]

def re_2pgf_updo_xph(wb,wf,wf1):
    if(wb == 20):
        return 	twopgf_updo_xph_l[(2*ffreq_xph)*(wf+ffreq_xph)+ (wf1+ffreq_xph), 3]

def im_2pgf_upup_xph(wb,wf,wf1):
    if(wb == 20):
	return  twopgf_updo_xph_l[(2*ffreq_xph)*(wf+ffreq_xph) + (wf1+ffreq_xph), 4]-twopgf_updo_ph_l[(2*ffreq_xph)*(wf+ffreq_xph) + (wf1+ffreq_xph), 4]

def im_2pgf_updo_xph(wb,wf,wf1):
    if(wb == 20):
        return 	 twopgf_updo_xph_l[(2*ffreq_xph)*(wf+ffreq_xph)+ (wf1+ffreq_xph), 4]

# ----------define arrays to store in hdf5 file

re_2pgf_arr_upup_xph = np.array([[[[[[[[[[re_2pgf_upup_xph(m,n,p) for i in range (0,QN_COUNT)] for j in range (0,QN_COUNT)]for ip in range(0,QN_COUNT)] for jp in range (0,QN_COUNT)] for q in range (0,PATCH_COUNT)] for k in range (0,PATCH_COUNT)] for kp in range (0,PATCH_COUNT)] for p in range (-ffreq_xph,ffreq_xph )] for n in range (-ffreq_xph,ffreq_xph)] for m in range (bfreq_xph,bfreq_xph+1)] )

im_2pgf_arr_upup_xph = np.array([[[[[[[[[[im_2pgf_upup_xph(m,n,p) for i in range (0,QN_COUNT)] for j in range (0,QN_COUNT)]for ip in range(0,QN_COUNT)] for jp in range (0,QN_COUNT)] for q in range (0,PATCH_COUNT)] for k in range (0,PATCH_COUNT)] for kp in range (0,PATCH_COUNT)] for p in range (-ffreq_xph,ffreq_xph )] for n in range (-ffreq_xph,ffreq_xph)] for m in range (bfreq_xph,bfreq_xph+1)] )

re_2pgf_arr_updo_xph = np.array([[[[[[[[[[re_2pgf_updo_xph(m,n,p) for i in range (0,QN_COUNT)] for j in range (0,QN_COUNT)]for ip in range(0,QN_COUNT)] for jp in range (0,QN_COUNT)] for q in range (0,PATCH_COUNT)] for k in range (0,PATCH_COUNT)] for kp in range (0,PATCH_COUNT)] for p in range (-ffreq_xph,ffreq_xph )] for n in range (-ffreq_xph,ffreq_xph)] for m in range (bfreq_xph,bfreq_xph+1)] )

im_2pgf_arr_updo_xph = np.array([[[[[[[[[[im_2pgf_updo_xph(m,n,p) for i in range (0,QN_COUNT)] for j in range (0,QN_COUNT)]for ip in range(0,QN_COUNT)] for jp in range (0,QN_COUNT)] for q in range (0,PATCH_COUNT)] for k in range (0,PATCH_COUNT)] for kp in range (0,PATCH_COUNT)] for p in range (-ffreq_xph,ffreq_xph )] for n in range (-ffreq_xph,ffreq_xph)] for m in range (bfreq_xph,bfreq_xph+1)] )

fgrid_arr_xph = np.array([(2*m +1)*pi/BETA for m in range (-ffreq_xph,ffreq_xph)])
bgrid_arr_xph = np.array([(2*n)*pi/BETA for n in range (bfreq_xph,bfreq_xph+1)])

#----------------------------------------Create HDF5 files-----------------------------------------

f.create_dataset('2PGF/XPH/RE_2PGF_UPUP', data=re_2pgf_arr_upup_xph, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('2PGF/XPH/IM_2PGF_UPUP', data=im_2pgf_arr_upup_xph, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('2PGF/XPH/RE_2PGF_UPDO', data=re_2pgf_arr_updo_xph, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('2PGF/XPH/IM_2PGF_UPDO', data=im_2pgf_arr_updo_xph, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('2PGF/XPH/fgrid', data=fgrid_arr_xph, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('2PGF/XPH/bgrid', data=bgrid_arr_xph, dtype='float64', compression="gzip", compression_opts=4)

#------------------------------------- GENERALIZED SUSCEPTIBILITY CHI XPH --------------------------------------------

# Cutting one disconnected diagram -> generalized susceptibility

def chi_upup_xph(wb,wf,wf1):
    if (wb == 0):
        return re_2pgf_upup_xph(wb,wf,wf1)+1j*im_2pgf_upup_xph(wb,wf,wf1)+BETA*G(wf)*G(wf1)
    else:
        return re_2pgf_upup_xph(wb,wf,wf1)+1j*im_2pgf_upup_xph(wb,wf,wf1)

def chi_updo_xph(wb,wf,wf1):
    return re_2pgf_updo_xph(wb,wf,wf1)+1j*im_2pgf_updo_xph(wb,wf,wf1)

# ----------define arrays to store in hdf5 file

re_chi_arr_upup_xph = np.array([[[[[[[[[[chi_upup_xph(m,n,p).real for i in range (0,QN_COUNT)] for j in range (0,QN_COUNT)]for ip in range(0,QN_COUNT)] for jp in range (0,QN_COUNT)] for q in range (0,PATCH_COUNT)] for k in range (0,PATCH_COUNT)] for kp in range (0,PATCH_COUNT)] for p in range (-ffreq_xph,ffreq_xph )] for n in range (-ffreq_xph,ffreq_xph)] for m in range (bfreq_xph,bfreq_xph+1)] )

im_chi_arr_upup_xph = np.array([[[[[[[[[[chi_upup_xph(m,n,p).imag for i in range (0,QN_COUNT)] for j in range (0,QN_COUNT)]for ip in range(0,QN_COUNT)] for jp in range (0,QN_COUNT)] for q in range (0,PATCH_COUNT)] for k in range (0,PATCH_COUNT)] for kp in range (0,PATCH_COUNT)] for p in range (-ffreq_xph,ffreq_xph )] for n in range (-ffreq_xph,ffreq_xph)] for m in range (bfreq_xph,bfreq_xph+1)] )

re_chi_arr_updo_xph = np.array([[[[[[[[[[chi_updo_xph(m,n,p).real for i in range (0,QN_COUNT)] for j in range (0,QN_COUNT)]for ip in range(0,QN_COUNT)] for jp in range (0,QN_COUNT)] for q in range (0,PATCH_COUNT)] for k in range (0,PATCH_COUNT)] for kp in range (0,PATCH_COUNT)] for p in range (-ffreq_xph,ffreq_xph )] for n in range (-ffreq_xph,ffreq_xph)] for m in range (bfreq_xph,bfreq_xph+1)] )

im_chi_arr_updo_xph = np.array([[[[[[[[[[chi_updo_xph(m,n,p).imag for i in range (0,QN_COUNT)] for j in range (0,QN_COUNT)]for ip in range(0,QN_COUNT)] for jp in range (0,QN_COUNT)] for q in range (0,PATCH_COUNT)] for k in range (0,PATCH_COUNT)] for kp in range (0,PATCH_COUNT)] for p in range (-ffreq_xph,ffreq_xph )] for n in range (-ffreq_xph,ffreq_xph)] for m in range (bfreq_xph,bfreq_xph+1)] )

#----------------------------------------Create HDF5 files-----------------------------------------
f.create_dataset('GENCHI/XPH/RE_GENCHI_UPUP', data=re_chi_arr_upup_xph, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('GENCHI/XPH/IM_GENCHI_UPUP', data=im_chi_arr_upup_xph, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('GENCHI/XPH/RE_GENCHI_UPDO', data=re_chi_arr_updo_xph, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('GENCHI/XPH/IM_GENCHI_UPDO', data=im_chi_arr_updo_xph, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('GENCHI/XPH/fgrid', data=fgrid_arr_xph, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('GENCHI/XPH/bgrid', data=bgrid_arr_xph, dtype='float64', compression="gzip", compression_opts=4)

#----------------------------------- FULL VERTEX XPH -----------------------------------------------------------

def chi_l_xph(wb,wf,wf1):
    return BETA*G(wf-myfloor_div2(wb))*G(wf1+myceil_div2(wb))

def chi_0_xph(wb,wf,wf1):
    if (wf == wf1):
        return BETA*G(wf-myfloor_div2(wb))*G(wf+myceil_div2(wb))
    else:
        return 0.0

def f_upup_xph(wb,wf,wf1):
    return BETA*BETA*(1.0/chi_l_xph(wb,wf,wf1))*(chi_upup_xph(wb,wf,wf1)-chi_0_xph(wb,wf,wf1))*(1.0/chi_l_xph(wb,wf1,wf))

def f_updo_xph(wb,wf,wf1):
    return BETA*BETA*(1.0/chi_l_xph(wb,wf,wf1))*(chi_updo_xph(wb,wf,wf1)-chi_0_xph(wb,wf,wf1))*(1.0/chi_l_xph(wb,wf1,wf))


# ----------define arrays to store in hdf5 file

re_f_arr_upup_xph = np.array([[[[[[[[[[f_upup_xph(m,n,p).real for i in range (0,QN_COUNT)] for j in range (0,QN_COUNT)]for ip in range(0,QN_COUNT)] for jp in range (0,QN_COUNT)] for q in range (0,PATCH_COUNT)] for k in range (0,PATCH_COUNT)] for kp in range (0,PATCH_COUNT)] for p in range (-ffreq_xph,ffreq_xph )] for n in range (-ffreq_xph,ffreq_xph)] for m in range (bfreq_xph,bfreq_xph+1)] )

im_f_arr_upup_xph = np.array([[[[[[[[[[f_upup_xph(m,n,p).imag for i in range (0,QN_COUNT)] for j in range (0,QN_COUNT)]for ip in range(0,QN_COUNT)] for jp in range (0,QN_COUNT)] for q in range (0,PATCH_COUNT)] for k in range (0,PATCH_COUNT)] for kp in range (0,PATCH_COUNT)] for p in range (-ffreq_xph,ffreq_xph )] for n in range (-ffreq_xph,ffreq_xph)] for m in range (bfreq_xph,bfreq_xph+1)] )

re_f_arr_updo_xph = np.array([[[[[[[[[[f_updo_xph(m,n,p).real for i in range (0,QN_COUNT)] for j in range (0,QN_COUNT)]for ip in range(0,QN_COUNT)] for jp in range (0,QN_COUNT)] for q in range (0,PATCH_COUNT)] for k in range (0,PATCH_COUNT)] for kp in range (0,PATCH_COUNT)] for p in range (-ffreq_xph,ffreq_xph )] for n in range (-ffreq_xph,ffreq_xph)] for m in range (bfreq_xph,bfreq_xph+1)] )

im_f_arr_updo_xph = np.array([[[[[[[[[[f_updo_xph(m,n,p).imag for i in range (0,QN_COUNT)] for j in range (0,QN_COUNT)]for ip in range(0,QN_COUNT)] for jp in range (0,QN_COUNT)] for q in range (0,PATCH_COUNT)] for k in range (0,PATCH_COUNT)] for kp in range (0,PATCH_COUNT)] for p in range (-ffreq_xph,ffreq_xph )] for n in range (-ffreq_xph,ffreq_xph)] for m in range (bfreq_xph,bfreq_xph+1)] )

#----------------------------------------Create HDF5 files-----------------------------------------

f.create_dataset('VERT/XPH/RE_F_UPUP', data=re_f_arr_upup_xph, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('VERT/XPH/IM_F_UPUP', data=im_f_arr_upup_xph, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('VERT/XPH/RE_F_UPDO', data=re_f_arr_updo_xph, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('VERT/XPH/IM_F_UPDO', data=im_f_arr_updo_xph, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('VERT/XPH/fgrid', data=fgrid_arr_xph, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('VERT/XPH/bgrid', data=bgrid_arr_xph, dtype='float64', compression="gzip", compression_opts=4)

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

#chi = np.loadtxt('/home/agnese/Coding/Vertex_structures/dat/pomerol/4SITES/U_'+str(U)+'_BETA_'+str(BETA)+'_FFREQ_'+ str(FFREQ)+'_BFREQ_'+ str(BFREQ)+'/chi_asympt')

bfreq_pp = CHI_BFREQ_ED
print bfreq_pp


def re_chi_updo_pp(wb):
    if(wb >= 0):
        #return U*U*chi[wb, 5]/BETA/BETA
        return 0.0; 
    else:
        #return U*U*chi[-wb, 5]/BETA/BETA 
        return 0.0; 

def im_chi_updo_pp(wb):
    if(wb >= 0):
        #return U*U*chi[wb, 6]/BETA/BETA
        return 0.0; 
    else:
        #return U*U*chi[-wb, 6]/BETA/BETA
        return 0.0; 

def re_chi_updo_ph(wb):
    if(wb >= 0):
        #return U*U*chi[wb, 3]/BETA/BETA
        return 0.0; 
    else:
        #return U*U*chi[-wb, 3]/BETA/BETA 
        return 0.0; 

def im_chi_updo_ph(wb):
    if(wb >= 0):
        #return U*U*chi[wb, 4]/BETA/BETA
        return 0.0; 
    else:
        #return U*U*chi[-wb, 4]/BETA/BETA
        return 0.0; 

def re_chi_upup_ph(wb):
    if(wb >= 0):
        #return U*U*chi[wb, 1]/BETA/BETA
        return 0.0; 
    else:
        #return U*U*chi[-wb, 1]/BETA/BETA 
        return 0.0; 

def im_chi_upup_ph(wb):
    if(wb >= 0):
        #return U*U*chi[wb, 2]/BETA/BETA
        return 0.0; 
    else:
        #return U*U*chi[-wb, 2]/BETA/BETA
        return 0.0; 

# ----------define arrays to store in hdf5 file -> FROM HERE OBJECT IN THE SHIFTED NOTATION!!

re_chi_arr_upup_pp = np.array([[[[[[0 for i in range (0,QN_COUNT)] for j in range (0,QN_COUNT)]for ip in range(0,QN_COUNT)] for jp in range (0,QN_COUNT)] for q in range (0,PATCH_COUNT)] for m in range (-CHI_BFREQ_ED,CHI_BFREQ_ED+1)] )

im_chi_arr_upup_pp = np.array([[[[[[0 for i in range (0,QN_COUNT)] for j in range (0,QN_COUNT)]for ip in range (0,QN_COUNT)] for jp in range (0,QN_COUNT)] for q in range (0,PATCH_COUNT)]for m in range (-CHI_BFREQ_ED,CHI_BFREQ_ED+1)] )

re_chi_arr_updo_pp = np.array([[[[[[re_chi_updo_pp(m) for i in range (0,QN_COUNT)] for j in range (0,QN_COUNT)] for ip in range (0,QN_COUNT)] for jp in range (0,QN_COUNT)] for q in range (0,PATCH_COUNT)] for m in range (-CHI_BFREQ_ED,CHI_BFREQ_ED+1)] )

im_chi_arr_updo_pp = np.array([[[[[[im_chi_updo_pp(m) for i in range (0,QN_COUNT)] for j in range (0,QN_COUNT)]for ip in range (0,QN_COUNT)] for jp in range (0,QN_COUNT)] for q in range (0,PATCH_COUNT)] for m in range (-CHI_BFREQ_ED,CHI_BFREQ_ED+1)] )

bgrid_arr_pp = np.array([(2*n)*pi/BETA for n in range (-CHI_BFREQ_ED,CHI_BFREQ_ED+1)])

# ----------define arrays to store in hdf5 file -> FROM HERE OBJECT IN THE SHIFTED NOTATION!!

re_chi_arr_upup_ph = np.array([[[[[[re_chi_upup_ph(m) for i in range (0,QN_COUNT)] for j in range (0,QN_COUNT)]for ip in range(0,QN_COUNT)] for jp in range (0,QN_COUNT)] for q in range (0,PATCH_COUNT)] for m in range (-CHI_BFREQ_ED,CHI_BFREQ_ED+1)] )

im_chi_arr_upup_ph = np.array([[[[[[im_chi_upup_ph(m) for i in range (0,QN_COUNT)] for j in range (0,QN_COUNT)]for ip in range (0,QN_COUNT)] for jp in range (0,QN_COUNT)] for q in range (0,PATCH_COUNT)]for m in range (-CHI_BFREQ_ED,CHI_BFREQ_ED+1)] )

re_chi_arr_updo_ph = np.array([[[[[[re_chi_updo_ph(m) for i in range (0,QN_COUNT)] for j in range (0,QN_COUNT)] for ip in range (0,QN_COUNT)] for jp in range (0,QN_COUNT)] for q in range (0,PATCH_COUNT)] for m in range (-CHI_BFREQ_ED,CHI_BFREQ_ED+1)] )

im_chi_arr_updo_ph = np.array([[[[[[im_chi_updo_ph(m) for i in range (0,QN_COUNT)] for j in range (0,QN_COUNT)]for ip in range (0,QN_COUNT)] for jp in range (0,QN_COUNT)] for q in range (0,PATCH_COUNT)] for m in range (-CHI_BFREQ_ED,CHI_BFREQ_ED+1)] )

bgrid_arr_ph = np.array([(2*n)*pi/BETA for n in range (-CHI_BFREQ_ED,CHI_BFREQ_ED+1)])

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

#trileg_pp = np.loadtxt('/home/agnese/Coding/Vertex_structures/dat/pomerol/4SITES/U_'+str(U)+'_BETA_'+str(BETA)+'_FFREQ_'+ str(FFREQ)+'_BFREQ_'+ str(BFREQ)+'/trileg_pp.dat')

ffreq_pp = TRI_FFREQ_ED
print ffreq_pp
bfreq_pp = TRI_BFREQ_ED
print bfreq_pp


def re_trileg_updo_pp(wb,wf):
	#return 	-U*(trileg_pp[(2*ffreq_pp)*(wb+bfreq_pp) + (wf+ffreq_pp), 2] - 1) - re_chi_updo_pp(wb)
        return 0.0; 

def re_trileg_upup_pp(wb,wf):
	return 	0.0

def im_trileg_updo_pp(wb,wf):
	#return 	-U*(trileg_pp[(2*ffreq_pp)*(wb+bfreq_pp)+ (wf+ffreq_pp), 3] - 1) - im_chi_updo_pp(wb) 
        return 0.0; 

def im_trileg_upup_pp(wb,wf):
	return 	0.0

# ----------define arrays to store in hdf5 file -> FROM HERE OBJECT IN THE SHIFTED NOTATION!!

re_trileg_arr_upup_pp = np.array([[[[[[[[re_trileg_upup_pp(m,n+myceil_div2(m)) for i in range (0,QN_COUNT)] for j in range (0,QN_COUNT)]for ip in range(0,QN_COUNT)] for jp in range (0,QN_COUNT)] for q in range (0,PATCH_COUNT)] for k in range (0,PATCH_COUNT)] for n in range (-TRI_FFREQ_SHIFT,TRI_FFREQ_SHIFT)] for m in range (-TRI_BFREQ_SHIFT,TRI_BFREQ_SHIFT+1)] )

im_trileg_arr_upup_pp = np.array([[[[[[[[im_trileg_upup_pp(m,n+myceil_div2(m)) for i in range (0,QN_COUNT)] for j in range (0,QN_COUNT)]for ip in range (0,QN_COUNT)] for jp in range (0,QN_COUNT)] for q in range (0,PATCH_COUNT)] for k in range (0,PATCH_COUNT)] for n in range (-TRI_FFREQ_SHIFT,TRI_FFREQ_SHIFT)] for m in range (-TRI_BFREQ_SHIFT,TRI_BFREQ_SHIFT+1)] )

re_trileg_arr_updo_pp = np.array([[[[[[[[re_trileg_updo_pp(m,n+myceil_div2(m)) for i in range (0,QN_COUNT)] for j in range (0,QN_COUNT)] for ip in range (0,QN_COUNT)] for jp in range (0,QN_COUNT)] for q in range (0,PATCH_COUNT)] for k in range (0,PATCH_COUNT)] for n in range (-TRI_FFREQ_SHIFT,TRI_FFREQ_SHIFT)] for m in range (-TRI_BFREQ_SHIFT,TRI_BFREQ_SHIFT+1)] )

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

#trileg_ph = np.loadtxt('/home/agnese/Coding/Vertex_structures/dat/pomerol/4SITES/U_'+str(U)+'_BETA_'+str(BETA)+'_FFREQ_'+ str(FFREQ)+'_BFREQ_'+ str(BFREQ)+'/trileg.dat')

ffreq_ph = TRI_FFREQ_ED
print ffreq_ph
bfreq_ph = TRI_BFREQ_ED
print bfreq_ph


def re_trileg_updo_ph(wb,wf):
	#return 	U*(trileg_ph[(2*ffreq_ph)*(wb+bfreq_ph) + (wf+ffreq_ph), 2] + 1) - re_chi_updo_ph(wb) 
        return 0.0; 

def re_trileg_upup_ph(wb,wf):
	#return 	U*trileg_ph[(2*ffreq_ph)*(wb+bfreq_ph) + (wf+ffreq_ph), 4] - re_chi_upup_ph(wb) 
        return 0.0; 

def im_trileg_updo_ph(wb,wf):
	#return 	U*(trileg_ph[(2*ffreq_ph)*(wb+bfreq_ph)+ (wf+ffreq_ph), 3] + 1) - im_chi_updo_ph(wb)
        return 0.0; 

def im_trileg_upup_ph(wb,wf):
	#return U*trileg_ph[(2*ffreq_ph)*(wb+bfreq_ph)+ (wf+ffreq_ph), 5] - im_chi_upup_ph(wb)
        return 0.0; 

# ----------define arrays to store in hdf5 file -> FROM HERE OBJECT IN THE SHIFTED NOTATION!!

re_trileg_arr_upup_ph = np.array([[[[[[[[re_trileg_upup_ph(m,n-myfloor_div2(m)) for i in range (0,QN_COUNT)] for j in range (0,QN_COUNT)]for ip in range(0,QN_COUNT)] for jp in range (0,QN_COUNT)] for q in range (0,PATCH_COUNT)] for k in range (0,PATCH_COUNT)] for n in range (-TRI_FFREQ_SHIFT,TRI_FFREQ_SHIFT)] for m in range (-TRI_BFREQ_SHIFT,TRI_BFREQ_SHIFT+1)] )

im_trileg_arr_upup_ph = np.array([[[[[[[[im_trileg_upup_ph(m,n-myfloor_div2(m)) for i in range (0,QN_COUNT)] for j in range (0,QN_COUNT)]for ip in range (0,QN_COUNT)] for jp in range (0,QN_COUNT)] for q in range (0,PATCH_COUNT)] for k in range (0,PATCH_COUNT)] for n in range (-TRI_FFREQ_SHIFT,TRI_FFREQ_SHIFT)] for m in range (-TRI_BFREQ_SHIFT,TRI_BFREQ_SHIFT+1)] )

re_trileg_arr_updo_ph = np.array([[[[[[[[re_trileg_updo_ph(m,n-myfloor_div2(m)) for i in range (0,QN_COUNT)] for j in range (0,QN_COUNT)] for ip in range (0,QN_COUNT)] for jp in range (0,QN_COUNT)] for q in range (0,PATCH_COUNT)] for k in range (0,PATCH_COUNT)] for n in range (-TRI_FFREQ_SHIFT,TRI_FFREQ_SHIFT)] for m in range (-TRI_BFREQ_SHIFT,TRI_BFREQ_SHIFT+1)] )

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



