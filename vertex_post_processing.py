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
from agneselib.mymath import *                          # mylibrary check in ~/usr/include/agneselib
from agneselib.inf_Msum import *
from _functools import partial

#---------------------------------------------------------------------------------

pi = math.pi

#---------------------------------------------------------------------------------

def run(command):
        output = subprocess.check_output(command, shell=True)
        return output

#-----------------------------------Read parameters----------------------------------------------
# ------U---------

a=raw_input('Enter the interaction value:') 
try:
    U=float(a)
except ValueError:
    sys.exit("Invalid interaction")

print ("Uhub value     " + str(U) )

#-------BETA-------

a = raw_input('Enter the value of beta:') 
try:
    beta=float(a) 
except ValueError:
    sys.exit("Invalid beta")

print ("beta value     " + str(beta) )

#----------------------------------------Create HDF5 files-----------------------------------------

if ('dat/dat_U'+ str(U)+'_beta'+ str(beta)+'_EDpomerol.h5'):
    os.system('rm dat/dat_U'+ str(U)+'_beta'+ str(beta)+'_EDpomerol.h5')

f = h5py.File('dat/dat_U'+ str(U)+'_beta'+ str(beta)+'_EDpomerol.h5', 'w')   # Create the datafile.h5

#---------------------------------------------------------------------------------
#GF

g_iw  = np.loadtxt('pomerol_output/U_0.25_beta_10.0_nbs2_32Wpos_80wfpos/gw_imag00.dat')
N_fermi_gf = g_iw.shape[0]
print ("Number of fermionic frequencies for the GF:")
print N_fermi_gf

def G(w):                           # imaginary part of the GF
 if (w >= 0):
     return g_iw[w,1]+1j*g_iw[w,2]
 else:
     return g_iw[-w-1,1]-1j*g_iw[-w-1,2]

g_arr_re = np.array([G(n).real for n in range (-N_fermi_gf,N_fermi_gf)])
g_arr_im = np.array([G(n).imag for n in range (-N_fermi_gf,N_fermi_gf)])
fgrid_arr = np.array([(2*m +1)*pi/beta for m in range (-N_fermi_gf,N_fermi_gf)])

#----------------------------------------Create HDF5 files-----------------------------------------

f.create_dataset('Giw/RE', data=g_arr_re, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('Giw/IM', data=g_arr_im, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('Giw/fgrid', data=fgrid_arr, dtype='float64', compression="gzip", compression_opts=4)


#---------------------------------------------------------------------------------------------------------
#
#                               PP
#
#---------------------------------------------------------------------------------------------------------
#--------------------------------------- 2PGF PP------------------------------------------------------------

vertex_pp = np.loadtxt("pomerol_output/U_0.25_beta_10.0_nbs2_32Wpos_80wfpos/2pgf_pp_shift.dat")
print vertex_pp.shape
ffreq_pp = int(np.transpose(vertex_pp)[1,:].max()+1)
print ffreq_pp
bfreq_pp = int(np.transpose(vertex_pp)[0,:].max())
print bfreq_pp


def re_2pgf_upup_pp(wb,wf,wf1):
	return 	vertex_pp[(2*ffreq_pp)*(2*ffreq_pp)*(wb+bfreq_pp)+(2*ffreq_pp)*(wf+ffreq_pp) + (wf1+ffreq_pp), 3]

def re_2pgf_updo_pp(wb,wf,wf1):
	return 	vertex_pp[(2*ffreq_pp)*(2*ffreq_pp)*(wb+bfreq_pp)+(2*ffreq_pp)*(wf+ffreq_pp)+ (wf1+ffreq_pp), 5]

def im_2pgf_upup_pp(wb,wf,wf1):
	return 	vertex_pp[(2*ffreq_pp)*(2*ffreq_pp)*(wb+bfreq_pp)+(2*ffreq_pp)*(wf+ffreq_pp) + (wf1+ffreq_pp), 4]

def im_2pgf_updo_pp(wb,wf,wf1):
	return 	vertex_pp[(2*ffreq_pp)*(2*ffreq_pp)*(wb+bfreq_pp)+(2*ffreq_pp)*(wf+ffreq_pp)+ (wf1+ffreq_pp), 6]

# ----------define arrays to store in hdf5 file

re_2pgf_arr_upup_pp = np.array([[[re_2pgf_upup_pp(m,n,p) for p in range (-ffreq_pp,ffreq_pp )] for n in range (-ffreq_pp,ffreq_pp)] for m in range (-bfreq_pp,bfreq_pp+1)] )

im_2pgf_arr_upup_pp = np.array([[[im_2pgf_upup_pp(m,n,p) for p in range (-ffreq_pp,ffreq_pp )] for n in range (-ffreq_pp,ffreq_pp)] for m in range (-bfreq_pp,bfreq_pp+1)] )

re_2pgf_arr_updo_pp = np.array([[[re_2pgf_updo_pp(m,n,p) for p in range (-ffreq_pp,ffreq_pp )] for n in range (-ffreq_pp,ffreq_pp)] for m in range (-bfreq_pp,bfreq_pp+1)] )

im_2pgf_arr_updo_pp = np.array([[[im_2pgf_updo_pp(m,n,p) for p in range (-ffreq_pp,ffreq_pp )] for n in range (-ffreq_pp,ffreq_pp)] for m in range (-bfreq_pp,bfreq_pp+1)] )

fgrid_arr_pp = np.array([(2*m +1)*pi/beta for m in range (-ffreq_pp,ffreq_pp)])
bgrid_arr_pp = np.array([(2*n)*pi/beta for n in range (-bfreq_pp,bfreq_pp+1)])

#----------------------------------------Create HDF5 files-----------------------------------------

f.create_dataset('2PGF/PP/RE_2PGF_UPUP', data=re_2pgf_arr_upup_pp, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('2PGF/PP/IM_2PGF_UPUP', data=im_2pgf_arr_upup_pp, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('2PGF/PP/RE_2PGF_UPDO', data=re_2pgf_arr_updo_pp, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('2PGF/PP/IM_2PGF_UPDO', data=im_2pgf_arr_updo_pp, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('2PGF/PP/fgrid', data=fgrid_arr_pp, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('2PGF/PP/bgrid', data=bgrid_arr_pp, dtype='float64', compression="gzip", compression_opts=4)

#------------------------------------- GENERALIZED SUSCEPTIBILITY CHI PP --------------------------------------------

# Cutting one disconnected diagram the crossed one -> generalized susceptibility


def chi_upup_pp(wb,wf,wf1):
    if (wf == -wf1-1-mymod_abs(wb)):
        return re_2pgf_upup_pp(wb,wf,wf1)+1j*im_2pgf_upup_pp(wb,wf,wf1)+beta*(G(wf+myceil_div2(wb))*G(myfloor_div2(wb)-wf-1))
    else:
        return re_2pgf_upup_pp(wb,wf,wf1)+1j*im_2pgf_upup_pp(wb,wf,wf1)

def chi_updo_pp(wb,wf,wf1):
    return re_2pgf_updo_pp(wb,wf,wf1)+1j*im_2pgf_updo_pp(wb,wf,wf1)

def chi_s_pp(wb,wf,wf1):
    return -chi_upup_pp(wb,wf,wf1) + 2*chi_updo_pp(wb,wf,wf1)

def chi_t_pp(wb,wf,wf1):
    return chi_upup_pp(wb,wf,wf1)


# ----------define arrays to store in hdf5 file

re_chi_arr_upup_pp = np.array([[[chi_upup_pp(m,n,p).real for p in range (-ffreq_pp,ffreq_pp )] for n in range (-ffreq_pp,ffreq_pp)] for m in range (-bfreq_pp,bfreq_pp+1)] )

im_chi_arr_upup_pp = np.array([[[chi_upup_pp(m,n,p).imag for p in range (-ffreq_pp,ffreq_pp )] for n in range (-ffreq_pp,ffreq_pp)] for m in range (-bfreq_pp,bfreq_pp+1)] )

re_chi_arr_updo_pp = np.array([[[chi_updo_pp(m,n,p).real for p in range (-ffreq_pp,ffreq_pp )] for n in range (-ffreq_pp,ffreq_pp)] for m in range (-bfreq_pp,bfreq_pp+1)] )

im_chi_arr_updo_pp = np.array([[[chi_updo_pp(m,n,p).imag for p in range (-ffreq_pp,ffreq_pp )] for n in range (-ffreq_pp,ffreq_pp)] for m in range (-bfreq_pp,bfreq_pp+1)] )

#----------------------------------------Create HDF5 files-----------------------------------------

f.create_dataset('GENCHI/PP/RE_GENCHI_UPUP', data=re_chi_arr_upup_pp, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('GENCHI/PP/IM_GENCHI_UPUP', data=im_chi_arr_upup_pp, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('GENCHI/PP/RE_GENCHI_UPDO', data=re_chi_arr_updo_pp, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('GENCHI/PP/IM_GENCHI_UPDO', data=im_chi_arr_updo_pp, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('GENCHI/PP/fgrid', data=fgrid_arr_pp, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('GENCHI/PP/bgrid', data=bgrid_arr_pp, dtype='float64', compression="gzip", compression_opts=4)


#----------------------------------- FULL VERTEX PP -----------------------------------------------------------

def chi_l_pp(wb,wf,wf1):  # Legs on one side of the diagram
    return beta*G(wf+myceil_div2(wb))*G(myfloor_div2(wb)-wf1-1)

def chi_0_pp(wb,wf,wf1):
    if (wf == wf1):
        return beta*G(wf+myceil_div2(wb))*G(myfloor_div2(wb)-wf-1)
    else:
        return 0.0

def f_upup_pp(wb,wf,wf1):
    return beta*beta*(1.0/(chi_l_pp(wb,wf,wf)))*(chi_upup_pp(wb,wf,wf1)-chi_0_pp(wb,wf,wf1))*(1.0/(chi_l_pp(wb,wf1,wf1)))

def f_updo_pp(wb,wf,wf1):
    return beta*beta*(1.0/(chi_l_pp(wb,wf,wf)))*(chi_updo_pp(wb,wf,wf1)-chi_0_pp(wb,wf,wf1))*(1.0/(chi_l_pp(wb,wf1,wf1)))

# ----------define arrays to store in hdf5 file

re_f_arr_upup_pp = np.array([[[f_upup_pp(m,n,p).real for p in range (-ffreq_pp,ffreq_pp )] for n in range (-ffreq_pp,ffreq_pp)] for m in range (-bfreq_pp,bfreq_pp+1)] )

im_f_arr_upup_pp = np.array([[[f_upup_pp(m,n,p).imag for p in range (-ffreq_pp,ffreq_pp )] for n in range (-ffreq_pp,ffreq_pp)] for m in range (-bfreq_pp,bfreq_pp+1)] )

re_f_arr_updo_pp = np.array([[[f_updo_pp(m,n,p).real for p in range (-ffreq_pp,ffreq_pp )] for n in range (-ffreq_pp,ffreq_pp)] for m in range (-bfreq_pp,bfreq_pp+1)] )

im_f_arr_updo_pp = np.array([[[f_updo_pp(m,n,p).imag for p in range (-ffreq_pp,ffreq_pp )] for n in range (-ffreq_pp,ffreq_pp)] for m in range (-bfreq_pp,bfreq_pp+1)] )

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

vertex_ph = np.loadtxt("pomerol_output/U_0.25_beta_10.0_nbs2_32Wpos_80wfpos/2pgf_ph_shift.dat")
ffreq_ph = int(np.transpose(vertex_ph)[1,:].max()+1)
print ffreq_ph
bfreq_ph = int(np.transpose(vertex_ph)[0,:].max())
print bfreq_ph

def re_2pgf_upup_ph(wb,wf,wf1):
	return 	vertex_ph[(2*ffreq_ph)*(2*ffreq_ph)*(wb+bfreq_ph)+(2*ffreq_ph)*(wf+ffreq_ph) + (wf1+ffreq_ph), 3]

def re_2pgf_updo_ph(wb,wf,wf1):
	return 	vertex_ph[(2*ffreq_ph)*(2*ffreq_ph)*(wb+bfreq_ph)+(2*ffreq_ph)*(wf+ffreq_ph)+ (wf1+ffreq_ph), 5]

def im_2pgf_upup_ph(wb,wf,wf1):
	return 	vertex_ph[(2*ffreq_ph)*(2*ffreq_ph)*(wb+bfreq_ph)+(2*ffreq_ph)*(wf+ffreq_ph) + (wf1+ffreq_ph), 4]

def im_2pgf_updo_ph(wb,wf,wf1):
	return 	vertex_ph[(2*ffreq_ph)*(2*ffreq_ph)*(wb+bfreq_ph)+(2*ffreq_ph)*(wf+ffreq_ph)+ (wf1+ffreq_ph), 6]

# ----------define arrays to store in hdf5 file

re_2pgf_arr_upup_ph = np.array([[[re_2pgf_upup_ph(m,n,p) for p in range (-ffreq_ph,ffreq_ph )] for n in range (-ffreq_ph,ffreq_ph)] for m in range (-bfreq_ph,bfreq_ph+1)] )

im_2pgf_arr_upup_ph = np.array([[[im_2pgf_upup_ph(m,n,p) for p in range (-ffreq_ph,ffreq_ph )] for n in range (-ffreq_ph,ffreq_ph)] for m in range (-bfreq_ph,bfreq_ph+1)] )

re_2pgf_arr_updo_ph = np.array([[[re_2pgf_updo_ph(m,n,p) for p in range (-ffreq_ph,ffreq_ph )] for n in range (-ffreq_ph,ffreq_ph)] for m in range (-bfreq_ph,bfreq_ph+1)] )

im_2pgf_arr_updo_ph = np.array([[[im_2pgf_updo_ph(m,n,p) for p in range (-ffreq_ph,ffreq_ph )] for n in range (-ffreq_ph,ffreq_ph)] for m in range (-bfreq_ph,bfreq_ph+1)] )

fgrid_arr_ph = np.array([(2*m +1)*pi/beta for m in range (-ffreq_ph,ffreq_ph)])
bgrid_arr_ph = np.array([(2*n)*pi/beta for n in range (-bfreq_ph,bfreq_ph+1)])

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
    if (wf == wf1):
        return re_2pgf_upup_ph(wb,wf,wf1)+1j*im_2pgf_upup_ph(wb,wf,wf1)+beta*(G(wf-myfloor_div2(wb))*G(wf1+myceil_div2(wb)))
    else:
        return re_2pgf_upup_ph(wb,wf,wf1)+1j*im_2pgf_upup_ph(wb,wf,wf1)

def chi_updo_ph(wb,wf,wf1):
    return re_2pgf_updo_ph(wb,wf,wf1)+1j*im_2pgf_updo_ph(wb,wf,wf1)

def tpgf_d_ph(wb,wf,wf1):
    return re_2pgf_upup_ph(wb,wf,wf1)+1j*im_2pgf_upup_ph(wb,wf,wf1)+re_2pgf_updo_ph(wb,wf,wf1)+1j*im_2pgf_updo_ph(wb,wf,wf1)

def tpgf_m_ph(wb,wf,wf1):
    return re_2pgf_upup_ph(wb,wf,wf1)+1j*im_2pgf_upup_ph(wb,wf,wf1) -re_2pgf_updo_ph(wb,wf,wf1)-1j*im_2pgf_updo_ph(wb,wf,wf1)

# ----------define arrays to store in hdf5 file

re_chi_arr_upup_ph = np.array([[[chi_upup_ph(m,n,p).real for p in range (-ffreq_ph,ffreq_ph )] for n in range (-ffreq_ph,ffreq_ph)] for m in range (-bfreq_ph,bfreq_ph+1)] )

im_chi_arr_upup_ph = np.array([[[chi_upup_ph(m,n,p).imag for p in range (-ffreq_ph,ffreq_ph )] for n in range (-ffreq_ph,ffreq_ph)] for m in range (-bfreq_ph,bfreq_ph+1)] )

re_chi_arr_updo_ph = np.array([[[chi_updo_ph(m,n,p).real for p in range (-ffreq_ph,ffreq_ph )] for n in range (-ffreq_ph,ffreq_ph)] for m in range (-bfreq_ph,bfreq_ph+1)] )

im_chi_arr_updo_ph = np.array([[[chi_updo_ph(m,n,p).imag for p in range (-ffreq_ph,ffreq_ph )] for n in range (-ffreq_ph,ffreq_ph)] for m in range (-bfreq_ph,bfreq_ph+1)] )

#----------------------------------------Create HDF5 files-----------------------------------------

f.create_dataset('GENCHI/PH/RE_GENCHI_UPUP', data=re_chi_arr_upup_ph, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('GENCHI/PH/IM_GENCHI_UPUP', data=im_chi_arr_upup_ph, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('GENCHI/PH/RE_GENCHI_UPDO', data=re_chi_arr_updo_ph, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('GENCHI/PH/IM_GENCHI_UPDO', data=im_chi_arr_updo_ph, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('GENCHI/PH/fgrid', data=fgrid_arr_ph, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('GENCHI/PH/bgrid', data=bgrid_arr_ph, dtype='float64', compression="gzip", compression_opts=4)

#----------------------------------- FULL VERTEX PH -----------------------------------------------------------

def chi_l_ph(wb,wf,wf1):
    return G(wf-myfloor_div2(wb))*G(wf1+myceil_div2(wb))

def chi_0_ph(wb,wf,wf1):
    if (wb == 0):
        return beta*G(wf)*G(wf1)
    else:
        return 0.0

def chi_x0_ph(wb,wf,wf1):
    if (wf == wf1):
        return -beta*G(wf-myfloor_div2(wb))*G(wf1+myceil_div2(wb))
    else:
        return 0.0

def f_upup_ph(wb,wf,wf1):
    return (1.0/chi_l_ph(wb,wf,wf1))*(chi_upup_ph(wb,wf,wf1)-chi_0_ph(wb,wf,wf1))*(1.0/chi_l_ph(wb,wf1,wf))

def f_updo_ph(wb,wf,wf1):
    return (1.0/chi_l_ph(wb,wf,wf1))*(chi_updo_ph(wb,wf,wf1)-chi_0_ph(wb,wf,wf1))*(1.0/chi_l_ph(wb,wf1,wf))


# ----------define arrays to store in hdf5 file

re_f_arr_upup_ph = np.array([[[f_upup_ph(m,n,p).real for p in range (-ffreq_ph,ffreq_ph )] for n in range (-ffreq_ph,ffreq_ph)] for m in range (-bfreq_ph,bfreq_ph+1)] )

im_f_arr_upup_ph = np.array([[[f_upup_ph(m,n,p).imag for p in range (-ffreq_ph,ffreq_ph )] for n in range (-ffreq_ph,ffreq_ph)] for m in range (-bfreq_ph,bfreq_ph+1)] )

re_f_arr_updo_ph = np.array([[[f_updo_ph(m,n,p).real for p in range (-ffreq_ph,ffreq_ph )] for n in range (-ffreq_ph,ffreq_ph)] for m in range (-bfreq_ph,bfreq_ph+1)] )

im_f_arr_updo_ph = np.array([[[f_updo_ph(m,n,p).imag for p in range (-ffreq_ph,ffreq_ph )] for n in range (-ffreq_ph,ffreq_ph)] for m in range (-bfreq_ph,bfreq_ph+1)] )

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

vertex_xph = np.loadtxt("pomerol_output/U_0.25_beta_10.0_nbs2_32Wpos_80wfpos/2pgf_xph_shift.dat")
ffreq_xph = int(np.transpose(vertex_xph)[1,:].max()+1)
print ffreq_xph
bfreq_xph = int(np.transpose(vertex_xph)[0,:].max())
print bfreq_xph

def re_2pgf_upup_xph(wb,wf,wf1):
	return 	vertex_xph[(2*ffreq_xph)*(2*ffreq_xph)*(wb+bfreq_xph)+(2*ffreq_xph)*(wf+ffreq_xph) + (wf1+ffreq_xph), 3]

def re_2pgf_updo_xph(wb,wf,wf1):
	return 	vertex_xph[(2*ffreq_xph)*(2*ffreq_xph)*(wb+bfreq_xph)+(2*ffreq_xph)*(wf+ffreq_xph)+ (wf1+ffreq_xph), 5]

def im_2pgf_upup_xph(wb,wf,wf1):
	return 	vertex_xph[(2*ffreq_xph)*(2*ffreq_xph)*(wb+bfreq_xph)+(2*ffreq_xph)*(wf+ffreq_xph) + (wf1+ffreq_xph), 4]

def im_2pgf_updo_xph(wb,wf,wf1):
	return 	vertex_xph[(2*ffreq_xph)*(2*ffreq_xph)*(wb+bfreq_xph)+(2*ffreq_xph)*(wf+ffreq_xph)+ (wf1+ffreq_xph), 6]

# ----------define arrays to store in hdf5 file

re_2pgf_arr_upup_xph = np.array([[[re_2pgf_upup_xph(m,n,p) for p in range (-ffreq_xph,ffreq_xph )] for n in range (-ffreq_xph,ffreq_xph)] for m in range (-bfreq_xph,bfreq_xph+1)] )

im_2pgf_arr_upup_xph = np.array([[[im_2pgf_upup_xph(m,n,p) for p in range (-ffreq_xph,ffreq_xph )] for n in range (-ffreq_xph,ffreq_xph)] for m in range (-bfreq_xph,bfreq_xph+1)] )

re_2pgf_arr_updo_xph = np.array([[[re_2pgf_updo_xph(m,n,p) for p in range (-ffreq_xph,ffreq_xph )] for n in range (-ffreq_xph,ffreq_xph)] for m in range (-bfreq_xph,bfreq_xph+1)] )

im_2pgf_arr_updo_xph = np.array([[[im_2pgf_updo_xph(m,n,p) for p in range (-ffreq_xph,ffreq_xph )] for n in range (-ffreq_xph,ffreq_xph)] for m in range (-bfreq_xph,bfreq_xph+1)] )

fgrid_arr_xph = np.array([(2*m +1)*pi/beta for m in range (-ffreq_xph,ffreq_xph)])
bgrid_arr_xph = np.array([(2*n)*pi/beta for n in range (-bfreq_xph,bfreq_xph+1)])

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
        return re_2pgf_upup_xph(wb,wf,wf1)+1j*im_2pgf_upup_xph(wb,wf,wf1)+beta*G(wf)*G(wf1)
    else:
        return re_2pgf_upup_xph(wb,wf,wf1)+1j*im_2pgf_upup_xph(wb,wf,wf1)

def chi_updo_xph(wb,wf,wf1):
    return re_2pgf_updo_xph(wb,wf,wf1)+1j*im_2pgf_updo_xph(wb,wf,wf1)

# ----------define arrays to store in hdf5 file

re_chi_arr_upup_xph = np.array([[[chi_upup_xph(m,n,p).real for p in range (-ffreq_xph,ffreq_xph )] for n in range (-ffreq_xph,ffreq_xph)] for m in range (-bfreq_xph,bfreq_xph+1)] )

im_chi_arr_upup_xph = np.array([[[chi_upup_xph(m,n,p).imag for p in range (-ffreq_xph,ffreq_xph )] for n in range (-ffreq_xph,ffreq_xph)] for m in range (-bfreq_xph,bfreq_xph+1)] )

re_chi_arr_updo_xph = np.array([[[chi_updo_xph(m,n,p).real for p in range (-ffreq_xph,ffreq_xph )] for n in range (-ffreq_xph,ffreq_xph)] for m in range (-bfreq_xph,bfreq_xph+1)] )

im_chi_arr_updo_xph = np.array([[[chi_updo_xph(m,n,p).imag for p in range (-ffreq_xph,ffreq_xph )] for n in range (-ffreq_xph,ffreq_xph)] for m in range (-bfreq_xph,bfreq_xph+1)] )

#----------------------------------------Create HDF5 files-----------------------------------------
f.create_dataset('GENCHI/XPH/RE_GENCHI_UPUP', data=re_chi_arr_upup_xph, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('GENCHI/XPH/IM_GENCHI_UPUP', data=im_chi_arr_upup_xph, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('GENCHI/XPH/RE_GENCHI_UPDO', data=re_chi_arr_updo_xph, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('GENCHI/XPH/IM_GENCHI_UPDO', data=im_chi_arr_updo_xph, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('GENCHI/XPH/fgrid', data=fgrid_arr_xph, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('GENCHI/XPH/bgrid', data=bgrid_arr_xph, dtype='float64', compression="gzip", compression_opts=4)

#----------------------------------- FULL VERTEX XPH -----------------------------------------------------------

def chi_l_xph(wb,wf,wf1):
    return G(wf-myfloor_div2(wb))*G(wf1+myceil_div2(wb))

def chi_0_xph(wb,wf,wf1):
    if (wf == wf1):
        return beta*G(wf-myfloor_div2(wb))*G(wf+myceil_div2(wb))
    else:
        return 0.0

def chi_x0_xph(wb,wf,wf1):
    if (wb == 0):
        return -beta*G(wf)*G(wf1)
    else:
        return 0.0

def f_upup_xph(wb,wf,wf1):
    return (1.0/chi_l_xph(wb,wf,wf1))*(chi_upup_xph(wb,wf,wf1)-chi_0_xph(wb,wf,wf1))*(1.0/chi_l_xph(wb,wf1,wf))

def f_updo_xph(wb,wf,wf1):
    return (1.0/chi_l_xph(wb,wf,wf1))*(chi_updo_xph(wb,wf,wf1)-chi_0_xph(wb,wf,wf1))*(1.0/chi_l_xph(wb,wf1,wf))


# ----------define arrays to store in hdf5 file

re_f_arr_upup_xph = np.array([[[f_upup_xph(m,n,p).real for p in range (-ffreq_xph,ffreq_xph )] for n in range (-ffreq_xph,ffreq_xph)] for m in range (-bfreq_xph,bfreq_xph+1)] )

im_f_arr_upup_xph = np.array([[[f_upup_xph(m,n,p).imag for p in range (-ffreq_xph,ffreq_xph )] for n in range (-ffreq_xph,ffreq_xph)] for m in range (-bfreq_xph,bfreq_xph+1)] )

re_f_arr_updo_xph = np.array([[[f_updo_xph(m,n,p).real for p in range (-ffreq_xph,ffreq_xph )] for n in range (-ffreq_xph,ffreq_xph)] for m in range (-bfreq_xph,bfreq_xph+1)] )

im_f_arr_updo_xph = np.array([[[f_updo_xph(m,n,p).imag for p in range (-ffreq_xph,ffreq_xph )] for n in range (-ffreq_xph,ffreq_xph)] for m in range (-bfreq_xph,bfreq_xph+1)] )

#----------------------------------------Create HDF5 files-----------------------------------------

f.create_dataset('VERT/XPH/RE_F_UPUP', data=re_f_arr_upup_xph, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('VERT/XPH/IM_F_UPUP', data=im_f_arr_upup_xph, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('VERT/XPH/RE_F_UPDO', data=re_f_arr_updo_xph, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('VERT/XPH/IM_F_UPDO', data=im_f_arr_updo_xph, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('VERT/XPH/fgrid', data=fgrid_arr_xph, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('VERT/XPH/bgrid', data=bgrid_arr_xph, dtype='float64', compression="gzip", compression_opts=4)


f.close()

