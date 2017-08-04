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
#from agneselibrary.mymath import *                          # mylibrary check in ~/usr/include/agneselib
from _functools import partial

#---------------------------------------------------------------------------------

U=1.0

beta=26.0

FFREQ = 159 #fermionic frequencies in the mixed notation
BFREQ = 160 #bosonic frequency transfer in the mixed notation

pi = math.pi

#---------------------------------------------------------------------------------

def run(command):
        output = subprocess.check_output(command, shell=True)
        return output

#----------------------------------------Create HDF5 files-----------------------------------------

if ('dat/U'+ str(U)+'_beta'+ str(beta)+'_FFREQ_'+ str(FFREQ)+'_BFREQ_' + str(BFREQ)+'.h5'):
    os.system('rm -r dat/U'+ str(U)+'_beta'+ str(beta)+'_FFREQ_'+ str(FFREQ)+'_BFREQ_' + str(BFREQ)+'.h5')

f = h5py.File('dat/U'+ str(U)+'_beta'+ str(beta)+'_FFREQ_'+ str(FFREQ)+'_BFREQ_' + str(BFREQ)+'.h5', 'w')   # Create the datafile.h5

##---------------------------------------------------------------------------------
##GF
#g_iw  = np.loadtxt('dat/pomerol/U_1.0_beta_26.0_FFREQ_20/BFREQ_10_NO_SHIFT/gw_imfreq_00.dat')
#N_fermi_gf = g_iw.shape[0]-20
#print ("Number of fermionic frequencies for the GF:")
#print N_fermi_gf
#
#def G(w):                           # imaginary part of the GF
# if (w >= 0):
#     return g_iw[w+20,2]+1j*g_iw[w+20,3]
# else:
#     return g_iw[-w+20-1,2]-1j*g_iw[-w+20-1,3]
#
#g_arr_re = np.array([G(n).real for n in range (-N_fermi_gf,N_fermi_gf)])
#g_arr_im = np.array([G(n).imag for n in range (-N_fermi_gf,N_fermi_gf)])
#fgrid_arr = np.array([(2*m +1)*pi/beta for m in range (-N_fermi_gf,N_fermi_gf)])
#
##----------------------------------------Create HDF5 files-----------------------------------------
#
#f.create_dataset('Giw/RE', data=g_arr_re, dtype='float64', compression="gzip", compression_opts=4)
#f.create_dataset('Giw/IM', data=g_arr_im, dtype='float64', compression="gzip", compression_opts=4)
#f.create_dataset('Giw/fgrid', data=fgrid_arr, dtype='float64', compression="gzip", compression_opts=4)
#

#---------------------------------------------------------------------------------------------------------
#
#                               PP
#
#---------------------------------------------------------------------------------------------------------
#--------------------------------------- 2PGF PP------------------------------------------------------------

vertex_updo = np.loadtxt('dat/vertex/vert_chi_pp_full')
print vertex_updo.shape
ffreq_pp = int(0.5*(np.transpose(vertex_updo)[1,:].max()*beta/pi-1))+1
print ffreq_pp
bfreq_pp = int(0.5*np.transpose(vertex_updo)[0,:].max()*beta/pi)
print bfreq_pp


def re_f_upup_pp(wb,wf,wf1):
	return 	vertex_updo[(2*ffreq_pp)*(2*ffreq_pp)*(wb+bfreq_pp)+(2*ffreq_pp)*(wf+ffreq_pp) + (wf1+ffreq_pp), 3]

def re_f_updo_pp(wb,wf,wf1):
	return 	vertex_updo[(2*ffreq_pp)*(2*ffreq_pp)*(wb+bfreq_pp)+(2*ffreq_pp)*(wf+ffreq_pp)+ (wf1+ffreq_pp), 5]

def im_f_upup_pp(wb,wf,wf1):

	return 	vertex_updo[(2*ffreq_pp)*(2*ffreq_pp)*(wb+bfreq_pp)+(2*ffreq_pp)*(wf+ffreq_pp) + (wf1+ffreq_pp), 4]
def im_f_updo_pp(wb,wf,wf1):
	return 	vertex_updo[(2*ffreq_pp)*(2*ffreq_pp)*(wb+bfreq_pp)+(2*ffreq_pp)*(wf+ffreq_pp)+ (wf1+ffreq_pp), 6]

# ----------define arrays to store in hdf5 file

re_f_arr_upup_pp = np.array([[[re_f_upup_pp(m,n,p) for p in range (-ffreq_pp,ffreq_pp )] for n in range (-ffreq_pp,ffreq_pp)] for m in range (-bfreq_pp,bfreq_pp+1)] )

im_f_arr_upup_pp = np.array([[[im_f_upup_pp(m,n,p) for p in range (-ffreq_pp,ffreq_pp )] for n in range (-ffreq_pp,ffreq_pp)] for m in range (-bfreq_pp,bfreq_pp+1)] )

re_f_arr_updo_pp = np.array([[[re_f_updo_pp(m,n,p) for p in range (-ffreq_pp,ffreq_pp )] for n in range (-ffreq_pp,ffreq_pp)] for m in range (-bfreq_pp,bfreq_pp+1)] )

im_f_arr_updo_pp = np.array([[[im_f_updo_pp(m,n,p) for p in range (-ffreq_pp,ffreq_pp )] for n in range (-ffreq_pp,ffreq_pp)] for m in range (-bfreq_pp,bfreq_pp+1)] )

fgrid_arr_pp = np.array([(2*m +1)*pi/beta for m in range (-ffreq_pp,ffreq_pp)])
bgrid_arr_pp = np.array([(2*n)*pi/beta for n in range (-bfreq_pp,bfreq_pp+1)])

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

vertex_updo = np.loadtxt('dat/vertex/vert_chi_full')
print vertex_updo.shape
ffreq_pp = int(0.5*(np.transpose(vertex_updo)[1,:].max()*beta/pi-1))+1
print ffreq_pp
bfreq_pp = int(0.5*np.transpose(vertex_updo)[0,:].max()*beta/pi)
print bfreq_pp


def re_f_upup_ph(wb,wf,wf1):
	return 	vertex_updo[(2*ffreq_pp)*(2*ffreq_pp)*(wb+bfreq_pp)+(2*ffreq_pp)*(wf+ffreq_pp) + (wf1+ffreq_pp), 3]

def re_f_updo_ph(wb,wf,wf1):
	return 	vertex_updo[(2*ffreq_pp)*(2*ffreq_pp)*(wb+bfreq_pp)+(2*ffreq_pp)*(wf+ffreq_pp)+ (wf1+ffreq_pp), 5]

def im_f_upup_ph(wb,wf,wf1):

	return 	vertex_updo[(2*ffreq_pp)*(2*ffreq_pp)*(wb+bfreq_pp)+(2*ffreq_pp)*(wf+ffreq_pp) + (wf1+ffreq_pp), 4]

def im_f_updo_ph(wb,wf,wf1):
	return 	vertex_updo[(2*ffreq_pp)*(2*ffreq_pp)*(wb+bfreq_pp)+(2*ffreq_pp)*(wf+ffreq_pp)+ (wf1+ffreq_pp), 6]

# ----------define arrays to store in hdf5 file

re_f_arr_upup_ph = np.array([[[re_f_upup_ph(m,n,p) for p in range (-ffreq_pp,ffreq_pp )] for n in range (-ffreq_pp,ffreq_pp)] for m in range (-bfreq_pp,bfreq_pp+1)] )

im_f_arr_upup_ph = np.array([[[im_f_upup_ph(m,n,p) for p in range (-ffreq_pp,ffreq_pp )] for n in range (-ffreq_pp,ffreq_pp)] for m in range (-bfreq_pp,bfreq_pp+1)] )

re_f_arr_updo_ph = np.array([[[re_f_updo_ph(m,n,p) for p in range (-ffreq_pp,ffreq_pp )] for n in range (-ffreq_pp,ffreq_pp)] for m in range (-bfreq_pp,bfreq_pp+1)] )

im_f_arr_updo_ph = np.array([[[im_f_updo_ph(m,n,p) for p in range (-ffreq_pp,ffreq_pp )] for n in range (-ffreq_pp,ffreq_pp)] for m in range (-bfreq_pp,bfreq_pp+1)] )

fgrid_arr_pp = np.array([(2*m +1)*pi/beta for m in range (-ffreq_pp,ffreq_pp)])
bgrid_arr_pp = np.array([(2*n)*pi/beta for n in range (-bfreq_pp,bfreq_pp+1)])

#----------------------------------------Create HDF5 files-----------------------------------------

f.create_dataset('VERT/PH/RE_F_UPUP', data=re_f_arr_upup_ph, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('VERT/PH/IM_F_UPUP', data=im_f_arr_upup_ph, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('VERT/PH/RE_F_UPDO', data=re_f_arr_updo_ph, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('VERT/PH/IM_F_UPDO', data=im_f_arr_updo_ph, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('VERT/PH/fgrid', data=fgrid_arr_pp, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('VERT/PH/bgrid', data=bgrid_arr_pp, dtype='float64', compression="gzip", compression_opts=4)

#---------------------------------------------------------------------------------------------------------
#
#                               XPH
#
#---------------------------------------------------------------------------------------------------------
#--------------------------------------- 2PGF XPH------------------------------------------------------------

def re_f_upup_xph(wb,wf,wf1):
	return 	-vertex_updo[(2*ffreq_pp)*(2*ffreq_pp)*(wb+bfreq_pp)+(2*ffreq_pp)*(wf+ffreq_pp) + (wf1+ffreq_pp), 3]

def re_f_updo_xph(wb,wf,wf1):
	return 	-vertex_updo[(2*ffreq_pp)*(2*ffreq_pp)*(wb+bfreq_pp)+(2*ffreq_pp)*(wf+ffreq_pp)+ (wf1+ffreq_pp), 3]+vertex_updo[(2*ffreq_pp)*(2*ffreq_pp)*(wb+bfreq_pp)+(2*ffreq_pp)*(wf+ffreq_pp) + (wf1+ffreq_pp), 5]

def im_f_upup_xph(wb,wf,wf1):

	return 	-vertex_updo[(2*ffreq_pp)*(2*ffreq_pp)*(wb+bfreq_pp)+(2*ffreq_pp)*(wf+ffreq_pp) + (wf1+ffreq_pp), 4]

def im_f_updo_xph(wb,wf,wf1):
	return 	-vertex_updo[(2*ffreq_pp)*(2*ffreq_pp)*(wb+bfreq_pp)+(2*ffreq_pp)*(wf+ffreq_pp)+ (wf1+ffreq_pp), 4]+vertex_updo[(2*ffreq_pp)*(2*ffreq_pp)*(wb+bfreq_pp)+(2*ffreq_pp)*(wf+ffreq_pp) + (wf1+ffreq_pp), 6]

# ----------define arrays to store in hdf5 file

re_f_arr_upup_xph = np.array([[[re_f_upup_xph(m,n,p) for p in range (-ffreq_pp,ffreq_pp )] for n in range (-ffreq_pp,ffreq_pp)] for m in range (-bfreq_pp,bfreq_pp+1)] )

im_f_arr_upup_xph = np.array([[[im_f_upup_xph(m,n,p) for p in range (-ffreq_pp,ffreq_pp )] for n in range (-ffreq_pp,ffreq_pp)] for m in range (-bfreq_pp,bfreq_pp+1)] )

re_f_arr_updo_xph = np.array([[[re_f_updo_xph(m,n,p) for p in range (-ffreq_pp,ffreq_pp )] for n in range (-ffreq_pp,ffreq_pp)] for m in range (-bfreq_pp,bfreq_pp+1)] )

im_f_arr_updo_xph = np.array([[[im_f_updo_xph(m,n,p) for p in range (-ffreq_pp,ffreq_pp )] for n in range (-ffreq_pp,ffreq_pp)] for m in range (-bfreq_pp,bfreq_pp+1)] )

fgrid_arr_pp = np.array([(2*m +1)*pi/beta for m in range (-ffreq_pp,ffreq_pp)])
bgrid_arr_pp = np.array([(2*n)*pi/beta for n in range (-bfreq_pp,bfreq_pp+1)])

#----------------------------------------Create HDF5 files-----------------------------------------

f.create_dataset('VERT/XPH/RE_F_UPUP', data=re_f_arr_upup_xph, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('VERT/XPH/IM_F_UPUP', data=im_f_arr_upup_xph, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('VERT/XPH/RE_F_UPDO', data=re_f_arr_updo_xph, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('VERT/XPH/IM_F_UPDO', data=im_f_arr_updo_xph, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('VERT/XPH/fgrid', data=fgrid_arr_pp, dtype='float64', compression="gzip", compression_opts=4)
f.create_dataset('VERT/XPH/bgrid', data=bgrid_arr_pp, dtype='float64', compression="gzip", compression_opts=4)



f.close()

