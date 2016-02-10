#!/usr/bin/python
#
#---------------------------------------------------------------------------------------
#
#NOTE:              This script plots Plus and Karrsch functions after the self-consistency
#                   No plots saved
#
#--------------------------------------IMPORTS ------------------------------------------

import h5py
import matplotlib.pyplot as pl
import numpy as np
import os
import sys
import subprocess
import math

#--------------------------------------MANAGING FILES ------------------------------------------

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

#----------------------------------------Read HDF5 files-----------------------------------------

if ('../dat'):
    f = h5py.File('../dat/dat_U'+ str(U)+'_beta'+ str(beta)+'_EDpomerol.h5', 'r')   # Read (and write) the hdf5 file in the directory "dat" if existing
else:
    sys.exit("No data file")

#---------  PARAMETERS  ---------------

pi = math.pi
shift = 0
#--------------------------------------VERTEX PLOTTING ------------------------------------------

print("Plotting P ...")


#--- Read
rep_upup_ph =  np.array(f["/P_func/PH/RE_P_UPUP"])
imp_upup_ph = np.array(f["/P_func/PH/IM_P_UPUP"])
rep_updo_ph =  np.array(f["/P_func/PH/RE_P_UPDO"])
imp_updo_ph = np.array(f["/P_func/PH/IM_P_UPDO"])
rep_updo_pp =  np.array(f["/P_func/PP/RE_P_UPDO"])
imp_updo_pp = np.array(f["/P_func/PP/IM_P_UPDO"])
rep_updo_xph =  np.array(f["/P_func/XPH/RE_P_UPDO"])
imp_updo_xph = np.array(f["/P_func/XPH/IM_P_UPDO"])
rep_upup_xph =  np.array(f["/P_func/XPH/RE_P_UPUP"])
imp_upup_xph = np.array(f["/P_func/XPH/IM_P_UPUP"])

#----Read Karrasch

rek_upup_ph =  np.array(f["/K_func/PH/RE_K_UPUP"])
imk_upup_ph = np.array(f["/K_func/PH/IM_K_UPUP"])
rek_updo_ph =  np.array(f["/K_func/PH/RE_K_UPDO"])
imk_updo_ph = np.array(f["/K_func/PH/IM_K_UPDO"])
rek_updo_pp =  np.array(f["/K_func/PP/RE_K_UPDO"])
imk_updo_pp = np.array(f["/K_func/PP/IM_K_UPDO"])
rek_updo_xph =  np.array(f["/K_func/XPH/RE_K_UPDO"])
imk_updo_xph = np.array(f["/K_func/XPH/IM_K_UPDO"])
rek_upup_xph =  np.array(f["/K_func/XPH/RE_K_UPUP"])
imk_upup_xph = np.array(f["/K_func/XPH/IM_K_UPUP"])

fdim = rep_upup_ph.shape[0]
fdim_pp = rep_updo_pp.shape[0]

if fdim <= shift:
    sys.exit("Error: Shift too large for vertex grid"); 

if fdim_pp <= shift:
    sys.exit("Error: Shift too large for the pp vertex grid"); 

#---  Helper functions --------------------------------------- 

def neg( w ):
    return fdim - w - 1

N_bose_p = (rep_upup_ph.shape[1]-1)/2
N_fermi_p = (rep_upup_ph.shape[0])/2


def plotP( use_pl, zarr, string ):
#    use_pl.set_aspect(1.0)
    pl.plot( np.arange(-N_fermi_p,N_fermi_p), zarr)
#    pl.xlim(-20,20)
#    pl.ylim(-0.005,0.005)
    use_pl.set_title( string , fontsize=10)
    return

def plotUpUpPRePH( use_pl ):
    title = r"$\operatorname{Re}P_{2,\uparrow\uparrow}(\Omega_{PH},\omega_n)$"
    zarr = rep_upup_ph[:,shift+N_bose_p]
    plotP( use_pl, zarr, title )
    return

def plotUpDoPRePH( use_pl ):
    title = r"$\operatorname{Re}P_{2,\uparrow\downarrow}(\Omega_{PH},\omega_n)$"
    zarr = rep_updo_ph[:,shift+N_bose_p]
    plotP( use_pl, zarr, title )
    return

def plotUpDoPRePP( use_pl ):
    title = r"$\operatorname{Re}P_{2,\uparrow\downarrow}(\Omega_{PP},\omega_n)$"
    zarr = rep_updo_pp[:,shift+N_bose_p]
    plotP( use_pl, zarr, title )
    return

def plotUpUpPReXPH( use_pl ):
    title = r"$\operatorname{Re}P_{2,\uparrow\uparrow}(\Omega_{XPH},\omega_n)$"
    zarr = rep_upup_xph[:,shift+N_bose_p]
    plotP( use_pl, zarr, title )
    return

def plotUpDoPReXPH( use_pl ):
    title = r"$\operatorname{Re}P_{2,\uparrow\downarrow}(\Omega_{XPH},\omega_n)$"
    zarr = rep_updo_xph[:,shift+N_bose_p]
    plotP( use_pl, zarr, title )
    return

pl.suptitle(r"$U=$" + str('{0:.3f}'.format(float(U))) + r"     $\beta=$" + str('{0:.3f}'.format(float(beta))) + r"     $\Omega_{\rm PH}= \Omega_{\rm PP}=$" + str(shift) + r"$*2\pi/\beta$")

#--- Plot Physical
plotUpUpPRePH( pl.subplot(2,3,2) )
pl.xlabel(r"$\omega_n$")
plotUpDoPRePH( pl.subplot(2,3,5) )
pl.xlabel(r"$\omega_m$")
plotUpDoPRePP( pl.subplot(2,3,4) )
pl.xlabel(r"$\omega_n$")
plotUpUpPReXPH( pl.subplot(2,3,3) )
pl.xlabel(r"$\omega_n$")
plotUpDoPReXPH( pl.subplot(2,3,6) )
pl.xlabel(r"$\omega_n$")
pl.tight_layout()

pl.show()


#--------------------------------------------------------------KARRASCH PLOTTING--------------------------------------------------------


print("Plotting K ...")


#--- Read
rek_upup_ph =  np.array(f["/K_func/PH/RE_K_UPUP"])
imk_upup_ph = np.array(f["/K_func/PH/IM_K_UPUP"])
rek_updo_ph =  np.array(f["/K_func/PH/RE_K_UPDO"])
imk_updo_ph = np.array(f["/K_func/PH/IM_K_UPDO"])
rek_updo_pp =  np.array(f["/K_func/PP/RE_K_UPDO"])
imk_updo_pp = np.array(f["/K_func/PP/IM_K_UPDO"])
rek_updo_xph =  np.array(f["/K_func/XPH/RE_K_UPDO"])
imk_updo_xph = np.array(f["/K_func/XPH/IM_K_UPDO"])
rek_upup_xph =  np.array(f["/K_func/XPH/RE_K_UPUP"])
imk_upup_xph = np.array(f["/K_func/XPH/IM_K_UPUP"])

bdim = rek_upup_ph.shape[0]

print bdim

N_bose_k = (rek_upup_ph.shape[0]-1)/2

print N_bose_k


def plotK( use_pl, zarr, string ):
#    use_pl.set_aspect(1.0)
    pl.plot( np.arange(-N_bose_k,N_bose_k+1), zarr)
#    pl.xlim(-20,20)
    use_pl.set_title( string , fontsize=10)
    return

def plotUpUpKRePH( use_pl ):
    title = r"$\operatorname{Re}K_{2,\uparrow\uparrow}(\Omega_{PH},\omega_n)$"
    zarr = rek_upup_ph[:]
    plotK( use_pl, zarr, title )
    return

def plotUpDoKRePH( use_pl ):
    title = r"$\operatorname{Re}K_{2,\uparrow\downarrow}(\Omega_{PH},\omega_n)$"
    zarr = rek_updo_ph[:]
    plotK( use_pl, zarr, title )
    return

def plotUpDoKRePP( use_pl ):
    title = r"$\operatorname{Re}K_{2,\uparrow\downarrow}(\Omega_{PP},\omega_n)$"
    zarr = rek_updo_pp[:]
    plotK( use_pl, zarr, title )
    return

def plotUpUpKReXPH( use_pl ):
    title = r"$\operatorname{Re}K_{2,\uparrow\uparrow}(\Omega_{XPH},\omega_n)$"
    zarr = rek_upup_xph[:]
    plotK( use_pl, zarr, title )
    return

def plotUpDoKReXPH( use_pl ):
    title = r"$\operatorname{Re}K_{2,\uparrow\downarrow}(\Omega_{XPH},\omega_n)$"
    zarr = rek_updo_xph[:]
    plotK( use_pl, zarr, title )
    return

pl.suptitle(r"$U=$" + str('{0:.3f}'.format(float(U))) + r"     $\beta=$" + str('{0:.3f}'.format(float(beta))) + r"     $\Omega_{\rm PH}= \Omega_{\rm PP}=$" + str(shift) + r"$*2\pi/\beta$")

#--- Plot Physical
plotUpUpKRePH( pl.subplot(2,3,2) )
pl.xlabel(r"$\omega_n$")
plotUpDoKRePH( pl.subplot(2,3,5) )
pl.xlabel(r"$\omega_m$")
plotUpDoKRePP( pl.subplot(2,3,4) )
pl.xlabel(r"$\omega_n$")
plotUpUpKReXPH( pl.subplot(2,3,3) )
pl.xlabel(r"$\omega_n$")
plotUpDoKReXPH( pl.subplot(2,3,6) )
pl.xlabel(r"$\omega_n$")
pl.tight_layout()

pl.show()

