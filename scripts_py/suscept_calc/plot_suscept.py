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
#from _funtools import partial


pi = math.pi
U= 4.0
beta=20.0
#Read data file

#Data fRG 2loop K1eff
if ('../dat/'):
#2loop
    f = h5py.File('../dat/fRG/karrasch_approx/2loop/dat_U4_Beta20_PFCB4096_EDB_OMFL_KAT_2LOOP_SU2_K1CenterMinus.h5', 'r')   # Read (and write) the hdf5 file in the directory "dat" if existing
#1loop    
#    f = h5py.File('../dat/fRG/karrasch_approx/2loop/dat_U1_Beta20_PFCB4096_EDB_OMFL_KAT_2LOOP_SU2_K1CenterMinus.h5', 'r')   # Read (and write) the hdf5 file in the directory "dat" if existing
else:
    sys.exit("No data file")

#Data fRG 2loop full treatment for K1
if ('../dat/'):
#2loop
    g = h5py.File('../dat/fRG/full/2loop/dat_U4_Beta20_PFCB4096_EDB_OMFL_KAT_2LOOP_SU2.h5', 'r')   # Read (and write) the hdf5 file in the directory "dat" if existing
#1loop
#    g = h5py.File('../dat/fRG/full/2loop/dat_U1_Beta20_PFCB4096_EDB_OMFL_KAT_2LOOP_SU2.h5', 'r')   # Read (and write) the hdf5 file in the directory "dat" if existing
else:
    sys.exit("No data file")

#Data ED full treatment
if ('../dat/'):
    h = h5py.File('../dat/ED/dat_U4_Beta20_PFCB94_EDB_PARQ_SU2.h5', 'r')   # Read (and write) the hdf5 file in the directory "dat" if existing
else:
    sys.exit("No data file")

#Data PA full treatment
if ('../dat/'):
    l = h5py.File('../dat/PARQ_APPR/dat_U4_Beta20_PFCB96_EDB_PARQ_APPR_SU2.h5', 'r')   # Read (and write) the hdf5 file in the directory "dat" if existing
else:
    sys.exit("No data file")

#Suscept fRG 2loop k1eff

rechi_updo_ph =  np.array(f["/suscept_fun/RE_PH"][:])
imchi_updo_ph = np.array(f["/suscept_fun/IM_PH"][:])
rechi_updo_pp =  np.array(f["/suscept_fun/RE_PP"][:])
imchi_updo_pp = np.array(f["/suscept_fun/IM_PP"][:])
rechi_updo_xph =  np.array(f["/suscept_fun/RE_XPH"][:])
imchi_updo_xph = np.array(f["/suscept_fun/IM_XPH"][:])

rechi_upup_ph =  np.array(f["/suscept_fun/RE_PH"][:]-f["/suscept_fun/RE_XPH"][:])
imchi_upup_ph =  np.array(f["/suscept_fun/IM_PH"][:]-f["/suscept_fun/IM_XPH"][:])
rechi_upup_xph = -rechi_upup_ph
imchi_upup_xph = -imchi_upup_ph

#K1 fRG 2loop k1eff

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


#K1/Suscept fRG 2loop full treatment

#K1
rek_updo_ph_full =  np.array(g["/chi_func/RE_PH"][:,0,0,0,0,0])
imk_updo_ph_full = np.array(g["/chi_func/IM_PH"][:,0,0,0,0,0])
rek_updo_pp_full =  np.array(g["/chi_func/RE_PP"][:,0,0,0,0,0])
imk_updo_pp_full = np.array(g["/chi_func/IM_PP"][:,0,0,0,0,0])
rek_updo_xph_full =  np.array(g["/chi_func/RE_XPH"][:,0,0,0,0,0])
imk_updo_xph_full = np.array(g["/chi_func/IM_XPH"][:,0,0,0,0,0])

rek_upup_ph_full =  np.array(g["/chi_func/RE_PH"][:,0,0,0,0,0]-g["/chi_func/RE_XPH"][:,0,0,0,0,0])
imk_upup_ph_full =  np.array(g["/chi_func/IM_PH"][:,0,0,0,0,0]-g["/chi_func/IM_XPH"][:,0,0,0,0,0])
rek_upup_xph_full = -rek_upup_ph_full
imk_upup_xph_full = -imk_upup_ph_full

#Suscept
rechi_updo_ph_full =  np.array(g["/suscept_fun/RE_PH"][:])
imchi_updo_ph_full = np.array(g["/suscept_fun/IM_PH"][:])
rechi_updo_pp_full =  np.array(g["/suscept_fun/RE_PP"][:])
imchi_updo_pp_full = np.array(g["/suscept_fun/IM_PP"][:])
rechi_updo_xph_full =  np.array(g["/suscept_fun/RE_XPH"][:])
imchi_updo_xph_full = np.array(g["/suscept_fun/IM_XPH"][:])

rechi_upup_ph_full =  np.array(g["/suscept_fun/RE_PH"][:]-g["/suscept_fun/RE_XPH"][:])
imchi_upup_ph_full =  np.array(g["/suscept_fun/IM_PH"][:]-g["/suscept_fun/IM_XPH"][:])
rechi_upup_xph_full = -rechi_upup_ph_full
imchi_upup_xph_full = -imchi_upup_ph_full

#K1 ED
rek_updo_ph_ED =  np.array(h["/chi_func/RE_PH"][:,0,0,0,0,0])
imk_updo_ph_ED = np.array(h["/chi_func/IM_PH"][:,0,0,0,0,0])
rek_updo_pp_ED =  np.array(h["/chi_func/RE_PP"][:,0,0,0,0,0])
imk_updo_pp_ED = np.array(h["/chi_func/IM_PP"][:,0,0,0,0,0])
rek_updo_xph_ED =  np.array(h["/chi_func/RE_XPH"][:,0,0,0,0,0])
imk_updo_xph_ED = np.array(h["/chi_func/IM_XPH"][:,0,0,0,0,0])

rek_upup_ph_ED =  np.array(h["/chi_func/RE_PH"][:,0,0,0,0,0]-h["/chi_func/RE_XPH"][:,0,0,0,0,0])
imk_upup_ph_ED =  np.array(h["/chi_func/IM_PH"][:,0,0,0,0,0]-h["/chi_func/IM_XPH"][:,0,0,0,0,0])
rek_upup_xph_ED = -rek_upup_ph_ED
imk_upup_xph_ED = -imk_upup_ph_ED

#K1 PA
rek_updo_ph_PA =  np.array(l["/chi_func/RE_PH"][:,0,0,0,0,0])
imk_updo_ph_PA = np.array(l["/chi_func/IM_PH"][:,0,0,0,0,0])
rek_updo_pp_PA =  np.array(l["/chi_func/RE_PP"][:,0,0,0,0,0])
imk_updo_pp_PA = np.array(l["/chi_func/IM_PP"][:,0,0,0,0,0])
rek_updo_xph_PA =  np.array(l["/chi_func/RE_XPH"][:,0,0,0,0,0])
imk_updo_xph_PA = np.array(l["/chi_func/IM_XPH"][:,0,0,0,0,0])

rek_upup_ph_PA =  np.array(l["/chi_func/RE_PH"][:,0,0,0,0,0]-l["/chi_func/RE_XPH"][:,0,0,0,0,0])
imk_upup_ph_PA =  np.array(l["/chi_func/IM_PH"][:,0,0,0,0,0]-l["/chi_func/IM_XPH"][:,0,0,0,0,0])
rek_upup_xph_PA = -rek_upup_ph_PA
imk_upup_xph_PA = -imk_upup_ph_PA


#Grids

N_bose_suscept = (rechi_updo_ph.shape[0]-1)/2 #We assume all channels to have the same bosonic range

N_bose_k = (rek_updo_ph.shape[0]-1)/2 #We assume all channels to have the same bosonic range

N_bose_suscept_full = (rechi_updo_ph_full.shape[0]-1)/2 #We assume all channels to have the same bosonic range
N_bose_k_full = (rek_updo_ph_full.shape[0]-1)/2 #We assume all channels to have the same bosonic range

N_bose_k_ED = (rek_updo_ph_ED.shape[0]-1)/2 #We assume all channels to have the same bosonic range
N_bose_k_PA = (rek_updo_ph_PA.shape[0]-1)/2 #We assume all channels to have the same bosonic range


#-----Plot Settings
pl.figure(figsize=(12,4))

N_bose_k_plot = min(N_bose_suscept, N_bose_k, N_bose_suscept_full, N_bose_k_full, N_bose_k_ED, N_bose_k_PA)

print N_bose_k_plot

pl.rc('xtick', labelsize=15) 
pl.rc('ytick', labelsize=15) 

Xpos = np.array([(2*n)*pi/beta for n in range(0,N_bose_k_plot+1)])
Xneg = np.array([(2*n)*pi/beta for n in range(-N_bose_k_plot,1)])
X = np.array([(2*n)*pi/beta for n in range(-N_bose_k_plot,N_bose_k_plot+1)])

#-----Plot functions

def plotUpDoKRePH( use_pl ):
    pl.plot( X, np.array([rechi_updo_ph[wb + N_bose_suscept].real for wb in range(-N_bose_k_plot,N_bose_k_plot+1)]), linestyle='-', color='g', ms=3, mew=0.2,linewidth=2, label= 'chi_k1eff')
   
#    pl.plot( X, np.array([1./U*1./U*rek_updo_ph[wb + N_bose_k].real for wb in range(-N_bose_k_plot,N_bose_k_plot+1)]), linestyle='-', color='r', ms=3, mew=0.2,linewidth=2, label= 'k1eff_2l')
#    pl.plot( X, np.array([rechi_updo_ph_1l[wb + N_bose_suscept].real for wb in range(-N_bose_k_plot,N_bose_k_plot+1)]), linestyle='-', color='r', ms=3, mew=0.2,linewidth=2, label= 'chi_k1eff_1l')
    
    pl.plot( X, np.array([rechi_updo_ph_full[wb + N_bose_suscept_full].real for wb in range(-N_bose_k_plot,N_bose_k_plot+1)]), linestyle='-', color='b', ms=3, mew=0.2,linewidth=2,label= 'chi_k1full')
    
    pl.plot( X, np.array([rek_updo_ph_full[wb + N_bose_k_full].real/U/U for wb in range(-N_bose_k_plot,N_bose_k_plot+1)]), linestyle='--', color='r', ms=3, mew=0.2,linewidth=2, label= 'k1_full')
   
    pl.plot( X, np.array([rek_updo_ph_ED[wb + N_bose_k_ED].real/U/U for wb in range(-N_bose_k_plot,N_bose_k_plot+1)]), linestyle='--', color='black', ms=2, mew=0.2,linewidth=2, label= 'k1_ED')
    
   
    pl.plot( X, np.array([rek_updo_ph_PA[wb + N_bose_k_PA].real/U/U for wb in range(-N_bose_k_plot,N_bose_k_plot+1)]), linestyle='-', color='violet', ms=2, mew=0.2,linewidth=2, label= 'k1_PA')
    
    return

def plotUpDoKRePP( use_pl ):
    pl.plot( X, np.array([rechi_updo_pp[wb + N_bose_suscept].real for wb in range(-N_bose_k_plot,N_bose_k_plot+1)]), linestyle='-', color='g', ms=3, mew=0.2,linewidth=2)
    
#    pl.plot( X, np.array([1./U*1./U*rek_updo_pp[wb + N_bose_k].real for wb in range(-N_bose_k_plot,N_bose_k_plot+1)]), linestyle='-', color='r', ms=3, mew=0.2,linewidth=2)
#    pl.plot( X, np.array([rechi_updo_pp_1l[wb + N_bose_suscept].real for wb in range(-N_bose_k_plot,N_bose_k_plot+1)]), linestyle='-', color='r', ms=3, mew=0.2,linewidth=2)
    
    pl.plot( X, np.array([rechi_updo_pp_full[wb + N_bose_suscept_full].real for wb in range(-N_bose_k_plot,N_bose_k_plot+1)]), linestyle='-', color='b', ms=3, mew=0.2,linewidth=2)
    
    pl.plot( X, np.array([1/U*1/U*rek_updo_pp_full[wb + N_bose_k_full].real for wb in range(-N_bose_k_plot,N_bose_k_plot+1)]), linestyle='--', color='r', ms=3, mew=0.2,linewidth=2)

    pl.plot( X, np.array([1./U*1./U*rek_updo_pp_ED[wb + N_bose_k_ED].real for wb in range(-N_bose_k_plot,N_bose_k_plot+1)]), linestyle='--', color='black', ms=2, mew=0.2,linewidth=2)
    
    pl.plot( X, np.array([rek_updo_pp_PA[wb + N_bose_k_PA].real/U/U for wb in range(-N_bose_k_plot,N_bose_k_plot+1)]), linestyle='-', color='violet', ms=2, mew=0.2,linewidth=2)

    return

def plotUpDoKReXPH( use_pl ):
    pl.plot( X, np.array([rechi_updo_xph[wb + N_bose_suscept].real for wb in range(-N_bose_k_plot,N_bose_k_plot+1)]), linestyle='-', color='g', ms=3, mew=0.2,linewidth=2)
    
#    pl.plot( X, np.array([1./U*1./U*rek_updo_xph[wb + N_bose_k].real for wb in range(-N_bose_k_plot,N_bose_k_plot+1)]), linestyle='-', color='r', ms=3, mew=0.2,linewidth=2)
#    pl.plot( X, np.array([rechi_updo_xph_1l[wb + N_bose_suscept].real for wb in range(-N_bose_k_plot,N_bose_k_plot+1)]), linestyle='-', color='r', ms=3, mew=0.2,linewidth=2)
    
    pl.plot( X, np.array([rechi_updo_xph_full[wb + N_bose_suscept_full].real for wb in range(-N_bose_k_plot,N_bose_k_plot+1)]), linestyle='-', color='b', ms=3, mew=0.2,linewidth=2)
    
    pl.plot( X, np.array([1./U*1./U*rek_updo_xph_full[wb + N_bose_k_full].real for wb in range(-N_bose_k_plot,N_bose_k_plot+1)]), linestyle='--', color='r', ms=3, mew=0.2,linewidth=2)
   
    pl.plot( X, np.array([1./U*1./U*rek_updo_xph_ED[wb + N_bose_k_ED].real for wb in range(-N_bose_k_plot,N_bose_k_plot+1)]), linestyle='--', color='black', ms=2, mew=0.2,linewidth=2)
    
    
    pl.plot( X, np.array([rek_updo_xph_PA[wb + N_bose_k_PA].real/U/U for wb in range(-N_bose_k_plot,N_bose_k_plot+1)]), linestyle='-', color='violet', ms=2, mew=0.2,linewidth=2)

    return

#--- Plot Physical
plotUpDoKRePH( pl.subplot(1,3,2) )
ax1 = pl.subplot(1,3,2)
pl.text(0.5,1.04,r"${\chi}_{ph,\uparrow\downarrow}^{\Omega}$",horizontalalignment='center', fontsize=20,transform = ax1.transAxes)
pl.xlabel(r"$\Omega$",fontsize=20)
pl.legend(loc=4,prop={'size':10})
pl.xlim([-3,3])
#pl.ylim([-20,0])
plotUpDoKRePP( pl.subplot(1,3,1) )
ax2=pl.subplot(1,3,1)
pl.text(0.5, 1.04,r"${\chi}_{pp,\uparrow\downarrow}^{\Omega}$",horizontalalignment='center', fontsize=20, transform = ax2.transAxes)
pl.xlabel(r"$\Omega$",fontsize=20)
pl.xlim([-3,3])
#pl.ylim([-0.0,2.0])
#pl.yticks([0.0, 0.5, 1.0, 1.5, 2.0])
pl.legend(loc=1)
plotUpDoKReXPH( pl.subplot(1,3,3) )
ax3=pl.subplot(1,3,3) 
pl.text(0.5, 1.04,r"${\chi}_{\overline{ph},\uparrow\downarrow}^{\Omega}$",horizontalalignment='center', fontsize=20, transform = ax3.transAxes)
pl.xlim([-3,3])
pl.xlabel(r"$\Omega$",fontsize=20)
pl.legend(loc=4)
#pl.yticks([-40, -30, -20, -10, 0.0]) 
pl.subplots_adjust(left=0.07, bottom=0.15, right=0.95, top=0.88, wspace=0.2, hspace=0.2)
#pl.tight_layout()

#--- Save to file
pl.savefig("Suscept_U1_2loop_comparison.png", dpi = 150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()

#-----Plot Settings_paper
#pl.figure(figsize=(8,8))

N_bose_k_plot = min(N_bose_suscept, N_bose_k, N_bose_suscept_full, N_bose_k_full, N_bose_k_ED, N_bose_k_PA)

print N_bose_k_plot

pl.rc('xtick', labelsize=20) 
pl.rc('ytick', labelsize=20) 

Xpos = np.array([(2*n)*pi/beta for n in range(0,N_bose_k_plot+1)])
Xneg = np.array([(2*n)*pi/beta for n in range(-N_bose_k_plot,1)])
X = np.array([(2*n)*pi/beta for n in range(-N_bose_k_plot,N_bose_k_plot+1)])

#-----Plot functions

def plotUpDoKRePH( use_pl ):
    pl.plot( X, np.array([-rechi_updo_ph[wb + N_bose_suscept].real for wb in range(-N_bose_k_plot,N_bose_k_plot+1)]), linestyle='-', color='limegreen', ms=10, mew=2,linewidth=2, label= '$\chi_{\mathrm{VC}, eff}$')
   
#    pl.plot( X, np.array([1./U*1./U*rek_updo_ph[wb + N_bose_k].real for wb in range(-N_bose_k_plot,N_bose_k_plot+1)]), linestyle='-', color='r', ms=3, mew=0.2,linewidth=2, label= 'k1eff_2l')
#    pl.plot( X, np.array([rechi_updo_ph_1l[wb + N_bose_suscept].real for wb in range(-N_bose_k_plot,N_bose_k_plot+1)]), linestyle='-', color='r', ms=3, mew=0.2,linewidth=2, label= 'chi_k1eff_1l')
    
    pl.plot( X, np.array([-rechi_updo_ph_full[wb + N_bose_suscept_full].real for wb in range(-N_bose_k_plot,N_bose_k_plot+1)]), linestyle='-', color='darkorchid', ms=10, mew=2,linewidth=2,label= '$\chi_{\mathrm{VC}}$')
    
    pl.plot( X, np.array([-rek_updo_ph_full[wb + N_bose_k_full].real/U/U for wb in range(-N_bose_k_plot,N_bose_k_plot+1)]), linestyle='-', color='g', ms=10, mew=2,linewidth=2, label= '$\chi$')
   
    pl.plot( X, np.array([-rek_updo_ph_ED[wb + N_bose_k_ED].real/U/U for wb in range(-N_bose_k_plot,N_bose_k_plot+1)]), linestyle='--', color='black', ms=10, mew=2,linewidth=2, label= '$\chi_{\mathrm{ED}}$')
    
   
    pl.plot( X, np.array([-rek_updo_ph_PA[wb + N_bose_k_PA].real/U/U for wb in range(-N_bose_k_plot,N_bose_k_plot+1)]), linestyle='-', color='darkturquoise', ms=10, mew=2,linewidth=2, label= '$\chi_{\mathrm{PA}}$')
    
    return

def plotUpDoKRePP( use_pl ):
    pl.plot( X, np.array([rechi_updo_pp[wb + N_bose_suscept].real for wb in range(-N_bose_k_plot,N_bose_k_plot+1)]), linestyle='-', color='g', ms=3, mew=0.2,linewidth=2)
    
#    pl.plot( X, np.array([1./U*1./U*rek_updo_pp[wb + N_bose_k].real for wb in range(-N_bose_k_plot,N_bose_k_plot+1)]), linestyle='-', color='r', ms=3, mew=0.2,linewidth=2)
#    pl.plot( X, np.array([rechi_updo_pp_1l[wb + N_bose_suscept].real for wb in range(-N_bose_k_plot,N_bose_k_plot+1)]), linestyle='-', color='r', ms=3, mew=0.2,linewidth=2)
    
    pl.plot( X, np.array([rechi_updo_pp_full[wb + N_bose_suscept_full].real for wb in range(-N_bose_k_plot,N_bose_k_plot+1)]), linestyle='-', color='b', ms=3, mew=0.2,linewidth=2)
    
    pl.plot( X, np.array([1/U*1/U*rek_updo_pp_full[wb + N_bose_k_full].real for wb in range(-N_bose_k_plot,N_bose_k_plot+1)]), linestyle='--', color='r', ms=3, mew=0.2,linewidth=2)

    pl.plot( X, np.array([1./U*1./U*rek_updo_pp_ED[wb + N_bose_k_ED].real for wb in range(-N_bose_k_plot,N_bose_k_plot+1)]), linestyle='--', color='black', ms=2, mew=0.2,linewidth=2)
    
    pl.plot( X, np.array([rek_updo_pp_PA[wb + N_bose_k_PA].real/U/U for wb in range(-N_bose_k_plot,N_bose_k_plot+1)]), linestyle='-', color='violet', ms=2, mew=0.2,linewidth=2)

    return

def plotUpDoKReXPH( use_pl ):
    pl.plot( X, np.array([rechi_updo_xph[wb + N_bose_suscept].real for wb in range(-N_bose_k_plot,N_bose_k_plot+1)]), linestyle='-', color='g', ms=3, mew=0.2,linewidth=2)
    
#    pl.plot( X, np.array([1./U*1./U*rek_updo_xph[wb + N_bose_k].real for wb in range(-N_bose_k_plot,N_bose_k_plot+1)]), linestyle='-', color='r', ms=3, mew=0.2,linewidth=2)
#    pl.plot( X, np.array([rechi_updo_xph_1l[wb + N_bose_suscept].real for wb in range(-N_bose_k_plot,N_bose_k_plot+1)]), linestyle='-', color='r', ms=3, mew=0.2,linewidth=2)
    
    pl.plot( X, np.array([rechi_updo_xph_full[wb + N_bose_suscept_full].real for wb in range(-N_bose_k_plot,N_bose_k_plot+1)]), linestyle='-', color='b', ms=3, mew=0.2,linewidth=2)
    
    pl.plot( X, np.array([1./U*1./U*rek_updo_xph_full[wb + N_bose_k_full].real for wb in range(-N_bose_k_plot,N_bose_k_plot+1)]), linestyle='--', color='r', ms=3, mew=0.2,linewidth=2)
   
    pl.plot( X, np.array([1./U*1./U*rek_updo_xph_ED[wb + N_bose_k_ED].real for wb in range(-N_bose_k_plot,N_bose_k_plot+1)]), linestyle='--', color='black', ms=2, mew=0.2,linewidth=2)
    
    
    pl.plot( X, np.array([rek_updo_xph_PA[wb + N_bose_k_PA].real/U/U for wb in range(-N_bose_k_plot,N_bose_k_plot+1)]), linestyle='-', color='violet', ms=2, mew=0.2,linewidth=2)

    return

#--- Plot Physical
plotUpDoKRePH( pl.subplot(1,1,1) )
#ax1 = pl.subplot(1,1,1)
#pl.text(0.5,1.04,r"${\chi}_{ph,\uparrow\downarrow}^{\Omega}$",horizontalalignment='center', fontsize=30,transform = ax1.transAxes)
pl.title(r"${\chi}_{ph,\uparrow\downarrow}^{\Omega}$",fontsize=30, y=1.04)
pl.xlabel(r"$\Omega$",fontsize=30)
pl.legend(loc=1,prop={'size':20})
pl.xlim([-3,3])
pl.ylim([0.0,1.5])
pl.yticks([0.0,0.5,1.0,1.5])
#pl.ylim([-20,0])
#plotUpDoKRePP( pl.subplot(1,3,1) )
#ax2=pl.subplot(1,3,1)
#pl.text(0.5, 1.04,r"${\chi}_{pp,\uparrow\downarrow}^{\Omega}$",horizontalalignment='center', fontsize=20, transform = ax2.transAxes)
#pl.xlabel(r"$\Omega$",fontsize=20)
#pl.xlim([-3,3])
##pl.ylim([-0.0,2.0])
##pl.yticks([0.0, 0.5, 1.0, 1.5, 2.0])
#pl.legend(loc=1)
#plotUpDoKReXPH( pl.subplot(1,3,3) )
#ax3=pl.subplot(1,3,3) 
#pl.text(0.5, 1.04,r"${\chi}_{\overline{ph},\uparrow\downarrow}^{\Omega}$",horizontalalignment='center', fontsize=20, transform = ax3.transAxes)
#pl.xlim([-3,3])
#pl.xlabel(r"$\Omega$",fontsize=20)
#pl.legend(loc=4)
#pl.yticks([-40, -30, -20, -10, 0.0]) 
#pl.subplots_adjust(left=0.07, bottom=0.15, right=0.95, top=0.88, wspace=0.2, hspace=0.2)
pl.tight_layout()

#--- Save to file
pl.savefig("Suscept_U4_2loop_comparison_PH.png", dpi = 150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()

#Plotting Notes on the susceptibility flow

#-----Plot Settings_paper
pl.figure(figsize=(8,8))

N_bose_k_plot = min(N_bose_suscept, N_bose_k, N_bose_suscept_full, N_bose_k_full, N_bose_k_ED, N_bose_k_PA)

print N_bose_k_plot

pl.rc('xtick', labelsize=18) 
pl.rc('ytick', labelsize=18) 

Xpos = np.array([(2*n)*pi/beta for n in range(0,N_bose_k_plot+1)])
Xneg = np.array([(2*n)*pi/beta for n in range(-N_bose_k_plot,1)])
X = np.array([(2*n)*pi/beta for n in range(-N_bose_k_plot,N_bose_k_plot+1)])

#-----Plot functions

def plotUpDoKRePH( use_pl ):
#    pl.plot( X, np.array([rechi_updo_ph[wb + N_bose_suscept].real for wb in range(-N_bose_k_plot,N_bose_k_plot+1)]), linestyle='-', color='g', ms=3, mew=0.2,linewidth=2, label= '$\chi_{\mathrm{VC}, eff}$')
   
#    pl.plot( X, np.array([1./U*1./U*rek_updo_ph[wb + N_bose_k].real for wb in range(-N_bose_k_plot,N_bose_k_plot+1)]), linestyle='-', color='r', ms=3, mew=0.2,linewidth=2, label= 'k1eff_2l')
#    pl.plot( X, np.array([rechi_updo_ph_1l[wb + N_bose_suscept].real for wb in range(-N_bose_k_plot,N_bose_k_plot+1)]), linestyle='-', color='r', ms=3, mew=0.2,linewidth=2, label= 'chi_k1eff_1l')
    
    pl.plot( X, np.array([rechi_updo_ph_full[wb + N_bose_suscept_full].real for wb in range(-N_bose_k_plot,N_bose_k_plot+1)]), linestyle='-', color='b', ms=3, mew=0.2,linewidth=2,label= '$\chi^{\Lambda_{f}}_{\mathrm{VC}}$')
    
    pl.plot( X, np.array([rek_updo_ph_full[wb + N_bose_k_full].real/U/U for wb in range(-N_bose_k_plot,N_bose_k_plot+1)]), linestyle='-', color='r', ms=3, mew=0.2,linewidth=2, label= '$\chi^{\Lambda_f}$')
   
    pl.plot( X, np.array([rek_updo_ph_ED[wb + N_bose_k_ED].real/U/U for wb in range(-N_bose_k_plot,N_bose_k_plot+1)]), linestyle='--', color='black', ms=3, mew=0.2,linewidth=2, label= '$\chi_{\mathrm{ED}}$')
    
   
#    pl.plot( X, np.array([rek_updo_ph_PA[wb + N_bose_k_PA].real/U/U for wb in range(-N_bose_k_plot,N_bose_k_plot+1)]), linestyle='-', color='violet', ms=2, mew=0.2,linewidth=2, label= '$\chi_{\mathrm{PA}}$')
    
    return


#--- Plot Physical
plotUpDoKRePH( pl.subplot(1,1,1) )
ax1 = pl.subplot(1,1,1)
pl.text(0.5,1.04,r"${\chi}_{ph,\uparrow\downarrow}^{\Omega}$",horizontalalignment='center', fontsize=24,transform = ax1.transAxes)
pl.xlabel(r"$\Omega$",fontsize=22)
pl.legend(loc=4,prop={'size':20})
pl.xlim([-3,3])
#pl.ylim([-20,0])
#plotUpDoKRePP( pl.subplot(1,3,1) )
#ax2=pl.subplot(1,3,1)
#pl.text(0.5, 1.04,r"${\chi}_{pp,\uparrow\downarrow}^{\Omega}$",horizontalalignment='center', fontsize=20, transform = ax2.transAxes)
#pl.xlabel(r"$\Omega$",fontsize=20)
#pl.xlim([-3,3])
##pl.ylim([-0.0,2.0])
##pl.yticks([0.0, 0.5, 1.0, 1.5, 2.0])
#pl.legend(loc=1)
#plotUpDoKReXPH( pl.subplot(1,3,3) )
#ax3=pl.subplot(1,3,3) 
#pl.text(0.5, 1.04,r"${\chi}_{\overline{ph},\uparrow\downarrow}^{\Omega}$",horizontalalignment='center', fontsize=20, transform = ax3.transAxes)
#pl.xlim([-3,3])
#pl.xlabel(r"$\Omega$",fontsize=20)
#pl.legend(loc=4)
#pl.yticks([-40, -30, -20, -10, 0.0]) 
#pl.subplots_adjust(left=0.07, bottom=0.15, right=0.95, top=0.88, wspace=0.2, hspace=0.2)
#pl.tight_layout()

#--- Save to file
pl.savefig("Suscept_U4_2loop_comparison_PH_notes.png", dpi = 150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()
