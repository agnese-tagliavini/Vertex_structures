#!/usr/bin/python

#===========================================================================================================================================
#
#                   NOTE: Plotting script to show the results of inverting the Bethe-Salpeter equations
#                         by using Agnese's corrections -> METHOD 2 in the main program
#
#                   CONTENTS:
#                   1) 2P irreducible vertex (Gamma) scanned nu'=0 and Omega=0
#                   2) 2P irreducible vertex (Gamma) scanned nu'= 0 and Omega \neq 0
#                   3) 2P irreducible vertex (Gamma) scanned nu'= nu and Omega=0
#                   4) 2P irreducible vertex (Gamma) scanned nu'= nu and Omega \neq 0
#                   5) 2P irreducible vertex differences for different sizes of the 00 block, no corrections (nu'= 0, Omega=0)
#   `               6) 2P irreducible vertex differences for different sizes of the 00 block and no corrections (nu'= 0, Omega \neq 0)
#                   7) 2P irreducible vertex differences for different sizes of the 00 block andno corrections (nu'= nu, Omega=0)
#                   8) 2P irreducible vertex differences for different sizes of the 00 block and no corrections (nu'= nu, Omega \neq 0)
#                   9) 2P irreducible vertex result convergence with respect to the inversion range
#
#===========================================================================================================================================
#--------------------------------------IMPORTS ------------------------------------------

import h5py
import matplotlib.pyplot as pl
import numpy as np
import os
import sys
import subprocess
import math
from agneselibrary.mymath import *

#--------------------------------------SETTINGS ------------------------------------------
#if -1, change overall vertex sign (switch definition)
vert_mul = -1
pi = math.pi

#--------------------------------------MANAGING FILES ------------------------------------------

def run(command):
    output = subprocess.check_output(command, shell=True)
    return output

most_recently_edited = run("ls -Art dat/ | tail -n 1")

#------COMPARISON METHOD 1 WITH AND WITHOUT CORRECTIONS -> CONVERGENCE WITH RESPECT TO THE INVERSION RANGE 

#fname2 = "../../dat/dat_U1_Beta50_PFCB37_PARQ_SU2_METH2_INVR1xCOUNT_ED_ONESHOT.h5" 
#fname2 = "../../dat/dat_U1_Beta50_PFCB150_PARQ_SU2_METH2_INVR1xCOUNT_CORR_SELFCON.h5"
#fname2 = "../../dat/dat_U1_Beta50_PFCB1000_PARQ_SU2_METH2_INVR1xCOUNT_ED.h5" 
#fname1 = "../../dat/dat_U1_Beta50_PFCB1000_PARQ_SU2_METH2_INVR1xCOUNT_CORR_ED.h5" 
fname2 = "../../dat/H5FILES/BETA50/4SITES/U1/POSTPROC/dat_U1_Beta50_PFCB1000_PARQ_SU2_METH2_INVR1xCOUNT_CORR_ED.h5" 
fname1 = "../../dat/dat_U1_Beta50_PFCB0_PARQ_SU2_METH2_INVR1xCOUNT_BIGINVR_CORR_ED_ONESHOT.h5" 


#------ DECIDE WHICH COMPARISON YOU WANT TO ANALYZE


if len(sys.argv) > 1:
    fname1 = str(sys.argv[1])
if len(sys.argv) > 1:
    fname2 = str(sys.argv[1])

fname1 = fname1.rstrip('\n') # strip newline of fname
f1 = h5py.File(fname1, "r")

fname2 = fname2.rstrip('\n') # strip newline of fname
f2 = h5py.File(fname2, "r")

os.system('mkdir -p log')
os.system('mkdir -p plots')


#--------------------------------------READ PARAMETERS FROM FILE ------------------------------------------

parVals = f1["/Params"].attrs.values()      #We are comparing different method with the same parameters

UINT =  parVals[0] # follows order in output.cpp
BETA =  parVals[1]
B =     parVals[2]
GAM_L = parVals[3]
DEL =   parVals[4]
EPS =   parVals[5]
PHI =   parVals[6]

shift=0

#--------------------------------------GENERAL PLOT SETTINGS------------------------------------------

pl.rc('xtick', labelsize=9) 
pl.rc('ytick', labelsize=9) 
pl.rc('text', usetex=True)
pl.rc('text.latex', preamble='\usepackage{amsmath}')

RE = r"$\operatorname{Re}"
IM = r"$\operatorname{Im}"


#-------------------------------------GAMMA PLOTTING ------------------------------------------

print("Plotting phi ...")

#--- Read gamma method 2 (INVR1xVERTR)
regamma_pp1 = vert_mul * np.array(f1["/phi_func/RE_PP"])
imgamma_pp1 = vert_mul * np.array(f1["/phi_func/IM_PP"])
regamma_ph1 = vert_mul * np.array(f1["/phi_func/RE_PH"])
imgamma_ph1 = vert_mul * np.array(f1["/phi_func/IM_PH"])
regamma_xph1 = vert_mul * np.array(f1["/phi_func/RE_XPH"])
imgamma_xph1 = vert_mul * np.array(f1["/phi_func/IM_XPH"])

#--- Read gamma method 2 (INVR2xVERTR)
regamma_pp2 = vert_mul * np.array(f2["/phi_func/RE_PP"])
imgamma_pp2 = vert_mul * np.array(f2["/phi_func/IM_PP"])
regamma_ph2 = vert_mul * np.array(f2["/phi_func/RE_PH"])
imgamma_ph2 = vert_mul * np.array(f2["/phi_func/IM_PH"])
regamma_xph2 = vert_mul * np.array(f2["/phi_func/RE_XPH"])
imgamma_xph2 = vert_mul * np.array(f2["/phi_func/IM_XPH"])


#--- RANGE FERMIONIC/BOSONIC FREQUENCIES

bdim = regamma_pp1.shape[0]
fdim = regamma_pp1.shape[1]
fdimo2 = fdim/2

#phigrid = np.array(f["/gamma_func/fgrid"])
gammagrid_plot = np.array([(2*n+1)*pi/BETA for n in range(-fdimo2,fdimo2)])

#---  Helper functions (include translation from  nambu to physical CAUTION : USING TIME REVERSAL SYMM) 
def plotgamma_axis( use_pl, arr1, arr2, string ):
    zarr1 = np.array([ arr1[shift + (bdim-1)/2,n,fdimo2,0,0,0,0,0,0,0] for n in range(fdim)])
    zarr2 = np.array([ arr2[shift + (bdim-1)/2,n,fdimo2,0,0,0,0,0,0,0] for n in range(fdim)])
    pl.plot( gammagrid_plot, zarr1, marker = 'o', linestyle='-', color='r', ms=3, mew=0.2, label=r"10x")
    pl.plot( gammagrid_plot, zarr2, marker = 'x',linestyle='-', color='b',ms=3, mew=0.2, label=r"20x")
    pl.xlim([0.3*min(gammagrid_plot),0.3*max(gammagrid_plot)])
    use_pl.set_title( string , fontsize=10 )
    pl.legend(prop={'size':5})

    return

def plotgamma_diag1( use_pl, arr1, arr2, string ):
    zarr1 = np.array([ arr1[shift + (bdim-1)/2,n,n,0,0,0,0,0,0,0] for n in range(fdim)])
    zarr2 = np.array([ arr2[shift + (bdim-1)/2,n,n,0,0,0,0,0,0,0] for n in range(fdim)])
    pl.plot( gammagrid_plot, zarr1, marker = 'o', linestyle='-', color='r', ms=3, mew=0.2, label=r"10x")
    pl.plot( gammagrid_plot, zarr2, marker = 'o',linestyle='-', color='b',ms=3, mew=0.2, label=r"20x")
    pl.xlim([0.3*min(gammagrid_plot),0.3*max(gammagrid_plot)])
    use_pl.set_title( string , fontsize=10 )
    pl.legend(prop={'size':5})

    return

#--- Plot along \nu'=0 

pl.suptitle(r"$U=$" + str(UINT) + r"     $\beta=$" + str(BETA) + r"     $\epsilon=$" + str(EPS) +  r"     $\Omega=$" + str(shift) + r"$*2\pi/\beta$")

plotgamma_axis( pl.subplot(2,2,1), regamma_pp1, regamma_pp2, RE + r"\Gamma^{\nu, \nu'=0, \omega =0}_{pp}$" ) # flip sign of w_out
pl.xlabel(r"$\nu$")
plotgamma_axis( pl.subplot(2,2,2), regamma_pp1 - regamma_pp1[:,:,::-1,:,:,:,:,:,:,:], regamma_pp2 - regamma_pp2[:,:,::-1,:,:,:,:,:,:,:], RE + r"\Gamma^{\nu, \nu'=0, \omega =0}_{t}$" ) # flip sign of w_out
pl.xlabel(r"$\nu$")
plotgamma_axis( pl.subplot(2,2,3), 2*regamma_ph1 - regamma_xph1, 2*regamma_ph2 - regamma_xph2, RE + r"\Gamma^{\nu, \nu'=0, \omega =0}_{d}$" )
pl.xlabel(r"$\nu$")
plotgamma_axis( pl.subplot(2,2,4), - regamma_xph1, - regamma_xph2, RE + r"\Gamma^{\nu, \nu'=0, \omega =0}_{m}$" )
pl.xlabel(r"$\nu$")

pl.tight_layout()

#--- Save to file
pl.savefig("plots/PHI_METHOD1_axis_nu'=0.png", dpi = 150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()

shift=20

pl.suptitle(r"$U=$" + str(UINT) + r"     $\beta=$" + str(BETA) + r"     $\epsilon=$" + str(EPS) +  r"     $\Omega=$" + str(shift) + r"$*2\pi/\beta$")

plotgamma_axis( pl.subplot(2,2,1), regamma_pp1, regamma_pp2, RE + r"\Gamma^{\nu, \nu'=0, \omega =0}_{pp}$" ) # flip sign of w_out
pl.xlabel(r"$\nu$")
plotgamma_axis( pl.subplot(2,2,2), regamma_pp1 - regamma_pp1[:,:,::-1,:,:,:,:,:,:,:], regamma_pp2 - regamma_pp2[:,:,::-1,:,:,:,:,:,:,:], RE + r"\Gamma^{\nu, \nu'=0, \omega =0}_{t}$" ) # flip sign of w_out
pl.xlabel(r"$\nu$")
plotgamma_axis( pl.subplot(2,2,3), 2*regamma_ph1 - regamma_xph1, 2*regamma_ph2 - regamma_xph2, RE + r"\Gamma^{\nu, \nu'=0, \omega =0}_{d}$" )
pl.xlabel(r"$\nu$")
plotgamma_axis( pl.subplot(2,2,4), - regamma_xph1, - regamma_xph2, RE + r"\Gamma^{\nu, \nu'=0, \omega =0}_{m}$" )
pl.xlabel(r"$\nu$")

pl.tight_layout()

#--- Save to file
pl.savefig("plots/PHI_METHOD1_Om="+str(shift)+"_axis_nu'=0.png", dpi = 150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()

#--- Plot along \nu'=\nu 

pl.suptitle(r"$U=$" + str(UINT) + r"     $\beta=$" + str(BETA) + r"     $\epsilon=$" + str(EPS) +  r"     $\Omega=$" + str(shift) + r"$*2\pi/\beta$")

plotgamma_diag1( pl.subplot(2,2,1), regamma_pp1, regamma_pp2, RE + r"\Gamma^{\nu, \nu'=0, \omega =0}_{pp}$" ) # flip sign of w_out
pl.xlabel(r"$\nu$")
plotgamma_diag1( pl.subplot(2,2,2), regamma_pp1 - regamma_pp1[:,:,::-1,:,:,:,:,:,:,:], regamma_pp2 - regamma_pp2[:,:,::-1,:,:,:,:,:,:,:], RE + r"\Gamma^{\nu, \nu'=0, \omega =0}_{t}$" ) # flip sign of w_out
pl.xlabel(r"$\nu$")
plotgamma_diag1( pl.subplot(2,2,3), 2*regamma_ph1 - regamma_xph1, 2*regamma_ph2 - regamma_xph2, RE + r"\Gamma^{\nu, \nu'=0, \omega =0}_{d}$" )
pl.xlabel(r"$\nu$")
plotgamma_diag1( pl.subplot(2,2,4), - regamma_xph1, - regamma_xph2, RE + r"\Gamma^{\nu, \nu'=0, \omega =0}_{m}$" )
pl.xlabel(r"$\nu$")

pl.tight_layout()

#--- Save to file
pl.savefig("plots/PHI_METHOD1_nu=nu'.png", dpi = 150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()

shift=20

pl.suptitle(r"$U=$" + str(UINT) + r"     $\beta=$" + str(BETA) + r"     $\epsilon=$" + str(EPS) +  r"     $\Omega=$" + str(shift) + r"$*2\pi/\beta$")

plotgamma_diag1( pl.subplot(2,2,1), regamma_pp1, regamma_pp2, RE + r"\Gamma^{\nu, \nu'=0, \omega =0}_{pp}$" ) # flip sign of w_out
pl.xlabel(r"$\nu$")
plotgamma_diag1( pl.subplot(2,2,2), regamma_pp1 - regamma_pp1[:,:,::-1,:,:,:,:,:,:,:], regamma_pp2 - regamma_pp2[:,:,::-1,:,:,:,:,:,:,:], RE + r"\Gamma^{\nu, \nu'=0, \omega =0}_{t}$" ) # flip sign of w_out
pl.xlabel(r"$\nu$")
plotgamma_diag1( pl.subplot(2,2,3), 2*regamma_ph1 - regamma_xph1, 2*regamma_ph2 - regamma_xph2, RE + r"\Gamma^{\nu, \nu'=0, \omega =0}_{d}$" )
pl.xlabel(r"$\nu$")
plotgamma_diag1( pl.subplot(2,2,4), - regamma_xph1, - regamma_xph2, RE + r"\Gamma^{\nu, \nu'=0, \omega =0}_{m}$" )
pl.xlabel(r"$\nu$")

pl.tight_layout()

#--- Save to file
pl.savefig("plots/PHI_METHOD1_Om="+str(shift)+"_nu=nu'.png", dpi = 150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()

#-------------------------------------GAMMA PLOTTING DIFFERENCES WITH RESPECT TO METHOD 1 INVR4xVERTR_ASYR10xINVR ------------------------------------------

print("Plotting delta phi ...")

#--- Read gamma method 2 (INVR1xVERTR)
diff_regamma_pp1 = vert_mul * np.subtract(f1["/phi_func/RE_PP"],f2["/phi_func/RE_PP"])
diff_imgamma_pp1 = vert_mul * np.subtract(f1["/phi_func/IM_PP"],f2["/phi_func/IM_PP"] )
diff_regamma_ph1 = vert_mul * np.subtract(f1["/phi_func/RE_PH"],f2["/phi_func/RE_PH"] )
diff_imgamma_ph1 = vert_mul * np.subtract(f1["/phi_func/IM_PH"],f2["/phi_func/IM_PH"] )
diff_regamma_xph1 = vert_mul * np.subtract(f1["/phi_func/RE_XPH"],f2["/phi_func/RE_XPH"])
diff_imgamma_xph1 = vert_mul * np.subtract(f1["/phi_func/IM_XPH"],f2["/phi_func/IM_XPH"])


#--- RANGE FERMIONIC/BOSONIC FREQUENCIES

bdim = regamma_pp1.shape[0]
fdim = regamma_pp1.shape[1]
fdimo2 = fdim/2

#phigrid = np.array(f["/gamma_func/fgrid"])
gammagrid_plot = np.array([(2*n+1)*pi/BETA for n in range(-fdimo2,fdimo2)])

#---  Helper functions (include translation from  nambu to physical CAUTION : USING TIME REVERSAL SYMM) 
def plotgamma_axis( use_pl, arr1,string ):
    zarr1 = np.array([ arr1[shift + (bdim-1)/2,n,fdimo2,0,0,0,0,0,0,0] for n in range(fdim)])
    pl.plot( gammagrid_plot, zarr1, marker = 'o', linestyle='-', color='r', ms=3, mew=0.2, label=r"10x-20x")
    pl.xlim([0.3*min(gammagrid_plot),0.3*max(gammagrid_plot)])
    use_pl.set_title( string , fontsize=10 )
    pl.legend(prop={'size':5})

    return

def plotgamma_diag1( use_pl, arr1, string ):
    zarr1 = np.array([ arr1[shift + (bdim-1)/2,n,n,0,0,0,0,0,0,0] for n in range(fdim)])
    pl.plot( gammagrid_plot, zarr1, marker = 'o', linestyle='-', color='r', ms=3, mew=0.2, label=r"10x-20x")
    pl.xlim([0.3*min(gammagrid_plot),0.3*max(gammagrid_plot)])
    use_pl.set_title( string , fontsize=10 )
    pl.legend(prop={'size':5})

    return

#--- Plot along axis \nu'=0

shift = 0
pl.suptitle(r"$U=$" + str(UINT) + r"     $\beta=$" + str(BETA) + r"     $\epsilon=$" + str(EPS) +  r"     $\Omega=$" + str(shift) + r"$*2\pi/\beta$" + r"     Notation: $\phi(\Omega,\omega,\omega')$")

plotgamma_axis( pl.subplot(2,2,1), diff_regamma_pp1, RE + r"\Gamma^{\nu, \nu'=0, \omega =0}_{s}$" ) # flip sign of w_out
pl.xlabel(r"$\nu$")
plotgamma_axis( pl.subplot(2,2,2), diff_regamma_pp1 - diff_regamma_pp1[:,:,::-1,:,:,:,:,:,:,:],  RE + r"\Gamma^{\nu, \nu'=0, \omega =0}_{t}$" ) # flip sign of w_out
pl.xlabel(r"$\nu$")
plotgamma_axis( pl.subplot(2,2,3), 2*diff_regamma_ph1 - diff_regamma_xph1, RE + r"\Gamma^{\nu, \nu'=0, \omega =0}_{d}$" )
pl.xlabel(r"$\nu$")
plotgamma_axis( pl.subplot(2,2,4), - diff_regamma_xph1, RE + r"\Gamma^{\nu, \nu'=0, \omega =0}_{m}$" )
pl.xlabel(r"$\nu$")

pl.tight_layout()

#--- Save to file
pl.savefig("plots/DeltaPHI_METHOD1_axis_nu'=0.png", dpi = 150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()

shift=20

pl.suptitle(r"$U=$" + str(UINT) + r"     $\beta=$" + str(BETA) + r"     $\epsilon=$" + str(EPS) +  r"     $\Omega=$" + str(shift) + r"$*2\pi/\beta$" + r"     Notation: $\phi(\Omega,\omega,\omega')$")

plotgamma_axis( pl.subplot(2,2,1), diff_regamma_pp1, RE + r"\Gamma^{\nu, \nu'=0, \omega =0}_{s}$" ) # flip sign of w_out
pl.xlabel(r"$\nu$")
plotgamma_axis( pl.subplot(2,2,2), diff_regamma_pp1 - diff_regamma_pp1[:,:,::-1,:,:,:,:,:,:,:],  RE + r"\Gamma^{\nu, \nu'=0, \omega =0}_{t}$" ) # flip sign of w_out
pl.xlabel(r"$\nu$")
plotgamma_axis( pl.subplot(2,2,3), 2*diff_regamma_ph1 - diff_regamma_xph1, RE + r"\Gamma^{\nu, \nu'=0, \omega =0}_{d}$" )
pl.xlabel(r"$\nu$")
plotgamma_axis( pl.subplot(2,2,4), - diff_regamma_xph1, RE + r"\Gamma^{\nu, \nu'=0, \omega =0}_{m}$" )
pl.xlabel(r"$\nu$")

pl.tight_layout()

#--- Save to file
pl.savefig("plots/DeltaPHI_METHOD1_Om="+str(shift)+"_axis_nu'=0.png", dpi = 150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()

#--- Plot along diag \nu'=\nu
shift = 0
pl.suptitle(r"$U=$" + str(UINT) + r"     $\beta=$" + str(BETA) + r"     $\epsilon=$" + str(EPS) +  r"     $\Omega=$" + str(shift) + r"$*2\pi/\beta$" + r"     Notation: $\phi(\Omega,\omega,\omega')$")

plotgamma_diag1( pl.subplot(2,2,1), diff_regamma_pp1, RE + r"\Gamma^{\nu, \nu'=0, \omega =0}_{s}$" ) # flip sign of w_out
pl.xlabel(r"$\nu$")
plotgamma_diag1( pl.subplot(2,2,2), diff_regamma_pp1 - diff_regamma_pp1[:,:,::-1,:,:,:,:,:,:,:],  RE + r"\Gamma^{\nu, \nu'=0, \omega =0}_{t}$" ) # flip sign of w_out
pl.xlabel(r"$\nu$")
plotgamma_diag1( pl.subplot(2,2,3), 2*diff_regamma_ph1 - diff_regamma_xph1, RE + r"\Gamma^{\nu, \nu'=0, \omega =0}_{d}$" )
pl.xlabel(r"$\nu$")
plotgamma_diag1( pl.subplot(2,2,4), - diff_regamma_xph1, RE + r"\Gamma^{\nu, \nu'=0, \omega =0}_{m}$" )
pl.xlabel(r"$\nu$")

pl.tight_layout()

#--- Save to file
pl.savefig("plots/DeltaPHI_METHOD1__nu=nu'.png", dpi = 150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()

shift=20

pl.suptitle(r"$U=$" + str(UINT) + r"     $\beta=$" + str(BETA) + r"     $\epsilon=$" + str(EPS) +  r"     $\Omega=$" + str(shift) + r"$*2\pi/\beta$" + r"     Notation: $\phi(\Omega,\omega,\omega')$")

plotgamma_diag1( pl.subplot(2,2,1), diff_regamma_pp1, RE + r"\Gamma^{\nu, \nu'=0, \omega =0}_{s}$" ) # flip sign of w_out
pl.xlabel(r"$\nu$")
plotgamma_diag1( pl.subplot(2,2,2), diff_regamma_pp1 - diff_regamma_pp1[:,:,::-1,:,:,:,:,:,:,:],  RE + r"\Gamma^{\nu, \nu'=0, \omega =0}_{t}$" ) # flip sign of w_out
pl.xlabel(r"$\nu$")
plotgamma_diag1( pl.subplot(2,2,3), 2*diff_regamma_ph1 - diff_regamma_xph1, RE + r"\Gamma^{\nu, \nu'=0, \omega =0}_{d}$" )
pl.xlabel(r"$\nu$")
plotgamma_diag1( pl.subplot(2,2,4), - diff_regamma_xph1, RE + r"\Gamma^{\nu, \nu'=0, \omega =0}_{m}$" )
pl.xlabel(r"$\nu$")

pl.tight_layout()

#--- Save to file
pl.savefig("plots/DeltaPHI_METHOD1_Om="+str(shift)+"_nu=nu'.png", dpi = 150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()


