#!/usr/bin/python

#===========================================================================================================================================
#
#                   NOTE: Plotting script to show the results of inverting the Bethe-Salpeter equations
#                         by using Stefan's corrections -> METHOD 1 in the main program
#
#                   CONTENTS:
#                   1) 2P irreducible vertex (Gamma) scanned nu'=0 and Omega=0
#                   2) 2P irreducible vertex (Gamma) scanned nu'= 0 and Omega \neq 0
#                   3) 2P irreducible vertex (Gamma) scanned nu'= nu and Omega=0
#                   4) 2P irreducible vertex (Gamma) scanned nu'= nu and Omega \neq 0
#                   5) 2P irreducible vertex differences for different sizes of the 00 block and 11 block, no corrections (nu'= 0, Omega=0)
#   `               6) 2P irreducible vertex differences for different sizes of the 00 block and 11 block, no corrections (nu'= 0, Omega \neq 0)
#                   7) 2P irreducible vertex differences for different sizes of the 00 block and 11 block, no corrections (nu'= nu, Omega=0)
#                   8) 2P irreducible vertex differences for different sizes of the 00 block and 11 block, no corrections (nu'= nu, Omega \neq 0)
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

fname0 = "../../dat/H5FILES/BETA50/4SITES/U1p75/POSTPROC/dat_U1p75_Beta50_PFCB1000_PARQ_SU2_METH2_INVR1xCOUNT_W0_CORR_ED.h5"
fname1 = "../../dat/H5FILES/BETA50/4SITES/U1p75/POSTPROC/dat_U1p75_Beta50_PFCB1000_PARQ_SU2_METH2_INVR1xCOUNT_CORR_ED.h5"
fname2 = "../../dat/H5FILES/BETA50/4SITES/U1p75/POSTPROC/dat_U1p75_Beta50_PFCB1000_PARQ_SU2_METH2_INVR1xCOUNT_ED.h5"
fname3 = "../../dat/H5FILES/BETA50/4SITES/U1p75/POSTPROC/dat_U1p75_Beta50_PFCB1000_PARQ_SU2_METH1_INVR1xCOUNT_ASYR10xINVR_CORR_ED.h5"

#fname0 = "../../dat/dat_U1p75_Beta50_PFCB1000_PARQ_SU2_METH2_INVR1xCOUNT_W20_ED.h5"
#fname1 = "../../dat/H5FILES/BETA50/4SITES/U1p75/POSTPROC/dat_U1p75_Beta50_PFCB1000_PARQ_SU2_METH2_INVR1xCOUNT_CORR_ED.h5"
#fname2 = "../../dat/H5FILES/BETA50/4SITES/U1p75/POSTPROC/dat_U1p75_Beta50_PFCB1000_PARQ_SU2_METH2_INVR1xCOUNT_ED.h5"
#fname3 = "../../dat/H5FILES/BETA50/4SITES/U1p75/POSTPROC/dat_U1p75_Beta50_PFCB1000_PARQ_SU2_METH1_INVR1xCOUNT_ASYR10xINVR_CORR_ED.h5"

#------ DECIDE WHICH COMPARISON YOU WANT TO ANALYZE


if len(sys.argv) > 1:
    fname0 = str(sys.argv[1])
if len(sys.argv) > 1:
    fname1 = str(sys.argv[1])
if len(sys.argv) > 1:
    fname2 = str(sys.argv[1])
if len(sys.argv) > 1:
    fname3 = str(sys.argv[1])

fname0 = fname0.rstrip('\n') # strip newline of fname
f0 = h5py.File(fname0, "r")

fname1 = fname1.rstrip('\n') # strip newline of fname
f1 = h5py.File(fname1, "r")

fname2 = fname2.rstrip('\n') # strip newline of fname
f2 = h5py.File(fname2, "r")

fname3 = fname3.rstrip('\n') # strip newline of fname
f3 = h5py.File(fname3, "r")

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

print("Plotting gamma ...")

#--- Read gamma method 2 (INVR1xVERTR)
regamma_pp0 = vert_mul * np.array(f0["/gamma_func/RE_PP"])
imgamma_pp0 = vert_mul * np.array(f0["/gamma_func/IM_PP"])
regamma_ph0 = vert_mul * np.array(f0["/gamma_func/RE_PH"])
imgamma_ph0 = vert_mul * np.array(f0["/gamma_func/IM_PH"])
regamma_xph0 = vert_mul * np.array(f0["/gamma_func/RE_XPH"])
imgamma_xph0 = vert_mul * np.array(f0["/gamma_func/IM_XPH"])

#--- Read gamma method 2 (INVR1xVERTR)
regamma_pp1 = vert_mul * np.array(f1["/gamma_func/RE_PP"])
imgamma_pp1 = vert_mul * np.array(f1["/gamma_func/IM_PP"])
regamma_ph1 = vert_mul * np.array(f1["/gamma_func/RE_PH"])
imgamma_ph1 = vert_mul * np.array(f1["/gamma_func/IM_PH"])
regamma_xph1 = vert_mul * np.array(f1["/gamma_func/RE_XPH"])
imgamma_xph1 = vert_mul * np.array(f1["/gamma_func/IM_XPH"])

#--- Read gamma method 2 (INVR1xVERTR NOCORR)
regamma_pp2 = vert_mul * np.array(f2["/gamma_func/RE_PP"])
imgamma_pp2 = vert_mul * np.array(f2["/gamma_func/IM_PP"])
regamma_ph2 = vert_mul * np.array(f2["/gamma_func/RE_PH"])
imgamma_ph2 = vert_mul * np.array(f2["/gamma_func/IM_PH"])
regamma_xph2 = vert_mul * np.array(f2["/gamma_func/RE_XPH"])
imgamma_xph2 = vert_mul * np.array(f2["/gamma_func/IM_XPH"])

#--- Read gamma method 2 (INVR1xVERTR-5)
regamma_pp3 = vert_mul * np.array(f3["/gamma_func/RE_PP"])
imgamma_pp3 = vert_mul * np.array(f3["/gamma_func/IM_PP"])
regamma_ph3 = vert_mul * np.array(f3["/gamma_func/RE_PH"])
imgamma_ph3 = vert_mul * np.array(f3["/gamma_func/IM_PH"])
regamma_xph3 = vert_mul * np.array(f3["/gamma_func/RE_XPH"])
imgamma_xph3 = vert_mul * np.array(f3["/gamma_func/IM_XPH"])

#--- RANGE FERMIONIC/BOSONIC FREQUENCIES

bdim0 = regamma_pp0.shape[0]
fdim0 = regamma_pp0.shape[1]
fdimo20 = fdim0/2
print bdim0, fdim0

bdim1 = regamma_pp1.shape[0]
fdim1 = regamma_pp1.shape[1]
fdimo21 = fdim1/2
print bdim1, fdim1

bdim2 = regamma_pp2.shape[0]
fdim2 = regamma_pp2.shape[1]
fdimo22 = fdim2/2
print bdim2, fdim2

bdim3= regamma_pp3.shape[0]
fdim3= regamma_pp3.shape[1]
fdimo23 = fdim3/2
print bdim3, fdim3

#
fdim_min = min([fdim0,fdim1,fdim2,fdim3])
print fdim_min
gammagrid_plot = np.array([(2*i+1)*pi/BETA for i in range(-fdim_min/2, fdim_min/2)])


#------------------------------------------- -PLOT DIRECT INV CORRECTIONS AND NO CORR ------------------------------------------------------------------------
shift= 0
def plotgamma_axis( use_pl, arr0, arr1, arr2,arr3, string ):
    zarr0 = np.array([ arr0[0, n+fdimo20-fdim_min/2,fdimo20,0,0,0,0,0,0,0] for n in range(fdim_min)])
    zarr1 = np.array([ arr1[shift + (bdim1-1)/2,n+fdimo21-fdim_min/2,fdimo21,0,0,0,0,0,0,0] for n in range(fdim_min)])
    zarr2 = np.array([ arr2[shift + (bdim2-1)/2,n+fdimo22-fdim_min/2,fdimo22,0,0,0,0,0,0,0] for n in range(fdim_min)])
    zarr3 = np.array([ arr3[shift + (bdim3-1)/2,n+fdimo23-fdim_min/2,fdimo23,0,0,0,0,0,0,0] for n in range(fdim_min)])
    pl.plot( gammagrid_plot, zarr0, marker = 'd', linestyle='-', color='k', ms=3, mew=0.2, label=r"INVR=120")
    pl.plot( gammagrid_plot, zarr1, marker = 'o', linestyle='-', color='r', ms=3, mew=0.2, label=r"INVR=60 METHOD2")
    pl.plot( gammagrid_plot, zarr3, marker = 's', linestyle='-', color='g', ms=3, mew=0.2, label=r"INVR=60 METHOD1")
    pl.plot( gammagrid_plot, zarr2, marker = 'x', linestyle='-', color='b', ms=3, mew=0.2, label=r"INVR=60 NOCORR")
    #pl.xlim([0.3*min(gammagrid_plot),0.3*max(gammagrid_plot)])
    use_pl.set_title( string , fontsize=10 )
    pl.legend(prop={'size':5})

    return

def plotgamma_diag1( use_pl, arr0, arr1, arr2, arr3, string ):
    zarr0 = np.array([ arr0[0, n+fdimo20-fdim_min/2,n+fdimo20-fdim_min/2,0,0,0,0,0,0,0] for n in range(fdim_min)])
    zarr1 = np.array([ arr1[shift + (bdim1-1)/2,n+fdimo21-fdim_min/2,n+fdimo21-fdim_min/2,0,0,0,0,0,0,0] for n in range(fdim_min)])
    zarr2 = np.array([ arr2[shift + (bdim2-1)/2,n+fdimo22-fdim_min/2,n+fdimo22-fdim_min/2,0,0,0,0,0,0,0] for n in range(fdim_min)])
    zarr3 = np.array([ arr3[shift + (bdim3-1)/2,n+fdimo23-fdim_min/2,n+fdimo23-fdim_min/2,0,0,0,0,0,0,0] for n in range(fdim_min)])
    pl.plot( gammagrid_plot, zarr0, marker = 'd', linestyle='-', color='k', ms=3, mew=0.2, label=r"INVR=120")
    pl.plot( gammagrid_plot, zarr1, marker = 'o', linestyle='-', color='r', ms=3, mew=0.2, label=r"INVR=60 METHOD2")
    pl.plot( gammagrid_plot, zarr3, marker = 's', linestyle='-', color='g', ms=3, mew=0.2, label=r"INVR=60 METHOD1")
    pl.plot( gammagrid_plot, zarr2, marker = 'x', linestyle='-', color='b', ms=3, mew=0.2, label=r"INVR=60 NOCORR")
    #pl.xlim([min(gammagrid_plot),max(gammagrid_plot)])
    use_pl.set_title( string , fontsize=10 )
    pl.legend(prop={'size':5})

    return

plotgamma_axis( pl.subplot(2,2,1), regamma_pp0 + regamma_pp0[:,:,::-1,:,:,:,:,:,:,:],regamma_pp1 + regamma_pp1[:,:,::-1,:,:,:,:,:,:,:], regamma_pp2 + regamma_pp2[:,:,::-1,:,:,:,:,:,:,:], regamma_pp3 + regamma_pp3[:,:,::-1,:,:,:,:,:,:,:],RE + r"\Gamma^{\nu, \nu'=0, \omega =0}_{s}$" ) # flip sign of w_out
pl.xlabel(r"$\nu$")
plotgamma_axis( pl.subplot(2,2,2), regamma_pp0 - regamma_pp0[:,:,::-1,:,:,:,:,:,:,:],regamma_pp1 - regamma_pp1[:,:,::-1,:,:,:,:,:,:,:], regamma_pp2 - regamma_pp2[:,:,::-1,:,:,:,:,:,:,:], regamma_pp3 - regamma_pp3[:,:,::-1,:,:,:,:,:,:,:],RE + r"\Gamma^{\nu, \nu'=0, \omega =0}_{t}$" ) # flip sign of w_out
pl.xlabel(r"$\nu$")
plotgamma_axis( pl.subplot(2,2,3), 2*regamma_ph0 - regamma_xph0,2*regamma_ph1 - regamma_xph1, 2*regamma_ph2 - regamma_xph2, 2*regamma_ph3 - regamma_xph3, RE + r"\Gamma^{\nu, \nu'=0, \omega =0}_{d}$" )
pl.xlabel(r"$\nu$")
plotgamma_axis( pl.subplot(2,2,4), - regamma_xph0,- regamma_xph1, - regamma_xph2, - regamma_xph3, RE + r"\Gamma^{\nu, \nu'=0, \omega =0}_{m}$" )
pl.xlabel(r"$\nu$")

pl.tight_layout()

#--- Save to file
pl.savefig("plots/Comparison_EXACTINV_METHOD2_METHOD1_axis_nu'=0.png", dpi = 150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()

plotgamma_diag1( pl.subplot(2,2,1), regamma_pp0 + regamma_pp0[:,:,::-1,:,:,:,:,:,:,:],regamma_pp1 + regamma_pp1[:,:,::-1,:,:,:,:,:,:,:], regamma_pp2 + regamma_pp2[:,:,::-1,:,:,:,:,:,:,:], regamma_pp3 + regamma_pp3[:,:,::-1,:,:,:,:,:,:,:], RE + r"\Gamma^{\nu=\nu', \omega =0}_{s}$" ) # flip sign of w_out
pl.xlabel(r"$\nu$")
plotgamma_diag1( pl.subplot(2,2,2), regamma_pp0 - regamma_pp0[:,:,::-1,:,:,:,:,:,:,:],regamma_pp1 - regamma_pp1[:,:,::-1,:,:,:,:,:,:,:], regamma_pp2 - regamma_pp2[:,:,::-1,:,:,:,:,:,:,:], regamma_pp3 - regamma_pp3[:,:,::-1,:,:,:,:,:,:,:], RE + r"\Gamma^{\nu= \nu', \omega =0}_{t}$" ) # flip sign of w_out
pl.xlabel(r"$\nu$")
plotgamma_diag1( pl.subplot(2,2,3), 2*regamma_ph0 - regamma_xph0,2*regamma_ph1 - regamma_xph1, 2*regamma_ph2 - regamma_xph2, 2*regamma_ph3 - regamma_xph3, RE + r"\Gamma^{\nu=\nu', \omega =0}_{d}$" )
pl.xlabel(r"$\nu$")
plotgamma_diag1( pl.subplot(2,2,4), - regamma_xph0,- regamma_xph1, - regamma_xph2,- regamma_xph3, RE + r"\Gamma^{\nu=\nu', \omega =0}_{m}$" )
pl.xlabel(r"$\nu$")

pl.tight_layout()

#--- Save to file
pl.savefig("plots/Comparison_EXACTINV_METHOD2_METHOD1_axis_nu=nu'.png", dpi = 150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()
#-------------------------------------------PLOT DIFFERENCE EXACTINV-INV_CORR EXACT_INV - INV_NOCORR --------------------------------------------------
shift= 0
def plotgamma_axis( use_pl, arr0, arr1, arr2, arr3, string ):
    zarr0 = np.array([ arr0[shift, n+fdimo20-fdim_min/2,fdimo20,0,0,0,0,0,0,0] - arr0[shift, n+fdimo20-fdim_min/2,fdimo20,0,0,0,0,0,0,0] for n in range(fdim_min)])
    zarr1 = np.array([ arr1[shift + (bdim1-1)/2,n+fdimo21-fdim_min/2,fdimo21,0,0,0,0,0,0,0]-arr0[shift, n+fdimo20-fdim_min/2,fdimo20,0,0,0,0,0,0,0] for n in range(fdim_min)])
    zarr2 = np.array([ arr2[shift + (bdim3-1)/2,n+fdimo23-fdim_min/2,fdimo23,0,0,0,0,0,0,0]-arr0[shift, n+fdimo20-fdim_min/2,fdimo20,0,0,0,0,0,0,0] for n in range(fdim_min)])
    zarr3= np.array([ arr3[shift + (bdim2-1)/2,n+fdimo22-fdim_min/2,fdimo22,0,0,0,0,0,0,0]-arr0[shift, n+fdimo20-fdim_min/2,fdimo20,0,0,0,0,0,0,0] for n in range(fdim_min)])
    pl.plot( gammagrid_plot, zarr0, marker = 'd', linestyle='-', color='k', ms=3, mew=0.2, label=r"$0$")
    pl.plot( gammagrid_plot, zarr1, marker = 'o', linestyle='-', color='r', ms=3, mew=0.2, label=r"INV METHOD2-EXACT")
    pl.plot( gammagrid_plot, zarr2, marker = 'x', linestyle='-', color='b', ms=3, mew=0.2, label=r"INV METHOD1-EXACT")
    pl.plot( gammagrid_plot, zarr3, marker = 's', linestyle='-', color='g', ms=3, mew=0.2, label=r"INV NOCORR-EXACT")
    pl.xlim([0.5*min(gammagrid_plot),0.5*max(gammagrid_plot)])
    use_pl.set_title( string , fontsize=10 )
    pl.legend(prop={'size':5})

    return

def plotgamma_diag1( use_pl, arr0, arr1, arr2, arr3, string ):
    zarr0 = np.array([ arr0[shift, n+fdimo20-fdim_min/2,n+fdimo20-fdim_min/2,0,0,0,0,0,0,0]-arr0[shift, n+fdimo20-fdim_min/2,n+fdimo20-fdim_min/2,0,0,0,0,0,0,0] for n in range(fdim_min)])
    zarr1 = np.array([ arr1[shift + (bdim1-1)/2,n+fdimo21-fdim_min/2,n+fdimo21-fdim_min/2,0,0,0,0,0,0,0]-arr0[shift, n+fdimo20-fdim_min/2,n+fdimo20-fdim_min/2,0,0,0,0,0,0,0] for n in range(fdim_min)])
    zarr2 = np.array([ arr2[shift + (bdim3-1)/2,n+fdimo23-fdim_min/2,n+fdimo23-fdim_min/2,0,0,0,0,0,0,0]-arr0[shift, n+fdimo20-fdim_min/2,n+fdimo20-fdim_min/2,0,0,0,0,0,0,0] for n in range(fdim_min)])
    zarr3 = np.array([ arr3[shift + (bdim2-1)/2,n+fdimo22-fdim_min/2,n+fdimo22-fdim_min/2,0,0,0,0,0,0,0]-arr0[shift, n+fdimo20-fdim_min/2,n+fdimo20-fdim_min/2,0,0,0,0,0,0,0] for n in range(fdim_min)])
    pl.plot( gammagrid_plot, zarr0, marker = 'd', linestyle='-', color='k', ms=3, mew=0.2, label=r"$0$")
    pl.plot( gammagrid_plot, zarr1, marker = 'o', linestyle='-', color='r', ms=3, mew=0.2, label=r"INV METHOD2-EXACT")
    pl.plot( gammagrid_plot, zarr2, marker = 'x', linestyle='-', color='b', ms=3, mew=0.2, label=r"INV METHOD1-EXACT")
    pl.plot( gammagrid_plot, zarr3, marker = 's', linestyle='-', color='g', ms=3, mew=0.2, label=r"INV NOCORR-EXACT")
    pl.xlim([0.5*min(gammagrid_plot),0.5*max(gammagrid_plot)])
    use_pl.set_title( string , fontsize=10 )
    pl.legend(prop={'size':5})

    return

plotgamma_diag1( pl.subplot(2,2,1), regamma_pp0 + regamma_pp0[:,:,::-1,:,:,:,:,:,:,:],regamma_pp1 + regamma_pp1[:,:,::-1,:,:,:,:,:,:,:], regamma_pp3 + regamma_pp3[:,:,::-1,:,:,:,:,:,:,:],regamma_pp2 + regamma_pp2[:,:,::-1,:,:,:,:,:,:,:], RE + r"\delta \Gamma^{\nu,\nu'=0, \omega =0}_{s}$" ) # flip sign of w_out
pl.xlabel(r"$\nu$")
plotgamma_diag1( pl.subplot(2,2,2), regamma_pp0 - regamma_pp0[:,:,::-1,:,:,:,:,:,:,:],regamma_pp1 - regamma_pp1[:,:,::-1,:,:,:,:,:,:,:], regamma_pp3 - regamma_pp3[:,:,::-1,:,:,:,:,:,:,:],regamma_pp2 - regamma_pp2[:,:,::-1,:,:,:,:,:,:,:], RE + r"\delta \Gamma^{\nu, \nu'=0, \omega =0}_{t}$" ) # flip sign of w_out
pl.xlabel(r"$\nu$")
plotgamma_diag1( pl.subplot(2,2,3), 2*regamma_ph0 - regamma_xph0,2*regamma_ph1 - regamma_xph1, 2*regamma_ph3 - regamma_xph3,2*regamma_ph2 - regamma_xph2,  RE + r"\delta \Gamma^{\nu,\nu'=0, \omega =0}_{d}$" )
pl.xlabel(r"$\nu$")
plotgamma_diag1( pl.subplot(2,2,4), - regamma_xph0,- regamma_xph1, - regamma_xph3, - regamma_xph2, RE + r"\delta \Gamma^{\nu,\nu'=0, \omega =0}_{m}$" )
pl.xlabel(r"$\nu$")

pl.tight_layout()

#--- Save to file
pl.savefig("plots/Delta_EXACTINV_METHOD2_METHOD1_axis_nu'=0.png", dpi = 150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()

plotgamma_diag1( pl.subplot(2,2,1), regamma_pp0 + regamma_pp0[:,:,::-1,:,:,:,:,:,:,:],regamma_pp1 + regamma_pp1[:,:,::-1,:,:,:,:,:,:,:], regamma_pp3 + regamma_pp3[:,:,::-1,:,:,:,:,:,:,:],regamma_pp2 + regamma_pp2[:,:,::-1,:,:,:,:,:,:,:], RE + r"\delta \Gamma^{\nu=\nu', \omega =0}_{s}$" ) # flip sign of w_out
pl.xlabel(r"$\nu$")
plotgamma_diag1( pl.subplot(2,2,2), regamma_pp0 - regamma_pp0[:,:,::-1,:,:,:,:,:,:,:],regamma_pp1 - regamma_pp1[:,:,::-1,:,:,:,:,:,:,:], regamma_pp3 - regamma_pp3[:,:,::-1,:,:,:,:,:,:,:],regamma_pp2 - regamma_pp2[:,:,::-1,:,:,:,:,:,:,:], RE + r"\delta \Gamma^{\nu= \nu', \omega =0}_{t}$" ) # flip sign of w_out
pl.xlabel(r"$\nu$")
plotgamma_diag1( pl.subplot(2,2,3), 2*regamma_ph0 - regamma_xph0,2*regamma_ph1 - regamma_xph1, 2*regamma_ph3 - regamma_xph3,2*regamma_ph2 - regamma_xph2,  RE + r"\delta \Gamma^{\nu=\nu', \omega =0}_{d}$" )
pl.xlabel(r"$\nu$")
plotgamma_diag1( pl.subplot(2,2,4), - regamma_xph0,- regamma_xph1, - regamma_xph3, - regamma_xph2, RE + r"\delta \Gamma^{\nu=\nu', \omega =0}_{m}$" )
pl.xlabel(r"$\nu$")

pl.tight_layout()

#--- Save to file
pl.savefig("plots/Delta_EXACTINV_METHOD2_METHOD1_axis_nu=nu'.png", dpi = 150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()
