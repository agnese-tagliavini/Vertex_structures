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

fname1 = "../../dat/dat_U1_Beta50_PFCB1000_PARQ_SU2_METH1_INVR1xCOUNT_ASYR5xINVR_ED.h5"
fname2 = "../../dat/dat_U1_Beta50_PFCB1000_PARQ_SU2_METH1_INVR1xCOUNT_ASYR10xINVR_ED.h5"
fname3 = "../../dat/dat_U1_Beta50_PFCB1000_PARQ_SU2_METH1_INVR1xCOUNT_ASYR15xINVR_ED.h5"
fname4 = "../../dat/dat_U1_Beta50_PFCB1000_PARQ_SU2_METH1_INVR1xCOUNT_ASYR20xINVR_ED.h5"
fname5 = "../../dat/dat_U1_Beta50_PFCB1000_PARQ_SU2_METH1_INVR1xCOUNT_ASYR25xINVR_ED.h5"
fname6 = "../../dat/dat_U1_Beta50_PFCB1000_PARQ_SU2_METH1_INVR1xCOUNT_ASYR30xINVR_ED.h5"

#------ DECIDE WHICH COMPARISON YOU WANT TO ANALYZE


if len(sys.argv) > 1:
    fname1 = str(sys.argv[1])
if len(sys.argv) > 1:
    fname2 = str(sys.argv[1])
if len(sys.argv) > 1:
    fname3 = str(sys.argv[1])
if len(sys.argv) > 1:
    fname4 = str(sys.argv[1])
if len(sys.argv) > 1:
    fname5 = str(sys.argv[1])
if len(sys.argv) > 1:
    fname6 = str(sys.argv[1])

fname1 = fname1.rstrip('\n') # strip newline of fname
f1 = h5py.File(fname1, "r")

fname2 = fname2.rstrip('\n') # strip newline of fname
f2 = h5py.File(fname2, "r")

fname3 = fname3.rstrip('\n') # strip newline of fname
f3 = h5py.File(fname3, "r")

fname4 = fname4.rstrip('\n') # strip newline of fname
f4 = h5py.File(fname4, "r")

fname5 = fname5.rstrip('\n') # strip newline of fname
f5 = h5py.File(fname5, "r")

fname6 = fname6.rstrip('\n') # strip newline of fname
f6 = h5py.File(fname6, "r")

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

#--- Read gamma method 2 (INVR1xVERTR-5 NOCORR)
regamma_pp4 = vert_mul * np.array(f4["/gamma_func/RE_PP"])
imgamma_pp4 = vert_mul * np.array(f4["/gamma_func/IM_PP"])
regamma_ph4 = vert_mul * np.array(f4["/gamma_func/RE_PH"])
imgamma_ph4 = vert_mul * np.array(f4["/gamma_func/IM_PH"])
regamma_xph4 = vert_mul * np.array(f4["/gamma_func/RE_XPH"])
imgamma_xph4 = vert_mul * np.array(f4["/gamma_func/IM_XPH"])

#--- Read gamma method 2 (INVR1xVERTR-10)
regamma_pp5 = vert_mul * np.array(f5["/gamma_func/RE_PP"])
imgamma_pp5 = vert_mul * np.array(f5["/gamma_func/IM_PP"])
regamma_ph5 = vert_mul * np.array(f5["/gamma_func/RE_PH"])
imgamma_ph5 = vert_mul * np.array(f5["/gamma_func/IM_PH"])
regamma_xph5 = vert_mul * np.array(f5["/gamma_func/RE_XPH"])
imgamma_xph5 = vert_mul * np.array(f5["/gamma_func/IM_XPH"])

#--- Read gamma method 2 (INVR1xVERTR-10 NOCORR)
regamma_pp6 = vert_mul * np.array(f6["/gamma_func/RE_PP"])
imgamma_pp6 = vert_mul * np.array(f6["/gamma_func/IM_PP"])
regamma_ph6 = vert_mul * np.array(f6["/gamma_func/RE_PH"])
imgamma_ph6 = vert_mul * np.array(f6["/gamma_func/IM_PH"])
regamma_xph6 = vert_mul * np.array(f6["/gamma_func/RE_XPH"])
imgamma_xph6 = vert_mul * np.array(f6["/gamma_func/IM_XPH"])

#--- RANGE FERMIONIC/BOSONIC FREQUENCIES

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

bdim4= regamma_pp4.shape[0]
fdim4= regamma_pp4.shape[1]
fdimo24 = fdim4/2
print bdim4, fdim4

bdim5= regamma_pp5.shape[0]
fdim5= regamma_pp5.shape[1]
fdimo25 = fdim5/2
print bdim5, fdim5

bdim6= regamma_pp6.shape[0]
fdim6= regamma_pp6.shape[1]
fdimo26 = fdim6/2
print bdim6, fdim6

#
fdim_min = min([fdim1,fdim3,fdim5,fdim7,fdim9,fdim11])
print fdim_min
gammagrid_plot = np.array([(2*i+1)*pi/BETA for i in range(-fdim_min/2, fdim_min/2)])

#---  Helper functions (include translation from  nambu to physical CAUTION : USING TIME REVERSAL SYMM) 
def plotgamma_axis( use_pl, arr1, arr2, arr3, arr4, arr5, arr6, string ):
    zarr1 = np.array([ arr1[shift + (bdim1-1)/2,n+fdimo21-fdim_min/2,fdimo21,0,0,0,0,0,0,0] for n in range(fdim_min)])
    zarr2 = np.array([ arr2[shift + (bdim2-1)/2,n+fdimo22-fdim_min/2,fdimo22,0,0,0,0,0,0,0] for n in range(fdim_min)])
    zarr3 = np.array([ arr3[shift + (bdim3-1)/2,n+fdimo23-fdim_min/2,fdimo23,0,0,0,0,0,0,0] for n in range(fdim_min)])
    zarr4 = np.array([ arr4[shift + (bdim4-1)/2,n+fdimo24-fdim_min/2,fdimo24,0,0,0,0,0,0,0] for n in range(fdim_min)])
    zarr5 = np.array([ arr5[shift + (bdim5-1)/2,n+fdimo25-fdim_min/2,fdimo25,0,0,0,0,0,0,0] for n in range(fdim_min)])
    zarr6 = np.array([ arr6[shift + (bdim6-1)/2,n+fdimo26-fdim_min/2,fdimo26,0,0,0,0,0,0,0] for n in range(fdim_min)])
    pl.plot( gammagrid_plot, zarr1, marker = 'o', linestyle='-', color='r', ms=3, mew=0.2, label=r"ASYR=5x")
    pl.plot( gammagrid_plot, zarr2, marker = 'x', linestyle='-', color='b', ms=3, mew=0.2, label=r"ASYR=10x")
    pl.plot( gammagrid_plot, zarr3, marker = '^', linestyle='-', color='g', ms=3, mew=0.2, label=r"ASYR=15x")
    pl.plot( gammagrid_plot, zarr4, marker = 'd', linestyle='-', color='c', ms=3, mew=0.2, label=r"ASYR=20x")
    pl.plot( gammagrid_plot, zarr5, marker = 's', linestyle='-', color='m', ms=3, mew=0.2, label=r"ASYR=25x")
    pl.plot( gammagrid_plot, zarr6, marker = '<', linestyle='-', color='k', ms=3, mew=0.2, label=r"ASYR=30x")
    #pl.xlim([0.3*min(gammagrid_plot),0.3*max(gammagrid_plot)])
    use_pl.set_title( string , fontsize=10 )
    pl.legend(prop={'size':5})

    return

def plotgamma_diag1( use_pl, arr1, arr2, arr3, arr4, arr5, arr6,string ):
    zarr1 = np.array([ arr1[shift + (bdim1-1)/2,n+fdimo21-fdim_min/2,n+fdimo21-fdim_min/2,0,0,0,0,0,0,0] for n in range(fdim_min)])
    zarr2 = np.array([ arr2[shift + (bdim2-1)/2,n+fdimo22-fdim_min/2,n+fdimo22-fdim_min/2,0,0,0,0,0,0,0] for n in range(fdim_min)])
    zarr3 = np.array([ arr3[shift + (bdim3-1)/2,n+fdimo23-fdim_min/2,n+fdimo23-fdim_min/2,0,0,0,0,0,0,0] for n in range(fdim_min)])
    zarr4 = np.array([ arr4[shift + (bdim4-1)/2,n+fdimo24-fdim_min/2,n+fdimo24-fdim_min/2,0,0,0,0,0,0,0] for n in range(fdim_min)])
    zarr5 = np.array([ arr5[shift + (bdim5-1)/2,n+fdimo25-fdim_min/2,n+fdimo25-fdim_min/2,0,0,0,0,0,0,0] for n in range(fdim_min)])
    zarr6 = np.array([ arr6[shift + (bdim6-1)/2,n+fdimo26-fdim_min/2,n+fdimo26-fdim_min/2,0,0,0,0,0,0,0] for n in range(fdim_min)])
    pl.plot( gammagrid_plot, zarr1, marker = 'o', linestyle='-', color='r', ms=3, mew=0.2, label=r"ASYR=5x")
    pl.plot( gammagrid_plot, zarr2, marker = 'x', linestyle='-', color='b', ms=3, mew=0.2, label=r"ASYR=10x")
    pl.plot( gammagrid_plot, zarr3, marker = '^', linestyle='-', color='g', ms=3, mew=0.2, label=r"ASYR=15x")
    pl.plot( gammagrid_plot, zarr4, marker = 'd', linestyle='-', color='c', ms=3, mew=0.2, label=r"ASYR=20x")
    pl.plot( gammagrid_plot, zarr5, marker = 's', linestyle='-', color='m', ms=3, mew=0.2, label=r"ASYR=25x")
    pl.plot( gammagrid_plot, zarr6, marker = '<', linestyle='-', color='k', ms=3, mew=0.2, label=r"ASYR=30x")
    #pl.xlim([min(gammagrid_plot),max(gammagrid_plot)])
    use_pl.set_title( string , fontsize=10 )
    pl.legend(prop={'size':5})

    return

#--- Plot along \nu'=0 

#pl.suptitle(r"$U=$" + str(UINT) + r"     $\beta=$" + str(BETA) + r"     $\epsilon=$" + str(EPS) +  r"     $\Omega=$" + str(shift) + r"$*2\pi/\beta$")

plotgamma_axis( pl.subplot(2,2,1), regamma_pp1 + regamma_pp1[:,:,::-1,:,:,:,:,:,:,:], regamma_pp2 + regamma_pp2[:,:,::-1,:,:,:,:,:,:,:], regamma_pp3 + regamma_pp3[:,:,::-1,:,:,:,:,:,:,:], regamma_pp4 + regamma_pp4[:,:,::-1,:,:,:,:,:,:,:],regamma_pp5 + regamma_pp5[:,:,::-1,:,:,:,:,:,:,:],regamma_pp6 + regamma_pp6[:,:,::-1,:,:,:,:,:,:,:], RE + r"\Gamma^{\nu, \nu'=0, \omega =0}_{s}$" ) # flip sign of w_out
pl.xlabel(r"$\nu$")
plotgamma_axis( pl.subplot(2,2,2), regamma_pp1 - regamma_pp1[:,:,::-1,:,:,:,:,:,:,:], regamma_pp2 - regamma_pp2[:,:,::-1,:,:,:,:,:,:,:], regamma_pp3 - regamma_pp3[:,:,::-1,:,:,:,:,:,:,:], regamma_pp4 - regamma_pp4[:,:,::-1,:,:,:,:,:,:,:], regamma_pp5 - regamma_pp5[:,:,::-1,:,:,:,:,:,:,:],regamma_pp6 - regamma_pp6[:,:,::-1,:,:,:,:,:,:,:], RE + r"\Gamma^{\nu, \nu'=0, \omega =0}_{t}$" ) # flip sign of w_out
pl.xlabel(r"$\nu$")
plotgamma_axis( pl.subplot(2,2,3), 2*regamma_ph1 - regamma_xph1, 2*regamma_ph2 - regamma_xph2, 2*regamma_ph3 - regamma_xph3, 2*regamma_ph4 - regamma_xph4,2*regamma_ph5 - regamma_xph5,2*regamma_ph6 - regamma_xph6, RE + r"\Gamma^{\nu, \nu'=0, \omega =0}_{d}$" )
pl.xlabel(r"$\nu$")
plotgamma_axis( pl.subplot(2,2,4), - regamma_xph1, - regamma_xph2, - regamma_xph3, - regamma_xph4,- regamma_xph5,- regamma_xph6, RE + r"\Gamma^{\nu, \nu'=0, \omega =0}_{m}$" )
pl.xlabel(r"$\nu$")

pl.tight_layout()

#--- Save to file
pl.savefig("plots/Gammas_METHOD1_axis_nu'=0.png", dpi = 150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()

shift=20

plotgamma_axis( pl.subplot(2,2,1), regamma_pp1 + regamma_pp1[:,:,::-1,:,:,:,:,:,:,:], regamma_pp2 + regamma_pp2[:,:,::-1,:,:,:,:,:,:,:], regamma_pp3 + regamma_pp3[:,:,::-1,:,:,:,:,:,:,:], regamma_pp4 + regamma_pp4[:,:,::-1,:,:,:,:,:,:,:],regamma_pp5 + regamma_pp5[:,:,::-1,:,:,:,:,:,:,:],regamma_pp6 + regamma_pp6[:,:,::-1,:,:,:,:,:,:,:], RE + r"\Gamma^{\nu, \nu'=0, \omega =20}_{s}$" ) # flip sign of w_out
pl.xlabel(r"$\nu$")
plotgamma_axis( pl.subplot(2,2,2), regamma_pp1 - regamma_pp1[:,:,::-1,:,:,:,:,:,:,:], regamma_pp2 - regamma_pp2[:,:,::-1,:,:,:,:,:,:,:], regamma_pp3 - regamma_pp3[:,:,::-1,:,:,:,:,:,:,:], regamma_pp4 - regamma_pp4[:,:,::-1,:,:,:,:,:,:,:], regamma_pp5 - regamma_pp5[:,:,::-1,:,:,:,:,:,:,:],regamma_pp6 - regamma_pp6[:,:,::-1,:,:,:,:,:,:,:], RE + r"\Gamma^{\nu, \nu'=0, \omega =20}_{t}$" ) # flip sign of w_out
pl.xlabel(r"$\nu$")
plotgamma_axis( pl.subplot(2,2,3), 2*regamma_ph1 - regamma_xph1, 2*regamma_ph2 - regamma_xph2, 2*regamma_ph3 - regamma_xph3, 2*regamma_ph4 - regamma_xph4,2*regamma_ph5 - regamma_xph5,2*regamma_ph6 - regamma_xph6, RE + r"\Gamma^{\nu, \nu'=0, \omega =20}_{d}$" )
pl.xlabel(r"$\nu$")
plotgamma_axis( pl.subplot(2,2,4), - regamma_xph1, - regamma_xph2, - regamma_xph3, - regamma_xph4,- regamma_xph5,- regamma_xph6, RE + r"\Gamma^{\nu, \nu'=0, \omega =20}_{m}$" )
pl.xlabel(r"$\nu$")

pl.tight_layout()

#--- Save to file
pl.savefig("plots/Gammas_METHOD1_Om="+str(shift)+"_axis_nu'=0.png", dpi = 150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()

shift=0
#--- Plot along \nu'=\nu 

plotgamma_diag1( pl.subplot(2,2,1), regamma_pp1 + regamma_pp1[:,:,::-1,:,:,:,:,:,:,:], regamma_pp2 + regamma_pp2[:,:,::-1,:,:,:,:,:,:,:], regamma_pp3 + regamma_pp3[:,:,::-1,:,:,:,:,:,:,:], regamma_pp4 + regamma_pp4[:,:,::-1,:,:,:,:,:,:,:],regamma_pp5 + regamma_pp5[:,:,::-1,:,:,:,:,:,:,:],regamma_pp6 + regamma_pp6[:,:,::-1,:,:,:,:,:,:,:], RE + r"\Gamma^{\nu, \nu', \omega =0}_{s}$" ) # flip sign of w_out
pl.xlabel(r"$\nu$")
plotgamma_diag1( pl.subplot(2,2,2), regamma_pp1 - regamma_pp1[:,:,::-1,:,:,:,:,:,:,:], regamma_pp2 - regamma_pp2[:,:,::-1,:,:,:,:,:,:,:], regamma_pp3 - regamma_pp3[:,:,::-1,:,:,:,:,:,:,:], regamma_pp4 - regamma_pp4[:,:,::-1,:,:,:,:,:,:,:], regamma_pp5 - regamma_pp5[:,:,::-1,:,:,:,:,:,:,:],regamma_pp6 - regamma_pp6[:,:,::-1,:,:,:,:,:,:,:], RE + r"\Gamma^{\nu, \nu', \omega =0}_{t}$" ) # flip sign of w_out
pl.xlabel(r"$\nu$")
plotgamma_diag1( pl.subplot(2,2,3), 2*regamma_ph1 - regamma_xph1, 2*regamma_ph2 - regamma_xph2, 2*regamma_ph3 - regamma_xph3, 2*regamma_ph4 - regamma_xph4,2*regamma_ph5 - regamma_xph5,2*regamma_ph6 - regamma_xph6, RE + r"\Gamma^{\nu, \nu', \omega =0}_{d}$" )
pl.xlabel(r"$\nu$")
plotgamma_diag1( pl.subplot(2,2,4), - regamma_xph1, - regamma_xph2, - regamma_xph3, - regamma_xph4,- regamma_xph5,- regamma_xph6, RE + r"\Gamma^{\nu, \nu', \omega =0}_{m}$" )
pl.xlabel(r"$\nu$")

pl.tight_layout()

#--- Save to file
pl.savefig("plots/Gammas_METHOD1_nu=nu'.png", dpi = 150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()

shift=20

plotgamma_diag1( pl.subplot(2,2,1), regamma_pp1 + regamma_pp1[:,:,::-1,:,:,:,:,:,:,:], regamma_pp2 + regamma_pp2[:,:,::-1,:,:,:,:,:,:,:], regamma_pp3 + regamma_pp3[:,:,::-1,:,:,:,:,:,:,:], regamma_pp4 + regamma_pp4[:,:,::-1,:,:,:,:,:,:,:],regamma_pp5 + regamma_pp5[:,:,::-1,:,:,:,:,:,:,:],regamma_pp6 + regamma_pp6[:,:,::-1,:,:,:,:,:,:,:], RE + r"\Gamma^{\nu, \nu', \omega =20}_{s}$" ) # flip sign of w_out
pl.xlabel(r"$\nu$")
plotgamma_diag1( pl.subplot(2,2,2), regamma_pp1 - regamma_pp1[:,:,::-1,:,:,:,:,:,:,:], regamma_pp2 - regamma_pp2[:,:,::-1,:,:,:,:,:,:,:], regamma_pp3 - regamma_pp3[:,:,::-1,:,:,:,:,:,:,:], regamma_pp4 - regamma_pp4[:,:,::-1,:,:,:,:,:,:,:], regamma_pp5 - regamma_pp5[:,:,::-1,:,:,:,:,:,:,:],regamma_pp6 - regamma_pp6[:,:,::-1,:,:,:,:,:,:,:], RE + r"\Gamma^{\nu, \nu', \omega =20}_{t}$" ) # flip sign of w_out
pl.xlabel(r"$\nu$")
plotgamma_diag1( pl.subplot(2,2,3), 2*regamma_ph1 - regamma_xph1, 2*regamma_ph2 - regamma_xph2, 2*regamma_ph3 - regamma_xph3, 2*regamma_ph4 - regamma_xph4,2*regamma_ph5 - regamma_xph5,2*regamma_ph6 - regamma_xph6, RE + r"\Gamma^{\nu, \nu', \omega =20}_{d}$" )
pl.xlabel(r"$\nu$")
plotgamma_diag1( pl.subplot(2,2,4), - regamma_xph1, - regamma_xph2, - regamma_xph3, - regamma_xph4,- regamma_xph5,- regamma_xph6, RE + r"\Gamma^{\nu, \nu', \omega =20}_{m}$" )
pl.xlabel(r"$\nu$")



pl.tight_layout()

#--- Save to file
pl.savefig("plots/Gammas_METHOD1_Om="+str(shift)+"_nu=nu'.png", dpi = 150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()

shift= 0
#-------------------------------------GAMMA PLOTTING DIFFERENCES NOCORRECTIONS MINUS CORRECTION RESULTS ------------------------------------------

print("Plotting delta gamma ...")

#--- INVR = COUNT 
diff_regamma_pp1 = vert_mul * np.subtract( f1["/gamma_func/RE_PP"], f6["/gamma_func/RE_PP"])
diff_imgamma_pp1 = vert_mul * np.subtract( f1["/gamma_func/IM_PP"], f6["/gamma_func/IM_PP"] )
diff_regamma_ph1 = vert_mul * np.subtract( f1["/gamma_func/RE_PH"], f6["/gamma_func/RE_PH"] )
diff_imgamma_ph1 = vert_mul * np.subtract( f1["/gamma_func/IM_PH"], f6["/gamma_func/IM_PH"] )
diff_regamma_xph1 = vert_mul * np.subtract(f1["/gamma_func/RE_XPH"],f6["/gamma_func/RE_XPH"])
diff_imgamma_xph1 = vert_mul * np.subtract(f1["/gamma_func/IM_XPH"],f6["/gamma_func/IM_XPH"])

#--- INVR=COUNT-5
diff_regamma_pp2 = vert_mul * np.subtract( f2["/gamma_func/RE_PP"], f6["/gamma_func/RE_PP"])
diff_imgamma_pp2 = vert_mul * np.subtract( f2["/gamma_func/IM_PP"], f6["/gamma_func/IM_PP"] )
diff_regamma_ph2 = vert_mul * np.subtract( f2["/gamma_func/RE_PH"], f6["/gamma_func/RE_PH"] )
diff_imgamma_ph2 = vert_mul * np.subtract( f2["/gamma_func/IM_PH"], f6["/gamma_func/IM_PH"] )
diff_regamma_xph2 = vert_mul * np.subtract(f2["/gamma_func/RE_XPH"],f6["/gamma_func/RE_XPH"])
diff_imgamma_xph2 = vert_mul * np.subtract(f2["/gamma_func/IM_XPH"],f6["/gamma_func/IM_XPH"])

#--- INVR=COUNT-10
diff_regamma_pp3 = vert_mul * np.subtract( f3["/gamma_func/RE_PP"], f6["/gamma_func/RE_PP"])
diff_imgamma_pp3 = vert_mul * np.subtract( f3["/gamma_func/IM_PP"], f6["/gamma_func/IM_PP"] )
diff_regamma_ph3 = vert_mul * np.subtract( f3["/gamma_func/RE_PH"], f6["/gamma_func/RE_PH"] )
diff_imgamma_ph3 = vert_mul * np.subtract( f3["/gamma_func/IM_PH"], f6["/gamma_func/IM_PH"] )
diff_regamma_xph3 = vert_mul * np.subtract(f3["/gamma_func/RE_XPH"],f6["/gamma_func/RE_XPH"])
diff_imgamma_xph3 = vert_mul * np.subtract(f3["/gamma_func/IM_XPH"],f6["/gamma_func/IM_XPH"])

#--- INVR=COUNT-15
diff_regamma_pp4 = vert_mul * np.subtract( f4["/gamma_func/RE_PP"], f6["/gamma_func/RE_PP"])
diff_imgamma_pp4 = vert_mul * np.subtract( f4["/gamma_func/IM_PP"], f6["/gamma_func/IM_PP"] )
diff_regamma_ph4 = vert_mul * np.subtract( f4["/gamma_func/RE_PH"], f6["/gamma_func/RE_PH"] )
diff_imgamma_ph4 = vert_mul * np.subtract( f4["/gamma_func/IM_PH"], f6["/gamma_func/IM_PH"] )
diff_regamma_xph4 = vert_mul * np.subtract(f4["/gamma_func/RE_XPH"],f6["/gamma_func/RE_XPH"])
diff_imgamma_xph4 = vert_mul * np.subtract(f4["/gamma_func/IM_XPH"],f6["/gamma_func/IM_XPH"])

#--- INVR=COUNT-20
diff_regamma_pp5 = vert_mul * np.subtract( f5["/gamma_func/RE_PP"], f6["/gamma_func/RE_PP"])
diff_imgamma_pp5 = vert_mul * np.subtract( f5["/gamma_func/IM_PP"], f6["/gamma_func/IM_PP"] )
diff_regamma_ph5 = vert_mul * np.subtract( f5["/gamma_func/RE_PH"], f6["/gamma_func/RE_PH"] )
diff_imgamma_ph5 = vert_mul * np.subtract( f5["/gamma_func/IM_PH"], f6["/gamma_func/IM_PH"] )
diff_regamma_xph5 = vert_mul * np.subtract(f5["/gamma_func/RE_XPH"],f6["/gamma_func/RE_XPH"])
diff_imgamma_xph5 = vert_mul * np.subtract(f5["/gamma_func/IM_XPH"],f6["/gamma_func/IM_XPH"])

diff_regamma_pp6 = vert_mul * np.subtract( f6["/gamma_func/RE_PP"], f6["/gamma_func/RE_PP"])
diff_imgamma_pp6 = vert_mul * np.subtract( f6["/gamma_func/IM_PP"], f6["/gamma_func/IM_PP"] )
diff_regamma_ph6 = vert_mul * np.subtract( f6["/gamma_func/RE_PH"], f6["/gamma_func/RE_PH"] )
diff_imgamma_ph6 = vert_mul * np.subtract( f6["/gamma_func/IM_PH"], f6["/gamma_func/IM_PH"] )
diff_regamma_xph6 = vert_mul * np.subtract(f6["/gamma_func/RE_XPH"],f6["/gamma_func/RE_XPH"])
diff_imgamma_xph6 = vert_mul * np.subtract(f6["/gamma_func/IM_XPH"],f6["/gamma_func/IM_XPH"])

#---  Helper functions (include translation from  nambu to physical CAUTION : USING TIME REVERSAL SYMM) 
def plotgamma_axis( use_pl, arr1, arr2, arr3, arr4, arr5,arr6, string ):
    zarr1 = np.array([ arr1[shift + (bdim1-1)/2,n+fdimo21-fdim_min/2,fdimo21,0,0,0,0,0,0,0] for n in range(fdim_min)])
    zarr2 = np.array([ arr2[shift + (bdim2-1)/2,n+fdimo22-fdim_min/2,fdimo22,0,0,0,0,0,0,0] for n in range(fdim_min)])
    zarr3 = np.array([ arr3[shift + (bdim3-1)/2,n+fdimo23-fdim_min/2,fdimo23,0,0,0,0,0,0,0] for n in range(fdim_min)])
    zarr4 = np.array([ arr4[shift + (bdim4-1)/2,n+fdimo24-fdim_min/2,fdimo24,0,0,0,0,0,0,0] for n in range(fdim_min)])
    zarr5 = np.array([ arr5[shift + (bdim5-1)/2,n+fdimo25-fdim_min/2,fdimo25,0,0,0,0,0,0,0] for n in range(fdim_min)])
    zarr6 = np.array([ arr6[shift + (bdim6-1)/2,n+fdimo26-fdim_min/2,fdimo26,0,0,0,0,0,0,0] for n in range(fdim_min)])
    pl.plot( gammagrid_plot, zarr1, marker = 'o', linestyle='-', color='r', ms=3, mew=0.2, label=r"ASYR=5x-30x")
    pl.plot( gammagrid_plot, zarr2, marker = 'x', linestyle='-', color='b', ms=3, mew=0.2, label=r"ASYR=10x-30x")
    pl.plot( gammagrid_plot, zarr3, marker = '^', linestyle='-', color='g', ms=3, mew=0.2, label=r"ASYR=15x-30x")
    pl.plot( gammagrid_plot, zarr4, marker = 'd', linestyle='-', color='c', ms=3, mew=0.2, label=r"ASYR=20x-30x")
    pl.plot( gammagrid_plot, zarr5, marker = 's', linestyle='-', color='m', ms=3, mew=0.2, label=r"ASYR=25x-30x")
    pl.plot( gammagrid_plot, zarr6, marker = '<', linestyle='-', color='k', ms=3, mew=0.2, label=r"ASYR=30x-30x")
    pl.xlim([0.3*min(gammagrid_plot),0.3*max(gammagrid_plot)])
    use_pl.set_title( string , fontsize=10 )
    pl.legend(prop={'size':5})

    return

def plotgamma_diag1( use_pl, arr1, arr2, arr3, arr4, arr5, arr6, string ):
    zarr1 = np.array([ arr1[shift + (bdim1-1)/2,n+fdimo21-fdim_min/2,n+fdimo21-fdim_min/2,0,0,0,0,0,0,0] for n in range(fdim_min)])
    zarr2 = np.array([ arr2[shift + (bdim2-1)/2,n+fdimo22-fdim_min/2,n+fdimo22-fdim_min/2,0,0,0,0,0,0,0] for n in range(fdim_min)])
    zarr3 = np.array([ arr3[shift + (bdim3-1)/2,n+fdimo23-fdim_min/2,n+fdimo23-fdim_min/2,0,0,0,0,0,0,0] for n in range(fdim_min)])
    zarr4 = np.array([ arr4[shift + (bdim4-1)/2,n+fdimo24-fdim_min/2,n+fdimo24-fdim_min/2,0,0,0,0,0,0,0] for n in range(fdim_min)])
    zarr5 = np.array([ arr5[shift + (bdim5-1)/2,n+fdimo25-fdim_min/2,n+fdimo25-fdim_min/2,0,0,0,0,0,0,0] for n in range(fdim_min)])
    zarr6 = np.array([ arr6[shift + (bdim6-1)/2,n+fdimo26-fdim_min/2,n+fdimo26-fdim_min/2,0,0,0,0,0,0,0] for n in range(fdim_min)])
    pl.plot( gammagrid_plot, zarr1, marker = 'o', linestyle='-', color='r', ms=3, mew=0.2, label=r"ASYR=5x-30x")
    pl.plot( gammagrid_plot, zarr2, marker = 'x', linestyle='-', color='b', ms=3, mew=0.2, label=r"ASYR=10x-30x")
    pl.plot( gammagrid_plot, zarr3, marker = '^', linestyle='-', color='g', ms=3, mew=0.2, label=r"ASYR=15x-30x")
    pl.plot( gammagrid_plot, zarr4, marker = 'd', linestyle='-', color='c', ms=3, mew=0.2, label=r"ASYR=20x-30x")
    pl.plot( gammagrid_plot, zarr5, marker = 's', linestyle='-', color='m', ms=3, mew=0.2, label=r"ASYR=25x-30x")
    pl.plot( gammagrid_plot, zarr6, marker = '<', linestyle='-', color='k', ms=3, mew=0.2, label=r"ASYR=30x-30x")
    pl.xlim([0.3*min(gammagrid_plot),0.3*max(gammagrid_plot)])
    use_pl.set_title( string , fontsize=10 )
    pl.legend(prop={'size':5})

    return

#--- Plot along axis \nu'=0

#pl.suptitle(r"$U=$" + str(UINT) + r"     $\beta=$" + str(BETA) + r"     $\epsilon=$" + str(EPS) +  r"     $\Omega=$" + str(shift) + r"$*2\pi/\beta$" + r"     Notation: $\phi(\Omega,\omega,\omega')$")

plotgamma_axis( pl.subplot(2,2,1), diff_regamma_pp1 + diff_regamma_pp1[:,:,::-1,:,:,:,:,:,:,:], diff_regamma_pp2 + diff_regamma_pp2[:,:,::-1,:,:,:,:,:,:,:], diff_regamma_pp3 + diff_regamma_pp3[:,:,::-1,:,:,:,:,:,:,:], diff_regamma_pp4 + diff_regamma_pp4[:,:,::-1,:,:,:,:,:,:,:],diff_regamma_pp5 + diff_regamma_pp5[:,:,::-1,:,:,:,:,:,:,:],diff_regamma_pp6 + diff_regamma_pp6[:,:,::-1,:,:,:,:,:,:,:],RE + r"\delta \Gamma^{\nu, \nu'=0, \omega =0}_{s}$" ) # flip sign of w_out
pl.xlabel(r"$\nu$")
plotgamma_axis( pl.subplot(2,2,2), diff_regamma_pp1 - diff_regamma_pp1[:,:,::-1,:,:,:,:,:,:,:], diff_regamma_pp2 - diff_regamma_pp2[:,:,::-1,:,:,:,:,:,:,:], diff_regamma_pp3 - diff_regamma_pp3[:,:,::-1,:,:,:,:,:,:,:], diff_regamma_pp4 - diff_regamma_pp4[:,:,::-1,:,:,:,:,:,:,:], diff_regamma_pp5 - diff_regamma_pp5[:,:,::-1,:,:,:,:,:,:,:],diff_regamma_pp6 - diff_regamma_pp6[:,:,::-1,:,:,:,:,:,:,:], RE + r"\delta \Gamma^{\nu, \nu'=0, \omega =0}_{t}$" ) # flip sign of w_out
pl.xlabel(r"$\nu$")
plotgamma_axis( pl.subplot(2,2,3), 2*diff_regamma_ph1 - diff_regamma_xph1, 2*diff_regamma_ph2 - diff_regamma_xph2, 2*diff_regamma_ph3 - diff_regamma_xph3, 2*diff_regamma_ph4 - diff_regamma_xph4,2*diff_regamma_ph5 - diff_regamma_xph5,2*diff_regamma_ph6 - diff_regamma_xph6, RE + r"\delta \Gamma^{\nu, \nu'=0, \omega =0}_{d}$" )
pl.xlabel(r"$\nu$")
plotgamma_axis( pl.subplot(2,2,4), - diff_regamma_xph1, - diff_regamma_xph2, - diff_regamma_xph3, - diff_regamma_xph4,- diff_regamma_xph5,- diff_regamma_xph6, RE + r"\delta \Gamma^{\nu, \nu'=0, \omega =0}_{m}$" )
pl.xlabel(r"$\nu$")

pl.tight_layout()

#--- Save to file
pl.savefig("plots/DeltaGammas_METHOD1_axis_nu'=0.png", dpi = 150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()

shift=20

#pl.suptitle(r"$U=$" + str(UINT) + r"     $\beta=$" + str(BETA) + r"     $\epsilon=$" + str(EPS) +  r"     $\Omega=$" + str(shift) + r"$*2\pi/\beta$" + r"     Notation: $\phi(\Omega,\omega,\omega')$")
plotgamma_axis( pl.subplot(2,2,1), diff_regamma_pp1 + diff_regamma_pp1[:,:,::-1,:,:,:,:,:,:,:], diff_regamma_pp2 + diff_regamma_pp2[:,:,::-1,:,:,:,:,:,:,:], diff_regamma_pp3 + diff_regamma_pp3[:,:,::-1,:,:,:,:,:,:,:], diff_regamma_pp4 + diff_regamma_pp4[:,:,::-1,:,:,:,:,:,:,:],diff_regamma_pp5 + diff_regamma_pp5[:,:,::-1,:,:,:,:,:,:,:],diff_regamma_pp6 + diff_regamma_pp6[:,:,::-1,:,:,:,:,:,:,:],RE + r"\delta \Gamma^{\nu, \nu'=0, \omega =20}_{s}$" ) # flip sign of w_out
pl.xlabel(r"$\nu$")
plotgamma_axis( pl.subplot(2,2,2), diff_regamma_pp1 - diff_regamma_pp1[:,:,::-1,:,:,:,:,:,:,:], diff_regamma_pp2 - diff_regamma_pp2[:,:,::-1,:,:,:,:,:,:,:], diff_regamma_pp3 - diff_regamma_pp3[:,:,::-1,:,:,:,:,:,:,:], diff_regamma_pp4 - diff_regamma_pp4[:,:,::-1,:,:,:,:,:,:,:], diff_regamma_pp5 - diff_regamma_pp5[:,:,::-1,:,:,:,:,:,:,:],diff_regamma_pp6 - diff_regamma_pp6[:,:,::-1,:,:,:,:,:,:,:], RE + r"\delta \Gamma^{\nu, \nu'=0, \omega =20}_{t}$" ) # flip sign of w_out
pl.xlabel(r"$\nu$")
plotgamma_axis( pl.subplot(2,2,3), 2*diff_regamma_ph1 - diff_regamma_xph1, 2*diff_regamma_ph2 - diff_regamma_xph2, 2*diff_regamma_ph3 - diff_regamma_xph3, 2*diff_regamma_ph4 - diff_regamma_xph4,2*diff_regamma_ph5 - diff_regamma_xph5,2*diff_regamma_ph6 - diff_regamma_xph6, RE + r"\delta \Gamma^{\nu, \nu'=0, \omega =20}_{d}$" )
pl.xlabel(r"$\nu$")
plotgamma_axis( pl.subplot(2,2,4), - diff_regamma_xph1, - diff_regamma_xph2, - diff_regamma_xph3, - diff_regamma_xph4,- diff_regamma_xph5,- diff_regamma_xph6, RE + r"\delta \Gamma^{\nu, \nu'=0, \omega =20}_{m}$" )
pl.xlabel(r"$\nu$")


pl.tight_layout()

#--- Save to file
pl.savefig("plots/DeltaGammas_METHOD1_Om="+str(shift)+"_axis_nu'=0.png", dpi = 150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()

shift=0
#--- Plot along diag \nu'=\nu

#pl.suptitle(r"$U=$" + str(UINT) + r"     $\beta=$" + str(BETA) + r"     $\epsilon=$" + str(EPS) +  r"     $\Omega=$" + str(shift) + r"$*2\pi/\beta$" + r"     Notation: $\phi(\Omega,\omega,\omega')$")
plotgamma_diag1( pl.subplot(2,2,1), diff_regamma_pp1 + diff_regamma_pp1[:,:,::-1,:,:,:,:,:,:,:], diff_regamma_pp2 + diff_regamma_pp2[:,:,::-1,:,:,:,:,:,:,:], diff_regamma_pp3 + diff_regamma_pp3[:,:,::-1,:,:,:,:,:,:,:], diff_regamma_pp4 + diff_regamma_pp4[:,:,::-1,:,:,:,:,:,:,:],diff_regamma_pp5 + diff_regamma_pp5[:,:,::-1,:,:,:,:,:,:,:],diff_regamma_pp6 + diff_regamma_pp6[:,:,::-1,:,:,:,:,:,:,:],RE + r"\delta \Gamma^{\nu=\nu', \omega =0}_{s}$" ) # flip sign of w_out
pl.xlabel(r"$\nu$")
plotgamma_diag1( pl.subplot(2,2,2), diff_regamma_pp1 - diff_regamma_pp1[:,:,::-1,:,:,:,:,:,:,:], diff_regamma_pp2 - diff_regamma_pp2[:,:,::-1,:,:,:,:,:,:,:], diff_regamma_pp3 - diff_regamma_pp3[:,:,::-1,:,:,:,:,:,:,:], diff_regamma_pp4 - diff_regamma_pp4[:,:,::-1,:,:,:,:,:,:,:], diff_regamma_pp5 - diff_regamma_pp5[:,:,::-1,:,:,:,:,:,:,:],diff_regamma_pp6 - diff_regamma_pp6[:,:,::-1,:,:,:,:,:,:,:], RE + r"\delta \Gamma^{\nu= \nu', \omega =0}_{t}$" ) # flip sign of w_out
pl.xlabel(r"$\nu$")
plotgamma_diag1( pl.subplot(2,2,3), 2*diff_regamma_ph1 - diff_regamma_xph1, 2*diff_regamma_ph2 - diff_regamma_xph2, 2*diff_regamma_ph3 - diff_regamma_xph3, 2*diff_regamma_ph4 - diff_regamma_xph4,2*diff_regamma_ph5 - diff_regamma_xph5,2*diff_regamma_ph6 - diff_regamma_xph6, RE + r"\delta \Gamma^{\nu=\nu', \omega =0}_{d}$" )
pl.xlabel(r"$\nu$")
plotgamma_diag1( pl.subplot(2,2,4), - diff_regamma_xph1, - diff_regamma_xph2, - diff_regamma_xph3, - diff_regamma_xph4,- diff_regamma_xph5,- diff_regamma_xph6, RE + r"\delta \Gamma^{\nu=\nu', \omega =0}_{m}$" )
pl.xlabel(r"$\nu$")

pl.tight_layout()

#--- Save to file
pl.savefig("plots/DeltaGammas_METHOD1__nu=nu'.png", dpi = 150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()

shift=20

#pl.suptitle(r"$U=$" + str(UINT) + r"     $\beta=$" + str(BETA) + r"     $\epsilon=$" + str(EPS) +  r"     $\Omega=$" + str(shift) + r"$*2\pi/\beta$" + r"     Notation: $\phi(\Omega,\omega,\omega')$")
plotgamma_diag1( pl.subplot(2,2,1), diff_regamma_pp1 + diff_regamma_pp1[:,:,::-1,:,:,:,:,:,:,:], diff_regamma_pp2 + diff_regamma_pp2[:,:,::-1,:,:,:,:,:,:,:], diff_regamma_pp3 + diff_regamma_pp3[:,:,::-1,:,:,:,:,:,:,:], diff_regamma_pp4 + diff_regamma_pp4[:,:,::-1,:,:,:,:,:,:,:],diff_regamma_pp5 + diff_regamma_pp5[:,:,::-1,:,:,:,:,:,:,:],diff_regamma_pp6 + diff_regamma_pp6[:,:,::-1,:,:,:,:,:,:,:],RE + r"\delta \Gamma^{\nu=\nu', \omega =20}_{s}$" ) # flip sign of w_out
pl.xlabel(r"$\nu$")
plotgamma_diag1( pl.subplot(2,2,2), diff_regamma_pp1 - diff_regamma_pp1[:,:,::-1,:,:,:,:,:,:,:], diff_regamma_pp2 - diff_regamma_pp2[:,:,::-1,:,:,:,:,:,:,:], diff_regamma_pp3 - diff_regamma_pp3[:,:,::-1,:,:,:,:,:,:,:], diff_regamma_pp4 - diff_regamma_pp4[:,:,::-1,:,:,:,:,:,:,:], diff_regamma_pp5 - diff_regamma_pp5[:,:,::-1,:,:,:,:,:,:,:],diff_regamma_pp6 - diff_regamma_pp6[:,:,::-1,:,:,:,:,:,:,:], RE + r"\delta \Gamma^{\nu= \nu', \omega =20}_{t}$" ) # flip sign of w_out
pl.xlabel(r"$\nu$")
plotgamma_diag1( pl.subplot(2,2,3), 2*diff_regamma_ph1 - diff_regamma_xph1, 2*diff_regamma_ph2 - diff_regamma_xph2, 2*diff_regamma_ph3 - diff_regamma_xph3, 2*diff_regamma_ph4 - diff_regamma_xph4,2*diff_regamma_ph5 - diff_regamma_xph5,2*diff_regamma_ph6 - diff_regamma_xph6, RE + r"\delta \Gamma^{\nu=\nu', \omega =20}_{d}$" )
pl.xlabel(r"$\nu$")
plotgamma_diag1( pl.subplot(2,2,4), - diff_regamma_xph1, - diff_regamma_xph2, - diff_regamma_xph3, - diff_regamma_xph4,- diff_regamma_xph5,- diff_regamma_xph6, RE + r"\delta \Gamma^{\nu=\nu', \omega =20}_{m}$" )
pl.xlabel(r"$\nu$")

pl.tight_layout()


#--- Save to file
pl.savefig("plots/DeltaGammas_METHOD1_Om="+str(shift)+"_nu=nu'.png", dpi = 150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()

shift=0
#--------------------------PLOT CONVERGENCE WITH RESPECT TO INVERSION RANGE ------------------------------------------

#INVROVR=[1,2,4]
INVROVR=[5,10,15,20,25,30]

#---  Helper functions (include translation from  nambu to physical CAUTION : USING TIME REVERSAL SYMM) 
def plotgamma_conv_axis( use_pl, arr1, arr2, arr3, arr4, arr5, arr6,string ):
    zarr = [ arr1[shift + (bdim1-1)/2,fdimo21,fdimo21,0,0,0,0,0,0,0], arr2[shift + (bdim12-1)/2,fdimo22,fdimo22,0,0,0,0,0,0,0],arr3[shift + (bdim3-1)/2,fdimo23,fdimo23,0,0,0,0,0,0,0],arr4[shift + (bdim4-1)/2,fdimo24,fdimo24,0,0,0,0,0,0,0], arr5[shift + (bdim5-1)/2,fdimo25,fdimo25,0,0,0,0,0,0,0], arr6[shift + (bdim6-1)/2,fdimo26,fdimo26,0,0,0,0,0,0,0]]
    pl.plot( INVROVR, zarr, marker = 'o', linestyle='-', color='b', ms=3, mew=0.2)
    use_pl.set_title( string , fontsize=10 )

    return


#--- Plot along axis \nu'=0

#pl.suptitle(r"$U=$" + str(UINT) + r"     $\beta=$" + str(BETA) + r"     $\epsilon=$" + str(EPS) +  r"     $\Omega=$" + str(shift) + r"$*2\pi/\beta$" + r"     Notation: $\phi(\Omega,\omega,\omega')$")

plotgamma_conv_axis( pl.subplot(2,2,1), regamma_pp1 + regamma_pp1[:,:,::-1,:,:,:,:,:,:,:], regamma_pp2 + regamma_pp2[:,:,::-1,:,:,:,:,:,:,:], regamma_pp3 + regamma_pp3[:,:,::-1,:,:,:,:,:,:,:], regamma_pp4 + regamma_pp4[:,:,::-1,:,:,:,:,:,:,:],regamma_pp5 + regamma_pp5[:,:,::-1,:,:,:,:,:,:,:],regamma_pp6 + regamma_pp6[:,:,::-1,:,:,:,:,:,:,:], RE + r"\Gamma^{0, 0, \omega =0}_{s}$" ) # flip sign of w_out
pl.xlabel(r"$ASYR/INVR$")
plotgamma_conv_axis( pl.subplot(2,2,2), regamma_pp1 - regamma_pp1[:,:,::-1,:,:,:,:,:,:,:], regamma_pp2 - regamma_pp2[:,:,::-1,:,:,:,:,:,:,:], regamma_pp3 - regamma_pp3[:,:,::-1,:,:,:,:,:,:,:], regamma_pp4 - regamma_pp4[:,:,::-1,:,:,:,:,:,:,:], regamma_pp5 - regamma_pp5[:,:,::-1,:,:,:,:,:,:,:],regamma_pp6 - regamma_pp6[:,:,::-1,:,:,:,:,:,:,:], RE + r"\Gamma^{\nu, \nu'=0, \omega =0}_{t}$" ) # flip sign of w_out
pl.xlabel(r"$ASYR/INVR$")
plotgamma_conv_axis( pl.subplot(2,2,3), 2*regamma_ph1 - regamma_xph1, 2*regamma_ph2 - regamma_xph2, 2*regamma_ph3 - regamma_xph3, 2*regamma_ph4 - regamma_xph4,2*regamma_ph5 - regamma_xph5,2*regamma_ph6 - regamma_xph6, RE + r"\Gamma^{\nu, \nu'=0, \omega =0}_{d}$" )
pl.xlabel(r"$ASYR/INVR$")
plotgamma_conv_axis( pl.subplot(2,2,4), - regamma_xph1, - regamma_xph2, - regamma_xph3, - regamma_xph4,- regamma_xph5,- regamma_xph6, RE + r"\Gamma^{\nu, \nu'=0, \omega =0}_{m}$" )
pl.xlabel(r"$ASYR/INVR$")

pl.tight_layout()

#--- Save to file
pl.savefig("plots/Gammas_METHOD1_conv.png", dpi = 150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()

shift=20

plotgamma_conv_axis( pl.subplot(2,2,1), regamma_pp1 + regamma_pp1[:,:,::-1,:,:,:,:,:,:,:], regamma_pp2 + regamma_pp2[:,:,::-1,:,:,:,:,:,:,:], regamma_pp3 + regamma_pp3[:,:,::-1,:,:,:,:,:,:,:], regamma_pp4 + regamma_pp4[:,:,::-1,:,:,:,:,:,:,:],regamma_pp5 + regamma_pp5[:,:,::-1,:,:,:,:,:,:,:],regamma_pp6 + regamma_pp6[:,:,::-1,:,:,:,:,:,:,:], RE + r"\Gamma^{0, 0, \omega =20}_{s}$" ) # flip sign of w_out
pl.xlabel(r"$ASYR/INVR$")
plotgamma_conv_axis( pl.subplot(2,2,2), regamma_pp1 - regamma_pp1[:,:,::-1,:,:,:,:,:,:,:], regamma_pp2 - regamma_pp2[:,:,::-1,:,:,:,:,:,:,:], regamma_pp3 - regamma_pp3[:,:,::-1,:,:,:,:,:,:,:], regamma_pp4 - regamma_pp4[:,:,::-1,:,:,:,:,:,:,:], regamma_pp5 - regamma_pp5[:,:,::-1,:,:,:,:,:,:,:],regamma_pp6 - regamma_pp6[:,:,::-1,:,:,:,:,:,:,:], RE + r"\Gamma^{\nu, \nu'=0, \omega =20}_{t}$" ) # flip sign of w_out
pl.xlabel(r"$ASYR/INVR$")
plotgamma_conv_axis( pl.subplot(2,2,3), 2*regamma_ph1 - regamma_xph1, 2*regamma_ph2 - regamma_xph2, 2*regamma_ph3 - regamma_xph3, 2*regamma_ph4 - regamma_xph4,2*regamma_ph5 - regamma_xph5,2*regamma_ph6 - regamma_xph6, RE + r"\Gamma^{\nu, \nu'=0, \omega =20}_{d}$" )
pl.xlabel(r"$ASYR/INVR$")
plotgamma_conv_axis( pl.subplot(2,2,4), - regamma_xph1, - regamma_xph2, - regamma_xph3, - regamma_xph4,- regamma_xph5,- regamma_xph6, RE + r"\Gamma^{\nu, \nu'=0, \omega =20}_{m}$" )
pl.xlabel(r"$ASYR/INVR$")

pl.tight_layout()

#--- Save to file
pl.savefig("plots/Gammas_METHOD1_Om="+str(shift)+"_conv.png", dpi = 150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()

