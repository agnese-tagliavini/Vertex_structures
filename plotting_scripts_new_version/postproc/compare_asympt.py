#!/usr/bin/python

#==========================================================================================================================
#
#               Plotting script for comparing the asymptotic functions calculated in the following ways:
#               
#               1) Full self-consistency
#               2) One shot calculation without fitting of the tails of the frequency summation
#               3) One shot calculation with fitting of the tails of the frequency summation
#
#               CONTENTS:
#               
#               - P func calculated with the three methods for different OMEGA transfer
#               - P differences with respect to the self-consistent result
#               - chi func calculated with the three methods
#               - chi differences with respect to the self-consistent result
#
#=========================================================================================================================
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

fname1 = "../../dat/H5FILES/BETA50/4SITES/U1/POSTPROC/dat_U1_Beta50_PFCB45_PARQ_SU2_METH2_INVR1xCOUNT_CORR_ED_ONESHOT.h5" 
fname2 = "../../dat/H5FILES/BETA50/4SITES/U1/POSTPROC/dat_U1_Beta50_PFCB1000_PARQ_SU2_METH2_INVR1xCOUNT_CORR_ED.h5" 
fname3 = "../../dat/H5FILES/BETA50/4SITES/U1/POSTPROC/dat_U1_Beta50_PFCB90_PARQ_SU2_METH2_INVR1xCOUNT_CORR_SELFCON.h5"


#------ DECIDE WHICH COMPARISON YOU WANT TO ANALYZE


if len(sys.argv) > 1:
    fname1 = str(sys.argv[1])
if len(sys.argv) > 1:
    fname2 = str(sys.argv[1])
if len(sys.argv) > 1:
    fname3 = str(sys.argv[1])

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

pl.rc('xtick', labelsize=10) 
pl.rc('ytick', labelsize=10) 
pl.rc('text', usetex=True)
pl.rc('text.latex', preamble='\usepackage{amsmath}')

RE = r"$\operatorname{Re}"
IM = r"$\operatorname{Im}"


def cm2inch(*tupl):
    inch = 2.54
    if isinstance(tupl[0], tuple):
        return tuple(i/inch for i in tupl[0])
    else:
        return tuple(i/inch for i in tupl)



#-------------------------------------P PLOTTING ------------------------------------------

print("Plotting P ...")

#--- Read P method 2 (INVR1xVERTR)
reP_pp1 = vert_mul * np.array(f1["/P_func/RE_PP"])
imP_pp1 = vert_mul * np.array(f1["/P_func/IM_PP"])
reP_ph1 = vert_mul * np.array(f1["/P_func/RE_PH"])
imP_ph1 = vert_mul * np.array(f1["/P_func/IM_PH"])
reP_xph1 = vert_mul * np.array(f1["/P_func/RE_XPH"])
imP_xph1 = vert_mul * np.array(f1["/P_func/IM_XPH"])

#--- Read P method 2 (INVR2xVERTR)
reP_pp2 = vert_mul * np.array(f2["/P_func/RE_PP"])
imP_pp2 = vert_mul * np.array(f2["/P_func/IM_PP"])
reP_ph2 = vert_mul * np.array(f2["/P_func/RE_PH"])
imP_ph2 = vert_mul * np.array(f2["/P_func/IM_PH"])
reP_xph2 = vert_mul * np.array(f2["/P_func/RE_XPH"])
imP_xph2 = vert_mul * np.array(f2["/P_func/IM_XPH"])

#--- Read P method 2 (INVR4xVERTR)
reP_pp3 = vert_mul * np.array(f3["/P_func/RE_PP"])
imP_pp3 = vert_mul * np.array(f3["/P_func/IM_PP"])
reP_ph3 = vert_mul * np.array(f3["/P_func/RE_PH"])
imP_ph3 = vert_mul * np.array(f3["/P_func/IM_PH"])
reP_xph3 = vert_mul * np.array(f3["/P_func/RE_XPH"])
imP_xph3 = vert_mul * np.array(f3["/P_func/IM_XPH"])


#--- RANGE FERMIONIC/BOSONIC FREQUENCIES

bdim1 = reP_pp1.shape[0]
fdim1 = reP_pp1.shape[1]
fdimo21 = fdim1/2

bdim2 = reP_pp2.shape[0]
fdim2 = reP_pp2.shape[1]
fdimo22 = fdim2/2

bdim3 = reP_pp3.shape[0]
fdim3 = reP_pp3.shape[1]
fdimo23 = fdim3/2
#phigrid = np.array(f["/P_func/fgrid"])

Pgrid_1 = np.array([(1*n+1)*pi/BETA for n in range(-fdimo21,fdimo21)])
Pgrid_2 = np.array([(1*n+1)*pi/BETA for n in range(-fdimo22,fdimo22)])
Pgrid_3 = np.array([(1*n+1)*pi/BETA for n in range(-fdimo23,fdimo23)])

Pgrid_plot = np.array([(2*n+1)*pi/BETA for n in range(-min(fdimo21,fdimo22,fdimo23),min(fdimo21,fdimo22,fdimo23))])

#---  Helper functions (include translation from  nambu to physical CAUTION : USING TIME REVERSAL SYMM) 
def plotP_axis( use_pl, arr1, arr2, arr3, string, legend ):
    zarr1 = np.array([ arr1[shift + (bdim1-1)/2,n,0,0,0,0,0,0] for n in range(fdim1)])
    zarr2 = np.array([ arr2[shift + (bdim2-1)/2,n,0,0,0,0,0,0] for n in range(fdim2)])
    zarr3 = np.array([ arr3[shift + (bdim3-1)/2,n,0,0,0,0,0,0] for n in range(fdim3)])
    pl.plot( Pgrid_1, zarr1,linestyle='None', marker = 'o',  color='r', ms=3, mew=0.3, label=r"Oneshot")
    pl.plot( Pgrid_2, zarr2,linestyle='None', marker = 'd',  color='b', ms=3, mew=0.3, label=r"ED")
    pl.plot( Pgrid_3, zarr3,linestyle='None', marker = 'x',  color='g', ms=3, mew=0.3, label=r"Selfcons")
    pl.xlim([min(Pgrid_plot),max(Pgrid_plot)])
    pl.xticks([ ])
    use_pl.set_title( string , fontsize=12 )
    if(legend):
        pl.legend(prop={'size':7})

    return


#--- Plot along \nu'=0 

#pl.suptitle(r"$U=$" + str(UINT) + r"     $\beta=$" + str(BETA) + r"     $\epsilon=$" + str(EPS) +  r"     $\Omega=$" + str(shift) + r"$*2\pi/\beta$")

fig = pl.figure(figsize=cm2inch(18.0,6.0))
plotP_axis( pl.subplot(1,3,1), reP_pp1 + reP_pp1, reP_pp2 + reP_pp2, reP_pp3 + reP_pp3, RE + r"\lambda^{\nu, \omega =0}_{s}$", False) # flip sign of w_out
#plotP_axis( pl.subplot1331,2), reP_pp1 - reP_pp1, reP_pp2 - reP_pp2, reP_pp3 - reP_pp3, RE + r"\lambda^{\nu, \omega =0}_{t}$", True) # flip sign of w_out
plotP_axis( pl.subplot(1,3,2), 2*reP_ph1 - reP_xph1, 2*reP_ph2 - reP_xph2, 2*reP_ph3 - reP_xph3, RE + r"\lambda^{\nu, \omega =0}_{d}$",False )
#pl.xlabel(r"$\nu$", fontsize=10)
plotP_axis( pl.subplot(1,3,3), - reP_xph1, - reP_xph2, - reP_xph3, RE + r"\lambda^{\nu,\omega =0}_{m}$",True )
#pl.xlabel(r"$\nu$", fontsize=10)

pl.tight_layout()

#--- Save to file
pl.savefig("plots/P_comparison_Om="+str(shift)+".png", dpi = 200)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()

#shift=1
#
##pl.suptitle(r"$U=$" + str(UINT) + r"     $\beta=$" + str(BETA) + r"     $\epsilon=$" + str(EPS) +  r"     $\Omega=$" + str(shift) + r"$*2\pi/\beta$")
#
#fig = pl.figure(figsize=cm2inch(12.0,12.0))
#plotP_axis( pl.subplot(2,2,1), reP_pp1 + reP_pp1, reP_pp2 + reP_pp2, reP_pp3 + reP_pp3, RE + r"\lambda^{\nu, \omega =0}_{s}$", False) # flip sign of w_out
#plotP_axis( pl.subplot(2,2,2), reP_pp1 - reP_pp1, reP_pp2 - reP_pp2, reP_pp3 - reP_pp3, RE + r"\lambda^{\nu, \omega =0}_{t}$", True) # flip sign of w_out
#plotP_axis( pl.subplot(2,2,3), 2*reP_ph1 - reP_xph1, 2*reP_ph2 - reP_xph2, 2*reP_ph3 - reP_xph3, RE + r"\lambda^{\nu, \omega =0}_{d}$",False )
#pl.xlabel(r"$\nu$", fontsize=10)
#plotP_axis( pl.subplot(2,2,4), - reP_xph1, - reP_xph2, - reP_xph3, RE + r"\lambda^{\nu,\omega =0}_{m}$",False )
#pl.xlabel(r"$\nu$", fontsize=10)
#
#pl.tight_layout()
#
##--- Save to file
#pl.savefig("plots/P_comparison_Om="+str(shift)+".png", dpi = 200)
#pl.figure(dpi=100) # Reset dpi to default
#pl.clf()
#
#shift=2
#
##pl.suptitle(r"$U=$" + str(UINT) + r"     $\beta=$" + str(BETA) + r"     $\epsilon=$" + str(EPS) +  r"     $\Omega=$" + str(shift) + r"$*2\pi/\beta$")
#
#fig = pl.figure(figsize=cm2inch(12.0,12.0))
#plotP_axis( pl.subplot(2,2,1), reP_pp1 + reP_pp1, reP_pp2 + reP_pp2, reP_pp3 + reP_pp3, RE + r"\lambda^{\nu, \omega =0}_{s}$", False) # flip sign of w_out
#plotP_axis( pl.subplot(2,2,2), reP_pp1 - reP_pp1, reP_pp2 - reP_pp2, reP_pp3 - reP_pp3, RE + r"\lambda^{\nu, \omega =0}_{t}$", True) # flip sign of w_out
#plotP_axis( pl.subplot(2,2,3), 2*reP_ph1 - reP_xph1, 2*reP_ph2 - reP_xph2, 2*reP_ph3 - reP_xph3, RE + r"\lambda^{\nu, \omega =0}_{d}$",False )
#pl.xlabel(r"$\nu$", fontsize=10)
#plotP_axis( pl.subplot(2,2,4), - reP_xph1, - reP_xph2, - reP_xph3, RE + r"\lambda^{\nu,\omega =0}_{m}$",False )
#pl.xlabel(r"$\nu$", fontsize=10)
#
#pl.tight_layout()
#
##--- Save to file
#pl.savefig("plots/P_comparison_Om="+str(shift)+".png", dpi = 200)
#pl.figure(dpi=100) # Reset dpi to default
#pl.clf()
#
#shift=4
#
##pl.suptitle(r"$U=$" + str(UINT) + r"     $\beta=$" + str(BETA) + r"     $\epsilon=$" + str(EPS) +  r"     $\Omega=$" + str(shift) + r"$*2\pi/\beta$")
#
#fig = pl.figure(figsize=cm2inch(12.0,12.0))
#plotP_axis( pl.subplot(2,2,1), reP_pp1 + reP_pp1, reP_pp2 + reP_pp2, reP_pp3 + reP_pp3, RE + r"\lambda^{\nu, \omega =0}_{s}$", False) # flip sign of w_out
#plotP_axis( pl.subplot(2,2,2), reP_pp1 - reP_pp1, reP_pp2 - reP_pp2, reP_pp3 - reP_pp3, RE + r"\lambda^{\nu, \omega =0}_{t}$", True) # flip sign of w_out
#plotP_axis( pl.subplot(2,2,3), 2*reP_ph1 - reP_xph1, 2*reP_ph2 - reP_xph2, 2*reP_ph3 - reP_xph3, RE + r"\lambda^{\nu, \omega =0}_{d}$",False )
#pl.xlabel(r"$\nu$", fontsize=10)
#plotP_axis( pl.subplot(2,2,4), - reP_xph1, - reP_xph2, - reP_xph3, RE + r"\lambda^{\nu,\omega =0}_{m}$",False )
#pl.xlabel(r"$\nu$", fontsize=10)
#
#pl.tight_layout()
#
##--- Save to file
#pl.savefig("plots/P_comparison_Om="+str(shift)+".png", dpi = 200)
#pl.figure(dpi=100) # Reset dpi to default
#pl.clf()
#
#shift=8
#
#pl.suptitle(r"$U=$" + str(UINT) + r"     $\beta=$" + str(BETA) + r"     $\epsilon=$" + str(EPS) +  r"     $\Omega=$" + str(shift) + r"$*2\pi/\beta$")
#
#plotP_axis( pl.subplot(2,2,1), reP_pp1 + reP_pp1, reP_pp2 + reP_pp2, reP_pp3 + reP_pp3, RE + r"\lambda^{\nu, \omega =8}_{s}$" , False) # flip sign of w_out
#pl.xlabel(r"$\nu$")
#plotP_axis( pl.subplot(2,2,2), reP_pp1 - reP_pp1, reP_pp2 - reP_pp2, reP_pp3 - reP_pp3, RE + r"\lambda^{\nu, \omega =8}_{t}$" ,True) # flip sign of w_out
#pl.xlabel(r"$\nu$")
#plotP_axis( pl.subplot(2,2,3), 2*reP_ph1 - reP_xph1, 2*reP_ph2 - reP_xph2, 2*reP_ph3 - reP_xph3, RE + r"\lambda^{\nu, \omega =8}_{d}$" , False)
#pl.xlabel(r"$\nu$")
#plotP_axis( pl.subplot(2,2,4), - reP_xph1, - reP_xph2, - reP_xph3, RE + r"\lambda^{\nu,\omega =8}_{m}$" , False)
#pl.xlabel(r"$\nu$")
#
#pl.tight_layout()
#
##--- Save to file
#pl.savefig("plots/P_comparison_Om="+str(shift)+".png", dpi = 150)
#pl.figure(dpi=100) # Reset dpi to default
#pl.clf()

shift=0
#-------------------------------------P PLOTTING DIFFERENCES WITH RESPECT TO THE FULL SELF-CONSISTENT CALCULATION ------------------------------------------

print("Plotting delta P ...")

#--- Helper functions for the subtractions of functions defined in different ranges


#-----DIFF 1-2
def diff_reP_pp1(W,w):
    return reP_pp1[W + (bdim1-1)/2,w+fdimo21,0,0,0,0,0,0]-reP_pp2[W + (bdim2-1)/2,w+fdimo22,0,0,0,0,0,0]

def diff_reP_ph1(W,w):
    return reP_ph1[W + (bdim1-1)/2,w+fdimo21,0,0,0,0,0,0]-reP_ph2[W + (bdim2-1)/2,w+fdimo22,0,0,0,0,0,0]

def diff_reP_xph1(W,w):
    return reP_xph1[W + (bdim1-1)/2,w+fdimo21,0,0,0,0,0,0]-reP_xph2[W + (bdim2-1)/2,w+fdimo22,0,0,0,0,0,0]


#---------0
def diff_reP_pp2(W,w):
    return reP_pp2[W + (bdim2-1)/2,w+fdimo22,0,0,0,0,0,0]-reP_pp2[W + (bdim2-1)/2,w+fdimo22,0,0,0,0,0,0]

def diff_reP_ph2(W,w):
    return reP_ph2[W + (bdim2-1)/2,w+fdimo22,0,0,0,0,0,0]-reP_ph2[W + (bdim2-1)/2,w+fdimo22,0,0,0,0,0,0]

def diff_reP_xph2(W,w):
    return reP_xph2[W + (bdim2-1)/2,w+fdimo22,0,0,0,0,0,0]-reP_xph2[W + (bdim2-1)/2,w+fdimo22,0,0,0,0,0,0]

#-------DIFF 3-2
def diff_reP_pp3(W,w):
    return reP_pp3[W + (bdim3-1)/2,w+fdimo23,0,0,0,0,0,0]-reP_pp2[W + (bdim2-1)/2,w+fdimo22,0,0,0,0,0,0]

def diff_reP_ph3(W,w):
    return reP_ph3[W + (bdim3-1)/2,w+fdimo23,0,0,0,0,0,0]-reP_ph2[W + (bdim2-1)/2,w+fdimo22,0,0,0,0,0,0]

def diff_reP_xph3(W,w):
    return reP_xph3[W + (bdim3-1)/2,w+fdimo23,0,0,0,0,0,0]-reP_xph2[W + (bdim2-1)/2,w+fdimo22,0,0,0,0,0,0]


#---  Helper functions (include translation from  nambu to physical CAUTION : USING TIME REVERSAL SYMM) 

def plotP_axis( use_pl, arr1, arr2, arr3, string,legend):
    pl.plot( Pgrid_plot, arr1,linestyle = 'None', marker = 'o', color='r', ms=3, mew=0.5, label=r"Oneshot")
    pl.plot( Pgrid_plot, arr2,linestyle = 'None', marker = 'd', color='b', ms=3, mew=0.5, label=r"ED")
    pl.plot( Pgrid_plot, arr3,linestyle = 'None', marker = 'x', color='g', ms=3, mew=0.5, label=r"Selfcons")
    pl.xlim([min(Pgrid_plot),max(Pgrid_plot)])
    use_pl.set_title( string , fontsize=12 )
    if(legend):
        pl.legend(prop={'size':7})

    return

#--- Plot along axis \nu'=0

#pl.suptitle(r"$U=$" + str(UINT) + r"     $\beta=$" + str(BETA) + r"     $\epsilon=$" + str(EPS) +  r"     $\Omega=$" + str(shift) + r"$*2\pi/\beta$" + r"     Notation: $\phi(\Omega,\omega,\omega')$")

fig = pl.figure(figsize=cm2inch(18.0,6.0))
plotP_axis( pl.subplot(1,3,1), np.array([diff_reP_pp1(shift,n) + diff_reP_pp1(shift,n) for n in range(-min(fdimo21,fdimo22,fdimo23),min(fdimo21,fdimo22,fdimo23))]), np.array([diff_reP_pp2(shift,n) + diff_reP_pp2(shift,n) for n in range(-min(fdimo21,fdimo22,fdimo23),min(fdimo21,fdimo22,fdimo23))]), np.array([diff_reP_pp3(shift,n) + diff_reP_pp3(shift,n) for n in range(-min(fdimo21,fdimo22,fdimo23),min(fdimo21,fdimo22,fdimo23))]), RE + r"\delta \lambda^{\nu, \omega =0}_{s}$", False ) # flip sign of w_out
pl.xlabel(r"$\nu$", fontsize=10)
#plotP_axis( pl.subplot(1,3,2), np.array([diff_reP_pp1(shift,n) - diff_reP_pp1(shift,n) for n in range(-min(fdimo21,fdimo22,fdimo23),min(fdimo21,fdimo22,fdimo23))]), np.array([diff_reP_pp2(shift,n) - diff_reP_pp2(shift,n) for n in range(-min(fdimo21,fdimo22,fdimo23),min(fdimo21,fdimo22,fdimo23))]), np.array([diff_reP_pp3(shift,n) - diff_reP_pp3(shift,n) for n in range(-min(fdimo21,fdimo22,fdimo23),min(fdimo21,fdimo22,fdimo23))]), RE + r"\delta \lambda^{\nu, \omega =0}_{t}$", True) # flip sign of w_out

plotP_axis( pl.subplot(1,3,2), np.array([2*diff_reP_ph1(shift,n) - diff_reP_xph1(shift,n) for n in range(-min(fdimo21,fdimo22,fdimo23),min(fdimo21,fdimo22,fdimo23))]), np.array([2*diff_reP_ph2(shift,n) - diff_reP_xph2(shift,n) for n in range(-min(fdimo21,fdimo22,fdimo23),min(fdimo21,fdimo22,fdimo23))]), np.array([2*diff_reP_ph3(shift,n) - diff_reP_xph3(shift,n) for n in range(-min(fdimo21,fdimo22,fdimo23),min(fdimo21,fdimo22,fdimo23))]), RE + r"\delta \lambda^{\nu, \omega =0}_{d}$" , False) # flip sign of w_out

pl.xlabel(r"$\nu$", fontsize=10)
plotP_axis( pl.subplot(1,3,3), np.array([- diff_reP_xph1(shift,n) for n in range(-min(fdimo21,fdimo22,fdimo23),min(fdimo21,fdimo22,fdimo23))]), np.array([- diff_reP_xph2(shift,n) for n in range(-min(fdimo21,fdimo22,fdimo23),min(fdimo21,fdimo22,fdimo23))]), np.array([ - diff_reP_xph3(shift,n) for n in range(-min(fdimo21,fdimo22,fdimo23),min(fdimo21,fdimo22,fdimo23))]), RE + r"\delta \lambda^{\nu, \omega =0}_{m}$" , False) # flip sign of w_out

pl.xlabel(r"$\nu$", fontsize=10)
pl.tight_layout()

#--- Save to file
pl.savefig("plots/DeltaP_comparison_Om="+str(shift)+".png", dpi = 200)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()


#-------------------------------------chi PLOTTING ------------------------------------------

print("Plotting chi ...")

#--- Read P method 2 (INVR1xVERTR)
rechi_pp1 = vert_mul * np.array(f1["/chi_func/RE_PP"])
imchi_pp1 = vert_mul * np.array(f1["/chi_func/IM_PP"])
rechi_ph1 = vert_mul * np.array(f1["/chi_func/RE_PH"])
imchi_ph1 = vert_mul * np.array(f1["/chi_func/IM_PH"])
rechi_xph1 = vert_mul * np.array(f1["/chi_func/RE_XPH"])
imchi_xph1 = vert_mul * np.array(f1["/chi_func/IM_XPH"])

#--- Read P method 2 (INVR2xVERTR)
rechi_pp2 = vert_mul * np.array(f2["/chi_func/RE_PP"])
imchi_pp2 = vert_mul * np.array(f2["/chi_func/IM_PP"])
rechi_ph2 = vert_mul * np.array(f2["/chi_func/RE_PH"])
imchi_ph2 = vert_mul * np.array(f2["/chi_func/IM_PH"])
rechi_xph2 = vert_mul * np.array(f2["/chi_func/RE_XPH"])
imchi_xph2 = vert_mul * np.array(f2["/chi_func/IM_XPH"])

#--- Read P method 2 (INVR4xVERTR)
rechi_pp3 = vert_mul * np.array(f3["/chi_func/RE_PP"])
imchi_pp3 = vert_mul * np.array(f3["/chi_func/IM_PP"])
rechi_ph3 = vert_mul * np.array(f3["/chi_func/RE_PH"])
imchi_ph3 = vert_mul * np.array(f3["/chi_func/IM_PH"])
rechi_xph3 = vert_mul * np.array(f3["/chi_func/RE_XPH"])
imchi_xph3 = vert_mul * np.array(f3["/chi_func/IM_XPH"])


#--- RANGE FERMIONIC/BOSONIC FREQUENCIES

bdim1 = rechi_pp1.shape[0]
bdimo21 = (bdim1-1)/2

bdim2 = rechi_pp2.shape[0]
bdimo22 = (bdim2-1)/2

bdim3 = rechi_pp3.shape[0]
bdimo23 = (bdim3-1)/2

#phigrid = np.array(f["/chi_func/fgrid"])

chigrid_1 = np.array([(2*n)*pi/BETA for n in range(-bdimo21,bdimo21+1)])
chigrid_2 = np.array([(2*n)*pi/BETA for n in range(-bdimo22,bdimo22+1)])
chigrid_3 = np.array([(2*n)*pi/BETA for n in range(-bdimo23,bdimo23+1)])

chigrid_plot = np.array([(2*n)*pi/BETA for n in range(-min(bdimo21,bdimo22,bdimo23),min(bdimo21,bdimo22,bdimo23)+1)])

#---  Helper functions (include translation from  nambu to physical CAUTION : USING TIME REVERSAL SYMM) 
def plotchi_axis( use_pl, arr1, arr2, arr3, string,legend ):
    zarr1 = np.array([ arr1[n ,0,0,0,0,0] for n in range(bdim1)])
    zarr2 = np.array([ arr2[n ,0,0,0,0,0] for n in range(bdim2)])
    zarr3 = np.array([ arr3[n ,0,0,0,0,0] for n in range(bdim3)])
    pl.plot( chigrid_1, zarr1, linestyle='None', marker = 'o',  color='r', ms=3, mew=0.5, label=r"Oneshot")
    pl.plot( chigrid_2, zarr2, linestyle='None', marker = 'd',  color='b', ms=3, mew=0.5, label=r"ED")
    pl.plot( chigrid_3, zarr3, linestyle='None', marker = 'x',  color='g', ms=3, mew=0.5, label=r"Selfcons")
    pl.xlim([min(chigrid_plot),max(chigrid_plot)])
    pl.xticks([ ])
    use_pl.set_title( string , fontsize=12 )
    if(legend):
        pl.legend(prop={'size':7})

    return


#--- Plot along \nu'=0 

#pl.suptitle(r"$U=$" + str(UINT) + r"     $\beta=$" + str(BETA) + r"     $\epsilon=$" + str(EPS) +  r"     $\Omega=$" + str(shift) + r"$*2\pi/\beta$")

fig = pl.figure(figsize=cm2inch(18.0,6.0))
plotchi_axis( pl.subplot(1,3,1), rechi_pp1 + rechi_pp1, rechi_pp2 + rechi_pp2, rechi_pp3 + rechi_pp3, RE + r"\chi^{ \omega}_{s}$", False ) # flip sign of w_out
#plotchi_axis( pl.subplot(1,3,2), rechi_pp1 - rechi_pp1, rechi_pp2 - rechi_pp2, rechi_pp3 - rechi_pp3, RE + r"\chi^{ \omega }_{t}$", True ) # flip sign of w_out
plotchi_axis( pl.subplot(1,3,2), 2*rechi_ph1 - rechi_xph1, 2*rechi_ph2 - rechi_xph2, 2*rechi_ph3 - rechi_xph3, RE + r"\chi^{ \omega }_{d}$" , False)
#pl.xlabel(r"$\omega$", fontsize=10)
plotchi_axis( pl.subplot(1,3,3), - rechi_xph1, - rechi_xph2, - rechi_xph3, RE + r"\chi^{\omega}_{m}$" , True)
#pl.xlabel(r"$\omega$", fontsize=10)

pl.tight_layout()

#--- Save to file
pl.savefig("plots/chi_comparison.png", dpi = 200)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()

#-------------------------------------chi PLOTTING DIFFERENCES WITH RESPECT TO THE FULL SELF-CONSISTENT CALCULATION ------------------------------------------

print("Plotting delta chi ...")

#--- Helper functions for the subtractions of functions defined in different ranges


#-----DIFF 1-2
def diff_rechi_pp1(W):
    return rechi_pp1[W + (bdim1-1)/2,0,0,0,0,0]-rechi_pp2[W + (bdim2-1)/2,0,0,0,0,0]

def diff_rechi_ph1(W):
    return rechi_ph1[W + (bdim1-1)/2,0,0,0,0,0]-rechi_ph2[W + (bdim2-1)/2,0,0,0,0,0]

def diff_rechi_xph1(W):
    return rechi_xph1[W + (bdim1-1)/2,0,0,0,0,0]-rechi_xph2[W + (bdim2-1)/2,0,0,0,0,0]


#---------0
def diff_rechi_pp2(W):
    return rechi_pp2[W + (bdim2-1)/2,0,0,0,0,0]-rechi_pp2[W + (bdim2-1)/2,0,0,0,0,0]

def diff_rechi_ph2(W):
    return rechi_ph2[W + (bdim2-1)/2,0,0,0,0,0]-rechi_ph2[W + (bdim2-1)/2,0,0,0,0,0]

def diff_rechi_xph2(W):
    return rechi_xph2[W + (bdim2-1)/2,0,0,0,0,0]-rechi_xph2[W + (bdim2-1)/2,0,0,0,0,0]

#-------DIFF 3-2
def diff_rechi_pp3(W):
    return rechi_pp3[W + (bdim3-1)/2,0,0,0,0,0]-rechi_pp2[W + (bdim2-1)/2,0,0,0,0,0]

def diff_rechi_ph3(W):
    return rechi_ph3[W + (bdim3-1)/2,0,0,0,0,0]-rechi_ph2[W + (bdim2-1)/2,0,0,0,0,0]

def diff_rechi_xph3(W):
    return rechi_xph3[W + (bdim3-1)/2,0,0,0,0,0]-rechi_xph2[W + (bdim2-1)/2,0,0,0,0,0]


#---  Helper functions (include translation from  nambu to physical CAUTION : USING TIME REVERSAL SYMM) 

def plotchi_axis( use_pl, arr1, arr2, arr3, string,legend ):
    pl.plot( chigrid_plot, arr1,linestyle='None', marker = 'o', color='r', ms=3, mew=0.5, label=r"Oneshot")
    pl.plot( chigrid_plot, arr2,linestyle='None', marker = 'd', color='b', ms=3, mew=0.5, label=r"ED")
    pl.plot( chigrid_plot, arr3,linestyle='None', marker = 'x', color='g', ms=3, mew=0.5, label=r"Selfcons")
    pl.xlim([min(chigrid_plot),max(chigrid_plot)])
    use_pl.set_title( string , fontsize=12 )
    if(legend):
        pl.legend(prop={'size':7})

    return

#--- Plot along axis \nu'=0

#pl.suptitle(r"$U=$" + str(UINT) + r"     $\beta=$" + str(BETA) + r"")

fig = pl.figure(figsize=cm2inch(18.0,6.0))
plotchi_axis( pl.subplot(1,3,1), np.array([diff_rechi_pp1(n) + diff_rechi_pp1(n) for n in range(-min(bdimo21,bdimo22,bdimo23),min(bdimo21,bdimo22,bdimo23)+1)]), np.array([diff_rechi_pp2(n) + diff_rechi_pp2(n) for n in range(-min(bdimo21,bdimo22,bdimo23),min(bdimo21,bdimo22,bdimo23)+1)]), np.array([diff_rechi_pp3(n) + diff_rechi_pp3(n) for n in range(-min(bdimo21,bdimo22,bdimo23),min(bdimo21,bdimo22,bdimo23)+1)]), RE + r"\delta \chi^{ \omega}_{s}$", False ) # flip sign of w_out
pl.xlabel(r"$\omega$", fontsize=10)
#plotchi_axis( pl.subplot(1,3,2), np.array([diff_rechi_pp1(n) - diff_rechi_pp1(n) for n in range(-min(bdimo21,bdimo22,bdimo23),min(bdimo21,bdimo22,bdimo23)+1)]), np.array([diff_rechi_pp2(n) - diff_rechi_pp2(n) for n in range(-min(bdimo21,bdimo22,bdimo23),min(bdimo21,bdimo22,bdimo23)+1)]), np.array([diff_rechi_pp3(n) - diff_rechi_pp3(n) for n in range(-min(bdimo21,bdimo22,bdimo23),min(bdimo21,bdimo22,bdimo23)+1)]), RE + r"\delta \chi^{ \omega}_{t}$", True) # flip sign of w_out
plotchi_axis( pl.subplot(1,3,2), np.array([2*diff_rechi_ph1(n) - diff_rechi_xph1(n) for n in range(-min(bdimo21,bdimo22,bdimo23),min(bdimo21,bdimo22,bdimo23)+1)]), np.array([2*diff_rechi_ph2(n) - diff_rechi_xph2(n) for n in range(-min(bdimo21,bdimo22,bdimo23),min(bdimo21,bdimo22,bdimo23)+1)]), np.array([2*diff_rechi_ph3(n) - diff_rechi_xph3(n) for n in range(-min(bdimo21,bdimo22,bdimo23),min(bdimo21,bdimo22,bdimo23)+1)]), RE + r"\delta \chi^{\omega}_{d}$" , False) # flip sign of w_out
pl.xlabel(r"$\omega$", fontsize=10)
plotchi_axis( pl.subplot(1,3,3), np.array([ - diff_rechi_xph1(n) for n in range(-min(bdimo21,bdimo22,bdimo23),min(bdimo21,bdimo22,bdimo23)+1)]), np.array([- diff_rechi_xph2(n) for n in range(-min(bdimo21,bdimo22,bdimo23),min(bdimo21,bdimo22,bdimo23)+1)]), np.array([ - diff_rechi_xph3(n) for n in range(-min(bdimo21,bdimo22,bdimo23),min(bdimo21,bdimo22,bdimo23)+1)]), RE + r"\delta \chi^{\omega}_{m}$" , False) # flip sign of w_out
pl.xlabel(r"$\omega$", fontsize=10)
pl.tight_layout()

#--- Save to file
pl.savefig("plots/Deltachi_comparison.png", dpi = 200)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()



