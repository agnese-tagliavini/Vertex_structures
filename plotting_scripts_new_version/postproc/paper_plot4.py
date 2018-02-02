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
fname1 = "../../dat/H5FILES/BETA50/4SITES/U1p75/POSTPROC/dat_U1p75_Beta50_PFCB1000_PARQ_SU2_METH2_INVR1xCOUNTm40_CORR_ED.h5"
fname2 = "../../dat/H5FILES/BETA50/4SITES/U1p75/POSTPROC/dat_U1p75_Beta50_PFCB1000_PARQ_SU2_METH1_INVR1xCOUNTm40_ASYR10xINVR_CORR_ED.h5"
fname3 = "../../dat/H5FILES/BETA50/4SITES/U1p75/POSTPROC/dat_U1p75_Beta50_PFCB1000_PARQ_SU2_METH2_INVR1xCOUNTm40_ED.h5"

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

shift= 0

#--------------------------------------GENERAL PLOT SETTINGS------------------------------------------

pl.rc('xtick', labelsize=9) 
pl.rc('ytick', labelsize=9) 
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

#=================================================================================
#
#                           READ F
#
#=================================================================================

#--- Read vert method 2 (INVR1xVERTR)
revert_pp0 = vert_mul * np.array(f0["/vert_func/RE_PP"])
imvert_pp0 = vert_mul * np.array(f0["/vert_func/IM_PP"])
revert_ph0 = vert_mul * np.array(f0["/vert_func/RE_PH"])
imvert_ph0 = vert_mul * np.array(f0["/vert_func/IM_PH"])
revert_xph0 = vert_mul * np.array(f0["/vert_func/RE_XPH"])
imvert_xph0 = vert_mul * np.array(f0["/vert_func/IM_XPH"])

#--- Read vert method 2 (INVR1xVERTR)
revert_pp1 = vert_mul * np.array(f1["/vert_func/RE_PP"])
imvert_pp1 = vert_mul * np.array(f1["/vert_func/IM_PP"])
revert_ph1 = vert_mul * np.array(f1["/vert_func/RE_PH"])
imvert_ph1 = vert_mul * np.array(f1["/vert_func/IM_PH"])
revert_xph1 = vert_mul * np.array(f1["/vert_func/RE_XPH"])
imvert_xph1 = vert_mul * np.array(f1["/vert_func/IM_XPH"])

#--- Read vert method 2 (INVR1xVERTR NOCORR)
revert_pp2 = vert_mul * np.array(f2["/vert_func/RE_PP"])
imvert_pp2 = vert_mul * np.array(f2["/vert_func/IM_PP"])
revert_ph2 = vert_mul * np.array(f2["/vert_func/RE_PH"])
imvert_ph2 = vert_mul * np.array(f2["/vert_func/IM_PH"])
revert_xph2 = vert_mul * np.array(f2["/vert_func/RE_XPH"])
imvert_xph2 = vert_mul * np.array(f2["/vert_func/IM_XPH"])

#--- Read vert method 2 (INVR1xVERTR-5)
revert_pp3 = vert_mul * np.array(f3["/vert_func/RE_PP"])
imvert_pp3 = vert_mul * np.array(f3["/vert_func/IM_PP"])
revert_ph3 = vert_mul * np.array(f3["/vert_func/RE_PH"])
imvert_ph3 = vert_mul * np.array(f3["/vert_func/IM_PH"])
revert_xph3 = vert_mul * np.array(f3["/vert_func/RE_XPH"])
imvert_xph3 = vert_mul * np.array(f3["/vert_func/IM_XPH"])


#=================================================================================
#
#                           READ Trileg
#
#=================================================================================

retrileg_pp = vert_mul * np.array(f1["/P_func/RE_PP"])
imtrileg_pp = vert_mul * np.array(f1["/P_func/IM_PP"])
retrileg_ph = vert_mul * np.array(f1["/P_func/RE_PH"])
imtrileg_ph = vert_mul * np.array(f1["/P_func/IM_PH"])
retrileg_xph = vert_mul * np.array(f1["/P_func/RE_XPH"])
imtrileg_xph = vert_mul * np.array(f1["/P_func/IM_XPH"])

#=================================================================================
#
#                           READ chi
#
#=================================================================================

rechi_pp = vert_mul * np.array(f1["/chi_func/RE_PP"])
imchi_pp = vert_mul * np.array(f1["/chi_func/IM_PP"])
rechi_ph = vert_mul * np.array(f1["/chi_func/RE_PH"])
imchi_ph = vert_mul * np.array(f1["/chi_func/IM_PH"])
rechi_xph = vert_mul * np.array(f1["/chi_func/RE_XPH"])
imchi_xph = vert_mul * np.array(f1["/chi_func/IM_XPH"])




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

bdim_tri= retrileg_pp.shape[0]
fdim_tri= retrileg_pp.shape[1]
fdim_trio2 = fdim_tri/2
print bdim_tri, fdim_tri

bdim_chi= rechi_pp.shape[0]
print bdim_chi

#=====================================================================================
#
#                       DEFINE ASYMPTOTICS VERTEX
#
#====================================================================================
#Helper fuinctions for asymptotics

def trileg_pp(wb,wf):
    if((wb in range(-(bdim_tri-1)/2, (bdim_tri-1)/2)) & (wf in range(-fdim_trio2, fdim_trio2))):
        return retrileg_pp[wb + (bdim_tri-1)/2, wf + fdim_tri/2, 0,0,0,0,0,0]
    else:
        return 0.0

def trileg_ph(wb,wf):
    if((wb in range(-(bdim_tri-1)/2, (bdim_tri-1)/2)) & (wf in range(-fdim_trio2, fdim_trio2))):
        return retrileg_ph[wb + (bdim_tri-1)/2, wf + fdim_tri/2, 0,0,0,0,0,0]
    else:
        return 0.0

def trileg_xph(wb,wf):
    if((wb in range(-(bdim_tri-1)/2, (bdim_tri-1)/2)) & (wf in range(-fdim_trio2, fdim_trio2))):
        return retrileg_xph[wb + (bdim_tri-1)/2, wf + fdim_tri/2, 0,0,0,0,0,0]
    else:
        return 0.0

def chi_pp(wb):
    if(wb in range(-(bdim_chi-1)/2, (bdim_chi-1)/2)):
        return rechi_pp[wb + (bdim_chi-1)/2,0,0,0,0,0]
    else:
        return 0.0

def chi_ph(wb):
    if(wb in range(-(bdim_chi-1)/2, (bdim_chi-1)/2)):
        return rechi_ph[wb + (bdim_chi-1)/2,0,0,0,0,0]
    else:
        return 0.0

def chi_xph(wb):
    if(wb in range(-(bdim_chi-1)/2, (bdim_chi-1)/2)):
        return rechi_xph[wb + (bdim_chi-1)/2,0,0,0,0,0]
    else:
        return 0.0

#===================================
#   F asy
#===================================
def vert_asy_pp(wb, wf, wf1):
    Wxph = wf1-wf
    Wph = -mymod_abs(wb) - wf1 -wf - 1
    return UINT + chi_pp(wb) + trileg_pp(wb, wf)+ trileg_pp(wb, wf1)+ chi_xph(Wxph) + trileg_xph(Wxph, wf +myceil_div2(wb) - myfloor_div2(Wxph)) + trileg_xph(Wxph, -wf1-1  +myfloor_div2(wb) + myfloor_div2(Wxph)) + chi_ph(Wph)+ trileg_ph(Wph, wf +myceil_div2(wb) + myfloor_div2(Wph)) +trileg_ph(Wph, wf1+ myceil_div2(wb) + myfloor_div2(Wph))  

def vert_asy_ph(wb, wf, wf1):
    Wxph = wf1-wf
    Wpp = mymod_abs(wb) + wf1 + wf + 1
    return UINT + chi_ph(wb) + trileg_ph(wb, wf) + trileg_ph(wb, wf1)+ chi_xph(Wxph) + trileg_xph(Wxph, wf - myfloor_div2(wb) + myfloor_div2(Wxph)) + trileg_xph(Wxph, wf1  +myceil_div2(wb) - myceil_div2(Wxph)) + chi_pp(Wpp)+ trileg_pp(Wpp, wf - myfloor_div2(wb) - myceil_div2(Wpp)) +trileg_pp(Wpp, wf1 - myfloor_div2(wb) - myceil_div2(Wpp))  

def vert_asy_xph(wb, wf, wf1):
    Wph = wf1-wf
    Wpp = mymod_abs(wb) + wf1 + wf + 1
    return UINT + chi_xph(wb) + trileg_xph(wb, wf) + trileg_xph(wb, wf1) + chi_ph(Wph) + trileg_ph(Wph, wf - myfloor_div2(wb) + myfloor_div2(Wph)) + trileg_ph(Wph, wf1  +myceil_div2(wb) - myceil_div2(Wph)) + chi_pp(Wpp)+ trileg_pp(Wpp, wf - myfloor_div2(wb) - myceil_div2(Wpp)) + trileg_pp(Wpp, -wf1-1 - myfloor_div2(wb) + myfloor_div2(Wpp))  

#==================================
#   GAMMA ASY   -> WITH K2
#==================================
def gamma_asy_pp(wb, wf, wf1):
    Wxph = wf1-wf
    Wph = -mymod_abs(wb) - wf1 -wf - 1
    return UINT + chi_xph(Wxph) + trileg_xph(Wxph, wf +myceil_div2(wb) - myfloor_div2(Wxph)) + trileg_xph(Wxph, -wf1-1  +myfloor_div2(wb) + myfloor_div2(Wxph)) + chi_ph(Wph)+ trileg_ph(Wph, wf +myceil_div2(wb) + myfloor_div2(Wph)) +trileg_ph(Wph, wf1+ myceil_div2(wb) + myfloor_div2(Wph))  

def gamma_asy_ph(wb, wf, wf1):
    Wxph = wf1-wf
    Wpp = mymod_abs(wb) + wf1 + wf + 1
    return UINT + trileg_ph(wb, wf1)+ chi_xph(Wxph) + trileg_xph(Wxph, wf - myfloor_div2(wb) + myfloor_div2(Wxph)) + trileg_xph(Wxph, wf1  +myceil_div2(wb) - myceil_div2(Wxph)) + chi_pp(Wpp)+ trileg_pp(Wpp, wf - myfloor_div2(wb) - myceil_div2(Wpp)) +trileg_pp(Wpp, wf1 - myfloor_div2(wb) - myceil_div2(Wpp))  

def gamma_asy_xph(wb, wf, wf1):
    Wph = wf1-wf
    Wpp = mymod_abs(wb) + wf1 + wf + 1
    return UINT + chi_ph(Wph) + trileg_ph(Wph, wf - myfloor_div2(wb) + myfloor_div2(Wph)) + trileg_ph(Wph, wf1  +myceil_div2(wb) - myceil_div2(Wph)) + chi_pp(Wpp)+ trileg_pp(Wpp, wf - myfloor_div2(wb) - myceil_div2(Wpp)) + trileg_pp(Wpp, -wf1-1 - myfloor_div2(wb) + myfloor_div2(Wpp))  

#==================================
#   GAMMA ASY   -> WITHOUT K2
#==================================

def gamma_asy_nok2_pp(wb, wf, wf1):
    Wxph = wf1-wf
    Wph = -mymod_abs(wb) - wf1 -wf - 1
    return UINT + chi_xph(Wxph) + chi_ph(Wph)  

def gamma_asy_nok2_ph(wb, wf, wf1):
    Wxph = wf1-wf
    Wpp = mymod_abs(wb) + wf1 + wf + 1
    return UINT + chi_xph(Wxph) + chi_pp(Wpp)  

def gamma_asy_nok2_xph(wb, wf, wf1):
    Wph = wf1-wf
    Wpp = mymod_abs(wb) + wf1 + wf + 1
    return UINT + chi_ph(Wph) + chi_pp(Wpp)  

#
fdim_min = min([fdim0,fdim1,fdim2,fdim3])-40
print "fdim_min:",fdim_min, fdimo20
gammagrid_plot = np.array([(2*i+1)*pi/BETA for i in range(0, fdimo20)])
gammagrid_plot2 = np.array([(2*i+1)*pi/BETA for i in range(-fdim_trio2, fdim_trio2)])
gammagrid_plotsmall = np.array([(2*i+1)*pi/BETA for i in range(0, fdim_min/2)])
gammagrid_plotsmall2 = np.array([(2*i+1)*pi/BETA for i in range(-fdim_min/2, fdim_min/2)])
extend_gammagrid_left = np.array([(2*i+1)*pi/BETA for i in range(-fdim_trio2, -fdim_min/2)])
extend_gammagrid_right = np.array([(2*i+1)*pi/BETA for i in range(fdim_min/2, fdim_trio2)])

#------------------------------------------- -PLOT DIRECT INV CORRECTIONS AND NO CORR ------------------------------------------------------------------------
shift= 0

def plotgamma_diag1( use_pl, arr0, arr1, arr2, arr3, string, legend, noxticks):
    zarr0 = np.array([ arr0[0, n+fdimo20,n+fdimo20,0,0,0,0,0,0,0] for n in range(0,fdimo20)])
    zarr1 = np.array([ arr1[shift + (bdim1-1)/2,n+fdimo21,n+fdimo21,0,0,0,0,0,0,0] for n in range(0,fdim_min/2)])
    zarr2 = np.array([ arr2[shift + (bdim2-1)/2,n+fdimo22,n+fdimo22,0,0,0,0,0,0,0] for n in range(0,fdim_min/2)])
    zarr3 = np.array([ arr3[shift + (bdim3-1)/2,n+fdimo23,n+fdimo23,0,0,0,0,0,0,0] for n in range(0,fdim_min/2)])
    pl.plot( gammagrid_plotsmall, zarr1, marker = 'o', linestyle='-',linewidth=1.2, color='b', ms=3, mew=0.5, label=r"M1")
    pl.plot( gammagrid_plotsmall, zarr2, marker = 's', linestyle='-',linewidth=1.2, color='g', ms=3, mew=0.5, label=r"M2")
    pl.plot( gammagrid_plotsmall, zarr3, marker = 'd', linestyle='-',linewidth=1.2, color='r', ms=3, mew=0.5, label=r"NC")
    pl.plot( gammagrid_plot, zarr0, marker = 'None', linestyle='-',linewidth=1.2, color='k', label=r"EX")
    pl.xlim([0.2*min(gammagrid_plot), 0.2*max(gammagrid_plot)])
    if(noxticks):
        pl.xticks([])
    use_pl.set_title( string , fontsize=12 )
    if(legend):
        pl.legend(prop={'size':7})

    return

fig = pl.figure(figsize=cm2inch(6.0,20.0))
ax = plotgamma_diag1( pl.subplot(4,1,1), regamma_pp0 + regamma_pp0[:,:,::-1,:,:,:,:,:,:,:],regamma_pp1 + regamma_pp1[:,:,::-1,:,:,:,:,:,:,:], regamma_pp2 + regamma_pp2[:,:,::-1,:,:,:,:,:,:,:], regamma_pp3 + regamma_pp3[:,:,::-1,:,:,:,:,:,:,:], RE + r"\Gamma^{\nu=\nu', \omega =0 \frac{\pi}{\beta}}_{s}$", True, True ) # flip sign of w_out
plotgamma_diag1( pl.subplot(4,1,2), regamma_pp0 - regamma_pp0[:,:,::-1,:,:,:,:,:,:,:],regamma_pp1 - regamma_pp1[:,:,::-1,:,:,:,:,:,:,:], regamma_pp2 - regamma_pp2[:,:,::-1,:,:,:,:,:,:,:], regamma_pp3 - regamma_pp3[:,:,::-1,:,:,:,:,:,:,:], RE + r"\Gamma^{\nu= \nu', \omega =0 \frac{\pi}{\beta}}_{t}$", False, True ) # flip sign of w_out
plotgamma_diag1( pl.subplot(4,1,3), 2*regamma_ph0 - regamma_xph0,2*regamma_ph1 - regamma_xph1, 2*regamma_ph2 - regamma_xph2, 2*regamma_ph3 - regamma_xph3, RE + r"\Gamma^{\nu=\nu', \omega =0 \frac{\pi}{\beta}}_{d}$", False, True )
#pl.xlabel(r"$\nu$", fontsize = 10)
plotgamma_diag1( pl.subplot(4,1,4), - regamma_xph0,- regamma_xph1, - regamma_xph2,- regamma_xph3, RE + r"\Gamma^{\nu=\nu', \omega =0 \frac{\pi}{\beta}}_{m}$", False, False )
pl.xlabel(r"$\nu$", fontsize = 10)

pl.tight_layout()

#--- Save to file
pl.savefig("plots/Comparison_EXACTINV_METHOD2_METHOD1_axis_nu=nu'_W0.png", dpi = 200)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()

shift= 0
def plotgamma_axis_s( use_pl,  arr1, string, legend ):
    zarr0 = np.array([ vert_asy_pp(shift,n,0)+vert_asy_pp(shift, -n-1, 0) for n in range(-fdim_trio2,fdim_trio2)])
    zarr1 = np.array([ arr1[0,n+fdimo20-fdim_min/2,fdimo20,0,0,0,0,0,0,0] for n in range(0,fdim_min)])
    zarr4 = np.array([ arr1[0,n+fdimo20,fdimo20,0,0,0,0,0,0,0] for n in range(-fdim_trio2, -fdim_min/2)])
    zarr5 = np.array([ arr1[0,n+fdimo20,fdimo20,0,0,0,0,0,0,0] for n in range(fdim_min/2, fdim_trio2)])
    pl.plot( gammagrid_plotsmall2, zarr1, marker = 'o', linestyle='-', linewidth=1.2, color='r', ms=3, mew=0.5, label=r"ED(small)")
    pl.plot( extend_gammagrid_left, zarr4, marker = 'o', linestyle='-', linewidth=1.2, color='#FA8072', ms=3, mew=0.5, markerfacecolor='None', label= "ED(large)")
    pl.plot( extend_gammagrid_right, zarr5, marker = 'o', linestyle='-', linewidth=1.2 ,color='#FA8072', ms=3, mew=0.5, markerfacecolor='None', label="_nolegend_")
    pl.plot( gammagrid_plot2, zarr0, marker = 'None', linestyle='-', linewidth=1.2, color='k', label=r"ASY")
    pl.xlim([-5.0, 5.0])
    use_pl.set_title( string , fontsize=12 )
    if(legend):
        pl.legend(prop={'size':6})

    return

def plotgamma_axis_t( use_pl, arr1, string, legend ):
    zarr0 = np.array([ vert_asy_pp(shift,n,0)-vert_asy_pp(shift, -n-1, 0) for n in range(-fdim_trio2,fdim_trio2)])
    zarr1 = np.array([ arr1[0,n+fdimo20-fdim_min/2,fdimo20,0,0,0,0,0,0,0] for n in range(0,fdim_min)])
    zarr4 = np.array([ arr1[0,n+fdimo20,fdimo20,0,0,0,0,0,0,0] for n in range(-fdim_trio2, -fdim_min/2)])
    zarr5 = np.array([ arr1[0,n+fdimo20,fdimo20,0,0,0,0,0,0,0] for n in range(fdim_min/2, fdim_trio2)])
    pl.plot( gammagrid_plotsmall2, zarr1, marker = 'o', linestyle='-', linewidth=1.2, color='r', ms=3, mew=0.5, label=r"ED(small)")
    pl.plot( extend_gammagrid_left, zarr4, marker = 'o', linestyle='-', linewidth=1.2, color='#FA8072', ms=3, mew=0.5, markerfacecolor='None', label= "ED(large)")
    pl.plot( extend_gammagrid_right, zarr5, marker = 'o', linestyle='-', linewidth=1.2 ,color='#FA8072', ms=3, mew=0.5, markerfacecolor='None', label="_nolegend_")
    pl.plot( gammagrid_plot2, zarr0, marker = 'None', linestyle='-', linewidth=1.2, color='k', label=r"ASY")
    pl.xlim([-5.0, 5.0])
    use_pl.set_title( string , fontsize=12 )
    if(legend):
        pl.legend(prop={'size':6})

    return
def plotgamma_axis_d( use_pl, arr1, string, legend ):
    zarr0 = np.array([ 2*vert_asy_ph(shift,n,0)-vert_asy_xph(shift, n, 0) for n in range(-fdim_trio2,fdim_trio2)])
    zarr1 = np.array([ arr1[0,n+fdimo20-fdim_min/2,fdimo20,0,0,0,0,0,0,0] for n in range(0,fdim_min)])
    zarr4 = np.array([ arr1[0,n+fdimo20,fdimo20,0,0,0,0,0,0,0] for n in range(-fdim_trio2, -fdim_min/2)])
    zarr5 = np.array([ arr1[0,n+fdimo20,fdimo20,0,0,0,0,0,0,0] for n in range(fdim_min/2, fdim_trio2)])
    pl.plot( gammagrid_plotsmall2, zarr1, marker = 'o', linestyle='-', linewidth=1.2, color='r', ms=3, mew=0.5, label=r"ED(small)")
    pl.plot( extend_gammagrid_left, zarr4, marker = 'o', linestyle='-', linewidth=1.2, color='#FA8072', ms=3, mew=0.5, markerfacecolor='None', label= "ED(large)")
    pl.plot( extend_gammagrid_right, zarr5, marker = 'o', linestyle='-', linewidth=1.2 ,color='#FA8072', ms=3, mew=0.5, markerfacecolor='None', label="_nolegend_")
    pl.plot( gammagrid_plot2, zarr0, marker = 'None', linestyle='-', linewidth=1.2, color='k', label=r"ASY")
    pl.xlim([-5.0, 5.0])
    use_pl.set_title( string , fontsize=12 )
    if(legend):
        pl.legend(prop={'size':6})

    return

def plotgamma_axis_m( use_pl, arr1, string, legend ):
    zarr0 = np.array([ -vert_asy_xph(shift, n, 0) for n in range(-fdim_trio2,fdim_trio2)])
    zarr1 = np.array([ arr1[0,n+fdimo20-fdim_min/2,fdimo20,0,0,0,0,0,0,0] for n in range(0,fdim_min)])
    zarr4 = np.array([ arr1[0,n+fdimo20,fdimo20,0,0,0,0,0,0,0] for n in range(-fdim_trio2, -fdim_min/2)])
    zarr5 = np.array([ arr1[0,n+fdimo20,fdimo20,0,0,0,0,0,0,0] for n in range(fdim_min/2, fdim_trio2)])
    pl.plot( gammagrid_plotsmall2, zarr1, marker = 'o', linestyle='-', linewidth=1.2, color='r', ms=3, mew=0.5, label=r"ED(small)")
    pl.plot( extend_gammagrid_left, zarr4, marker = 'o', linestyle='-', linewidth=1.2, color='#FA8072', ms=3, mew=0.5, markerfacecolor='None', label= "ED(large)")
    pl.plot( extend_gammagrid_right, zarr5, marker = 'o', linestyle='-', linewidth=1.2 ,color='#FA8072', ms=3, mew=0.5, markerfacecolor='None', label="_nolegend_")
    pl.plot( gammagrid_plot2, zarr0, marker = 'None', linestyle='-', linewidth=1.2, color='k', label=r"ASY")
    pl.xlim([-5.0, 5.0])
    use_pl.set_title( string , fontsize=12 )
    if(legend):
        pl.legend(prop={'size':6})

    return

fig = pl.figure(figsize=cm2inch(12.0,12.0))
plotgamma_axis_s( pl.subplot(2,2,1), revert_pp0 + revert_pp0[:,:,::-1,:,:,:,:,:,:,:], RE + r"F^{\nu,\nu'=\pi/\beta, \omega = 0 }_{s}$", False ) # flip sign of w_out
plotgamma_axis_t( pl.subplot(2,2,2), revert_pp0 - revert_pp0[:,:,::-1,:,:,:,:,:,:,:], RE + r"F^{\nu,\nu'=\pi/\beta, \omega = 0 }_{t}$", True ) # flip sign of w_out
plotgamma_axis_d( pl.subplot(2,2,3), 2*revert_ph0 - revert_xph0, RE + r"F^{\nu,\nu'=\pi/\beta, \omega =0 }_{d}$", False )
pl.xlabel(r"$\nu$", fontsize = 10)
plotgamma_axis_m( pl.subplot(2,2,4), -revert_xph0, RE + r"F^{\nu,\nu'=\pi/\beta, \omega =0}_{m}$", False )
pl.xlabel(r"$\nu$", fontsize = 10)

pl.tight_layout()

#--- Save to file
pl.savefig("plots/Comparison_F_EXACTINV_METHOD2_METHOD1_axis_nu'=0_W0.png", dpi = 200)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()


#======================================================================================================
#
#                               GAMMA
#
#======================================================================================================

shift= 0

def plotgamma_axis_s( use_pl,  arr1, arr2,arr3, string, legend ):
    zarr0 = np.array([ gamma_asy_pp(shift,n,n)+gamma_asy_pp(shift, n, -n-1) for n in range(-fdim_trio2,fdim_trio2)])
    zarr2 = np.array([ gamma_asy_nok2_pp(shift,n,n)+gamma_asy_nok2_pp(shift, n, -n-1) for n in range(-fdim_trio2,fdim_trio2)])
    zarr1 = np.array([ arr1[0,n+fdimo20-fdim_min/2,n+fdimo20-fdim_min/2,0,0,0,0,0,0,0] for n in range(0,fdim_min)])
    zarr4 = np.array([ arr1[0,n+fdimo20,n+fdimo20,0,0,0,0,0,0,0] for n in range(-fdim_trio2, -fdim_min/2)])
    zarr5 = np.array([ arr1[0,n+fdimo20,n+fdimo20,0,0,0,0,0,0,0] for n in range(fdim_min/2, fdim_trio2)])
    pl.plot( gammagrid_plotsmall2, zarr1, marker = 'o', linestyle='-', linewidth=1.2,color='r', ms=3, mew=0.5, label=r"EX(small)")
    pl.plot( extend_gammagrid_left, zarr4, marker = 'o', linestyle='-', linewidth=1.2, color='#FA8072', ms=3, mew=0.5, markerfacecolor='None', label= "EX(large)")
    pl.plot( extend_gammagrid_right, zarr5, marker = 'o', linestyle='-', linewidth=1.2 ,color='#FA8072', ms=3, mew=0.5, markerfacecolor='None', label="_nolegend_")
    pl.plot( gammagrid_plot2, zarr2, marker = 'None', linestyle='--', linewidth=1.2,color='k', label=r"ASY ")
    pl.plot( gammagrid_plot2, zarr0, marker = 'None', linestyle='-', linewidth=1.2,color='b', label=r"ASY+$\lambda$")
    pl.xlim([-7.0, 7.0])
    use_pl.set_title( string , fontsize=12 )
    if(legend):
        pl.legend(prop={'size':7})

    return

def plotgamma_axis_t( use_pl, arr1, arr2,arr3, string, legend ):
    zarr0 = np.array([ gamma_asy_pp(shift,n,n)-gamma_asy_pp(shift, n, -n-1) for n in range(-fdim_trio2,fdim_trio2)])
    zarr2 = np.array([ gamma_asy_nok2_pp(shift,n,n)-gamma_asy_nok2_pp(shift, n, -n-1) for n in range(-fdim_trio2,fdim_trio2)])
    zarr1 = np.array([ arr1[0,n+fdimo20-fdim_min/2,n+fdimo20-fdim_min/2,0,0,0,0,0,0,0] for n in range(0,fdim_min)])
    zarr4 = np.array([ arr1[0,n+fdimo20,n+fdimo20,0,0,0,0,0,0,0] for n in range(-fdim_trio2, -fdim_min/2)])
    zarr5 = np.array([ arr1[0,n+fdimo20,n+fdimo20,0,0,0,0,0,0,0] for n in range(fdim_min/2, fdim_trio2)])
    pl.plot( gammagrid_plotsmall2, zarr1, marker = 'o', linestyle='-', linewidth=1.2,color='r', ms=3, mew=0.5, label=r"EX(small)")
    pl.plot( extend_gammagrid_left, zarr4, marker = 'o', linestyle='-', linewidth=1.2, color='#FA8072', ms=3, mew=0.5, markerfacecolor='None', label= "EX(large)")
    pl.plot( extend_gammagrid_right, zarr5, marker = 'o', linestyle='-', linewidth=1.2 ,color='#FA8072', ms=3, mew=0.5, markerfacecolor='None', label="_nolegend_")
    pl.plot( gammagrid_plot2, zarr2, marker = 'None', linestyle='--', linewidth=1.2,color='k', label=r"ASY ")
    pl.plot( gammagrid_plot2, zarr0, marker = 'None', linestyle='-', linewidth=1.2,color='b', label=r"ASY+$\lambda$")
    pl.xlim([-7.0, 7.0])
    use_pl.set_title( string , fontsize=12 )
    if(legend):
        pl.legend(prop={'size':7})

    return
def plotgamma_axis_d( use_pl, arr1, arr2,arr3, string, legend ):
    zarr0 = np.array([ 2*gamma_asy_ph(shift,n,n)-gamma_asy_xph(shift, n, n) for n in range(-fdim_trio2,fdim_trio2)])
    zarr2 = np.array([ 2*gamma_asy_nok2_ph(shift,n,n)-gamma_asy_nok2_xph(shift, n, n) for n in range(-fdim_trio2,fdim_trio2)])
    zarr1 = np.array([ arr1[0,n+fdimo20-fdim_min/2,n+fdimo20-fdim_min/2,0,0,0,0,0,0,0] for n in range(0,fdim_min)])
    zarr4 = np.array([ arr1[0,n+fdimo20,n+fdimo20,0,0,0,0,0,0,0] for n in range(-fdim_trio2, -fdim_min/2)])
    zarr5 = np.array([ arr1[0,n+fdimo20,n+fdimo20,0,0,0,0,0,0,0] for n in range(fdim_min/2, fdim_trio2)])
    pl.plot( gammagrid_plotsmall2, zarr1, marker = 'o', linestyle='-',linewidth=1.2, color='r', ms=3, mew=0.5, label=r"EX(small)")
    pl.plot( extend_gammagrid_left, zarr4, marker = 'o', linestyle='-', linewidth=1.2, color='#FA8072', ms=3, mew=0.5, markerfacecolor='None', label= "EX(large)")
    pl.plot( extend_gammagrid_right, zarr5, marker = 'o', linestyle='-', linewidth=1.2 ,color='#FA8072', ms=3, mew=0.5, markerfacecolor='None', label="_nolegend_")
    pl.plot( gammagrid_plot2, zarr2, marker = 'None', linestyle='--', linewidth=1.2,color='k', label=r"ASY ")
    pl.plot( gammagrid_plot2, zarr0, marker = 'None', linestyle='-', linewidth=1.2,color='b', label=r"ASY+$\lambda$")
    pl.xlim([-7.0, 7.0])
    use_pl.set_title( string , fontsize=12 )
    if(legend):
        pl.legend(prop={'size':7})

    return

def plotgamma_axis_m( use_pl, arr1, arr2,arr3, string, legend ):
    zarr0 = np.array([ -gamma_asy_xph(shift, n, n) for n in range(-fdim_trio2,fdim_trio2)])
    zarr2 = np.array([ -gamma_asy_nok2_xph(shift, n, n) for n in range(-fdim_trio2,fdim_trio2)])
    zarr1 = np.array([ arr1[0,n+fdimo20-fdim_min/2,n+fdimo20-fdim_min/2,0,0,0,0,0,0,0] for n in range(0,fdim_min)])
    zarr4 = np.array([ arr1[0,n+fdimo20,n+fdimo20,0,0,0,0,0,0,0] for n in range(-fdim_trio2, -fdim_min/2)])
    zarr5 = np.array([ arr1[0,n+fdimo20,n+fdimo20,0,0,0,0,0,0,0] for n in range(fdim_min/2, fdim_trio2)])
    pl.plot( gammagrid_plotsmall2, zarr1, marker = 'o', linestyle='-', linewidth=1.2,color='r', ms=3, mew=0.5, label=r"EX(small)")
    pl.plot( extend_gammagrid_left, zarr4, marker = 'o', linestyle='-', linewidth=1.2, color='#FA8072', ms=3, mew=0.5, markerfacecolor='None', label= "ED(large)")
    pl.plot( extend_gammagrid_right, zarr5, marker = 'o', linestyle='-', linewidth=1.2 ,color='#FA8072', ms=3, mew=0.5, markerfacecolor='None', label="_nolegend_")
    pl.plot( gammagrid_plot2, zarr2, marker = 'None', linestyle='--', linewidth=1.2,color='k', label=r"ASY ")
    pl.plot( gammagrid_plot2, zarr0, marker = 'None', linestyle='-', linewidth=1.2,color='b', label=r"ASY+$\lambda$")
    pl.xlim([-7.0, 7.0])
    use_pl.set_title( string , fontsize=12 )
    if(legend):
        pl.legend(prop={'size':7})

    return

fig = pl.figure(figsize=cm2inch(12.0,12.0)) 
plotgamma_axis_s( pl.subplot(2,2,1), regamma_pp0 + regamma_pp0[:,:,::-1,:,:,:,:,:,:,:], regamma_pp2 + regamma_pp2[:,:,::-1,:,:,:,:,:,:,:], regamma_pp3 + regamma_pp3[:,:,::-1,:,:,:,:,:,:,:], RE + r"\Gamma^{\nu=\nu', \omega=0 }_{s}$", False ) # flip sign of w_out
plotgamma_axis_t( pl.subplot(2,2,2), regamma_pp0 - regamma_pp0[:,:,::-1,:,:,:,:,:,:,:], regamma_pp2 - regamma_pp2[:,:,::-1,:,:,:,:,:,:,:], regamma_pp3 - regamma_pp3[:,:,::-1,:,:,:,:,:,:,:], RE + r"\Gamma^{\nu=\nu', \omega =0 }_{t}$", True ) # flip sign of w_out
plotgamma_axis_d( pl.subplot(2,2,3),2*regamma_ph0 - regamma_xph0, 2*regamma_ph2 - regamma_xph2, 2*regamma_ph3 - regamma_xph3, RE + r"\Gamma^{\nu=\nu', \omega =0 }_{d}$", False )
pl.xlabel(r"$\nu$", fontsize = 10)
plotgamma_axis_m( pl.subplot(2,2,4),- regamma_xph0, - regamma_xph0,- regamma_xph3, RE + r"\Gamma^{\nu=\nu', \omega =0 }_{m}$", False )
pl.xlabel(r"$\nu$", fontsize = 10)

pl.tight_layout()

#--- Save to file
pl.savefig("plots/GAMMA_axis_nu'=nu_W0.png", dpi = 200)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()
