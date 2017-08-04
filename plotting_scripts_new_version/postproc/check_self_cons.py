#!/usr/bin/python

#=============================================================================================
#
#                       NOTE: Plots of the comparison between the asymptotic functions calculate exactly
#                             by means of Vienna's TU ED code and the the ones obtain by means of 
#                             the self-consistency
#
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
vert_mul = 1
pi = math.pi

#--------------------------------------MANAGING FILES ------------------------------------------

def run(command):
    output = subprocess.check_output(command, shell=True)
    return output

most_recently_edited = run("ls -Art dat/ | tail -n 1")

fname1 = "dat/dat_U1_Beta50_PFCB1000_PARQ_SU2_ED_DIVERGENT.h5"
#fname2 = "dat/H5FILES/BETA50/4SITES/U1/POSTPROC/dat_U1_Beta50_PFCB800_PARQ_SU2_SELFCON.h5"
fname2 = "dat/dat_U1_Beta50_PFCB800_PARQ_SU2_SELFCON_new_version.h5"

if len(sys.argv) > 1:
    fname1 = str(sys.argv[1])

fname1 = fname1.rstrip('\n') # strip newline of fname
f1 = h5py.File(fname1, "r")

if len(sys.argv) > 1:
    fname2 = str(sys.argv[1])

fname2 = fname2.rstrip('\n') # strip newline of fname
f2 = h5py.File(fname2, "r")

os.system('mkdir -p log')
os.system('mkdir -p plots')

#os.system('rm plots/prev_Vert.png 2> log/plot.log')
#os.system('rm plots/prev_phi.png 2> log/plot.log')
#os.system('rm plots/prev_P_func.png 2> log/plot.log')
#os.system('rm plots/prev_R_func.png 2> log/plot.log')

#os.system('mv plots/Vert.png plots/prev_Vert.png 2> log/plot.log')
#os.system('mv plots/phi.png plots/prev_phi.png 2> log/plot.log')
#os.system('mv plots/P_func.png plots/prev_P_func.png 2> log/plot.log')
#os.system('mv plots/R_func.png plots/prev_R_func.png 2> log/plot.log')


#--------------------------------------READ PARAMETERS FROM FILE ------------------------------------------

parVals = f1["/Params"].attrs.values()

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
#pl.rc('text', usetex=True)
#pl.rc('text.latex', preamble='\usepackage{amsmath}')

RE = r"$\operatorname{Re}"
IM = r"$\operatorname{Im}"

#--------------------------------------Chi PLOTTING ------------------------------------------

print("Plotting chi functions ...")


#--- Read ED K1
bosgrid1 = np.array(f1["/chi_func/bgrid"])
rechi_pp1 = vert_mul * np.array(f1["/chi_func/RE_PP"])
imchi_pp1 = vert_mul * np.array(f1["/chi_func/IM_PP"])
rechi_ph1 = vert_mul * np.array(f1["/chi_func/RE_PH"])
imchi_ph1 = vert_mul * np.array(f1["/chi_func/IM_PH"])
rechi_xph1 = vert_mul * np.array(f1["/chi_func/RE_XPH"])
imchi_xph1 = vert_mul * np.array(f1["/chi_func/IM_XPH"])
fdim_bos1 = bosgrid1.shape[0]


#--- Read SELF-CONSISTENT K1
bosgrid2 = np.array(f2["/chi_func/bgrid"])
rechi_pp2 = vert_mul * np.array(f2["/chi_func/RE_PP"])
imchi_pp2 = vert_mul * np.array(f2["/chi_func/IM_PP"])
rechi_ph2 = vert_mul * np.array(f2["/chi_func/RE_PH"])
imchi_ph2 = vert_mul * np.array(f2["/chi_func/IM_PH"])
rechi_xph2 = vert_mul * np.array(f2["/chi_func/RE_XPH"])
imchi_xph2 = vert_mul * np.array(f2["/chi_func/IM_XPH"])
fdim_bos2 = bosgrid2.shape[0]

fdim_bos_plot = min([fdim_bos2,fdim_bos1])

bosgrid_plot = np.array([(2*n)*pi/BETA for n in range(-(fdim_bos_plot-1)/2, (fdim_bos_plot-1)/2 + 1)])

x_range_fact = 1.0

#--- Helper functions
def plotchi( use_pl, arr1,arr2, title ):
    pl.plot(bosgrid_plot, np.array([arr1[i+(fdim_bos1-1)/2,0,0,0,0,0]-arr2[i+(fdim_bos2-1)/2,0,0,0,0,0] for i in range(-(fdim_bos_plot-1)/2, (fdim_bos_plot-1)/2 + 1)]), 'b-', ms=3, mew=0.2)
    pl.xlim([x_range_fact*min(bosgrid_plot),x_range_fact*max(bosgrid_plot)])
    use_pl.set_title(title)
    return

pl.suptitle(r"$U=$" + str(UINT) + r"     $\beta=$" + str(BETA) + r"     $\epsilon=$" + str(EPS) + "  BEPs ")

plotchi( pl.subplot(2,3,1), rechi_pp1 - rechi_pp1, rechi_pp2 - rechi_pp2, RE + r"\delta \chi^{PP}_{\uparrow\uparrow}$" )
plotchi( pl.subplot(2,3,2), rechi_ph1 - rechi_xph1, rechi_ph2 - rechi_xph2, RE + r"\delta \chi^{PH}_{\uparrow\uparrow}$" )
plotchi( pl.subplot(2,3,3), rechi_xph1 - rechi_ph1, rechi_xph2 - rechi_ph2, RE + r"\delta \chi^{XPH}_{\uparrow\uparrow}$" )

plotchi( pl.subplot(2,3,4), rechi_pp1,rechi_pp2, RE + r"\delta \chi^{PP}_{\uparrow\downarrow}$" )
pl.xlabel(r"$\Omega$")
plotchi( pl.subplot(2,3,5), rechi_ph1,rechi_ph2, RE + r"\delta \chi^{PH}_{\uparrow\downarrow}$" )
pl.xlabel(r"$\Omega$")
plotchi( pl.subplot(2,3,6), rechi_xph1,rechi_xph2, RE + r"\delta \chi^{XPH}_{\uparrow\downarrow}$" )
pl.xlabel(r"$\Omega$")

pl.tight_layout()

pl.savefig("plots/delta_chi.png", dpi=150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()

#--------------------------------------P PLOTTING ------------------------------------------

print("Plotting P ...")

#--- Read
reP_pp1 = vert_mul * np.array(f1["/P_func/RE_PP"])
imP_pp1 = vert_mul * np.array(f1["/P_func/IM_PP"])
reP_ph1 = vert_mul * np.array(f1["/P_func/RE_PH"])
imP_ph1 = vert_mul * np.array(f1["/P_func/IM_PH"])
reP_xph1 = vert_mul * np.array(f1["/P_func/RE_XPH"])
imP_xph1 = vert_mul * np.array(f1["/P_func/IM_XPH"])

bdim1 = reP_pp1.shape[0]
fdim1 = reP_pp1.shape[1]

reP_pp2 = vert_mul * np.array(f2["/P_func/RE_PP"])
imP_pp2 = vert_mul * np.array(f2["/P_func/IM_PP"])
reP_ph2 = vert_mul * np.array(f2["/P_func/RE_PH"])
imP_ph2 = vert_mul * np.array(f2["/P_func/IM_PH"])
reP_xph2 = vert_mul * np.array(f2["/P_func/RE_XPH"])
imP_xph2 = vert_mul * np.array(f2["/P_func/IM_XPH"])

bdim2 = reP_pp2.shape[0]
fdim2 = reP_pp2.shape[1]

Pgrid2 = np.array(f2["/P_func/fgrid"])
Pgrid2 = np.array(f2["/P_func/fgrid"])

fdim_plot = min([fdim1,fdim2])

Pgrid_plot = np.array([(2*n)*pi/BETA for n in range(-fdim_plot/2, fdim_plot/2)])

x_range_fact = 1.0

#--- Helper functions
def plotP_Omega( W, arr1,arr2, color ):
    pl.plot(Pgrid_plot + (W % 2) * math.pi/BETA, np.array([arr1[W + (bdim1-1)/2,n + fdim1/2,0,0,0,0,0,0]- arr2[W + (bdim2-1)/2,n + fdim2/2,0,0,0,0,0,0] for n in range(-fdim_plot/2, fdim_plot/2)]), color, ms=3, mew=0.2, label=r"$\Omega=" + str(W*2) + r"\pi/\beta$")

def plotP( use_pl, arr1,arr2, title, legend ):
    plotP_Omega( 0, arr1,arr2, 'r-' ) 
    plotP_Omega( 1, arr1,arr2, 'g-' ) 
    plotP_Omega( 2, arr1,arr2, 'b-' ) 
    plotP_Omega( 8, arr1,arr2, 'c-' ) 
    plotP_Omega( 16, arr1,arr2, 'm-' ) 
    #plotP_Omega( 32, arr, 'y-' ) 
    pl.xlim([x_range_fact*min(Pgrid_plot),x_range_fact*max(Pgrid_plot)])
    if ( legend ):
        pl.legend(prop={'size':7})

    use_pl.set_title(title)
    return

pl.suptitle(r"$U=$" + str(UINT) + r"     $\beta=$" + str(BETA) + r"     $\epsilon=$" + str(EPS) + "  BEPs ")

plotP( pl.subplot(2,3,1), reP_pp1 - reP_pp1, reP_pp2 - reP_pp2, RE + r"\delta P^{PP}_{\uparrow\uparrow}$", True )
plotP( pl.subplot(2,3,2), reP_ph1 - reP_xph1,reP_ph2 - reP_xph2, RE + r"\delta P^{PH}_{\uparrow\uparrow}$", False )
plotP( pl.subplot(2,3,3), reP_xph1- reP_ph1,reP_xph2- reP_ph2, RE + r"\delta P^{XPH}_{\uparrow\uparrow}$", False )

plotP( pl.subplot(2,3,4), reP_pp1, reP_pp2, RE + r"\delta P^{PP}_{\uparrow\downarrow}$", False )
pl.xlabel(r"$\omega$")
plotP( pl.subplot(2,3,5), reP_ph1, reP_ph2, RE + r"\delta P^{PH}_{\uparrow\downarrow}$", False )
pl.xlabel(r"$\omega$")
plotP( pl.subplot(2,3,6), reP_xph1, reP_xph2, RE + r"\delta P^{XPH}_{\uparrow\downarrow}$", False )
pl.xlabel(r"$\omega$")

pl.tight_layout()

pl.savefig("plots/delta_P.png", dpi=150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()

