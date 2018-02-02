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

#fname0 = "../../dat/dat_U1_Beta50_PFCB0_PARQ_SU2_METH2_INVR1xCOUNT_BIGINVR_CORR_ED_ONESHOT.h5"
#fname1 = "../../dat/H5FILES/BETA50/4SITES/U1p75/POSTPROC/dat_U1p75_Beta50_PFCB1000_PARQ_SU2_METH2_INVR1xCOUNT_CORR_ED.h5"
#fname2 = "../../dat/H5FILES/BETA50/4SITES/U1p75/POSTPROC/dat_U1p75_Beta50_PFCB1000_PARQ_SU2_METH2_INVR1xCOUNT_ED.h5"
#fname3 = "../../dat/H5FILES/BETA50/4SITES/U1p75/POSTPROC/dat_U1p75_Beta50_PFCB1000_PARQ_SU2_METH2_INVR1xCOUNTm5_CORR_ED.h5"
#fname4 = "../../dat/H5FILES/BETA50/4SITES/U1p75/POSTPROC/dat_U1p75_Beta50_PFCB1000_PARQ_SU2_METH2_INVR1xCOUNTm5_ED.h5"
#fname5 = "../../dat/H5FILES/BETA50/4SITES/U1p75/POSTPROC/dat_U1p75_Beta50_PFCB1000_PARQ_SU2_METH2_INVR1xCOUNTm10_CORR_ED.h5"
#fname6 = "../../dat/H5FILES/BETA50/4SITES/U1p75/POSTPROC/dat_U1p75_Beta50_PFCB1000_PARQ_SU2_METH2_INVR1xCOUNTm10_ED.h5"
#fname7 = "../../dat/H5FILES/BETA50/4SITES/U1p75/POSTPROC/dat_U1p75_Beta50_PFCB1000_PARQ_SU2_METH2_INVR1xCOUNTm15_CORR_ED.h5"
#fname8 = "../../dat/H5FILES/BETA50/4SITES/U1p75/POSTPROC/dat_U1p75_Beta50_PFCB1000_PARQ_SU2_METH2_INVR1xCOUNTm15_ED.h5"
#fname9 = "../../dat/H5FILES/BETA50/4SITES/U1p75/POSTPROC/dat_U1p75_Beta50_PFCB1000_PARQ_SU2_METH2_INVR1xCOUNTm20_CORR_ED.h5"
#fname10 = "../../dat/H5FILES/BETA50/4SITES/U1p75/POSTPROC/dat_U1p75_Beta50_PFCB1000_PARQ_SU2_METH2_INVR1xCOUNTm20_ED.h5"
#fname11 = "../../dat/H5FILES/BETA50/4SITES/U1p75/POSTPROC/dat_U1p75_Beta50_PFCB1000_PARQ_SU2_METH2_INVR1xCOUNTm30_CORR_ED.h5"
#fname12 = "../../dat/H5FILES/BETA50/4SITES/U1p75/POSTPROC/dat_U1p75_Beta50_PFCB1000_PARQ_SU2_METH2_INVR1xCOUNTm30_ED.h5"
#fname13 = "../../dat/H5FILES/BETA50/4SITES/U1p75/POSTPROC/dat_U1p75_Beta50_PFCB1000_PARQ_SU2_METH2_INVR1xCOUNTm40_CORR_ED.h5"
#fname14 = "../../dat/H5FILES/BETA50/4SITES/U1p75/POSTPROC/dat_U1p75_Beta50_PFCB1000_PARQ_SU2_METH2_INVR1xCOUNTm40_ED.h5"
#fname15 = "../../dat/H5FILES/BETA50/4SITES/U1p75/POSTPROC/dat_U1p75_Beta50_PFCB1000_PARQ_SU2_METH2_INVR1xCOUNTm50_CORR_ED.h5"
#fname16 = "../../dat/H5FILES/BETA50/4SITES/U1p75/POSTPROC/dat_U1p75_Beta50_PFCB1000_PARQ_SU2_METH2_INVR1xCOUNTm50_ED.h5"
#
#fname17 = "../../dat/H5FILES/BETA50/4SITES/U1p75/POSTPROC/dat_U1p75_Beta50_PFCB1000_PARQ_SU2_METH2_INVR1xCOUNT_W20_CORR_ED.h5"
#fname18 = "../../dat/H5FILES/BETA50/4SITES/U1p75/POSTPROC/dat_U1p75_Beta50_PFCB1000_PARQ_SU2_METH2_INVR1xCOUNT_W20_ED.h5"
#fname19 = "../../dat/H5FILES/BETA50/4SITES/U1p75/POSTPROC/dat_U1p75_Beta50_PFCB1000_PARQ_SU2_METH2_INVR1xCOUNTm20_W20_CORR_ED.h5"
#fname20 = "../../dat/H5FILES/BETA50/4SITES/U1p75/POSTPROC/dat_U1p75_Beta50_PFCB1000_PARQ_SU2_METH2_INVR1xCOUNTm20_W20_ED.h5"
#fname21 = "../../dat/H5FILES/BETA50/4SITES/U1p75/POSTPROC/dat_U1p75_Beta50_PFCB1000_PARQ_SU2_METH2_INVR1xCOUNTm40_W20_CORR_ED.h5"
#fname22 = "../../dat/H5FILES/BETA50/4SITES/U1p75/POSTPROC/dat_U1p75_Beta50_PFCB1000_PARQ_SU2_METH2_INVR1xCOUNTm40_W20_ED.h5"

fname17b = "../../dat/H5FILES/BETA50/4SITES/U1p75/POSTPROC/dat_U1p75_Beta50_PFCB1000_PARQ_SU2_METH2_INVR1xCOUNT_W20_CORR_ED.h5"
fname18b = "../../dat/H5FILES/BETA50/4SITES/U1p75/POSTPROC/dat_U1p75_Beta50_PFCB1000_PARQ_SU2_METH2_INVR1xCOUNTm_W20_ED.h5"
fname19b = "../../dat/H5FILES/BETA50/4SITES/U1p75/POSTPROC/dat_U1p75_Beta50_PFCB1000_PARQ_SU2_METH2_INVR1xCOUNTm20_W20_CORR_ED.h5"
fname20b = "../../dat/H5FILES/BETA50/4SITES/U1p75/POSTPROC/dat_U1p75_Beta50_PFCB1000_PARQ_SU2_METH2_INVR1xCOUNTm20_W20_ED.h5"
fname21b = "../../dat/H5FILES/BETA50/4SITES/U1p75/POSTPROC/dat_U1p75_Beta50_PFCB1000_PARQ_SU2_METH2_INVR1xCOUNTm40_W20_CORR_ED.h5"
fname22b = "../../dat/H5FILES/BETA50/4SITES/U1p75/POSTPROC/dat_U1p75_Beta50_PFCB1000_PARQ_SU2_METH2_INVR1xCOUNTm40_W20_ED.h5"

fname1 = "../../dat/H5FILES/BETA50/4SITES/U1p75/POSTPROC/dat_U1p75_Beta50_PFCB1000_PARQ_SU2_METH1_INVR1xCOUNT_ASYR10xINVR_CORR_ED.h5"
fname2 = "../../dat/H5FILES/BETA50/4SITES/U1p75/POSTPROC/dat_U1p75_Beta50_PFCB1000_PARQ_SU2_METH2_INVR1xCOUNT_ED.h5"
fname3 = "../../dat/H5FILES/BETA50/4SITES/U1p75/POSTPROC/dat_U1p75_Beta50_PFCB1000_PARQ_SU2_METH1_INVR1xCOUNTm5_ASYR10xINVR_CORR_ED.h5"
fname4 = "../../dat/H5FILES/BETA50/4SITES/U1p75/POSTPROC/dat_U1p75_Beta50_PFCB1000_PARQ_SU2_METH2_INVR1xCOUNTm5_ED.h5"
fname5 = "../../dat/H5FILES/BETA50/4SITES/U1p75/POSTPROC/dat_U1p75_Beta50_PFCB1000_PARQ_SU2_METH1_INVR1xCOUNTm10_ASYR10xINVR_CORR_ED.h5"
fname6 = "../../dat/H5FILES/BETA50/4SITES/U1p75/POSTPROC/dat_U1p75_Beta50_PFCB1000_PARQ_SU2_METH2_INVR1xCOUNTm10_ED.h5"
fname7 = "../../dat/H5FILES/BETA50/4SITES/U1p75/POSTPROC/dat_U1p75_Beta50_PFCB1000_PARQ_SU2_METH1_INVR1xCOUNTm15_ASYR10xINVR_CORR_ED.h5"
fname8 = "../../dat/H5FILES/BETA50/4SITES/U1p75/POSTPROC/dat_U1p75_Beta50_PFCB1000_PARQ_SU2_METH2_INVR1xCOUNTm15_ED.h5"
fname9 = "../../dat/H5FILES/BETA50/4SITES/U1p75/POSTPROC/dat_U1p75_Beta50_PFCB1000_PARQ_SU2_METH1_INVR1xCOUNTm20_ASYR10xINVR_CORR_ED.h5"
fname10 = "../../dat/H5FILES/BETA50/4SITES/U1p75/POSTPROC/dat_U1p75_Beta50_PFCB1000_PARQ_SU2_METH2_INVR1xCOUNTm20_ED.h5"
fname11 = "../../dat/H5FILES/BETA50/4SITES/U1p75/POSTPROC/dat_U1p75_Beta50_PFCB1000_PARQ_SU2_METH1_INVR1xCOUNTm30_ASYR10xINVR_CORR_ED.h5"
fname12 = "../../dat/H5FILES/BETA50/4SITES/U1p75/POSTPROC/dat_U1p75_Beta50_PFCB1000_PARQ_SU2_METH2_INVR1xCOUNTm30_ED.h5"
fname13 = "../../dat/H5FILES/BETA50/4SITES/U1p75/POSTPROC/dat_U1p75_Beta50_PFCB1000_PARQ_SU2_METH1_INVR1xCOUNTm40_ASYR10xINVR_CORR_ED.h5"
fname14 = "../../dat/H5FILES/BETA50/4SITES/U1p75/POSTPROC/dat_U1p75_Beta50_PFCB1000_PARQ_SU2_METH2_INVR1xCOUNTm40_ED.h5"
fname15 = "../../dat/H5FILES/BETA50/4SITES/U1p75/POSTPROC/dat_U1p75_Beta50_PFCB1000_PARQ_SU2_METH1_INVR1xCOUNTm50_ASYR10xINVR_CORR_ED.h5"
fname16 = "../../dat/H5FILES/BETA50/4SITES/U1p75/POSTPROC/dat_U1p75_Beta50_PFCB1000_PARQ_SU2_METH2_INVR1xCOUNTm50_ED.h5"

fname17 = "../../dat/H5FILES/BETA50/4SITES/U1p75/POSTPROC/dat_U1p75_Beta50_PFCB1000_PARQ_SU2_METH1_INVR1xCOUNT_ASYR10xINVR_W0_CORR_ED.h5"
fname18 = "../../dat/H5FILES/BETA50/4SITES/U1p75/POSTPROC/dat_U1p75_Beta50_PFCB1000_PARQ_SU2_METH2_INVR1xCOUNT_W0_ED.h5"
fname19 = "../../dat/H5FILES/BETA50/4SITES/U1p75/POSTPROC/dat_U1p75_Beta50_PFCB1000_PARQ_SU2_METH1_INVR1xCOUNTm20_ASYR10xINVR_W0_CORR_ED.h5"
fname20 = "../../dat/H5FILES/BETA50/4SITES/U1p75/POSTPROC/dat_U1p75_Beta50_PFCB1000_PARQ_SU2_METH2_INVR1xCOUNTm20_W0_ED.h5"
fname21 = "../../dat/H5FILES/BETA50/4SITES/U1p75/POSTPROC/dat_U1p75_Beta50_PFCB1000_PARQ_SU2_METH1_INVR1xCOUNTm40_ASYR10xINVR_W0_CORR_ED.h5"
fname22 = "../../dat/H5FILES/BETA50/4SITES/U1p75/POSTPROC/dat_U1p75_Beta50_PFCB1000_PARQ_SU2_METH2_INVR1xCOUNTm40_W0_ED.h5"

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
if len(sys.argv) > 1:
    fname7 = str(sys.argv[1])
if len(sys.argv) > 1:
    fname8 = str(sys.argv[1])
if len(sys.argv) > 1:
    fname9 = str(sys.argv[1])
if len(sys.argv) > 1:
    fname10 = str(sys.argv[1])
if len(sys.argv) > 1:
    fname11 = str(sys.argv[1])
if len(sys.argv) > 1:
    fname12 = str(sys.argv[1])
if len(sys.argv) > 1:
    fname13 = str(sys.argv[1])
if len(sys.argv) > 1:
    fname14 = str(sys.argv[1])
if len(sys.argv) > 1:
    fname15 = str(sys.argv[1])
if len(sys.argv) > 1:
    fname16 = str(sys.argv[1])

if len(sys.argv) > 1:
    fname17 = str(sys.argv[1])
if len(sys.argv) > 1:
    fname18 = str(sys.argv[1])
if len(sys.argv) > 1:
    fname19 = str(sys.argv[1])
if len(sys.argv) > 1:
    fname20 = str(sys.argv[1])
if len(sys.argv) > 1:
    fname21 = str(sys.argv[1])
if len(sys.argv) > 1:
    fname22 = str(sys.argv[1])

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

fname7 = fname7.rstrip('\n') # strip newline of fname
f7 = h5py.File(fname7, "r")

fname8 = fname8.rstrip('\n') # strip newline of fname
f8 = h5py.File(fname8, "r")

fname9 = fname9.rstrip('\n') # strip newline of fname
f9 = h5py.File(fname9, "r")

fname10 = fname10.rstrip('\n') # strip newline of fname
f10 = h5py.File(fname10, "r")

fname11 = fname11.rstrip('\n') # strip newline of fname
f11 = h5py.File(fname11, "r")

fname12 = fname12.rstrip('\n') # strip newline of fname
f12 = h5py.File(fname12, "r")

fname13 = fname13.rstrip('\n') # strip newline of fname
f13 = h5py.File(fname13, "r")

fname14 = fname14.rstrip('\n') # strip newline of fname
f14 = h5py.File(fname14, "r")

fname15 = fname15.rstrip('\n') # strip newline of fname
f15 = h5py.File(fname15, "r")

fname16 = fname16.rstrip('\n') # strip newline of fname
f16 = h5py.File(fname16, "r")

fname17 = fname17.rstrip('\n') # strip newline of fname
f17 = h5py.File(fname17, "r")

fname18 = fname18.rstrip('\n') # strip newline of fname
f18 = h5py.File(fname18, "r")

fname19 = fname19.rstrip('\n') # strip newline of fname
f19 = h5py.File(fname19, "r")

fname20 = fname20.rstrip('\n') # strip newline of fname
f20 = h5py.File(fname20, "r")

fname21 = fname21.rstrip('\n') # strip newline of fname
f21 = h5py.File(fname21, "r")

fname22 = fname22.rstrip('\n') # strip newline of fname
f22 = h5py.File(fname22, "r")

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

#--- Read gamma method 2 (INVR1xVERTR-15)
regamma_pp7 = vert_mul * np.array(f7["/gamma_func/RE_PP"])
imgamma_pp7 = vert_mul * np.array(f7["/gamma_func/IM_PP"])
regamma_ph7 = vert_mul * np.array(f7["/gamma_func/RE_PH"])
imgamma_ph7 = vert_mul * np.array(f7["/gamma_func/IM_PH"])
regamma_xph7 = vert_mul * np.array(f7["/gamma_func/RE_XPH"])
imgamma_xph7 = vert_mul * np.array(f7["/gamma_func/IM_XPH"])

#--- Read gamma method 2 (INVR1xVERTR-15 NOCORR)
regamma_pp8 = vert_mul * np.array(f8["/gamma_func/RE_PP"])
imgamma_pp8 = vert_mul * np.array(f8["/gamma_func/IM_PP"])
regamma_ph8 = vert_mul * np.array(f8["/gamma_func/RE_PH"])
imgamma_ph8 = vert_mul * np.array(f8["/gamma_func/IM_PH"])
regamma_xph8 = vert_mul * np.array(f8["/gamma_func/RE_XPH"])
imgamma_xph8 = vert_mul * np.array(f8["/gamma_func/IM_XPH"])

#--- Read gamma method 2 (INVR1xVERTR-20)
regamma_pp9 = vert_mul * np.array(f9["/gamma_func/RE_PP"])
imgamma_pp9 = vert_mul * np.array(f9["/gamma_func/IM_PP"])
regamma_ph9 = vert_mul * np.array(f9["/gamma_func/RE_PH"])
imgamma_ph9 = vert_mul * np.array(f9["/gamma_func/IM_PH"])
regamma_xph9 = vert_mul * np.array(f9["/gamma_func/RE_XPH"])
imgamma_xph9 = vert_mul * np.array(f9["/gamma_func/IM_XPH"])

#--- Read gamma method 2 (INVR1xVERTR-15 NOCORR)
regamma_pp10 = vert_mul * np.array(f10["/gamma_func/RE_PP"])
imgamma_pp10 = vert_mul * np.array(f10["/gamma_func/IM_PP"])
regamma_ph10 = vert_mul * np.array(f10["/gamma_func/RE_PH"])
imgamma_ph10 = vert_mul * np.array(f10["/gamma_func/IM_PH"])
regamma_xph10 = vert_mul * np.array(f10["/gamma_func/RE_XPH"])
imgamma_xph10 = vert_mul * np.array(f10["/gamma_func/IM_XPH"])

#--- Read gamma method 2 (INVR1xVERTR-20)
regamma_pp11 = vert_mul * np.array(f11["/gamma_func/RE_PP"])
imgamma_pp11 = vert_mul * np.array(f11["/gamma_func/IM_PP"])
regamma_ph11 = vert_mul * np.array(f11["/gamma_func/RE_PH"])
imgamma_ph11 = vert_mul * np.array(f11["/gamma_func/IM_PH"])
regamma_xph11 = vert_mul * np.array(f11["/gamma_func/RE_XPH"])
imgamma_xph11 = vert_mul * np.array(f11["/gamma_func/IM_XPH"])

#--- Read gamma method 2 (INVR1xVERTR-15 NOCORR)
regamma_pp12 = vert_mul * np.array(f12["/gamma_func/RE_PP"])
imgamma_pp12 = vert_mul * np.array(f12["/gamma_func/IM_PP"])
regamma_ph12 = vert_mul * np.array(f12["/gamma_func/RE_PH"])
imgamma_ph12 = vert_mul * np.array(f12["/gamma_func/IM_PH"])
regamma_xph12 = vert_mul * np.array(f12["/gamma_func/RE_XPH"])
imgamma_xph12 = vert_mul * np.array(f12["/gamma_func/IM_XPH"])

#--- Read gamma method 2 (INVR1xVERTR-20)
regamma_pp13 = vert_mul * np.array(f13["/gamma_func/RE_PP"])
imgamma_pp13 = vert_mul * np.array(f13["/gamma_func/IM_PP"])
regamma_ph13 = vert_mul * np.array(f13["/gamma_func/RE_PH"])
imgamma_ph13 = vert_mul * np.array(f13["/gamma_func/IM_PH"])
regamma_xph13 = vert_mul * np.array(f13["/gamma_func/RE_XPH"])
imgamma_xph13 = vert_mul * np.array(f13["/gamma_func/IM_XPH"])

#--- Read gamma method 2 (INVR1xVERTR-15 NOCORR)
regamma_pp14 = vert_mul * np.array(f14["/gamma_func/RE_PP"])
imgamma_pp14 = vert_mul * np.array(f14["/gamma_func/IM_PP"])
regamma_ph14 = vert_mul * np.array(f14["/gamma_func/RE_PH"])
imgamma_ph14 = vert_mul * np.array(f14["/gamma_func/IM_PH"])
regamma_xph14 = vert_mul * np.array(f14["/gamma_func/RE_XPH"])
imgamma_xph14 = vert_mul * np.array(f14["/gamma_func/IM_XPH"])

#--- Read gamma method 2 (INVR1xVERTR-20)
regamma_pp15 = vert_mul * np.array(f15["/gamma_func/RE_PP"])
imgamma_pp15 = vert_mul * np.array(f15["/gamma_func/IM_PP"])
regamma_ph15 = vert_mul * np.array(f15["/gamma_func/RE_PH"])
imgamma_ph15 = vert_mul * np.array(f15["/gamma_func/IM_PH"])
regamma_xph15 = vert_mul * np.array(f15["/gamma_func/RE_XPH"])
imgamma_xph15 = vert_mul * np.array(f15["/gamma_func/IM_XPH"])

#--- Read gamma method 2 (INVR1xVERTR-15 NOCORR)
regamma_pp16 = vert_mul * np.array(f16["/gamma_func/RE_PP"])
imgamma_pp16 = vert_mul * np.array(f16["/gamma_func/IM_PP"])
regamma_ph16 = vert_mul * np.array(f16["/gamma_func/RE_PH"])
imgamma_ph16 = vert_mul * np.array(f16["/gamma_func/IM_PH"])
regamma_xph16 = vert_mul * np.array(f16["/gamma_func/RE_XPH"])
imgamma_xph16 = vert_mul * np.array(f16["/gamma_func/IM_XPH"])

#--- Read gamma method 2 (INVR1xVERTR-20)
regamma_pp17 = vert_mul * np.array(f17["/gamma_func/RE_PP"])
imgamma_pp17 = vert_mul * np.array(f17["/gamma_func/IM_PP"])
regamma_ph17 = vert_mul * np.array(f17["/gamma_func/RE_PH"])
imgamma_ph17 = vert_mul * np.array(f17["/gamma_func/IM_PH"])
regamma_xph17 = vert_mul * np.array(f17["/gamma_func/RE_XPH"])
imgamma_xph17 = vert_mul * np.array(f17["/gamma_func/IM_XPH"])

#--- Read gamma method 2 (INVR1xVERTR-15 NOCORR)
regamma_pp18 = vert_mul * np.array(f18["/gamma_func/RE_PP"])
imgamma_pp18 = vert_mul * np.array(f18["/gamma_func/IM_PP"])
regamma_ph18 = vert_mul * np.array(f18["/gamma_func/RE_PH"])
imgamma_ph18 = vert_mul * np.array(f18["/gamma_func/IM_PH"])
regamma_xph18 = vert_mul * np.array(f18["/gamma_func/RE_XPH"])
imgamma_xph18 = vert_mul * np.array(f18["/gamma_func/IM_XPH"])

#--- Read gamma method 2 (INVR1xVERTR-20)
regamma_pp19 = vert_mul * np.array(f19["/gamma_func/RE_PP"])
imgamma_pp19 = vert_mul * np.array(f19["/gamma_func/IM_PP"])
regamma_ph19 = vert_mul * np.array(f19["/gamma_func/RE_PH"])
imgamma_ph19 = vert_mul * np.array(f19["/gamma_func/IM_PH"])
regamma_xph19 = vert_mul * np.array(f19["/gamma_func/RE_XPH"])
imgamma_xph19 = vert_mul * np.array(f19["/gamma_func/IM_XPH"])

#--- Read gamma method 2 (INVR1xVERTR-15 NOCORR)
regamma_pp20 = vert_mul * np.array(f20["/gamma_func/RE_PP"])
imgamma_pp20 = vert_mul * np.array(f20["/gamma_func/IM_PP"])
regamma_ph20 = vert_mul * np.array(f20["/gamma_func/RE_PH"])
imgamma_ph20 = vert_mul * np.array(f20["/gamma_func/IM_PH"])
regamma_xph20 = vert_mul * np.array(f20["/gamma_func/RE_XPH"])
imgamma_xph20 = vert_mul * np.array(f20["/gamma_func/IM_XPH"])

#--- Read gamma method 2 (INVR1xVERTR-20)
regamma_pp21 = vert_mul * np.array(f21["/gamma_func/RE_PP"])
imgamma_pp21 = vert_mul * np.array(f21["/gamma_func/IM_PP"])
regamma_ph21 = vert_mul * np.array(f21["/gamma_func/RE_PH"])
imgamma_ph21 = vert_mul * np.array(f21["/gamma_func/IM_PH"])
regamma_xph21 = vert_mul * np.array(f21["/gamma_func/RE_XPH"])
imgamma_xph21 = vert_mul * np.array(f21["/gamma_func/IM_XPH"])

#--- Read gamma method 2 (INVR1xVERTR-15 NOCORR)
regamma_pp22 = vert_mul * np.array(f22["/gamma_func/RE_PP"])
imgamma_pp22 = vert_mul * np.array(f22["/gamma_func/IM_PP"])
regamma_ph22 = vert_mul * np.array(f22["/gamma_func/RE_PH"])
imgamma_ph22 = vert_mul * np.array(f22["/gamma_func/IM_PH"])
regamma_xph22 = vert_mul * np.array(f22["/gamma_func/RE_XPH"])
imgamma_xph22 = vert_mul * np.array(f22["/gamma_func/IM_XPH"])

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

bdim7= regamma_pp7.shape[0]
fdim7= regamma_pp7.shape[1]
fdimo27 = fdim7/2
print bdim7, fdim7

bdim8= regamma_pp8.shape[0]
fdim8= regamma_pp8.shape[1]
fdimo28 = fdim8/2
print bdim8, fdim8

bdim9= regamma_pp9.shape[0]
fdim9= regamma_pp9.shape[1]
fdimo29 = fdim9/2
print bdim9, fdim9

bdim10= regamma_pp10.shape[0]
fdim10= regamma_pp10.shape[1]
fdimo210 = fdim10/2
print bdim10, fdim10

bdim11= regamma_pp11.shape[0]
fdim11= regamma_pp11.shape[1]
fdimo211 = fdim11/2
print bdim11, fdim11

bdim12= regamma_pp12.shape[0]
fdim12= regamma_pp12.shape[1]
fdimo212 = fdim12/2
print bdim12, fdim12

bdim13= regamma_pp13.shape[0]
fdim13= regamma_pp13.shape[1]
fdimo213 = fdim13/2
print bdim13, fdim13

bdim14= regamma_pp14.shape[0]
fdim14= regamma_pp14.shape[1]
fdimo214 = fdim14/2
print bdim14, fdim14

bdim15= regamma_pp15.shape[0]
fdim15= regamma_pp15.shape[1]
fdimo215 = fdim15/2
print bdim15, fdim15

bdim16= regamma_pp16.shape[0]
fdim16= regamma_pp16.shape[1]
fdimo216 = fdim16/2
print bdim16, fdim16

bdim17= regamma_pp17.shape[0]
fdim17= regamma_pp17.shape[1]
fdimo217 = fdim17/2
print bdim17, fdim17

bdim18= regamma_pp18.shape[0]
fdim18= regamma_pp18.shape[1]
fdimo218 = fdim18/2
print bdim18, fdim18

bdim19= regamma_pp19.shape[0]
fdim19= regamma_pp19.shape[1]
fdimo219 = fdim19/2
print bdim19, fdim19

bdim20= regamma_pp20.shape[0]
fdim20= regamma_pp20.shape[1]
fdimo220 = fdim20/2
print bdim20, fdim20

bdim21= regamma_pp21.shape[0]
fdim21= regamma_pp21.shape[1]
fdimo221 = fdim21/2
print bdim21, fdim21

bdim22= regamma_pp22.shape[0]
fdim22= regamma_pp22.shape[1]
fdimo222 = fdim22/2
print bdim22, fdim22
#
fdim_min = min([fdim1,fdim3,fdim5,fdim7,fdim9,fdim11,fdim13,fdim15,fdim17,fdim19,fdim21])
print fdim_min
gammagrid_plot = np.array([(2*i+1)*pi/BETA for i in range(-fdim_min/2, fdim_min/2)])


shift=0
#--------------------------PLOT CONVERGENCE WITH RESPECT TO INVERSION RANGE ------------------------------------------

#INVROVR=[1,2,4]
INVROVR=[10,20,30,40,45,50,55,60,80,100,120]

#---  Helper functions (include translation from  nambu to physical CAUTION : USING TIME REVERSAL SYMM) 
def plotgamma_conv_axis( use_pl, arr1, arr2, arr3, arr4, arr5, arr6,arr7,arr8, arr9, arr10, arr11,string ):
    zarr = [ arr11[shift + (bdim15-1)/2,fdimo215,fdimo215,0,0,0,0,0,0,0], arr10[shift + (bdim13-1)/2,fdimo213,fdimo213,0,0,0,0,0,0,0],arr9[shift + (bdim11-1)/2,fdimo211,fdimo211,0,0,0,0,0,0,0],arr8[shift + (bdim9-1)/2,fdimo29,fdimo29,0,0,0,0,0,0,0], arr7[shift + (bdim7-1)/2,fdimo27,fdimo27,0,0,0,0,0,0,0], arr6[shift + (bdim5-1)/2,fdimo25,fdimo25,0,0,0,0,0,0,0], arr5[shift + (bdim3-1)/2,fdimo23,fdimo23,0,0,0,0,0,0,0], arr4[shift + (bdim1-1)/2,fdimo21,fdimo21,0,0,0,0,0,0,0], arr3[0,fdimo221, fdimo221, 0,0,0,0,0,0,0], arr2[0,fdimo219, fdimo219, 0,0,0,0,0,0,0],arr1[0,fdimo217, fdimo217, 0,0,0,0,0,0,0]]
    pl.plot( INVROVR, zarr, marker = 'o', linestyle='-', color='b', ms=3, mew=0.2)
    use_pl.set_title( string , fontsize=10 )

    return


#--- Plot along axis \nu'=0

#pl.suptitle(r"$U=$" + str(UINT) + r"     $\beta=$" + str(BETA) + r"     $\epsilon=$" + str(EPS) +  r"     $\Omega=$" + str(shift) + r"$*2\pi/\beta$" + r"     Notation: $\phi(\Omega,\omega,\omega')$")

plotgamma_conv_axis( pl.subplot(2,2,1), regamma_pp17 + regamma_pp17[:,:,::-1,:,:,:,:,:,:,:],regamma_pp19 + regamma_pp19[:,:,::-1,:,:,:,:,:,:,:],regamma_pp21 + regamma_pp21[:,:,::-1,:,:,:,:,:,:,:],  regamma_pp1 + regamma_pp1[:,:,::-1,:,:,:,:,:,:,:], regamma_pp3 + regamma_pp3[:,:,::-1,:,:,:,:,:,:,:], regamma_pp5 + regamma_pp5[:,:,::-1,:,:,:,:,:,:,:], regamma_pp7 + regamma_pp7[:,:,::-1,:,:,:,:,:,:,:],regamma_pp9 + regamma_pp9[:,:,::-1,:,:,:,:,:,:,:],regamma_pp11 + regamma_pp11[:,:,::-1,:,:,:,:,:,:,:], regamma_pp13 + regamma_pp13[:,:,::-1,:,:,:,:,:,:,:],regamma_pp15 + regamma_pp15[:,:,::-1,:,:,:,:,:,:,:], RE + r"\Gamma^{\nu=\nu'=0, \omega =0}_{s}$" ) # flip sign of w_out
pl.xlabel(r"$INVR$")
plotgamma_conv_axis( pl.subplot(2,2,2), regamma_pp17 - regamma_pp17[:,:,::-1,:,:,:,:,:,:,:],regamma_pp19 - regamma_pp19[:,:,::-1,:,:,:,:,:,:,:],regamma_pp21 - regamma_pp21[:,:,::-1,:,:,:,:,:,:,:],regamma_pp1 - regamma_pp1[:,:,::-1,:,:,:,:,:,:,:], regamma_pp3 - regamma_pp3[:,:,::-1,:,:,:,:,:,:,:], regamma_pp5 - regamma_pp5[:,:,::-1,:,:,:,:,:,:,:], regamma_pp7 - regamma_pp7[:,:,::-1,:,:,:,:,:,:,:], regamma_pp9 - regamma_pp9[:,:,::-1,:,:,:,:,:,:,:],regamma_pp11 - regamma_pp11[:,:,::-1,:,:,:,:,:,:,:],regamma_pp13 - regamma_pp13[:,:,::-1,:,:,:,:,:,:,:],regamma_pp15 - regamma_pp15[:,:,::-1,:,:,:,:,:,:,:], RE + r"\Gamma^{\nu=\nu'=0, \omega =0}_{t}$" ) # flip sign of w_out
pl.xlabel(r"$INVR$")
plotgamma_conv_axis( pl.subplot(2,2,3), 2*regamma_ph17 - regamma_xph17,2*regamma_ph19 - regamma_xph19,2*regamma_ph21 - regamma_xph21,2*regamma_ph1 - regamma_xph1, 2*regamma_ph3 - regamma_xph3, 2*regamma_ph5 - regamma_xph5, 2*regamma_ph7 - regamma_xph7,2*regamma_ph9 - regamma_xph9,2*regamma_ph11 - regamma_xph11,2*regamma_ph13 - regamma_xph13,2*regamma_ph15 - regamma_xph15, RE + r"\Gamma^{\nu=\nu'=0, \omega =0}_{d}$" )
pl.xlabel(r"$INVR$")
plotgamma_conv_axis( pl.subplot(2,2,4), - regamma_xph17,- regamma_xph19,- regamma_xph21,- regamma_xph1, - regamma_xph3, - regamma_xph5, - regamma_xph7,- regamma_xph9,- regamma_xph11,- regamma_xph13,- regamma_xph15, RE + r"\Gamma^{\nu = \nu'=0, \omega =0}_{m}$" )
pl.xlabel(r"$INVR$")

pl.tight_layout()

#--- Save to file
pl.savefig("plots/Gammas_METHOD1_conv_W0.png", dpi = 150)
#pl.savefig("plots/Gammas_METHOD1_conv.png", dpi = 150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()

