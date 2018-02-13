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

fname1 = "../../dat/H5FILES/BETA50/4SITES/U1/POSTPROC/dat_U1_Beta50_PFCB1000_PARQ_SU2_METH2_INVR1xCOUNT_CORR_ED_01_FULL.h5"
fname1b = "../../dat/H5FILES/BETA50/4SITES/U1/POSTPROC/dat_U1_Beta50_PFCB1000_PARQ_SU2_METH1_INVR1xCOUNT_ASYR10xINVR_CORR_ED_01_FULL.h5"
fname2 = "../../dat/H5FILES/BETA50/4SITES/U1/POSTPROC/dat_U1_Beta50_PFCB1000_PARQ_SU2_METH2_INVR1xCOUNT_ED.h5"
fname3 = "../../dat/H5FILES/BETA50/4SITES/U1/POSTPROC/dat_U1_Beta50_PFCB1000_PARQ_SU2_METH2_INVR1xCOUNTm5_CORR_ED_01_FULL.h5"
fname3b = "../../dat/H5FILES/BETA50/4SITES/U1/POSTPROC/dat_U1_Beta50_PFCB1000_PARQ_SU2_METH1_INVR1xCOUNTm5_ASYR10xINVR_CORR_ED_01_FULL.h5"
fname4 = "../../dat/H5FILES/BETA50/4SITES/U1/POSTPROC/dat_U1_Beta50_PFCB1000_PARQ_SU2_METH2_INVR1xCOUNTm5_ED.h5"
fname5 = "../../dat/H5FILES/BETA50/4SITES/U1/POSTPROC/dat_U1_Beta50_PFCB1000_PARQ_SU2_METH2_INVR1xCOUNTm10_CORR_ED_01_FULL.h5"
fname5b = "../../dat/H5FILES/BETA50/4SITES/U1/POSTPROC/dat_U1_Beta50_PFCB1000_PARQ_SU2_METH1_INVR1xCOUNTm10_ASYR10xINVR_CORR_ED_01_FULL.h5"
fname6 = "../../dat/H5FILES/BETA50/4SITES/U1/POSTPROC/dat_U1_Beta50_PFCB1000_PARQ_SU2_METH2_INVR1xCOUNTm10_ED.h5"
fname7 = "../../dat/H5FILES/BETA50/4SITES/U1/POSTPROC/dat_U1_Beta50_PFCB1000_PARQ_SU2_METH2_INVR1xCOUNTm15_CORR_ED_01_FULL.h5"
fname7b = "../../dat/H5FILES/BETA50/4SITES/U1/POSTPROC/dat_U1_Beta50_PFCB1000_PARQ_SU2_METH1_INVR1xCOUNTm15_ASYR10xINVR_CORR_ED_01_FULL.h5"
fname8 = "../../dat/H5FILES/BETA50/4SITES/U1/POSTPROC/dat_U1_Beta50_PFCB1000_PARQ_SU2_METH2_INVR1xCOUNTm15_ED.h5"
fname9 = "../../dat/H5FILES/BETA50/4SITES/U1/POSTPROC/dat_U1_Beta50_PFCB1000_PARQ_SU2_METH2_INVR1xCOUNTm20_CORR_ED_01_FULL.h5"
fname9b = "../../dat/H5FILES/BETA50/4SITES/U1/POSTPROC/dat_U1_Beta50_PFCB1000_PARQ_SU2_METH1_INVR1xCOUNTm20_ASYR10xINVR_CORR_ED_01_FULL.h5"
fname10 = "../../dat/H5FILES/BETA50/4SITES/U1/POSTPROC/dat_U1_Beta50_PFCB1000_PARQ_SU2_METH2_INVR1xCOUNTm20_ED.h5"
fname11 = "../../dat/H5FILES/BETA50/4SITES/U1/POSTPROC/dat_U1_Beta50_PFCB1000_PARQ_SU2_METH2_INVR1xCOUNTm30_CORR_ED_01_FULL.h5"
fname11b = "../../dat/H5FILES/BETA50/4SITES/U1/POSTPROC/dat_U1_Beta50_PFCB1000_PARQ_SU2_METH1_INVR1xCOUNTm30_ASYR10xINVR_CORR_ED_01_FULL.h5"
fname12 = "../../dat/H5FILES/BETA50/4SITES/U1/POSTPROC/dat_U1_Beta50_PFCB1000_PARQ_SU2_METH2_INVR1xCOUNTm30_ED.h5"
fname13 = "../../dat/H5FILES/BETA50/4SITES/U1/POSTPROC/dat_U1_Beta50_PFCB1000_PARQ_SU2_METH2_INVR1xCOUNTm40_CORR_ED_01_FULL.h5"
fname13b = "../../dat/H5FILES/BETA50/4SITES/U1/POSTPROC/dat_U1_Beta50_PFCB1000_PARQ_SU2_METH1_INVR1xCOUNTm40_ASYR10xINVR_CORR_ED_01_FULL.h5"
fname14 = "../../dat/H5FILES/BETA50/4SITES/U1/POSTPROC/dat_U1_Beta50_PFCB1000_PARQ_SU2_METH2_INVR1xCOUNTm40_ED.h5"
fname15 = "../../dat/H5FILES/BETA50/4SITES/U1/POSTPROC/dat_U1_Beta50_PFCB1000_PARQ_SU2_METH2_INVR1xCOUNTm50_CORR_ED_01_FULL.h5"
fname15b = "../../dat/H5FILES/BETA50/4SITES/U1/POSTPROC/dat_U1_Beta50_PFCB1000_PARQ_SU2_METH1_INVR1xCOUNTm50_ASYR10xINVR_CORR_ED_01_FULL.h5"
fname16 = "../../dat/H5FILES/BETA50/4SITES/U1/POSTPROC/dat_U1_Beta50_PFCB1000_PARQ_SU2_METH2_INVR1xCOUNTm50_ED.h5"

fname17 = "../../dat/H5FILES/BETA50/4SITES/U1/POSTPROC/dat_U1_Beta50_PFCB1000_PARQ_SU2_METH2_INVR1xCOUNT_W20_CORR_ED.h5"
fname17b = "../../dat/H5FILES/BETA50/4SITES/U1/POSTPROC/dat_U1_Beta50_PFCB1000_PARQ_SU2_METH1_INVR1xCOUNT_ASYR10xINVR_W20_CORR_ED.h5"
fname18 = "../../dat/H5FILES/BETA50/4SITES/U1/POSTPROC/dat_U1_Beta50_PFCB1000_PARQ_SU2_METH2_INVR1xCOUNT_W20_ED.h5"
fname19 = "../../dat/H5FILES/BETA50/4SITES/U1/POSTPROC/dat_U1_Beta50_PFCB1000_PARQ_SU2_METH2_INVR1xCOUNTm20_W20_CORR_ED.h5"
fname19b = "../../dat/H5FILES/BETA50/4SITES/U1/POSTPROC/dat_U1_Beta50_PFCB1000_PARQ_SU2_METH1_INVR1xCOUNTm20_ASYR10xINVR_W20_CORR_ED.h5"
fname20 = "../../dat/H5FILES/BETA50/4SITES/U1/POSTPROC/dat_U1_Beta50_PFCB1000_PARQ_SU2_METH2_INVR1xCOUNTm20_W20_ED.h5"
fname21 = "../../dat/H5FILES/BETA50/4SITES/U1/POSTPROC/dat_U1_Beta50_PFCB1000_PARQ_SU2_METH2_INVR1xCOUNTm40_W20_CORR_ED.h5"
fname21b = "../../dat/H5FILES/BETA50/4SITES/U1/POSTPROC/dat_U1_Beta50_PFCB1000_PARQ_SU2_METH1_INVR1xCOUNTm40_ASYR10xINVR_W20_CORR_ED.h5"
fname22 = "../../dat/H5FILES/BETA50/4SITES/U1/POSTPROC/dat_U1_Beta50_PFCB1000_PARQ_SU2_METH2_INVR1xCOUNTm40_W20_ED.h5"


#------ DECIDE WHICH COMPARISON YOU WANT TO ANALYZE


if len(sys.argv) > 1:
    fname1 = str(sys.argv[1])
if len(sys.argv) > 1:
    fname1b = str(sys.argv[1])
if len(sys.argv) > 1:
    fname2 = str(sys.argv[1])
if len(sys.argv) > 1:
    fname3 = str(sys.argv[1])
if len(sys.argv) > 1:
    fname3b = str(sys.argv[1])
if len(sys.argv) > 1:
    fname4 = str(sys.argv[1])
if len(sys.argv) > 1:
    fname5 = str(sys.argv[1])
if len(sys.argv) > 1:
    fname5b = str(sys.argv[1])
if len(sys.argv) > 1:
    fname6 = str(sys.argv[1])
if len(sys.argv) > 1:
    fname7 = str(sys.argv[1])
if len(sys.argv) > 1:
    fname7b = str(sys.argv[1])
if len(sys.argv) > 1:
    fname8 = str(sys.argv[1])
if len(sys.argv) > 1:
    fname9 = str(sys.argv[1])
if len(sys.argv) > 1:
    fname9b = str(sys.argv[1])
if len(sys.argv) > 1:
    fname10 = str(sys.argv[1])
if len(sys.argv) > 1:
    fname11 = str(sys.argv[1])
if len(sys.argv) > 1:
    fname11b = str(sys.argv[1])
if len(sys.argv) > 1:
    fname12 = str(sys.argv[1])
if len(sys.argv) > 1:
    fname13 = str(sys.argv[1])
if len(sys.argv) > 1:
    fname13b = str(sys.argv[1])
if len(sys.argv) > 1:
    fname14 = str(sys.argv[1])
if len(sys.argv) > 1:
    fname15 = str(sys.argv[1])
if len(sys.argv) > 1:
    fname15b = str(sys.argv[1])
if len(sys.argv) > 1:
    fname16 = str(sys.argv[1])

if len(sys.argv) > 1:
    fname17 = str(sys.argv[1])
if len(sys.argv) > 1:
    fname17b = str(sys.argv[1])
if len(sys.argv) > 1:
    fname18 = str(sys.argv[1])
if len(sys.argv) > 1:
    fname19 = str(sys.argv[1])
if len(sys.argv) > 1:
    fname19b = str(sys.argv[1])
if len(sys.argv) > 1:
    fname20 = str(sys.argv[1])
if len(sys.argv) > 1:
    fname21 = str(sys.argv[1])
if len(sys.argv) > 1:
    fname21b = str(sys.argv[1])
if len(sys.argv) > 1:
    fname22 = str(sys.argv[1])

fname1 = fname1.rstrip('\n') # strip newline of fname
f1 = h5py.File(fname1, "r")
fname1b = fname1b.rstrip('\n') # strip newline of fname
f1b = h5py.File(fname1b, "r")

fname2 = fname2.rstrip('\n') # strip newline of fname
f2 = h5py.File(fname2, "r")

fname3 = fname3.rstrip('\n') # strip newline of fname
f3 = h5py.File(fname3, "r")
fname3b = fname3b.rstrip('\n') # strip newline of fname
f3b = h5py.File(fname3b, "r")

fname4 = fname4.rstrip('\n') # strip newline of fname
f4 = h5py.File(fname4, "r")

fname5 = fname5.rstrip('\n') # strip newline of fname
f5 = h5py.File(fname5, "r")
fname5b = fname5b.rstrip('\n') # strip newline of fname
f5b = h5py.File(fname5b, "r")

fname6 = fname6.rstrip('\n') # strip newline of fname
f6 = h5py.File(fname6, "r")

fname7 = fname7.rstrip('\n') # strip newline of fname
f7 = h5py.File(fname7, "r")
fname7b = fname7b.rstrip('\n') # strip newline of fname
f7b = h5py.File(fname7b, "r")

fname8 = fname8.rstrip('\n') # strip newline of fname
f8 = h5py.File(fname8, "r")

fname9 = fname9.rstrip('\n') # strip newline of fname
f9 = h5py.File(fname9, "r")
fname9b = fname9b.rstrip('\n') # strip newline of fname
f9b = h5py.File(fname9b, "r")

fname10 = fname10.rstrip('\n') # strip newline of fname
f10 = h5py.File(fname10, "r")

fname11 = fname11.rstrip('\n') # strip newline of fname
f11 = h5py.File(fname11, "r")
fname11b = fname11b.rstrip('\n') # strip newline of fname
f11b = h5py.File(fname11b, "r")

fname12 = fname12.rstrip('\n') # strip newline of fname
f12 = h5py.File(fname12, "r")

fname13 = fname13.rstrip('\n') # strip newline of fname
f13 = h5py.File(fname13, "r")
fname13b = fname13b.rstrip('\n') # strip newline of fname
f13b = h5py.File(fname13b, "r")

fname14 = fname14.rstrip('\n') # strip newline of fname
f14 = h5py.File(fname14, "r")

fname15 = fname15.rstrip('\n') # strip newline of fname
f15 = h5py.File(fname15, "r")
fname15b = fname15b.rstrip('\n') # strip newline of fname
f15b = h5py.File(fname15b, "r")

fname16 = fname16.rstrip('\n') # strip newline of fname
f16 = h5py.File(fname16, "r")

fname17 = fname17.rstrip('\n') # strip newline of fname
f17 = h5py.File(fname17, "r")
fname17b = fname17b.rstrip('\n') # strip newline of fname
f17b = h5py.File(fname17b, "r")

fname18 = fname18.rstrip('\n') # strip newline of fname
f18 = h5py.File(fname18, "r")

fname19 = fname19.rstrip('\n') # strip newline of fname
f19 = h5py.File(fname19, "r")
fname19b = fname19b.rstrip('\n') # strip newline of fname
f19b = h5py.File(fname19b, "r")

fname20 = fname20.rstrip('\n') # strip newline of fname
f20 = h5py.File(fname20, "r")

fname21 = fname21.rstrip('\n') # strip newline of fname
f21 = h5py.File(fname21, "r")
fname21b = fname21b.rstrip('\n') # strip newline of fname
f21b = h5py.File(fname21b, "r")

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

shift= 20

#--------------------------------------GENERAL PLOT SETTINGS------------------------------------------

pl.rc('xtick', labelsize=10) 
pl.rc('ytick', labelsize=10) 
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

regamma_pp1b = vert_mul * np.array(f1b["/gamma_func/RE_PP"])
imgamma_pp1b = vert_mul * np.array(f1b["/gamma_func/IM_PP"])
regamma_ph1b = vert_mul * np.array(f1b["/gamma_func/RE_PH"])
imgamma_ph1b = vert_mul * np.array(f1b["/gamma_func/IM_PH"])
regamma_xph1b = vert_mul * np.array(f1b["/gamma_func/RE_XPH"])
imgamma_xph1b = vert_mul * np.array(f1b["/gamma_func/IM_XPH"])

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

regamma_pp3b = vert_mul * np.array(f3b["/gamma_func/RE_PP"])
imgamma_pp3b = vert_mul * np.array(f3b["/gamma_func/IM_PP"])
regamma_ph3b = vert_mul * np.array(f3b["/gamma_func/RE_PH"])
imgamma_ph3b = vert_mul * np.array(f3b["/gamma_func/IM_PH"])
regamma_xph3b = vert_mul * np.array(f3b["/gamma_func/RE_XPH"])
imgamma_xph3b = vert_mul * np.array(f3b["/gamma_func/IM_XPH"])

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

regamma_pp5b = vert_mul * np.array(f5b["/gamma_func/RE_PP"])
imgamma_pp5b = vert_mul * np.array(f5b["/gamma_func/IM_PP"])
regamma_ph5b = vert_mul * np.array(f5b["/gamma_func/RE_PH"])
imgamma_ph5b = vert_mul * np.array(f5b["/gamma_func/IM_PH"])
regamma_xph5b = vert_mul * np.array(f5b["/gamma_func/RE_XPH"])
imgamma_xph5b = vert_mul * np.array(f5b["/gamma_func/IM_XPH"])

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

regamma_pp7b = vert_mul * np.array(f7b["/gamma_func/RE_PP"])
imgamma_pp7b = vert_mul * np.array(f7b["/gamma_func/IM_PP"])
regamma_ph7b = vert_mul * np.array(f7b["/gamma_func/RE_PH"])
imgamma_ph7b = vert_mul * np.array(f7b["/gamma_func/IM_PH"])
regamma_xph7b = vert_mul * np.array(f7b["/gamma_func/RE_XPH"])
imgamma_xph7b = vert_mul * np.array(f7b["/gamma_func/IM_XPH"])

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

regamma_pp9b = vert_mul * np.array(f9b["/gamma_func/RE_PP"])
imgamma_pp9b = vert_mul * np.array(f9b["/gamma_func/IM_PP"])
regamma_ph9b = vert_mul * np.array(f9b["/gamma_func/RE_PH"])
imgamma_ph9b = vert_mul * np.array(f9b["/gamma_func/IM_PH"])
regamma_xph9b = vert_mul * np.array(f9b["/gamma_func/RE_XPH"])
imgamma_xph9b = vert_mul * np.array(f9b["/gamma_func/IM_XPH"])

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

regamma_pp11b = vert_mul * np.array(f11b["/gamma_func/RE_PP"])
imgamma_pp11b = vert_mul * np.array(f11b["/gamma_func/IM_PP"])
regamma_ph11b = vert_mul * np.array(f11b["/gamma_func/RE_PH"])
imgamma_ph11b = vert_mul * np.array(f11b["/gamma_func/IM_PH"])
regamma_xph11b = vert_mul * np.array(f11b["/gamma_func/RE_XPH"])
imgamma_xph11b = vert_mul * np.array(f11b["/gamma_func/IM_XPH"])

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

regamma_pp13b = vert_mul * np.array(f13b["/gamma_func/RE_PP"])
imgamma_pp13b = vert_mul * np.array(f13b["/gamma_func/IM_PP"])
regamma_ph13b = vert_mul * np.array(f13b["/gamma_func/RE_PH"])
imgamma_ph13b = vert_mul * np.array(f13b["/gamma_func/IM_PH"])
regamma_xph13b = vert_mul * np.array(f13b["/gamma_func/RE_XPH"])
imgamma_xph13b = vert_mul * np.array(f13b["/gamma_func/IM_XPH"])

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

regamma_pp15b = vert_mul * np.array(f15b["/gamma_func/RE_PP"])
imgamma_pp15b = vert_mul * np.array(f15b["/gamma_func/IM_PP"])
regamma_ph15b = vert_mul * np.array(f15b["/gamma_func/RE_PH"])
imgamma_ph15b = vert_mul * np.array(f15b["/gamma_func/IM_PH"])
regamma_xph15b = vert_mul * np.array(f15b["/gamma_func/RE_XPH"])
imgamma_xph15b = vert_mul * np.array(f15b["/gamma_func/IM_XPH"])

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

regamma_pp17b = vert_mul * np.array(f17b["/gamma_func/RE_PP"])
imgamma_pp17b = vert_mul * np.array(f17b["/gamma_func/IM_PP"])
regamma_ph17b = vert_mul * np.array(f17b["/gamma_func/RE_PH"])
imgamma_ph17b = vert_mul * np.array(f17b["/gamma_func/IM_PH"])
regamma_xph17b = vert_mul * np.array(f17b["/gamma_func/RE_XPH"])
imgamma_xph17b = vert_mul * np.array(f17b["/gamma_func/IM_XPH"])

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

regamma_pp19b = vert_mul * np.array(f19b["/gamma_func/RE_PP"])
imgamma_pp19b = vert_mul * np.array(f19b["/gamma_func/IM_PP"])
regamma_ph19b = vert_mul * np.array(f19b["/gamma_func/RE_PH"])
imgamma_ph19b = vert_mul * np.array(f19b["/gamma_func/IM_PH"])
regamma_xph19b = vert_mul * np.array(f19b["/gamma_func/RE_XPH"])
imgamma_xph19b = vert_mul * np.array(f19b["/gamma_func/IM_XPH"])

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

regamma_pp21b = vert_mul * np.array(f21b["/gamma_func/RE_PP"])
imgamma_pp21b = vert_mul * np.array(f21b["/gamma_func/IM_PP"])
regamma_ph21b = vert_mul * np.array(f21b["/gamma_func/RE_PH"])
imgamma_ph21b = vert_mul * np.array(f21b["/gamma_func/IM_PH"])
regamma_xph21b = vert_mul * np.array(f21b["/gamma_func/RE_XPH"])
imgamma_xph21b = vert_mul * np.array(f21b["/gamma_func/IM_XPH"])

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

#-------------------------------------GAMMA PLOTTING DIFFERENCES NOCORRECTIONS MINUS CORRECTION RESULTS ------------------------------------------

print("Plotting delta gamma ...")

#---------------------------------------------------------------------------------------------
#
#                   METHOD 2
#
#---------------------------------------------------------------------------------------------

#--- INVR = 120 
diff_regamma_pp1 = vert_mul * np.subtract( f18["/gamma_func/RE_PP"], f17["/gamma_func/RE_PP"])
diff_imgamma_pp1 = vert_mul * np.subtract( f18["/gamma_func/IM_PP"], f17["/gamma_func/IM_PP"] )
diff_regamma_ph1 = vert_mul * np.subtract( f18["/gamma_func/RE_PH"], f17["/gamma_func/RE_PH"] )
diff_imgamma_ph1 = vert_mul * np.subtract( f18["/gamma_func/IM_PH"], f17["/gamma_func/IM_PH"] )
diff_regamma_xph1 = vert_mul * np.subtract(f18["/gamma_func/RE_XPH"],f17["/gamma_func/RE_XPH"])
diff_imgamma_xph1 = vert_mul * np.subtract(f18["/gamma_func/IM_XPH"],f17["/gamma_func/IM_XPH"])

#--- INVR = 100 
diff_regamma_pp2 = vert_mul * np.subtract( f20["/gamma_func/RE_PP"], f19["/gamma_func/RE_PP"])
diff_imgamma_pp2 = vert_mul * np.subtract( f20["/gamma_func/IM_PP"], f19["/gamma_func/IM_PP"] )
diff_regamma_ph2 = vert_mul * np.subtract( f20["/gamma_func/RE_PH"], f19["/gamma_func/RE_PH"] )
diff_imgamma_ph2 = vert_mul * np.subtract( f20["/gamma_func/IM_PH"], f19["/gamma_func/IM_PH"] )
diff_regamma_xph2 = vert_mul * np.subtract(f20["/gamma_func/RE_XPH"],f19["/gamma_func/RE_XPH"])
diff_imgamma_xph2 = vert_mul * np.subtract(f20["/gamma_func/IM_XPH"],f19["/gamma_func/IM_XPH"])

#--- INVR = 80
diff_regamma_pp3 = vert_mul * np.subtract( f22["/gamma_func/RE_PP"], f21["/gamma_func/RE_PP"])
diff_imgamma_pp3 = vert_mul * np.subtract( f22["/gamma_func/IM_PP"], f21["/gamma_func/IM_PP"] )
diff_regamma_ph3 = vert_mul * np.subtract( f22["/gamma_func/RE_PH"], f21["/gamma_func/RE_PH"] )
diff_imgamma_ph3 = vert_mul * np.subtract( f22["/gamma_func/IM_PH"], f21["/gamma_func/IM_PH"] )
diff_regamma_xph3 = vert_mul * np.subtract(f22["/gamma_func/RE_XPH"],f21["/gamma_func/RE_XPH"])
diff_imgamma_xph3 = vert_mul * np.subtract(f22["/gamma_func/IM_XPH"],f21["/gamma_func/IM_XPH"])

#--- INVR = 60 
diff_regamma_pp4 = vert_mul * np.subtract( f2["/gamma_func/RE_PP"], f1["/gamma_func/RE_PP"])
diff_imgamma_pp4 = vert_mul * np.subtract( f2["/gamma_func/IM_PP"], f1["/gamma_func/IM_PP"] )
diff_regamma_ph4 = vert_mul * np.subtract( f2["/gamma_func/RE_PH"], f1["/gamma_func/RE_PH"] )
diff_imgamma_ph4 = vert_mul * np.subtract( f2["/gamma_func/IM_PH"], f1["/gamma_func/IM_PH"] )
diff_regamma_xph4 = vert_mul * np.subtract(f2["/gamma_func/RE_XPH"],f1["/gamma_func/RE_XPH"])
diff_imgamma_xph4 = vert_mul * np.subtract(f2["/gamma_func/IM_XPH"],f1["/gamma_func/IM_XPH"])

#--- INVR = 55
diff_regamma_pp5 = vert_mul * np.subtract( f4["/gamma_func/RE_PP"], f3["/gamma_func/RE_PP"])
diff_imgamma_pp5 = vert_mul * np.subtract( f4["/gamma_func/IM_PP"], f3["/gamma_func/IM_PP"] )
diff_regamma_ph5 = vert_mul * np.subtract( f4["/gamma_func/RE_PH"], f3["/gamma_func/RE_PH"] )
diff_imgamma_ph5 = vert_mul * np.subtract( f4["/gamma_func/IM_PH"], f3["/gamma_func/IM_PH"] )
diff_regamma_xph5 = vert_mul * np.subtract(f4["/gamma_func/RE_XPH"],f3["/gamma_func/RE_XPH"])
diff_imgamma_xph5 = vert_mul * np.subtract(f4["/gamma_func/IM_XPH"],f3["/gamma_func/IM_XPH"])

#--- INVR= 50
diff_regamma_pp6 = vert_mul * np.subtract( f6["/gamma_func/RE_PP"], f5["/gamma_func/RE_PP"])
diff_imgamma_pp6 = vert_mul * np.subtract( f6["/gamma_func/IM_PP"], f5["/gamma_func/IM_PP"] )
diff_regamma_ph6 = vert_mul * np.subtract( f6["/gamma_func/RE_PH"], f5["/gamma_func/RE_PH"] )
diff_imgamma_ph6 = vert_mul * np.subtract( f6["/gamma_func/IM_PH"], f5["/gamma_func/IM_PH"] )
diff_regamma_xph6 = vert_mul * np.subtract(f6["/gamma_func/RE_XPH"],f5["/gamma_func/RE_XPH"])
diff_imgamma_xph6 = vert_mul * np.subtract(f6["/gamma_func/IM_XPH"],f5["/gamma_func/IM_XPH"])

#--- INVR= 45
diff_regamma_pp7 = vert_mul * np.subtract( f8["/gamma_func/RE_PP"], f7["/gamma_func/RE_PP"])
diff_imgamma_pp7 = vert_mul * np.subtract( f8["/gamma_func/IM_PP"], f7["/gamma_func/IM_PP"] )
diff_regamma_ph7 = vert_mul * np.subtract( f8["/gamma_func/RE_PH"], f7["/gamma_func/RE_PH"] )
diff_imgamma_ph7 = vert_mul * np.subtract( f8["/gamma_func/IM_PH"], f7["/gamma_func/IM_PH"] )
diff_regamma_xph7 = vert_mul * np.subtract(f8["/gamma_func/RE_XPH"],f7["/gamma_func/RE_XPH"])
diff_imgamma_xph7 = vert_mul * np.subtract(f8["/gamma_func/IM_XPH"],f7["/gamma_func/IM_XPH"])

#--- INVR=40
diff_regamma_pp8 = vert_mul * np.subtract( f10["/gamma_func/RE_PP"], f9["/gamma_func/RE_PP"])
diff_imgamma_pp8 = vert_mul * np.subtract( f10["/gamma_func/IM_PP"], f9["/gamma_func/IM_PP"] )
diff_regamma_ph8 = vert_mul * np.subtract( f10["/gamma_func/RE_PH"], f9["/gamma_func/RE_PH"] )
diff_imgamma_ph8 = vert_mul * np.subtract( f10["/gamma_func/IM_PH"], f9["/gamma_func/IM_PH"] )
diff_regamma_xph8 = vert_mul * np.subtract(f10["/gamma_func/RE_XPH"],f9["/gamma_func/RE_XPH"])
diff_imgamma_xph8 = vert_mul * np.subtract(f10["/gamma_func/IM_XPH"],f9["/gamma_func/IM_XPH"])

#--- INVR=30
diff_regamma_pp9 = vert_mul * np.subtract( f12["/gamma_func/RE_PP"], f11["/gamma_func/RE_PP"])
diff_imgamma_pp9 = vert_mul * np.subtract( f12["/gamma_func/IM_PP"], f11["/gamma_func/IM_PP"] )
diff_regamma_ph9 = vert_mul * np.subtract( f12["/gamma_func/RE_PH"], f11["/gamma_func/RE_PH"] )
diff_imgamma_ph9 = vert_mul * np.subtract( f12["/gamma_func/IM_PH"], f11["/gamma_func/IM_PH"] )
diff_regamma_xph9 = vert_mul * np.subtract(f12["/gamma_func/RE_XPH"],f11["/gamma_func/RE_XPH"])
diff_imgamma_xph9 = vert_mul * np.subtract(f12["/gamma_func/IM_XPH"],f11["/gamma_func/IM_XPH"])

#--- INVR=20
diff_regamma_pp10 = vert_mul * np.subtract( f14["/gamma_func/RE_PP"], f13["/gamma_func/RE_PP"])
diff_imgamma_pp10 = vert_mul * np.subtract( f14["/gamma_func/IM_PP"], f13["/gamma_func/IM_PP"] )
diff_regamma_ph10 = vert_mul * np.subtract( f14["/gamma_func/RE_PH"], f13["/gamma_func/RE_PH"] )
diff_imgamma_ph10 = vert_mul * np.subtract( f14["/gamma_func/IM_PH"], f13["/gamma_func/IM_PH"] )
diff_regamma_xph10 = vert_mul * np.subtract(f14["/gamma_func/RE_XPH"],f13["/gamma_func/RE_XPH"])
diff_imgamma_xph10 = vert_mul * np.subtract(f14["/gamma_func/IM_XPH"],f13["/gamma_func/IM_XPH"])

#--- INVR=10
diff_regamma_pp11 = vert_mul * np.subtract( f16["/gamma_func/RE_PP"], f15["/gamma_func/RE_PP"])
diff_imgamma_pp11 = vert_mul * np.subtract( f16["/gamma_func/IM_PP"], f15["/gamma_func/IM_PP"] )
diff_regamma_ph11 = vert_mul * np.subtract( f16["/gamma_func/RE_PH"], f15["/gamma_func/RE_PH"] )
diff_imgamma_ph11 = vert_mul * np.subtract( f16["/gamma_func/IM_PH"], f15["/gamma_func/IM_PH"] )
diff_regamma_xph11 = vert_mul * np.subtract(f16["/gamma_func/RE_XPH"],f15["/gamma_func/RE_XPH"])
diff_imgamma_xph11 = vert_mul * np.subtract(f16["/gamma_func/IM_XPH"],f15["/gamma_func/IM_XPH"])

#---------------------------------------------------------------------------------------------
#
#                   METHOD 1
#
#---------------------------------------------------------------------------------------------

#--- INVR = 120 
diff_regamma_pp1b = vert_mul * np.subtract( f18["/gamma_func/RE_PP"], f17b["/gamma_func/RE_PP"])
diff_imgamma_pp1b = vert_mul * np.subtract( f18["/gamma_func/IM_PP"], f17b["/gamma_func/IM_PP"] )
diff_regamma_ph1b = vert_mul * np.subtract( f18["/gamma_func/RE_PH"], f17b["/gamma_func/RE_PH"] )
diff_imgamma_ph1b = vert_mul * np.subtract( f18["/gamma_func/IM_PH"], f17b["/gamma_func/IM_PH"] )
diff_regamma_xph1b = vert_mul * np.subtract(f18["/gamma_func/RE_XPH"],f17b["/gamma_func/RE_XPH"])
diff_imgamma_xph1b = vert_mul * np.subtract(f18["/gamma_func/IM_XPH"],f17b["/gamma_func/IM_XPH"])

#--- INVR = 100 
diff_regamma_pp2b = vert_mul * np.subtract( f20["/gamma_func/RE_PP"], f19b["/gamma_func/RE_PP"])
diff_imgamma_pp2b = vert_mul * np.subtract( f20["/gamma_func/IM_PP"], f19b["/gamma_func/IM_PP"] )
diff_regamma_ph2b = vert_mul * np.subtract( f20["/gamma_func/RE_PH"], f19b["/gamma_func/RE_PH"] )
diff_imgamma_ph2b = vert_mul * np.subtract( f20["/gamma_func/IM_PH"], f19b["/gamma_func/IM_PH"] )
diff_regamma_xph2b = vert_mul * np.subtract(f20["/gamma_func/RE_XPH"],f19b["/gamma_func/RE_XPH"])
diff_imgamma_xph2b = vert_mul * np.subtract(f20["/gamma_func/IM_XPH"],f19b["/gamma_func/IM_XPH"])

#--- INVR = 80
diff_regamma_pp3b = vert_mul * np.subtract( f22["/gamma_func/RE_PP"], f21b["/gamma_func/RE_PP"])
diff_imgamma_pp3b = vert_mul * np.subtract( f22["/gamma_func/IM_PP"], f21b["/gamma_func/IM_PP"] )
diff_regamma_ph3b = vert_mul * np.subtract( f22["/gamma_func/RE_PH"], f21b["/gamma_func/RE_PH"] )
diff_imgamma_ph3b = vert_mul * np.subtract( f22["/gamma_func/IM_PH"], f21b["/gamma_func/IM_PH"] )
diff_regamma_xph3b = vert_mul * np.subtract(f22["/gamma_func/RE_XPH"],f21b["/gamma_func/RE_XPH"])
diff_imgamma_xph3b = vert_mul * np.subtract(f22["/gamma_func/IM_XPH"],f21b["/gamma_func/IM_XPH"])

#--- INVR = 60 
diff_regamma_pp4b = vert_mul * np.subtract( f2["/gamma_func/RE_PP"], f1b["/gamma_func/RE_PP"])
diff_imgamma_pp4b = vert_mul * np.subtract( f2["/gamma_func/IM_PP"], f1b["/gamma_func/IM_PP"] )
diff_regamma_ph4b = vert_mul * np.subtract( f2["/gamma_func/RE_PH"], f1b["/gamma_func/RE_PH"] )
diff_imgamma_ph4b = vert_mul * np.subtract( f2["/gamma_func/IM_PH"], f1b["/gamma_func/IM_PH"] )
diff_regamma_xph4b = vert_mul * np.subtract(f2["/gamma_func/RE_XPH"],f1b["/gamma_func/RE_XPH"])
diff_imgamma_xph4b = vert_mul * np.subtract(f2["/gamma_func/IM_XPH"],f1b["/gamma_func/IM_XPH"])

#--- INVR = 55
diff_regamma_pp5b = vert_mul * np.subtract( f4["/gamma_func/RE_PP"], f3b["/gamma_func/RE_PP"])
diff_imgamma_pp5b = vert_mul * np.subtract( f4["/gamma_func/IM_PP"], f3b["/gamma_func/IM_PP"] )
diff_regamma_ph5b = vert_mul * np.subtract( f4["/gamma_func/RE_PH"], f3b["/gamma_func/RE_PH"] )
diff_imgamma_ph5b = vert_mul * np.subtract( f4["/gamma_func/IM_PH"], f3b["/gamma_func/IM_PH"] )
diff_regamma_xph5b = vert_mul * np.subtract(f4["/gamma_func/RE_XPH"],f3b["/gamma_func/RE_XPH"])
diff_imgamma_xph5b = vert_mul * np.subtract(f4["/gamma_func/IM_XPH"],f3b["/gamma_func/IM_XPH"])

#--- INVR= 50
diff_regamma_pp6b = vert_mul * np.subtract( f6["/gamma_func/RE_PP"], f5b["/gamma_func/RE_PP"])
diff_imgamma_pp6b = vert_mul * np.subtract( f6["/gamma_func/IM_PP"], f5b["/gamma_func/IM_PP"] )
diff_regamma_ph6b = vert_mul * np.subtract( f6["/gamma_func/RE_PH"], f5b["/gamma_func/RE_PH"] )
diff_imgamma_ph6b = vert_mul * np.subtract( f6["/gamma_func/IM_PH"], f5b["/gamma_func/IM_PH"] )
diff_regamma_xph6b = vert_mul * np.subtract(f6["/gamma_func/RE_XPH"],f5b["/gamma_func/RE_XPH"])
diff_imgamma_xph6b = vert_mul * np.subtract(f6["/gamma_func/IM_XPH"],f5b["/gamma_func/IM_XPH"])

#--- INVR= 45
diff_regamma_pp7b = vert_mul * np.subtract( f8["/gamma_func/RE_PP"], f7b["/gamma_func/RE_PP"])
diff_imgamma_pp7b = vert_mul * np.subtract( f8["/gamma_func/IM_PP"], f7b["/gamma_func/IM_PP"] )
diff_regamma_ph7b = vert_mul * np.subtract( f8["/gamma_func/RE_PH"], f7b["/gamma_func/RE_PH"] )
diff_imgamma_ph7b = vert_mul * np.subtract( f8["/gamma_func/IM_PH"], f7b["/gamma_func/IM_PH"] )
diff_regamma_xph7b = vert_mul * np.subtract(f8["/gamma_func/RE_XPH"],f7b["/gamma_func/RE_XPH"])
diff_imgamma_xph7b = vert_mul * np.subtract(f8["/gamma_func/IM_XPH"],f7b["/gamma_func/IM_XPH"])

#--- INVR=40
diff_regamma_pp8b = vert_mul * np.subtract( f10["/gamma_func/RE_PP"], f9b["/gamma_func/RE_PP"])
diff_imgamma_pp8b = vert_mul * np.subtract( f10["/gamma_func/IM_PP"], f9b["/gamma_func/IM_PP"] )
diff_regamma_ph8b = vert_mul * np.subtract( f10["/gamma_func/RE_PH"], f9b["/gamma_func/RE_PH"] )
diff_imgamma_ph8b = vert_mul * np.subtract( f10["/gamma_func/IM_PH"], f9b["/gamma_func/IM_PH"] )
diff_regamma_xph8b = vert_mul * np.subtract(f10["/gamma_func/RE_XPH"],f9b["/gamma_func/RE_XPH"])
diff_imgamma_xph8b = vert_mul * np.subtract(f10["/gamma_func/IM_XPH"],f9b["/gamma_func/IM_XPH"])

#--- INVR=30
diff_regamma_pp9b = vert_mul * np.subtract( f12["/gamma_func/RE_PP"], f11b["/gamma_func/RE_PP"])
diff_imgamma_pp9b = vert_mul * np.subtract( f12["/gamma_func/IM_PP"], f11b["/gamma_func/IM_PP"] )
diff_regamma_ph9b = vert_mul * np.subtract( f12["/gamma_func/RE_PH"], f11b["/gamma_func/RE_PH"] )
diff_imgamma_ph9b = vert_mul * np.subtract( f12["/gamma_func/IM_PH"], f11b["/gamma_func/IM_PH"] )
diff_regamma_xph9b = vert_mul * np.subtract(f12["/gamma_func/RE_XPH"],f11b["/gamma_func/RE_XPH"])
diff_imgamma_xph9b = vert_mul * np.subtract(f12["/gamma_func/IM_XPH"],f11b["/gamma_func/IM_XPH"])

#--- INVR=20
diff_regamma_pp10b = vert_mul * np.subtract( f14["/gamma_func/RE_PP"], f13b["/gamma_func/RE_PP"])
diff_imgamma_pp10b = vert_mul * np.subtract( f14["/gamma_func/IM_PP"], f13b["/gamma_func/IM_PP"] )
diff_regamma_ph10b = vert_mul * np.subtract( f14["/gamma_func/RE_PH"], f13b["/gamma_func/RE_PH"] )
diff_imgamma_ph10b = vert_mul * np.subtract( f14["/gamma_func/IM_PH"], f13b["/gamma_func/IM_PH"] )
diff_regamma_xph10b = vert_mul * np.subtract(f14["/gamma_func/RE_XPH"],f13b["/gamma_func/RE_XPH"])
diff_imgamma_xph10b = vert_mul * np.subtract(f14["/gamma_func/IM_XPH"],f13b["/gamma_func/IM_XPH"])

#--- INVR=10
diff_regamma_pp11 = vert_mul * np.subtract( f16["/gamma_func/RE_PP"], f15["/gamma_func/RE_PP"])
diff_imgamma_pp11 = vert_mul * np.subtract( f16["/gamma_func/IM_PP"], f15["/gamma_func/IM_PP"] )
diff_regamma_ph11 = vert_mul * np.subtract( f16["/gamma_func/RE_PH"], f15["/gamma_func/RE_PH"] )
diff_imgamma_ph11 = vert_mul * np.subtract( f16["/gamma_func/IM_PH"], f15["/gamma_func/IM_PH"] )
diff_regamma_xph11 = vert_mul * np.subtract(f16["/gamma_func/RE_XPH"],f15["/gamma_func/RE_XPH"])
diff_imgamma_xph11 = vert_mul * np.subtract(f16["/gamma_func/IM_XPH"],f15["/gamma_func/IM_XPH"])

#---  Helper functions (include translation from  nambu to physical CAUTION : USING TIME REVERSAL SYMM) 
def plotgamma_axis( use_pl, arr1, arr2, arr3, arr4, arr5,arr6, string, legend ):
    zarr1 = np.array([ arr1[0,n+fdimo217-fdim_min/2,fdimo217,0,0,0,0,0,0,0] for n in range(fdim_min)])
    zarr2 = np.array([ arr2[0,n+fdimo219-fdim_min/2,fdimo219,0,0,0,0,0,0,0] for n in range(fdim_min)])
    zarr3 = np.array([ arr3[0,n+fdimo221-fdim_min/2,fdimo221,0,0,0,0,0,0,0] for n in range(fdim_min)])
    zarr4 = np.array([ arr4[shift + (bdim1-1)/2,n+fdimo21-fdim_min/2,fdimo21,0,0,0,0,0,0,0] for n in range(fdim_min)])
    zarr5 = np.array([ arr5[shift + (bdim9-1)/2,n+fdimo29-fdim_min/2,fdimo29,0,0,0,0,0,0,0] for n in range(fdim_min)])
    zarr6 = np.array([ arr6[shift + (bdim13-1)/2,n+fdimo213-fdim_min/2,fdimo213,0,0,0,0,0,0,0] for n in range(fdim_min)])
    pl.plot( gammagrid_plot, zarr1, marker = 'None', linestyle='-', color='#000000',  label=r"N$_{\text{inv}}=240$")
    print string
    print "240 ... 40", zarr1[0],zarr2[0],zarr3[0],zarr4[0],zarr5[0],zarr6[0]
    pl.plot( gammagrid_plot, zarr2, marker = 'None', linestyle='-', color='#00008B',  label=r"N$_{\text{inv}}=200$")
    pl.plot( gammagrid_plot, zarr3, marker = 'None', linestyle='-', color='#1E90FF',  label=r"N$_{\text{inv}}=160$")
    pl.plot( gammagrid_plot, zarr4, marker = 'None', linestyle='-', color='#32CD32',  label=r"N$_{\text{inv}}=120$")
    pl.plot( gammagrid_plot, zarr5, marker = 'None', linestyle='-', color='#FF0000',  label=r"N$_{\text{inv}}=80$")
    pl.plot( gammagrid_plot, zarr6, marker = 'None', linestyle='-', color='#FFA500',  label=r"N$_{\text{inv}}=40$")
    pl.axhline(y=0, linewidth=2,linestyle = '--', color = 'k')
    pl.xlim([min(gammagrid_plot),max(gammagrid_plot)])
    #ax.set(xlim=[0.8*min(gammagrid_plot),0.8*max(gammagrid_plot)], ylim=[zarr6[0]-0.01,0.01], aspect=1)
    use_pl.set_title( string , fontsize=13 )
    if (legend):
        pl.legend(prop={'size':9})

    return

def plotgamma_diag1( use_pl, arr1, arr2, arr3, arr4, arr5, arr6, string, legend ):
    zarr1 = np.array([ arr1[0,n+fdimo217-fdim_min/2,n+fdimo217-fdim_min/2,0,0,0,0,0,0,0] for n in range(fdim_min)])
    zarr2 = np.array([ arr2[0,n+fdimo219-fdim_min/2,n+fdimo219-fdim_min/2,0,0,0,0,0,0,0] for n in range(fdim_min)])
    zarr3 = np.array([ arr3[0,n+fdimo221-fdim_min/2,n+fdimo221-fdim_min/2,0,0,0,0,0,0,0] for n in range(fdim_min)])
    zarr4 = np.array([ arr4[shift + (bdim1-1)/2,n+fdimo21-fdim_min/2,n+fdimo21-fdim_min/2,0,0,0,0,0,0,0] for n in range(fdim_min)])
    zarr5 = np.array([ arr5[shift + (bdim9-1)/2,n+fdimo29-fdim_min/2,n+fdimo29-fdim_min/2,0,0,0,0,0,0,0] for n in range(fdim_min)])
    zarr6 = np.array([ arr6[shift + (bdim10-1)/2,n+fdimo213-fdim_min/2,n+fdimo213-fdim_min/2,0,0,0,0,0,0,0] for n in range(fdim_min)])
    print string
    print "240 ... 40", zarr1[0],zarr2[0],zarr3[0],zarr4[0],zarr5[0],zarr6[0]
    pl.plot( gammagrid_plot, zarr1, marker = 'None', linestyle='-', color='#000000',  label=r"N$_{\text{inv}}=240$")
    pl.plot( gammagrid_plot, zarr2, marker = 'None', linestyle='-', color='#00008B',  label=r"N$_{\text{inv}}=200$")
    pl.plot( gammagrid_plot, zarr3, marker = 'None', linestyle='-', color='#1E90FF',  label=r"N$_{\text{inv}}=160$")
    pl.plot( gammagrid_plot, zarr4, marker = 'None', linestyle='-', color='#32CD32',  label=r"N$_{\text{inv}}=120$")
    pl.plot( gammagrid_plot, zarr5, marker = 'None', linestyle='-', color='#FF0000',  label=r"N$_{\text{inv}}=80$")
    pl.plot( gammagrid_plot, zarr6, marker = 'None', linestyle='-', color='#FFA500',  label=r"N$_{\text{inv}}=40$")
    pl.axhline(y=0, linewidth=2,linestyle = '--', color = 'k')
    pl.xlim([min(gammagrid_plot),max(gammagrid_plot)])
    #ax.set(xlim=[0.8*min(gammagrid_plot),0.8*max(gammagrid_plot)], ylim=[zarr6[0]-0.01,0.01], aspect=1)
    use_pl.set_title( string , fontsize=13 )
    if(legend):
        pl.legend(prop={'size':9})

    return

#=============================================================================
#
#                   METHOD_2
#=============================================================================
fig, ax = pl.subplots()

pl.figure(figsize=(12, 5))

#-----Print for table 

plotgamma_axis( pl.subplot(1,3,1), diff_regamma_pp1 + diff_regamma_pp1[:,:,::-1,:,:,:,:,:,:,:], diff_regamma_pp2 + diff_regamma_pp2[:,:,::-1,:,:,:,:,:,:,:], diff_regamma_pp3 + diff_regamma_pp3[:,:,::-1,:,:,:,:,:,:,:], diff_regamma_pp4 + diff_regamma_pp4[:,:,::-1,:,:,:,:,:,:,:],diff_regamma_pp8 + diff_regamma_pp8[:,:,::-1,:,:,:,:,:,:,:],diff_regamma_pp10 + diff_regamma_pp10[:,:,::-1,:,:,:,:,:,:,:],RE + r"\delta \Gamma^{\nu, \nu'=\pi/\beta, \omega =40 \frac{\pi}{\beta}}_{s}$", True ) # flip sign of w_out
pl.xlabel(r"$\nu$", fontsize = 12)
plotgamma_axis( pl.subplot(1,3,2), 2*diff_regamma_ph1 - diff_regamma_xph1, 2*diff_regamma_ph2 - diff_regamma_xph2, 2*diff_regamma_ph3 - diff_regamma_xph3, 2*diff_regamma_ph4- diff_regamma_xph4,2*diff_regamma_ph8- diff_regamma_xph8,2*diff_regamma_ph10 - diff_regamma_xph10, RE + r"\delta \Gamma^{\nu, \nu'=\pi/\beta, \omega =40 \frac{\pi}{\beta}}_{d}$", False )
pl.xlabel(r"$\nu$", fontsize = 12)
plotgamma_axis( pl.subplot(1,3,3), - diff_regamma_xph1, - diff_regamma_xph2, - diff_regamma_xph3, - diff_regamma_xph4,- diff_regamma_xph8,- diff_regamma_xph10, RE + r"\delta \Gamma^{\nu, \nu'=\pi/\beta, \omega =40 \frac{\pi}{\beta}}_{m}$", False )
pl.xlabel(r"$\nu$", fontsize = 12)

pl.tight_layout()

#--- Save to file
pl.savefig("plots/DeltaGammas_METHOD2_axis_nu'=0_W20.png", dpi = 150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()

pl.figure(figsize=(12, 5))

plotgamma_axis( pl.subplot(1,3,1), diff_regamma_pp1 + diff_regamma_pp1[:,:,::-1,:,:,:,:,:,:,:], diff_regamma_pp2 + diff_regamma_pp2[:,:,::-1,:,:,:,:,:,:,:], diff_regamma_pp3 + diff_regamma_pp3[:,:,::-1,:,:,:,:,:,:,:], diff_regamma_pp4 + diff_regamma_pp4[:,:,::-1,:,:,:,:,:,:,:],diff_regamma_pp8 + diff_regamma_pp8[:,:,::-1,:,:,:,:,:,:,:],diff_regamma_pp10 + diff_regamma_pp10[:,:,::-1,:,:,:,:,:,:,:],RE + r"\delta \Gamma^{\nu=\nu'=, \omega =40 \frac{\pi}{\beta}}_{s}$", True ) # flip sign of w_out
pl.xlabel(r"$\nu$", fontsize = 12)
plotgamma_axis( pl.subplot(1,3,2), 2*diff_regamma_ph1 - diff_regamma_xph1, 2*diff_regamma_ph2 - diff_regamma_xph2, 2*diff_regamma_ph3 - diff_regamma_xph3, 2*diff_regamma_ph4 - diff_regamma_xph4,2*diff_regamma_ph8- diff_regamma_xph8,2*diff_regamma_ph10 - diff_regamma_xph10, RE + r"\delta \Gamma^{\nu=\nu', \omega =40 \frac{\pi}{\beta}}_{d}$", False )
pl.xlabel(r"$\nu$", fontsize = 12)
plotgamma_axis( pl.subplot(1,3,3), - diff_regamma_xph1, - diff_regamma_xph2, - diff_regamma_xph3, - diff_regamma_xph4,- diff_regamma_xph8,- diff_regamma_xph10, RE + r"\delta \Gamma^{\nu=\nu', \omega =40 \frac{\pi}{\beta}}_{m}$", False )
pl.xlabel(r"$\nu$", fontsize = 12)
pl.tight_layout()

#--- Save to file
pl.savefig("plots/DeltaGammas_METHOD2_nu=nu'_W20.png", dpi = 150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()

pl.figure(figsize=(12, 5))
#=============================================================================
#
#                   METHOD_1
#=============================================================================
#--- Plot along axis \nu'=0
plotgamma_axis( pl.subplot(1,3,1), diff_regamma_pp1b + diff_regamma_pp1b[:,:,::-1,:,:,:,:,:,:,:], diff_regamma_pp2b + diff_regamma_pp2b[:,:,::-1,:,:,:,:,:,:,:], diff_regamma_pp3b + diff_regamma_pp3b[:,:,::-1,:,:,:,:,:,:,:], diff_regamma_pp4b + diff_regamma_pp4b[:,:,::-1,:,:,:,:,:,:,:],diff_regamma_pp8b + diff_regamma_pp8b[:,:,::-1,:,:,:,:,:,:,:],diff_regamma_pp10b + diff_regamma_pp10b[:,:,::-1,:,:,:,:,:,:,:],RE + r"\delta \Gamma^{\nu, \nu'=\pi/\beta, \omega =40 \frac{\pi}{\beta}}_{s}$", True ) # flip sign of w_out
pl.xlabel(r"$\nu$", fontsize = 12)
plotgamma_axis( pl.subplot(1,3,2), 2*diff_regamma_ph1b - diff_regamma_xph1b, 2*diff_regamma_ph2b - diff_regamma_xph2b, 2*diff_regamma_ph3b - diff_regamma_xph3b, 2*diff_regamma_ph4b- diff_regamma_xph4b,2*diff_regamma_ph8b- diff_regamma_xph8b,2*diff_regamma_ph10b - diff_regamma_xph10b, RE + r"\delta \Gamma^{\nu, \nu'=\pi/\beta, \omega =40 \frac{\pi}{\beta}}_{d}$", False )
pl.xlabel(r"$\nu$", fontsize = 12)
plotgamma_axis( pl.subplot(1,3,3), - diff_regamma_xph1b, - diff_regamma_xph2b, - diff_regamma_xph3b, - diff_regamma_xph4b,- diff_regamma_xph8b,- diff_regamma_xph10b, RE + r"\delta \Gamma^{\nu, \nu'=\pi/\beta, \omega =40 \frac{\pi}{\beta}}_{m}$", False )
pl.xlabel(r"$\nu$", fontsize = 12)

pl.tight_layout()

#--- Save to file
pl.savefig("plots/DeltaGammas_METHOD1_axis_nu'=0_W20.png", dpi = 150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()

pl.figure(figsize=(12, 5))

plotgamma_axis( pl.subplot(1,3,1), diff_regamma_pp1b + diff_regamma_pp1b[:,:,::-1,:,:,:,:,:,:,:], diff_regamma_pp2b + diff_regamma_pp2b[:,:,::-1,:,:,:,:,:,:,:], diff_regamma_pp3b + diff_regamma_pp3b[:,:,::-1,:,:,:,:,:,:,:], diff_regamma_pp4b + diff_regamma_pp4b[:,:,::-1,:,:,:,:,:,:,:],diff_regamma_pp8b + diff_regamma_pp8b[:,:,::-1,:,:,:,:,:,:,:],diff_regamma_pp10b + diff_regamma_pp10b[:,:,::-1,:,:,:,:,:,:,:],RE + r"\delta \Gamma^{\nu=\nu'=, \omega =40 \frac{\pi}{\beta}}_{s}$", True ) # flip sign of w_out
pl.xlabel(r"$\nu$", fontsize = 12)
plotgamma_axis( pl.subplot(1,3,2), 2*diff_regamma_ph1b - diff_regamma_xph1b, 2*diff_regamma_ph2b - diff_regamma_xph2b, 2*diff_regamma_ph3b - diff_regamma_xph3b, 2*diff_regamma_ph4b - diff_regamma_xph4b,2*diff_regamma_ph8b- diff_regamma_xph8b,2*diff_regamma_ph10b - diff_regamma_xph10b, RE + r"\delta \Gamma^{\nu=\nu', \omega =40 \frac{\pi}{\beta}}_{d}$", False )
pl.xlabel(r"$\nu$", fontsize = 12)
plotgamma_axis( pl.subplot(1,3,3), - diff_regamma_xph1b, - diff_regamma_xph2b, - diff_regamma_xph3b, - diff_regamma_xph4b,- diff_regamma_xph8b,- diff_regamma_xph10b, RE + r"\delta \Gamma^{\nu=\nu', \omega =40 \frac{\pi}{\beta}}_{m}$", False )
pl.xlabel(r"$\nu$", fontsize = 12)
pl.tight_layout()

#--- Save to file
pl.savefig("plots/DeltaGammas_METHOD1_nu=nu'_W20.png", dpi = 150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()

