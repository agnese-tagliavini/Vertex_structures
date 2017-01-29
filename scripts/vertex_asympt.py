#!/usr/bin/python
#
# This script provides the calculation of the following objects from the ED vertex:
#
# - Plus functions in the pp and in the ph-notations for the upup and the updown cases
# - Karrasch functions in the pp and ph-notations for the upup ans updown cases
# - Three-leg object in the three channels : charge, spin and pairing
# - Susceptibilities in the three channels : charge, spin and pairing
#
# Everything will be stored in the HDF5 FILE which already contains all the info about the ED calculation of the vertex
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
from agneselibrary.mymath import *
from agneselibrary.translate_notation import *
from agneselibrary.inf_Msum import *
from _functools import partial

#-----------------------------------Read parameters----------------------------------------------
# ------U---------

#a=raw_input('Enter the interaction value:') 
#try:
#    U=float(a)
#except ValueError:
#    sys.exit("Invalid interaction")
#
#print ("Uhub value     " + str(U) )
#
##-------BETA-------
#
#a = raw_input('Enter the value of beta:') 
#try:
#    beta=float(a) 
#except ValueError:
#    sys.exit("Invalid beta")
#
#print ("beta value     " + str(beta) )

U = 1.0
beta = 26.0

#----------------------------------------Read HDF5 files-----------------------------------------

if ('../dat'):
    f = h5py.File('../dat/U'+ str(U)+'_beta'+ str(beta)+'_FFREQ_20_BFREQ_20.h5', 'r+')   # Read (and write) the hdf5 file in the directory "dat" if existing
else:
    sys.exit("No data file")

#---------  PARAMETERS  ---------------

sign = -1.0
pi = math.pi
N_iter_max = 20
accuracy = 10**(-12)

#---------------------- Read objects involved in the calculation from HDF5 file -------------------------------

#--------------------------------- READING ----------------------------------------------------------

#GF

re_g_iw = f["Giw/RE"][:]
im_g_iw = f["Giw/IM"][:]           
print im_g_iw.shape[0] 
fgrid_gf = im_g_iw.shape[0]
N_fermi_gf = fgrid_gf/2


#VERTEX

re_f_upup_ph = f["VERT/PH/RE_F_UPUP"][:]
re_f_updo_ph = f["VERT/PH/RE_F_UPDO"][:]
re_f_upup_xph = f["VERT/XPH/RE_F_UPUP"][:]
re_f_updo_xph = f["VERT/XPH/RE_F_UPDO"][:]
re_f_upup_pp = f["VERT/PP/RE_F_UPUP"][:]
re_f_updo_pp = f["VERT/PP/RE_F_UPDO"][:]


im_f_upup_ph = f["VERT/PH/IM_F_UPUP"][:]
im_f_updo_ph = f["VERT/PH/IM_F_UPDO"][:]
im_f_upup_xph = f["VERT/XPH/IM_F_UPUP"][:]
im_f_updo_xph = f["VERT/XPH/IM_F_UPDO"][:]
im_f_upup_pp = f["VERT/PP/IM_F_UPUP"][:]
im_f_updo_pp = f["VERT/PP/IM_F_UPDO"][:]

print re_f_upup_ph.shape
#--------------------------------------- FREQUENCY RANGES AND VARIABLES ------------------------

# We assume all the channels to have the same B/F grids

fgrid = f["VERT/PH/fgrid"][:].shape[0] 
bgrid = f["VERT/PH/bgrid"][:].shape[0]

print fgrid

N_bose = (bgrid-1)/2 # to create a bosonic frequency grid from -N_bose to N_bose
N_fermi= (fgrid)/2   # to create a fermionic frequency grid from -N_fermi to N_fermi

print bgrid

# FREQUENCY RANGES OF THE ASYMPTOTIC STRUCTURES
# THE FFREQUENCY RANGE AND THE SUMMATION RANGE HAS TO BE >> THAN THE BFREQUENCY RANGE IN ORDER TO RESOLVE THE STRUCTURE

#N_bose_big = 2*N_bose

#Karrasch function box:
N_bose_big_k = 2*N_bose
#N_bose_big_k = 100
bgrid_big_k = 2*N_bose_big_k + 1

#Plus function box
N_bose_big_p = 2*N_bose
#N_bose_big_p = 100
bgrid_big_p = 2*N_bose_big_p + 1

N_fermi_big = 2*N_fermi
#N_fermi_big = 90
fgrid_big = 2*(N_fermi_big)

N_fermi_sum = 2*max(N_bose_big_k,N_fermi_big,N_bose_big_p)

print N_fermi_sum, N_bose_big_k, N_bose_big_p, N_fermi_big

wNum = int(N_fermi_sum/100.0*10.0) # 10% of the the range/2 of my finite Matsubara summation
#wNum = 10
iMin = N_fermi_sum-wNum
Nl = 5  #Order of my fitting function (see library inf_Msum)

print wNum

sum_fit = generate_sum_func(iMin,wNum,Nl) #It returns the Matsubara summation function {-infty,infty}
dblsum_fit = generate_dblsum_func(iMin,wNum,Nl) #It returns double Matsubara summation function {-infty,Infty}

#
# Initial condition for asymptotics
#

if (('P_func' in f) or ('K_func' in f)):
    print "P and K already exist!!"
    K_upup_ph_loop  = np.array(f["/K_func/PH/RE_K_UPUP"])
    K_upup_xph_loop = np.array(f["/K_func/XPH/RE_K_UPUP"])
    K_updo_pp_loop  = np.array(f["/K_func/PP/RE_K_UPDO"])
    K_updo_ph_loop  = np.array(f["/K_func/PH/RE_K_UPDO"])
    K_updo_xph_loop = np.array(f["/K_func/XPH/RE_K_UPDO"])
    P_upup_ph_loop    = np.array(f["/P_func/PH/RE_P_UPUP"])
    P_upup_xph_loop   = np.array(f["/P_func/XPH/RE_P_UPUP"])
    P_updo_pp_loop    = np.array(f["/P_func/PP/RE_P_UPDO"])
    P_updo_ph_loop    = np.array(f["/P_func/PH/RE_P_UPDO"])
    P_updo_xph_loop   = np.array(f["/P_func/XPH/RE_P_UPDO"])
    
else:
    print " P and K do not exist yet!!"
    K_upup_ph_loop  = np.zeros ( bgrid_big_k , dtype='float64' )
    K_upup_xph_loop = np.zeros ( bgrid_big_k , dtype='float64' )
    K_updo_pp_loop  = np.zeros ( bgrid_big_k , dtype='float64' )
    K_updo_ph_loop  = np.zeros ( bgrid_big_k , dtype='float64' )
    K_updo_xph_loop = np.zeros ( bgrid_big_k , dtype='float64' )
    P_upup_ph_loop    = np.zeros ( (fgrid_big, bgrid_big_p) , dtype='float64' )
    P_upup_xph_loop   = np.zeros ( (fgrid_big,bgrid_big_p) , dtype='float64' )
    P_updo_pp_loop    = np.zeros ( (fgrid_big,bgrid_big_p) , dtype='float64' )
    P_updo_ph_loop    = np.zeros ( (fgrid_big,bgrid_big_p) , dtype='float64' )
    P_updo_xph_loop   = np.zeros ( (fgrid_big,bgrid_big_p) , dtype='float64' )


#-----------GF------------------------------------------------

def G(wf):
    if (wf >= -N_fermi_gf and wf < N_fermi_gf):
        return re_g_iw[wf+N_fermi_gf]+1j*im_g_iw[wf+N_fermi_gf]
    else:
        wMat = 1j*(2*wf +1)*pi/beta     #GF asymtotic behavior
        return 1.0/wMat


# Construction of the inner vertex box

def isInside(i,j,k):
    return abs(i) <= N_bose and j >= -N_fermi and j < N_fermi and k >= -N_fermi and k < N_fermi


# -------------------------------- Main loop--------------------------------------------------
#                                                                                    

for ind in range(N_iter_max):

    print "Inside the loop!"
     
    # Update the old asymptotic structures (from the previous loop)
                                                                                    
    def K_upup_ph_old(wb):
        if (abs(wb) <= N_bose_big_k):
            return K_upup_ph_loop[wb+N_bose_big_k]
        else:
            return 0.0
    def K_updo_ph_old(wb):
        if (abs(wb) <= N_bose_big_k):
            return K_updo_ph_loop[wb+N_bose_big_k]
        else:
            return 0.0
    def K_updo_pp_old(wb):
        if (abs(wb) <= N_bose_big_k):
            return K_updo_pp_loop[wb+N_bose_big_k]
        else:
            return 0.0
    def K_upup_xph_old(wb):
        if (abs(wb) <= N_bose_big_k):
            return K_upup_xph_loop[wb+N_bose_big_k]
        else:
            return 0.0
    def K_updo_xph_old(wb):
        if (abs(wb) <= N_bose_big_k):
            return K_updo_xph_loop[wb+N_bose_big_k]
        else:
            return 0.0
    
    def P_upup_ph_old(wb,wf):
        if (abs(wb) <= N_bose_big_p and wf >= -N_fermi_big and wf < N_fermi_big):
            return P_upup_ph_loop[wf+N_fermi_big,wb+N_bose_big_p]
        else:
            return 0.0
    def P_updo_ph_old(wb,wf):
        if (abs(wb) <= N_bose_big_p and wf >= -N_fermi_big and wf < N_fermi_big):
            return P_updo_ph_loop[wf+N_fermi_big,wb+N_bose_big_p]
        else:
            return 0.0
    def P_updo_pp_old(wb,wf):
        if (abs(wb) <= N_bose_big_p and wf >= -N_fermi_big and wf < N_fermi_big):
            return P_updo_pp_loop[wf+N_fermi_big,wb+N_bose_big_p]
        else:
            return 0.0
    def P_upup_xph_old(wb,wf):
        if (abs(wb) <= N_bose_big_p and wf >= -N_fermi_big and wf < N_fermi_big):
            return P_upup_xph_loop[wf+N_fermi_big,wb+N_bose_big_p]
        else:
            return 0.0
    def P_updo_xph_old(wb,wf):
        if (abs(wb) <= N_bose_big_p and wf >= -N_fermi_big and wf < N_fermi_big):
            return P_updo_xph_loop[wf+N_fermi_big,wb+N_bose_big_p]
        else:
            return 0.0
    
    # Update the asymptotic structures for the VERTEX IN ALL CHANNELS
    
    print "Full vertex with asymptotics"

    def f_upup_fun_ph(i,j,k):
        if isInside(i,j,k):
            return re_f_upup_ph[i + N_bose, j+N_fermi, k + N_fermi]+1j*im_f_upup_ph[i + N_bose, j+N_fermi, k + N_fermi]
        else:
            return K_upup_ph_old(i) + P_upup_ph_old(i,j)+ P_upup_ph_old(i,k) + K_upup_xph_old(PHtoXPH((i,j,k))[0]) + P_upup_xph_old(PHtoXPH((i,j,k))[0],PHtoXPH((i,j,k))[1])+P_upup_xph_old(PHtoXPH((i,j,k))[0],PHtoXPH((i,j,k))[2])
               
    def f_updo_fun_ph(i,j,k):
        if isInside(i,j,k):
            return re_f_updo_ph[i + N_bose, j+N_fermi, k + N_fermi]+1j*im_f_updo_ph[i + N_bose, j+N_fermi, k + N_fermi]
        else:
            return - U + K_updo_ph_old(i) + P_updo_ph_old(i,j)+P_updo_ph_old(i,k) + K_updo_xph_old(PHtoXPH((i,j,k))[0]) + P_updo_xph_old(PHtoXPH((i,j,k))[0],PHtoXPH((i,j,k))[1]) + P_updo_xph_old(PHtoXPH((i,j,k))[0],PHtoXPH((i,j,k))[2])+ K_updo_pp_old(PHtoPP((i,j,k))[0]) + P_updo_pp_old(PHtoPP((i,j,k))[0],PHtoPP((i,j,k))[1]) + P_updo_pp_old(PHtoPP((i,j,k))[0],PHtoPP((i,j,k))[2])
    
    def f_updo_fun_pp(i,j,k):
        if isInside(i,j,k):
            return re_f_updo_pp[i + N_bose, j+N_fermi, k + N_fermi]+1j*im_f_updo_pp[i + N_bose, j+N_fermi, k + N_fermi]
        else:
            return - U + K_updo_pp_old(i) + P_updo_pp_old(i,j) + P_updo_pp_old(i,k) + K_updo_ph_old(PPtoPH((i,j,k))[0]) + P_updo_ph_old(PPtoPH((i,j,k))[0],PPtoPH((i,j,k))[1]) + P_updo_ph_old(PPtoPH((i,j,k))[0],PPtoPH((i,j,k))[2])+ K_updo_xph_old(PPtoXPH((i,j,k))[0]) + P_updo_xph_old(PPtoXPH((i,j,k))[0],PPtoXPH((i,j,k))[1]) + P_updo_xph_old(PPtoXPH((i,j,k))[0],PPtoXPH((i,j,k))[2])
    
    def f_xupdo_fun_pp(i,j,k):
        if isInside(i,j,k):
            return re_f_upup_pp[i + N_bose, j+N_fermi, k + N_fermi]+1j*im_f_upup_pp[i + N_bose, j+N_fermi, k + N_fermi]-re_f_updo_pp[i + N_bose, j+N_fermi, k + N_fermi]-1j*im_f_updo_pp[i + N_bose, j+N_fermi, k + N_fermi]
        else:
            return U - K_updo_pp_old(i) - P_updo_pp_old(i,j) - P_updo_pp_old(i,k) - K_updo_xph_old(PPtoXPH((i,j,-k-mymod_abs(i)-1))[0]) - P_updo_xph_old(PPtoXPH((i,j,-k-mymod_abs(i)-1))[0],PPtoXPH((i,j,-k-mymod_abs(i)-1))[1]) - P_updo_xph_old(PPtoXPH((i,j,-k-mymod_abs(i)-1))[0],PPtoXPH((i,j,-k-mymod_abs(i)-1))[2])- K_updo_ph_old(PPtoPH((i,j,-k-mymod_abs(i)-1))[0]) - P_updo_ph_old(PPtoPH((i,j,-k-mymod_abs(i)-1))[0],PPtoPH((i,j,-k-mymod_abs(i)-1))[1]) - P_updo_ph_old(PPtoPH((i,j,-k-mymod_abs(i)-1))[0],PPtoPH((i,j,-k-mymod_abs(i)-1))[2])    

    def f_upup_fun_xph(i,j,k):
        if isInside(i,j,k):
            return re_f_upup_xph[i + N_bose, j+N_fermi, k + N_fermi]+1j*im_f_upup_xph[i + N_bose, j+N_fermi, k + N_fermi]
        else:
            return K_upup_xph_old(i) + P_upup_xph_old(i,j)+P_upup_xph_old(i,k) + K_upup_ph_old(XPHtoPH((i,j,k))[0]) + P_upup_ph_old(XPHtoPH((i,j,k))[0],XPHtoPH((i,j,k))[1])+ P_upup_ph_old(XPHtoPH((i,j,k))[0],XPHtoPH((i,j,k))[2]) 
    
    def f_updo_fun_xph(i,j,k):
        if isInside(i,j,k):
            return re_f_updo_xph[i + N_bose, j+N_fermi, k + N_fermi]+1j*im_f_updo_xph[i + N_bose, j+N_fermi, k + N_fermi]
        else:
            return - U + K_updo_xph_old(i) + P_updo_xph_old(i,j)+P_updo_xph_old(i,k) + K_updo_ph_old(XPHtoPH((i,j,k))[0]) + P_updo_ph_old(XPHtoPH((i,j,k))[0],XPHtoPH((i,j,k))[1]) + P_updo_ph_old(XPHtoPH((i,j,k))[0],XPHtoPH((i,j,k))[2])+ K_updo_pp_old(XPHtoPP((i,j,k))[0]) + P_updo_pp_old(XPHtoPP((i,j,k))[0],XPHtoPP((i,j,k))[1]) + P_updo_pp_old(XPHtoPP((i,j,k))[0],XPHtoPP((i,j,k))[2])

    # Update of the new asymptotic structures (go to "END")
    
    #-----------------------------------------------------------------------------------------------------------------------------
    #                                           Karrasch functions
    #------------------------------------------------------------------------------------------------------------------------------
    #
    #                                               PH
    #                                                                                    
    #
    #------------------------------------------------------------------------------------------------------------------------------
    print "PH channel diagrams"
    #-------------------------------Karrasch functions and susceptibilities in the ph notation----------------------------------------
    
    def int_chi_0_ph(i,k):
        return -G(k-myfloor_div2(i))*G(k+myceil_div2(i))*1.0/beta   # (-) for internal loop
    
    def int_chi_upup_ph(i,j,k):
        return G(j-myfloor_div2(i))*G(j+myceil_div2(i))*G(k-myfloor_div2(i))*G(k+myceil_div2(i))*f_upup_fun_ph(i,j,k)*1.0/beta*1.0/beta #(-)^2 two internal loops
    
    def int_chi_updo_ph(i,j,k):
        return G(j-myfloor_div2(i))*G(j+myceil_div2(i))*G(k-myfloor_div2(i))*G(k+myceil_div2(i))*f_updo_fun_ph(i,j,k)*1.0/beta*1.0/beta 
    
    # Arrays creation 

    chi_0_ph_arr = np.array([sum_fit(partial(int_chi_0_ph,i)) for i in range(-N_bose_big_k,N_bose_big_k+1)])
    print "first chi ph:" 
    chi_upup_ph_arr = np.array([dblsum_fit(partial(int_chi_upup_ph,i)) for i in range(-N_bose_big_k,N_bose_big_k+1)])
    print "chi_upup ph:"
    chi_updo_ph_arr = np.array([dblsum_fit(partial(int_chi_updo_ph,i)) for i in range(-N_bose_big_k,N_bose_big_k+1)])
    print "chi_updo ph:"
    
    #-------------------------- Karrasch and Plus function arrays PH-------------------------------------------
    
    K_upup_ph_new = U*U*(chi_0_ph_arr + chi_upup_ph_arr)
    print "K_upup_ph:"
    K_updo_ph_new = U*U*chi_updo_ph_arr
    print "K_updo_ph:"
    
    #------------------------------------------------------------------------------------------------------------------------------
    #
    #                                               PP
    #                                                                                    
    #
    #------------------------------------------------------------------------------------------------------------------------------
    print "PP channel diagrams"
    #-------------------------------Karrasch functions and susceptibilities in the ph notation----------------------------------------
    
    def int_chi_0_pp(i,k):
        return G(myfloor_div2(i)-k-1)*G(k+myceil_div2(i))*1.0/beta  # 2 for the spin, 1/2 for the equivalent lines
    
    def int_chi_updo_pp(i,j,k):
        return 0.5*G(myfloor_div2(i)-j-1)*G(j+myceil_div2(i))*G(myfloor_div2(i)-k-1)*G(k+myceil_div2(i))*f_updo_fun_pp(i,j,k)*1.0/beta*1.0/beta # 2 for the spin, 1/4 for 2 pairs of equivalent lines
    
    def int_chi_xupdo_pp(i,j,k):
        return -0.5*G(myfloor_div2(i)-j-1)*G(j+myceil_div2(i))*G(myfloor_div2(i)-k-1)*G(k+myceil_div2(i))*f_xupdo_fun_pp(i,j,k)*1.0/beta*1.0/beta #2 for the spin, 1/4 for 2 pairs of equivalent lines, (-) from the swapping of one of the two bare vertices 
    
    # Arrays creation 

    chi_0_pp_arr = np.array([sum_fit(partial(int_chi_0_pp,i)) for i in range(-N_bose_big_k,N_bose_big_k+1)])
    print "first chi pp:"
    chi_updo_pp_arr = np.array([dblsum_fit(partial(int_chi_updo_pp,i)) for i in range(-N_bose_big_k,N_bose_big_k+1)])
    print "chi_updo pp:"
    chi_xupdo_pp_arr = np.array([dblsum_fit(partial(int_chi_xupdo_pp,i)) for i in range(-N_bose_big_k,N_bose_big_k+1)])
    print "chi_xupdo pp:"
    #-------------------------Definition of Plus function, Karrasch --------------------------
    
    K_updo_pp_new = U*U*(chi_0_pp_arr+chi_updo_pp_arr+chi_xupdo_pp_arr)
    print "K_updo_pp:"
    #------------------------------------------------------------------------------------------------------------------------------
    #
    #                                               XPH
    #                                                                                    
    #
    #------------------------------------------------------------------------------------------------------------------------------
    print "XPH channel diagrams"
    
    #-------------------------------Karrasch functions and susceptibilities in the ph notation----------------------------------------
    
    def int_chi_0_xph(i,k):
        return G(k-myfloor_div2(i))*G(k+myceil_div2(i))*1.0/beta #(-) from loop * (-) from swapping-> UPUP ; no loops no swapping-> UPDO 
    
    def int_chi_upup_xph(i,j,k):
        return -G(j-myfloor_div2(i))*G(j+myceil_div2(i))*G(k-myfloor_div2(i))*G(k+myceil_div2(i))*f_upup_fun_ph(i,k,j)*1.0/(beta*beta) #(-) from swapping* (-)^2 from 2 internal loops 
    
    def int_chi_updo_xph(i,j,k):
        return G(j-myfloor_div2(i))*G(j+myceil_div2(i))*G(k-myfloor_div2(i))*G(k+myceil_div2(i))*f_updo_fun_xph(i,k,j)*1.0/(beta*beta) # no loops no swapping 
    
    # Arrays creation 

    chi_0_xph_arr = np.array([sum_fit(partial(int_chi_0_xph,i)) for i in range(-N_bose_big_k,N_bose_big_k+1)])
    print "first chi_xph:"
    chi_upup_xph_arr = np.array([dblsum_fit(partial(int_chi_upup_xph,i)) for i in range(-N_bose_big_k,N_bose_big_k+1)])
    print "chi_upup xph:"
    chi_updo_xph_arr = np.array([dblsum_fit(partial(int_chi_updo_xph,i)) for i in range(-N_bose_big_k,N_bose_big_k+1)])
    print "chi_updo xph:"
    
    #-------------------------Definition of Plus function, Karrasch PH--------------------------
    
    K_upup_xph_new = U*U*(chi_0_xph_arr + chi_upup_xph_arr) 
    print "K_upup_xph:"
    K_updo_xph_new = U*U*(chi_0_xph_arr + chi_updo_xph_arr)
    print "K_updo_xph:"

    # Compute variance

    epsik = np.max(map(mynorm, [K_upup_ph_new - K_upup_ph_loop, K_updo_ph_new - K_updo_ph_loop, K_updo_pp_new - K_updo_pp_loop, K_upup_xph_new - K_upup_xph_loop, K_updo_xph_new - K_updo_xph_loop]))
    print "Epsi Karraschl:", epsik 

    #--------Update karrasch functions--------
    print "Update Karrasch functions"

    K_upup_ph_loop = K_upup_ph_new
    K_updo_ph_loop = K_updo_ph_new
    K_updo_pp_loop = K_updo_pp_new
    K_upup_xph_loop = K_upup_xph_new
    K_updo_xph_loop = K_updo_xph_new

#    break
    #-----------------------------------------------------------------------------------------------------------------------------
    #                                           Plus functions
    #------------------------------------------------------------------------------------------------------------------------------
    #
    #                                               PH
    #                                                                                    
    #
    #------------------------------------------------------------------------------------------------------------------------------
    print "PH channel diagrams"
    #----------------------------------Plus functions and three-leg vertex in the ph notation---------------------------------------
    
    def int_gamma_3_upup_ph(i,j,k):
         return -G(k-myfloor_div2(i))*G(k+myceil_div2(i))*f_updo_fun_ph(i,j,k)*1.0/beta # (-) from the internal loop
    
    def int_gamma_3_updo_ph(i,j,k):
        return -G(k-myfloor_div2(i))*G(k+myceil_div2(i))*f_upup_fun_ph(i,j,k)*1.0/beta # (-) from the internal loop
    
    # Arrays creation

    gamma_3_upup_ph_arr = np.array([[ sum_fit(partial(int_gamma_3_upup_ph,i,j)) for i in range(-N_bose_big_p,N_bose_big_p+1)] for j in range(-N_fermi_big, N_fermi_big)])
    gamma_3_updo_ph_arr = np.array([[ sum_fit(partial(int_gamma_3_updo_ph,i,j)) for i in range(-N_bose_big_p,N_bose_big_p+1)] for j in range(-N_fermi_big, N_fermi_big)])
    print "gammas3 ph"
    
    #-------------------------Definition of Plus function--------------------------------------------------------------------------

    P_upup_ph_new = (sign*U)*gamma_3_upup_ph_arr - K_upup_ph_new[N_bose_big_k - N_bose_big_p:N_bose_big_k + N_bose_big_p+1]
    print "Plus upup ph"
    P_updo_ph_new =  (sign*U)*gamma_3_updo_ph_arr - K_updo_ph_new[N_bose_big_k - N_bose_big_p:N_bose_big_k + N_bose_big_p+1]
    print "Plus updo ph"   
    
    #------------------------------------------------------------------------------------------------------------------------------
    #
    #                                               PP
    #                                                                                    
    #
    #------------------------------------------------------------------------------------------------------------------------------
    print "PP channel diagrams"
    #----------------------------------Plus functions and three-leg vertex in the pp notation---------------------------------------
    
    def int_gamma_3_updo_pp(i,j,k):
        return 0.5*G(myfloor_div2(i)-k-1)*G(k+myceil_div2(i))*f_updo_fun_pp(i,j,k)*1.0/beta # 1/2 from two equivalent lines

    def int_gamma_3_xupdo_pp(i,j,k):
        return -0.5*G(myfloor_div2(i)-k-1)*G(k+myceil_div2(i))*f_xupdo_fun_pp(i,j,k)*1.0/beta # 1/2 from two equivalent lines and (-) sign from the swapping of two lines for the bare vertex
    
    # Arrays creation

    gamma_3_updo_pp_arr = np.array([[ sum_fit(partial(int_gamma_3_updo_pp,i,j)) for i in range(-N_bose_big_p,N_bose_big_p+1)] for j in range(-N_fermi_big, N_fermi_big)])
    gamma_3_xupdo_pp_arr = np.array([[ sum_fit(partial(int_gamma_3_xupdo_pp,i,j)) for i in range(-N_bose_big_p,N_bose_big_p+1)] for j in range(-N_fermi_big, N_fermi_big)])
    print "gammas3 pp"
    
    #-------------------------Definition of Plus function--------------------------
    
    P_updo_pp_new = (sign*U)*(gamma_3_updo_pp_arr+gamma_3_xupdo_pp_arr)- K_updo_pp_new[N_bose_big_k - N_bose_big_p:N_bose_big_k + N_bose_big_p+1]
    print "Plus updo pp"

    #------------------------------------------------------------------------------------------------------------------------------
    #
    #                                               XPH
    #                                                                                    
    #
    #------------------------------------------------------------------------------------------------------------------------------
    print "XPH channel diagrams"
    
    #----------------------------------Plus functions and three-leg vertex in the xph notation---------------------------------------
    
    def int_gamma_3_xupdo_xph(i,j,k):
        return G(k-myfloor_div2(i))*G(k+myceil_div2(i))*f_updo_fun_ph(i,k,j)*1.0/beta #(-) from loop * (-) from swapping -> UPUP 
    
    def int_gamma_3_updo_xph(i,j,k):
        return G(k-myfloor_div2(i))*G(k+myceil_div2(i))*f_updo_fun_xph(i,k,j)*1.0/beta # no loops no swapping 
    
    
    # Arrays creation
    
    gamma_3_xupdo_xph_arr = np.array([[ sum_fit(partial(int_gamma_3_xupdo_xph,i,j)) for i in range(-N_bose_big_p,N_bose_big_p+1)] for j in range(-N_fermi_big, N_fermi_big)])
    gamma_3_updo_xph_arr = np.array([[ sum_fit(partial(int_gamma_3_updo_xph,i,j)) for i in range(-N_bose_big_p,N_bose_big_p+1)] for j in range(-N_fermi_big, N_fermi_big)])
    print "gammas3 xph"
    
    #-------------------------Definition of Plus function XPH--------------------------
    
    P_upup_xph_new = (sign*U)*gamma_3_xupdo_xph_arr - K_upup_xph_new[N_bose_big_k - N_bose_big_p:N_bose_big_k + N_bose_big_p+1]
    print "Plus upup xph"
    P_updo_xph_new = (sign*U)*gamma_3_updo_xph_arr - K_updo_xph_new[N_bose_big_k - N_bose_big_p:N_bose_big_k + N_bose_big_p+1] 
    print "Plus updo xph"
    #------------------------------------- END FUNCTION UPDATE ---------------------------------------------------------------------------
    
    # Compute variance
    epsip = np.max(map(mynorm, [P_upup_ph_new - P_upup_ph_loop, P_updo_ph_new - P_updo_ph_loop, P_updo_pp_new - P_updo_pp_loop, P_upup_xph_new - P_upup_xph_loop, P_updo_xph_new - P_updo_xph_loop]))
    print "Epsi plus:",epsip
    epsi = max(epsik,epsip)
    
    print "Convergency parameter epsi:",epsi
    
    # Update
    P_upup_ph_loop = P_upup_ph_new
    P_updo_ph_loop = P_updo_ph_new
    P_updo_pp_loop = P_updo_pp_new
    P_upup_xph_loop = P_upup_xph_new
    P_updo_xph_loop = P_updo_xph_new
    
    
    # Check after update
    if (epsi <= accuracy):
        break                           # Loop exit

# Now you are out of the loop!

# HDF5 NEW GROUP CREATION

if (('P_func' in f) or ('K_func' in f)):
    del f['P_func']
    del f['K_func']

p_grp = f.require_group("P_func")
p_subgrp_ph = p_grp.require_group("PH")
p_subgrp_pp = p_grp.require_group("PP")
p_subgrp_xph = p_grp.require_group("XPH")

# HDF5 NEW SUBGROUP CREATION

k_grp = f.require_group("K_func")
k_subgrp_ph = k_grp.require_group("PH")
k_subgrp_pp = k_grp.require_group("PP")
k_subgrp_xph = k_grp.require_group("XPH")

# HDF5 NEW DATASET CREATION

p_subgrp_ph.create_dataset("RE_P_UPUP", data= P_upup_ph_loop.real, dtype='float64', compression="gzip", compression_opts=4)
p_subgrp_ph.create_dataset("IM_P_UPUP", data= P_upup_ph_loop.imag, dtype='float64', compression="gzip", compression_opts=4)
p_subgrp_ph.create_dataset("RE_P_UPDO", data= P_updo_ph_loop.real, dtype='float64', compression="gzip", compression_opts=4)
p_subgrp_ph.create_dataset("IM_P_UPDO", data= P_updo_ph_loop.imag, dtype='float64', compression="gzip", compression_opts=4)
p_subgrp_pp.create_dataset("RE_P_UPDO", data= P_updo_pp_loop.real, dtype='float64', compression="gzip", compression_opts=4)
p_subgrp_pp.create_dataset("IM_P_UPDO", data= P_updo_pp_loop.imag, dtype='float64', compression="gzip", compression_opts=4)
p_subgrp_xph.create_dataset("RE_P_UPUP", data= P_upup_xph_loop.real, dtype='float64', compression="gzip", compression_opts=4)
p_subgrp_xph.create_dataset("IM_P_UPUP", data= P_upup_xph_loop.imag, dtype='float64', compression="gzip", compression_opts=4)
p_subgrp_xph.create_dataset("RE_P_UPDO", data= P_updo_xph_loop.real, dtype='float64', compression="gzip", compression_opts=4)
p_subgrp_xph.create_dataset("IM_P_UPDO", data= P_updo_xph_loop.imag, dtype='float64', compression="gzip", compression_opts=4)

k_subgrp_ph.create_dataset("RE_K_UPUP", data= K_upup_ph_loop.real, dtype='float64', compression="gzip", compression_opts=4)
k_subgrp_ph.create_dataset("IM_K_UPUP", data= K_upup_ph_loop.imag, dtype='float64', compression="gzip", compression_opts=4)
k_subgrp_ph.create_dataset("RE_K_UPDO", data= K_updo_ph_loop.real, dtype='float64', compression="gzip", compression_opts=4)
k_subgrp_ph.create_dataset("IM_K_UPDO", data= K_updo_ph_loop.imag, dtype='float64', compression="gzip", compression_opts=4)

k_subgrp_pp.create_dataset("RE_K_UPDO", data= K_updo_pp_loop.real, dtype='float64', compression="gzip", compression_opts=4)
k_subgrp_pp.create_dataset("IM_K_UPDO", data= K_updo_pp_loop.imag, dtype='float64', compression="gzip", compression_opts=4)
k_subgrp_xph.create_dataset("RE_K_UPUP", data= K_upup_xph_loop.real, dtype='float64', compression="gzip", compression_opts=4)
k_subgrp_xph.create_dataset("IM_K_UPUP", data= K_upup_xph_loop.imag, dtype='float64', compression="gzip", compression_opts=4)
k_subgrp_xph.create_dataset("RE_K_UPDO", data= K_updo_xph_loop.real, dtype='float64', compression="gzip", compression_opts=4)
k_subgrp_xph.create_dataset("IM_K_UPDO", data= K_updo_xph_loop.imag, dtype='float64', compression="gzip", compression_opts=4)

f.close()


#-------------------------------------CHECK PLOT FULL VERTEX WITH ASYMPTOTICS----------------
#------------------------------- PLOTTING FULL VERTEX -----------------------------------

print ("Plotting full Vetex...")

pl.rc('xtick', labelsize=9) 
pl.rc('ytick', labelsize=9) 
cmap = pl.get_cmap('jet') # Bianconiglio

pl.figsize=(13, 7)

bos = 0

def plotVert_ED( use_pl, zarr, string):
    use_pl.set_aspect(1.0)
    pl.pcolormesh(np.array([(2*i+1)*pi/beta for i in range(-N_fermi_sum,N_fermi_sum)]), np.array([(2*i+1)*pi/beta for i in range(-N_fermi_sum,N_fermi_sum)]),np.ma.masked_where( np.isnan(zarr), zarr ))
    use_pl.set_title( string , fontsize=10)    
    pl.colorbar(shrink=0.6)
    return

def plotre_f_upup_ph( use_pl ):
    title=r"$ReF^{PH}_{\uparrow \uparrow}$"
    zarr = np.array([[ f_upup_fun_ph(bos,n,m).real  for n in range(-N_fermi_sum,N_fermi_sum)] for m in range(-N_fermi_sum,N_fermi_sum)])
    plotVert_ED( use_pl, zarr, title ) 
    return

def plotre_f_updo_ph( use_pl ):
    title=r"$ReF^{PH}_{\uparrow \downarrow}$"
    zarr = np.array([[ f_updo_fun_ph(bos,n,m).real  for n in range(-N_fermi_sum,N_fermi_sum)] for m in range(-N_fermi_sum,N_fermi_sum)])
    plotVert_ED( use_pl, zarr, title ) 
    return

def plotre_f_xupdo_pp( use_pl ):
    title=r"$ReF^{PP}_{\uparrow \uparrow}$"
    zarr = np.array([[ f_xupdo_fun_pp(bos,n,m).real for n in range(-N_fermi_sum,N_fermi_sum)] for m in range(-N_fermi_sum,N_fermi_sum)])
    plotVert_ED( use_pl, zarr, title ) 
    return

def plotre_f_updo_pp( use_pl ):
    title=r"$ReF^{PP}_{\uparrow \downarrow}$"
    zarr = np.array([[ f_updo_fun_pp(bos,n,m).real for n in range(-N_fermi_sum,N_fermi_sum)] for m in range(-N_fermi_sum,N_fermi_sum)])
    plotVert_ED( use_pl, zarr, title ) 
    return

def plotre_f_upup_xph( use_pl ):
    title=r"$ReF^{XPH}_{2\uparrow \uparrow}$"
    zarr = np.array([[ (f_upup_fun_xph(bos,n,m)).real  for n in range(-N_fermi_sum,N_fermi_sum)] for m in range(-N_fermi_sum,N_fermi_sum)])
    plotVert_ED( use_pl, zarr, title ) 
    return

def plotre_f_updo_xph( use_pl ):
    title=r"$ReF^{XPH}_{2\uparrow \downarrow}$"
    zarr = np.array([[ (f_updo_fun_xph(bos,n,m)).real  for n in range(-N_fermi_sum,N_fermi_sum)] for m in range(-N_fermi_sum,N_fermi_sum)])
    plotVert_ED( use_pl, zarr, title ) 
    return

pl.suptitle(r"$U=$" + str('{0:.3f}'.format(float(U))) + r"     $\beta=$" + str('{0:.3f}'.format(float(beta))) + r"     $\Omega_{\rm PP}=$" + str(bos) + r"$*2\pi/\beta$")

plotre_f_xupdo_pp(pl.subplot(2,3,1))
pl.xlabel(r"$\omega_n$", fontsize=20)
pl.ylabel(r"$\omega_m$", fontsize=20)
pl.grid()
plotre_f_upup_ph(pl.subplot(2,3,2))
pl.xlabel(r"$\omega_n$", fontsize=20)
pl.ylabel(r"$\omega_m$", fontsize=20)
pl.grid()
plotre_f_upup_xph(pl.subplot(2,3,3))
pl.xlabel(r"$\omega_n$", fontsize=20)
pl.ylabel(r"$\omega_m$", fontsize=20)
pl.grid()
plotre_f_updo_pp(pl.subplot(2,3,4))
pl.xlabel(r"$\omega_n$", fontsize=20)
pl.ylabel(r"$\omega_m$", fontsize=20)
pl.grid()
plotre_f_updo_ph(pl.subplot(2,3,5))
pl.xlabel(r"$\omega_n$", fontsize=20)
pl.ylabel(r"$\omega_m$", fontsize=20)
pl.grid()
plotre_f_updo_xph(pl.subplot(2,3,6))
pl.xlabel(r"$\omega_n$", fontsize=20)
pl.ylabel(r"$\omega_m$", fontsize=20)
pl.grid()

pl.tight_layout()
#pl.show()

#--- Save to file
pl.savefig("Vertex_extended.png")

pl.clf()




