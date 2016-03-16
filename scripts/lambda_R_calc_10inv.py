#!/usr/bin/python
#
# This script provides the calculation of the following objects from the ED vertex:
#   - The "Rest" function subtracting from the phi the corresponding plus and karrasch functions
#   - The Fully Irreducible vertex Lam 
#
# The Gamma and the Phi exploit the asymptotic structures already obtainend from the self-consistency!
# In particular the Gamma is calculated using the extended Full vertex (with the asymptotics taken from the plus and the karrash functions)
# and taken the inversion
# Everything will be stored in the HDF5 FILE which already contains all the info about the ED calculation of the vertex
#
# plotting using plot.py
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
from _functools import partial

#-----------------------------------Read parameters----------------------------------------------
# ------U---------
U = 1.0
#-------BETA-------
beta = 10.0

#----------------------------------------Read HDF5 files-----------------------------------------

if ('../dat'):
    f = h5py.File('../dat/10inv/dat_U'+ str(U)+'_beta'+ str(beta)+'_EDpomerol.h5', 'r+')   # Read (and write) the hdf5 file in the directory "dat" if existing
else:
    sys.exit("No data file")

#---------  PARAMETERS  ---------------

pi = math.pi

#---------------------- Read objects involved in the calculation from HDF5 file -------------------------------

#--------------------------------- READING ----------------------------------------------------------

#GF

re_g_iw = f["Giw/RE"][:]
im_g_iw = f["Giw/IM"][:]           

print im_g_iw.shape[0] 
fgrid_gf = im_g_iw.shape[0]
N_fermi_gf = fgrid_gf/2

#2PGF

re_2pgf_upup_ph = f["2PGF/PH/RE_2PGF_UPUP"][:]
re_2pgf_updo_ph = f["2PGF/PH/RE_2PGF_UPDO"][:]
re_2pgf_upup_xph = f["2PGF/XPH/RE_2PGF_UPUP"][:]
re_2pgf_updo_xph = f["2PGF/XPH/RE_2PGF_UPDO"][:]
re_2pgf_upup_pp = f["2PGF/PP/RE_2PGF_UPUP"][:]
re_2pgf_updo_pp = f["2PGF/PP/RE_2PGF_UPDO"][:]


im_2pgf_upup_ph = f["2PGF/PH/IM_2PGF_UPUP"][:]
im_2pgf_updo_ph = f["2PGF/PH/IM_2PGF_UPDO"][:]
im_2pgf_upup_xph = f["2PGF/XPH/IM_2PGF_UPUP"][:]
im_2pgf_updo_xph = f["2PGF/XPH/IM_2PGF_UPDO"][:]
im_2pgf_upup_pp = f["2PGF/PP/IM_2PGF_UPUP"][:]
im_2pgf_updo_pp = f["2PGF/PP/IM_2PGF_UPDO"][:]

#GENERALIZED SUSCEPTIBILITY

re_chi_upup_ph = f["GENCHI/PH/RE_GENCHI_UPUP"][:]
re_chi_updo_ph = f["GENCHI/PH/RE_GENCHI_UPDO"][:]
re_chi_upup_xph = f["GENCHI/XPH/RE_GENCHI_UPUP"][:]
re_chi_updo_xph = f["GENCHI/XPH/RE_GENCHI_UPDO"][:]
re_chi_upup_pp = f["GENCHI/PP/RE_GENCHI_UPUP"][:]
re_chi_updo_pp = f["GENCHI/PP/RE_GENCHI_UPDO"][:]


im_chi_upup_ph = f["GENCHI/PH/IM_GENCHI_UPUP"][:]
im_chi_updo_ph = f["GENCHI/PH/IM_GENCHI_UPDO"][:]
im_chi_upup_xph = f["GENCHI/XPH/IM_GENCHI_UPUP"][:]
im_chi_updo_xph = f["GENCHI/XPH/IM_GENCHI_UPDO"][:]
im_chi_upup_pp = f["GENCHI/PP/IM_GENCHI_UPUP"][:]
im_chi_updo_pp = f["GENCHI/PP/IM_GENCHI_UPDO"][:]

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

# We assume all the channels to have the same B/F grids

fgrid = f["VERT/PH/fgrid"][:].shape[0] 
bgrid = f["VERT/PH/bgrid"][:].shape[0]

N_bose = (bgrid-1)/2 # to create a bosonic frequency grid from -N_bose to N_bose
N_fermi= (fgrid)/2   # to create a fermionic frequency grid from -N_fermi to N_fermi

#Plus

rep_upup_ph =  np.array(f["/P_func/PH/RE_P_UPUP"])
imp_upup_ph = np.array(f["/P_func/PH/IM_P_UPUP"])
rep_updo_ph =  np.array(f["/P_func/PH/RE_P_UPDO"])
imp_updo_ph = np.array(f["/P_func/PH/IM_P_UPDO"])
rep_updo_pp =  np.array(f["/P_func/PP/RE_P_UPDO"])
imp_updo_pp = np.array(f["/P_func/PP/IM_P_UPDO"])
rep_updo_xph =  np.array(f["/P_func/XPH/RE_P_UPDO"])
imp_updo_xph = np.array(f["/P_func/XPH/IM_P_UPDO"])
rep_upup_xph =  np.array(f["/P_func/XPH/RE_P_UPUP"])
imp_upup_xph = np.array(f["/P_func/XPH/IM_P_UPUP"])

#Karrasch

rek_upup_ph =  np.array(f["/K_func/PH/RE_K_UPUP"])
imk_upup_ph = np.array(f["/K_func/PH/IM_K_UPUP"])
rek_updo_ph =  np.array(f["/K_func/PH/RE_K_UPDO"])
imk_updo_ph = np.array(f["/K_func/PH/IM_K_UPDO"])
rek_updo_pp =  np.array(f["/K_func/PP/RE_K_UPDO"])
imk_updo_pp = np.array(f["/K_func/PP/IM_K_UPDO"])
rek_updo_xph =  np.array(f["/K_func/XPH/RE_K_UPDO"])
imk_updo_xph = np.array(f["/K_func/XPH/IM_K_UPDO"])
rek_upup_xph =  np.array(f["/K_func/XPH/RE_K_UPUP"])
imk_upup_xph = np.array(f["/K_func/XPH/IM_K_UPUP"])

N_bose_k = (rek_upup_ph.shape[0] -1)/2
N_bose_p = (rep_upup_ph.shape[1]-1)/2
N_fermi_p = (rep_upup_ph.shape[0])/2

#GRIDS

fgrid_arr = np.array([(2*m +1)*pi/beta for m in range (-N_fermi, N_fermi)])
bgrid_arr = np.array([(2*n)*pi/beta for n in range (-N_bose, N_bose+1)])

#Inside box

def isInside(i,j,k):
    return abs(i) <= N_bose and j >= -N_fermi and j < N_fermi and k >= -N_fermi and k < N_fermi

#-------------------------------- FUNCTION DEFINITION--------------------------------------

#-----------GF------------------------------------------------

def G(wf):
    if (wf >= -N_fermi_gf and wf < N_fermi_gf):
        return re_g_iw[wf+N_fermi_gf]+1j*im_g_iw[wf+N_fermi_gf]
    else:
        wMat = 1j*(2*wf +1)*pi/beta     #GF asymtotic behavior
        return 1.0/wMat


#----------------- BARE BUBBLES IN ALL CHANNELS -------------------
#PP
def chi_0_pp(wb,wf,wf1):
    if (wf == wf1):
        return beta*G(wf+myceil_div2(wb))*G(myfloor_div2(wb)-wf-1)
    else:
        return 0.0
#PH
def chi_0_ph(wb,wf,wf1):
    if (wb == 0):
        return beta*G(wf)*G(wf1)
    else:
        return 0.0

def chi_x0_ph(wb,wf,wf1):
    if (wf == wf1):
        return -beta*G(wf-myfloor_div2(wb))*G(wf1+myceil_div2(wb))
    else:
        return 0.0

#XPH
def chi_0_xph(wb,wf,wf1):
    if (wf == wf1):
        return beta*G(wf-myfloor_div2(wb))*G(wf+myceil_div2(wb))
    else:
        return 0.0

#---------------------------KARRASCH--------------------------

def K_upup_ph(wb):
    if (abs(wb) <= N_bose_k):
        return rek_upup_ph[wb+N_bose_k]+1j*imk_upup_ph[wb+N_bose_k]
    else:
        return 0.0

def K_updo_ph(wb):
    if (abs(wb) <= N_bose_k):
        return rek_updo_ph[wb+N_bose_k]+1j*imk_updo_ph[wb+N_bose_k]
    else:
        return 0.0

def K_updo_pp(wb):
    if (abs(wb) <= N_bose_k):
        return rek_updo_pp[wb+N_bose_k]+1j*imk_updo_pp[wb+N_bose_k]
    else:
        return 0.0

def K_upup_xph(wb):
    if (abs(wb) <= N_bose_k):
        return rek_upup_xph[wb+N_bose_k]+1j*imk_upup_xph[wb+N_bose_k]
    else:
        return 0.0

def K_updo_xph(wb):
    if (abs(wb) <= N_bose_k):
        return rek_updo_xph[wb+N_bose_k]+1j*imk_updo_xph[wb+N_bose_k]
    else:
        return 0.0

#---------------------------P func--------------------------------

def P_upup_ph(wb,wf):
    if (abs(wb) <= N_bose_p and wf >= -N_fermi_p and wf < N_fermi_p):
        return rep_upup_ph[wf+N_fermi_p,wb+N_bose_p]+1j*imp_upup_ph[wf+N_fermi_p,wb+N_bose_p]
    else:
#        print "P_upup_ph out"
        return 0.0
def P_updo_ph(wb,wf):
    if (abs(wb) <= N_bose_p and wf >= -N_fermi_p and wf < N_fermi_p):
        return rep_updo_ph[wf+N_fermi_p,wb+N_bose_p]+1j*imp_updo_ph[wf+N_fermi_p,wb+N_bose_p]
    else:
#        print "P_updo_ph out"
        return 0.0
def P_updo_pp(wb,wf):
    if (abs(wb) <= N_bose_p and wf >= -N_fermi_p and wf < N_fermi_p):
        return rep_updo_pp[wf+N_fermi_p,wb+N_bose_p]+1j*imp_updo_pp[wf+N_fermi_p,wb+N_bose_p]
    else:
#        print "P_updo_pp out"
        return 0.0
def P_upup_xph(wb,wf):
    if (abs(wb) <= N_bose_p and wf >= -N_fermi_p and wf < N_fermi_p):
        return rep_upup_xph[wf+N_fermi_p,wb+N_bose_p]+1j*imp_upup_xph[wf+N_fermi_p,wb+N_bose_p]
    else:
#        print "P_upup_xph out"
        return 0.0
def P_updo_xph(wb,wf):
    if (abs(wb) <= N_bose_p and wf >= -N_fermi_p and wf < N_fermi_p):
        return rep_updo_xph[wf+N_fermi_p,wb+N_bose_p]+1j*imp_updo_xph[wf+N_fermi_p,wb+N_bose_p]
    else:
#        print "K_updo_xph out"
        return 0.0

# Update the asymptotic structures for the VERTEX IN ALL CHANNELS

def f_upup_fun_ph(i,j,k):
    if isInside(i,j,k):
        return re_f_upup_ph[i + N_bose, j+N_fermi, k + N_fermi]+1j*im_f_upup_ph[i + N_bose, j+N_fermi, k + N_fermi]
    else:
        return K_upup_ph(i) + P_upup_ph(i,j)+ P_upup_ph(i,k) + K_upup_xph(PHtoXPH((i,j,k))[0]) + P_upup_xph(PHtoXPH((i,j,k))[0],PHtoXPH((i,j,k))[1])+P_upup_xph(PHtoXPH((i,j,k))[0],PHtoXPH((i,j,k))[2])
           
def f_updo_fun_ph(i,j,k):
    if isInside(i,j,k):
        return re_f_updo_ph[i + N_bose, j+N_fermi, k + N_fermi]+1j*im_f_updo_ph[i + N_bose, j+N_fermi, k + N_fermi]
    else:
        return - U + K_updo_ph(i) + P_updo_ph(i,j)+P_updo_ph(i,k) + K_updo_xph(PHtoXPH((i,j,k))[0]) + P_updo_xph(PHtoXPH((i,j,k))[0],PHtoXPH((i,j,k))[1]) + P_updo_xph(PHtoXPH((i,j,k))[0],PHtoXPH((i,j,k))[2])+ K_updo_pp(PHtoPP((i,j,k))[0]) + P_updo_pp(PHtoPP((i,j,k))[0],PHtoPP((i,j,k))[1]) + P_updo_pp(PHtoPP((i,j,k))[0],PHtoPP((i,j,k))[2])

def f_upup_fun_pp(i,j,k):
    if isInside(i,j,k):
        return re_f_upup_pp[i + N_bose, j+N_fermi, k + N_fermi]+1j*im_f_updo_pp[i + N_bose, j+N_fermi, k + N_fermi]
    else:
        return  K_upup_ph(PPtoPH((i,j,k))[0]) + P_upup_ph(PPtoPH((i,j,k))[0],PPtoPH((i,j,k))[1]) + P_upup_ph(PPtoPH((i,j,k))[0],PPtoPH((i,j,k))[2])+ K_upup_xph(PPtoXPH((i,j,k))[0]) + P_upup_xph(PPtoXPH((i,j,k))[0],PPtoXPH((i,j,k))[1]) + P_upup_xph(PPtoXPH((i,j,k))[0],PPtoXPH((i,j,k))[2])

def f_updo_fun_pp(i,j,k):
    if isInside(i,j,k):
        return re_f_updo_pp[i + N_bose, j+N_fermi, k + N_fermi]+1j*im_f_updo_pp[i + N_bose, j+N_fermi, k + N_fermi]
    else:
        return - U + K_updo_pp(i) + P_updo_pp(i,j) + P_updo_pp(i,k) + K_updo_ph(PPtoPH((i,j,k))[0]) + P_updo_ph(PPtoPH((i,j,k))[0],PPtoPH((i,j,k))[1]) + P_updo_ph(PPtoPH((i,j,k))[0],PPtoPH((i,j,k))[2])+ K_updo_xph(PPtoXPH((i,j,k))[0]) + P_updo_xph(PPtoXPH((i,j,k))[0],PPtoXPH((i,j,k))[1]) + P_updo_xph(PPtoXPH((i,j,k))[0],PPtoXPH((i,j,k))[2])

def f_xupdo_fun_pp(i,j,k):
    if isInside(i,j,k):
        return re_f_upup_pp[i + N_bose, j+N_fermi, k + N_fermi]+1j*im_f_upup_pp[i + N_bose, j+N_fermi, k + N_fermi]-re_f_updo_pp[i + N_bose, j+N_fermi, k + N_fermi]-1j*im_f_updo_pp[i + N_bose, j+N_fermi, k + N_fermi]
    else:
        return U - K_updo_pp(i) - P_updo_pp(i,j) - P_updo_pp(i,k) - K_updo_xph(PPtoXPH((i,j,-k-mymod_abs(i)-1))[0]) - P_updo_xph(PPtoXPH((i,j,-k-mymod_abs(i)-1))[0],PPtoXPH((i,j,-k-mymod_abs(i)-1))[1]) - P_updo_xph(PPtoXPH((i,j,-k-mymod_abs(i)-1))[0],PPtoXPH((i,j,-k-mymod_abs(i)-1))[2])- K_updo_ph(PPtoPH((i,j,-k-mymod_abs(i)-1))[0]) - P_updo_ph(PPtoPH((i,j,-k-mymod_abs(i)-1))[0],PPtoPH((i,j,-k-mymod_abs(i)-1))[1]) - P_updo_ph(PPtoPH((i,j,-k-mymod_abs(i)-1))[0],PPtoPH((i,j,-k-mymod_abs(i)-1))[2])    

def f_upup_fun_xph(i,j,k):
    if isInside(i,j,k):
        return re_f_upup_xph[i + N_bose, j+N_fermi, k + N_fermi]+1j*im_f_upup_xph[i + N_bose, j+N_fermi, k + N_fermi]
    else:
        return K_upup_xph(i) + P_upup_xph(i,j)+P_upup_xph(i,k) + K_upup_ph(XPHtoPH((i,j,k))[0]) + P_upup_ph(XPHtoPH((i,j,k))[0],XPHtoPH((i,j,k))[1])+ P_upup_ph(XPHtoPH((i,j,k))[0],XPHtoPH((i,j,k))[2]) 

def f_updo_fun_xph(i,j,k):
    if isInside(i,j,k):
        return re_f_updo_xph[i + N_bose, j+N_fermi, k + N_fermi]+1j*im_f_updo_xph[i + N_bose, j+N_fermi, k + N_fermi]
    else:
        return - U + K_updo_xph(i) + P_updo_xph(i,j)+P_updo_xph(i,k) + K_updo_ph(XPHtoPH((i,j,k))[0]) + P_updo_ph(XPHtoPH((i,j,k))[0],XPHtoPH((i,j,k))[1]) + P_updo_ph(XPHtoPH((i,j,k))[0],XPHtoPH((i,j,k))[2])+ K_updo_pp(XPHtoPP((i,j,k))[0]) + P_updo_pp(XPHtoPP((i,j,k))[0],XPHtoPP((i,j,k))[1]) + P_updo_pp(XPHtoPP((i,j,k))[0],XPHtoPP((i,j,k))[2])

#--------------------------------GENCHI FUNC--------------------------------------
#PP -> chi (singlet, triplet)

def chi_s_pp(wb,wf,wf1):
    if isInside(wb,wf,wf1):
        return -re_chi_upup_pp[wb+N_bose,wf+N_fermi,wf1+N_fermi]-1j*im_chi_upup_pp[wb+N_bose,wf+N_fermi,wf1+N_fermi]+ 2*(re_chi_updo_pp[wb+N_bose,wf+N_fermi,wf1+N_fermi] +1j*im_chi_updo_pp[wb+N_bose,wf+N_fermi,wf1+N_fermi])
    else:
        return G(wf+myceil_div2(wb))*G(-wf-1+myfloor_div2(wb))*(-f_upup_fun_pp(wb,wf,wf1)+2*f_updo_fun_pp(wb,wf,wf1))*G(wf1+myceil_div2(wb))*G(-wf1-1+myfloor_div2(wb))+chi_0_pp(wb,wf,wf1)
def chi_t_pp(wb,wf,wf1):
    if isInside(wb,wf,wf1):
        return re_chi_upup_pp[wb+N_bose,wf+N_fermi,wf1+N_fermi] + 1j*im_chi_upup_pp[wb+N_bose,wf+N_fermi,wf1+N_fermi] 
    else:
        return G(wf+myceil_div2(wb))*G(-wf-1+myfloor_div2(wb))*(f_upup_fun_pp(wb,wf,wf1))*G(wf1+myceil_div2(wb))*G(-wf1-1+myfloor_div2(wb))+chi_0_pp(wb,wf,wf1)

#PH -> 2pgf(magnetic, density)

def tpgf_d_ph(wb,wf,wf1):
    if isInside(wb,wf,wf1):
        return re_2pgf_upup_ph[wb+N_bose,wf+N_fermi,wf1+N_fermi]+1j*im_2pgf_upup_ph[wb+N_bose,wf+N_fermi,wf1+N_fermi]+re_2pgf_updo_ph[wb+N_bose,wf+N_fermi,wf1+N_fermi]+1j*im_2pgf_updo_ph[wb+N_bose,wf+N_fermi,wf1+N_fermi]
    else:
        return G(wf+myceil_div2(wb))*G(wf-myfloor_div2(wb))*(f_upup_fun_ph(wb,wf,wf1)+f_updo_fun_ph(wb,wf,wf1))*G(wf1+myceil_div2(wb))*G(wf1-myfloor_div2(wb))+2*chi_0_ph(wb,wf,wf1)+chi_x0_ph(wb,wf,wf1)
def tpgf_m_ph(wb,wf,wf1):
    if isInside(wb,wf,wf1):
        return re_2pgf_upup_ph[wb+N_bose,wf+N_fermi,wf1+N_fermi]+1j*im_2pgf_upup_ph[wb+N_bose,wf+N_fermi,wf1+N_fermi] -re_2pgf_updo_ph[wb+N_bose,wf+N_fermi,wf1+N_fermi]-1j*im_2pgf_updo_ph[wb+N_bose,wf+N_fermi,wf1+N_fermi]
    else:
        return G(wf+myceil_div2(wb))*G(wf-myfloor_div2(wb))*(f_upup_fun_ph(wb,wf,wf1)-f_updo_fun_ph(wb,wf,wf1))*G(wf1+myceil_div2(wb))*G(wf1-myfloor_div2(wb))+chi_x0_ph(wb,wf,wf1)

#-----------------------------GAMMA BETHE SALPETER INVERSION

N_fermi_inv = 10*N_fermi

def chis_chi0_arr(wb):
    return np.array([[chi_s_pp(wb,wf,wf1) + chi_0_pp(wb,wf,wf1) for wf in range(-N_fermi_inv,N_fermi_inv)] for wf1 in range(-N_fermi_inv,N_fermi_inv)])

def chit_chi0_arr(wb):
    return np.array([[chi_t_pp(wb,wf,wf1) + chi_0_pp(wb,wf,wf1) for wf in range(-N_fermi_inv,N_fermi_inv)] for wf1 in range(-N_fermi_inv,N_fermi_inv)])

def chi_0_pp_arr(wb):
    return np.array([[chi_0_pp(wb,wf,wf1) for wf in range(-N_fermi_inv,N_fermi_inv)] for wf1 in range(-N_fermi_inv,N_fermi_inv)])

gamma_t_arr = np.array([beta*beta*(2*inv(chi_0_pp_arr(wb))- 4*inv(chit_chi0_arr(wb))) for wb in range (-N_bose,N_bose+1)])
gamma_s_arr = np.array([beta*beta*(2*inv(chi_0_pp_arr(wb)) - 4*inv(chis_chi0_arr(wb))) for wb in range (-N_bose,N_bose+1)])

gamma_upup_pp_arr = gamma_t_arr[:,N_fermi_inv-N_fermi:N_fermi_inv+N_fermi,N_fermi_inv-N_fermi:N_fermi_inv+N_fermi] 
gamma_updo_pp_arr = 0.5*(gamma_t_arr[:,N_fermi_inv-N_fermi:N_fermi_inv+N_fermi,N_fermi_inv-N_fermi:N_fermi_inv+N_fermi] + gamma_s_arr[:,N_fermi_inv-N_fermi:N_fermi_inv+N_fermi,N_fermi_inv-N_fermi:N_fermi_inv+N_fermi]) 

#PH

def m_d_arr_ph(wb):
    return np.array([[tpgf_d_ph(wb,wf,wf1)-2*chi_0_ph(wb,wf,wf1) for wf in range(-N_fermi_inv,N_fermi_inv)] for wf1 in range(-N_fermi_inv,N_fermi_inv)])

def m_m_arr_ph(wb):
    return np.array([[tpgf_m_ph(wb,wf,wf1) for wf in range(-N_fermi_inv,N_fermi_inv)] for wf1 in range(-N_fermi_inv,N_fermi_inv)])

def chi_x0_ph_arr(wb):
    return np.array([[chi_x0_ph(wb,wf,wf1) for wf in range(-N_fermi_inv,N_fermi_inv)] for wf1 in range(-N_fermi_inv,N_fermi_inv)])


gamma_d_arr = np.array([beta*beta*(inv(chi_x0_ph_arr(wb))-inv(m_d_arr_ph(wb))) for wb in range (-N_bose,N_bose+1)])
gamma_m_arr = np.array([beta*beta*(inv(chi_x0_ph_arr(wb))-inv(m_m_arr_ph(wb))) for wb in range (-N_bose,N_bose+1)])

gamma_upup_ph_arr = 0.5*(gamma_d_arr[:,N_fermi_inv-N_fermi:N_fermi_inv+N_fermi,N_fermi_inv-N_fermi:N_fermi_inv+N_fermi] + gamma_m_arr[:,N_fermi_inv-N_fermi:N_fermi_inv+N_fermi,N_fermi_inv-N_fermi:N_fermi_inv+N_fermi])
gamma_updo_ph_arr = 0.5*(gamma_d_arr[:,N_fermi_inv-N_fermi:N_fermi_inv+N_fermi,N_fermi_inv-N_fermi:N_fermi_inv+N_fermi] - gamma_m_arr[:,N_fermi_inv-N_fermi:N_fermi_inv+N_fermi,N_fermi_inv-N_fermi:N_fermi_inv+N_fermi])

#XPH

def m_updo_fun_xph(wb,wf,wf1):
    if isInside(wb,wf,wf1):
        return re_2pgf_updo_xph[wb+N_bose,wf+N_fermi,wf1+N_fermi]+1j*im_2pgf_updo_xph[wb+N_bose,wf+N_fermi,wf1+N_fermi]
    else:
        return G(wf+myceil_div2(wb))*G(wf-myfloor_div2(wb))*(f_updo_fun_xph(wb,wf,wf1))*G(wf1+myceil_div2(wb))*G(wf1-myfloor_div2(wb))+chi_0_xph(wb,wf,wf1)

def m_updo_arr_xph(wb):
    return np.array([[m_updo_fun_xph(wb,wf,wf1) for wf in range(-N_fermi_inv,N_fermi_inv)] for wf1 in range(-N_fermi_inv,N_fermi_inv)])


def chi_0_xph_arr(wb):
    return np.array([[chi_0_xph(wb,wf,wf1) for wf in range(-N_fermi_inv,N_fermi_inv)] for wf1 in range(-N_fermi_inv,N_fermi_inv)])

gamma_upup_xph_arr = -gamma_upup_ph_arr
gamma_updo_xph_arr_big = np.array([beta*beta*(+inv(chi_0_xph_arr(wb))-inv(m_updo_arr_xph(wb))) for wb in range (-N_bose,N_bose+1)])
gamma_updo_xph_arr = gamma_updo_xph_arr_big[:,N_fermi_inv-N_fermi:N_fermi_inv+N_fermi,N_fermi_inv-N_fermi:N_fermi_inv+N_fermi]


#---------------------------Store in hdf5 file---------------------------------------

if (('GAMMA' in f)):
    del f['GAMMA']

if ('PHI' in f):
    del f['PHI']

gamma_grp = f.require_group("GAMMA")
gamma_subgrp_ph = gamma_grp.require_group("PH")
gamma_subgrp_pp = gamma_grp.require_group("PP")
gamma_subgrp_xph = gamma_grp.require_group("XPH")
#PP

gamma_subgrp_pp.create_dataset('RE_GAMMA_UPUP', data= gamma_upup_pp_arr.real, dtype='float64', compression="gzip", compression_opts=4)
gamma_subgrp_pp.create_dataset('IM_GAMMA_UPUP', data= gamma_upup_pp_arr.imag, dtype='float64', compression="gzip", compression_opts=4)
gamma_subgrp_pp.create_dataset('RE_GAMMA_UPDO', data= gamma_updo_pp_arr.real, dtype='float64', compression="gzip", compression_opts=4)
gamma_subgrp_pp.create_dataset('IM_GAMMA_UPDO', data= gamma_updo_pp_arr.imag, dtype='float64', compression="gzip", compression_opts=4)
gamma_subgrp_pp.create_dataset('fgrid', data=fgrid_arr, dtype='float64', compression="gzip", compression_opts=4)
gamma_subgrp_pp.create_dataset('bgrid', data=bgrid_arr, dtype='float64', compression="gzip", compression_opts=4)

#PH
gamma_subgrp_ph.create_dataset('RE_GAMMA_UPUP', data= gamma_upup_ph_arr.real, dtype='float64', compression="gzip", compression_opts=4)
gamma_subgrp_ph.create_dataset('IM_GAMMA_UPUP', data= gamma_upup_ph_arr.imag, dtype='float64', compression="gzip", compression_opts=4)
gamma_subgrp_ph.create_dataset('RE_GAMMA_UPDO', data= gamma_updo_ph_arr.real, dtype='float64', compression="gzip", compression_opts=4)
gamma_subgrp_ph.create_dataset('IM_GAMMA_UPDO', data= gamma_updo_ph_arr.imag, dtype='float64', compression="gzip", compression_opts=4)
gamma_subgrp_ph.create_dataset('fgrid', data=fgrid_arr, dtype='float64', compression="gzip", compression_opts=4)
gamma_subgrp_ph.create_dataset('bgrid', data=bgrid_arr, dtype='float64', compression="gzip", compression_opts=4)

#XPH

gamma_subgrp_xph.create_dataset('RE_GAMMA_UPUP', data= gamma_upup_xph_arr.real, dtype='float64', compression="gzip", compression_opts=4)
gamma_subgrp_xph.create_dataset('IM_GAMMA_UPUP', data= gamma_upup_xph_arr.imag, dtype='float64', compression="gzip", compression_opts=4)
gamma_subgrp_xph.create_dataset('RE_GAMMA_UPDO', data= gamma_updo_xph_arr.real, dtype='float64', compression="gzip", compression_opts=4)
gamma_subgrp_xph.create_dataset('IM_GAMMA_UPDO', data= gamma_updo_xph_arr.imag, dtype='float64', compression="gzip", compression_opts=4)
gamma_subgrp_xph.create_dataset('fgrid', data=fgrid_arr, dtype='float64', compression="gzip", compression_opts=4)
gamma_subgrp_xph.create_dataset('bgrid', data=bgrid_arr, dtype='float64', compression="gzip", compression_opts=4)


#PHIS FUNCTIONS


def phi_upup_fun_ph(i,j,k):
    if isInside(i,j,k):
        return f_upup_fun_ph(i,j,k)-gamma_upup_ph_arr[i + N_bose, j+N_fermi, k + N_fermi] 
    else:
       return K_upup_ph(i) + P_upup_ph(i,j) + P_upup_ph(i,k) 
    
def phi_updo_fun_ph(i,j,k):
    if isInside(i,j,k):
        return f_updo_fun_ph(i,j,k)-gamma_updo_ph_arr[i + N_bose, j+N_fermi, k + N_fermi]
    else:
        return K_updo_ph(i) + P_updo_ph(i,j)+P_updo_ph(i,k)
    
def phi_upup_fun_pp(i,j,k):
    if isInside(i,j,k):
        return f_upup_fun_pp(i,j,k)-gamma_upup_pp_arr[i + N_bose, j+N_fermi, k + N_fermi]
    else:
        return 0.0 
    
def phi_updo_fun_pp(i,j,k):
    if isInside(i,j,k):
        return f_updo_fun_pp(i,j,k)-gamma_updo_pp_arr[i + N_bose, j+N_fermi, k + N_fermi] 
    else:
        return K_updo_pp(i) + P_updo_pp(i,j) + P_updo_pp(i,k)
    
def phi_upup_fun_xph(i,j,k):
    if isInside(i,j,k):
        return f_upup_fun_xph(i,j,k)-gamma_upup_xph_arr[i + N_bose, j+N_fermi, k + N_fermi]
    else:
        return K_upup_xph(i) + P_upup_xph(i,j)  + P_upup_xph(i,k)
    
def phi_updo_fun_xph(i,j,k):
    if isInside(i,j,k):
        return f_updo_fun_xph(i,j,k)-gamma_updo_xph_arr[i + N_bose, j+N_fermi, k + N_fermi]
    else:
        return K_updo_xph(i) + P_updo_xph(i,j)+ P_updo_xph(i,k)

#---------------------Store phi on hdf5 file-----------------------------------------
phi_upup_ph_arr = np.array([[[phi_upup_fun_ph(i,j,k) for j in range(-N_fermi,N_fermi)] for k in range(-N_fermi,N_fermi)] for i in range(-N_bose, N_bose+1)])
phi_updo_ph_arr = np.array([[[phi_updo_fun_ph(i,j,k) for j in range(-N_fermi,N_fermi)] for k in range(-N_fermi,N_fermi)] for i in range(-N_bose, N_bose+1)])
phi_upup_pp_arr = np.array([[[phi_upup_fun_pp(i,j,k) for j in range(-N_fermi,N_fermi)] for k in range(-N_fermi,N_fermi)] for i in range(-N_bose, N_bose+1)])
phi_updo_pp_arr = np.array([[[phi_updo_fun_pp(i,j,k) for j in range(-N_fermi,N_fermi)] for k in range(-N_fermi,N_fermi)] for i in range(-N_bose, N_bose+1)])
phi_upup_xph_arr = np.array([[[phi_upup_fun_xph(i,j,k) for j in range(-N_fermi,N_fermi)] for k in range(-N_fermi,N_fermi)] for i in range(-N_bose, N_bose+1)])
phi_updo_xph_arr = np.array([[[phi_updo_fun_xph(i,j,k) for j in range(-N_fermi,N_fermi)] for k in range(-N_fermi,N_fermi)] for i in range(-N_bose, N_bose+1)])

phi_grp = f.require_group("PHI")
phi_subgrp_ph = phi_grp.require_group("PH")
phi_subgrp_pp = phi_grp.require_group("PP")
phi_subgrp_xph = phi_grp.require_group("XPH")
#PP
phi_subgrp_pp.create_dataset('RE_PHI_UPUP', data= phi_upup_pp_arr.real, dtype='float64', compression="gzip", compression_opts=4)
phi_subgrp_pp.create_dataset('IM_PHI_UPUP', data= phi_upup_pp_arr.imag, dtype='float64', compression="gzip", compression_opts=4)
phi_subgrp_pp.create_dataset('RE_PHI_UPDO', data= phi_updo_pp_arr.real, dtype='float64', compression="gzip", compression_opts=4)
phi_subgrp_pp.create_dataset('IM_PHI_UPDO', data= phi_updo_pp_arr.imag, dtype='float64', compression="gzip", compression_opts=4)
phi_subgrp_pp.create_dataset('fgrid', data=fgrid_arr, dtype='float64', compression="gzip", compression_opts=4)
phi_subgrp_pp.create_dataset('bgrid', data=bgrid_arr, dtype='float64', compression="gzip", compression_opts=4)

#PH

phi_subgrp_ph.create_dataset('RE_PHI_UPUP', data= phi_upup_ph_arr.real, dtype='float64', compression="gzip", compression_opts=4)
phi_subgrp_ph.create_dataset('IM_PHI_UPUP', data= phi_upup_ph_arr.imag, dtype='float64', compression="gzip", compression_opts=4)
phi_subgrp_ph.create_dataset('RE_PHI_UPDO', data= phi_updo_ph_arr.real, dtype='float64', compression="gzip", compression_opts=4)
phi_subgrp_ph.create_dataset('IM_PHI_UPDO', data= phi_updo_ph_arr.imag, dtype='float64', compression="gzip", compression_opts=4)
phi_subgrp_ph.create_dataset('fgrid', data=fgrid_arr, dtype='float64', compression="gzip", compression_opts=4)
phi_subgrp_ph.create_dataset('bgrid', data=bgrid_arr, dtype='float64', compression="gzip", compression_opts=4)

#XPH

phi_subgrp_xph.create_dataset('RE_PHI_UPUP', data= phi_upup_xph_arr.real, dtype='float64', compression="gzip", compression_opts=4)
phi_subgrp_xph.create_dataset('IM_PHI_UPUP', data= phi_upup_xph_arr.imag, dtype='float64', compression="gzip", compression_opts=4)
phi_subgrp_xph.create_dataset('RE_PHI_UPDO', data= phi_updo_xph_arr.real, dtype='float64', compression="gzip", compression_opts=4)
phi_subgrp_xph.create_dataset('IM_PHI_UPDO', data= phi_updo_xph_arr.imag, dtype='float64', compression="gzip", compression_opts=4)
phi_subgrp_xph.create_dataset('fgrid', data=fgrid_arr, dtype='float64', compression="gzip", compression_opts=4)
phi_subgrp_xph.create_dataset('bgrid', data=bgrid_arr, dtype='float64', compression="gzip", compression_opts=4)

#-----------------------------------REST FUNCTION----------------------------------------

def R_upup_ph(i,j,k):
    return phi_upup_fun_ph(i,j,k)-P_upup_ph(i,j)-P_upup_ph(i,k)-K_upup_ph(i)

def R_updo_ph(i,j,k):
    return phi_updo_fun_ph(i,j,k)-P_updo_ph(i,j)-P_updo_ph(i,k)-K_updo_ph(i)

def R_upup_pp(i,j,k):
    return phi_upup_fun_pp(i,j,k)

def R_updo_pp(i,j,k):
    return phi_updo_fun_pp(i,j,k)-P_updo_pp(i,j)-P_updo_pp(i,k)-K_updo_pp(i)

def R_upup_xph(i,j,k):
    return phi_upup_fun_xph(i,j,k)-P_upup_xph(i,j)-P_upup_xph(i,k)-K_upup_xph(i)

def R_updo_xph(i,j,k):
    return phi_updo_fun_xph(i,j,k)-P_updo_xph(i,j)-P_updo_xph(i,k)-K_updo_xph(i)

#---------------------------------------LAMBDA--------------------------------------------

def Lam_upup_ph(wb,wf,wf1):
    return f_upup_fun_ph(wb,wf,wf1) - phi_upup_fun_ph(wb,wf,wf1) - phi_upup_fun_pp(PHtoPP((wb,wf,wf1))[0],PHtoPP((wb,wf,wf1))[1],PHtoPP((wb,wf,wf1))[2]) - phi_upup_fun_xph(PHtoXPH((wb,wf,wf1))[0], PHtoXPH((wb,wf,wf1))[1], PHtoXPH((wb,wf,wf1))[2])

def Lam_updo_ph(wb,wf,wf1):
    return f_updo_fun_ph(wb,wf,wf1) - phi_updo_fun_ph(wb,wf,wf1) - phi_updo_fun_pp(PHtoPP((wb,wf,wf1))[0],PHtoPP((wb,wf,wf1))[1],PHtoPP((wb,wf,wf1))[2]) - phi_updo_fun_xph(PHtoXPH((wb,wf,wf1))[0], PHtoXPH((wb,wf,wf1))[1], PHtoXPH((wb,wf,wf1))[2])

def Lam_upup_pp(wb,wf,wf1):
    return f_upup_fun_pp(wb,wf,wf1) - phi_upup_fun_pp(wb,wf,wf1) - phi_upup_fun_ph(PPtoPH((wb,wf,wf1))[0],PPtoPH((wb,wf,wf1))[1],PPtoPH((wb,wf,wf1))[2]) - phi_upup_fun_xph(PPtoXPH((wb,wf,wf1))[0], PPtoXPH((wb,wf,wf1))[1], PPtoXPH((wb,wf,wf1))[2])

def Lam_updo_pp(wb,wf,wf1):
    return f_updo_fun_pp(wb,wf,wf1) - phi_updo_fun_pp(wb,wf,wf1) - phi_updo_fun_ph(PPtoPH((wb,wf,wf1))[0],PPtoPH((wb,wf,wf1))[1],PPtoPH((wb,wf,wf1))[2]) - phi_updo_fun_xph(PPtoXPH((wb,wf,wf1))[0], PPtoXPH((wb,wf,wf1))[1], PPtoXPH((wb,wf,wf1))[2])

def Lam_upup_xph(wb,wf,wf1):
    return f_upup_fun_xph(wb,wf,wf1) - phi_upup_fun_xph(wb,wf,wf1) - phi_upup_fun_ph(XPHtoPH((wb,wf,wf1))[0],XPHtoPH((wb,wf,wf1))[1],XPHtoPH((wb,wf,wf1))[2]) - phi_upup_fun_pp(XPHtoPP((wb,wf,wf1))[0], XPHtoPP((wb,wf,wf1))[1], XPHtoPP((wb,wf,wf1))[2])

def Lam_updo_xph(wb,wf,wf1):
    return f_updo_fun_xph(wb,wf,wf1) - phi_updo_fun_xph(wb,wf,wf1) - phi_updo_fun_ph(XPHtoPH((wb,wf,wf1))[0],XPHtoPH((wb,wf,wf1))[1],XPHtoPH((wb,wf,wf1))[2]) - phi_updo_fun_pp(XPHtoPP((wb,wf,wf1))[0], XPHtoPP((wb,wf,wf1))[1], XPHtoPP((wb,wf,wf1))[2])

#----------------------------------------ARRAY CREATION------------------------------------

R_upup_ph_arr = np.array([[[R_upup_ph(wb,wf,wf1) for wb in range(-N_bose, N_bose+1)] for wf in range(-N_fermi, N_fermi)] for wf1 in range(-N_fermi, N_fermi)])
R_updo_ph_arr = np.array([[[R_updo_ph(wb,wf,wf1) for wb in range(-N_bose, N_bose+1)] for wf in range(-N_fermi, N_fermi)] for wf1 in range(-N_fermi, N_fermi)])
R_upup_pp_arr = np.array([[[R_upup_pp(wb,wf,wf1) for wb in range(-N_bose, N_bose+1)] for wf in range(-N_fermi, N_fermi)] for wf1 in range(-N_fermi, N_fermi)])
R_updo_pp_arr = np.array([[[R_updo_pp(wb,wf,wf1) for wb in range(-N_bose, N_bose+1)] for wf in range(-N_fermi, N_fermi)] for wf1 in range(-N_fermi, N_fermi)])
R_upup_xph_arr = np.array([[[R_upup_xph(wb,wf,wf1) for wb in range(-N_bose, N_bose+1)] for wf in range(-N_fermi, N_fermi)] for wf1 in range(-N_fermi, N_fermi)])
R_updo_xph_arr = np.array([[[R_updo_xph(wb,wf,wf1) for wb in range(-N_bose, N_bose+1)] for wf in range(-N_fermi, N_fermi)] for wf1 in range(-N_fermi, N_fermi)])

Lam_upup_ph_arr = np.array([[[Lam_upup_ph(wb,wf,wf1) for wb in range(-N_bose, N_bose+1)] for wf in range(-N_fermi, N_fermi)] for wf1 in range(-N_fermi, N_fermi)])
Lam_updo_ph_arr = np.array([[[Lam_updo_ph(wb,wf,wf1) for wb in range(-N_bose, N_bose+1)] for wf in range(-N_fermi, N_fermi)] for wf1 in range(-N_fermi, N_fermi)])
Lam_upup_pp_arr = np.array([[[Lam_upup_pp(wb,wf,wf1) for wb in range(-N_bose, N_bose+1)] for wf in range(-N_fermi, N_fermi)] for wf1 in range(-N_fermi, N_fermi)])
Lam_updo_pp_arr = np.array([[[Lam_updo_pp(wb,wf,wf1) for wb in range(-N_bose, N_bose+1)] for wf in range(-N_fermi, N_fermi)] for wf1 in range(-N_fermi, N_fermi)])
Lam_upup_xph_arr = np.array([[[Lam_upup_xph(wb,wf,wf1) for wb in range(-N_bose, N_bose+1)] for wf in range(-N_fermi, N_fermi)] for wf1 in range(-N_fermi, N_fermi)])
Lam_updo_xph_arr = np.array([[[Lam_updo_xph(wb,wf,wf1) for wb in range(-N_bose, N_bose+1)] for wf in range(-N_fermi, N_fermi)] for wf1 in range(-N_fermi, N_fermi)])

# HDF5 NEW GROUP CREATION

if (('Lambda' in f) or ('R_func' in f)):
    del f['Lambda']
    del f['R_func']

Lam_grp = f.require_group("Lambda")
Lam_subgrp_ph = Lam_grp.require_group("PH")
Lam_subgrp_pp = Lam_grp.require_group("PP")
Lam_subgrp_xph = Lam_grp.require_group("XPH")

# HDF5 NEW DATASET CREATION

Lam_subgrp_ph.create_dataset("RE_LAMBDA_UPUP", data= Lam_upup_ph_arr.real, dtype='float64', compression="gzip", compression_opts=4)
Lam_subgrp_ph.create_dataset("IM_LAMBDA_UPUP", data= Lam_upup_ph_arr.imag, dtype='float64', compression="gzip", compression_opts=4)
Lam_subgrp_ph.create_dataset("RE_LAMBDA_UPDO", data= Lam_updo_ph_arr.real, dtype='float64', compression="gzip", compression_opts=4)
Lam_subgrp_ph.create_dataset("IM_LAMBDA_UPDO", data= Lam_updo_ph_arr.imag, dtype='float64', compression="gzip", compression_opts=4)
Lam_subgrp_pp.create_dataset("RE_LAMBDA_UPUP", data= Lam_upup_pp_arr.real, dtype='float64', compression="gzip", compression_opts=4)
Lam_subgrp_pp.create_dataset("IM_LAMBDA_UPUP", data= Lam_upup_pp_arr.imag, dtype='float64', compression="gzip", compression_opts=4)
Lam_subgrp_pp.create_dataset("RE_LAMBDA_UPDO", data= Lam_updo_pp_arr.real, dtype='float64', compression="gzip", compression_opts=4)
Lam_subgrp_pp.create_dataset("IM_LAMBDA_UPDO", data= Lam_updo_pp_arr.imag, dtype='float64', compression="gzip", compression_opts=4)
Lam_subgrp_xph.create_dataset("RE_LAMBDA_UPUP", data= Lam_upup_xph_arr.real, dtype='float64', compression="gzip", compression_opts=4)
Lam_subgrp_xph.create_dataset("IM_LAMBDA_UPUP", data= Lam_upup_xph_arr.imag, dtype='float64', compression="gzip", compression_opts=4)
Lam_subgrp_xph.create_dataset("RE_LAMBDA_UPDO", data= Lam_updo_xph_arr.real, dtype='float64', compression="gzip", compression_opts=4)
Lam_subgrp_xph.create_dataset("IM_LAMBDA_UPDO", data= Lam_updo_xph_arr.imag, dtype='float64', compression="gzip", compression_opts=4)


R_grp = f.require_group("R_func")
R_subgrp_ph = R_grp.require_group("PH")
R_subgrp_pp = R_grp.require_group("PP")
R_subgrp_xph = R_grp.require_group("XPH")

# HDF5 NEW DATASET CREATION

R_subgrp_ph.create_dataset("RE_R_UPUP", data= R_upup_ph_arr.real, dtype='float64', compression="gzip", compression_opts=4)
R_subgrp_ph.create_dataset("IM_R_UPUP", data= R_upup_ph_arr.imag, dtype='float64', compression="gzip", compression_opts=4)
R_subgrp_ph.create_dataset("RE_R_UPDO", data= R_updo_ph_arr.real, dtype='float64', compression="gzip", compression_opts=4)
R_subgrp_ph.create_dataset("IM_R_UPDO", data= R_updo_ph_arr.imag, dtype='float64', compression="gzip", compression_opts=4)
R_subgrp_pp.create_dataset("RE_R_UPUP", data= R_upup_pp_arr.real, dtype='float64', compression="gzip", compression_opts=4)
R_subgrp_pp.create_dataset("IM_R_UPUP", data= R_upup_pp_arr.imag, dtype='float64', compression="gzip", compression_opts=4)
R_subgrp_pp.create_dataset("RE_R_UPDO", data= R_updo_pp_arr.real, dtype='float64', compression="gzip", compression_opts=4)
R_subgrp_pp.create_dataset("IM_R_UPDO", data= R_updo_pp_arr.imag, dtype='float64', compression="gzip", compression_opts=4)
R_subgrp_xph.create_dataset("RE_R_UPUP", data= R_upup_xph_arr.real, dtype='float64', compression="gzip", compression_opts=4)
R_subgrp_xph.create_dataset("IM_R_UPUP", data= R_upup_xph_arr.imag, dtype='float64', compression="gzip", compression_opts=4)
R_subgrp_xph.create_dataset("RE_R_UPDO", data= R_updo_xph_arr.real, dtype='float64', compression="gzip", compression_opts=4)
R_subgrp_xph.create_dataset("IM_R_UPDO", data= R_updo_xph_arr.imag, dtype='float64', compression="gzip", compression_opts=4)

f.close()
