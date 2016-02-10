#!/usr/bin/python
#
# This script provides the calculation of the following objects from the ED vertex:
#   - The "Rest" function subtracting from the phi the corresponding plus and karrasch functions
#   - The Fully Irreducible vertex Lam 
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
from agneselib.mymath import *
from agneselib.translate_notation import *
from agneselib.inf_Msum import *
from _functools import partial

#-----------------------------------Read parameters----------------------------------------------
# ------U---------

a=raw_input('Enter the interaction value:') 
try:
    U=float(a)
except ValueError:
    sys.exit("Invalid interaction")

print ("Uhub value     " + str(U) )

#-------BETA-------

a = raw_input('Enter the value of beta:') 
try:
    beta=float(a) 
except ValueError:
    sys.exit("Invalid beta")

print ("beta value     " + str(beta) )

#----------------------------------------Read HDF5 files-----------------------------------------

if ('dat'):
    f = h5py.File('dat/dat_U'+ str(U)+'_beta'+ str(beta)+'_EDpomerol.h5', 'r+')   # Read (and write) the hdf5 file in the directory "dat" if existing
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

# Update the asymptotic structures for the VERTEX IN ALL CHANNELS

def f_upup_fun_ph(i,j,k):
    return re_f_upup_ph[i + N_bose, j+N_fermi, k + N_fermi]+1j*im_f_upup_ph[i + N_bose, j+N_fermi, k + N_fermi]
           
def f_updo_fun_ph(i,j,k):
    return re_f_updo_ph[i + N_bose, j+N_fermi, k + N_fermi]+1j*im_f_updo_ph[i + N_bose, j+N_fermi, k + N_fermi]

def f_upup_fun_pp(i,j,k):
    return re_f_upup_pp[i + N_bose, j+N_fermi, k + N_fermi]+1j*im_f_updo_pp[i + N_bose, j+N_fermi, k + N_fermi]

def f_updo_fun_pp(i,j,k):
    return re_f_updo_pp[i + N_bose, j+N_fermi, k + N_fermi]+1j*im_f_updo_pp[i + N_bose, j+N_fermi, k + N_fermi]

def f_xupdo_fun_pp(i,j,k):
    return re_f_upup_pp[i + N_bose, j+N_fermi, k + N_fermi]+1j*im_f_upup_pp[i + N_bose, j+N_fermi, k + N_fermi]-re_f_updo_pp[i + N_bose, j+N_fermi, k + N_fermi]-1j*im_f_updo_pp[i + N_bose, j+N_fermi, k + N_fermi]

def f_upup_fun_xph(i,j,k):
    return re_f_upup_xph[i + N_bose, j+N_fermi, k + N_fermi]+1j*im_f_upup_xph[i + N_bose, j+N_fermi, k + N_fermi]

def f_updo_fun_xph(i,j,k):
    return re_f_updo_xph[i + N_bose, j+N_fermi, k + N_fermi]+1j*im_f_updo_xph[i + N_bose, j+N_fermi, k + N_fermi]

#--------------------------------GENCHI FUNC--------------------------------------
#PP -> chi (singlet, triplet)

def chi_s_pp(wb,wf,wf1):
    return -re_chi_upup_pp[wb+N_bose,wf+N_fermi,wf1+N_fermi]-1j*im_chi_upup_pp[wb+N_bose,wf+N_fermi,wf1+N_fermi]+ 2*(re_chi_updo_pp[wb+N_bose,wf+N_fermi,wf1+N_fermi] +1j*im_chi_updo_pp[wb+N_bose,wf+N_fermi,wf1+N_fermi])
def chi_t_pp(wb,wf,wf1):
    return re_chi_upup_pp[wb+N_bose,wf+N_fermi,wf1+N_fermi] + 1j*im_chi_upup_pp[wb+N_bose,wf+N_fermi,wf1+N_fermi] 

#PH -> 2pgf(magnetic, density)

def tpgf_d_ph(wb,wf,wf1):
    return re_2pgf_upup_ph[wb+N_bose,wf+N_fermi,wf1+N_fermi]+1j*im_2pgf_upup_ph[wb+N_bose,wf+N_fermi,wf1+N_fermi]+re_2pgf_updo_ph[wb+N_bose,wf+N_fermi,wf1+N_fermi]+1j*im_2pgf_updo_ph[wb+N_bose,wf+N_fermi,wf1+N_fermi]

def tpgf_m_ph(wb,wf,wf1):
    return re_2pgf_upup_ph[wb+N_bose,wf+N_fermi,wf1+N_fermi]+1j*im_2pgf_upup_ph[wb+N_bose,wf+N_fermi,wf1+N_fermi] -re_2pgf_updo_ph[wb+N_bose,wf+N_fermi,wf1+N_fermi]-1j*im_2pgf_updo_ph[wb+N_bose,wf+N_fermi,wf1+N_fermi]

#-----------------------------GAMMA BETHE SALPETER INVERSION


def chis_chi0_arr(wb):
    return np.array([[chi_s_pp(wb,wf,wf1) + chi_0_pp(wb,wf,wf1) for wf in range(-N_fermi,N_fermi)] for wf1 in range(-N_fermi,N_fermi)])

def chit_chi0_arr(wb):
    return np.array([[chi_t_pp(wb,wf,wf1) + chi_0_pp(wb,wf,wf1) for wf in range(-N_fermi,N_fermi)] for wf1 in range(-N_fermi,N_fermi)])

def chi_0_pp_arr(wb):
    return np.array([[chi_0_pp(wb,wf,wf1) for wf in range(-N_fermi,N_fermi)] for wf1 in range(-N_fermi,N_fermi)])

gamma_t_arr = np.array([beta*beta*(2*inv(chi_0_pp_arr(wb))- 4*inv(chit_chi0_arr(wb))) for wb in range (-N_bose,N_bose+1)])
gamma_s_arr = np.array([beta*beta*(2*inv(chi_0_pp_arr(wb)) - 4*inv(chis_chi0_arr(wb))) for wb in range (-N_bose,N_bose+1)])

gamma_upup_pp_arr = gamma_t_arr
gamma_updo_pp_arr = 0.5*(gamma_t_arr + gamma_s_arr) 

#PH

def m_d_arr_ph(wb):
    return np.array([[tpgf_d_ph(wb,wf,wf1)-2*chi_0_ph(wb,wf,wf1) for wf in range(-N_fermi,N_fermi)] for wf1 in range(-N_fermi,N_fermi)])

def m_m_arr_ph(wb):
    return np.array([[tpgf_m_ph(wb,wf,wf1) for wf in range(-N_fermi,N_fermi)] for wf1 in range(-N_fermi,N_fermi)])

def chi_x0_ph_arr(wb):
    return np.array([[chi_x0_ph(wb,wf,wf1) for wf in range(-N_fermi,N_fermi)] for wf1 in range(-N_fermi,N_fermi)])


gamma_d_arr = np.array([beta*beta*(inv(chi_x0_ph_arr(wb))-inv(m_d_arr_ph(wb))) for wb in range (-N_bose,N_bose+1)])
gamma_m_arr = np.array([beta*beta*(inv(chi_x0_ph_arr(wb))-inv(m_m_arr_ph(wb))) for wb in range (-N_bose,N_bose+1)])

gamma_upup_ph_arr = 0.5*(gamma_d_arr + gamma_m_arr)
gamma_updo_ph_arr = 0.5*(gamma_d_arr - gamma_m_arr)

#XPH

def m_updo_fun_xph(wb,wf,wf1):
    return re_2pgf_updo_xph[wb+N_bose,wf+N_fermi,wf1+N_fermi]+1j*im_2pgf_updo_xph[wb+N_bose,wf+N_fermi,wf1+N_fermi]

def m_updo_arr_xph(wb):
    return np.array([[m_updo_fun_xph(wb,wf,wf1) for wf in range(-N_fermi,N_fermi)] for wf1 in range(-N_fermi,N_fermi)])


def chi_0_xph_arr(wb):
    return np.array([[chi_0_xph(wb,wf,wf1) for wf in range(-N_fermi,N_fermi)] for wf1 in range(-N_fermi,N_fermi)])


gamma_upup_xph_arr = - gamma_upup_ph_arr
gamma_updo_xph_arr_big = np.array([beta*beta*(+inv(chi_0_xph_arr(wb))-inv(m_updo_arr_xph(wb))) for wb in range (-N_bose,N_bose+1)])
gamma_updo_xph_arr = gamma_updo_xph_arr_big

#-------------------------------Plotting Vertex-----------------

print "Plotting Vertex extended"

N_fermi_plot = N_fermi

pl.rc('xtick', labelsize=9) 
pl.rc('ytick', labelsize=9) 
cmap = pl.get_cmap('jet') # Bianconiglio

pl.figsize=(13, 7)

shift = 0

def plotVert_ED( use_pl, zarr, string):
    use_pl.set_aspect(1.0)
    pl.pcolormesh(np.array([(2*i+1)*pi/beta for i in range(-N_fermi_plot,N_fermi_plot)]), np.array([(2*i+1)*pi/beta for i in range(-N_fermi_plot,N_fermi_plot)]),np.ma.masked_where( np.isnan(zarr), zarr ))
    use_pl.set_title( string , fontsize=10)    
    pl.colorbar(shrink=0.6)
    return

def plotre_f_upup_ph( use_pl ):
    title=r"$\operatorname{Re}F^{PH}_{\uparrow \uparrow}$"
    zarr = np.array([[ f_upup_fun_ph(shift,n,m).real  for n in range(-N_fermi_plot,N_fermi_plot)] for m in range(-N_fermi_plot,N_fermi_plot)])
    plotVert_ED( use_pl, zarr, title ) 
    return

def plotre_f_updo_ph( use_pl ):
    title=r"$\operatorname{Re}F^{PH}_{\uparrow \downarrow}$"
    zarr = np.array([[ f_updo_fun_ph(shift,n,m).real  for n in range(-N_fermi_plot,N_fermi_plot)] for m in range(-N_fermi_plot,N_fermi_plot)])
    plotVert_ED( use_pl, zarr, title ) 
    return

def plotre_f_upup_pp( use_pl ):
    title=r"$\operatorname{Re}F^{PP}_{\uparrow \uparrow}$"
    zarr = np.array([[ f_upup_fun_pp(shift,n,m).real for n in range(-N_fermi_plot,N_fermi_plot)] for m in range(-N_fermi_plot,N_fermi_plot)])
    plotVert_ED( use_pl, zarr, title ) 
    return

def plotre_f_updo_pp( use_pl ):
    title=r"$\operatorname{Re}F^{PP}_{\uparrow \downarrow}$"
    zarr = np.array([[ f_updo_fun_pp(shift,n,m).real for n in range(-N_fermi_plot,N_fermi_plot)] for m in range(-N_fermi_plot,N_fermi_plot)])
    plotVert_ED( use_pl, zarr, title ) 
    return

def plotre_f_upup_xph( use_pl ):
    title=r"$\operatorname{Re}F^{XPH}_{2\uparrow \uparrow}$"
    zarr = np.array([[ (f_upup_fun_xph(shift,n,m)).real  for n in range(-N_fermi_plot,N_fermi_plot)] for m in range(-N_fermi_plot,N_fermi_plot)])
    plotVert_ED( use_pl, zarr, title ) 
    return

def plotre_f_updo_xph( use_pl ):
    title=r"$\operatorname{Re}F^{XPH}_{2\uparrow \downarrow}$"
    zarr = np.array([[ (f_updo_fun_xph(shift,n,m)).real  for n in range(-N_fermi_plot,N_fermi_plot)] for m in range(-N_fermi_plot,N_fermi_plot)])
    plotVert_ED( use_pl, zarr, title ) 
    return

pl.suptitle(r"$U=$" + str('{0:.3f}'.format(float(U))) + r"     $\beta=$" + str('{0:.3f}'.format(float(beta))) + r"     $\Omega_{\rm PP}=$" + str(shift) + r"$*2\pi/\beta$")

plotre_f_upup_pp(pl.subplot(2,3,1))
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

pl.show()
pl.clf()
#------------------------------- PLOTTING GAMMA PP-----------------------------------

print ("Plotting Gamma PP...")

pl.rc('xtick', labelsize=9) 
pl.rc('ytick', labelsize=9) 
cmap = pl.get_cmap('jet') # Bianconiglio

pl.figsize=(13, 7)

bos = 0

def plotVert_ED( use_pl, zarr, string):
    use_pl.set_aspect(1.0)
    pl.pcolormesh(np.array([(2*i+1)*pi/beta for i in range(-N_fermi,N_fermi)]), np.array([(2*i+1)*pi/beta for i in range(-N_fermi,N_fermi)]),np.ma.masked_where( np.isnan(zarr), zarr ))
    use_pl.set_title( string , fontsize=10)    
    pl.colorbar(shrink=0.6)
    return

def plotre_gamma_upup_pp( use_pl ):
    title=r"$ReGamma^{PP}_{\uparrow \uparrow}$"
    zarr = gamma_upup_pp_arr[N_bose+bos,:,:].real
    plotVert_ED( use_pl, zarr, title ) 
    return

def plotim_gamma_upup_pp( use_pl ):
    title=r"$ImGamma^{PP}_{\uparrow \uparrow}$"
    zarr = gamma_upup_pp_arr[N_bose+bos,:,:].imag
    plotVert_ED( use_pl, zarr, title ) 
    return

def plotre_gamma_updo_pp( use_pl ):
    title=r"$ReGamma^{PP}_{\uparrow \downarrow}$"
    zarr = gamma_updo_pp_arr[N_bose+bos,:,:].real
    plotVert_ED( use_pl, zarr, title ) 
    return

def plotim_gamma_updo_pp( use_pl ):
    title=r"$ImGamma^{PP}_{2\uparrow \downarrow}$"
    zarr = gamma_updo_pp_arr[N_bose+bos,:,:].imag
    plotVert_ED( use_pl, zarr, title ) 
    return

pl.suptitle(r"$U=$" + str('{0:.3f}'.format(float(U))) + r"     $\beta=$" + str('{0:.3f}'.format(float(beta))) + r"     $\Omega_{\rm PP}=$" + str(bos) + r"$*2\pi/\beta$")

plotre_gamma_upup_pp(pl.subplot(2,2,1))
pl.xlabel(r"$\omega_n$", fontsize=20)
pl.ylabel(r"$\omega_m$", fontsize=20)
pl.grid()
plotim_gamma_upup_pp(pl.subplot(2,2,2) )
pl.grid()
plotre_gamma_updo_pp(pl.subplot(2,2,3))
pl.grid()
plotim_gamma_updo_pp(pl.subplot(2,2,4) )
pl.grid()

pl.show()
pl.clf()


#------------------------------- PLOTTING GAMMA PH -----------------------------------

print ("Plotting Gamma PH...")

pl.rc('xtick', labelsize=9) 
pl.rc('ytick', labelsize=9) 
cmap = pl.get_cmap('jet') # Bianconiglio

pl.figsize=(13, 7)

bos = 0

def plotVert_ED( use_pl, zarr, string):
    use_pl.set_aspect(1.0)
    pl.pcolormesh(np.array([(2*i+1)*pi/beta for i in range(-N_fermi,N_fermi)]), np.array([(2*i+1)*pi/beta for i in range(-N_fermi,N_fermi)]),np.ma.masked_where( np.isnan(zarr), zarr ))
    use_pl.set_title( string , fontsize=10)    
    pl.colorbar(shrink=0.6)
    return

def plotre_gamma_upup_ph( use_pl ):
    title=r"$ReGamma^{PH}_{\uparrow \uparrow}$"
    zarr = gamma_upup_ph_arr[N_bose+bos,:,:].real
    plotVert_ED( use_pl, zarr, title ) 
    return

def plotim_gamma_upup_ph( use_pl ):
    title=r"$ImGamma^{PH}_{\uparrow \uparrow}$"
    zarr = gamma_upup_ph_arr[N_bose+bos,:,:].imag
    plotVert_ED( use_pl, zarr, title ) 
    return

def plotre_gamma_updo_ph( use_pl ):
    title=r"$ReGamma^{PH}_{\uparrow \downarrow}$"
    zarr = gamma_updo_ph_arr[N_bose+bos,:,:].real
    plotVert_ED( use_pl, zarr, title ) 
    return

def plotim_gamma_updo_ph( use_pl ):
    title=r"$ImGamma^{PH}_{2\uparrow \downarrow}$"
    zarr = gamma_updo_ph_arr[N_bose+bos,:,:].imag
    plotVert_ED( use_pl, zarr, title ) 
    return


pl.suptitle(r"$U=$" + str('{0:.3f}'.format(float(U))) + r"     $\beta=$" + str('{0:.3f}'.format(float(beta))) + r"     $\Omega_{\rm PH}=$" + str(bos) + r"$*2\pi/\beta$")

plotre_gamma_upup_ph(pl.subplot(2,2,1))
pl.xlabel(r"$\omega_n$", fontsize=20)
pl.ylabel(r"$\omega_m$", fontsize=20)
pl.grid()
plotim_gamma_upup_ph(pl.subplot(2,2,2) )
pl.grid()
plotre_gamma_updo_ph(pl.subplot(2,2,3))
pl.grid()
plotim_gamma_updo_ph(pl.subplot(2,2,4) )
pl.grid()

pl.show()
pl.clf()


#------------------------------- PLOTTING GAMMA XPH -----------------------------------

print ("Plotting Gamma XPH...")

pl.rc('xtick', labelsize=9) 
pl.rc('ytick', labelsize=9) 
cmap = pl.get_cmap('jet') # Bianconiglio

pl.figsize=(13, 7)

bos = 0

def plotVert_ED( use_pl, zarr, string):
    use_pl.set_aspect(1.0)
    pl.pcolormesh(np.array([(2*i+1)*pi/beta for i in range(-N_fermi,N_fermi)]), np.array([(2*i+1)*pi/beta for i in range(-N_fermi,N_fermi)]),np.ma.masked_where( np.isnan(zarr), zarr ))
    use_pl.set_title( string , fontsize=10)    
    pl.colorbar(shrink=0.6)
    return

def plotre_gamma_upup_xph( use_pl ):
    title=r"$ReGamma^{XPH}_{\uparrow \uparrow}$"
    zarr = gamma_upup_xph_arr[N_bose+bos,:,:].real
    plotVert_ED( use_pl, zarr, title ) 
    return

def plotim_gamma_upup_xph( use_pl ):
    title=r"$ImGamma^{XPH}_{\uparrow \uparrow}$"
    zarr = gamma_upup_xph_arr[N_bose+bos,:,:].imag
    plotVert_ED( use_pl, zarr, title ) 
    return

def plotre_gamma_updo_xph( use_pl ):
    title=r"$ReGamma^{XPH}_{\uparrow \downarrow}$"
    zarr = gamma_updo_xph_arr[N_bose+bos,:,:].real
    plotVert_ED( use_pl, zarr, title ) 
    return

def plotim_gamma_updo_xph( use_pl ):
    title=r"$ImGamma^{XPH}_{2\uparrow \downarrow}$"
    zarr = gamma_updo_xph_arr[N_bose+bos,:,:].imag
    plotVert_ED( use_pl, zarr, title ) 
    return

pl.suptitle(r"$U=$" + str('{0:.3f}'.format(float(U))) + r"     $\beta=$" + str('{0:.3f}'.format(float(beta))) + r"     $\Omega_{\rm XPH}=$" + str(bos) + r"$*2\pi/\beta$")

plotre_gamma_upup_xph(pl.subplot(2,2,1))
pl.xlabel(r"$\omega_n$", fontsize=20)
pl.ylabel(r"$\omega_m$", fontsize=20)
pl.grid()
plotim_gamma_upup_xph(pl.subplot(2,2,2) )
pl.grid()
plotre_gamma_updo_xph(pl.subplot(2,2,3))
pl.grid()
plotim_gamma_updo_xph(pl.subplot(2,2,4) )
pl.grid()

pl.show()
pl.clf()


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
    return f_upup_fun_ph(i,j,k)-gamma_upup_ph_arr[i + N_bose, j+N_fermi, k + N_fermi] 
    
def phi_updo_fun_ph(i,j,k):
    return f_updo_fun_ph(i,j,k)-gamma_updo_ph_arr[i + N_bose, j+N_fermi, k + N_fermi]
    
def phi_upup_fun_pp(i,j,k):
    return f_upup_fun_pp(i,j,k)-gamma_upup_pp_arr[i + N_bose, j+N_fermi, k + N_fermi]
    
def phi_updo_fun_pp(i,j,k):
    return f_updo_fun_pp(i,j,k)-gamma_updo_pp_arr[i + N_bose, j+N_fermi, k + N_fermi] 
    
def phi_upup_fun_xph(i,j,k):
    return f_upup_fun_xph(i,j,k)-gamma_upup_xph_arr[i + N_bose, j+N_fermi, k + N_fermi]
    
def phi_updo_fun_xph(i,j,k):
    return f_updo_fun_xph(i,j,k)-gamma_updo_xph_arr[i + N_bose, j+N_fermi, k + N_fermi]

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

#---------------------------KARRASCH--------------------------

def K_upup_ph(wb):
    if (abs(wb) <= N_bose):
        return phi_upup_ph_arr[wb+N_bose,1,1]
    else:
        return 0.0

def K_updo_ph(wb):
    if (abs(wb) <= N_bose):
        return phi_updo_ph_arr[wb+N_bose,1,1]
    else:
        return 0.0

def K_upup_pp(wb):
    if (abs(wb) <= N_bose):
        return phi_upup_pp_arr[wb+N_bose,1,1]
    else:
        return 0.0

def K_updo_pp(wb):
    if (abs(wb) <= N_bose):
        return phi_updo_pp_arr[wb+N_bose,1,1]
    else:
        return 0.0

def K_upup_xph(wb):
    if (abs(wb) <= N_bose):
        return phi_upup_xph_arr[wb+N_bose,1,1]
    else:
        return 0.0

def K_updo_xph(wb):
    if (abs(wb) <= N_bose):
        return phi_updo_xph_arr[wb+N_bose,1,1]
    else:
        return 0.0

#---------------------------P func--------------------------------

def P_upup_ph(wb,wf):
    if (abs(wb) <= N_bose and wf >= -N_fermi and wf < N_fermi):
        return  phi_upup_fun_ph(wb,wf,-N_fermi+1)- K_upup_ph(wb)
    else:
#        print "P_upup_ph out"
        return 0.0
def P_updo_ph(wb,wf):
    if (abs(wb) <= N_bose and wf >= -N_fermi and wf < N_fermi):
        return phi_updo_fun_ph(wb,wf,-N_fermi+1)- K_updo_ph(wb)
    else:
#        print "P_updo_ph out"
        return 0.0
def P_upup_pp(wb,wf):
    if (abs(wb) <= N_bose and wf >= -N_fermi and wf < N_fermi):
        return phi_upup_fun_pp(wb,wf,-N_fermi+1)- K_upup_pp(wb)
    else:
#        print "P_updo_pp out"
        return 0.0
def P_updo_pp(wb,wf):
    if (abs(wb) <= N_bose and wf >= -N_fermi and wf < N_fermi):
        return phi_updo_fun_pp(wb,wf,-N_fermi+1)- K_updo_pp(wb)
    else:
#        print "P_updo_pp out"
        return 0.0
def P_upup_xph(wb,wf):
    if (abs(wb) <= N_bose and wf >= -N_fermi and wf < N_fermi):
        return phi_upup_fun_xph(wb,wf,-N_fermi+1)- K_upup_xph(wb)
    else:
#        print "P_upup_xph out"
        return 0.0

def P_updo_xph(wb,wf):
    if (abs(wb) <= N_bose and wf >= -N_fermi and wf < N_fermi):
        return phi_updo_fun_xph(wb,wf,-N_fermi+1)- K_updo_xph(wb)
    else:
#        print "K_updo_xph out"
        return 0.0
 
 #Save Karrasch and Plus extarcted via scanning

K_upup_ph_arr = np.array([K_upup_ph(i) for i in range(-N_bose, N_bose+1)])
K_updo_ph_arr = np.array([K_updo_ph(i) for i in range(-N_bose, N_bose+1)])
K_upup_pp_arr = np.array([K_upup_pp(i) for i in range(-N_bose, N_bose+1)])
K_updo_pp_arr = np.array([K_updo_pp(i) for i in range(-N_bose, N_bose+1)])
K_upup_xph_arr = np.array([K_upup_xph(i) for i in range(-N_bose, N_bose+1)])
K_updo_xph_arr = np.array([K_updo_xph(i) for i in range(-N_bose, N_bose+1)])

P_upup_ph_arr = np.array([[P_upup_ph(i,j) for j in range (-N_fermi, N_fermi)] for i in range(-N_bose, N_bose+1)])
P_updo_ph_arr = np.array([[P_updo_ph(i,j) for j in range (-N_fermi, N_fermi)] for i in range(-N_bose, N_bose+1)])
P_upup_pp_arr = np.array([[P_upup_pp(i,j) for j in range (-N_fermi, N_fermi)] for i in range(-N_bose, N_bose+1)])
P_updo_pp_arr = np.array([[P_updo_pp(i,j) for j in range (-N_fermi, N_fermi)] for i in range(-N_bose, N_bose+1)])
P_upup_xph_arr = np.array([[P_upup_xph(i,j) for j in range (-N_fermi, N_fermi)] for i in range(-N_bose, N_bose+1)])
P_updo_xph_arr = np.array([[P_updo_xph(i,j) for j in range (-N_fermi, N_fermi)] for i in range(-N_bose, N_bose+1)])

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

p_subgrp_ph.create_dataset("RE_P_UPUP", data= P_upup_ph_arr.real, dtype='float64', compression="gzip", compression_opts=4)
p_subgrp_ph.create_dataset("IM_P_UPUP", data= P_upup_ph_arr.imag, dtype='float64', compression="gzip", compression_opts=4)
p_subgrp_ph.create_dataset("RE_P_UPDO", data= P_updo_ph_arr.real, dtype='float64', compression="gzip", compression_opts=4)
p_subgrp_ph.create_dataset("IM_P_UPDO", data= P_updo_ph_arr.imag, dtype='float64', compression="gzip", compression_opts=4)
p_subgrp_pp.create_dataset("RE_P_UPUP", data= P_upup_pp_arr.real, dtype='float64', compression="gzip", compression_opts=4)
p_subgrp_pp.create_dataset("IM_P_UPUP", data= P_upup_pp_arr.imag, dtype='float64', compression="gzip", compression_opts=4)
p_subgrp_pp.create_dataset("RE_P_UPDO", data= P_updo_pp_arr.real, dtype='float64', compression="gzip", compression_opts=4)
p_subgrp_pp.create_dataset("IM_P_UPDO", data= P_updo_pp_arr.imag, dtype='float64', compression="gzip", compression_opts=4)
p_subgrp_xph.create_dataset("RE_P_UPUP", data= P_upup_xph_arr.real, dtype='float64', compression="gzip", compression_opts=4)
p_subgrp_xph.create_dataset("IM_P_UPUP", data= P_upup_xph_arr.imag, dtype='float64', compression="gzip", compression_opts=4)
p_subgrp_xph.create_dataset("RE_P_UPDO", data= P_updo_xph_arr.real, dtype='float64', compression="gzip", compression_opts=4)
p_subgrp_xph.create_dataset("IM_P_UPDO", data= P_updo_xph_arr.imag, dtype='float64', compression="gzip", compression_opts=4)

k_subgrp_ph.create_dataset("RE_K_UPUP", data= K_upup_ph_arr.real, dtype='float64', compression="gzip", compression_opts=4)
k_subgrp_ph.create_dataset("IM_K_UPUP", data= K_upup_ph_arr.imag, dtype='float64', compression="gzip", compression_opts=4)
k_subgrp_ph.create_dataset("RE_K_UPDO", data= K_updo_ph_arr.real, dtype='float64', compression="gzip", compression_opts=4)
k_subgrp_ph.create_dataset("IM_K_UPDO", data= K_updo_ph_arr.imag, dtype='float64', compression="gzip", compression_opts=4)
k_subgrp_pp.create_dataset("RE_K_UPUP", data= K_upup_pp_arr.real, dtype='float64', compression="gzip", compression_opts=4)
k_subgrp_pp.create_dataset("IM_K_UPUP", data= K_upup_pp_arr.imag, dtype='float64', compression="gzip", compression_opts=4)
k_subgrp_pp.create_dataset("RE_K_UPDO", data= K_updo_pp_arr.real, dtype='float64', compression="gzip", compression_opts=4)
k_subgrp_pp.create_dataset("IM_K_UPDO", data= K_updo_pp_arr.imag, dtype='float64', compression="gzip", compression_opts=4)
k_subgrp_xph.create_dataset("RE_K_UPUP", data= K_upup_xph_arr.real, dtype='float64', compression="gzip", compression_opts=4)
k_subgrp_xph.create_dataset("IM_K_UPUP", data= K_upup_xph_arr.imag, dtype='float64', compression="gzip", compression_opts=4)
k_subgrp_xph.create_dataset("RE_K_UPDO", data= K_updo_xph_arr.real, dtype='float64', compression="gzip", compression_opts=4)
k_subgrp_xph.create_dataset("IM_K_UPDO", data= K_updo_xph_arr.imag, dtype='float64', compression="gzip", compression_opts=4)

#---------------------------------------NOW EXTEND PLUS AND F WITH THE EXTRACTED ASYPTOTICS -------------------------------

def f_ext_upup_ph(i,j,k):
    if IsInside(i,j,k):
        return f_upup_ph(i,j,k)
    else:
        return K_upup_ph(i) + P_upup_ph(i,j)+ P_upup_ph(i,k) + K_upup_xph(PHtoXPH((i,j,k))[0]) + P_upup_xph(PHtoXPH((i,j,k))[0],PHtoXPH((i,j,k))[1])+P_upup_xph(PHtoXPH((i,j,k))[0],PHtoXPH((i,j,k))[2])
           
def f_ext_updo_ph(i,j,k):
    if IsInside(i,j,k):
        return f_updo_ph(i,j,k)
    else:
        return - U + K_updo_ph(i) + P_updo_ph(i,j)+P_updo_ph(i,k) + K_updo_xph(PHtoXPH((i,j,k))[0]) + P_updo_xph(PHtoXPH((i,j,k))[0],PHtoXPH((i,j,k))[1]) + P_updo_xph(PHtoXPH((i,j,k))[0],PHtoXPH((i,j,k))[2])+ K_updo_pp(PHtoPP((i,j,k))[0]) + P_updo_pp(PHtoPP((i,j,k))[0],PHtoPP((i,j,k))[1]) + P_updo_pp(PHtoPP((i,j,k))[0],PHtoPP((i,j,k))[2])

def f_ext_upup_pp(i,j,k):
    if IsInside(i,j,k):
        return f_upup_pp(i,j,k)
    else:
        return  K_upup_pp(i) + P_upup_pp(i,j)+ P_upup_pp(i,k)+K_upup_ph(PPtoPH((i,j,k))[0]) + P_upup_ph(PPtoPH((i,j,k))[0],PPtoPH((i,j,k))[1]) + P_upup_ph(PPtoPH((i,j,k))[0],PPtoPH((i,j,k))[2])+ K_upup_xph(PPtoXPH((i,j,k))[0]) + P_upup_xph(PPtoXPH((i,j,k))[0],PPtoXPH((i,j,k))[1]) + P_upup_xph(PPtoXPH((i,j,k))[0],PPtoXPH((i,j,k))[2])

def f_ext_updo_pp(i,j,k):
    if IsInside(i,j,k):
        return f_updo_pp(i,j,k)
    else:
        return - U + K_updo_pp(i) + P_updo_pp(i,j) + P_updo_pp(i,k) + K_updo_ph(PPtoPH((i,j,k))[0]) + P_updo_ph(PPtoPH((i,j,k))[0],PPtoPH((i,j,k))[1]) + P_updo_ph(PPtoPH((i,j,k))[0],PPtoPH((i,j,k))[2])+ K_updo_xph(PPtoXPH((i,j,k))[0]) + P_updo_xph(PPtoXPH((i,j,k))[0],PPtoXPH((i,j,k))[1]) + P_updo_xph(PPtoXPH((i,j,k))[0],PPtoXPH((i,j,k))[2])

def f_ext_upup_xph(i,j,k):
    if IsInside(i,j,k):
        return f_upup_xph(i,j,k)
    else:
        return K_upup_xph(i) + P_upup_xph(i,j)+P_upup_xph(i,k) + K_upup_ph(XPHtoPH((i,j,k))[0]) + P_upup_ph(XPHtoPH((i,j,k))[0],XPHtoPH((i,j,k))[1])+ P_upup_ph(XPHtoPH((i,j,k))[0],XPHtoPH((i,j,k))[2]) 

def f_ext_updo_xph(i,j,k):
    if IsInside(i,j,k):
        return f_updo_xph(i,j,k)
    else:
        return - U + K_updo_xph(i) + P_updo_xph(i,j)+P_updo_xph(i,k) + K_updo_ph(XPHtoPH((i,j,k))[0]) + P_updo_ph(XPHtoPH((i,j,k))[0],XPHtoPH((i,j,k))[1]) + P_updo_ph(XPHtoPH((i,j,k))[0],XPHtoPH((i,j,k))[2])+ K_updo_pp(XPHtoPP((i,j,k))[0]) + P_updo_pp(XPHtoPP((i,j,k))[0],XPHtoPP((i,j,k))[1]) + P_updo_pp(XPHtoPP((i,j,k))[0],XPHtoPP((i,j,k))[2])

#PHIs

def phi_ext_upup_fun_ph(i,j,k):
    if isInside(i,j,k):
        return phi_upup_fun_ph(i,j,k) 
    else:
       return K_upup_ph(i) + P_upup_ph(i,j) + P_upup_ph(i,k) 
    
def phi_ext_updo_fun_ph(i,j,k):
    if isInside(i,j,k):
        return phi_updo_fun_ph(i,j,k)
    else:
        return K_updo_ph(i) + P_updo_ph(i,j)+P_updo_ph(i,k)
    
def phi_ext_upup_fun_pp(i,j,k):
    if isInside(i,j,k):
        return phi_upup_fun_pp(i,j,k)
    else:
        return K_upup_pp(i) + P_upup_pp(i,j)+P_upup_pp(i,k)
    
def phi_ext_updo_fun_pp(i,j,k):
    if isInside(i,j,k):
        return phi_updo_fun_pp(i,j,k) 
    else:
        return K_updo_pp(i) + P_updo_pp(i,j) + P_updo_pp(i,k)
    
def phi_ext_upup_fun_xph(i,j,k):
    if isInside(i,j,k):
        return phi_upup_fun_xph(i,j,k)
    else:
        return K_upup_xph(i) + P_upup_xph(i,j)  + P_upup_xph(i,k)
    
def phi_ext_updo_fun_xph(i,j,k):
    if isInside(i,j,k):
        return phi_updo_fun_xph(i,j,k) 
    else:
        return K_updo_xph(i) + P_updo_xph(i,j)+ P_updo_xph(i,k)

#-----------------------------------REST FUNCTION----------------------------------------

def R_upup_ph(i,j,k):
    return phi_ext_upup_fun_ph(i,j,k)-P_upup_ph(i,j)-P_upup_ph(i,k)-K_upup_ph(i)

def R_updo_ph(i,j,k):
    return phi_ext_updo_fun_ph(i,j,k)-P_updo_ph(i,j)-P_updo_ph(i,k)-K_updo_ph(i)

def R_upup_pp(i,j,k):
    return phi_ext_upup_fun_pp(i,j,k)-P_upup_pp(i,j)-P_upup_pp(i,k)-K_upup_pp(i)

def R_updo_pp(i,j,k):
    return phi_ext_updo_fun_pp(i,j,k)-P_updo_pp(i,j)-P_updo_pp(i,k)-K_updo_pp(i)

def R_upup_xph(i,j,k):
    return phi_ext_upup_fun_xph(i,j,k)-P_upup_xph(i,j)-P_upup_xph(i,k)-K_upup_xph(i)

def R_updo_xph(i,j,k):
    return phi_ext_updo_fun_xph(i,j,k)-P_updo_xph(i,j)-P_updo_xph(i,k)-K_updo_xph(i)

#---------------------------------------LAMBDA--------------------------------------------

def Lam_upup_ph(wb,wf,wf1):
    return f_upup_fun_ph(wb,wf,wf1) - phi_ext_upup_fun_ph(wb,wf,wf1) - phi_ext_upup_fun_pp(PHtoPP((wb,wf,wf1))[0],PHtoPP((wb,wf,wf1))[1],PHtoPP((wb,wf,wf1))[2]) - phi_ext_upup_fun_xph(PHtoXPH((wb,wf,wf1))[0], PHtoXPH((wb,wf,wf1))[1], PHtoXPH((wb,wf,wf1))[2])

def Lam_updo_ph(wb,wf,wf1):
    return f_updo_fun_ph(wb,wf,wf1) - phi_ext_updo_fun_ph(wb,wf,wf1) - phi_ext_updo_fun_pp(PHtoPP((wb,wf,wf1))[0],PHtoPP((wb,wf,wf1))[1],PHtoPP((wb,wf,wf1))[2]) - phi_ext_updo_fun_xph(PHtoXPH((wb,wf,wf1))[0], PHtoXPH((wb,wf,wf1))[1], PHtoXPH((wb,wf,wf1))[2])

def Lam_upup_pp(wb,wf,wf1):
    return f_upup_fun_pp(wb,wf,wf1) - phi_ext_upup_fun_pp(wb,wf,wf1) - phi_ext_upup_fun_ph(PPtoPH((wb,wf,wf1))[0],PPtoPH((wb,wf,wf1))[1],PPtoPH((wb,wf,wf1))[2]) - phi_ext_upup_fun_xph(PPtoXPH((wb,wf,wf1))[0], PPtoXPH((wb,wf,wf1))[1], PPtoXPH((wb,wf,wf1))[2])

def Lam_updo_pp(wb,wf,wf1):
    return f_updo_fun_pp(wb,wf,wf1) - phi_ext_updo_fun_pp(wb,wf,wf1) - phi_ext_updo_fun_ph(PPtoPH((wb,wf,wf1))[0],PPtoPH((wb,wf,wf1))[1],PPtoPH((wb,wf,wf1))[2]) - phi_ext_updo_fun_xph(PPtoXPH((wb,wf,wf1))[0], PPtoXPH((wb,wf,wf1))[1], PPtoXPH((wb,wf,wf1))[2])

def Lam_upup_xph(wb,wf,wf1):
    return f_upup_fun_xph(wb,wf,wf1) - phi_ext_upup_fun_xph(wb,wf,wf1) - phi_ext_upup_fun_ph(XPHtoPH((wb,wf,wf1))[0],XPHtoPH((wb,wf,wf1))[1],XPHtoPH((wb,wf,wf1))[2]) - phi_ext_upup_fun_pp(XPHtoPP((wb,wf,wf1))[0], XPHtoPP((wb,wf,wf1))[1], XPHtoPP((wb,wf,wf1))[2])

def Lam_updo_xph(wb,wf,wf1):
    return f_updo_fun_xph(wb,wf,wf1) - phi_ext_updo_fun_xph(wb,wf,wf1) - phi_ext_updo_fun_ph(XPHtoPH((wb,wf,wf1))[0],XPHtoPH((wb,wf,wf1))[1],XPHtoPH((wb,wf,wf1))[2]) - phi_ext_updo_fun_pp(XPHtoPP((wb,wf,wf1))[0], XPHtoPP((wb,wf,wf1))[1], XPHtoPP((wb,wf,wf1))[2])

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
