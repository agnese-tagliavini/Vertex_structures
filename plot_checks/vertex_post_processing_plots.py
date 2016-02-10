#!/usr/bin/python
#------------------------------------------------------------------------------------
# CONTENTS:
# The script loads and shiftes the vert_centralized back in order to compare with the original one.
# From a vimdiff comparison the two vertices should be the same
#
#------------------------------------IMPORTS------------------------------------

import numpy as np
import matplotlib.pyplot as pl
import os
import sys
import subprocess
import cmath
import math
from os.path import isfile, join
import matplotlib.cm as cm
from numpy.linalg import inv
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from agneselib.mymath import *                          # mylibrary check in ~/usr/include/agneselib
from agneselib.inf_Msum import * 
from _functools import partial
print myfloor_div2(3)
#---------------------------------------------------------------------------------

U=1.0

beta=10.0

pi = math.pi

print math.log10(1.0)
#---------------------------------------------------------------------------------

def run(command):
        output = subprocess.check_output(command, shell=True)
        return output

#---------------------------------------------------------------------------------
#GF

g_iw  = np.loadtxt('../pomerol_output/U_1.0_beta_10.0_2bs/gw_imag00.dat')
N_fermi_gf = g_iw.shape[0]
print ("Number of fermionic frequencies for the GF:")
print N_fermi_gf

def G(w):                           # imaginary part of the GF
 if (w >= 0):
     return g_iw[w,1]+1j*g_iw[w,2]
 else:
     return g_iw[-w-1,1]-1j*g_iw[-w-1,2]

#--------------------------------------- 2PGF PP------------------------------------------------------------

vertex_pp = np.loadtxt("../pomerol_output/U_1.0_beta_10.0_2bs/2pgf_pp_shift.dat")
print vertex_pp.shape
ffreq_original_pp = int(np.transpose(vertex_pp)[1,:].max()+1)
print ffreq_original_pp
bfreq_original_pp = int(np.transpose(vertex_pp)[0,:].max())
print bfreq_original_pp

Nl =3
#wNum = int(ffreq_original_pp/100.0*10.0)
wNum = 1
iMin = ffreq_original_pp-wNum

inv_fit = generate_invfit_func(iMin,wNum,Nl)

def re_2pgf_upup_pp(wb,wf,wf1):
	return 	vertex_pp[(2*ffreq_original_pp)*(2*ffreq_original_pp)*(wb+bfreq_original_pp)+(2*ffreq_original_pp)*(wf+ffreq_original_pp) + (wf1+ffreq_original_pp), 3]

def re_2pgf_updo_pp(wb,wf,wf1):
	return 	vertex_pp[(2*ffreq_original_pp)*(2*ffreq_original_pp)*(wb+bfreq_original_pp)+(2*ffreq_original_pp)*(wf+ffreq_original_pp)+ (wf1+ffreq_original_pp), 5]

def im_2pgf_upup_pp(wb,wf,wf1):
	return 	vertex_pp[(2*ffreq_original_pp)*(2*ffreq_original_pp)*(wb+bfreq_original_pp)+(2*ffreq_original_pp)*(wf+ffreq_original_pp) + (wf1+ffreq_original_pp), 4]

def im_2pgf_updo_pp(wb,wf,wf1):
	return 	vertex_pp[(2*ffreq_original_pp)*(2*ffreq_original_pp)*(wb+bfreq_original_pp)+(2*ffreq_original_pp)*(wf+ffreq_original_pp)+ (wf1+ffreq_original_pp), 6]

# ----------define arrays to store in hdf5 file
#We need to throw away wNum frequencies from the original 2pgf (because of the gamma fitting) -> we store all quantities in this new range 

ffreq_pp = ffreq_original_pp-wNum
bfreq_pp = bfreq_original_pp

# Cutting one disconnected diagram the crossed one -> generalized susceptibility

def chi_upup_pp(wb,wf,wf1):
    if (wf == -wf1-1-mymod_abs(wb)):
        return re_2pgf_upup_pp(wb,wf,wf1)+1j*im_2pgf_upup_pp(wb,wf,wf1)+beta*(G(wf+myceil_div2(wb))*G(myfloor_div2(wb)-wf-1))
    else:
        return re_2pgf_upup_pp(wb,wf,wf1)+1j*im_2pgf_upup_pp(wb,wf,wf1)

def chi_updo_pp(wb,wf,wf1):
    return re_2pgf_updo_pp(wb,wf,wf1)+1j*im_2pgf_updo_pp(wb,wf,wf1)

def chi_s_pp(wb,wf,wf1):
    return -chi_upup_pp(wb,wf,wf1) + 2*chi_updo_pp(wb,wf,wf1)

def chi_t_pp(wb,wf,wf1):
    return chi_upup_pp(wb,wf,wf1)


#----------------------------------- FULL VERTEX PP -----------------------------------------------------------


def chi_l_pp(wb,wf,wf1):  # Legs on one side of the diagram
    return beta*G(wf+myceil_div2(wb))*G(myfloor_div2(wb)-wf1-1)

def chi_0_pp(wb,wf,wf1):
    if (wf == wf1):
        return beta*G(wf+myceil_div2(wb))*G(myfloor_div2(wb)-wf-1)
    else:
        return 0.0

def f_upup_pp(wb,wf,wf1):
    return beta*beta*(1.0/(chi_l_pp(wb,wf,wf)))*(chi_upup_pp(wb,wf,wf1)-chi_0_pp(wb,wf,wf1))*(1.0/(chi_l_pp(wb,wf1,wf1)))

def f_updo_pp(wb,wf,wf1):
    return beta*beta*(1.0/(chi_l_pp(wb,wf,wf)))*(chi_updo_pp(wb,wf,wf1)-chi_0_pp(wb,wf,wf1))*(1.0/(chi_l_pp(wb,wf1,wf1)))

#---------------------------------- 2PI VERTEX GAMMA PP ----------------------------------------------------------

def chis_chi0_arr(wb,wf,wf1):
    return chi_s_pp(wb,wf,wf1) + chi_0_pp(wb,wf,wf1)


def chit_chi0_arr(wb,wf,wf1):
    return chi_t_pp(wb,wf,wf1) + chi_0_pp(wb,wf,wf1)


def chi_0_pp_arr(wb,wf,wf1):
    return chi_0_pp(wb,wf,wf1)

gamma_t_arr = np.array([beta*beta*(2*inv_fit(partial(chi_0_pp_arr,wb)) - 4*inv_fit(partial(chit_chi0_arr,wb))) for wb in range (-bfreq_pp,bfreq_pp+1)])
gamma_s_arr = np.array([beta*beta*(2*inv_fit(partial(chi_0_pp_arr,wb)) - 4*inv_fit(partial(chis_chi0_arr,wb))) for wb in range (-bfreq_pp,bfreq_pp+1)])

gamma_upup_pp_arr = gamma_t_arr 
gamma_updo_pp_arr = 0.5*(gamma_t_arr + gamma_s_arr) 
print gamma_upup_pp_arr.shape
print  gamma_updo_pp_arr[0+bfreq_pp,:,:].max()

#----------------------------------2PR VERTEX PHI----------------------------------------------------------------------------------------

f_upup_pp_arr = np.array([[[f_upup_pp(wb,wf,wf1) for wf in range (-ffreq_pp,ffreq_pp)] for wf1 in range(-ffreq_pp,ffreq_pp)] for wb in range (-bfreq_pp,bfreq_pp+1)])
f_updo_pp_arr = np.array([[[f_updo_pp(wb,wf,wf1) for wf in range (-ffreq_pp,ffreq_pp)] for wf1 in range(-ffreq_pp,ffreq_pp)] for wb in range (-bfreq_pp,bfreq_pp+1)])

print f_upup_pp_arr.shape


phi_upup_pp_arr = f_upup_pp_arr - gamma_upup_pp_arr
phi_updo_pp_arr = f_updo_pp_arr - gamma_updo_pp_arr

#------------------------------- PLOTTING 2PGF -----------------------------------

print ("Plotting generalized 2PGF...")

pl.rc('xtick', labelsize=9) 
pl.rc('ytick', labelsize=9) 
cmap = pl.get_cmap('jet') # Bianconiglio

pl.figsize=(13, 7)

bos = 0

def plotVert_ED( use_pl, zarr, string):
    use_pl.set_aspect(1.0)
    pl.pcolormesh(np.array([(2*i+1)*pi/beta for i in range(-ffreq_pp,ffreq_pp)]), np.array([(2*i+1)*pi/beta for i in range(-ffreq_pp,ffreq_pp)]),np.ma.masked_where( np.isnan(zarr), zarr ))
    use_pl.set_title( string , fontsize=10)    
    pl.colorbar(shrink=0.6)
    return

def plotVert_ED1( use_pl, zarr, string):
    use_pl.set_aspect(1.0)
    pl.pcolormesh(np.array([(2*i+1)*pi/beta for i in range(-ffreq,ffreq)]), np.array([(2*i+1)*pi/beta for i in range(-ffreq,ffreq)]),np.ma.masked_where( np.isnan(zarr), zarr ))
    use_pl.set_title( string , fontsize=10)    
    pl.colorbar(shrink=0.6)
    return

def plotre_2pgf_upup_pp( use_pl ):
    title=r"$ReG^{PP}_{2\uparrow \uparrow}$"
    zarr = np.array([[ re_2pgf_upup_pp(bos,n,m)  for n in range(-ffreq_pp,ffreq_pp)] for m in range(-ffreq_pp,ffreq_pp)])
    plotVert_ED( use_pl, zarr, title ) 
    return

def plotim_2pgf_upup_pp( use_pl ):
    title=r"$ImG^{PP}_{2\uparrow \uparrow}$"
    zarr = np.array([[ im_2pgf_upup_pp(bos,n,m)  for n in range(-ffreq_pp,ffreq_pp)] for m in range(-ffreq_pp,ffreq_pp)])
    plotVert_ED( use_pl, zarr, title ) 
    return

def plotre_2pgf_updo_pp( use_pl ):
    title=r"$ReG^{PP}_{2\uparrow \downarrow}$"
    zarr = np.array([[ re_2pgf_updo_pp(bos,n,m)  for n in range(-ffreq_pp,ffreq_pp)] for m in range(-ffreq_pp,ffreq_pp)])
    plotVert_ED( use_pl, zarr, title ) 
    return

def plotim_2pgf_updo_pp( use_pl ):
    title=r"$ImG^{PP}_{2\uparrow \downarrow}$"
    zarr = np.array([[ im_2pgf_updo_pp(bos,n,m)  for n in range(-ffreq_pp,ffreq_pp)] for m in range(-ffreq_pp,ffreq_pp)])
    plotVert_ED( use_pl, zarr, title ) 
    return

pl.suptitle(r"$U=$" + str('{0:.3f}'.format(float(U))) + r"     $\beta=$" + str('{0:.3f}'.format(float(beta))) + r"     $\Omega_{\rm PP}=$" + str(bos) + r"$*2\pi/\beta$")

plotre_2pgf_upup_pp(pl.subplot(2,2,1))
pl.xlabel(r"$\omega_n$", fontsize=20)
pl.ylabel(r"$\omega_m$", fontsize=20)
pl.grid()
plotim_2pgf_upup_pp(pl.subplot(2,2,2) )
pl.grid()
plotre_2pgf_updo_pp(pl.subplot(2,2,3))
pl.grid()
plotim_2pgf_updo_pp(pl.subplot(2,2,4) )
pl.grid()

pl.show()
pl.clf()


#------------------------------- PLOTTING GENERALIZED SUSCEPTIBILITY -----------------------------------

print ("Plotting generalized Chi...")

pl.rc('xtick', labelsize=9) 
pl.rc('ytick', labelsize=9) 
cmap = pl.get_cmap('jet') # Bianconiglio

pl.figsize=(13, 7)

bos = 0

def plotVert_ED( use_pl, zarr, string):
    use_pl.set_aspect(1.0)
    pl.pcolormesh(np.array([(2*i+1)*pi/beta for i in range(-ffreq_pp,ffreq_pp)]), np.array([(2*i+1)*pi/beta for i in range(-ffreq_pp,ffreq_pp)]),np.ma.masked_where( np.isnan(zarr), zarr ))
    use_pl.set_title( string , fontsize=10)    
    pl.colorbar(shrink=0.6)
    return

def plotVert_ED1( use_pl, zarr, string):
    use_pl.set_aspect(1.0)
    pl.pcolormesh(np.array([(2*i+1)*pi/beta for i in range(-ffreq,ffreq)]), np.array([(2*i+1)*pi/beta for i in range(-ffreq,ffreq)]),np.ma.masked_where( np.isnan(zarr), zarr ))
    use_pl.set_title( string , fontsize=10)    
    pl.colorbar(shrink=0.6)
    return

def plotre_chi_upup_pp( use_pl ):
    title=r"$Re\chi^{PP}_{\uparrow \uparrow}$"
    zarr = np.array([[ chi_upup_pp(bos,n,m).real  for n in range(-ffreq_pp,ffreq_pp)] for m in range(-ffreq_pp,ffreq_pp)])
    plotVert_ED( use_pl, zarr, title ) 
    return

def plotim_chi_upup_pp( use_pl ):
    title=r"$Im\chi^{PP}_{\uparrow \uparrow}$"
    zarr = np.array([[ chi_upup_pp(bos,n,m).imag  for n in range(-ffreq_pp,ffreq_pp)] for m in range(-ffreq_pp,ffreq_pp)])
    plotVert_ED( use_pl, zarr, title ) 
    return

def plotre_chi_updo_pp( use_pl ):
    title=r"$Re\chi^{PP}_{\uparrow \downarrow}$"
    zarr = np.array([[ chi_updo_pp(bos,n,m).real for n in range(-ffreq_pp,ffreq_pp)] for m in range(-ffreq_pp,ffreq_pp)])
    plotVert_ED( use_pl, zarr, title ) 
    return

def plotim_chi_updo_pp( use_pl ):
    title=r"$Im\chi^{PP}_{2\uparrow \downarrow}$"
    zarr = np.array([[ (chi_updo_pp(bos,n,m)).imag  for n in range(-ffreq_pp,ffreq_pp)] for m in range(-ffreq_pp,ffreq_pp)])
    plotVert_ED( use_pl, zarr, title ) 
    return

pl.suptitle(r"$U=$" + str('{0:.3f}'.format(float(U))) + r"     $\beta=$" + str('{0:.3f}'.format(float(beta))) + r"     $\Omega_{\rm PP}=$" + str(bos) + r"$*2\pi/\beta$")

plotre_chi_upup_pp(pl.subplot(2,2,1))
pl.xlabel(r"$\omega_n$", fontsize=20)
pl.ylabel(r"$\omega_m$", fontsize=20)
pl.grid()
plotim_chi_upup_pp(pl.subplot(2,2,2) )
pl.grid()
plotre_chi_updo_pp(pl.subplot(2,2,3))
pl.grid()
plotim_chi_updo_pp(pl.subplot(2,2,4) )
pl.grid()

pl.show()
pl.clf()


#------------------------------- PLOTTING FULL VERTEX -----------------------------------

print ("Plotting full Vetex...")

pl.rc('xtick', labelsize=9) 
pl.rc('ytick', labelsize=9) 
cmap = pl.get_cmap('jet') # Bianconiglio

pl.figsize=(13, 7)

bos = 0

def plotVert_ED( use_pl, zarr, string):
    use_pl.set_aspect(1.0)
    pl.pcolormesh(np.array([(2*i+1)*pi/beta for i in range(-ffreq_pp,ffreq_pp)]), np.array([(2*i+1)*pi/beta for i in range(-ffreq_pp,ffreq_pp)]),np.ma.masked_where( np.isnan(zarr), zarr ))
    use_pl.set_title( string , fontsize=10)    
    pl.colorbar(shrink=0.6)
    return

def plotVert_ED1( use_pl, zarr, string):
    use_pl.set_aspect(1.0)
    pl.pcolormesh(np.array([(2*i+1)*pi/beta for i in range(-ffreq,ffreq)]), np.array([(2*i+1)*pi/beta for i in range(-ffreq,ffreq)]),np.ma.masked_where( np.isnan(zarr), zarr ))
    use_pl.set_title( string , fontsize=10)    
    pl.colorbar(shrink=0.6)
    return

def plotre_f_upup_pp( use_pl ):
    title=r"$ReF^{PP}_{\uparrow \uparrow}$"
    zarr = np.array([[ f_upup_pp(bos,n,m).real  for n in range(-ffreq_pp,ffreq_pp)] for m in range(-ffreq_pp,ffreq_pp)])
    plotVert_ED( use_pl, zarr, title ) 
    return

def plotim_f_upup_pp( use_pl ):
    title=r"$ImF^{PP}_{\uparrow \uparrow}$"
    zarr = np.array([[ f_upup_pp(bos,n,m).imag  for n in range(-ffreq_pp,ffreq_pp)] for m in range(-ffreq_pp,ffreq_pp)])
    plotVert_ED( use_pl, zarr, title ) 
    return

def plotre_f_updo_pp( use_pl ):
    title=r"$ReF^{PP}_{\uparrow \downarrow}$"
    zarr = np.array([[ f_updo_pp(bos,n,m).real for n in range(-ffreq_pp,ffreq_pp)] for m in range(-ffreq_pp,ffreq_pp)])
    plotVert_ED( use_pl, zarr, title ) 
    return

def plotim_f_updo_pp( use_pl ):
    title=r"$ImF^{PP}_{2\uparrow \downarrow}$"
    zarr = np.array([[ (f_updo_pp(bos,n,m)).imag  for n in range(-ffreq_pp,ffreq_pp)] for m in range(-ffreq_pp,ffreq_pp)])
    plotVert_ED( use_pl, zarr, title ) 
    return

pl.suptitle(r"$U=$" + str('{0:.3f}'.format(float(U))) + r"     $\beta=$" + str('{0:.3f}'.format(float(beta))) + r"     $\Omega_{\rm PP}=$" + str(bos) + r"$*2\pi/\beta$")

plotre_f_upup_pp(pl.subplot(2,2,1))
pl.xlabel(r"$\omega_n$", fontsize=20)
pl.ylabel(r"$\omega_m$", fontsize=20)
pl.grid()
plotim_f_upup_pp(pl.subplot(2,2,2) )
pl.grid()
plotre_f_updo_pp(pl.subplot(2,2,3))
pl.grid()
plotim_f_updo_pp(pl.subplot(2,2,4) )
pl.grid()

pl.show()
pl.clf()


#------------------------------- PLOTTING GAMMA  -----------------------------------

print ("Plotting Gamma...")

pl.rc('xtick', labelsize=9) 
pl.rc('ytick', labelsize=9) 
cmap = pl.get_cmap('jet') # Bianconiglio

pl.figsize=(13, 7)

bos = 0

def plotVert_ED( use_pl, zarr, string):
    use_pl.set_aspect(1.0)
    pl.pcolormesh(np.array([(2*i+1)*pi/beta for i in range(-ffreq_pp,ffreq_pp)]), np.array([(2*i+1)*pi/beta for i in range(-ffreq_pp,ffreq_pp)]),np.ma.masked_where( np.isnan(zarr), zarr ))
    use_pl.set_title( string , fontsize=10)    
    pl.colorbar(shrink=0.6)
    return

def plotVert_ED1( use_pl, zarr, string):
    use_pl.set_aspect(1.0)
    pl.pcolormesh(np.array([(2*i+1)*pi/beta for i in range(-ffreq,ffreq)]), np.array([(2*i+1)*pi/beta for i in range(-ffreq,ffreq)]),np.ma.masked_where( np.isnan(zarr), zarr ))
    use_pl.set_title( string , fontsize=10)    
    pl.colorbar(shrink=0.6)
    return

def plotre_gamma_upup_pp( use_pl ):
    title=r"$ReGamma^{PP}_{\uparrow \uparrow}$"
    zarr = gamma_upup_pp_arr[bfreq_pp+bos,:,:].real
#    zarr = np.array([[ gamma_upup_pp(bos,n,m).real  for n in range(-ffreq_pp,ffreq_pp)] for m in range(-ffreq_pp,ffreq_pp)])
    plotVert_ED( use_pl, zarr, title ) 
    return

def plotim_gamma_upup_pp( use_pl ):
    title=r"$ImGamma^{PP}_{\uparrow \uparrow}$"
    zarr = gamma_upup_pp_arr[bfreq_pp+bos,:,:].imag
    plotVert_ED( use_pl, zarr, title ) 
    return

def plotre_gamma_updo_pp( use_pl ):
    title=r"$ReGamma^{PP}_{\uparrow \downarrow}$"
    zarr = gamma_updo_pp_arr[bfreq_pp+bos,:,:].real
    plotVert_ED( use_pl, zarr, title ) 
    return

def plotim_gamma_updo_pp( use_pl ):
    title=r"$ImGamma^{PP}_{2\uparrow \downarrow}$"
    zarr = gamma_updo_pp_arr[bfreq_pp+bos,:,:].imag
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


#------------------------------- PLOTTING PHI PP  -----------------------------------

print ("Plotting Phi PP...")

pl.rc('xtick', labelsize=9) 
pl.rc('ytick', labelsize=9) 
cmap = pl.get_cmap('jet') # Bianconiglio

pl.figsize=(13, 7)

bos = 0

def plotVert_ED( use_pl, zarr, string):
    use_pl.set_aspect(1.0)
    pl.pcolormesh(np.array([(2*i+1)*pi/beta for i in range(-ffreq_pp,ffreq_pp)]), np.array([(2*i+1)*pi/beta for i in range(-ffreq_pp,ffreq_pp)]),np.ma.masked_where( np.isnan(zarr), zarr ))
    use_pl.set_title( string , fontsize=10)    
    pl.colorbar(shrink=0.6)
    return

def plotVert_ED1( use_pl, zarr, string):
    use_pl.set_aspect(1.0)
    pl.pcolormesh(np.array([(2*i+1)*pi/beta for i in range(-ffreq,ffreq)]), np.array([(2*i+1)*pi/beta for i in range(-ffreq,ffreq)]),np.ma.masked_where( np.isnan(zarr), zarr ))
    use_pl.set_title( string , fontsize=10)    
    pl.colorbar(shrink=0.6)
    return

def plotre_phi_upup_pp( use_pl ):
    title=r"$RePhi^{PP}_{\uparrow \uparrow}$"
    zarr = phi_upup_pp_arr[bfreq_pp+bos,:,:].real
#    zarr = np.array([[ phi_upup_pp(bos,n,m).real  for n in range(-ffreq_pp,ffreq_pp)] for m in range(-ffreq_pp,ffreq_pp)])
    plotVert_ED( use_pl, zarr, title ) 
    return

def plotim_phi_upup_pp( use_pl ):
    title=r"$ImPhi^{PP}_{\uparrow \uparrow}$"
    zarr = phi_upup_pp_arr[bfreq_pp+bos,:,:].imag
    plotVert_ED( use_pl, zarr, title ) 
    return

def plotre_phi_updo_pp( use_pl ):
    title=r"$RePhi^{PP}_{\uparrow \downarrow}$"
    zarr = phi_updo_pp_arr[bfreq_pp+bos,:,:].real
    plotVert_ED( use_pl, zarr, title ) 
    return

def plotim_phi_updo_pp( use_pl ):
    title=r"$ImPhi^{PP}_{2\uparrow \downarrow}$"
    zarr = phi_updo_pp_arr[bfreq_pp+bos,:,:].imag
    plotVert_ED( use_pl, zarr, title ) 
    return

pl.suptitle(r"$U=$" + str('{0:.3f}'.format(float(U))) + r"     $\beta=$" + str('{0:.3f}'.format(float(beta))) + r"     $\Omega_{\rm PP}=$" + str(bos) + r"$*2\pi/\beta$")

plotre_phi_upup_pp(pl.subplot(2,2,1))
pl.xlabel(r"$\omega_n$", fontsize=20)
pl.ylabel(r"$\omega_m$", fontsize=20)
pl.grid()
plotim_phi_upup_pp(pl.subplot(2,2,2) )
pl.grid()
plotre_phi_updo_pp(pl.subplot(2,2,3))
pl.grid()
plotim_phi_updo_pp(pl.subplot(2,2,4) )
pl.grid()

pl.show()
pl.clf()


#--------------------------------------- 2PGF PH------------------------------------------------------------

vertex_ph = np.loadtxt("../pomerol_output/U_1.0_beta_10.0_2bs/2pgf_ph_shift.dat")
ffreq_original_ph = int(np.transpose(vertex_ph)[1,:].max()+1)
print ffreq_original_ph
bfreq_original_ph = int(np.transpose(vertex_ph)[0,:].max())
print bfreq_original_ph

Nl =3
#wNum = int(ffreq_original_ph/100.0*10.0)
wNum = 1
iMin = ffreq_original_ph-wNum

inv_fit = generate_invfit_func(iMin,wNum,Nl)
#g_arr = np.array([G(i) for i in range(-N_fermi_gf,N_fermi_gf)])

def re_2pgf_upup_ph(wb,wf,wf1):
	return 	vertex_ph[(2*ffreq_original_ph)*(2*ffreq_original_ph)*(wb+bfreq_original_ph)+(2*ffreq_original_ph)*(wf+ffreq_original_ph) + (wf1+ffreq_original_ph), 3]

def re_2pgf_updo_ph(wb,wf,wf1):
	return 	vertex_ph[(2*ffreq_original_ph)*(2*ffreq_original_ph)*(wb+bfreq_original_ph)+(2*ffreq_original_ph)*(wf+ffreq_original_ph)+ (wf1+ffreq_original_ph), 5]

def im_2pgf_upup_ph(wb,wf,wf1):
	return 	vertex_ph[(2*ffreq_original_ph)*(2*ffreq_original_ph)*(wb+bfreq_original_ph)+(2*ffreq_original_ph)*(wf+ffreq_original_ph) + (wf1+ffreq_original_ph), 4]

def im_2pgf_updo_ph(wb,wf,wf1):
	return 	vertex_ph[(2*ffreq_original_ph)*(2*ffreq_original_ph)*(wb+bfreq_original_ph)+(2*ffreq_original_ph)*(wf+ffreq_original_ph)+ (wf1+ffreq_original_ph), 6]

#------------------------------------- GENERALIZED SUSCEPTIBILITY CHI PH --------------------------------------------

#We need to throw away wNum frequencies from the original 2pgf (because of the gamma fitting) -> we store all quantities in this new range 

ffreq_ph = ffreq_original_ph - wNum 
bfreq_ph = bfreq_original_ph

# Cutting one disconnected diagram -> generalized susceptibility

def chi_upup_ph(wb,wf,wf1):
    if (wf == wf1):
        return re_2pgf_upup_ph(wb,wf,wf1)+1j*im_2pgf_upup_ph(wb,wf,wf1)+beta*(G(wf-myfloor_div2(wb))*G(wf1+myceil_div2(wb)))
    else:
        return re_2pgf_upup_ph(wb,wf,wf1)+1j*im_2pgf_upup_ph(wb,wf,wf1)

def chi_updo_ph(wb,wf,wf1):
    return re_2pgf_updo_ph(wb,wf,wf1)+1j*im_2pgf_updo_ph(wb,wf,wf1)

def tpgf_d_ph(wb,wf,wf1):
    return re_2pgf_upup_ph(wb,wf,wf1)+1j*im_2pgf_upup_ph(wb,wf,wf1)+re_2pgf_updo_ph(wb,wf,wf1)+1j*im_2pgf_updo_ph(wb,wf,wf1)

def tpgf_m_ph(wb,wf,wf1):
    return re_2pgf_upup_ph(wb,wf,wf1)+1j*im_2pgf_upup_ph(wb,wf,wf1) -re_2pgf_updo_ph(wb,wf,wf1)+1j*im_2pgf_updo_ph(wb,wf,wf1)

#----------------------------------- FULL VERTEX PH -----------------------------------------------------------

def chi_l_ph(wb,wf,wf1):
    return G(wf-myfloor_div2(wb))*G(wf1+myceil_div2(wb))

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

def f_upup_ph(wb,wf,wf1):
    return (1.0/chi_l_ph(wb,wf,wf1))*(chi_upup_ph(wb,wf,wf1)-chi_0_ph(wb,wf,wf1))*(1.0/chi_l_ph(wb,wf1,wf))

def f_updo_ph(wb,wf,wf1):
    return (1.0/chi_l_ph(wb,wf,wf1))*(chi_updo_ph(wb,wf,wf1)-chi_0_ph(wb,wf,wf1))*(1.0/chi_l_ph(wb,wf1,wf))

#---------------------------------- 2PI VERTEX GAMMA PH ----------------------------------------------------------

def m_d_arr_ph(wb,wf,wf1):
    return tpgf_d_ph(wb,wf,wf1)-2*chi_0_ph(wb,wf,wf1)

def m_m_arr_ph(wb,wf,wf1):
    return tpgf_m_ph(wb,wf,wf1)

gamma_d_arr = np.array([beta*beta*(inv_fit(partial(chi_x0_ph,wb))-inv_fit(partial(m_d_arr_ph,wb))) for wb in range (-bfreq_ph,bfreq_ph+1)])
gamma_m_arr = np.array([beta*beta*(inv_fit(partial(chi_x0_ph,wb))-inv_fit(partial(m_m_arr_ph,wb))) for wb in range (-bfreq_ph,bfreq_ph+1)])

gamma_upup_ph_arr = 0.5*(gamma_d_arr + gamma_m_arr)
gamma_updo_ph_arr = 0.5*(gamma_d_arr - gamma_m_arr)

print gamma_upup_ph_arr.shape
print gamma_upup_ph_arr[0+bfreq_pp,:,:].max() ,gamma_updo_ph_arr[0+bfreq_pp,0+ffreq_pp,2*ffreq_pp-1]
#----------------------------------2PR VERTEX PHI----------------------------------------------------------------------------------------

f_upup_ph_arr = np.array([[[f_upup_ph(wb,wf,wf1) for wf in range (-ffreq_ph,ffreq_ph)] for wf1 in range(-ffreq_ph,ffreq_ph)] for wb in range (-bfreq_ph,bfreq_ph+1)])
f_updo_ph_arr = np.array([[[f_updo_ph(wb,wf,wf1) for wf in range (-ffreq_ph,ffreq_ph)] for wf1 in range(-ffreq_ph,ffreq_ph)] for wb in range (-bfreq_ph,bfreq_ph+1)])

phi_upup_ph_arr = f_upup_ph_arr - gamma_upup_ph_arr
phi_updo_ph_arr = f_updo_ph_arr - gamma_updo_ph_arr

#------------------------------- PLOTTING 2PGF -----------------------------------

print ("Plotting generalized 2PGF...")

pl.rc('xtick', labelsize=9) 
pl.rc('ytick', labelsize=9) 
cmap = pl.get_cmap('jet') # Bianconiglio

pl.figsize=(13, 7)

bos = 0

def plotVert_ED( use_pl, zarr, string):
    use_pl.set_aspect(1.0)
    pl.pcolormesh(np.array([(2*i+1)*pi/beta for i in range(-ffreq_ph,ffreq_ph)]), np.array([(2*i+1)*pi/beta for i in range(-ffreq_ph,ffreq_ph)]),np.ma.masked_where( np.isnan(zarr), zarr ))
#    ax.set_yscale('log')
#    ax.set_xscale('log')
    use_pl.set_title( string , fontsize=10)    
    pl.colorbar(shrink=0.6)
    return

def plotVert_ED1( use_pl, zarr, string):
    use_pl.set_aspect(1.0)
    pl.pcolormesh(np.array([(2*i+1)*pi/beta for i in range(-ffreq,ffreq)]), np.array([(2*i+1)*pi/beta for i in range(-ffreq,ffreq)]),np.ma.masked_where( np.isnan(zarr), zarr ))
    use_pl.set_title( string , fontsize=10)    
    pl.colorbar(shrink=0.6)
    return

def plotre_2pgf_upup_ph( use_pl ):
    title=r"$ReG^{PH}_{2\uparrow \uparrow}$"
    zarr = np.array([[ re_2pgf_upup_ph(bos,n,m)  for n in range(-ffreq_ph,ffreq_ph)] for m in range(-ffreq_ph,ffreq_ph)])
    plotVert_ED( use_pl, zarr, title ) 
    return

def plotim_2pgf_upup_ph( use_pl ):
    title=r"$ImG^{PH}_{2\uparrow \uparrow}$"
    zarr = np.array([[im_2pgf_upup_ph(bos,n,m)  for n in range(-ffreq_ph,ffreq_ph)] for m in range(-ffreq_ph,ffreq_ph)])
    plotVert_ED( use_pl, zarr, title ) 
    return

def plotre_2pgf_updo_ph( use_pl ):
    title=r"$ReG^{PH}_{2\uparrow \downarrow}$"
    zarr = np.array([[ re_2pgf_updo_ph(bos,n,m)  for n in range(-ffreq_ph,ffreq_ph)] for m in range(-ffreq_ph,ffreq_ph)])
    plotVert_ED( use_pl, zarr, title ) 
    return

def plotim_2pgf_updo_ph( use_pl ):
    title=r"$ImG^{PH}_{2\uparrow \downarrow}$"
    zarr = np.array([[im_2pgf_updo_ph(bos,n,m)  for n in range(-ffreq_ph,ffreq_ph)] for m in range(-ffreq_ph,ffreq_ph)])
    plotVert_ED( use_pl, zarr, title ) 
    return

pl.suptitle(r"$U=$" + str('{0:.3f}'.format(float(U))) + r"     $\beta=$" + str('{0:.3f}'.format(float(beta))) + r"     $\Omega_{\rm PH}=$" + str(bos) + r"$*2\pi/\beta$")

plotre_2pgf_upup_ph(pl.subplot(2,2,1))
pl.xlabel(r"$\omega_n$", fontsize=20)
pl.ylabel(r"$\omega_m$", fontsize=20)
pl.grid()
plotim_2pgf_upup_ph(pl.subplot(2,2,2) )
pl.grid()
plotre_2pgf_updo_ph(pl.subplot(2,2,3))
pl.grid()
plotim_2pgf_updo_ph(pl.subplot(2,2,4) )
pl.grid()

pl.show()
pl.clf()


#------------------------------- PLOTTING GENERALIZED SUSCEPTIBILITY -----------------------------------

print ("Plotting generalized Chi...")

pl.rc('xtick', labelsize=9) 
pl.rc('ytick', labelsize=9) 
cmap = pl.get_cmap('jet') # Bianconiglio

pl.figsize=(13, 7)

bos = 0

def plotVert_ED( use_pl, zarr, string):
    use_pl.set_aspect(1.0)
    pl.pcolormesh(np.array([(2*i+1)*pi/beta for i in range(-ffreq_ph,ffreq_ph)]), np.array([(2*i+1)*pi/beta for i in range(-ffreq_ph,ffreq_ph)]),np.ma.masked_where( np.isnan(zarr), zarr ))
    use_pl.set_title( string , fontsize=10)    
    pl.colorbar(shrink=0.6)
    return

def plotVert_ED1( use_pl, zarr, string):
    use_pl.set_aspect(1.0)
    pl.pcolormesh(np.array([(2*i+1)*pi/beta for i in range(-ffreq,ffreq)]), np.array([(2*i+1)*pi/beta for i in range(-ffreq,ffreq)]),np.ma.masked_where( np.isnan(zarr), zarr ))
    use_pl.set_title( string , fontsize=10)    
    pl.colorbar(shrink=0.6)
    return

def plotre_chi_upup_ph( use_pl ):
    title=r"$Re\chi^{PH}_{\uparrow \uparrow}$"
    zarr = np.array([[ chi_upup_ph(bos,n,m).real  for n in range(-ffreq_ph,ffreq_ph)] for m in range(-ffreq_ph,ffreq_ph)])
    plotVert_ED( use_pl, zarr, title ) 
    return

def plotim_chi_upup_ph( use_pl ):
    title=r"$Im\chi^{PH}_{\uparrow \uparrow}$"
    zarr = np.array([[ chi_upup_ph(bos,n,m).imag  for n in range(-ffreq_ph,ffreq_ph)] for m in range(-ffreq_ph,ffreq_ph)])
    plotVert_ED( use_pl, zarr, title ) 
    return

def plotre_chi_updo_ph( use_pl ):
    title=r"$Re\chi^{PH}_{\uparrow \downarrow}$"
    zarr = np.array([[ chi_updo_ph(bos,n,m).real for n in range(-ffreq_ph,ffreq_ph)] for m in range(-ffreq_ph,ffreq_ph)])
    plotVert_ED( use_pl, zarr, title ) 
    return

def plotim_chi_updo_ph( use_pl ):
    title=r"$Im\chi^{PH}_{2\uparrow \downarrow}$"
    zarr = np.array([[ (chi_updo_ph(bos,n,m)).imag  for n in range(-ffreq_ph,ffreq_ph)] for m in range(-ffreq_ph,ffreq_ph)])
    plotVert_ED( use_pl, zarr, title ) 
    return

pl.suptitle(r"$U=$" + str('{0:.3f}'.format(float(U))) + r"     $\beta=$" + str('{0:.3f}'.format(float(beta))) + r"     $\Omega_{\rm PH}=$" + str(bos) + r"$*2\pi/\beta$")

plotre_chi_upup_ph(pl.subplot(2,2,1))
pl.xlabel(r"$\omega_n$", fontsize=20)
pl.ylabel(r"$\omega_m$", fontsize=20)
pl.grid()
plotim_chi_upup_ph(pl.subplot(2,2,2) )
pl.grid()
plotre_chi_updo_ph(pl.subplot(2,2,3))
pl.grid()
plotim_chi_updo_ph(pl.subplot(2,2,4) )
pl.grid()

pl.show()
pl.clf()


#------------------------------- PLOTTING FULL VERTEX -----------------------------------

print ("Plotting full Vetex...")

pl.rc('xtick', labelsize=9) 
pl.rc('ytick', labelsize=9) 
cmap = pl.get_cmap('jet') # Bianconiglio

pl.figsize=(13, 7)

bos = 0

def plotVert_ED( use_pl, zarr, string):
    use_pl.set_aspect(1.0)
    pl.pcolormesh(np.array([(2*i+1)*pi/beta for i in range(-ffreq_ph,ffreq_ph)]), np.array([(2*i+1)*pi/beta for i in range(-ffreq_ph,ffreq_ph)]),np.ma.masked_where( np.isnan(zarr), zarr ))
    use_pl.set_title( string , fontsize=10)    
    pl.colorbar(shrink=0.6)
    return

def plotVert_ED1( use_pl, zarr, string):
    use_pl.set_aspect(1.0)
    pl.pcolormesh(np.array([(2*i+1)*pi/beta for i in range(-ffreq,ffreq)]), np.array([(2*i+1)*pi/beta for i in range(-ffreq,ffreq)]),np.ma.masked_where( np.isnan(zarr), zarr ))
    use_pl.set_title( string , fontsize=10)    
    pl.colorbar(shrink=0.6)
    return

def plotre_f_upup_ph( use_pl ):
    title=r"$ReF^{PH}_{\uparrow \uparrow}$"
    zarr = np.array([[ f_upup_ph(bos,n,m).real  for n in range(-ffreq_ph,ffreq_ph)] for m in range(-ffreq_ph,ffreq_ph)])
    plotVert_ED( use_pl, zarr, title ) 
    return

def plotim_f_upup_ph( use_pl ):
    title=r"$ImF^{PH}_{\uparrow \uparrow}$"
    zarr = np.array([[ f_upup_ph(bos,n,m).imag  for n in range(-ffreq_ph,ffreq_ph)] for m in range(-ffreq_ph,ffreq_ph)])
    plotVert_ED( use_pl, zarr, title ) 
    return

def plotre_f_updo_ph( use_pl ):
    title=r"$ReF^{PH}_{\uparrow \downarrow}$"
    zarr = np.array([[ f_updo_ph(bos,n,m).real for n in range(-ffreq_ph,ffreq_ph)] for m in range(-ffreq_ph,ffreq_ph)])
    plotVert_ED( use_pl, zarr, title ) 
    return

def plotim_f_updo_ph( use_pl ):
    title=r"$ImF^{PH}_{2\uparrow \downarrow}$"
    zarr = np.array([[ (f_updo_ph(bos,n,m)).imag  for n in range(-ffreq_ph,ffreq_ph)] for m in range(-ffreq_ph,ffreq_ph)])
    plotVert_ED( use_pl, zarr, title ) 
    return

pl.suptitle(r"$U=$" + str('{0:.3f}'.format(float(U))) + r"     $\beta=$" + str('{0:.3f}'.format(float(beta))) + r"     $\Omega_{\rm PH}=$" + str(bos) + r"$*2\pi/\beta$")

plotre_f_upup_ph(pl.subplot(2,2,1))
pl.xlabel(r"$\omega_n$", fontsize=20)
pl.ylabel(r"$\omega_m$", fontsize=20)
pl.grid()
plotim_f_upup_ph(pl.subplot(2,2,2) )
pl.grid()
plotre_f_updo_ph(pl.subplot(2,2,3))
pl.grid()
plotim_f_updo_ph(pl.subplot(2,2,4) )
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
    pl.pcolormesh(np.array([(2*i+1)*pi/beta for i in range(-ffreq_pp,ffreq_pp)]), np.array([(2*i+1)*pi/beta for i in range(-ffreq_pp,ffreq_pp)]),np.ma.masked_where( np.isnan(zarr), zarr ))
    use_pl.set_title( string , fontsize=10)    
    pl.colorbar(shrink=0.6)
    return

def plotVert_ED1( use_pl, zarr, string):
    use_pl.set_aspect(1.0)
    pl.pcolormesh(np.array([(2*i+1)*pi/beta for i in range(-ffreq,ffreq)]), np.array([(2*i+1)*pi/beta for i in range(-ffreq,ffreq)]),np.ma.masked_where( np.isnan(zarr), zarr ))
    use_pl.set_title( string , fontsize=10)    
    pl.colorbar(shrink=0.6)
    return

def plotre_gamma_upup_ph( use_pl ):
    title=r"$ReGamma^{PH}_{\uparrow \uparrow}$"
    zarr = gamma_upup_ph_arr[bfreq_ph+bos,:,:].real
#    zarr = np.array([[ gamma_upup_ph(bos,n,m).real  for n in range(-ffreq_ph,ffreq_ph)] for m in range(-ffreq_ph,ffreq_ph)])
    plotVert_ED( use_pl, zarr, title ) 
    return

def plotim_gamma_upup_ph( use_pl ):
    title=r"$ImGamma^{PH}_{\uparrow \uparrow}$"
    zarr = gamma_upup_ph_arr[bfreq_ph+bos,:,:].imag
    plotVert_ED( use_pl, zarr, title ) 
    return

def plotre_gamma_updo_ph( use_pl ):
    title=r"$ReGamma^{PH}_{\uparrow \downarrow}$"
    zarr = gamma_updo_ph_arr[bfreq_ph+bos,:,:].real
    plotVert_ED( use_pl, zarr, title ) 
    return

def plotim_gamma_updo_ph( use_pl ):
    title=r"$ImGamma^{PH}_{2\uparrow \downarrow}$"
    zarr = gamma_updo_ph_arr[bfreq_ph+bos,:,:].imag
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


#------------------------------- PLOTTING PHI PH  -----------------------------------

print ("Plotting Phi PH...")

pl.rc('xtick', labelsize=9) 
pl.rc('ytick', labelsize=9) 
cmap = pl.get_cmap('jet') # Bianconiglio

pl.figsize=(13, 7)

bos = 0

def plotVert_ED( use_pl, zarr, string):
    use_pl.set_aspect(1.0)
    pl.pcolormesh(np.array([(2*i+1)*pi/beta for i in range(-ffreq_ph,ffreq_ph)]), np.array([(2*i+1)*pi/beta for i in range(-ffreq_ph,ffreq_ph)]),np.ma.masked_where( np.isnan(zarr), zarr ))
    use_pl.set_title( string , fontsize=10)    
    pl.colorbar(shrink=0.6)
    return

def plotVert_ED1( use_pl, zarr, string):
    use_pl.set_aspect(1.0)
    pl.pcolormesh(np.array([(2*i+1)*pi/beta for i in range(-ffreq,ffreq)]), np.array([(2*i+1)*pi/beta for i in range(-ffreq,ffreq)]),np.ma.masked_where( np.isnan(zarr), zarr ))
    use_pl.set_title( string , fontsize=10)    
    pl.colorbar(shrink=0.6)
    return

def plotre_phi_upup_ph( use_pl ):
    title=r"$RePhi^{PH}_{\uparrow \uparrow}$"
    zarr = phi_upup_ph_arr[bfreq_ph+bos,:,:].real
#    zarr = np.array([[ phi_upup_ph(bos,n,m).real  for n in range(-ffreq_ph,ffreq_ph)] for m in range(-ffreq_ph,ffreq_ph)])
    plotVert_ED( use_pl, zarr, title ) 
    return

def plotim_phi_upup_ph( use_pl ):
    title=r"$ImPhi^{PH}_{\uparrow \uparrow}$"
    zarr = phi_upup_ph_arr[bfreq_ph+bos,:,:].imag
    plotVert_ED( use_pl, zarr, title ) 
    return

def plotre_phi_updo_ph( use_pl ):
    title=r"$RePhi^{PH}_{\uparrow \downarrow}$"
    zarr = phi_updo_ph_arr[bfreq_ph+bos,:,:].real
    plotVert_ED( use_pl, zarr, title ) 
    return

def plotim_phi_updo_ph( use_pl ):
    title=r"$ImPhi^{PH}_{2\uparrow \downarrow}$"
    zarr = phi_updo_ph_arr[bfreq_ph+bos,:,:].imag
    plotVert_ED( use_pl, zarr, title ) 
    return


pl.suptitle(r"$U=$" + str('{0:.3f}'.format(float(U))) + r"     $\beta=$" + str('{0:.3f}'.format(float(beta))) + r"     $\Omega_{\rm PH}=$" + str(bos) + r"$*2\pi/\beta$")

plotre_phi_upup_ph(pl.subplot(2,2,1))
pl.xlabel(r"$\omega_n$", fontsize=20)
pl.ylabel(r"$\omega_m$", fontsize=20)
pl.grid()
plotim_phi_upup_ph(pl.subplot(2,2,2) )
pl.grid()
plotre_phi_updo_ph(pl.subplot(2,2,3))
pl.grid()
plotim_phi_updo_ph(pl.subplot(2,2,4) )
pl.grid()

pl.show()
pl.clf()


#--------------------------------------- 2PGF XPH------------------------------------------------------------

vertex_xph = np.loadtxt("../pomerol_output/U_1.0_beta_10.0_2bs/2pgf_xph_shift.dat")
ffreq_original_xph = int(np.transpose(vertex_xph)[1,:].max()+1)
print ffreq_original_xph
bfreq_original_xph = int(np.transpose(vertex_xph)[0,:].max())
print bfreq_original_xph

Nl =3
#wNum = int(ffreq_original_pp/100.0*10.0)
wNum = 1
iMin = ffreq_original_xph-wNum
#g_arr = np.array([G(i) for i in range(-N_fermi_gf,N_fermi_gf)])
inv_fit = generate_invfit_func(iMin,wNum,Nl)


def re_2pgf_upup_xph(wb,wf,wf1):
	return 	vertex_xph[(2*ffreq_original_xph)*(2*ffreq_original_xph)*(wb+bfreq_original_xph)+(2*ffreq_original_xph)*(wf+ffreq_original_xph) + (wf1+ffreq_original_xph), 3]

def re_2pgf_updo_xph(wb,wf,wf1):
	return 	vertex_xph[(2*ffreq_original_xph)*(2*ffreq_original_xph)*(wb+bfreq_original_xph)+(2*ffreq_original_xph)*(wf+ffreq_original_xph)+ (wf1+ffreq_original_xph), 5]

def im_2pgf_upup_xph(wb,wf,wf1):
	return 	vertex_xph[(2*ffreq_original_xph)*(2*ffreq_original_xph)*(wb+bfreq_original_xph)+(2*ffreq_original_xph)*(wf+ffreq_original_xph) + (wf1+ffreq_original_xph), 4]

def im_2pgf_updo_xph(wb,wf,wf1):
	return 	vertex_xph[(2*ffreq_original_xph)*(2*ffreq_original_xph)*(wb+bfreq_original_xph)+(2*ffreq_original_xph)*(wf+ffreq_original_xph)+ (wf1+ffreq_original_xph), 6]

#------------------------------------- GENERALIZED SUSCEPTIBILITY CHI XPH --------------------------------------------

#We need to throw away wNum frequencies from the original 2pgf (because of the gamma fitting) -> we store all quantities in this new range 

ffreq_xph = ffreq_original_xph - wNum 
bfreq_xph = bfreq_original_xph

# Cutting one disconnected diagram -> generalized susceptibility

def chi_upup_xph(wb,wf,wf1):
    if (wb == 0):
        return re_2pgf_upup_xph(wb,wf,wf1)+1j*im_2pgf_upup_xph(wb,wf,wf1)+beta*G(wf)*G(wf1)
    else:
        return re_2pgf_upup_xph(wb,wf,wf1)+1j*im_2pgf_upup_xph(wb,wf,wf1)

def chi_updo_xph(wb,wf,wf1):
    return re_2pgf_updo_xph(wb,wf,wf1)+1j*im_2pgf_updo_xph(wb,wf,wf1)

#----------------------------------- FULL VERTEX XPH -----------------------------------------------------------

def chi_l_xph(wb,wf,wf1):
    return G(wf-myfloor_div2(wb))*G(wf1+myceil_div2(wb))

def chi_0_xph(wb,wf,wf1):
    if (wf == wf1):
        return beta*G(wf-myfloor_div2(wb))*G(wf+myceil_div2(wb))
    else:
        return 0.0

def chi_x0_xph(wb,wf,wf1):
    if (wb == 0):
        return -beta*G(wf)*G(wf1)
    else:
        return 0.0

def f_upup_xph(wb,wf,wf1):
    return (1.0/chi_l_xph(wb,wf,wf1))*(chi_upup_xph(wb,wf,wf1)-chi_0_xph(wb,wf,wf1))*(1.0/chi_l_xph(wb,wf1,wf))

def f_updo_xph(wb,wf,wf1):
    return (1.0/chi_l_xph(wb,wf,wf1))*(chi_updo_xph(wb,wf,wf1)-chi_0_xph(wb,wf,wf1))*(1.0/chi_l_xph(wb,wf1,wf))


#---------------------------------- 2PI VERTEX GAMMA XPH ----------------------------------------------------------

def m_updo_arr_xph(wb,wf,wf1):
    return re_2pgf_updo_xph(wb,wf,wf1)+1j*im_2pgf_updo_xph(wb,wf,wf1)

gamma_updo_xph_arr = np.array([beta*beta*(+inv_fit(partial(chi_0_xph,wb))-inv_fit(partial(m_updo_arr_xph,wb))) for wb in range (-bfreq_xph,bfreq_xph+1)])

gamma_upup_xph_arr = -gamma_upup_ph_arr

print gamma_upup_xph_arr.shape
print gamma_updo_xph_arr[0+bfreq_xph,0+ffreq_xph,2*ffreq_xph-1]
#----------------------------------2PR VERTEX XPHI----------------------------------------------------------------------------------------

f_upup_xph_arr = np.array([[[f_upup_xph(wb,wf,wf1) for wf in range (-ffreq_xph,ffreq_xph)] for wf1 in range(-ffreq_xph,ffreq_xph)] for wb in range (-bfreq_xph,bfreq_xph+1)])
f_updo_xph_arr = np.array([[[f_updo_xph(wb,wf,wf1) for wf in range (-ffreq_xph,ffreq_xph)] for wf1 in range(-ffreq_xph,ffreq_xph)] for wb in range (-bfreq_xph,bfreq_xph+1)])

phi_upup_xph_arr = f_upup_xph_arr - gamma_upup_xph_arr
phi_updo_xph_arr = f_updo_xph_arr - gamma_updo_xph_arr

#------------------------------- PLOTTING 2PGF -----------------------------------

print ("Plotting generalized 2PGF...")

pl.rc('xtick', labelsize=9) 
pl.rc('ytick', labelsize=9) 
cmap = pl.get_cmap('jet') # Bianconiglio

pl.figsize=(13, 7)

bos = 0

def plotVert_ED( use_pl, zarr, string):
    use_pl.set_aspect(1.0)
    pl.pcolormesh(np.array([(2*i+1)*pi/beta for i in range(-ffreq_xph,ffreq_xph)]), np.array([(2*i+1)*pi/beta for i in range(-ffreq_xph,ffreq_xph)]),np.ma.masked_where( np.isnan(zarr), zarr ))
    use_pl.set_title( string , fontsize=10)    
    pl.colorbar(shrink=0.6)
    return

def plotVert_ED1( use_pl, zarr, string):
    use_pl.set_aspect(1.0)
    pl.pcolormesh(np.array([(2*i+1)*pi/beta for i in range(-ffreq,ffreq)]), np.array([(2*i+1)*pi/beta for i in range(-ffreq,ffreq)]),np.ma.masked_where( np.isnan(zarr), zarr ))
    use_pl.set_title( string , fontsize=10)    
    pl.colorbar(shrink=0.6)
    return

def plotre_2pgf_upup_xph( use_pl ):
    title=r"$ReG^{XPH}_{2\uparrow \uparrow}$"
    zarr = np.array([[ re_2pgf_upup_xph(bos,n,m)  for n in range(-ffreq_xph,ffreq_xph)] for m in range(-ffreq_xph,ffreq_xph)])
    plotVert_ED( use_pl, zarr, title ) 
    return

def plotim_2pgf_upup_xph( use_pl ):
    title=r"$ImG^{XPH}_{2\uparrow \uparrow}$"
    zarr = np.array([[ im_2pgf_upup_xph(bos,n,m)  for n in range(-ffreq_xph,ffreq_xph)] for m in range(-ffreq_xph,ffreq_xph)])
    plotVert_ED( use_pl, zarr, title ) 
    return

def plotre_2pgf_updo_xph( use_pl ):
    title=r"$ReG^{XPH}_{2\uparrow \downarrow}$"
    zarr = np.array([[ re_2pgf_updo_xph(bos,n,m)  for n in range(-ffreq_xph,ffreq_xph)] for m in range(-ffreq_xph,ffreq_xph)])
    plotVert_ED( use_pl, zarr, title ) 
    return

def plotim_2pgf_updo_xph( use_pl ):
    title=r"$ImG^{XPH}_{2\uparrow \downarrow}$"
    zarr = np.array([[ im_2pgf_updo_xph(bos,n,m)  for n in range(-ffreq_xph,ffreq_xph)] for m in range(-ffreq_xph,ffreq_xph)])
    plotVert_ED( use_pl, zarr, title ) 
    return

pl.suptitle(r"$U=$" + str('{0:.3f}'.format(float(U))) + r"     $\beta=$" + str('{0:.3f}'.format(float(beta))) + r"     $\Omega_{\rm XPH}=$" + str(bos) + r"$*2\pi/\beta$")

plotre_2pgf_upup_xph(pl.subplot(2,2,1))
pl.xlabel(r"$\omega_n$", fontsize=20)
pl.ylabel(r"$\omega_m$", fontsize=20)
pl.grid()
plotim_2pgf_upup_xph(pl.subplot(2,2,2) )
pl.grid()
plotre_2pgf_updo_xph(pl.subplot(2,2,3))
pl.grid()
plotim_2pgf_updo_xph(pl.subplot(2,2,4) )
pl.grid()

pl.show()
pl.clf()


#------------------------------- PLOTTING GENERALIZED SUSCEPTIBILITY -----------------------------------

print ("Plotting generalized Chi...")

pl.rc('xtick', labelsize=9) 
pl.rc('ytick', labelsize=9) 
cmap = pl.get_cmap('jet') # Bianconiglio

pl.figsize=(13, 7)

bos = 0

def plotVert_ED( use_pl, zarr, string):
    use_pl.set_aspect(1.0)
    pl.pcolormesh(np.array([(2*i+1)*pi/beta for i in range(-ffreq_xph,ffreq_xph)]), np.array([(2*i+1)*pi/beta for i in range(-ffreq_xph,ffreq_xph)]),np.ma.masked_where( np.isnan(zarr), zarr ))
    use_pl.set_title( string , fontsize=10)    
    pl.colorbar(shrink=0.6)
    return

def plotVert_ED1( use_pl, zarr, string):
    use_pl.set_aspect(1.0)
    pl.pcolormesh(np.array([(2*i+1)*pi/beta for i in range(-ffreq,ffreq)]), np.array([(2*i+1)*pi/beta for i in range(-ffreq,ffreq)]),np.ma.masked_where( np.isnan(zarr), zarr ))
    use_pl.set_title( string , fontsize=10)    
    pl.colorbar(shrink=0.6)
    return

def plotre_chi_upup_xph( use_pl ):
    title=r"$Re\chi^{XPH}_{\uparrow \uparrow}$"
    zarr = np.array([[ chi_upup_xph(bos,n,m).real  for n in range(-ffreq_xph,ffreq_xph)] for m in range(-ffreq_xph,ffreq_xph)])
    plotVert_ED( use_pl, zarr, title ) 
    return

def plotim_chi_upup_xph( use_pl ):
    title=r"$Im\chi^{XPH}_{\uparrow \uparrow}$"
    zarr = np.array([[ chi_upup_xph(bos,n,m).imag  for n in range(-ffreq_xph,ffreq_xph)] for m in range(-ffreq_xph,ffreq_xph)])
    plotVert_ED( use_pl, zarr, title ) 
    return

def plotre_chi_updo_xph( use_pl ):
    title=r"$Re\chi^{XPH}_{\uparrow \downarrow}$"
    zarr = np.array([[ chi_updo_xph(bos,n,m).real for n in range(-ffreq_xph,ffreq_xph)] for m in range(-ffreq_xph,ffreq_xph)])
    plotVert_ED( use_pl, zarr, title ) 
    return

def plotim_chi_updo_xph( use_pl ):
    title=r"$Im\chi^{XPH}_{2\uparrow \downarrow}$"
    zarr = np.array([[ (chi_updo_xph(bos,n,m)).imag  for n in range(-ffreq_xph,ffreq_xph)] for m in range(-ffreq_xph,ffreq_xph)])
    plotVert_ED( use_pl, zarr, title ) 
    return

pl.suptitle(r"$U=$" + str('{0:.3f}'.format(float(U))) + r"     $\beta=$" + str('{0:.3f}'.format(float(beta))) + r"     $\Omega_{\rm XPH}=$" + str(bos) + r"$*2\pi/\beta$")

plotre_chi_upup_xph(pl.subplot(2,2,1))
pl.xlabel(r"$\omega_n$", fontsize=20)
pl.ylabel(r"$\omega_m$", fontsize=20)
pl.grid()
plotim_chi_upup_xph(pl.subplot(2,2,2) )
pl.grid()
plotre_chi_updo_xph(pl.subplot(2,2,3))
pl.grid()
plotim_chi_updo_xph(pl.subplot(2,2,4) )
pl.grid()

pl.show()
pl.clf()


#------------------------------- PLOTTING FULL VERTEX -----------------------------------

print ("Plotting full Vetex...")

pl.rc('xtick', labelsize=9) 
pl.rc('ytick', labelsize=9) 
cmap = pl.get_cmap('jet') # Bianconiglio

pl.figsize=(13, 7)

bos = 0

def plotVert_ED( use_pl, zarr, string):
    use_pl.set_aspect(1.0)
    pl.pcolormesh(np.array([(2*i+1)*pi/beta for i in range(-ffreq_xph,ffreq_xph)]), np.array([(2*i+1)*pi/beta for i in range(-ffreq_xph,ffreq_xph)]),np.ma.masked_where( np.isnan(zarr), zarr ))
    use_pl.set_title( string , fontsize=10)    
    pl.colorbar(shrink=0.6)
    return

def plotVert_ED1( use_pl, zarr, string):
    use_pl.set_aspect(1.0)
    pl.pcolormesh(np.array([(2*i+1)*pi/beta for i in range(-ffreq,ffreq)]), np.array([(2*i+1)*pi/beta for i in range(-ffreq,ffreq)]),np.ma.masked_where( np.isnan(zarr), zarr ))
    use_pl.set_title( string , fontsize=10)    
    pl.colorbar(shrink=0.6)
    return

def plotre_f_upup_xph( use_pl ):
    title=r"$ReF^{XPH}_{\uparrow \uparrow}$"
    zarr = np.array([[ f_upup_xph(bos,n,m).real  for n in range(-ffreq_xph,ffreq_xph)] for m in range(-ffreq_xph,ffreq_xph)])
    plotVert_ED( use_pl, zarr, title ) 
    return

def plotim_f_upup_xph( use_pl ):
    title=r"$ImF^{XPH}_{\uparrow \uparrow}$"
    zarr = np.array([[ f_upup_xph(bos,n,m).imag  for n in range(-ffreq_xph,ffreq_xph)] for m in range(-ffreq_xph,ffreq_xph)])
    plotVert_ED( use_pl, zarr, title ) 
    return

def plotre_f_updo_xph( use_pl ):
    title=r"$ReF^{XPH}_{\uparrow \downarrow}$"
    zarr = np.array([[ f_updo_xph(bos,n,m).real for n in range(-ffreq_xph,ffreq_xph)] for m in range(-ffreq_xph,ffreq_xph)])
    plotVert_ED( use_pl, zarr, title ) 
    return

def plotim_f_updo_xph( use_pl ):
    title=r"$ImF^{XPH}_{2\uparrow \downarrow}$"
    zarr = np.array([[ (f_updo_xph(bos,n,m)).imag  for n in range(-ffreq_xph,ffreq_xph)] for m in range(-ffreq_xph,ffreq_xph)])
    plotVert_ED( use_pl, zarr, title ) 
    return

pl.suptitle(r"$U=$" + str('{0:.3f}'.format(float(U))) + r"     $\beta=$" + str('{0:.3f}'.format(float(beta))) + r"     $\Omega_{\rm XPH}=$" + str(bos) + r"$*2\pi/\beta$")

plotre_f_upup_xph(pl.subplot(2,2,1))
pl.xlabel(r"$\omega_n$", fontsize=20)
pl.ylabel(r"$\omega_m$", fontsize=20)
pl.grid()
plotim_f_upup_xph(pl.subplot(2,2,2) )
pl.grid()
plotre_f_updo_xph(pl.subplot(2,2,3))
pl.grid()
plotim_f_updo_xph(pl.subplot(2,2,4) )
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
    pl.pcolormesh(np.array([(2*i+1)*pi/beta for i in range(-ffreq_pp,ffreq_pp)]), np.array([(2*i+1)*pi/beta for i in range(-ffreq_pp,ffreq_pp)]),np.ma.masked_where( np.isnan(zarr), zarr ))
    use_pl.set_title( string , fontsize=10)    
    pl.colorbar(shrink=0.6)
    return

def plotVert_ED1( use_pl, zarr, string):
    use_pl.set_aspect(1.0)
    pl.pcolormesh(np.array([(2*i+1)*pi/beta for i in range(-ffreq,ffreq)]), np.array([(2*i+1)*pi/beta for i in range(-ffreq,ffreq)]),np.ma.masked_where( np.isnan(zarr), zarr ))
    use_pl.set_title( string , fontsize=10)    
    pl.colorbar(shrink=0.6)
    return

def plotre_gamma_upup_xph( use_pl ):
    title=r"$ReGamma^{XPH}_{\uparrow \uparrow}$"
    zarr = gamma_upup_xph_arr[bfreq_xph+bos,:,:].real
#    zarr = np.array([[ gamma_upup_xph(bos,n,m).real  for n in range(-ffreq_xph,ffreq_xph)] for m in range(-ffreq_xph,ffreq_xph)])
    plotVert_ED( use_pl, zarr, title ) 
    return

def plotim_gamma_upup_xph( use_pl ):
    title=r"$ImGamma^{XPH}_{\uparrow \uparrow}$"
    zarr = gamma_upup_xph_arr[bfreq_xph+bos,:,:].imag
    plotVert_ED( use_pl, zarr, title ) 
    return

def plotre_gamma_updo_xph( use_pl ):
    title=r"$ReGamma^{XPH}_{\uparrow \downarrow}$"
    zarr = gamma_updo_xph_arr[bfreq_xph+bos,:,:].real
    plotVert_ED( use_pl, zarr, title ) 
    return

def plotim_gamma_updo_xph( use_pl ):
    title=r"$ImGamma^{XPH}_{2\uparrow \downarrow}$"
    zarr = gamma_updo_xph_arr[bfreq_xph+bos,:,:].imag
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


#------------------------------- PLOTTING XPHI XPH  -----------------------------------

print ("Plotting Phi XPH...")

pl.rc('xtick', labelsize=9) 
pl.rc('ytick', labelsize=9) 
cmap = pl.get_cmap('jet') # Bianconiglio

pl.figsize=(13, 7)

bos = 0

def plotVert_ED( use_pl, zarr, string):
    use_pl.set_aspect(1.0)
    pl.pcolormesh(np.array([(2*i+1)*pi/beta for i in range(-ffreq_xph,ffreq_xph)]), np.array([(2*i+1)*pi/beta for i in range(-ffreq_xph,ffreq_xph)]),np.ma.masked_where( np.isnan(zarr), zarr ))
    use_pl.set_title( string , fontsize=10)    
    pl.colorbar(shrink=0.6)
    return

def plotVert_ED1( use_pl, zarr, string):
    use_pl.set_aspect(1.0)
    pl.pcolormesh(np.array([(2*i+1)*pi/beta for i in range(-ffreq,ffreq)]), np.array([(2*i+1)*pi/beta for i in range(-ffreq,ffreq)]),np.ma.masked_where( np.isnan(zarr), zarr ))
    use_pl.set_title( string , fontsize=10)    
    pl.colorbar(shrink=0.6)
    return

def plotre_xphi_upup_xph( use_pl ):
    title=r"$RePhi^{XPH}_{\uparrow \uparrow}$"
    zarr = phi_upup_xph_arr[bfreq_xph+bos,:,:].real
#    zarr = np.array([[ phi_upup_xph(bos,n,m).real  for n in range(-ffreq_xph,ffreq_xph)] for m in range(-ffreq_xph,ffreq_xph)])
    plotVert_ED( use_pl, zarr, title ) 
    return

def plotim_xphi_upup_xph( use_pl ):
    title=r"$ImPhi^{XPH}_{\uparrow \uparrow}$"
    zarr = phi_upup_xph_arr[bfreq_xph+bos,:,:].imag
    plotVert_ED( use_pl, zarr, title ) 
    return

def plotre_xphi_updo_xph( use_pl ):
    title=r"$RePhi^{XPH}_{\uparrow \downarrow}$"
    zarr = phi_updo_xph_arr[bfreq_xph+bos,:,:].real
    plotVert_ED( use_pl, zarr, title ) 
    return

def plotim_xphi_updo_xph( use_pl ):
    title=r"$ImPhi^{XPH}_{2\uparrow \downarrow}$"
    zarr = phi_updo_xph_arr[bfreq_xph+bos,:,:].imag
    plotVert_ED( use_pl, zarr, title ) 
    return

pl.suptitle(r"$U=$" + str('{0:.3f}'.format(float(U))) + r"     $\beta=$" + str('{0:.3f}'.format(float(beta))) + r"     $\Omega_{\rm XPH}=$" + str(bos) + r"$*2\pi/\beta$")

plotre_xphi_upup_xph(pl.subplot(2,2,1))
pl.xlabel(r"$\omega_n$", fontsize=20)
pl.ylabel(r"$\omega_m$", fontsize=20)
pl.grid()
plotim_xphi_upup_xph(pl.subplot(2,2,2) )
pl.grid()
plotre_xphi_updo_xph(pl.subplot(2,2,3))
pl.grid()
plotim_xphi_updo_xph(pl.subplot(2,2,4) )
pl.grid()

pl.show()
pl.clf()

