#!/usr/bin/python

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

fname = "dat/BETA26/U1/POMEROL/dat_U1_Beta26_PFCB150_PARQ_SU2_METH2_INVR4xVERTR.h5"


if len(sys.argv) > 1:
    fname = str(sys.argv[1])

fname = fname.rstrip('\n') # strip newline of fname
f = h5py.File(fname, "r")

os.system('mkdir -p log')
os.system('mkdir -p plots')


#--------------------------------------READ PARAMETERS FROM FILE ------------------------------------------

parVals = f["/Params"].attrs.values()

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

#---------------------------------------------------------------------------------------------------------
#
#                              GREENS FUNCTION STEFAN
#
#---------------------------------------------------------------------------------------------------------
#--------------------------------------- GF------------------------------------------------------------

gf_s = np.loadtxt('dat/BETA26/U1/STEFAN/gm_wim')

print gf_s.shape

def GF_s(w):
    if (w >= 0):
        return gf_s[w,1]+1.j*gf_s[w,2]
    else:
        return gf_s[-w-1,1]-1.j*gf_s[-w-1,2]

print gf_s.shape[0]
Giwgrid_s = np.arange(-gf_s.shape[0],gf_s.shape[0])

#--------------------------------------G(iw) PLOTTING ------------------------------------------

print("Plotting Green function ...")


#--- Read
Giwgrid_p = np.array(f["/Giw/fgrid"])
reGiw_p = f["/Giw/RE"]
imGiw_p = f["/Giw/IM"]
fdim_p = reGiw_p.shape[0]

def GF_p(w):
    return reGiw_p[w+fdim_p/2,0,0,0]+1.j*imGiw_p[w+fdim_p/2,0,0,0]

if (max(Giwgrid_p) > max(Giwgrid_s)):
    Giwgrid_plot = Giwgrid_s
else:
    Giwgrid_plot = Giwgrid_p

#--- Helper functions
def plotGiw( use_pl, arr, string ):
    pl.plot( Giwgrid_plot, arr, 'bx', ms=3, mew=0.2)
    pl.xlim([min(Giwgrid_plot),max(Giwgrid_plot)])
    use_pl.set_title(string)
    return


pl.suptitle(r"$U=$" + str(UINT) + r"     $\beta=$" + str(BETA) + r"     $\epsilon=$" + str(EPS))

#--- Plot physical
plotGiw( pl.subplot(2,2,1), np.array([GF_p(i).real for i in Giwgrid_plot]), RE + "G(i\omega)_{{pom}}$" ) 
pl.xlabel(r"$\omega_n$")
plotGiw( pl.subplot(2,2,2), np.array([GF_p(i).imag for i in Giwgrid_plot]), IM + "G(i\omega)_{{pom}}$" ) 
pl.xlabel(r"$\omega_n$")
plotGiw( pl.subplot(2,2,3), np.array([GF_s(i).real for i in Giwgrid_plot]), RE + "G(i\omega)_{{stef}}$" ) 
pl.xlabel(r"$\omega_n$")
plotGiw( pl.subplot(2,2,4), np.array([GF_s(i).imag for i in Giwgrid_plot]), IM + "G(i\omega)_{{stef}}$" ) 
pl.xlabel(r"$\omega_n$")

pl.tight_layout()

#--- Save to file
pl.savefig("plots/Giw_comp.png", dpi=150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()


shift=0



#---------------------------------------------------------------------------------------------------------
#
#                              GENCHI PP STEFAN
#
#---------------------------------------------------------------------------------------------------------
#--------------------------------------- GENCHI PP------------------------------------------------------------

genchi_pp_s = np.loadtxt('dat/BETA26/U1/STEFAN/vert_chi_pp')
print genchi_pp_s.shape
ffreq_pp = int(0.5*(np.transpose(genchi_pp_s)[1,:].max()*BETA/pi-1))+1
print ffreq_pp
bfreq_pp = int(0.5*np.transpose(genchi_pp_s)[0,:].max()*BETA/pi)
print bfreq_pp


def re_genchi_s_upup_pp(wb,wf,wf1):
	return 	genchi_pp_s[(2*ffreq_pp)*(2*ffreq_pp)*(wb+bfreq_pp)+(2*ffreq_pp)*(wf+ffreq_pp) + (wf1+ffreq_pp), 3]

def re_genchi_s_updo_pp(wb,wf,wf1):
	return 	genchi_pp_s[(2*ffreq_pp)*(2*ffreq_pp)*(wb+bfreq_pp)+(2*ffreq_pp)*(wf+ffreq_pp)+ (wf1+ffreq_pp), 5]

def im_genchi_s_upup_pp(wb,wf,wf1):
	return 	genchi_pp_s[(2*ffreq_pp)*(2*ffreq_pp)*(wb+bfreq_pp)+(2*ffreq_pp)*(wf+ffreq_pp) + (wf1+ffreq_pp), 4]

def im_genchi_s_updo_pp(wb,wf,wf1):
	return 	genchi_pp_s[(2*ffreq_pp)*(2*ffreq_pp)*(wb+bfreq_pp)+(2*ffreq_pp)*(wf+ffreq_pp)+ (wf1+ffreq_pp), 6]

#---------------------------------------------------------------------------------------------------------
#
#                               GENCHI PH STEFAN
#
#---------------------------------------------------------------------------------------------------------
#--------------------------------------- GENCHI PH------------------------------------------------------------

genchi_ph_s = np.loadtxt('dat/BETA26/U1/STEFAN/vert_chi')
print genchi_ph_s.shape
ffreq_ph = int(0.5*(np.transpose(genchi_ph_s)[1,:].max()*BETA/pi-1))+1
print ffreq_pp
bfreq_ph = int(0.5*np.transpose(genchi_ph_s)[0,:].max()*BETA/pi)
print bfreq_pp


def re_genchi_s_upup_ph(wb,wf,wf1):
	return 	genchi_ph_s[(2*ffreq_ph)*(2*ffreq_ph)*(wb+bfreq_ph)+(2*ffreq_ph)*(wf+ffreq_ph) + (wf1+ffreq_ph), 3]

def re_genchi_s_updo_ph(wb,wf,wf1):
	return 	genchi_ph_s[(2*ffreq_ph)*(2*ffreq_ph)*(wb+bfreq_ph)+(2*ffreq_ph)*(wf+ffreq_ph)+ (wf1+ffreq_ph), 5]

def im_genchi_s_upup_ph(wb,wf,wf1):

	return 	genchi_ph_s[(2*ffreq_ph)*(2*ffreq_ph)*(wb+bfreq_ph)+(2*ffreq_ph)*(wf+ffreq_ph) + (wf1+ffreq_ph), 4]
def im_genchi_s_updo_ph(wb,wf,wf1):
	return 	genchi_ph_s[(2*ffreq_ph)*(2*ffreq_ph)*(wb+bfreq_ph)+(2*ffreq_ph)*(wf+ffreq_ph)+ (wf1+ffreq_ph), 6]

#---------------------------------------------------------------------------------------------------------
#
#                               GENCHI XPH STEFAN
#
#---------------------------------------------------------------------------------------------------------
#--------------------------------------- GENCHI XPH------------------------------------------------------------

def re_genchi_s_upup_xph(wb,wf,wf1):
	return 	-genchi_ph_s[(2*ffreq_ph)*(2*ffreq_ph)*(wb+bfreq_ph)+(2*ffreq_ph)*(wf+ffreq_ph) + (wf1+ffreq_ph), 3]

def re_genchi_s_updo_xph(wb,wf,wf1):
	return 	-genchi_ph_s[(2*ffreq_ph)*(2*ffreq_ph)*(wb+bfreq_ph)+(2*ffreq_ph)*(wf+ffreq_ph)+ (wf1+ffreq_ph), 3]+genchi_ph_s[(2*ffreq_ph)*(2*ffreq_ph)*(wb+bfreq_ph)+(2*ffreq_ph)*(wf+ffreq_ph) + (wf1+ffreq_ph), 5]

def im_genchi_s_upup_xph(wb,wf,wf1):
	return 	-genchi_ph_s[(2*ffreq_ph)*(2*ffreq_ph)*(wb+bfreq_ph)+(2*ffreq_ph)*(wf+ffreq_ph) + (wf1+ffreq_ph), 4]

def im_genchi_s_updo_xph(wb,wf,wf1):
	return 	-genchi_ph_s[(2*ffreq_ph)*(2*ffreq_ph)*(wb+bfreq_ph)+(2*ffreq_ph)*(wf+ffreq_ph)+ (wf1+ffreq_ph), 4]+genchi_ph_s[(2*ffreq_ph)*(2*ffreq_ph)*(wb+bfreq_ph)+(2*ffreq_ph)*(wf+ffreq_ph) + (wf1+ffreq_ph), 6]


#--------------------------------------GENCHI PLOTTING POMEROL VS STEFAN ------------------------------------------


#--- Read
regenchi_pp = vert_mul * np.array(f["/genchi_func/RE_PP"])
imgenchi_pp = vert_mul * np.array(f["/genchi_func/IM_PP"])
regenchi_ph = vert_mul * np.array(f["/genchi_func/RE_PH"])
imgenchi_ph = vert_mul * np.array(f["/genchi_func/IM_PH"])
regenchi_xph = vert_mul * np.array(f["/genchi_func/RE_XPH"])
imgenchi_xph = vert_mul * np.array(f["/genchi_func/IM_XPH"])

bdim = regenchi_pp.shape[0]
fdim = regenchi_pp.shape[1]

genchigrid = np.array(f["/genchi_func/fgrid"])

N_fermi_plot = min([fdim/2, ffreq_pp])

print N_fermi_plot

fdim_plot = 2*N_fermi_plot


#---  Helper functions (include translation from  nambu to physical CAUTION : USING TIME REVERSAL SYMM) 
def plotgenchi_p( use_pl, arr, string ):
    use_pl.set_aspect(1.0)
    zarr = np.array([[ arr[shift + (bdim-1)/2,n+fdim/2,m+fdim/2,0,0,0,0,0,0,0] for n in range(-N_fermi_plot, N_fermi_plot)] for m in range(-N_fermi_plot, N_fermi_plot)])
    pl.pcolormesh( np.arange(-N_fermi_plot, N_fermi_plot), np.arange(-N_fermi_plot, N_fermi_plot), zarr )
    pl.ylim([-N_fermi_plot, N_fermi_plot])
    pl.xlim([-N_fermi_plot,N_fermi_plot])
    use_pl.set_title( string , fontsize=10 )
    pl.colorbar(shrink=0.6) 
    return

def plotgenchi_s( use_pl, func, string ):
    use_pl.set_aspect(1.0)
    zarr = np.array([[func(shift,n,m) for n in range(-N_fermi_plot,N_fermi_plot)] for m in range(-N_fermi_plot,N_fermi_plot)])
    pl.pcolormesh( np.arange(-N_fermi_plot, N_fermi_plot), np.arange(-N_fermi_plot, N_fermi_plot), zarr )
    pl.ylim([-N_fermi_plot,N_fermi_plot])
    pl.xlim([-N_fermi_plot,N_fermi_plot])
    use_pl.set_title( string , fontsize=10 )
    pl.colorbar(shrink=0.6) 
    return

print("Plotting genchi PP RE ...")

#--- Plot 
pl.suptitle(r"$U=$" + str(UINT) + r"     $\beta=$" + str(BETA) + r"     $\epsilon=$" + str(EPS) +  r"     $\Omega=$" + str(shift) + r"$*2\pi/\beta$" + r"     Notation: $\chi(\Omega,\omega,\omega')$")

plotgenchi_p( pl.subplot(2,2,1), regenchi_pp - regenchi_pp[:,:,::-1,:,:,:,:,:,:,:], RE + r"\chi^{PP,pom}_{\uparrow\uparrow}(\Omega_{PP},\omega,\omega')$" ) # flip sign of w_out
pl.ylabel(r"$\omega'$")
pl.xlabel(r"$\omega$")
plotgenchi_p( pl.subplot(2,2,2), regenchi_pp , RE + r"\chi^{PP,pom}_{\uparrow\downarrow}(\Omega_{PP},\omega,\omega')$" ) # flip sign of w_out
plotgenchi_s( pl.subplot(2,2,3), re_genchi_s_upup_pp, RE + r"\chi^{PP,stef}_{\uparrow\uparrow}(\Omega_{PP},\omega,\omega')$" ) # flip sign of w_out
pl.ylabel(r"$\omega'$")
pl.xlabel(r"$\omega$")
plotgenchi_s( pl.subplot(2,2,4), re_genchi_s_updo_pp , RE + r"\chi^{PP,pom}_{\uparrow\downarrow}(\Omega_{PP},\omega,\omega')$" ) # flip sign of w_out
pl.ylabel(r"$\omega'$")
pl.xlabel(r"$\omega$")

pl.tight_layout()

#--- Save to file
pl.savefig("plots/genchi_re_PP.png", dpi = 150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()

print("Plotting genchi PP IM ...")
#--- Plot 
pl.suptitle(r"$U=$" + str(UINT) + r"     $\beta=$" + str(BETA) + r"     $\epsilon=$" + str(EPS) +  r"     $\Omega=$" + str(shift) + r"$*2\pi/\beta$" + r"     Notation: $\chi(\Omega,\omega,\omega')$")

plotgenchi_p( pl.subplot(2,2,1), imgenchi_pp - imgenchi_pp[:,:,::-1,:,:,:,:,:,:,:], IM + r"\chi^{PP,pom}_{\uparrow\uparrow}(\Omega_{PP},\omega,\omega')$" ) # flip sign of w_out
pl.ylabel(r"$\omega'$")
pl.xlabel(r"$\omega$")
plotgenchi_p( pl.subplot(2,2,2), imgenchi_pp , IM + r"\chi^{PP,pom}_{\uparrow\downarrow}(\Omega_{PP},\omega,\omega')$" ) # flip sign of w_out
plotgenchi_s( pl.subplot(2,2,3), im_genchi_s_upup_pp, IM + r"\chi^{PP,stef}_{\uparrow\uparrow}(\Omega_{PP},\omega,\omega')$" ) # flip sign of w_out
pl.ylabel(r"$\omega'$")
pl.xlabel(r"$\omega$")
plotgenchi_s( pl.subplot(2,2,4), im_genchi_s_updo_pp , IM + r"\chi^{PP,pom}_{\uparrow\downarrow}(\Omega_{PP},\omega,\omega')$" ) # flip sign of w_out
pl.ylabel(r"$\omega'$")
pl.xlabel(r"$\omega$")

pl.tight_layout()

#--- Save to file
pl.savefig("plots/genchi_im_PP.png", dpi = 150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()


print("Plotting genchi PH RE ...")

#--- Plot 
pl.suptitle(r"$U=$" + str(UINT) + r"     $\beta=$" + str(BETA) + r"     $\epsilon=$" + str(EPS) +  r"     $\Omega=$" + str(shift) + r"$*2\pi/\beta$" + r"     Notation: $\chi(\Omega,\omega,\omega')$")

plotgenchi_p( pl.subplot(2,2,1), regenchi_ph - regenchi_xph, RE + r"\chi^{PH,pom}_{\uparrow\uparrow}(\Omega_{PH},\omega,\omega')$" ) # flip sign of w_out
pl.ylabel(r"$\omega'$")
pl.xlabel(r"$\omega$")
plotgenchi_p( pl.subplot(2,2,2), regenchi_ph , RE + r"\chi^{PH,pom}_{\uparrow\downarrow}(\Omega_{PH},\omega,\omega')$" ) # flip sign of w_out
plotgenchi_s( pl.subplot(2,2,3), re_genchi_s_upup_ph, RE + r"\chi^{PH,stef}_{\uparrow\uparrow}(\Omega_{PH},\omega,\omega')$" ) # flip sign of w_out
pl.ylabel(r"$\omega'$")
pl.xlabel(r"$\omega$")
plotgenchi_s( pl.subplot(2,2,4), re_genchi_s_updo_ph , RE + r"\chi^{PH,pom}_{\uparrow\downarrow}(\Omega_{PH},\omega,\omega')$" ) # flip sign of w_out
pl.ylabel(r"$\omega'$")
pl.xlabel(r"$\omega$")

pl.tight_layout()

#--- Save to file
pl.savefig("plots/genchi_re_PH.png", dpi = 150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()

print("Plotting genchi PH IM ...")
#--- Plot 
pl.suptitle(r"$U=$" + str(UINT) + r"     $\beta=$" + str(BETA) + r"     $\epsilon=$" + str(EPS) +  r"     $\Omega=$" + str(shift) + r"$*2\pi/\beta$" + r"     Notation: $\chi(\Omega,\omega,\omega')$")

plotgenchi_p( pl.subplot(2,2,1), imgenchi_ph - imgenchi_xph, IM + r"\chi^{PH,pom}_{\uparrow\uparrow}(\Omega_{PH},\omega,\omega')$" ) # flip sign of w_out
pl.ylabel(r"$\omega'$")
pl.xlabel(r"$\omega$")
plotgenchi_p( pl.subplot(2,2,2), imgenchi_ph , IM + r"\chi^{XPH,pom}_{\uparrow\downarrow}(\Omega_PH},\omega,\omega')$" ) # flip sign of w_out
plotgenchi_s( pl.subplot(2,2,3), im_genchi_s_upup_ph, IM + r"\chi^{PH,stef}_{\uparrow\uparrow}(\Omega_{PH},\omega,\omega')$" ) # flip sign of w_out
pl.ylabel(r"$\omega'$")
pl.xlabel(r"$\omega$")
plotgenchi_s( pl.subplot(2,2,4), im_genchi_s_updo_ph , IM + r"\chi^{PH,pom}_{\uparrow\downarrow}(\Omega_{PH},\omega,\omega')$" ) # flip sign of w_out
pl.ylabel(r"$\omega'$")
pl.xlabel(r"$\omega$")

pl.tight_layout()

#--- Save to file
pl.savefig("plots/genchi_im_PH.png", dpi = 150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()

print("Plotting genchi XPH RE ...")

#--- Plot 
pl.suptitle(r"$U=$" + str(UINT) + r"     $\beta=$" + str(BETA) + r"     $\epsilon=$" + str(EPS) +  r"     $\Omega=$" + str(shift) + r"$*2\pi/\beta$" + r"     Notation: $\chi(\Omega,\omega,\omega')$")

plotgenchi_p( pl.subplot(2,2,1), regenchi_xph - regenchi_ph, RE + r"\chi^{XPH,pom}_{\uparrow\uparrow}(\Omega_{XPH},\omega,\omega')$" ) # flip sign of w_out
pl.ylabel(r"$\omega'$")
pl.xlabel(r"$\omega$")
plotgenchi_p( pl.subplot(2,2,2), regenchi_xph , RE + r"\chi^{XPH,pom}_{\uparrow\downarrow}(\Omega_{XPH},\omega,\omega')$" ) # flip sign of w_out
plotgenchi_s( pl.subplot(2,2,3), re_genchi_s_upup_xph, RE + r"\chi^{XPH,stef}_{\uparrow\uparrow}(\Omega_{XPH},\omega,\omega')$" ) # flip sign of w_out
pl.ylabel(r"$\omega'$")
pl.xlabel(r"$\omega$")
plotgenchi_s( pl.subplot(2,2,4), re_genchi_s_updo_xph , RE + r"\chi^{XPH,pom}_{\uparrow\downarrow}(\Omega_{XPH},\omega,\omega')$" ) # flip sign of w_out
pl.ylabel(r"$\omega'$")
pl.xlabel(r"$\omega$")

pl.tight_layout()

#--- Save to file
pl.savefig("plots/genchi_re_XPH.png", dpi = 150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()

print("Plotting genchi XPH IM ...")
#--- Plot 
pl.suptitle(r"$U=$" + str(UINT) + r"     $\beta=$" + str(BETA) + r"     $\epsilon=$" + str(EPS) +  r"     $\Omega=$" + str(shift) + r"$*2\pi/\beta$" + r"     Notation: $\chi(\Omega,\omega,\omega')$")

plotgenchi_p( pl.subplot(2,2,1), imgenchi_xph - imgenchi_ph, IM + r"\chi^{XPH,pom}_{\uparrow\uparrow}(\Omega_{XPH},\omega,\omega')$" ) # flip sign of w_out
pl.ylabel(r"$\omega'$")
pl.xlabel(r"$\omega$")
plotgenchi_p( pl.subplot(2,2,2), imgenchi_xph , IM + r"\chi^{XPH,pom}_{\uparrow\downarrow}(\Omega_XPH},\omega,\omega')$" ) # flip sign of w_out
plotgenchi_s( pl.subplot(2,2,3), im_genchi_s_upup_xph, IM + r"\chi^{XPH,stef}_{\uparrow\uparrow}(\Omega_{XPH},\omega,\omega')$" ) # flip sign of w_out
pl.ylabel(r"$\omega'$")
pl.xlabel(r"$\omega$")
plotgenchi_s( pl.subplot(2,2,4), im_genchi_s_updo_xph , IM + r"\chi^{XPH,pom}_{\uparrow\downarrow}(\Omega_{XPH},\omega,\omega')$" ) # flip sign of w_out
pl.ylabel(r"$\omega'$")
pl.xlabel(r"$\omega$")

pl.tight_layout()

#--- Save to file
pl.savefig("plots/genchi_im_XPH.png", dpi = 150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()
