#!/usr/bin/python

#=============================================================================================================
#
#                       Plottings script for Pomerol results before extracting the asymptotics
#                       and inverting the Bethe-Salpeter equations
#
#                       CONTENTS:
#                       - Self-energy
#                       - Vertex in all channels
#                       - Generalized susceptibility in all channels
#                       - 2P Green's function in all channels
#
#============================================================================================================

#--------------------------------------IMPORTS ------------------------------------------
import h5py
import matplotlib.pyplot as pl
import numpy as np
import os
import sys
import subprocess
import math
from agneselibrary.mymath import *
from agneselibrary.plot_tools import *
from agneselibrary.translate_notation import *
from matplotlib.colors import LinearSegmentedColormap

#--------------------------------------MANAGING FILES ------------------------------------------

def run(command):
    output = subprocess.check_output(command, shell=True)
    return output

most_recently_edited = run("ls -Art dat/ | tail -n 1")

fname = "/home/agnese/Coding/Vertex_structures/dat/H5FILES/BETA26/4SITES/U1/PREPROC/U1.0_BETA_26.0_FFREQ_50_BFREQ_75_noQN_idx.h5"

if len(sys.argv) > 1:
    fname = str(sys.argv[1])

fname = fname.rstrip('\n') # strip newline of fname
f = h5py.File(fname, "r")

os.system('mkdir -p log')
os.system('mkdir -p plots')

os.system('rm plots/prev_Sig.png 2> log/plot.log')
os.system('rm plots/prev_Giw.png 2> log/plot.log')
os.system('rm plots/prev_Vert.png 2> log/plot.log')
os.system('rm plots/prev_Phi.png 2> log/plot.log')
os.system('rm plots/prev_P.png 2> log/plot.log')
os.system('rm plots/prev_K.png 2> log/plot.log')
os.system('rm plots/prev_Trileg.png 2> log/plot.log')
os.system('rm plots/prev_Chi.png 2> log/plot.log')
os.system('rm plots/prev_R.png 2> log/plot.log')
os.system('mv plots/Sig.png plots/prev_Sig.png 2> log/plot.log')
os.system('mv plots/Giw.png plots/prev_Giw.png 2> log/plot.log') 
os.system('mv plots/Vert.png plots/prev_Vert.png 2> log/plot.log')
os.system('mv plots/Phi.png plots/prev_Phi.png 2> log/plot.log')
os.system('mv plots/P.png plots/prev_P.png 2> log/plot.log')
os.system('mv plots/K.png plots/prev_K.png 2> log/plot.log')
os.system('mv plots/Trileg.png plots/prev_Trileg.png 2> log/plot.log')
os.system('mv plots/Chi.png plots/prev_Chi.png 2> log/plot.log')
os.system('mv plots/R.png plots/prev_R.png 2> log/plot.log')

#--------------------------------------READ PARAMETERS FROM FILE ------------------------------------------
print "Check UINT and BETA in plot.py!"

U =  1.0          # follows order in script_conversion_hdf5_demetrio.py
beta = 26.0
pi = math.pi

shift=30

#-----------------------------------GREEN'S FUNCTION--------------------------

print("Plotting Green function ...")


#--- Read
Giwgrid = np.array(f["/Giw/fgrid"])
reGiw = f["/Giw/RE"]
imGiw = f["/Giw/IM"]
fdim = reGiw.shape[0]


#--- Helper functions

def plotGiw( use_pl, arr, string ):
    pl.plot( Giwgrid, arr, 'bx', ms=3, mew=0.2)
#    pl.xlim([1.2*min(vertgrid),1.2*max(vertgrid)])
    use_pl.set_title(string)
    return

def plotGiwRe( use_pl ):
    plotGiw( use_pl, reGiw[:], r"$\operatorname{Re}G(i\omega)$")
    return

def plotGiwIm( use_pl ):
    plotGiw( use_pl, imGiw[:], r"$\operatorname{Im}G(i\omega)$")
    return

pl.suptitle(r"$U=$" + str('{0:.3f}'.format(float(U))) + r"     $\beta=$" + str('{0:.3f}'.format(float(beta))))

plotGiwRe( pl.subplot(1,2,1) ) 
pl.xlabel(r"$\omega_n$",fontsize=8)
plotGiwIm( pl.subplot(1,2,2) ) 
pl.xlabel(r"$\omega_n$",fontsize=8)

pl.tight_layout()

#--- Save to file
pl.savefig("plots/Giw.png", dpi=150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()


#-----------------------------------SE FUNCTION--------------------------

print("Plotting SE ...")


#--- Read
Siggrid = np.array(f["/Sig/fgrid"])
reSig = f["/Sig/RE"]
imSig = f["/Sig/IM"]
fdim = reSig.shape[0]


#--- Helper functions

def plotSig( use_pl, arr, string ):
    pl.plot( Siggrid, arr, 'bx', ms=3, mew=0.2)
#    pl.xlim([1.2*min(vertgrid),1.2*max(vertgrid)])
    use_pl.set_title(string)
    return

def plotSigRe( use_pl ):
    plotGiw( use_pl, reSig[:], r"$\operatorname{Re}\Sigma(i\omega)$")
    return

def plotSigIm( use_pl ):
    plotGiw( use_pl, imSig[:], r"$\operatorname{Im}\Sigma(i\omega)$")
    return

pl.suptitle(r"$U=$" + str('{0:.3f}'.format(float(U))) + r"     $\beta=$" + str('{0:.3f}'.format(float(beta))))

plotSigRe( pl.subplot(1,2,1) ) 
pl.xlabel(r"$\omega_n$",fontsize=8)
plotSigIm( pl.subplot(1,2,2) ) 
pl.xlabel(r"$\omega_n$",fontsize=8)

pl.tight_layout()

#--- Save to file
pl.savefig("plots/Sig.png", dpi=150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()


#--------------------------------------VERTEX PLOTTING ------------------------------------------

print("Plotting vertex ...")


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
print N_bose
print N_fermi

def isInside(i,j,k):
    return abs(i) <= N_bose and j >= -N_fermi and j < N_fermi and k >= -N_fermi and k < N_fermi


#-------------------------------- FUNCTION DEFINITION--------------------------------------

# 
# Update the asymptotic structures (for the moment firt oder diagrams outside) for the VERTEX IN ALL CHANNELS

def f_upup_fun_ph(i,j,k):
    if isInside(i,j,k):
        return re_f_upup_ph[i + N_bose, j+N_fermi, k + N_fermi]+1j*im_f_upup_ph[i + N_bose, j+N_fermi, k + N_fermi]
    else:
          return 0.0 

def f_updo_fun_ph(i,j,k):
    if isInside(i,j,k):
        return re_f_updo_ph[i + N_bose, j+N_fermi, k + N_fermi]+1j*im_f_updo_ph[i + N_bose, j+N_fermi, k + N_fermi]
    else:
        return -U

def f_upup_fun_pp(i,j,k):
    if isInside(i,j,k):
        return re_f_upup_pp[i + N_bose, j+N_fermi, k + N_fermi]+1j*im_f_updo_pp[i + N_bose, j+N_fermi, k + N_fermi]
    else:
        return 0.0

def f_updo_fun_pp(i,j,k):
    if isInside(i,j,k):
        return re_f_updo_pp[i + N_bose, j+N_fermi, k + N_fermi]+1j*im_f_updo_pp[i + N_bose, j+N_fermi, k + N_fermi]
    else:
        return -U
def f_upup_fun_xph(i,j,k):
    if isInside(i,j,k):
        return re_f_upup_xph[i + N_bose, j+N_fermi, k + N_fermi]+1j*im_f_upup_xph[i + N_bose, j+N_fermi, k + N_fermi]
    else:
        return 0.0
def f_updo_fun_xph(i,j,k):
    if isInside(i,j,k):
        return re_f_updo_xph[i + N_bose, j+N_fermi, k + N_fermi]+1j*im_f_updo_xph[i + N_bose, j+N_fermi, k + N_fermi]
    else:
        return -U

#-------------------------------Plotting Vertex-----------------


N_fermi_plot = int(1.2*N_fermi)

pl.rc('xtick', labelsize=9) 
pl.rc('ytick', labelsize=9) 
cmap = pl.get_cmap('jet') 

X = np.array([(2*i+1)*pi/beta for i in range(-N_fermi_plot,N_fermi_plot)])
#pl.figsize=(13, 7)

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
    title=r"$\operatorname{Re}F^{XPH}_{\uparrow \uparrow}$"
    zarr = np.array([[ (f_upup_fun_xph(shift,n,m)).real  for n in range(-N_fermi_plot,N_fermi_plot)] for m in range(-N_fermi_plot,N_fermi_plot)])
    plotVert_ED( use_pl, zarr, title ) 
    return

def plotre_f_updo_xph( use_pl ):
    title=r"$\operatorname{Re}F^{XPH}_{\uparrow \downarrow}$"
    zarr = np.array([[ (f_updo_fun_xph(shift,n,m)).real  for n in range(-N_fermi_plot,N_fermi_plot)] for m in range(-N_fermi_plot,N_fermi_plot)])
    plotVert_ED( use_pl, zarr, title ) 
    return

pl.suptitle(r"$U=$" + str('{0:.3f}'.format(float(U))) + r"     $\beta=$" + str('{0:.3f}'.format(float(beta))) + r"     $\Omega_{\rm PP}= \Omega_{\rm PH}= \Omega_{\rm XPH}=$" + str(shift) + r"$*2\pi/\beta$")

plotre_f_upup_pp(pl.subplot(2,3,1))
pl.ylabel(r"$\omega_m$",fontsize=10)
#pl.xticks([-40,-20,0,20,40])
#pl.yticks([-40,-20,0,20,40])
plotre_f_upup_ph(pl.subplot(2,3,2))
#pl.xticks([-40,-20,0,20,40])
#pl.yticks([-40,-20,0,20,40])
plotre_f_upup_xph(pl.subplot(2,3,3))
#pl.xticks([-40,-20,0,20,40])
#pl.yticks([-40,-20,0,20,40])
plotre_f_updo_pp(pl.subplot(2,3,4))
pl.xlabel(r"$\omega_n$",fontsize=10)
pl.ylabel(r"$\omega_m$",fontsize=10)
#pl.xticks([-40,-20,0,20,40])
#pl.yticks([-40,-20,0,20,40])
plotre_f_updo_ph(pl.subplot(2,3,5))
#pl.xticks([-40,-20,0,20,40])
#pl.yticks([-40,-20,0,20,40])
pl.xlabel(r"$\omega_n$",fontsize=10)
plotre_f_updo_xph(pl.subplot(2,3,6))
pl.xlabel(r"$\omega_n$",fontsize=10)
#pl.xticks([-40,-20,0,20,40])
#pl.yticks([-40,-20,0,20,40])

pl.tight_layout()

#--- Save to file
pl.savefig("plots/Vert.png", dpi = 150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()

#-------------------------------Plotting Generalized Susceptibility-----------------

print("Plotting Generalized Susceptibility ...")

#GENCHIs

re_genchi_upup_ph = f["GENCHI/PH/RE_GENCHI_UPUP"][:]
re_genchi_updo_ph = f["GENCHI/PH/RE_GENCHI_UPDO"][:]
re_genchi_upup_xph = f["GENCHI/XPH/RE_GENCHI_UPUP"][:]
re_genchi_updo_xph = f["GENCHI/XPH/RE_GENCHI_UPDO"][:]
re_genchi_upup_pp = f["GENCHI/PP/RE_GENCHI_UPUP"][:]
re_genchi_updo_pp = f["GENCHI/PP/RE_GENCHI_UPDO"][:]


im_genchi_upup_ph = f["GENCHI/PH/IM_GENCHI_UPUP"][:]
im_genchi_updo_ph = f["GENCHI/PH/IM_GENCHI_UPDO"][:]
im_genchi_upup_xph = f["GENCHI/XPH/IM_GENCHI_UPUP"][:]
im_genchi_updo_xph = f["GENCHI/XPH/IM_GENCHI_UPDO"][:]
im_genchi_upup_pp = f["GENCHI/PP/IM_GENCHI_UPUP"][:]
im_genchi_updo_pp = f["GENCHI/PP/IM_GENCHI_UPDO"][:]

N_fermi_plot = N_fermi

pl.rc('xtick', labelsize=9) 
pl.rc('ytick', labelsize=9) 
cmap = pl.get_cmap('jet') # Bianconiglio

pl.figsize=(13, 7)

def plotGamma_ED( use_pl, zarr, string):
    use_pl.set_aspect(1.0)
    pl.pcolormesh(np.array([(2*i+1)*pi/beta for i in range(-N_fermi_plot,N_fermi_plot)]), np.array([(2*i+1)*pi/beta for i in range(-N_fermi_plot,N_fermi_plot)]),np.ma.masked_where( np.isnan(zarr), zarr ))
    use_pl.set_title( string , fontsize=10)    
    pl.colorbar(shrink=0.6)
    return

def plotre_genchi_upup_ph( use_pl ):
    title=r"$\operatorname{Re}\chi^{PH}_{\uparrow \uparrow}$"
    zarr = re_genchi_upup_ph[shift+N_bose,:,:].real
    plotGamma_ED( use_pl, zarr, title ) 
    return

def plotre_genchi_updo_ph( use_pl ):
    title=r"$\operatorname{Re}\chi^{PH}_{\uparrow \downarrow}$"
    zarr = re_genchi_updo_ph[shift+N_bose,:,:].real
    plotGamma_ED( use_pl, zarr, title ) 
    return

def plotre_genchi_upup_pp( use_pl ):
    title=r"$\operatorname{Re}\chi^{PP}_{\uparrow \uparrow}$"
    zarr = re_genchi_upup_pp[shift+N_bose,:,:].real
    plotGamma_ED( use_pl, zarr, title ) 
    return

def plotre_genchi_updo_pp( use_pl ):
    title=r"$\operatorname{Re}\chi^{PP}_{\uparrow \downarrow}$"
    zarr = re_genchi_updo_pp[shift+N_bose,:,:].real
    plotGamma_ED( use_pl, zarr, title ) 
    return

def plotre_genchi_upup_xph( use_pl ):
    title=r"$\operatorname{Re}\chi^{XPH}_{\uparrow \uparrow}$"
    zarr = re_genchi_upup_xph[shift+N_bose,:,:].real
    plotGamma_ED( use_pl, zarr, title ) 
    return

def plotre_genchi_updo_xph( use_pl ):
    title=r"$\operatorname{Re}\chi^{XPH}_{\uparrow \downarrow}$"
    zarr = re_genchi_updo_xph[shift+N_bose,:,:].real
    plotGamma_ED( use_pl, zarr, title ) 
    return


pl.suptitle(r"$U=$" + str('{0:.3f}'.format(float(U))) + r"     $\beta=$" + str('{0:.3f}'.format(float(beta))) + r"     $\Omega_{\rm PP}= \Omega_{\rm PH}= \Omega_{\rm XPH}=$" + str(shift) + r"$*2\pi/\beta$")

plotre_genchi_upup_pp(pl.subplot(2,3,1))
pl.ylabel(r"$\omega_m$", fontsize=10)
#pl.xticks([-20,-10,0,10,20])
#pl.yticks([-20,-10,0,10,20])
plotre_genchi_upup_ph(pl.subplot(2,3,2))
#pl.xticks([-20,-10,0,10,20])
#pl.yticks([-20,-10,0,10,20])
plotre_genchi_upup_xph(pl.subplot(2,3,3))
#pl.xticks([-20,-10,0,10,20])
#pl.yticks([-20,-10,0,10,20])
plotre_genchi_updo_pp(pl.subplot(2,3,4))
#pl.xticks([-20,-10,0,10,20])
#pl.yticks([-20,-10,0,10,20])
pl.xlabel(r"$\omega_n$",fontsize=10)
pl.ylabel(r"$\omega_m$", fontsize=10)
plotre_genchi_updo_ph(pl.subplot(2,3,5))
#pl.xticks([-20,-10,0,10,20])
#pl.yticks([-20,-10,0,10,20])
pl.xlabel(r"$\omega_n$",fontsize=10)
plotre_genchi_updo_xph(pl.subplot(2,3,6))
#pl.xticks([-20,-10,0,10,20])
#pl.yticks([-20,-10,0,10,20])
pl.xlabel(r"$\omega_n$",fontsize=10)

pl.tight_layout()

#--- Save to file
pl.savefig("plots/Genchi.png", dpi = 150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()

#-------------------------------------------------2PGF--------------------------------------------

print("Plotting 2PGF ...")

#RS

re_2PGF_upup_ph = f["2PGF/PH/RE_2PGF_UPUP"][:]
re_2PGF_updo_ph = f["2PGF/PH/RE_2PGF_UPDO"][:]
re_2PGF_upup_xph = f["2PGF/XPH/RE_2PGF_UPUP"][:]
re_2PGF_updo_xph = f["2PGF/XPH/RE_2PGF_UPDO"][:]
re_2PGF_upup_pp = f["2PGF/PP/RE_2PGF_UPUP"][:]
re_2PGF_updo_pp = f["2PGF/PP/RE_2PGF_UPDO"][:]


im_2PGF_upup_ph = f["2PGF/PH/IM_2PGF_UPUP"][:]
im_2PGF_updo_ph = f["2PGF/PH/IM_2PGF_UPDO"][:]
im_2PGF_upup_xph = f["2PGF/XPH/IM_2PGF_UPUP"][:]
im_2PGF_updo_xph = f["2PGF/XPH/IM_2PGF_UPDO"][:]
im_2PGF_upup_pp = f["2PGF/PP/IM_2PGF_UPUP"][:]
im_2PGF_updo_pp = f["2PGF/PP/IM_2PGF_UPDO"][:]


#-------------------------------Plotting 2PGF-----------------

N_fermi_plot = N_fermi

pl.rc('xtick', labelsize=9) 
pl.rc('ytick', labelsize=9) 
cmap = pl.get_cmap('jet') # Bianconiglio

#pl.figsize=(13, 7)

def plot2PGF_ED( use_pl, zarr, string):
    use_pl.set_aspect(1.0)
    pl.pcolormesh(np.array([(2*i+1)*pi/beta for i in range(-N_fermi_plot,N_fermi_plot)]), np.array([(2*i+1)*pi/beta for i in range(-N_fermi_plot,N_fermi_plot)]),np.ma.masked_where( np.isnan(zarr), zarr ))
    use_pl.set_title( string , fontsize=10)    
    pl.colorbar(shrink=0.6)
    return

def plotre_2PGF_upup_ph( use_pl ):
    title=r"$\operatorname{Re}G^{PH}_{2,\uparrow \uparrow}$"
    zarr =  re_2PGF_upup_ph[shift+N_bose,:,:] 
    plot2PGF_ED( use_pl, zarr, title ) 
    return

def plotre_2PGF_updo_ph( use_pl ):
    title=r"$\operatorname{Re}G^{PH}_{2,\uparrow \downarrow}$"
    zarr =  re_2PGF_updo_ph[shift+N_bose,:,:]
    plot2PGF_ED( use_pl, zarr, title ) 
    return

def plotre_2PGF_upup_pp( use_pl ):
    title=r"$\operatorname{Re}G^{PP}_{2,\uparrow \uparrow}$"
    zarr =  re_2PGF_upup_pp[shift+N_bose,:,:]
    plot2PGF_ED( use_pl, zarr, title ) 
    return

def plotre_2PGF_updo_pp( use_pl ):
    title=r"$\operatorname{Re}G^{PP}_{2,\uparrow \downarrow}}}$"
    zarr =  re_2PGF_updo_pp[shift+N_bose,:,:]    
    plot2PGF_ED( use_pl, zarr, title ) 
    return

def plotre_2PGF_upup_xph( use_pl ):
    title=r"$\operatorname{Re}G^{XPH}_{2,\uparrow \uparrow}$"
    zarr =  re_2PGF_upup_xph[shift+N_bose,:,:]
    plot2PGF_ED( use_pl, zarr, title ) 
    return

def plotre_2PGF_updo_xph( use_pl ):
    title=r"$\operatorname{Re}G^{XPH}_{2,\uparrow \downarrow}$"
    zarr =  re_2PGF_updo_xph[shift+N_bose,:,:]
    plot2PGF_ED( use_pl, zarr, title ) 
    return

pl.suptitle(r"$U=$" + str('{0:.3f}'.format(float(U))) + r"     $\beta=$" + str('{0:.3f}'.format(float(beta))) + r"     $\Omega_{\rm PP}= \Omega_{\rm PH}= \Omega_{\rm XPH}=$" + str(shift) + r"$*2\pi/\beta$")

plotre_2PGF_upup_pp(pl.subplot(2,3,1))
pl.ylabel(r"$\omega_m$",fontsize=10)
#pl.xticks([-20,-10,0,10,20])
#pl.yticks([-20,-10,0,10,20])
plotre_2PGF_upup_ph(pl.subplot(2,3,2))
#pl.xticks([-20,-10,0,10,20])
#pl.yticks([-20,-10,0,10,20])
plotre_2PGF_upup_xph(pl.subplot(2,3,3))
#pl.xticks([-20,-10,0,10,20])
#pl.yticks([-20,-10,0,10,20])
plotre_2PGF_updo_pp(pl.subplot(2,3,4))
pl.xlabel(r"$\omega_n$",fontsize=10)
pl.ylabel(r"$\omega_m$",fontsize=10)
#pl.xticks([-20,-10,0,10,20])
#pl.yticks([-20,-10,0,10,20])
plotre_2PGF_updo_ph(pl.subplot(2,3,5))
pl.xlabel(r"$\omega_n$",fontsize=10)
#pl.xticks([-20,-10,0,10,20])
#pl.yticks([-20,-10,0,10,20])
plotre_2PGF_updo_xph(pl.subplot(2,3,6))
pl.xlabel(r"$\omega_n$",fontsize=10)
#pl.xticks([-20,-10,0,10,20])
#pl.yticks([-20,-10,0,10,20])

pl.tight_layout()

#--- Save to file
pl.savefig("plots/2PGF.png", dpi = 150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()

