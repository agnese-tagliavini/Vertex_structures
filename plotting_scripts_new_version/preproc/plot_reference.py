#!/usr/bin/python

#===================================================================================================================
#
#                   NOTE: Plotting script for the outputs of Vienna's TU ED program
#                   
#                   CONTENTS:
#                   - Full veretx in PP PH and XPH channels
#                   - Generalized susceptibility in all channels
#
#===================================================================================================================
#--------------------------------------IMPORTS ------------------------------------------

import h5py
import matplotlib.pyplot as pl
import numpy as np
import os
import sys
import subprocess
import math
#from agneselibrary.mymath import *
#from agneselibrary.plot_tools import *
#from agneselibrary.translate_notation import *
from matplotlib.colors import LinearSegmentedColormap

#--------------------------------------MANAGING FILES ------------------------------------------

def run(command):
    output = subprocess.check_output(command, shell=True)
    return output

#most_recently_edited = run("ls -Art dat/ | tail -n 1")
#
#fname = "/home/gprimes/Coding/Vertex_structures/dat/U1.0_beta26.0_FFREQ_159_BFREQ_160.h5"
#
#if len(sys.argv) > 1:
#    fname = str(sys.argv[1])
#
#fname = fname.rstrip('\n') # strip newline of fname
#f = h5py.File(fname, "r")

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

shift=0

#---------------------------------------------------------------------------------------------------------
#
#                               PP
#
#---------------------------------------------------------------------------------------------------------
#--------------------------------------- VERTX PP------------------------------------------------------------

vertex_updo_pp = np.loadtxt('/home/gprimes/Coding/Vertex_structures/dat/vertex/vert_chi_pp_full')
print vertex_updo_pp.shape
ffreq_pp = int(0.5*(np.transpose(vertex_updo_pp)[1,:].max()*beta/pi-1))+1
print ffreq_pp
bfreq_pp = int(0.5*np.transpose(vertex_updo_pp)[0,:].max()*beta/pi)
print bfreq_pp


def re_f_upup_pp(wb,wf,wf1):
	return 	vertex_updo_pp[(2*ffreq_pp)*(2*ffreq_pp)*(wb+bfreq_pp)+(2*ffreq_pp)*(wf+ffreq_pp) + (wf1+ffreq_pp), 3]

def re_f_updo_pp(wb,wf,wf1):
	return 	vertex_updo_pp[(2*ffreq_pp)*(2*ffreq_pp)*(wb+bfreq_pp)+(2*ffreq_pp)*(wf+ffreq_pp)+ (wf1+ffreq_pp), 5]

def im_f_upup_pp(wb,wf,wf1):
	return 	vertex_updo_pp[(2*ffreq_pp)*(2*ffreq_pp)*(wb+bfreq_pp)+(2*ffreq_pp)*(wf+ffreq_pp) + (wf1+ffreq_pp), 4]

def im_f_updo_pp(wb,wf,wf1):
	return 	vertex_updo_pp[(2*ffreq_pp)*(2*ffreq_pp)*(wb+bfreq_pp)+(2*ffreq_pp)*(wf+ffreq_pp)+ (wf1+ffreq_pp), 6]

#---------------------------------------------------------------------------------------------------------
#
#                               PH
#
#---------------------------------------------------------------------------------------------------------
#--------------------------------------- VERTX PH------------------------------------------------------------

vertex_updo = np.loadtxt('/home/gprimes/Coding/Vertex_structures/dat/vertex/vert_chi_full')
print vertex_updo.shape
ffreq_pp = int(0.5*(np.transpose(vertex_updo)[1,:].max()*beta/pi-1))+1
print ffreq_pp
bfreq_pp = int(0.5*np.transpose(vertex_updo)[0,:].max()*beta/pi)
print bfreq_pp


def re_f_upup_ph(wb,wf,wf1):
	return 	vertex_updo[(2*ffreq_pp)*(2*ffreq_pp)*(wb+bfreq_pp)+(2*ffreq_pp)*(wf+ffreq_pp) + (wf1+ffreq_pp), 3]

def re_f_updo_ph(wb,wf,wf1):
	return 	vertex_updo[(2*ffreq_pp)*(2*ffreq_pp)*(wb+bfreq_pp)+(2*ffreq_pp)*(wf+ffreq_pp)+ (wf1+ffreq_pp), 5]

def im_f_upup_ph(wb,wf,wf1):

	return 	vertex_updo[(2*ffreq_pp)*(2*ffreq_pp)*(wb+bfreq_pp)+(2*ffreq_pp)*(wf+ffreq_pp) + (wf1+ffreq_pp), 4]
def im_f_updo_ph(wb,wf,wf1):
	return 	vertex_updo[(2*ffreq_pp)*(2*ffreq_pp)*(wb+bfreq_pp)+(2*ffreq_pp)*(wf+ffreq_pp)+ (wf1+ffreq_pp), 6]

#---------------------------------------------------------------------------------------------------------
#
#                               XPH
#
#---------------------------------------------------------------------------------------------------------
#--------------------------------------- VERTX XPH------------------------------------------------------------

def re_f_upup_xph(wb,wf,wf1):
	return 	-vertex_updo[(2*ffreq_pp)*(2*ffreq_pp)*(wb+bfreq_pp)+(2*ffreq_pp)*(wf+ffreq_pp) + (wf1+ffreq_pp), 3]

def re_f_updo_xph(wb,wf,wf1):
	return 	-vertex_updo[(2*ffreq_pp)*(2*ffreq_pp)*(wb+bfreq_pp)+(2*ffreq_pp)*(wf+ffreq_pp)+ (wf1+ffreq_pp), 3]+vertex_updo[(2*ffreq_pp)*(2*ffreq_pp)*(wb+bfreq_pp)+(2*ffreq_pp)*(wf+ffreq_pp) + (wf1+ffreq_pp), 5]

def im_f_upup_xph(wb,wf,wf1):
	return 	-vertex_updo[(2*ffreq_pp)*(2*ffreq_pp)*(wb+bfreq_pp)+(2*ffreq_pp)*(wf+ffreq_pp) + (wf1+ffreq_pp), 4]

def im_f_updo_xph(wb,wf,wf1):
	return 	-vertex_updo[(2*ffreq_pp)*(2*ffreq_pp)*(wb+bfreq_pp)+(2*ffreq_pp)*(wf+ffreq_pp)+ (wf1+ffreq_pp), 4]+vertex_updo[(2*ffreq_pp)*(2*ffreq_pp)*(wb+bfreq_pp)+(2*ffreq_pp)*(wf+ffreq_pp) + (wf1+ffreq_pp), 6]

#-------------------------------Plotting Vertex-----------------


N_fermi_plot = 20

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
    zarr = np.array([[ re_f_upup_ph(shift,n,m)  for n in range(-N_fermi_plot,N_fermi_plot)] for m in range(-N_fermi_plot,N_fermi_plot)])
    plotVert_ED( use_pl, zarr, title ) 
    return

def plotre_f_updo_ph( use_pl ):
    title=r"$\operatorname{Re}F^{PH}_{\uparrow \downarrow}$"
    zarr = np.array([[ re_f_updo_ph(shift,n,m)  for n in range(-N_fermi_plot,N_fermi_plot)] for m in range(-N_fermi_plot,N_fermi_plot)])
    plotVert_ED( use_pl, zarr, title ) 
    return

def plotre_f_upup_pp( use_pl ):
    title=r"$\operatorname{Re}F^{PP}_{\uparrow \uparrow}$"
    zarr = np.array([[ re_f_upup_pp(shift,n,m) for n in range(-N_fermi_plot,N_fermi_plot)] for m in range(-N_fermi_plot,N_fermi_plot)])
    plotVert_ED( use_pl, zarr, title ) 
    return

def plotre_f_updo_pp( use_pl ):
    title=r"$\operatorname{Re}F^{PP}_{\uparrow \downarrow}$"
    zarr = np.array([[ re_f_updo_pp(shift,n,m) for n in range(-N_fermi_plot,N_fermi_plot)] for m in range(-N_fermi_plot,N_fermi_plot)])
    plotVert_ED( use_pl, zarr, title ) 
    return

def plotre_f_upup_xph( use_pl ):
    title=r"$\operatorname{Re}F^{XPH}_{\uparrow \uparrow}$"
    zarr = np.array([[ (re_f_upup_xph(shift,n,m))  for n in range(-N_fermi_plot,N_fermi_plot)] for m in range(-N_fermi_plot,N_fermi_plot)])
    plotVert_ED( use_pl, zarr, title ) 
    return

def plotre_f_updo_xph( use_pl ):
    title=r"$\operatorname{Re}F^{XPH}_{\uparrow \downarrow}$"
    zarr = np.array([[ (re_f_updo_xph(shift,n,m))  for n in range(-N_fermi_plot,N_fermi_plot)] for m in range(-N_fermi_plot,N_fermi_plot)])
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

#---------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------
#
#                               PP
#
#---------------------------------------------------------------------------------------------------------
#--------------------------------------- GENCHI PP------------------------------------------------------------

genchi_updo_pp = np.loadtxt('/home/gprimes/Coding/Vertex_structures/dat/vertex/vert_chi_pp')
print genchi_updo_pp.shape
ffreq_pp = int(0.5*(np.transpose(genchi_updo_pp)[1,:].max()*beta/pi-1))+1
print ffreq_pp
bfreq_pp = int(0.5*np.transpose(genchi_updo_pp)[0,:].max()*beta/pi)
print bfreq_pp


def re_genchi_upup_pp(wb,wf,wf1):
	return 	genchi_updo_pp[(2*ffreq_pp)*(2*ffreq_pp)*(wb+bfreq_pp)+(2*ffreq_pp)*(wf+ffreq_pp) + (wf1+ffreq_pp), 3]

def re_genchi_updo_pp(wb,wf,wf1):
	return 	genchi_updo_pp[(2*ffreq_pp)*(2*ffreq_pp)*(wb+bfreq_pp)+(2*ffreq_pp)*(wf+ffreq_pp)+ (wf1+ffreq_pp), 5]

def im_genchi_upup_pp(wb,wf,wf1):
	return 	genchi_updo_pp[(2*ffreq_pp)*(2*ffreq_pp)*(wb+bfreq_pp)+(2*ffreq_pp)*(wf+ffreq_pp) + (wf1+ffreq_pp), 4]

def im_genchi_updo_pp(wb,wf,wf1):
	return 	genchi_updo_pp[(2*ffreq_pp)*(2*ffreq_pp)*(wb+bfreq_pp)+(2*ffreq_pp)*(wf+ffreq_pp)+ (wf1+ffreq_pp), 6]

#---------------------------------------------------------------------------------------------------------
#
#                               PH
#
#---------------------------------------------------------------------------------------------------------
#--------------------------------------- GENCHI PH------------------------------------------------------------

genchi_updo = np.loadtxt('/home/gprimes/Coding/Vertex_structures/dat/vertex/vert_chi')
print genchi_updo.shape
ffreq_pp = int(0.5*(np.transpose(genchi_updo)[1,:].max()*beta/pi-1))+1
print ffreq_pp
bfreq_pp = int(0.5*np.transpose(genchi_updo)[0,:].max()*beta/pi)
print bfreq_pp


def re_genchi_upup_ph(wb,wf,wf1):
	return 	genchi_updo[(2*ffreq_pp)*(2*ffreq_pp)*(wb+bfreq_pp)+(2*ffreq_pp)*(wf+ffreq_pp) + (wf1+ffreq_pp), 3]

def re_genchi_updo_ph(wb,wf,wf1):
	return 	genchi_updo[(2*ffreq_pp)*(2*ffreq_pp)*(wb+bfreq_pp)+(2*ffreq_pp)*(wf+ffreq_pp)+ (wf1+ffreq_pp), 5]

def im_genchi_upup_ph(wb,wf,wf1):

	return 	genchi_updo[(2*ffreq_pp)*(2*ffreq_pp)*(wb+bfreq_pp)+(2*ffreq_pp)*(wf+ffreq_pp) + (wf1+ffreq_pp), 4]
def im_genchi_updo_ph(wb,wf,wf1):
	return 	genchi_updo[(2*ffreq_pp)*(2*ffreq_pp)*(wb+bfreq_pp)+(2*ffreq_pp)*(wf+ffreq_pp)+ (wf1+ffreq_pp), 6]

#---------------------------------------------------------------------------------------------------------
#
#                               XPH
#
#---------------------------------------------------------------------------------------------------------
#--------------------------------------- GENCHI XPH------------------------------------------------------------

def re_genchi_upup_xph(wb,wf,wf1):
	return 	-genchi_updo[(2*ffreq_pp)*(2*ffreq_pp)*(wb+bfreq_pp)+(2*ffreq_pp)*(wf+ffreq_pp) + (wf1+ffreq_pp), 3]

def re_genchi_updo_xph(wb,wf,wf1):
	return 	-genchi_updo[(2*ffreq_pp)*(2*ffreq_pp)*(wb+bfreq_pp)+(2*ffreq_pp)*(wf+ffreq_pp)+ (wf1+ffreq_pp), 3]+genchi_updo[(2*ffreq_pp)*(2*ffreq_pp)*(wb+bfreq_pp)+(2*ffreq_pp)*(wf+ffreq_pp) + (wf1+ffreq_pp), 5]

def im_genchi_upup_xph(wb,wf,wf1):
	return 	-genchi_updo[(2*ffreq_pp)*(2*ffreq_pp)*(wb+bfreq_pp)+(2*ffreq_pp)*(wf+ffreq_pp) + (wf1+ffreq_pp), 4]

def im_genchi_updo_xph(wb,wf,wf1):
	return 	-genchi_updo[(2*ffreq_pp)*(2*ffreq_pp)*(wb+bfreq_pp)+(2*ffreq_pp)*(wf+ffreq_pp)+ (wf1+ffreq_pp), 4]+genchi_updo[(2*ffreq_pp)*(2*ffreq_pp)*(wb+bfreq_pp)+(2*ffreq_pp)*(wf+ffreq_pp) + (wf1+ffreq_pp), 6]

#-------------------------------Plotting GENCHI-----------------


N_fermi_plot = ffreq_pp

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

def plotre_chi_upup_ph( use_pl ):
    title=r"$\operatorname{Re}\chi^{PH}_{\uparrow \uparrow}$"
    zarr = np.array([[ re_genchi_upup_ph(shift,n,m)  for n in range(-N_fermi_plot,N_fermi_plot)] for m in range(-N_fermi_plot,N_fermi_plot)])
    plotVert_ED( use_pl, zarr, title ) 
    return

def plotre_chi_updo_ph( use_pl ):
    title=r"$\operatorname{Re}\chi^{PH}_{\uparrow \downarrow}$"
    zarr = np.array([[ re_genchi_updo_ph(shift,n,m)  for n in range(-N_fermi_plot,N_fermi_plot)] for m in range(-N_fermi_plot,N_fermi_plot)])
    plotVert_ED( use_pl, zarr, title ) 
    return

def plotre_chi_upup_pp( use_pl ):
    title=r"$\operatorname{Re}\chi^{PP}_{\uparrow \uparrow}$"
    zarr = np.array([[ re_genchi_upup_pp(shift,n,m) for n in range(-N_fermi_plot,N_fermi_plot)] for m in range(-N_fermi_plot,N_fermi_plot)])
    plotVert_ED( use_pl, zarr, title ) 
    return

def plotre_chi_updo_pp( use_pl ):
    title=r"$\operatorname{Re}\chi^{PP}_{\uparrow \downarrow}$"
    zarr = np.array([[ re_genchi_updo_pp(shift,n,m) for n in range(-N_fermi_plot,N_fermi_plot)] for m in range(-N_fermi_plot,N_fermi_plot)])
    plotVert_ED( use_pl, zarr, title ) 
    return

def plotre_chi_upup_xph( use_pl ):
    title=r"$\operatorname{Re}\chi^{XPH}_{\uparrow \uparrow}$"
    zarr = np.array([[ (re_genchi_upup_xph(shift,n,m))  for n in range(-N_fermi_plot,N_fermi_plot)] for m in range(-N_fermi_plot,N_fermi_plot)])
    plotVert_ED( use_pl, zarr, title ) 
    return

def plotre_chi_updo_xph( use_pl ):
    title=r"$\operatorname{Re}\chi^{XPH}_{\uparrow \downarrow}$"
    zarr = np.array([[ (re_genchi_updo_xph(shift,n,m))  for n in range(-N_fermi_plot,N_fermi_plot)] for m in range(-N_fermi_plot,N_fermi_plot)])
    plotVert_ED( use_pl, zarr, title ) 
    return

pl.suptitle(r"$U=$" + str('{0:.3f}'.format(float(U))) + r"     $\beta=$" + str('{0:.3f}'.format(float(beta))) + r"     $\Omega_{\rm PP}= \Omega_{\rm PH}= \Omega_{\rm XPH}=$" + str(shift) + r"$*2\pi/\beta$")

plotre_chi_upup_pp(pl.subplot(2,3,1))
pl.ylabel(r"$\omega_m$",fontsize=10)
#pl.xticks([-40,-20,0,20,40])
#pl.yticks([-40,-20,0,20,40])
plotre_chi_upup_ph(pl.subplot(2,3,2))
#pl.xticks([-40,-20,0,20,40])
#pl.yticks([-40,-20,0,20,40])
plotre_chi_upup_xph(pl.subplot(2,3,3))
#pl.xticks([-40,-20,0,20,40])
#pl.yticks([-40,-20,0,20,40])
plotre_chi_updo_pp(pl.subplot(2,3,4))
pl.xlabel(r"$\omega_n$",fontsize=10)
pl.ylabel(r"$\omega_m$",fontsize=10)
#pl.xticks([-40,-20,0,20,40])
#pl.yticks([-40,-20,0,20,40])
plotre_chi_updo_ph(pl.subplot(2,3,5))
#pl.xticks([-40,-20,0,20,40])
#pl.yticks([-40,-20,0,20,40])
pl.xlabel(r"$\omega_n$",fontsize=10)
plotre_chi_updo_xph(pl.subplot(2,3,6))
pl.xlabel(r"$\omega_n$",fontsize=10)
#pl.xticks([-40,-20,0,20,40])
#pl.yticks([-40,-20,0,20,40])

pl.tight_layout()

#--- Save to file
pl.savefig("plots/GENCHI.png", dpi = 150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()