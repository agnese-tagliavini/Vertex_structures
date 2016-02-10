#!/usr/bin/python

#--------------------------------------IMPORTS ------------------------------------------

import h5py
import matplotlib.pyplot as pl
import numpy as np
import os
import sys
import subprocess
import math
from agneselib.mymath import *
from agneselib.translate_notation import *

#--------------------------------------MANAGING FILES ------------------------------------------

def run(command):
    output = subprocess.check_output(command, shell=True)
    return output

most_recently_edited = run("ls -Art dat/ | tail -n 1")

#fname = "dat/" + most_recently_edited
fname = "dat/dat_U2.0_beta10.0_EDpomerol.h5"
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

U =  2.0          # follows order in script_conversion_hdf5_demetrio.py
beta = 10.0
pi = math.pi

shift=0

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
pl.xlabel(r"$\omega_n$")
plotGiwIm( pl.subplot(1,2,2) ) 
pl.xlabel(r"$\omega_n$")

pl.tight_layout()

#--- Save to file
pl.savefig("plots/Giw.png", dpi=150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()


#--------------------------------------PLUS  PLOTTING ------------------------------------------

print("Plotting P ...")


#--- Read
rep_upup_ph =  np.array(f["/P_func/PH/RE_P_UPUP"])
imp_upup_ph = np.array(f["/P_func/PH/IM_P_UPUP"])
rep_updo_ph =  np.array(f["/P_func/PH/RE_P_UPDO"])
imp_updo_ph = np.array(f["/P_func/PH/IM_P_UPDO"])
bgrid_p = rep_upup_ph.shape[0]
fgrid_p = rep_upup_ph.shape[1]
rep_upup_pp = np.array(f["/P_func/PP/RE_P_UPUP"])
imp_upup_pp = np.array(f["/P_func/PP/IM_P_UPUP"])
rep_updo_pp =  np.array(f["/P_func/PP/RE_P_UPDO"])
imp_updo_pp = np.array(f["/P_func/PP/IM_P_UPDO"])
rep_updo_xph =  np.array(f["/P_func/XPH/RE_P_UPDO"])
imp_updo_xph = np.array(f["/P_func/XPH/IM_P_UPDO"])
rep_upup_xph =  np.array(f["/P_func/XPH/RE_P_UPUP"])
imp_upup_xph = np.array(f["/P_func/XPH/IM_P_UPUP"])

if fgrid_p <= shift:
    sys.exit("Error: Shift too large for vertex grid"); 

N_bose_p = (rep_upup_ph.shape[0]-1)/2
N_fermi_p = (rep_upup_ph.shape[1])/2
N_bose_k = (rek_upup_ph.shape[0]-1)/2 

print N_bose_p, N_fermi_p, N_bose_k

def plotP( use_pl, zarr, string ):
    pl.plot( np.arange(-N_fermi_p,N_fermi_p), zarr)
    use_pl.set_title( string , fontsize=10)
    return

def plotUpUpPRePH( use_pl ):
    title = r"$\operatorname{Re}P_{2,\uparrow\uparrow}(\Omega_{PH},\omega_n)$"
    zarr = rep_upup_ph[shift+N_bose_p,:]
    plotP( use_pl, zarr, title )
    return

def plotUpDoPRePH( use_pl ):
    title = r"$\operatorname{Re}P_{2,\uparrow\downarrow}(\Omega_{PH},\omega_n)$"
    zarr = rep_updo_ph[shift+N_bose_p,:]
    plotP( use_pl, zarr, title )
    return

def plotUpUpPRePP( use_pl ):
    title = r"$\operatorname{Re}P_{2,\uparrow\uparrow}(\Omega_{PP},\omega_n)$"
    zarr = rep_upup_pp[shift+N_bose_p,:]
    plotP( use_pl, zarr, title )
    return

def plotUpDoPRePP( use_pl ):
    title = r"$\operatorname{Re}P_{2,\uparrow\downarrow}(\Omega_{PP},\omega_n)$"
    zarr = rep_updo_pp[shift+N_bose_p,:]
    plotP( use_pl, zarr, title )
    return

def plotUpUpPReXPH( use_pl ):
    title = r"$\operatorname{Re}P_{2,\uparrow\uparrow}(\Omega_{XPH},\omega_n)$"
    zarr = rep_upup_xph[shift+N_bose_p,:]
    plotP( use_pl, zarr, title )
    return

def plotUpDoPReXPH( use_pl ):
    title = r"$\operatorname{Re}P_{2,\uparrow\downarrow}(\Omega_{XPH},\omega_n)$"
    zarr = rep_updo_xph[shift+N_bose_p,:]
    plotP( use_pl, zarr, title )
    return

pl.suptitle(r"$U=$" + str('{0:.3f}'.format(float(U))) + r"     $\beta=$" + str('{0:.3f}'.format(float(beta))) + r"     $\Omega_{\rm PH}= \Omega_{\rm PP}=$" + str(shift) + r"$*2\pi/\beta$")

#--- Plot Physical
plotUpUpPRePH( pl.subplot(2,3,2) )
pl.xlabel(r"$\omega_n$")
plotUpDoPRePH( pl.subplot(2,3,5) )
pl.xlabel(r"$\omega_m$")
plotUpUpPRePP( pl.subplot(2,3,1) )
pl.xlabel(r"$\omega_n$")
plotUpDoPRePP( pl.subplot(2,3,4) )
pl.xlabel(r"$\omega_n$")
plotUpUpPReXPH( pl.subplot(2,3,3) )
pl.xlabel(r"$\omega_n$")
plotUpDoPReXPH( pl.subplot(2,3,6) )
pl.xlabel(r"$\omega_n$")
pl.tight_layout()

#--- Save to file
pl.savefig("plots/P_func.png", dpi = 150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()


#--------------------------------------------------------------KARRASCH PLOTTING--------------------------------------------------------


print("Plotting K ...")

#----Read Karrasch

rek_upup_ph =  np.array(f["/K_func/PH/RE_K_UPUP"])
imk_upup_ph = np.array(f["/K_func/PH/IM_K_UPUP"])
rek_updo_ph =  np.array(f["/K_func/PH/RE_K_UPDO"])
imk_updo_ph = np.array(f["/K_func/PH/IM_K_UPDO"])
bgrid_k = rek_upup_ph.shape[0]
rek_upup_pp = np.array(f["/K_func/PP/RE_K_UPUP"])
imk_upup_pp = np.array(f["/K_func/PP/IM_K_UPUP"])
rek_updo_pp =  np.array(f["/K_func/PP/RE_K_UPDO"])
imk_updo_pp = np.array(f["/K_func/PP/IM_K_UPDO"])
rek_updo_xph =  np.array(f["/K_func/XPH/RE_K_UPDO"])
imk_updo_xph = np.array(f["/K_func/XPH/IM_K_UPDO"])
rek_upup_xph =  np.array(f["/K_func/XPH/RE_K_UPUP"])
imk_upup_xph = np.array(f["/K_func/XPH/IM_K_UPUP"])


N_bose_k = (rek_upup_ph.shape[0]-1)/2

print N_bose_k


def plotK( use_pl, zarr, string ):
    pl.plot( np.arange(-N_bose_k,N_bose_k+1), zarr)
    use_pl.set_title( string , fontsize=10)
    return

def plotUpUpKRePH( use_pl ):
    title = r"$\operatorname{Re}K_{2,\uparrow\uparrow}(\Omega_{PH},\omega_n)$"
    zarr = rek_upup_ph[:]
    plotK( use_pl, zarr, title )
    return

def plotUpDoKRePH( use_pl ):
    title = r"$\operatorname{Re}K_{2,\uparrow\downarrow}(\Omega_{PH},\omega_n)$"
    zarr = rek_updo_ph[:]
    plotK( use_pl, zarr, title )
    return

def plotUpUpKRePP( use_pl ):
    title = r"$\operatorname{Re}K_{2,\uparrow\uparrow}(\Omega_{PP},\omega_n)$"
    zarr = rek_updo_pp[:]
    plotK( use_pl, zarr, title )
    return

def plotUpDoKRePP( use_pl ):
    title = r"$\operatorname{Re}K_{2,\uparrow\downarrow}(\Omega_{PP},\omega_n)$"
    zarr = rek_updo_pp[:]
    plotK( use_pl, zarr, title )
    return

def plotUpUpKReXPH( use_pl ):
    title = r"$\operatorname{Re}K_{2,\uparrow\uparrow}(\Omega_{XPH},\omega_n)$"
    zarr = rek_upup_xph[:]
    plotK( use_pl, zarr, title )
    return

def plotUpDoKReXPH( use_pl ):
    title = r"$\operatorname{Re}K_{2,\uparrow\downarrow}(\Omega_{XPH},\omega_n)$"
    zarr = rek_updo_xph[:]
    plotK( use_pl, zarr, title )
    return

pl.suptitle(r"$U=$" + str('{0:.3f}'.format(float(U))) + r"     $\beta=$" + str('{0:.3f}'.format(float(beta))) + r"     $\Omega_{\rm PH}= \Omega_{\rm PP}=$" + str(shift) + r"$*2\pi/\beta$")

#--- Plot Physical
plotUpUpKRePH( pl.subplot(2,3,2) )
pl.xlabel(r"$\omega_n$")
plotUpDoKRePH( pl.subplot(2,3,5) )
pl.xlabel(r"$\omega_m$")
plotUpUpKRePP( pl.subplot(2,3,1) )
pl.xlabel(r"$\omega_n$")
plotUpDoKRePP( pl.subplot(2,3,4) )
pl.xlabel(r"$\omega_n$")
plotUpUpKReXPH( pl.subplot(2,3,3) )
pl.xlabel(r"$\omega_n$")
plotUpDoKReXPH( pl.subplot(2,3,6) )
pl.xlabel(r"$\omega_n$")
pl.tight_layout()

#--- Save to file
pl.savefig("plots/K_func.png", dpi = 150)
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

print re_f_upup_ph.shape

# We assume all the channels to have the same B/F grids

fgrid = f["VERT/PH/fgrid"][:].shape[0] 
bgrid = f["VERT/PH/bgrid"][:].shape[0]

print fgrid

N_bose = (bgrid-1)/2 # to create a bosonic frequency grid from -N_bose to N_bose
N_fermi= (fgrid)/2   # to create a fermionic frequency grid from -N_fermi to N_fermi

def isInside(i,j,k):
    return abs(i) <= N_bose and j >= -N_fermi and j < N_fermi and k >= -N_fermi and k < N_fermi


#-------------------------------- FUNCTION DEFINITION--------------------------------------

# 
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

def K_upup_pp(wb):
    if (abs(wb) <= N_bose_k):
        return rek_upup_pp[wb+N_bose_k]+1j*imk_upup_pp[wb+N_bose_k]
    else:

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

def P_upup_ph(wb,wf):
    if (abs(wb) <= N_bose_p and wf >= -N_fermi_p and wf < N_fermi_p):
        return rep_upup_ph[wb+N_bose_p,wf+N_fermi_p]+1j*imp_upup_ph[wb+N_bose_p,wf+N_fermi_p]
    else:
        return 0.0
def P_updo_ph(wb,wf):
    if (abs(wb) <= N_bose_p and wf >= -N_fermi_p and wf < N_fermi_p):
        return rep_updo_ph[wb+N_bose_p,wf+N_fermi_p]+1j*imp_updo_ph[wb+N_bose_p,wf+N_fermi_p]
    else:
        return 0.0

def P_upup_pp(wb,wf):
    if (abs(wb) <= N_bose_p and wf >= -N_fermi_p and wf < N_fermi_p):
        return rep_upup_pp[wb+N_bose_p,wf+N_fermi_p]+1j*imp_upup_pp[wb+N_bose_p,wf+N_fermi_p]
    else:
        return 0.0
def P_updo_pp(wb,wf):
    if (abs(wb) <= N_bose_p and wf >= -N_fermi_p and wf < N_fermi_p):
        return rep_updo_pp[wb+N_bose_p,wf+N_fermi_p]+1j*imp_updo_pp[wb+N_bose_p,wf+N_fermi_p]
    else:
        return 0.0
def P_upup_xph(wb,wf):
    if (abs(wb) <= N_bose_p and wf >= -N_fermi_p and wf < N_fermi_p):
        return rep_upup_xph[wb+N_bose_p,wf+N_fermi_p]+1j*imp_upup_xph[wb+N_bose_p,wf+N_fermi_p]
    else:
        return 0.0
def P_updo_xph(wb,wf):
    if (abs(wb) <= N_bose_p and wf >= -N_fermi_p and wf < N_fermi_p):
        return rep_updo_xph[wb+N_bose_p,wf+N_fermi_p]+1j*imp_updo_xph[wb+N_bose_p,wf+N_fermi_p]
    else:
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
        return  K_upup_pp(i) + P_upup_pp(i,j)+ P_upup_pp(i,k)+K_upup_ph(PPtoPH((i,j,k))[0]) + P_upup_ph(PPtoPH((i,j,k))[0],PPtoPH((i,j,k))[1]) + P_upup_ph(PPtoPH((i,j,k))[0],PPtoPH((i,j,k))[2])+ K_upup_xph(PPtoXPH((i,j,k))[0]) + P_upup_xph(PPtoXPH((i,j,k))[0],PPtoXPH((i,j,k))[1]) + P_upup_xph(PPtoXPH((i,j,k))[0],PPtoXPH((i,j,k))[2])

def f_updo_fun_pp(i,j,k):
    if isInside(i,j,k):
        return re_f_updo_pp[i + N_bose, j+N_fermi, k + N_fermi]+1j*im_f_updo_pp[i + N_bose, j+N_fermi, k + N_fermi]
    else:
        return - U + K_updo_pp(i) + P_updo_pp(i,j) + P_updo_pp(i,k) + K_updo_ph(PPtoPH((i,j,k))[0]) + P_updo_ph(PPtoPH((i,j,k))[0],PPtoPH((i,j,k))[1]) + P_updo_ph(PPtoPH((i,j,k))[0],PPtoPH((i,j,k))[2])+ K_updo_xph(PPtoXPH((i,j,k))[0]) + P_updo_xph(PPtoXPH((i,j,k))[0],PPtoXPH((i,j,k))[1]) + P_updo_xph(PPtoXPH((i,j,k))[0],PPtoXPH((i,j,k))[2])

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
#-------------------------------Plotting Vertex-----------------


N_fermi_plot = N_fermi

pl.rc('xtick', labelsize=9) 
pl.rc('ytick', labelsize=9) 
cmap = pl.get_cmap('jet') # Bianconiglio

pl.figsize=(13, 7)

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

pl.tight_layout()

#--- Save to file
pl.savefig("plots/Vert_original.png", dpi = 150)
#pl.savefig("plots/Vert.png", dpi = 150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()

#-------------------------------Plotting Cut Vertex-----------------


N_fermi_plot = N_fermi

def plotVert_ED( use_pl, zarr, string):
    pl.plot(np.array([(2*i+1)*pi/beta for i in range(-N_fermi_plot,N_fermi_plot)]), zarr, 'bx', ms=3, mew=0.2) 
    use_pl.set_title( string , fontsize=10)    
    return

def plotre_f_upup_ph_cut_1diag( use_pl ):
    title=r"$\operatorname{Re}F^{PH}_{\uparrow \uparrow}(\Omega_{PH},\nu_n, \nu_n)$"
    zarr = np.array([ f_upup_fun_ph(shift,n,n).real  for n in range(-N_fermi_plot,N_fermi_plot)])
    plotVert_ED( use_pl, zarr, title ) 
    return

def plotre_f_updo_ph_cut_1diag( use_pl ):
    title=r"$\operatorname{Re}F^{PH}_{\uparrow \downarrow}(\Omega_{PH},\nu_n, \nu_n)$"
    zarr = np.array([ f_updo_fun_ph(shift,n,n).real  for n in range(-N_fermi_plot,N_fermi_plot)])
    plotVert_ED( use_pl, zarr, title ) 
    return

def plotre_f_upup_ph_cut_2diag( use_pl ):
    title=r"$\operatorname{Re}F^{PH}_{\uparrow \uparrow}(\Omega_{PH},\nu_n,-\nu_n)$"
    zarr = np.array([ f_upup_fun_ph(shift,n,-n-1).real  for n in range(-N_fermi_plot,N_fermi_plot)])
    plotVert_ED( use_pl, zarr, title ) 
    return
                     
def plotre_f_updo_ph_cut_2diag( use_pl ):
    title=r"$\operatorname{Re}F^{PH}_{\uparrow \downarrow}(\Omega_{PH},\nu_n,-\nu_n)$"
    zarr = np.array([ f_updo_fun_ph(shift,n,-n-1).real  for n in range(-N_fermi_plot,N_fermi_plot)])
    plotVert_ED( use_pl, zarr, title ) 
    return

pl.suptitle(r"$U=$" + str('{0:.3f}'.format(float(U))) + r"     $\beta=$" + str('{0:.3f}'.format(float(beta))) + r"     $\Omega_{\rm PH}=$" + str(shift) + r"$*2\pi/\beta$")

plotre_f_upup_ph_cut_1diag(pl.subplot(2,2,1))
pl.xlabel(r"$\omega_n$", fontsize=20)
pl.grid()
plotre_f_upup_ph_cut_2diag(pl.subplot(2,2,2))
pl.xlabel(r"$\omega_n$", fontsize=20)
pl.grid()
plotre_f_updo_ph_cut_1diag(pl.subplot(2,2,3))
pl.xlabel(r"$\omega_n$", fontsize=20)
pl.grid()
plotre_f_updo_ph_cut_2diag(pl.subplot(2,2,4))
pl.xlabel(r"$\omega_n$", fontsize=20)
pl.grid()

#pl.tight_layout()

#--- Save to file
pl.savefig("plots/Vert_cut_PH.png", dpi = 150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()

#-------------------------------Plotting Vertex-----------------

print("Plotting Gamma ...")

#PHIs

re_gamma_upup_ph = f["GAMMA/PH/RE_GAMMA_UPUP"][:]
re_gamma_updo_ph = f["GAMMA/PH/RE_GAMMA_UPDO"][:]
re_gamma_upup_xph = f["GAMMA/XPH/RE_GAMMA_UPUP"][:]
re_gamma_updo_xph = f["GAMMA/XPH/RE_GAMMA_UPDO"][:]
re_gamma_upup_pp = f["GAMMA/PP/RE_GAMMA_UPUP"][:]
re_gamma_updo_pp = f["GAMMA/PP/RE_GAMMA_UPDO"][:]


im_gamma_upup_ph = f["GAMMA/PH/IM_GAMMA_UPUP"][:]
im_gamma_updo_ph = f["GAMMA/PH/IM_GAMMA_UPDO"][:]
im_gamma_upup_xph = f["GAMMA/XPH/IM_GAMMA_UPUP"][:]
im_gamma_updo_xph = f["GAMMA/XPH/IM_GAMMA_UPDO"][:]
im_gamma_upup_pp = f["GAMMA/PP/IM_GAMMA_UPUP"][:]
im_gamma_updo_pp = f["GAMMA/PP/IM_GAMMA_UPDO"][:]

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

def plotre_gamma_upup_ph( use_pl ):
    title=r"$\operatorname{Re}\Gamma^{PH}_{\uparrow \uparrow}$"
    zarr = re_gamma_upup_ph[shift+N_bose,:,:].real
    plotGamma_ED( use_pl, zarr, title ) 
    return

def plotre_gamma_updo_ph( use_pl ):
    title=r"$\operatorname{Re}\Gamma^{PH}_{\uparrow \downarrow}$"
    zarr = re_gamma_updo_ph[shift+N_bose,:,:].real
    plotGamma_ED( use_pl, zarr, title ) 
    return

def plotre_gamma_upup_pp( use_pl ):
    title=r"$\operatorname{Re}\Gamma^{PP}_{\uparrow \uparrow}$"
    zarr = re_gamma_upup_pp[shift+N_bose,:,:].real
    plotGamma_ED( use_pl, zarr, title ) 
    return

def plotre_gamma_updo_pp( use_pl ):
    title=r"$\operatorname{Re}\Gamma^{PP}_{\uparrow \downarrow}$"
    zarr = re_gamma_updo_pp[shift+N_bose,:,:].real
    plotGamma_ED( use_pl, zarr, title ) 
    return

def plotre_gamma_upup_xph( use_pl ):
    title=r"$\operatorname{Re}\Gamma^{XPH}_{\uparrow \uparrow}$"
    zarr = re_gamma_upup_xph[shift+N_bose,:,:].real
    plotGamma_ED( use_pl, zarr, title ) 
    return

def plotre_gamma_updo_xph( use_pl ):
    title=r"$\operatorname{Re}\Gamma^{XPH}_{\uparrow \downarrow}$"
    zarr = re_gamma_updo_xph[shift+N_bose,:,:].real
    plotGamma_ED( use_pl, zarr, title ) 
    return


pl.suptitle(r"$U=$" + str('{0:.3f}'.format(float(U))) + r"     $\beta=$" + str('{0:.3f}'.format(float(beta))) + r"     $\Omega_{\rm PP}=$" + str(shift) + r"$*2\pi/\beta$")

plotre_gamma_upup_pp(pl.subplot(2,3,1))
pl.xlabel(r"$\omega_n$", fontsize=20)
pl.ylabel(r"$\omega_m$", fontsize=20)
pl.grid()
plotre_gamma_upup_ph(pl.subplot(2,3,2))
pl.xlabel(r"$\omega_n$", fontsize=20)
pl.ylabel(r"$\omega_m$", fontsize=20)
pl.grid()
plotre_gamma_upup_xph(pl.subplot(2,3,3))
pl.xlabel(r"$\omega_n$", fontsize=20)
pl.ylabel(r"$\omega_m$", fontsize=20)
pl.grid()
plotre_gamma_updo_pp(pl.subplot(2,3,4))
pl.xlabel(r"$\omega_n$", fontsize=20)
pl.ylabel(r"$\omega_m$", fontsize=20)
pl.grid()
plotre_gamma_updo_ph(pl.subplot(2,3,5))
pl.xlabel(r"$\omega_n$", fontsize=20)
pl.ylabel(r"$\omega_m$", fontsize=20)
pl.grid()
plotre_gamma_updo_xph(pl.subplot(2,3,6))
pl.xlabel(r"$\omega_n$", fontsize=20)
pl.ylabel(r"$\omega_m$", fontsize=20)
pl.grid()

pl.tight_layout()

#--- Save to file
pl.savefig("plots/Gamma.png", dpi = 150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()

#--------------------------------------------PHI PLOTTING------------------------------------------
print("Plotting Phi ...")

#PHIs

re_phi_upup_ph = f["PHI/PH/RE_PHI_UPUP"][:]
re_phi_updo_ph = f["PHI/PH/RE_PHI_UPDO"][:]
re_phi_upup_xph = f["PHI/XPH/RE_PHI_UPUP"][:]
re_phi_updo_xph = f["PHI/XPH/RE_PHI_UPDO"][:]
re_phi_upup_pp = f["PHI/PP/RE_PHI_UPUP"][:]
re_phi_updo_pp = f["PHI/PP/RE_PHI_UPDO"][:]


im_phi_upup_ph = f["PHI/PH/IM_PHI_UPUP"][:]
im_phi_updo_ph = f["PHI/PH/IM_PHI_UPDO"][:]
im_phi_upup_xph = f["PHI/XPH/IM_PHI_UPUP"][:]
im_phi_updo_xph = f["PHI/XPH/IM_PHI_UPDO"][:]
im_phi_upup_pp = f["PHI/PP/IM_PHI_UPUP"][:]
im_phi_updo_pp = f["PHI/PP/IM_PHI_UPDO"][:]

# We assume the phi fucntions to have the same b and f grids
#PHIS FUNCTIONS


def phi_upup_fun_ph(i,j,k):
    if isInside(i,j,k):
        return re_phi_upup_ph[i + N_bose, j+N_fermi, k + N_fermi]+1j*im_phi_upup_ph[i + N_bose, j+N_fermi, k + N_fermi]
    else:
        return K_upup_ph(i) + P_upup_ph(i,j) +P_upup_ph(i,k) 
    
def phi_updo_fun_ph(i,j,k):
    if isInside(i,j,k):
        return re_phi_updo_ph[i + N_bose, j+N_fermi, k + N_fermi]+1j*im_phi_updo_ph[i + N_bose, j+N_fermi, k + N_fermi]
    else:
        return K_updo_ph(i) + P_updo_ph(i,j)+P_updo_ph(i,k)
    
def phi_upup_fun_pp(i,j,k):
    if isInside(i,j,k):
        return re_phi_upup_pp[i + N_bose, j+N_fermi, k + N_fermi]+1j*im_phi_updo_pp[i + N_bose, j+N_fermi, k + N_fermi]
    else:
        return K_upup_pp(i) + P_upup_pp(i,j) +P_upup_pp(i,k) 
    
def phi_updo_fun_pp(i,j,k):
    if isInside(i,j,k):
        return re_phi_updo_pp[i + N_bose, j+N_fermi, k + N_fermi]+1j*im_phi_updo_pp[i + N_bose, j+N_fermi, k + N_fermi]
    else:
        return K_updo_pp(i) + P_updo_pp(i,j) + P_updo_pp(i,k)
    
def phi_upup_fun_xph(i,j,k):
    if isInside(i,j,k):
        return re_phi_upup_xph[i + N_bose, j+N_fermi, k + N_fermi]+1j*im_phi_upup_xph[i + N_bose, j+N_fermi, k + N_fermi]
    else:
        return K_upup_xph(i) + P_upup_xph(i,j)  + P_upup_xph(i,k)
    
def phi_updo_fun_xph(i,j,k):
    if isInside(i,j,k):
        return re_phi_updo_xph[i + N_bose, j+N_fermi, k + N_fermi]+1j*im_phi_updo_xph[i + N_bose, j+N_fermi, k + N_fermi]
    else:
        return K_updo_xph(i) + P_updo_xph(i,j)+ P_updo_xph(i,k)
    
#-------------------------------Plotting Phi-----------------

N_fermi_plot = 2*N_fermi

pl.rc('xtick', labelsize=9) 
pl.rc('ytick', labelsize=9) 
cmap = pl.get_cmap('jet') # Bianconiglio

pl.figsize=(13, 7)

def plotPhi_ED( use_pl, zarr, string):
    use_pl.set_aspect(1.0)
    pl.pcolormesh(np.array([(2*i+1)*pi/beta for i in range(-N_fermi_plot,N_fermi_plot)]), np.array([(2*i+1)*pi/beta for i in range(-N_fermi_plot,N_fermi_plot)]),np.ma.masked_where( np.isnan(zarr), zarr ))
    use_pl.set_title( string , fontsize=10)    
    pl.colorbar(shrink=0.6)
    return

def plotre_phi_upup_ph( use_pl ):
    title=r"$\operatorname{Re}\Phi^{PH}_{\uparrow \uparrow}$"
    zarr = np.array([[ phi_upup_fun_ph(shift,n,m).real  for n in range(-N_fermi_plot,N_fermi_plot)] for m in range(-N_fermi_plot,N_fermi_plot)])
    plotPhi_ED( use_pl, zarr, title ) 
    return

def plotre_phi_updo_ph( use_pl ):
    title=r"$\operatorname{Re}\Phi^{PH}_{\uparrow \downarrow}$"
    zarr = np.array([[ phi_updo_fun_ph(shift,n,m).real  for n in range(-N_fermi_plot,N_fermi_plot)] for m in range(-N_fermi_plot,N_fermi_plot)])
    plotPhi_ED( use_pl, zarr, title ) 
    return

def plotre_phi_upup_pp( use_pl ):
    title=r"$\operatorname{Re}\Phi^{PP}_{\uparrow \uparrow}$"
    zarr = np.array([[ phi_upup_fun_pp(shift,n,m).real for n in range(-N_fermi_plot,N_fermi_plot)] for m in range(-N_fermi_plot,N_fermi_plot)])
    plotPhi_ED( use_pl, zarr, title ) 
    return

def plotre_phi_updo_pp( use_pl ):
    title=r"$\operatorname{Re}\Phi^{PP}_{\uparrow \downarrow}$"
    zarr = np.array([[ phi_updo_fun_pp(shift,n,m).real for n in range(-N_fermi_plot,N_fermi_plot)] for m in range(-N_fermi_plot,N_fermi_plot)])
    plotPhi_ED( use_pl, zarr, title ) 
    return

def plotre_phi_upup_xph( use_pl ):
    title=r"$\operatorname{Re}\Phi^{XPH}_{2\uparrow \uparrow}$"
    zarr = np.array([[ (phi_upup_fun_xph(shift,n,m)).real  for n in range(-N_fermi_plot,N_fermi_plot)] for m in range(-N_fermi_plot,N_fermi_plot)])
    plotPhi_ED( use_pl, zarr, title ) 
    return

def plotre_phi_updo_xph( use_pl ):
    title=r"$\operatorname{Re}\Phi^{XPH}_{2\uparrow \downarrow}$"
    zarr = np.array([[ (phi_updo_fun_xph(shift,n,m)).real  for n in range(-N_fermi_plot,N_fermi_plot)] for m in range(-N_fermi_plot,N_fermi_plot)])
    plotPhi_ED( use_pl, zarr, title ) 
    return

pl.suptitle(r"$U=$" + str('{0:.3f}'.format(float(U))) + r"     $\beta=$" + str('{0:.3f}'.format(float(beta))) + r"     $\Omega_{\rm PP}=$" + str(shift) + r"$*2\pi/\beta$")

plotre_phi_upup_pp(pl.subplot(2,3,1))
pl.xlabel(r"$\omega_n$", fontsize=20)
pl.ylabel(r"$\omega_m$", fontsize=20)
pl.grid()
plotre_phi_upup_ph(pl.subplot(2,3,2))
pl.xlabel(r"$\omega_n$", fontsize=20)
pl.ylabel(r"$\omega_m$", fontsize=20)
pl.grid()
plotre_phi_upup_xph(pl.subplot(2,3,3))
pl.xlabel(r"$\omega_n$", fontsize=20)
pl.ylabel(r"$\omega_m$", fontsize=20)
pl.grid()
plotre_phi_updo_pp(pl.subplot(2,3,4))
pl.xlabel(r"$\omega_n$", fontsize=20)
pl.ylabel(r"$\omega_m$", fontsize=20)
pl.grid()
plotre_phi_updo_ph(pl.subplot(2,3,5))
pl.xlabel(r"$\omega_n$", fontsize=20)
pl.ylabel(r"$\omega_m$", fontsize=20)
pl.grid()
plotre_phi_updo_xph(pl.subplot(2,3,6))
pl.xlabel(r"$\omega_n$", fontsize=20)
pl.ylabel(r"$\omega_m$", fontsize=20)
pl.grid()

pl.tight_layout()

#--- Save to file
pl.savefig("plots/Phi.png", dpi = 150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()

#-------------------------------------------------LAMBDA--------------------------------------------

print("Plotting Lambda ...")

#LAMBDAS

re_lambda_upup_ph = f["Lambda/PH/RE_LAMBDA_UPUP"][:]
re_lambda_updo_ph = f["Lambda/PH/RE_LAMBDA_UPDO"][:]
re_lambda_upup_xph = f["Lambda/XPH/RE_LAMBDA_UPUP"][:]
re_lambda_updo_xph = f["Lambda/XPH/RE_LAMBDA_UPDO"][:]
re_lambda_upup_pp = f["Lambda/PP/RE_LAMBDA_UPUP"][:]
re_lambda_updo_pp = f["Lambda/PP/RE_LAMBDA_UPDO"][:]


im_lambda_upup_ph = f["Lambda/PH/IM_LAMBDA_UPUP"][:]
im_lambda_updo_ph = f["Lambda/PH/IM_LAMBDA_UPDO"][:]
im_lambda_upup_xph = f["Lambda/XPH/IM_LAMBDA_UPUP"][:]
im_lambda_updo_xph = f["Lambda/XPH/IM_LAMBDA_UPDO"][:]
im_lambda_upup_pp = f["Lambda/PP/IM_LAMBDA_UPUP"][:]
im_lambda_updo_pp = f["Lambda/PP/IM_LAMBDA_UPDO"][:]


#-------------------------------Plotting Lambda-----------------

N_fermi_plot = N_fermi

pl.rc('xtick', labelsize=9) 
pl.rc('ytick', labelsize=9) 
cmap = pl.get_cmap('jet') # Bianconiglio

pl.figsize=(13, 7)

def plotLambda_ED( use_pl, zarr, string):
    use_pl.set_aspect(1.0)
    pl.pcolormesh(np.array([(2*i+1)*pi/beta for i in range(-N_fermi_plot,N_fermi_plot)]), np.array([(2*i+1)*pi/beta for i in range(-N_fermi_plot,N_fermi_plot)]),np.ma.masked_where( np.isnan(zarr), zarr ))
    use_pl.set_title( string , fontsize=10)    
    pl.colorbar(shrink=0.6)
    return

def plotre_lambda_upup_ph( use_pl ):
    title=r"$\operatorname{Re}\Lambda^{PH}_{\uparrow \uparrow}$"
    zarr =  re_lambda_upup_ph[:,:,shift+N_bose]
    plotLambda_ED( use_pl, zarr, title ) 
    return

def plotre_lambda_updo_ph( use_pl ):
    title=r"$\operatorname{Re}\Lambda^{PH}_{\uparrow \downarrow}$"
    zarr =  re_lambda_updo_ph[:,:,shift+N_bose]
    plotLambda_ED( use_pl, zarr, title ) 
    return

def plotre_lambda_upup_pp( use_pl ):
    title=r"$\operatorname{Re}\Lambda^{PP}_{\uparrow \uparrow}$"
    zarr =  re_lambda_upup_pp[:,:,shift+N_bose]
    plotLambda_ED( use_pl, zarr, title ) 
    return

def plotre_lambda_updo_pp( use_pl ):
    title=r"$\operatorname{Re}\Lambda^{PP}_{\uparrow \downarrow}$"
    zarr =  re_lambda_updo_pp[:,:,shift+N_bose]
    plotLambda_ED( use_pl, zarr, title ) 
    return

def plotre_lambda_upup_xph( use_pl ):
    title=r"$\operatorname{Re}\Lambda^{XPH}_{2\uparrow \uparrow}$"
    zarr =  re_lambda_upup_xph[:,:,shift+N_bose]
    plotLambda_ED( use_pl, zarr, title ) 
    return

def plotre_lambda_updo_xph( use_pl ):
    title=r"$\operatorname{Re}\Lambda^{XPH}_{2\uparrow \downarrow}$"
    zarr =  re_lambda_updo_xph[:,:,shift+N_bose] 
    plotLambda_ED( use_pl, zarr, title ) 
    return

pl.suptitle(r"$U=$" + str('{0:.3f}'.format(float(U))) + r"     $\beta=$" + str('{0:.3f}'.format(float(beta))) + r"     $\Omega_{\rm PP}=$" + str(shift) + r"$*2\pi/\beta$")

plotre_lambda_upup_pp(pl.subplot(2,3,1))
pl.xlabel(r"$\omega_n$", fontsize=20)
pl.ylabel(r"$\omega_m$", fontsize=20)
pl.grid()
plotre_lambda_upup_ph(pl.subplot(2,3,2))
pl.xlabel(r"$\omega_n$", fontsize=20)
pl.ylabel(r"$\omega_m$", fontsize=20)
pl.grid()
plotre_lambda_upup_xph(pl.subplot(2,3,3))
pl.xlabel(r"$\omega_n$", fontsize=20)
pl.ylabel(r"$\omega_m$", fontsize=20)
pl.grid()
plotre_lambda_updo_pp(pl.subplot(2,3,4))
pl.xlabel(r"$\omega_n$", fontsize=20)
pl.ylabel(r"$\omega_m$", fontsize=20)
pl.grid()
plotre_lambda_updo_ph(pl.subplot(2,3,5))
pl.xlabel(r"$\omega_n$", fontsize=20)
pl.ylabel(r"$\omega_m$", fontsize=20)
pl.grid()
plotre_lambda_updo_xph(pl.subplot(2,3,6))
pl.xlabel(r"$\omega_n$", fontsize=20)
pl.ylabel(r"$\omega_m$", fontsize=20)
pl.grid()

pl.tight_layout()

#--- Save to file
pl.savefig("plots/Lambda.png", dpi = 150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()


#-------------------------------------------------REST FUNCTION--------------------------------------------

print("Plotting R_func ...")

#RS

re_R_upup_ph = f["R_func/PH/RE_R_UPUP"][:]
re_R_updo_ph = f["R_func/PH/RE_R_UPDO"][:]
re_R_upup_xph = f["R_func/XPH/RE_R_UPUP"][:]
re_R_updo_xph = f["R_func/XPH/RE_R_UPDO"][:]
re_R_upup_pp = f["R_func/PP/RE_R_UPUP"][:]
re_R_updo_pp = f["R_func/PP/RE_R_UPDO"][:]


im_R_upup_ph = f["R_func/PH/IM_R_UPUP"][:]
im_R_updo_ph = f["R_func/PH/IM_R_UPDO"][:]
im_R_upup_xph = f["R_func/XPH/IM_R_UPUP"][:]
im_R_updo_xph = f["R_func/XPH/IM_R_UPDO"][:]
im_R_upup_pp = f["R_func/PP/IM_R_UPUP"][:]
im_R_updo_pp = f["R_func/PP/IM_R_UPDO"][:]


#-------------------------------Plotting R_func-----------------

N_fermi_plot = N_fermi

pl.rc('xtick', labelsize=9) 
pl.rc('ytick', labelsize=9) 
cmap = pl.get_cmap('jet') # Bianconiglio

pl.figsize=(13, 7)

def plotR_func_ED( use_pl, zarr, string):
    use_pl.set_aspect(1.0)
    pl.pcolormesh(np.array([(2*i+1)*pi/beta for i in range(-N_fermi_plot,N_fermi_plot)]), np.array([(2*i+1)*pi/beta for i in range(-N_fermi_plot,N_fermi_plot)]),np.ma.masked_where( np.isnan(zarr), zarr ))
    use_pl.set_title( string , fontsize=10)    
    pl.colorbar(shrink=0.6)
    return

def plotre_R_upup_ph( use_pl ):
    title=r"$\operatorname{Re}R^{PH}_{\uparrow \uparrow}$"
    zarr =  re_R_upup_ph[:,:,shift+N_bose] 
    plotR_func_ED( use_pl, zarr, title ) 
    return

def plotre_R_updo_ph( use_pl ):
    title=r"$\operatorname{Re}R^{PH}_{\uparrow \downarrow}$"
    zarr =  re_R_updo_ph[:,:,shift+N_bose]
    plotR_func_ED( use_pl, zarr, title ) 
    return

def plotre_R_upup_pp( use_pl ):
    title=r"$\operatorname{Re}R^{PP}_{\uparrow \uparrow}$"
    zarr =  re_R_upup_pp[:,:,shift+N_bose]
    plotR_func_ED( use_pl, zarr, title ) 
    return

def plotre_R_updo_pp( use_pl ):
    title=r"$\operatorname{Re}R^{PP}_{\uparrow \downarrow}$"
    zarr =  re_R_updo_pp[:,:,shift+N_bose]    
    plotR_func_ED( use_pl, zarr, title ) 
    return

def plotre_R_upup_xph( use_pl ):
    title=r"$\operatorname{Re}R^{XPH}_{2\uparrow \uparrow}$"
    zarr =  re_R_upup_xph[:,:,shift+N_bose]
    plotR_func_ED( use_pl, zarr, title ) 
    return

def plotre_R_updo_xph( use_pl ):
    title=r"$\operatorname{Re}R^{XPH}_{2\uparrow \downarrow}$"
    zarr =  re_R_updo_xph[:,:,shift+N_bose]
    plotR_func_ED( use_pl, zarr, title ) 
    return

pl.suptitle(r"$U=$" + str('{0:.3f}'.format(float(U))) + r"     $\beta=$" + str('{0:.3f}'.format(float(beta))) + r"     $\Omega_{\rm PP}=$" + str(shift) + r"$*2\pi/\beta$")

plotre_R_upup_pp(pl.subplot(2,3,1))
pl.xlabel(r"$\omega_n$", fontsize=20)
pl.ylabel(r"$\omega_m$", fontsize=20)
pl.grid()
plotre_R_upup_ph(pl.subplot(2,3,2))
pl.xlabel(r"$\omega_n$", fontsize=20)
pl.ylabel(r"$\omega_m$", fontsize=20)
pl.grid()
plotre_R_upup_xph(pl.subplot(2,3,3))
pl.xlabel(r"$\omega_n$", fontsize=20)
pl.ylabel(r"$\omega_m$", fontsize=20)
pl.grid()
plotre_R_updo_pp(pl.subplot(2,3,4))
pl.xlabel(r"$\omega_n$", fontsize=20)
pl.ylabel(r"$\omega_m$", fontsize=20)
pl.grid()
plotre_R_updo_ph(pl.subplot(2,3,5))
pl.xlabel(r"$\omega_n$", fontsize=20)
pl.ylabel(r"$\omega_m$", fontsize=20)
pl.grid()
plotre_R_updo_xph(pl.subplot(2,3,6))
pl.xlabel(r"$\omega_n$", fontsize=20)
pl.ylabel(r"$\omega_m$", fontsize=20)
pl.grid()

pl.tight_layout()

#--- Save to file
pl.savefig("plots/R_func.png", dpi = 150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()


