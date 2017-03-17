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

fname = "dat/dat_U1_Beta20_PFCB120_PARQ_SU2_METH2.h5"

if len(sys.argv) > 1:
    fname = str(sys.argv[1])

fname = fname.rstrip('\n') # strip newline of fname
f = h5py.File(fname, "r")

os.system('mkdir -p log')
os.system('mkdir -p plots')

#os.system('rm plots/prev_Vert.png 2> log/plot.log')
#os.system('rm plots/prev_phi.png 2> log/plot.log')
#os.system('rm plots/prev_P_func.png 2> log/plot.log')
#os.system('rm plots/prev_R_func.png 2> log/plot.log')

#os.system('mv plots/Vert.png plots/prev_Vert.png 2> log/plot.log')
#os.system('mv plots/phi.png plots/prev_phi.png 2> log/plot.log')
#os.system('mv plots/P_func.png plots/prev_P_func.png 2> log/plot.log')
#os.system('mv plots/R_func.png plots/prev_R_func.png 2> log/plot.log')


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

#--------------------------------------SELF ENERGY PLOTTING ------------------------------------------

print("Plotting self-energy ...")


#--- Read
siggrid = np.array(f["/Sig/fgrid"])
vertgrid = np.array(f["/Vert/fgrid"])
resig = f["/Sig/RE"]
imsig = f["/Sig/IM"]
fdim = resig.shape[0]


#--- Helper functions
def plotSig( use_pl, arr, string ):
    pl.plot( siggrid, arr[:,0,0,0], 'bx', ms=3, mew=0.2)
    pl.xlim([min(siggrid),max(siggrid)])
    use_pl.set_title(string)
    return

pl.suptitle(r"$U=$" + str(UINT) + r"     $\beta=$" + str(BETA) + r"     $\epsilon=$" + str(EPS))

#--- Plot physical
plotSig( pl.subplot(2,2,1), resig, RE + "\Sigma(i\omega)$" ) 
plotSig( pl.subplot(2,2,2), imsig, IM + "\Sigma(i\omega)$" ) 

pl.tight_layout()

#--- Save to file
pl.savefig("plots/Sig.png", dpi=150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()

#--------------------------------------G(iw) PLOTTING ------------------------------------------

print("Plotting Green function ...")


#--- Read
Giwgrid = np.array(f["/Giw/fgrid"])
reGiw = f["/Giw/RE"]
imGiw = f["/Giw/IM"]
fdim = reGiw.shape[0]


#--- Helper functions
def plotGiw( use_pl, arr, string ):
    pl.plot( Giwgrid, arr[:,0,0,0], 'bx', ms=3, mew=0.2)
    pl.xlim([min(Giwgrid),max(Giwgrid)])
    use_pl.set_title(string)
    return

pl.suptitle(r"$U=$" + str(UINT) + r"     $\beta=$" + str(BETA) + r"     $\epsilon=$" + str(EPS))

#--- Plot physical
plotGiw( pl.subplot(2,2,1), reGiw, RE + "G(i\omega)$" ) 
pl.xlabel(r"$\omega_n$")
plotGiw( pl.subplot(2,2,2), imGiw, IM + "G(i\omega)$" ) 
pl.xlabel(r"$\omega_n$")

pl.tight_layout()

#--- Save to file
pl.savefig("plots/Giw.png", dpi=150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()

#--------------------------------------FLOW OBSERVABLES ------------------------------------------

#print("Plotting flowing observables ...")

##--- Read
#Lam_arr = np.array(f["/Flow_obs/LAM"])
#meff_arr = np.array(f["/Flow_obs/EFF_MASS"])
#meff_err_arr = np.array(f["/Flow_obs/ERR_EFF_MASS"])
#maxCpl_arr = np.array(f["/Flow_obs/ABS_MAX_CPL"])

#def TK( U ):
    #return np.sqrt(U/2.0)*np.exp(-math.pi*U/8.0)

#nrg_x = [ 0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20 ]
#nrg_y = [ 1, 1.01, 1.06, 1.14, 1.26, 1.42, 1.62, 1.88, 2.18, 2.56, 3.03, 4.32, 6.14, 8.74, 12.47, 17.78, 25.5, 36.5, 52.3, 74.9, 109.4, 158, 228, 328, 474, 681]

#def plotEffM( use_pl ):
    #pl.plot( Lam_arr, meff_arr, 'bx', ms=3, mew=0.2)
    #use_pl.set_title( r"$m^*$" )
    ##pl.plot( Lam_arr**2 * UINT, meff_arr, 'bx', ms=3, mew=0.2)
    ##pl.plot( nrg_x, nrg_y, 'r--', ms=3, mew=0.2 )
    ##pl.xlim([0,UINT])
    ##pl.ylim([1.0,1.0/TK(UINT)])
    ##pl.yscale('log')
    #return

#def plotObs( use_pl, y_arr, string ):
    #pl.plot( Lam_arr, y_arr, 'bx', ms=3, mew=0.2)
    #use_pl.set_title( string )
    ##pl.plot( Lam_arr**2 * UINT, y_arr, 'bx', ms=3, mew=0.2)
    ##pl.yscale('log')
    #return

#plotEffM( pl.subplot(2,1,1) )
#pl.xlabel(r"$\Lambda$")
#plotObs( pl.subplot(2,1,2), maxCpl_arr, r"Max Cpl")
#pl.xlabel(r"$\Lambda$")

#pl.tight_layout()

#pl.savefig("plots/flow_obs.png", dpi=150)
#pl.figure(dpi=100) # Reset dpi to default
#pl.clf()

#--------------------------------------VERTEX PLOTTING ------------------------------------------

print("Plotting vertex ...")

#--- Read
revert = vert_mul * np.array(f["/Vert/RE"])
imvert = vert_mul * np.array(f["/Vert/IM"])
fdim = revert.shape[0]

if fdim <= shift:
    sys.exit("Error: Shift to large for vertex grid"); 

#---  Helper functions
def neg( w ):
    return fdim - w - 1

def check_bounds( w1, w2, w1p ):
    if ( w1 < 0 or w1 > fdim - 1 or w2 < 0 or w2 > fdim - 1 or w1p < 0 or w1p > fdim - 1 ):
        return False
    return True

def ReVert( w1, w2, w1p, i, j, k, l ):
    if ( not check_bounds( w1, w2, w1p ) ):
            return float('nan')
    return revert[w1,w2,w1p,0,0,0,0,0,0,0]

def ImVert( w1, w2, w1p, i, j, k, l ):
    if ( not check_bounds( w1, w2, w1p ) ):
            return float('nan')
    return imvert[w1,w2,w1p,0,0,0,0,0,0,0]

def ReVertUpDown( w1, w2, w1p ):
    if ( not check_bounds( w1, neg(w2), w1p ) ):
            return float('nan')
    return revert[w1, w2, w1p, 0, 0, 0,0,0,0,0]

def ImVertUpDown( w1, w2, w1p ):
    if ( not check_bounds( w1, neg(w2), w1p ) ):
            return float('nan')
    return imvert[w1, w2, w1p, 0, 0, 0,0,0,0,0]

def ReVertUpUp( w1, w2, w1p ):
    if ( not check_bounds( w1, w2, w1p ) ):
            return float('nan')
    return ReVertUpDown( w1, w2, w1p ) - ReVertUpDown( w1, w2, w1+w2-w1p-1 )

def ImVertUpUp( w1, w2, w1p ):
    if ( not check_bounds( w1, w2, w1p ) ):
            return float('nan')
    return ImVertUpDown( w1, w2, w1p ) - ImVertUpDown( w1, w2, w1+w2-w1p-1 )


def plotVert( use_pl, zarr, string ):
    use_pl.set_aspect(1.0)
    pl.pcolormesh( vertgrid, vertgrid, np.ma.masked_where( np.isnan(zarr), zarr ) )
    pl.ylim([min(vertgrid),max(vertgrid)])
    pl.xlim([min(vertgrid),max(vertgrid)])
    use_pl.set_title( string , fontsize=10)
    pl.colorbar(shrink=0.6) 
    return

def plotUpUpVertRePP( use_pl ):
    title = RE + "\gamma_{2,\uparrow\uparrow}(\omega_n,\Omega_{PP}-\omega_n,\omega_m)$"
    zarr = np.array([[ ReVertUpUp(n,shift+neg(n),shift+neg(m)) for n in range(fdim)] for m in range(fdim)])
    plotVert( use_pl, zarr, title )
    return

def plotUpDownVertRePP( use_pl ):
    title = RE + "\gamma_{2,\uparrow\downarrow}(\omega_n,\Omega_{PP}-\omega_n,\omega_m)$"
    zarr = np.array([[ ReVertUpDown(n,shift+neg(n),shift+neg(m)) for n in range(fdim)] for m in range(fdim)])
    plotVert( use_pl, zarr, title )
    return

def plotUpUpVertRePH( use_pl ):
    title = RE + "\gamma_{2,\uparrow\uparrow}(\omega_n,\Omega_{PH}+\omega_m,\Omega_{PH}+\omega_n)$"
    zarr = np.array([[ ReVertUpUp(n,m+shift,n+shift) for n in range(fdim)] for m in range(fdim)])
    plotVert( use_pl, zarr, title )
    return

def plotUpDownVertRePH( use_pl ):
    title = RE + "\gamma_{2,\uparrow\downarrow}(\omega_n,\Omega_{PH}+\omega_m,\Omega_{PH}+\omega_n)$"
    zarr = np.array([[ ReVertUpDown(n,m+shift,n+shift) for n in range(fdim)] for m in range(fdim)])
    plotVert( use_pl, zarr, title )
    return

def plotUpUpVertReXPH( use_pl ):
    title = RE + "\gamma_{2,\uparrow\uparrow}(\omega_n,\Omega_{XPH}+\omega_m,\omega_m)$"
    zarr = np.array([[ ReVertUpUp(n,m+shift,m) for n in range(fdim)] for m in range(fdim)])
    plotVert( use_pl, zarr, title )
    return

def plotUpDownVertReXPH( use_pl ):
    title = RE + "\gamma_{2,\uparrow\downarrow}(\omega_n,\Omega_{XPH}+\omega_m,\omega_m)$"
    zarr = np.array([[ ReVertUpDown(n,m+shift,m) for n in range(fdim)] for m in range(fdim)])
    plotVert( use_pl, zarr, title )
    return

#--- Plot
pl.suptitle(r"$U=$" + str(UINT) + r"     $\beta=$" + str(BETA) + r"     $\epsilon=$" + str(EPS) +  r"     $\Omega_{\rm PP}=\Omega_{\rm PH}=\Omega_{\rm xPH}=$" + str(shift) + r"$*2\pi/\beta$" + r"     Notation: $\gamma_{2}(\omega_1,\omega_2,\omega_1')$")

plotUpUpVertRePP( pl.subplot(2,3,1) )
pl.ylabel(r"$\omega_m$")
plotUpUpVertRePH( pl.subplot(2,3,2) )
plotUpUpVertReXPH( pl.subplot(2,3,3) )
plotUpDownVertRePP( pl.subplot(2,3,4) )
pl.xlabel(r"$\omega_n$")
pl.ylabel(r"$\omega_m$")
plotUpDownVertRePH( pl.subplot(2,3,5) )
pl.xlabel(r"$\omega_n$")
plotUpDownVertReXPH( pl.subplot(2,3,6) )
pl.xlabel(r"$\omega_n$")
pl.tight_layout()

#--- Save to file
pl.savefig("plots/Vert.png", dpi = 150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()

shift=6

#--- Plot
pl.suptitle(r"$U=$" + str(UINT) + r"     $\beta=$" + str(BETA) + r"     $\epsilon=$" + str(EPS) +  r"     $\Omega_{\rm PP}=\Omega_{\rm PH}=\Omega_{\rm xPH}=$" + str(shift) + r"$*2\pi/\beta$" + r"     Notation: $\gamma_{2}(\omega_1,\omega_2,\omega_1')$")

#--- Plot Physical
plotUpUpVertRePP( pl.subplot(2,3,1) )
pl.ylabel(r"$\omega_m$")
plotUpUpVertRePH( pl.subplot(2,3,2) )
plotUpUpVertReXPH( pl.subplot(2,3,3) )

plotUpDownVertRePP( pl.subplot(2,3,4) )
pl.xlabel(r"$\omega_n$")
pl.ylabel(r"$\omega_m$")
plotUpDownVertRePH( pl.subplot(2,3,5) )
pl.xlabel(r"$\omega_n$")
plotUpDownVertReXPH( pl.subplot(2,3,6) )
pl.xlabel(r"$\omega_n$")

pl.tight_layout()

#--- Save to file
pl.savefig("plots/Vert_shift.png", dpi = 150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()

shift=0

#--------------------------------------LAMBDA PLOTTING ------------------------------------------

print("Plotting Lambda ...")

#--- Read
relambda = vert_mul * np.array(f["/Lambda/RE"])
imlambda = vert_mul * np.array(f["/Lambda/IM"])
fdim = relambda.shape[0]
lambdagrid = np.array(f["/Lambda/fgrid"])

if fdim <= shift:
    sys.exit("Error: Shift to large for vertex grid"); 

#---  Helper functions
def neg( w ):
    return fdim - w - 1

def check_bounds_mix( W, w, wp ):
    if ( W < 0 or W > bdim or w < 0 or w > fdim - 1 or wp < 0 or wp > fdim - 1 ):
        return False
    return True

def check_bounds( w1, w2, w1p ):
    if ( w1 < 0 or w1 > fdim - 1 or w2 < 0 or w2 > fdim - 1 or w1p < 0 or w1p > fdim - 1 ):
        return False
    return True

def ReLambda( w1, w2, w1p, i, j, k, l ):
    if ( not check_bounds( w1, w2, w1p ) ):
            return float('nan')
    return relambda[w1,w2,w1p,0,0,0,0,0,0,0]

def ImLambda( w1, w2, w1p, i, j, k, l ):
    if ( not check_bounds( w1, w2, w1p ) ):
            return float('nan')
    return imlambda[w1,w2,w1p,0,0,0,0,0,0,0]

def ReLambdaUpDown( w1, w2, w1p ):
    if ( not check_bounds( w1, neg(w2), w1p ) ):
            return float('nan')
    return relambda[w1, w2, w1p, 0, 0, 0,0,0,0,0]

def ImLambdaUpDown( w1, w2, w1p ):
    if ( not check_bounds( w1, neg(w2), w1p ) ):
            return float('nan')
    return imlambda[w1, w2, w1p, 0, 0, 0,0,0,0,0]

def ReLambdaUpUp( w1, w2, w1p ):
    if ( not check_bounds( w1, w2, w1p ) ):
            return float('nan')
    return ReLambdaUpDown( w1, w2, w1p ) - ReLambdaUpDown( w2, w1, w1p )

def ImLambdaUpUp( w1, w2, w1p ):
    if ( not check_bounds( w1, w2, w1p ) ):
            return float('nan')
    return ImLambdaUpDown( w1, w2, w1p ) - ImLambdaUpDown( w2, w1, w1p )


def plotLambda( use_pl, zarr, string ):
    use_pl.set_aspect(1.0)
    pl.pcolormesh( lambdagrid, lambdagrid, np.ma.masked_where( np.isnan(zarr), zarr ) )
    pl.ylim([min(lambdagrid),max(lambdagrid)])
    pl.xlim([min(lambdagrid),max(lambdagrid)])
    use_pl.set_title( string , fontsize=10)
    pl.colorbar(shrink=0.6) 
    return

def plotUpUpLambdaRePP( use_pl ):
    title = r"$\operatorname{Re}\Lambda_{\uparrow\uparrow}(\omega_n,\Omega_{PP}-\omega_n,\omega_m)$"
    zarr = np.array([[ ReLambdaUpUp(n,shift+neg(n),shift+neg(m)) for n in range(fdim)] for m in range(fdim)])
    plotLambda( use_pl, zarr, title )
    return

def plotUpDownLambdaRePP( use_pl ):
    title = r"$\operatorname{Re}\Lambda_{\uparrow\downarrow}(\omega_n,\Omega_{PP}-\omega_n,\omega_m)$"
    zarr = np.array([[ ReLambdaUpDown(n,shift+neg(n),shift+neg(m)) for n in range(fdim)] for m in range(fdim)])
    plotLambda( use_pl, zarr, title )
    return

def plotUpUpLambdaRePH( use_pl ):
    title = r"$\operatorname{Re}\Lambda_{\uparrow\uparrow}(\omega_n,\Omega_{PH}+\omega_m,\Omega_{PH}+\omega_n)$"
    zarr = np.array([[ ReLambdaUpUp(n,m+shift,n+shift) for n in range(fdim)] for m in range(fdim)])
    plotLambda( use_pl, zarr, title )
    return

def plotUpDownLambdaRePH( use_pl ):
    title = r"$\operatorname{Re}\Lambda_{\uparrow\downarrow}(\omega_n,\Omega_{PH}+\omega_m,\Omega_{PH}+\omega_n)$"
    zarr = np.array([[ ReLambdaUpDown(n,m+shift,n+shift) for n in range(fdim)] for m in range(fdim)])
    plotLambda( use_pl, zarr, title )
    return

def plotUpUpLambdaRePHX( use_pl ):
    title = r"$\operatorname{Re}\Lambda_{\uparrow\uparrow}(\omega_n,\Omega_{PHX}+\omega_m,\omega_m)$"
    zarr = np.array([[ ReLambdaUpUp(n,m+shift,m) for n in range(fdim)] for m in range(fdim)])
    plotLambda( use_pl, zarr, title )
    return

def plotUpDownLambdaRePHX( use_pl ):
    title = r"$\operatorname{Re}\Lambda_{\uparrow\downarrow}(\omega_n,\Omega_{PHX}+\omega_m,\omega_m)$"
    zarr = np.array([[ ReLambdaUpDown(n,m+shift,m) for n in range(fdim)] for m in range(fdim)])
    plotLambda( use_pl, zarr, title )
    return


#--- plot
pl.suptitle(r"$U=$" + str(UINT) + r"     $\beta=$" + str(BETA) + r"     $\epsilon=$" + str(EPS) +  r"     $\Omega_{\rm PP}=\Omega_{\rm PH}=\Omega_{\rm xPH}=$" + str(shift) + r"$*2\pi/\beta$" + r"     Notation: $\Lambda(\omega_1,\omega_2,\omega_1')$")

plotUpUpLambdaRePP( pl.subplot(2,3,1) )
pl.ylabel(r"$\omega_m$")
plotUpUpLambdaRePH( pl.subplot(2,3,2) )
plotUpUpLambdaRePHX( pl.subplot(2,3,3) )
plotUpDownLambdaRePP( pl.subplot(2,3,4) )
pl.xlabel(r"$\omega_n$")
pl.ylabel(r"$\omega_m$")
plotUpDownLambdaRePH( pl.subplot(2,3,5) )
pl.xlabel(r"$\omega_n$")
plotUpDownLambdaRePHX( pl.subplot(2,3,6) )
pl.xlabel(r"$\omega_n$")
pl.tight_layout()

shift=6

#--- Save to file
pl.savefig("plots/Lambda.png", dpi = 150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()

#--- plot
pl.suptitle(r"$U=$" + str(UINT) + r"     $\beta=$" + str(BETA) + r"     $\epsilon=$" + str(EPS) +  r"     $\Omega_{\rm PP}=\Omega_{\rm PH}=\Omega_{\rm xPH}=$" + str(shift) + r"$*2\pi/\beta$" + r"     Notation: $\Lambda(\omega_1,\omega_2,\omega_1')$")

plotUpUpLambdaRePP( pl.subplot(2,3,1) )
pl.ylabel(r"$\omega_m$")
plotUpUpLambdaRePH( pl.subplot(2,3,2) )
plotUpUpLambdaRePHX( pl.subplot(2,3,3) )
plotUpDownLambdaRePP( pl.subplot(2,3,4) )
pl.xlabel(r"$\omega_n$")
pl.ylabel(r"$\omega_m$")
plotUpDownLambdaRePH( pl.subplot(2,3,5) )
pl.xlabel(r"$\omega_n$")
plotUpDownLambdaRePHX( pl.subplot(2,3,6) )
pl.xlabel(r"$\omega_n$")
pl.tight_layout()

#--- Save to file
pl.savefig("plots/Lambda_shift.png", dpi = 150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()

shift=0

#--------------------------------------GENCHI PLOTTING ------------------------------------------

print("Plotting genchi ...")

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

#---  Helper functions (include translation from  nambu to physical CAUTION : USING TIME REVERSAL SYMM) 
def plotgenchi( use_pl, arr, string ):
    use_pl.set_aspect(1.0)
    zarr = np.array([[ arr[shift + (bdim-1)/2,n,m,0,0,0,0,0,0,0] for n in range(fdim)] for m in range(fdim)])
    pl.pcolormesh( genchigrid, genchigrid, zarr )
    pl.ylim([min(genchigrid),max(genchigrid)])
    pl.xlim([min(genchigrid),max(genchigrid)])
    use_pl.set_title( string , fontsize=10 )
    pl.colorbar(shrink=0.6) 
    return

#--- Plot 
pl.suptitle(r"$U=$" + str(UINT) + r"     $\beta=$" + str(BETA) + r"     $\epsilon=$" + str(EPS) +  r"     $\Omega=$" + str(shift) + r"$*2\pi/\beta$" + r"     Notation: $F(\Omega,\omega,\omega')$")

plotgenchi( pl.subplot(2,3,1), regenchi_pp - regenchi_pp[:,:,::-1,:,:,:,:,:,:,:], RE + r"\chi^{PP}_{\uparrow\uparrow}(\Omega_{PP},\omega,\omega')$" ) # flip sign of w_out
pl.ylabel(r"$\omega'$")
pl.xlabel(r"$\omega$")
plotgenchi( pl.subplot(2,3,2), regenchi_ph - regenchi_xph, RE + r"\chi^{PH}_{\uparrow\uparrow}(\Omega_{PH},\omega,\omega')$" )
pl.xlabel(r"$\omega$")
plotgenchi( pl.subplot(2,3,3), regenchi_xph - regenchi_ph, RE + r"\chi^{XPH}_{\uparrow\uparrow}(\Omega_{XPH},\omega,\omega')$" )
pl.xlabel(r"$\omega$")

plotgenchi( pl.subplot(2,3,4), regenchi_pp, RE + r"\chi^{PP}_{\uparrow\downarrow}(\Omega_{PP},\omega,\omega')$" )
pl.ylabel(r"$\omega'$")
pl.xlabel(r"$\omega$")
plotgenchi( pl.subplot(2,3,5), regenchi_ph, RE + r"\chi^{PH}_{\uparrow\downarrow}(\Omega_{PH},\omega,\omega')$" )
pl.xlabel(r"$\omega$")
plotgenchi( pl.subplot(2,3,6), regenchi_xph, RE + r"\chi^{XPH}_{\uparrow\downarrow}(\Omega_{XPH},\omega,\omega')$" )
pl.xlabel(r"$\omega$")

pl.tight_layout()

#--- Save to file
pl.savefig("plots/genchi.png", dpi = 150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()

shift=6

#--- Plot
pl.suptitle(r"$U=$" + str(UINT) + r"     $\beta=$" + str(BETA) + r"     $\epsilon=$" + str(EPS) +  r"     $\Omega=$" + str(shift) + r"$*2\pi/\beta$" + r"     Notation: $F(\Omega,\omega,\omega')$")

plotgenchi( pl.subplot(2,3,1), regenchi_pp - regenchi_pp[:,:,::-1,:,:,:,:,:,:,:], RE + r"F^{PP}_{\uparrow\uparrow}(\Omega_{PP},\omega,\omega')$" ) # flip sign of w_out
pl.ylabel(r"$\omega'$")
pl.xlabel(r"$\omega$")
plotgenchi( pl.subplot(2,3,2), regenchi_ph - regenchi_xph, RE + r"\chi^{PH}_{\uparrow\uparrow}(\Omega_{PH},\omega,\omega')$" )
pl.xlabel(r"$\omega$")
plotgenchi( pl.subplot(2,3,3), regenchi_xph - regenchi_ph, RE + r"\chi^{XPH}_{\uparrow\uparrow}(\Omega_{XPH},\omega,\omega')$" )
pl.xlabel(r"$\omega$")

plotgenchi( pl.subplot(2,3,4), regenchi_pp, RE + r"\chi^{PP}_{\uparrow\downarrow}(\Omega_{PP},\omega,\omega')$" )
pl.ylabel(r"$\omega'$")
pl.xlabel(r"$\omega$")
plotgenchi( pl.subplot(2,3,5), regenchi_ph, RE + r"\chi^{PH}_{\uparrow\downarrow}(\Omega_{PH},\omega,\omega')$" )
pl.xlabel(r"$\omega$")
plotgenchi( pl.subplot(2,3,6), regenchi_xph, RE + r"\chi^{XPH}_{\uparrow\downarrow}(\Omega_{XPH},\omega,\omega')$" )
pl.xlabel(r"$\omega$")

pl.tight_layout()

#--- Save to file
pl.savefig("plots/genchi_shift.png", dpi = 150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()

shift=0

#--------------------------------------VERT PLOTTING ------------------------------------------

print("Plotting vert ...")

#--- Read
revert_pp = vert_mul * np.array(f["/vert_func/RE_PP"])
imvert_pp = vert_mul * np.array(f["/vert_func/IM_PP"])
revert_ph = vert_mul * np.array(f["/vert_func/RE_PH"])
imvert_ph = vert_mul * np.array(f["/vert_func/IM_PH"])
revert_xph = vert_mul * np.array(f["/vert_func/RE_XPH"])
imvert_xph = vert_mul * np.array(f["/vert_func/IM_XPH"])

bdim = revert_pp.shape[0]
fdim = revert_pp.shape[1]

vertgrid = np.array(f["/vert_func/fgrid"])

#---  Helper functions (include translation from  nambu to physical CAUTION : USING TIME REVERSAL SYMM) 
def plotvert( use_pl, arr, string ):
    use_pl.set_aspect(1.0)
    zarr = np.array([[ arr[shift + (bdim-1)/2,n,m,0,0,0,0,0,0,0] for n in range(fdim)] for m in range(fdim)])
    pl.pcolormesh( vertgrid, vertgrid, zarr )
    pl.ylim([min(vertgrid),max(vertgrid)])
    pl.xlim([min(vertgrid),max(vertgrid)])
    use_pl.set_title( string , fontsize=10 )
    pl.colorbar(shrink=0.6) 
    return

#--- Plot 
pl.suptitle(r"$U=$" + str(UINT) + r"     $\beta=$" + str(BETA) + r"     $\epsilon=$" + str(EPS) +  r"     $\Omega=$" + str(shift) + r"$*2\pi/\beta$" + r"     Notation: $F(\Omega,\omega,\omega')$")

plotvert( pl.subplot(2,3,1), revert_pp - revert_pp[:,:,::-1,:,:,:,:,:,:,:], RE + r"F^{PP}_{\uparrow\uparrow}(\Omega_{PP},\omega,\omega')$" ) # flip sign of w_out
pl.ylabel(r"$\omega'$")
pl.xlabel(r"$\omega$")
plotvert( pl.subplot(2,3,2), revert_ph - revert_xph, RE + r"F^{PH}_{\uparrow\uparrow}(\Omega_{PH},\omega,\omega')$" )
pl.xlabel(r"$\omega$")
plotvert( pl.subplot(2,3,3), revert_xph - revert_ph, RE + r"F^{XPH}_{\uparrow\uparrow}(\Omega_{XPH},\omega,\omega')$" )
pl.xlabel(r"$\omega$")

plotvert( pl.subplot(2,3,4), revert_pp, RE + r"F^{PP}_{\uparrow\downarrow}(\Omega_{PP},\omega,\omega')$" )
pl.ylabel(r"$\omega'$")
pl.xlabel(r"$\omega$")
plotvert( pl.subplot(2,3,5), revert_ph, RE + r"F^{PH}_{\uparrow\downarrow}(\Omega_{PH},\omega,\omega')$" )
pl.xlabel(r"$\omega$")
plotvert( pl.subplot(2,3,6), revert_xph, RE + r"F^{XPH}_{\uparrow\downarrow}(\Omega_{XPH},\omega,\omega')$" )
pl.xlabel(r"$\omega$")

pl.tight_layout()

#--- Save to file
pl.savefig("plots/vert.png", dpi = 150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()

shift=6

#--- Plot
pl.suptitle(r"$U=$" + str(UINT) + r"     $\beta=$" + str(BETA) + r"     $\epsilon=$" + str(EPS) +  r"     $\Omega=$" + str(shift) + r"$*2\pi/\beta$" + r"     Notation: $F(\Omega,\omega,\omega')$")

plotvert( pl.subplot(2,3,1), revert_pp - revert_pp[:,:,::-1,:,:,:,:,:,:,:], RE + r"F^{PP}_{\uparrow\uparrow}(\Omega_{PP},\omega,\omega')$" ) # flip sign of w_out
pl.ylabel(r"$\omega'$")
pl.xlabel(r"$\omega$")
plotvert( pl.subplot(2,3,2), revert_ph - revert_xph, RE + r"F^{PH}_{\uparrow\uparrow}(\Omega_{PH},\omega,\omega')$" )
pl.xlabel(r"$\omega$")
plotvert( pl.subplot(2,3,3), revert_xph - revert_ph, RE + r"F^{XPH}_{\uparrow\uparrow}(\Omega_{XPH},\omega,\omega')$" )
pl.xlabel(r"$\omega$")

plotvert( pl.subplot(2,3,4), revert_pp, RE + r"F^{PP}_{\uparrow\downarrow}(\Omega_{PP},\omega,\omega')$" )
pl.ylabel(r"$\omega'$")
pl.xlabel(r"$\omega$")
plotvert( pl.subplot(2,3,5), revert_ph, RE + r"F^{PH}_{\uparrow\downarrow}(\Omega_{PH},\omega,\omega')$" )
pl.xlabel(r"$\omega$")
plotvert( pl.subplot(2,3,6), revert_xph, RE + r"F^{XPH}_{\uparrow\downarrow}(\Omega_{XPH},\omega,\omega')$" )
pl.xlabel(r"$\omega$")

pl.tight_layout()

#--- Save to file
pl.savefig("plots/vert_shift.png", dpi = 150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()

shift=0

#--------------------------------------PHI PLOTTING ------------------------------------------

print("Plotting phi ...")

#--- Read
rephi_pp = vert_mul * np.array(f["/phi_func/RE_PP"])
imphi_pp = vert_mul * np.array(f["/phi_func/IM_PP"])
rephi_ph = vert_mul * np.array(f["/phi_func/RE_PH"])
imphi_ph = vert_mul * np.array(f["/phi_func/IM_PH"])
rephi_xph = vert_mul * np.array(f["/phi_func/RE_XPH"])
imphi_xph = vert_mul * np.array(f["/phi_func/IM_XPH"])

bdim = rephi_pp.shape[0]
fdim = rephi_pp.shape[1]
fdimo2 = fdim/2
phigrid = np.array(f["/phi_func/fgrid"])
phigrid_plot = np.array([(2*n+1)*pi/BETA for n in range(-fdimo2+1,fdimo2-1)])
#rephi_pp_upup = np.array([[[[[[[[[[ rephi_pp[W,m,n,K,k,kp,s1,s2,s3,s3] - rephi_pp[W,m,-n-1+fdim-mymod_abs(W),K,k,kp,s1,s2,s3,s4] for s4 in range(0,1)] for s3 in range(0,1)]for s2 in range(0,1)]for s1 in range(0,1)] for kp in range(0,1)] for k in range(0,1)] for K in range(0,1)] for n in range(fdim)] for m in range(fdim)] for W in range(bdim+1)])

def rephi_pp_updo( W, w, wp, K, k, kp , s1, s2, s3, s4 ):
    if ( not check_bounds_mix( W, w, wp ) ):
            return 0.0
    return rephi_pp[W,w,wp,K,k,kp,s1,s2,s3,s4]

#---  Helper functions (include translation from  nambu to physical CAUTION : USING TIME REVERSAL SYMM) 
def plotphi( use_pl, arr, string ):
    use_pl.set_aspect(1.0)
    zarr = np.array([[ arr[shift + (bdim-1)/2,n,m,0,0,0,0,0,0,0] for n in range(fdim)] for m in range(fdim)])
    pl.pcolormesh( phigrid, phigrid, zarr )
    pl.ylim([min(phigrid),max(phigrid)])
    pl.xlim([min(phigrid),max(phigrid)])
    use_pl.set_title( string , fontsize=10 )
    pl.colorbar(shrink=0.6) 
    return

def plotphiupup( use_pl, string ):
    use_pl.set_aspect(1.0)
    zarr = np.array([[ rephi_pp_updo(shift + (bdim-1)/2,n+fdim/2,m+fdim/2,0,0,0,0,0,0,0)-rephi_pp_updo(shift + (bdim-1)/2,n+fdim/2,-m-1-mymod_abs(shift + (bdim-1)/2)+fdim/2,0,0,0,0,0,0,0) for n in range(-fdimo2,fdimo2)] for m in range(-fdimo2,fdimo2)])
    pl.pcolormesh( phigrid, phigrid, zarr )
    pl.ylim([min(phigrid),max(phigrid)])
    pl.xlim([min(phigrid),max(phigrid)])
    use_pl.set_title( string , fontsize=10 )
    pl.colorbar(shrink=0.6) 
    return
#--- Plot 
pl.suptitle(r"$U=$" + str(UINT) + r"     $\beta=$" + str(BETA) + r"     $\epsilon=$" + str(EPS) +  r"     $\Omega=$" + str(shift) + r"$*2\pi/\beta$" + r"     Notation: $\phi(\Omega,\omega,\omega')$")

plotphiupup( pl.subplot(2,3,1), RE + r"\phi^{PP}_{\uparrow\uparrow}(\Omega_{PP},\omega,\omega')$" ) # flip sign of w_out
pl.ylabel(r"$\omega'$")
pl.xlabel(r"$\omega$")
plotphi( pl.subplot(2,3,2), rephi_ph - rephi_xph, RE + r"\phi^{PH}_{\uparrow\uparrow}(\Omega_{PH},\omega,\omega')$" )
pl.xlabel(r"$\omega$")
plotphi( pl.subplot(2,3,3), rephi_xph - rephi_ph, RE + r"\phi^{XPH}_{\uparrow\uparrow}(\Omega_{XPH},\omega,\omega')$" )
pl.xlabel(r"$\omega$")

plotphi( pl.subplot(2,3,4), rephi_pp, RE + r"\phi^{PP}_{\uparrow\downarrow}(\Omega_{PP},\omega,\omega')$" )
pl.ylabel(r"$\omega'$")
pl.xlabel(r"$\omega$")
plotphi( pl.subplot(2,3,5), rephi_ph, RE + r"\phi^{PH}_{\uparrow\downarrow}(\Omega_{PH},\omega,\omega')$" )
pl.xlabel(r"$\omega$")
plotphi( pl.subplot(2,3,6), rephi_xph, RE + r"\phi^{XPH}_{\uparrow\downarrow}(\Omega_{XPH},\omega,\omega')$" )
pl.xlabel(r"$\omega$")

pl.tight_layout()

#--- Save to file
pl.savefig("plots/phi.png", dpi = 150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()

shift=6

#--- Plot
pl.suptitle(r"$U=$" + str(UINT) + r"     $\beta=$" + str(BETA) + r"     $\epsilon=$" + str(EPS) +  r"     $\Omega=$" + str(shift) + r"$*2\pi/\beta$" + r"     Notation: $\phi(\Omega,\omega,\omega')$")

plotphiupup( pl.subplot(2,3,1),  RE + r"\phi^{PP}_{\uparrow\uparrow}(\Omega_{PP},\omega,\omega')$" ) # flip sign of w_out
pl.ylabel(r"$\omega'$")
pl.xlabel(r"$\omega$")
plotphi( pl.subplot(2,3,2), rephi_ph - rephi_xph, RE + r"\phi^{PH}_{\uparrow\uparrow}(\Omega_{PH},\omega,\omega')$" )
pl.xlabel(r"$\omega$")
plotphi( pl.subplot(2,3,3), rephi_xph - rephi_ph, RE + r"\phi^{XPH}_{\uparrow\uparrow}(\Omega_{XPH},\omega,\omega')$" )
pl.xlabel(r"$\omega$")

plotphi( pl.subplot(2,3,4), rephi_pp, RE + r"\phi^{PP}_{\uparrow\downarrow}(\Omega_{PP},\omega,\omega')$" )
pl.ylabel(r"$\omega'$")
pl.xlabel(r"$\omega$")
plotphi( pl.subplot(2,3,5), rephi_ph, RE + r"\phi^{PH}_{\uparrow\downarrow}(\Omega_{PH},\omega,\omega')$" )
pl.xlabel(r"$\omega$")
plotphi( pl.subplot(2,3,6), rephi_xph, RE + r"\phi^{XPH}_{\uparrow\downarrow}(\Omega_{XPH},\omega,\omega')$" )
pl.xlabel(r"$\omega$")

pl.tight_layout()

#--- Save to file
pl.savefig("plots/phi_shift.png", dpi = 150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()

shift=0

#--------------------------------------Chi PLOTTING ------------------------------------------

print("Plotting chi functions ...")


#--- Read
bosgrid = np.array(f["/chi_func/bgrid"])
rechi_pp = vert_mul * np.array(f["/chi_func/RE_PP"])
imchi_pp = vert_mul * np.array(f["/chi_func/IM_PP"])
rechi_ph = vert_mul * np.array(f["/chi_func/RE_PH"])
imchi_ph = vert_mul * np.array(f["/chi_func/IM_PH"])
rechi_xph = vert_mul * np.array(f["/chi_func/RE_XPH"])
imchi_xph = vert_mul * np.array(f["/chi_func/IM_XPH"])
fdim_bos = bosgrid.shape[0]

x_range_fact = 0.3

#--- Helper functions
def plotchi( use_pl, arr, title ):
    pl.plot(bosgrid, arr[:,0,0,0,0,0], 'b-', ms=3, mew=0.2)
    pl.xlim([x_range_fact*min(bosgrid),x_range_fact*max(bosgrid)])
    use_pl.set_title(title)
    return

pl.suptitle(r"$U=$" + str(UINT) + r"     $\beta=$" + str(BETA) + r"     $\epsilon=$" + str(EPS) + "  BEPs ")

plotchi( pl.subplot(2,3,1), rechi_pp - rechi_pp, RE + r"\chi^{PP}_{\uparrow\uparrow}$" )
plotchi( pl.subplot(2,3,2), rechi_ph - rechi_xph, RE + r"\chi^{PH}_{\uparrow\uparrow}$" )
plotchi( pl.subplot(2,3,3), rechi_xph - rechi_ph, RE + r"\chi^{XPH}_{\uparrow\uparrow}$" )

plotchi( pl.subplot(2,3,4), rechi_pp, RE + r"\chi^{PP}_{\uparrow\downarrow}$" )
pl.xlabel(r"$\Omega$")
plotchi( pl.subplot(2,3,5), rechi_ph, RE + r"\chi^{PH}_{\uparrow\downarrow}$" )
pl.xlabel(r"$\Omega$")
plotchi( pl.subplot(2,3,6), rechi_xph, RE + r"\chi^{XPH}_{\uparrow\downarrow}$" )
pl.xlabel(r"$\Omega$")

pl.tight_layout()

pl.savefig("plots/chi.png", dpi=150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()

#--------------------------------------P PLOTTING ------------------------------------------

print("Plotting P ...")

#--- Read
reP_pp = vert_mul * np.array(f["/P_func/RE_PP"])
imP_pp = vert_mul * np.array(f["/P_func/IM_PP"])
reP_ph = vert_mul * np.array(f["/P_func/RE_PH"])
imP_ph = vert_mul * np.array(f["/P_func/IM_PH"])
reP_xph = vert_mul * np.array(f["/P_func/RE_XPH"])
imP_xph = vert_mul * np.array(f["/P_func/IM_XPH"])

bdim = reP_pp.shape[0]
fdim = reP_pp.shape[1]

Pgrid = np.array(f["/P_func/fgrid"])

x_range_fact = 1.0

#--- Helper functions
def plotP_Omega( W, arr, color ):
    pl.plot(Pgrid + (W % 2) * math.pi/BETA, arr[W + (bdim-1)/2,:,0,0,0,0,0,0], color, ms=3, mew=0.2, label=r"$\Omega=" + str(W*2) + r"\pi/\beta$")

def plotP( use_pl, arr, title, legend ):
    plotP_Omega( 0, arr, 'r-' ) 
    plotP_Omega( 1, arr, 'g-' ) 
    plotP_Omega( 2, arr, 'b-' ) 
    plotP_Omega( 8, arr, 'c-' ) 
    plotP_Omega( 16, arr, 'm-' ) 
    #plotP_Omega( 32, arr, 'y-' ) 
    pl.xlim([x_range_fact*min(Pgrid),x_range_fact*max(Pgrid)])
    if ( legend ):
        pl.legend(prop={'size':7})

    use_pl.set_title(title)
    return

pl.suptitle(r"$U=$" + str(UINT) + r"     $\beta=$" + str(BETA) + r"     $\epsilon=$" + str(EPS) + "  BEPs ")

plotP( pl.subplot(2,3,1), reP_pp - reP_pp, RE + r"P^{PP}_{\uparrow\uparrow}$", True )
plotP( pl.subplot(2,3,2), reP_ph - reP_xph, RE + r"P^{PH}_{\uparrow\uparrow}$", False )
plotP( pl.subplot(2,3,3), reP_xph - reP_ph, RE + r"P^{XPH}_{\uparrow\uparrow}$", False )

plotP( pl.subplot(2,3,4), reP_pp, RE + r"P^{PP}_{\uparrow\downarrow}$", False )
pl.xlabel(r"$\omega$")
plotP( pl.subplot(2,3,5), reP_ph, RE + r"P^{PH}_{\uparrow\downarrow}$", False )
pl.xlabel(r"$\omega$")
plotP( pl.subplot(2,3,6), reP_xph, RE + r"P^{XPH}_{\uparrow\downarrow}$", False )
pl.xlabel(r"$\omega$")

pl.tight_layout()

pl.savefig("plots/P.png", dpi=150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()

#--------------------------------------R PLOTTING ------------------------------------------

print("Plotting R ...")

#--- Read
reR_pp = vert_mul * np.array(f["/R_func/RE_PP"])
imR_pp = vert_mul * np.array(f["/R_func/IM_PP"])
reR_ph = vert_mul * np.array(f["/R_func/RE_PH"])
imR_ph = vert_mul * np.array(f["/R_func/IM_PH"])
reR_xph = vert_mul * np.array(f["/R_func/RE_XPH"])
imR_xph = vert_mul * np.array(f["/R_func/IM_XPH"])

bdim = reR_pp.shape[0]
fdim = reR_pp.shape[1]

Rgrid = np.array(f["/R_func/fgrid"])

#---  Helper functions 
def plotR( use_pl, arr, string ):
    use_pl.set_aspect(1.0)
    zarr = np.array([[ arr[shift + (bdim-1)/2,n,m,0,0,0,0,0,0,0] for n in range(fdim)] for m in range(fdim)])
    pl.pcolormesh( Rgrid, Rgrid, zarr )
    pl.ylim([min(Rgrid),max(Rgrid)])
    pl.xlim([min(Rgrid),max(Rgrid)])
    use_pl.set_title( string , fontsize=10 )
    pl.colorbar(shrink=0.6) 
    return


#--- Plot 
pl.suptitle(r"$U=$" + str(UINT) + r"     $\beta=$" + str(BETA) + r"     $\epsilon=$" + str(EPS) +  r"     $\Omega=$" + str(shift) + r"$*2\pi/\beta$" + r"     Notation: $R(\Omega,\omega,\omega')$")

plotR( pl.subplot(2,3,1), reR_pp - reR_pp[:,:,::-1,:,:,:,:,:,:,:], RE + r"R^{PP}_{\uparrow\uparrow}(\Omega_{PP},\omega,\omega')$" ) # flip sign of w_out
pl.ylabel(r"$\omega'$")
pl.xlabel(r"$\omega$")
plotR( pl.subplot(2,3,2), reR_ph - reR_xph, RE + r"R^{PH}_{\uparrow\uparrow}(\Omega_{PH},\omega,\omega')$" )
pl.xlabel(r"$\omega$")
plotR( pl.subplot(2,3,3), reR_xph - reR_ph, RE + r"R^{XPH}_{\uparrow\uparrow}(\Omega_{XPH},\omega,\omega')$" )
pl.xlabel(r"$\omega$")

plotR( pl.subplot(2,3,4), reR_pp, RE + r"R^{PP}_{\uparrow\downarrow}(\Omega_{PP},\omega,\omega')$" )
pl.ylabel(r"$\omega'$")
pl.xlabel(r"$\omega$")
plotR( pl.subplot(2,3,5), reR_ph, RE + r"R^{PH}_{\uparrow\downarrow}(\Omega_{PH},\omega,\omega')$" )
pl.xlabel(r"$\omega$")
plotR( pl.subplot(2,3,6), reR_xph, RE + r"R^{XPH}_{\uparrow\downarrow}(\Omega_{XPH},\omega,\omega')$" )
pl.xlabel(r"$\omega$")

pl.tight_layout()

#--- Save to file
pl.savefig("plots/R.png", dpi = 150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()

shift=6

#--- Plot
pl.suptitle(r"$U=$" + str(UINT) + r"     $\beta=$" + str(BETA) + r"     $\epsilon=$" + str(EPS) +  r"     $\Omega=$" + str(shift) + r"$*2\pi/\beta$" + r"     Notation: $R(\Omega,\omega,\omega')$")


plotR( pl.subplot(2,3,1), reR_pp - reR_pp[:,:,::-1,:,:,:,:,:,:,:], RE + r"R^{PP}_{\uparrow\uparrow}(\Omega_{PP},\omega,\omega')$" ) # flip sign of w_out
pl.ylabel(r"$\omega'$")
pl.xlabel(r"$\omega$")
plotR( pl.subplot(2,3,2), reR_ph - reR_xph, RE + r"R^{PH}_{\uparrow\uparrow}(\Omega_{PH},\omega,\omega')$" )
pl.xlabel(r"$\omega$")
plotR( pl.subplot(2,3,3), reR_xph - reR_ph, RE + r"R^{XPH}_{\uparrow\uparrow}(\Omega_{XPH},\omega,\omega')$" )
pl.xlabel(r"$\omega$")

plotR( pl.subplot(2,3,4), reR_pp, RE + r"R^{PP}_{\uparrow\downarrow}(\Omega_{PP},\omega,\omega')$" )
pl.ylabel(r"$\omega'$")
pl.xlabel(r"$\omega$")
plotR( pl.subplot(2,3,5), reR_ph, RE + r"R^{PH}_{\uparrow\downarrow}(\Omega_{PH},\omega,\omega')$" )
pl.xlabel(r"$\omega$")
plotR( pl.subplot(2,3,6), reR_xph, RE + r"R^{XPH}_{\uparrow\downarrow}(\Omega_{XPH},\omega,\omega')$" )
pl.xlabel(r"$\omega$")

pl.tight_layout()

#--- Save to file
pl.savefig("plots/R_shift.png", dpi = 150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()

#---- Open files
#os.system('feh plots/phi.png')
#os.system('feh plots/Giw.png')
#os.system('feh plots/karr_func.png')
#os.system('feh plots/flow_obs.png')
#os.system('feh plots/chi.png')
#os.system('feh plots/trileg.png')
#os.system('feh plots/Sig.png')
