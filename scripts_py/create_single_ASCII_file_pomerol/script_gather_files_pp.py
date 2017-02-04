#!/usr/bin/python

#-----------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------
#                                   CONTENTS: 
#    
#   Merge 2P Green's Function obtained from Pomerol for different bosonic frequency transferred
#   in a single data file with the following form:
#   bos_freq  first_ferm_freq second_ferm_freq Re_chi1111 Im_chi1111 Re_chi1010 Im_chi1010
#
#------------------------------------------------------------------------------------------
#--------------------------------------IMPORTS ------------------------------------------

import h5py
import matplotlib.pyplot as pl
import numpy as np
import os
import sys
import subprocess
import matplotlib.cm as cm
from numpy.linalg import inv
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from os.path import isfile, join
import math

#--------------------------------------MANAGING FILES ------------------------------------------

def run(command):
    output = subprocess.check_output(command, shell=True)
    return output
pos_bfreq = 32
num_bfreq = 2*(pos_bfreq)-1
beta = 10.0
pi = math.pi

def twopgfupup(i):
    arr = np.loadtxt("../pomerol_output/datadirname/chi_pp_shift_1111_W%s.dat"%i)
    return arr

print twopgfupup(0)[:,0]

def twopgfupdo(i):
    arr1 = np.loadtxt("../pomerol_output/datadirname/chi_pp_shift_1001_W%s.dat"%i)
    return arr1

def filename(i):
    return '../pomerol_output/datadirname/2pgf_pp_shift_W%s.dat'%i

index_list = np.array([i for i in range(0,num_bfreq)])
print index_list

for index in index_list:
    np.savetxt('../pomerol_output/datadirname/2pgf_pp_shift_W%s.dat'%index, np.column_stack((twopgfupup(index)[:,0], twopgfupup(index)[:,1], twopgfupup(index)[:,2], twopgfupup(index)[:,3], twopgfupup(index)[:,4],twopgfupdo(index)[:,3], twopgfupdo(index)[:,4])))

filenames=[filename(i) for i in range(0,num_bfreq)] 

with open('../pomerol_output/datadirname/2pgf_pp_shift.dat', 'w') as outfile:
     for fname in filenames:
         with open(fname) as infile:
              for line in infile:
                  outfile.write(line)                                                        
