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


pos_bfreq = 90
num_bfreq = 2*(pos_bfreq)+1
beta = 50.0
pi = math.pi


dataloc="/home/agnese/Coding/Vertex_structures/scripts_py/create_single_ASCII_file_EDtriamp/U_1.0_beta_50_hf/"

dirname = dataloc + "ALL_IN_ONE_PH/" 


if os.path.exists(dirname):
    os.system('rm -r ALL_IN_ONE_PH')

os.mkdir(dirname)
def trilegamp(i):
    arr1 = np.loadtxt("triamp%s"%i)
    return arr1

def filename_ph(i):
    return 'ALL_IN_ONE_PH/trileg_ph_W%s.dat'%i

index_list = np.array([i for i in range(0,num_bfreq)])

print index_list

for index in index_list:
    np.savetxt('ALL_IN_ONE_PH/trileg_ph_W%s.dat'%index, np.column_stack((2*(index-pos_bfreq)*pi/beta, trilegamp(index)[:,0], trilegamp(index)[:,1], trilegamp(index)[:,2],trilegamp(index)[:,3], trilegamp(index)[:,4]))

filenames=[filename_pp(i) for i in range(0,num_bfreq)] 

with open('ALL_IN_ONE_PH/trileg.dat', 'w') as outfile:
     for fname in filenames:
         with open(fname) as infile:
              for line in infile:
                  outfile.write(line)                                                        
os.system('rm ALL_IN_ONE_PH/trileg_ph_W*.dat')                                                        
