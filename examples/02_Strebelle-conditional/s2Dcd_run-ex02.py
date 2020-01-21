#!/usr/bin/env python3
'''
Example 2: Conditioning data
************************************************************************

:file name: `s2Dcd_run-ex02.py`

:directory: `/examples/02_Strebelle-conditional`.

Same as `Example 1`, but with conditioning points.

'''

# Import "standard" modules
import os
import time
import sys
import numpy as np
import random
import copy

# Import the "s2Dcd" modules
import s2Dcd.s2Dcd as s2Dcd
import s2Dcd.deesse as mpds_interface
import s2Dcd.utili as utili
import s2Dcd.gslibnumpy as gslibnumpy

# If you can use a parallel version of the MPS simulation core,
# define here the number of threads.
# 
s2Dcd.nb_threads = 8

# Print header and start counting time
time_start = utili.print_start() 

# The random seed, useful to get repeatable results
seed = np.random.RandomState(456833)

# Load the default simulation parameters from a template file
template_in = "template.in"
par_template = mpds_interface.Param(file_name=template_in)

# Create with the same files some templates for the diffent directions.
# If there is no TI along a given directgion, use the keyword "None".
par_Xnorm = mpds_interface.Param(file_name=template_in)
par_Xnorm.tis[0].file_name = None
par_Ynorm = mpds_interface.Param(file_name=template_in)
par_Ynorm.tis[0].file_name = os.path.join("..","data","strebelle", 
                                           "ti_250x1x250.gslib")
par_Znorm = mpds_interface.Param(file_name=template_in)
par_Znorm.tis[0].file_name = os.path.join("..","data","strebelle", 
                                           "ti_250x250x1.gslib")

# The max number of simulation steps to be performed. You can use it
# if you want to stop the simulation before the 3D domain is
# completed.
step_max = 3000

# Instantiate an object that contains the dimensions of the simulation
# grid
simODS = par_template.grid

s2Dcd.print_sim_info(simODS, par_Xnorm, par_Ynorm, par_Znorm)

# Array to store hard data and the result of the simulations. By
# default, the cells that don't contain a simulated (or a default
# conditioning data) are set to the value of the variable
# `s2Dcd.no_data`.
hard_data = s2Dcd.no_data*np.ones(
    (simODS.nx, simODS.ny, simODS.nz), 'float')

# Add some conditioning data from a GSLIB point data file
s2Dcd.add_gslib_pointdata(["data_points.gslib"], hard_data, simODS)

# Load an additional conditioning section and add it to the hard data
cond_sect_file = os.path.join(
    "..","data", "strebelle", "ti_250x1x250.gslib")
hd_section = gslibnumpy.gslib2numpy_onevar(cond_sect_file)[0]

# Assign to the hard data only a sub-section of the loaded
# conditioning section.
hard_data[:,4,:] = hd_section[0:simODS.nx,0,0:simODS.nz]

# Create the simulation sequence
seq = s2Dcd.create_seq(simODS, par_Xnorm, par_Ynorm, par_Znorm)
    
#
# Simulation
#
s2Dcd.sim_run(seq, step_max, hard_data, simODS, par_template, seed)

# Stop counting time
utili.print_stop(time_start)


        

    
    
 
