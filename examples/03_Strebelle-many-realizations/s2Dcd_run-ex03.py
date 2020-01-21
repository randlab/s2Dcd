#!/usr/bin/env python3
'''
Example 3: Multiple realizations
***********************************************************************

:file name: `s2Dcd_run-ex03.py`

:directory: `/examples/03_Strebelle-many-realizations`.

In this example, the same simulation framework of `Example 1` is applied
many times in order to obtain different simulations with different 
random seeds.

'''

# Import "standard" modules
import os
import numpy as np

# Import non-standard modules
import s2Dcd.s2Dcd as s2Dcd
import s2Dcd.deesse as mpds_interface
import s2Dcd.utili as utili 

# Print header and start counting time
time_start = utili.print_start() 

# Load the default simulation parameters from a template file
template_in = "template.in"
par_template = mpds_interface.Param(file_name=template_in)

# Create with the same files some templates for the diffent directions.
# If there is no TI along a given direction, use the keyword "None".
par_Xnorm = mpds_interface.Param(file_name=template_in)
par_Xnorm.tis[0].file_name = None
par_Ynorm = mpds_interface.Param(file_name=template_in)
par_Ynorm.tis[0].file_name =  os.path.join("..","data","strebelle", 
                                           "ti_250x1x250.gslib")
par_Znorm = mpds_interface.Param(file_name=template_in)
par_Znorm.tis[0].file_name =  os.path.join("..","data","strebelle",
                                           "ti_250x250x1.gslib")

# The max number of simulation steps to be performed. You can use it if you
# want to stop the simulation before the 3D domain is completed.
step_max = 3000

# Instantiate an object that contains the dimensions of the simulation
# grid
simODS = par_template.grid

s2Dcd.print_sim_info(simODS, par_Xnorm, par_Ynorm, par_Znorm)

# Create the simulation sequence
seq = s2Dcd.create_seq(simODS, par_Xnorm, par_Ynorm, par_Znorm)

nb_realizations = 3

for i in range(nb_realizations):
    # The random seed, useful to get repeatable results
    seed = np.random.RandomState(np.random.randint(10e+6, 10e+8))
    
    # Array to store hard data and the result of the simulations. By
    # default, the cells that don't contain a simulated (or a default
    # conditioning data) are set to the value of the variable
    # `s2Dcd.no_data`.
    hard_data = s2Dcd.no_data*np.ones(
        (simODS.nx, simODS.ny, simODS.nz), 'float')

    print('\n    *** Run realization #', i, '***')
    # Simulation
    s2Dcd.sim_run(seq, step_max, hard_data, simODS, par_template, seed)

    # Rename the result and store it with a name related to the
    # realization number
    os.rename("result.gslib","result%05i.gslib"%(i))

# Stop counting time
utili.print_stop(time_start)


        

    
    
 
