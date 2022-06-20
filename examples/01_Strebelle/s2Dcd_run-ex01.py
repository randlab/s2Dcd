#!/usr/bin/env python3

'''
Example1: Strebelle's TI along two directions
*************************************************************************

:file name: `s2Dcd_run-ex01_geone.py`

:directory: `/examples/deesse/01_Strebelle`.

In this example, the `s2Dcd` method is applied to the celebrated
Strebelle's 2D training image to obtain a 3D simulation.
Here we use the same training images along the directions normal to the axis
*z* and *y*.
   
No conditioning data are considered.

.. note::
    Another possibility could be, for example, to use the two TI normal to
    axis *x* and *y*. In this case the hypothesis about the 3D geometry are
    of course completely different.

'''

# Import "standard" modules
import os
import numpy as np
import geone as gn
import json
import copy

# Import non-standard modules
import s2Dcd.s2Dcd as s2Dcd
import s2Dcd.utili  as utili

# Print header and start counting time
time_start = utili.print_start() 

# If you can use a parallel version of the MPS simulation core,
# define here the number of threads.
nthreads = 8

# Defalut information that 
ti3Dini_file = "ti3Dini.json"

# Read an external file containing the dimensions of the output grid.
with open(ti3Dini_file, "r") as json_in:
    ti3Dini = json.load(json_in)
# The initial content is NaN
ti3Dini["val"] = np.full([ti3Dini["nx"], ti3Dini["ny"], ti3Dini["nz"]], np.nan)

# Use the parameters to create an empty 3D TI
res3D = gn.img.Img(**ti3Dini)

# Create a dictionary containing the parameters that will be used to set up
# the default input parameters of the DeeSse
ds3Din_file = "ds3Din.json"
with open(ds3Din_file, "r") as json_in:
    ds3Din_dict = json.load(json_in)
ds3Din_dict["TI"] = res3D

# Set up the general parameters to be used for the target 3D final result
ds3Din = gn.deesseinterface.DeesseInput(**ds3Din_dict)

# Set up the main parameters for the simulation normal to axis *x*, starting
# from the 3D template
dsXin = None

# Set up the main parameters for the simulation normal to axis *y*, starting
# from the 3D template
dsYin = copy.deepcopy(ds3Din)
dsYin.ny = 1

# Read the training image normal to the *y* axis
tiY = gn.img.readImageGslib(os.path.join("..","data","strebelle", 
                                           "ti_250x1x250.gslib"))
dsYin.TI = tiY

# Set up the main parameters for the simulation normal to axis *z*, starting
# from the 3D template
dsZin = copy.deepcopy(ds3Din)
dsZin.nz = 1
# Read the training image normal to the *z* axis
tiZ = gn.img.readImageGslib(os.path.join("..","data","strebelle", 
                                           "ti_250x250x1.gslib"))
dsZin.TI = tiZ

# The max number of simulation steps to be performed. You can use it
# if you want to stop the simulation before the 3D domain is
# completed.
step_max = 300000

# Print some simulation info
s2Dcd.print_sim_info(ds3Din, dsXin, dsYin, dsZin, nthreads)

# Create the simulation sequence
seq = s2Dcd.create_seq(ds3Din, dsXin, dsYin, dsZin, nthreads)

#
# Simulation
#
s2Dcd.sim_run(seq, step_max, res3D, ds3Din, nthreads)

# Stop counting time
utili.print_stop(time_start)

gn.imgplot3d.drawImage3D_surface(res3D, text='TI', scalar_bar_kwargs={'vertical':True})

gn.img.writeImageVtk(res3D, "res3D.vtk", missing_value=-9999999)
