#!/usr/bin/env python3

'''
Example1: Strebelle's TI along two directions
*************************************************************************

:file name: `s2Dcd_run-ex01.py`

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
    
.. remark::
    This is only a toy example, useful to demonstrate the capabilities of the
    `s2Dcd` approach, and it is not pretended to reproduce any geological
    "reality" with the current setting used here.

'''

# Import "standard" modules
import os
import numpy as np
import geone as gn
import json
import copy

# Import non-standard modules
import s2Dcd.s2Dcd as s2Dcd
import s2Dcd.utili as utili

# Print header and start counting time
time_start = utili.print_start() 

# If you can use a parallel version of the MPS simulation core,
# define here the number of threads.
nthreads = 8

# %% Read some info about the TI

# Read from an external file the dimensions that would have a 3D TI
# (when composed by the 2D slices considered here)
ti3Ddict_file = "ti3Ddict.json"
with open(ti3Ddict_file, "r") as json_in:
    ti3Ddict = json.load(json_in)

# Use the parameters to create an empty 3D TI
ti3Ddict["val"] = np.full([ti3Ddict["nx"], ti3Ddict["ny"], ti3Ddict["nz"]],
                          np.nan)
ti3D = gn.img.Img(**ti3Ddict)

# %% Read the simulation grid and the main simulation parameters    
#

# Create a dictionary containing the parameters that will be used to set up
# the default input parameters of the DeeSse
ds3Ddict_file = "ds3Ddict.json"
with open(ds3Ddict_file, "r") as json_in:
    ds3Din_dict = json.load(json_in)
nx = ds3Din_dict["nx"]    
ny = ds3Din_dict["ny"]
nz = ds3Din_dict["nz"]
sx = ds3Din_dict["sx"]
sy = ds3Din_dict["sy"]
sz = ds3Din_dict["sz"]
ox = ds3Din_dict["ox"]
oy = ds3Din_dict["oy"]
oz = ds3Din_dict["oz"]
nv = ds3Din_dict["nv"]
varname = ds3Din_dict["varname"]
name = "res3D.gslib"
# Add a 3D TI to the parameter file
ds3Din_dict["TI"] = np.array([ti3D])

# %% Create an empty simulation grid to be filled with the sequential 2D simulations

# Create an empty simulation grid (Img) to be filled with the 2D simulations
val = np.full([ds3Din_dict["nx"], ds3Din_dict["ny"], ds3Din_dict["nz"]], np.nan)
res3D = gn.img.Img(nx, ny, nz, sx, sy, sz, ox, oy, oz, nv, val, varname, name)

# %% Set up the parameters and the training images

# Set up the general parameters to be used for the target 3D final result
ds3Din = gn.deesseinterface.DeesseInput(**ds3Din_dict)

# Set up the main parameters for the simulation along the plane *yz*
# from the 3D template (along this direction we do not have a TI in this example)
ds_yz_in = None

# Set up the main parameters for the simulation along the plane *xz*, starting
# from the 3D template
ds_xz_in = copy.deepcopy(ds3Din)
# For the moment we only set up the size of the SG. The origin will be defined
# within the simulation sequence.
ds_xz_in.ny = 1

# Read the training image along the plane xz
ti_xz = gn.img.readImageGslib(os.path.join("..","data","strebelle", 
                                           "ti_250x1x250.gslib"))
ds_xz_in.TI = np.array([ti_xz])

# Set up the main parameters for the simulation along the plane *xy*, starting
# from the 3D template
ds_xy_in = copy.deepcopy(ds3Din)
ds_xy_in.nz = 1
# Read the training image along the plane *xy*
ti_xy = gn.img.readImageGslib(os.path.join("..","data","strebelle", 
                                           "ti_250x250x1.gslib"))
ds_xy_in.TI = np.array([ti_xy])

# %% perform the simulation

# The max number of simulation steps to be performed. You can use it
# if you want to stop the simulation before the 3D domain is
# completed (i.e. for quick testing purposes...).
step_max = 1000000

# Print some simulation info
s2Dcd.print_sim_info(ds3Din, ds_yz_in, ds_xz_in, ds_xy_in, nthreads)

# Create the simulation sequence
seq = s2Dcd.create_seq(ds3Din, ds_yz_in, ds_xz_in, ds_xy_in, nthreads)

#
# Simulation
#
s2Dcd.sim_run(seq, step_max, res3D, ds3Din, nthreads)

# Stop counting time
utili.print_stop(time_start)


gn.img.writeImageVtk(res3D, "res3D.vtk", missing_value=-9999999)

# %% Print the result in 3D with PyVista
gn.imgplot3d.drawImage3D_surface(res3D, text='TI', scalar_bar_kwargs={'vertical':True})
