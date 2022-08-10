'''
:license:
    This file is part of s2Dcd.

    s2Dcd is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
 
    s2Dcd is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with s2Dcd.  If not, see <https://www.gnu.org/licenses/>.

:this_file: `ext.py`

:Purpose:
    A library containing the definition for the extensions for a number
    of formats...

:Version:

* 0.0.2, 2013-02-26 :
    Added some formats and improved the readibility.

* 0.0.1, 2012-05-15 :
    The first version.


'''


# --------------------------------------------------------------------------
# Generic formats
# --------------------------------------------------------------------------
TXT = ".txt"
LOG = ".log"
EXE = ".exe"
PDF = ".pdf"

# --------------------------------------------------------------------------
# ParaView Specific format
# --------------------------------------------------------------------------
VTK = ".vtk"
VTM = ".vtm"


# --------------------------------------------------------------------------
# Formats related to various geostatistics tools
# --------------------------------------------------------------------------
GSLIB = ".gslib"
SGEMS = ".sgems"
IMPALA_IN = ".in"
MPDS_IN = ".in"
# This is a little "generic", but...
DAT = ".dat"
# The output files of FLUVSIM have this extension
OUT = ".out"


# --------------------------------------------------------------------------
# Python related formats
# --------------------------------------------------------------------------
NPZ = ".npz"


# --------------------------------------------------------------------------
# Images related formats
# --------------------------------------------------------------------------
PNG = ".png"


# --------------------------------------------------------------------------
# Parflow related formats
# --------------------------------------------------------------------------
# Parflow binary files
PFB = ".pfb"




