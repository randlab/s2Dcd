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

:this file:
    `utili.py`

:Purpose:
    A collection of small utility functions

:Version:
    * 0.6 , 2017-12-15 :
        * Improved the way computing time is computed by using
          ``time.perf_counter()``
    * 0.5 , 2015-12-03 :
        * Some minor cleanup
    * 0.4 , 2013-03-06 :
        * Converted to Python3.X
    * 0.3 , 2013-03-06 :
        * Last version before the conversion to python3.X
        * Added a function to skip the C-like comments.
    * 0.2 , 2012-05-02 :
        * Added a function usefult for the output print...
    * 0.1 , 2012-03-09 :
        * First version.

:Authors:
    Alessandro Comunian

'''

import sys
import time
import inspect
import ntpath

AUTHOR = 'Alessandro Comunian'

# Define the length of the integers in "add_file_id"
NB_DGT = 6

# Number of significant figures for the "e" format
SIG_FIG = 6

def print_start(title=None):
    '''
    Print the "begin" header and start counting time...

    Parameters:
        title: string, optional
            A message to be printed in a "formatted stype". If a
            string is not provided, then use the name of the calling
            script.

    Returns:
        The output of the function `time.clock()`.

    '''

    if title:
        message = title
    else:
        message = ntpath.basename(inspect.stack()[1][1])

    decor = "    **********" + "*"*len(message)
        
    print("\n")
    print(decor)
    print("    *** ", message, " ***")
    print(decor)
    print("    START: ", time.ctime())

    return time.perf_counter()  # Start counting time


def print_stop(time_start=None):
    '''
    Print the "end" message and the computing time.

    Parameters:
        time_start: ouput of the function :py:func:`print_start`, optional
            The function should be called in conjunction with the
            function `print_start`. If a value for time_start is not
            provided, then simply print out the current time and the
            stop message.

    Returns:
        Print out the computing time
        
    .. note:
        True CPU time may differ. This way of computing performance
        gives only a rough estimate of the performances.
    '''

    print("\n    STOP:  ", time.ctime())
    if time_start: 
        print("    Ellapsed time [sec]: %.2f" % (time.perf_counter() - time_start))
    print("\n")
    

def dtype_fmt( x):
    '''
    Returns a printing format according to the type of the input numpy
    array.

    Parameters:
        x: numpy array

    Returns:
        A string that defines the output format that should be used for the
        formatted output, that is '%d' is the case of an integer, and
        '%.4e' in case of a float.
        In case of an error, returs -1.
    '''

    if x.dtype.name in ['int','int8','int32','int64']:
        dtypefmt = '%d'
    elif x.dtype.name in  ['float','float8','float32','float64']:
        dtypefmt = '%.{0}e'.format(SIG_FIG)
    else:
        print("    Error in `dtype_fmt`, unknown data type"
              " ({0})".format(x.dtype.name))
        return -1
    
    return dtypefmt


def skip_ccomments(par_file):
    '''
    Read the connection to a file and extract a list containing all
    the fields, excluding line commenten in a *C-like* style, that is
    starting with "//" or enclosed by "/*...*/".
    
    Parameters
        file: file object
            The file to be read.
    Returns
        A list containing all the strings (the parameters), separated
        by a space.  Of course, the structure of the list depends on
        the structure of the paramters file.

    .. note::
        For the moment only the space characters are considered as field
        separators. All the fields in the list are returned as strings,
        therefore you will have to "cast" them according to your needs.
    
    .. warning::
        The comments at the end of a line containing parameters are
        not deleted. Therefore, used should be aware about what they
        are reading on each line...
    '''

    multiline_comm = 0 # A flag to check if we are in a /*...*/ comment
    
    par_list = [] # This list will contains all the important lines
    
    for line in par_file:
        x=line.split() # Create a list with the content of the line
        if x: # (the string is not empty)
            # Skip all the comments
            if x[0] == '//' or line.startswith('//'):
                continue
            elif ('/*' in line) or multiline_comm==1:
                if '*/' in line:
                    multiline_comm = 0
                    continue
                else:
                    multiline_comm = 1
                    continue
            
        # Add the important lines to a list
        if x:
            par_list.append(x)

    return par_list



def savepdf(filename, title=''):
    """
    Save a matplotlib figure as PDF and set up some useful information as
    metadata in the PDF file.

    Parameters:
        filename : string
            The name of the PDF file.
        title : string, optional
            The title of the PDF document.

    Returns:
        A PDF file containing the figure with the selected metadata.
   
    """
    import sys
    import pylab as pl
    from matplotlib.backends.backend_pdf import PdfPages

    pdffig = PdfPages(filename)

    pl.savefig(pdffig, format="pdf")

    d = pdffig.infodict()
    d['Title'] = title
    d['Author'] = AUTHOR
    d['Subject'] = "Created by: "+"  ".join(sys.argv)
#    d['Keywords'] = 'PdfPages example'
    
    pdffig.close()


    

def add_file_id(file_name, file_id):
    """
    Add a file ID to a file name, given an int ID.

    Parameters:
        file_name: string
            The name of the file (with extension...)
        file_id: int
            The ID which is associated to the file name
    """
    
    import os
    
    (file_root, file_ext) = os.path.splitext(file_name)
    id_format = '{{0}}_{{1:0{0:d}d}}'.format(NB_DGT)
    file_root_new = id_format.format(file_root, file_id)
    
    return file_root_new + file_ext
    
    
def progress_bar(i, nb_steps):
    """
    A simple progress bar

    Parameters:
        i: integer
            The current index

        nb_steps:
            The total number of steps
    """
    if i%(nb_steps//20)==0:
        sys.stdout.write('\r')
        sys.stdout.write(\
            "    [%-20s] %3d%%" % ('='*(i*20//nb_steps), 5*i*20//nb_steps))
        sys.stdout.flush()

def progress_bar_end():
    """
    Print the end of the progress bar
    """
    sys.stdout.write('\r')
    sys.stdout.write("    [%-20s] %3d%%" % ('='*20, 100))
    sys.stdout.flush()
