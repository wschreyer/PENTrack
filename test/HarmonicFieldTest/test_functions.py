# Some standard import statements
# statements begining with !pip are ones you might need to run once
# just to install these packages.

import ROOT
import os
import shutil
from rootpy.io import root_open, DoesNotExist
# from rootpy.plotting import Hist, Hist2D
from rootpy import testdata
from rootpy import asrootpy

from rootpy.extern.six.moves import range
# from rootpy.plotting import Hist, Hist2D, Hist3D, HistStack, Legend, Canvas
# from rootpy.interactive import wait

import sys
import os
import logging
import time
import decimal
import datetime as dt
import numpy as np
from scipy.optimize import curve_fit
from collections import OrderedDict

# !pip install uncertainties
import uncertainties
from uncertainties import unumpy
from uncertainties import *
# !pip install lmfit
from lmfit import Model

from IPython import get_ipython
ipython = get_ipython()
# if you use your own separate scripts with function definitions
# these commands make your notebook grab updates from those script
# files every time you run a code cell in the notebook. saves time.
ipython.magic("load_ext autoreload")
ipython.magic("autoreload 2")
logging.basicConfig(level=logging.INFO, stream=sys.stdout)

# for plotting
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
ipython.magic("matplotlib inline")
plt.rcParams['figure.dpi'] = 200
# ipython.magic("matplotlib inline")
# plt.rcParams['figure.dpi'] = 200

from mpl_toolkits.mplot3d import Axes3D

###############################################################################
###############################################################################

def load_root(B_dict, filename):
    """A function that loads data from root files into a data dictionary
    
    Arguments:
        filename {string} -- name of root file to load
        B_dict {dict} -- dictionary containing data from the PENTrack BField
        Cut
    
    Returns:
        dict -- the same dictionary, with new data added, is returned
    """
    # this might need to be adjusted depending on path to the file
    # it remains convenient for dictionary access to keep the actual
    # filename string short
    f = root_open(filename)

    x     = np.empty((0,1), float)
    y     = np.empty((0,1), float)
    z     = np.empty((0,1), float)
    Bx    = np.empty((0,1), float)
    By    = np.empty((0,1), float)
    Bz    = np.empty((0,1), float)
    dBxdx = np.empty((0,1), float)
    dBydy = np.empty((0,1), float)
    dBzdz = np.empty((0,1), float)

    for evt in f.mytree:
        x  = np.append(x, evt.x)
        y  = np.append(y, evt.y)
        z  = np.append(z, evt.z)
        Bx = np.append(Bx, evt.Bx)
        By = np.append(By, evt.By)
        Bz = np.append(Bz, evt.Bz)
        dBxdx = np.append(dBxdx, evt.dBxdx)
        dBydy = np.append(dBydy, evt.dBydy)
        dBzdz = np.append(dBzdz, evt.dBzdz)
    
    B_dict[filename, 'x'] = x
    B_dict[filename, 'y'] = y
    B_dict[filename, 'z'] = z
    B_dict[filename, 'Bx'] = Bx
    B_dict[filename, 'By'] = By
    B_dict[filename, 'Bz'] = Bz
    B_dict[filename, 'dBxdx'] = dBxdx
    B_dict[filename, 'dBydy'] = dBydy
    B_dict[filename, 'dBzdz'] = dBzdz
    
    B_dict[filename, 'x', 'label'] = r'$x$'
    B_dict[filename, 'y', 'label'] = r'$y$'
    B_dict[filename, 'z', 'label'] = r'$z$'
    B_dict[filename, 'Bx', 'label'] = r'$B_x$'
    B_dict[filename, 'By', 'label'] = r'$B_y$'
    B_dict[filename, 'Bz', 'label'] = r'$B_z$'
    B_dict[filename, 'dBxdx', 'label'] = r'$\partial_xB_x$'
    B_dict[filename, 'dBydy', 'label'] = r'$\partial_yB_y$'
    B_dict[filename, 'dBzdz', 'label'] = r'$\partial_zB_z$'
    
    return B_dict
   
###############################################################################
###############################################################################

def plot_3D_surf(B_dict, filename, p0, p1, p2):
    
    a = B_dict[filename, p0]
    b = B_dict[filename, p1]
    c = B_dict[filename, p2]
    
    fig, ax = plt.subplots()
    ax = plt.axes(projection='3d')
    ax.scatter(a, b, c, c=c, cmap='viridis', linewidth=0.5);
    fig.tight_layout()
    
    ax.set_xlabel(B_dict[filename, p0, 'label'])
    ax.set_ylabel(B_dict[filename, p1, 'label'])
    ax.set_zlabel(B_dict[filename, p2, 'label']);
    
    
    