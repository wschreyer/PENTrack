# Some standard import statements
# statements begining with !pip are ones you might need to run once
# just to install these packages.

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
ipython.magic("matplotlib inline")
plt.rcParams['figure.dpi'] = 200

from mpl_toolkits.mplot3d import Axes3D