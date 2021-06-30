import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import cm
import matplotlib.colors as colors
import numpy.fft as fft
import corner
import h5py
import sys
import copy 
import matplotlib as mpl

import tools
import map_cosmo
import power_spectrum
from tflib import get_1D_TF, get_2D_TF, PS_plotter

from sklearn.decomposition import PCA

simpath = "some path"
mappath = "some path"
noisepath = "some path"


# Plot PS and TF:

outputplot = "some file name"

PS_plotter(simpath, mappath, noisepath, outname = outputplot)


# Generate HDF5 file with TF(k) in 1D:

TF, k = get_1D_TF(simpath, mappath, noisepath)

dfile = h5py.File("TF_1d.h5", "a")

dfile.create_dataset("TF", data = TF)
dfile.create_dataset("k", data = k)
dfile.close()

# Generate HDF5 file with TF(k) in 1D:

TF, k = get_2D_TF(simpath, mappath, noisepath)

dfile = h5py.File("TF_2d.h5", "a")

dfile.create_dataset("TF", data = TF)
dfile.create_dataset("k", data = k)
dfile.close()