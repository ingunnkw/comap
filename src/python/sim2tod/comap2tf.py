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
from tflib import get_1D_TF, get_2D_TF, PS_plotter, plot_two_TFs_and_diff

from sklearn.decomposition import PCA

simpath_2   = "/mn/stornext/d22/cmbco/comap/protodir/maps/simulations/cube_maps/co2_cube_slow_az.h5"
mappath_2   = "/mn/stornext/d22/cmbco/comap/protodir/maps/simulations/co2_map_slow_az_sim.h5"
noisepath_2 = "/mn/stornext/d22/cmbco/comap/protodir/maps/co2_map_slow_az.h5"

# simpath = "/mn/stornext/d22/cmbco/comap/protodir/maps/simulations/cube_maps/co2_cube_slow_az_norm_0_003.h5"
# mappath = "/mn/stornext/d22/cmbco/comap/protodir/maps/simulations/co2_map_slow_az_sim_norm_0_003.h5"
# noisepath = "/mn/stornext/d22/cmbco/comap/protodir/maps/co2_map_slow_az_norm_0_003.h5"

# simpath_1   = "/mn/stornext/d22/cmbco/comap/protodir/maps/simulations/cube_maps/co2_cube_new_clock.h5"
# mappath_1   = "/mn/stornext/d22/cmbco/comap/protodir/maps/simulations/co2_map_new_clock_sim.h5"
# noisepath_1 = "/mn/stornext/d22/cmbco/comap/protodir/maps/co2_map_new_clock.h5"

# simpath = "/mn/stornext/d22/cmbco/comap/nils/COMAP_general/data/maps/test_tf/co2_cube_test_tf.h5"
# mappath = "/mn/stornext/d22/cmbco/comap/nils/COMAP_general/data/maps/test_tf/co2_test_tf_sim.h5"
# noisepath = "/mn/stornext/d22/cmbco/comap/nils/COMAP_general/data/maps/test_tf/co2_test_tf.h5"

simpath_1 = "/mn/stornext/d22/cmbco/comap/nils/COMAP_general/data/maps/test_tf/co2_cube_test_tf_selection.h5"
mappath_1 = "/mn/stornext/d22/cmbco/comap/nils/COMAP_general/data/maps/test_tf/co2_test_tf_selection_sim.h5"
noisepath_1 = "/mn/stornext/d22/cmbco/comap/nils/COMAP_general/data/maps/test_tf/co2_test_tf_selection.h5"


# Plot PS and TF:

outputplot = "/mn/stornext/d22/cmbco/comap/nils/COMAP_general/figs/tf_new_az_speed_before_after_diff.pdf"



#PS_plotter(simpath, mappath, noisepath, outname = outputplot)


plot_two_TFs_and_diff(simpath_1, mappath_1, noisepath_1, 
                      simpath_2, mappath_2, noisepath_2, 
                      outputplot)

"""
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
"""