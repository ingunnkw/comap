import numpy as np
import matplotlib.pyplot as plt
import numpy.fft as fft
from mpi4py import MPI
import corner
import h5py

import tools
import master
import comap2ps

deg2mpc = 76.22  # at redshift 2.9
dz2mpc = 699.62  # redshift 2.4 to 3.4
#freq = np.linspace(26, 34, 257)
dz = (1+2.9) ** 2 * 32.2e-3 / 115
redshift = np.linspace(2.9 - 128*dz, 2.9 + 128*dz, 257)
n = 20
# p = 10
n_x = n + 1
n_y = n + 1
n_z = 65  # n

n_k = 10

x = np.linspace(0, 2, n_x)
y = np.linspace(0, 2, n_y)
z = np.linspace(0, 10, n_z)
dx = x[1] - x[0]
dy = y[1] - y[0]
dz = z[1] - z[0]

signal_map = np.random.randn(n, n, 64)

x_ind, y_ind, z_ind = np.indices(signal_map.shape)

r = np.hypot(x[x_ind] - 1, y[y_ind] - 1, z[z_ind] - 5)

rms_map = (r / np.max(r.flatten()) + 0.05) * np.std(signal_map.flatten()) ** 2.5 / 5.0
datamap = rms_map * np.random.randn(*rms_map.shape)  # + signal_map
my_mask = np.zeros_like(datamap)
my_mask[(rms_map < 0.21)] = 1.0
# print(datamap)
outname = 'map.h5'
f2 = h5py.File(outname, 'w')
f2.create_dataset('rms', data=rms_map)
f2.create_dataset('map', data=datamap)
f2.create_dataset('mask', data=my_mask)
f2.create_dataset('x', data=tools.edge2cent(x))
f2.create_dataset('y', data=tools.edge2cent(y))
#f2.create_dataset('', data=rms_map)
f2.close()