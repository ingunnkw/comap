import numpy as np
import matplotlib.pyplot as plt
import numpy.fft as fft
import corner
import h5py
import sys

import tools
import map_cosmo
import power_spectrum


try:
    mappath = sys.argv[1]
except IndexError:
    print('Missing filename!')
    print('Usage: python ps_2d.py mappath')
    sys.exit(1)

my_map = map_cosmo.MapCosmo(mappath)

my_ps = power_spectrum.PowerSpectrum(my_map)

ps_2d, k, nmodes = my_ps.calculate_ps(do_2d=True)

#my_ps.make_h5(outname = "/mn/stornext/d16/cmbco/comap/nils/COMAP_general/data/maps/PS/out.h5")

fig, ax = plt.subplots(1,1)
img = ax.imshow(np.log10(ps_2d), interpolation='none', origin='lower', extent=[0,1,0,1])
#plt.imshow(np.log10(nmodes), interpolation='none', origin='lower')
cbar = fig.colorbar(img)
cbar.set_label(r'$\log_{10}(\tilde{P}_{\parallel, \bot}(k))$ [$\mu$K${}^2$ (Mpc)${}^3$]')

def log2lin(x, k_edges):
    loglen = np.log10(k_edges[-1]) - np.log10(k_edges[0])
    logx = np.log10(x) - np.log10(k_edges[0])
    return logx / loglen


# ax.set_xscale('log')
minorticks = [0.0002, 0.0003, 0.0004, 0.0005, 0.0006, 0.0007, 0.0008, 0.0009,
              0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009,
              0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09,
              0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
              2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0,
              20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0,
              200.0, 300.0, 400.0, 500.0, 600.0, 700.0, 800.0, 900.0]

majorticks = [1.e-04, 1.e-03, 1.e-02, 1.e-01, 1.e+00, 1.e+01, 1.e+02]
majorlabels = ['$10^{-4}$', '$10^{-3}$', '$10^{-2}$', '$10^{-1}$', '$10^{0}$', '$10^{1}$', '$10^{2}$']

xbins = my_ps.k_bin_edges_par

ticklist_x = log2lin(minorticks, xbins)
majorlist_x = log2lin(majorticks, xbins)

ybins = my_ps.k_bin_edges_perp

ticklist_y = log2lin(minorticks, ybins)
majorlist_y = log2lin(majorticks, ybins)


ax.set_xticks(ticklist_x, minor=True)
ax.set_xticks(majorlist_x, minor=False)
ax.set_xticklabels(majorlabels, minor=False)
ax.set_yticks(ticklist_y, minor=True)
ax.set_yticks(majorlist_y, minor=False)
ax.set_yticklabels(majorlabels, minor=False)

plt.xlabel(r'$k_{\parallel}$')
plt.ylabel(r'$k_{\bot}$')
plt.xlim(0, 1)
plt.ylim(0, 1)
#plt.savefig('ps_par_vs_perp_nmodes.png')
plt.savefig('ps_par_vs_perp.png')
plt.show()

