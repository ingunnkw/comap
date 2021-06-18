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
    print('Missing filepath!')
    print('Usage: python ps_script.py mappath')
    sys.exit(1)

my_map = map_cosmo.MapCosmo(mappath)

my_ps = power_spectrum.PowerSpectrum(my_map)

ps, k, nmodes = my_ps.calculate_ps()

rms_mean, rms_sig = my_ps.run_noise_sims(10)

my_ps.make_h5()

fig = plt.figure()

ax1 = fig.add_subplot(211)
ax1.errorbar(k,  ps, rms_sig, fmt='o', label=r'$\tilde{P}_{data}(k)$')
ax1.plot(k, rms_mean, 'k', label=r'$\tilde{P}_{noise}(k)$', alpha=0.4)
# ax1.plot(k_th, ps_th * 10, 'r--', label=r'$10 * P_{Theory}(k)$')
ax1.set_ylabel(r'$\tilde{P}(k)$ [$\mu$K${}^2$ (Mpc)${}^3$]')
ax1.set_ylim(1e5, 1e8)  # ax1.set_ylim(0, 0.1)
ax1.set_yscale('log')
ax1.set_xscale('log')
ax1.grid()
plt.legend()

ax2 = fig.add_subplot(212)
ax2.errorbar(k, (ps - rms_mean) / rms_sig, rms_sig / rms_sig, fmt='o', label=r'$\tilde{P}_{data}(k) - \tilde{P}_{noise}(k)$')
ax2.plot(k, 0 * rms_mean, 'k', alpha=0.4)
ax2.set_ylabel(r'$\tilde{P}(k) / \sigma_\tilde{P}$')
ax2.set_xlabel(r'$k$ [Mpc${}^{-1}$]')
#ax2.set_ylim(-7, 20)
ax2.set_ylim(-7, 100)
ax2.set_xscale('log')
ax2.grid()
plt.legend()


plt.savefig('ps' + my_map.save_string + '.png', bbox_inches='tight')
plt.show()
