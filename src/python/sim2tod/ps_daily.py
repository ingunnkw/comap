import numpy as np
import matplotlib.pyplot as plt
import numpy.fft as fft
import h5py
import sys
import multiprocessing

import tools
import map_cosmo
import power_spectrum

try:
    mapname = sys.argv[1]
except IndexError:
    print('Missing filename!')
    print('Usage: python ps_script.py mapname')
    sys.exit(1)

save_folder = '/mn/stornext/u3/haavarti/www_docs/diag/ps/'    

prefix = mapname[:-6].rpartition("/")[-1]

my_map = map_cosmo.MapCosmo(mapname)

my_ps = power_spectrum.PowerSpectrum(my_map)

ps, k, _ = my_ps.calculate_ps()

rms_mean, rms_sig = my_ps.run_noise_sims(10)

my_ps.make_h5()

fig = plt.figure()

ax1 = fig.add_subplot(211)
ax1.errorbar(k, ps, rms_sig, fmt='o', label=r'$\tilde{P}_{data}(k)$')
ax1.plot(k, rms_mean, 'k', label=r'$\tilde{P}_{noise}(k)$', alpha=0.4)
ax1.set_ylabel(r'$\tilde{P}(k)$ [$\mu$K${}^2$ Mpc${}^3$]')
# ax1.set_ylim(0, 0.012)#rms_mean.mean() * 3)
ax1.set_yscale('log')
ax1.set_xscale('log')
ax1.grid()
ax1.legend()

ax2 = fig.add_subplot(212)

ax2.errorbar(k, (ps - rms_mean) / rms_sig, rms_sig / rms_sig,
             fmt='o', label=r'$\tilde{P}_{data}(k) - \tilde{P}_{noise}(k)$')
ax2.plot(k, 0 * rms_mean, 'k', alpha=0.4)
ax2.set_ylabel(r'$\tilde{P}(k) / \sigma_\tilde{P}$')
ax2.set_xlabel(r'$k$ [Mpc${}^{-1}$]')
ax2.set_ylim(-7, 20)
ax2.set_xscale('log')
ax2.grid()
plt.legend()


plt.savefig(save_folder + prefix + 'ps.pdf', bbox_inches='tight')

cols = ['b', 'g', 'r', 'c', 'm', 'y', 'k']

linetype = 'o'

fig = plt.figure(figsize=(12, 6))
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212)


def get_feed_ps(feed):
    my_map = map_cosmo.MapCosmo(mapname, feed=feed)

    my_ps = power_spectrum.PowerSpectrum(my_map)

    ps, k, nmodes = my_ps.calculate_ps()

    rms_mean, rms_sig = my_ps.run_noise_sims(10)
    my_ps.make_h5()
    return ps, k, rms_mean, rms_sig


pool = multiprocessing.Pool(20)
ps_arr, k_arr, rms_mean_arr, rms_sig_arr = zip(
    *pool.map(get_feed_ps, range(1, 20)))

for feed in range(1, 20):
    if feed == 8:
        linetype = 'x'
    if feed == 15:
        linetype = 'v'
    ps = ps_arr[feed - 1]
    k = k_arr[feed - 1]
    rms_mean = rms_mean_arr[feed - 1]
    rms_sig = rms_sig_arr[feed - 1]

    col = cols[(feed-1) % 7]

    ax1.errorbar(k, ps, rms_sig, color=col, fmt=linetype,
                 label='Feed %02i' % feed)  # label=r'$P_{data}(k)$')
    # label=r'$P_{noise}(k)$', alpha=0.4)
    ax1.plot(k, rms_mean, color=col, alpha=0.4)
    ax1.set_ylabel(r'$\tilde{P}(k)$ [$\mu$K${}^2$ Mpc${}^3$]')
#    ax1.set_ylim(0, 0.012)#rms_mean.mean() * 3)

    ax2.errorbar(k, (ps - rms_mean) / rms_sig, rms_sig / rms_sig, color=col,
                 fmt=linetype, label=r'$\tilde{P}_{data}(k) - \tilde{P}_{noise}(k)$')
    ax2.set_ylabel(r'$\tilde{P}(k) / \sigma_\tilde{P}$')
    ax2.set_xlabel(r'$k$ [Mpc${}^{-1}$]')
    ax2.set_ylim(-7, 20)

    # print('Done with feed ', feed)
ax1.set_yscale('log')
ax1.set_xscale('log')
ax1.grid()
ax2.plot(k, 0 * rms_mean, 'k', alpha=0.4)
ax2.set_xscale('log')
ax2.grid()
plt.figlegend(*ax1.get_legend_handles_labels(), loc='upper right')
# ax1.legend(loc='upper right')
# ax2.legend()
plt.savefig(save_folder + prefix + 'allpix_ps.pdf', bbox_inches='tight')
# plt.show()
