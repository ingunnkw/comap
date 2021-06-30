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
    mapname = sys.argv[1]
    mapname2 = sys.argv[2]
except IndexError:
    print('Missing filename!')
    print('Usage: python ps_script.py mapname mapname2')
    sys.exit(1)

feeds = [1, 2, 3, 5, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]

for i in feeds:
    for j in feeds:
        my_map = map_cosmo.MapCosmo(mapname, feed=i)
        my_map2 = map_cosmo.MapCosmo(mapname2, feed=j)

        tot_rms = np.zeros_like(my_map.rms)
        where = np.where((my_map.rms * my_map2.rms) != 0)
        tot_rms[where] = np.sqrt(my_map.rms[where] ** 2 + my_map2.rms[where] ** 2) 

        summap = np.zeros_like(my_map.map)
        summap[where] = (my_map.map[where] + my_map2.map[where])
        diffmap = np.zeros_like(my_map.map)
        diffmap[where] = (my_map.map[where] - my_map2.map[where])

        my_map.rms = tot_rms
        my_map.map = summap

        sum_ps = power_spectrum.PowerSpectrum(my_map)

        ps, k, nmodes = sum_ps.calculate_ps()

        rms_mean, rms_sig = sum_ps.run_noise_sims(10)

        outname = 'spectra/ps_sum' + my_map.save_string + '_' + my_map2.map_string + '_%02i.h5' % j   

        sum_ps.make_h5(outname=outname)

        ps_sum, k, _ = sum_ps.calculate_ps()

        my_map.map = diffmap
        diff_ps = power_spectrum.PowerSpectrum(my_map)
        ps_diff, k, _ = diff_ps.calculate_ps()

        diff_ps.rms_ps_mean = rms_mean
        diff_ps.rms_ps_std = rms_sig

        outname = 'spectra/ps_diff' + my_map.save_string + '_' + my_map2.map_string + '_%02i.h5' % j            

        diff_ps.make_h5(outname=outname)

        lim = np.mean(np.abs(ps_diff[6:])) * 3

        fig = plt.figure()
        fig.suptitle('feed %02i vs feed %02i' % (i, j))
        ax1 = fig.add_subplot(211)
        ax1.errorbar(k, ps_sum, rms_sig, fmt='o', label=r'$\tilde{P}_{sum}(k)$')
        ax1.errorbar(k, ps_diff, rms_sig, fmt='o', label=r'$\tilde{P}_{diff}(k)$')
        ax1.plot(k, rms_mean, 'k', label=r'$\tilde{P}_{noise}(k)$', alpha=0.4)
        # ax1.plot(k_th / h, ps_th * h ** 3 * 10, 'r--', label=r'$10 * P_{Theory}(k)$')
        ax1.set_ylabel(r'$\tilde{P}(k)$ [$\mu$K${}^2$ Mpc${}^3$]')
        #ax1.set_ylim(0, lim)  # ax1.set_ylim(0, 0.1)
        ax1.set_yscale('log')
        ax1.set_xscale('log')

        plt.legend()

        ax2 = fig.add_subplot(212)
        ax2.errorbar(k, (ps_sum - rms_mean) / rms_sig, rms_sig / rms_sig, fmt='o', label=r'$\tilde{P}_{sum}(k) - \tilde{P}_{noise}(k)$')
        ax2.errorbar(k, (ps_diff - rms_mean) / rms_sig, rms_sig / rms_sig, fmt='o', label=r'$\tilde{P}_{diff}(k) - \tilde{P}_{noise}(k)$')
        ax2.plot(k, 0 * rms_mean, 'k', alpha=0.4)
        ax2.set_ylabel(r'$\tilde{P}(k) / \sigma_\tilde{P}$')
        ax2.set_xlabel(r'$k$ [Mpc${}^{-1}$]')
        ax2.set_ylim(-7, 20)
        ax2.set_xscale('log')
        plt.legend()

        folder = '/mn/stornext/u3/haavarti/www_docs/files/diff/'
        plt.savefig(folder + 'ps_diff' + my_map.save_string + '_' + my_map2.map_string + '_%02i.png' % j, bbox_inches='tight')
        print('Done with %02i, %02i!' % (i, j))
#plt.show()
