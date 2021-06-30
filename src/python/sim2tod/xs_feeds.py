import numpy as np
import matplotlib.pyplot as plt
import numpy.fft as fft
import corner
import h5py
import sys

import tools
import map_cosmo
import power_spectrum
import multiprocessing

try:
    mapname = sys.argv[1]
    # mapname2 = sys.argv[2]
except IndexError:
    print('Missing filename!')
    print('Usage: python ps_script.py mapname')  # mapname2')
    sys.exit(1)
#feeds = [1, 2, 3, 5, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]

nums = list(range(19 * 19))
#for i in feeds:
#    for j in feeds:
def calc_feed_feed_cross_spectrum(p):
    
    i = p // 19 + 1
    j = p % 19 + 1
    if i == 4 or i == 6 or i == 7:
        return p
    if j == 4 or j == 6 or j == 7:
        return p
    
    my_map = map_cosmo.MapCosmo(mapname, feed=i, jk='half', split=0)
    my_map2 = map_cosmo.MapCosmo(mapname, feed=j, jk='half', split=1)
    print(my_map.rms.shape)
    print(my_map.rms[28, 68, :])
    my_xs = power_spectrum.CrossSpectrum(my_map, my_map2)
    
    xs, k, nmodes = my_xs.calculate_xs()
    
    rms_mean, rms_sig = my_xs.run_noise_sims(100, seed=42)

    outname = 'spectra/xs_feeds_par' + my_map.save_string + '_' + my_map2.map_string + '_%02i.h5' % j   
    my_xs.make_h5(outname=outname, save_noise_realizations=True)

    lim = np.mean(np.abs(xs[4:])) * 4
    #print(xs)
    fig = plt.figure()
    fig.suptitle('feed %02i x feed %02i' % (i, j))
    ax1 = fig.add_subplot(211)
    ax1.errorbar(k, xs, rms_sig, fmt='o', label=r'$\tilde{C}_{data}(k)$')
    ax1.plot(k, 0 * rms_mean, 'k', label=r'$\tilde{C}_{noise}(k)$', alpha=0.4)
    #ax1.plot(k, theory[1:] * 100, 'r--', label=r'$100 * P_{Theory}(k)$')
    #ax1.plot(k_th, ps_th * 10, 'r--', label=r'$10 * P_{Theory}(k)$')
    ax1.set_ylabel(r'$\tilde{C}(k)$ [$\mu$K${}^2$ Mpc${}^3$]')
    ax1.set_ylim(-lim, lim)              # ax1.set_ylim(0, 0.1)
    #ax1.set_yscale('log')
    ax1.set_xscale('log')
    ax1.grid()
    plt.legend()

    ax2 = fig.add_subplot(212)
    ax2.errorbar(k, xs / rms_sig, rms_sig / rms_sig, fmt='o', label=r'$\tilde{C}_{data}(k)$')
    ax2.plot(k, 0 * rms_mean, 'k', alpha=0.4)
    ax2.set_ylabel(r'$\tilde{C}(k) / \sigma_\tilde{C}$')
    ax2.set_xlabel(r'$k$ [Mpc${}^{-1}$]')
    ax2.set_ylim(-12, 12)
    ax2.set_xscale('log')
    ax2.grid()
    plt.legend()

    folder = '/mn/stornext/d16/www_cmb/comap/xs/'
    plt.savefig(folder + 'xs_co7_par_%02i_%02i' % (i,j) + my_map.map_string + '_' + my_map2.map_string + '.png', bbox_inches='tight')
    print('Done with %02i, %02i!' % (i, j))
    return p

pool = multiprocessing.Pool(64)
np.array(pool.map(calc_feed_feed_cross_spectrum, nums))
        

