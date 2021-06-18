import numpy as np
import matplotlib.pyplot as plt
import numpy.fft as fft
import corner
import h5py
import sys
import scipy.interpolate

import tools
import map_cosmo
import power_spectrum


n_sim = 100
n_k = 14
n_feed = 19

xs = np.zeros((n_feed, n_feed, n_k))
# rms_xs_mean = np.zeros_like(xs)
rms_xs_std = np.zeros_like(xs)

rms_xs = np.zeros((n_feed, n_feed, n_k, n_sim))

chi2 = np.zeros((n_feed, n_feed))
noise = np.zeros_like(chi2)

n_sum = 0

k = np.zeros(n_k)
xs_sum = np.zeros(n_k)
rms_xs_sum = np.zeros((n_k, n_sim))
xs_div = np.zeros(n_k)
for i in range(n_feed):
    for j in range(n_feed):
        try:
            #filepath = 'spectra/xs_feeds_co2_half_0_map_full_%02i_full_%02i.h5' % (i+1, j+1)
#            filepath = 'spectra/xs_feeds_co6_half_0_good_hp_map_%02i_map_%02i.h5' % (i+1, j+1)
            filepath = 'spectra/xs_feeds_par_co7_map_complete_%02i_complete_%02i.h5' % (i+1, j+1)
            with h5py.File(filepath, mode="r") as my_file:
                xs[i, j] = np.array(my_file['xs'][:])
                rms_xs[i, j] = np.array(my_file['rms_xs'][:])
                rms_xs_std[i, j] = np.array(my_file['rms_xs_std'][:])
                k[:] = np.array(my_file['k'][:])
            
            
        except:
            xs[i, j] = np.nan
            rms_xs[i, j] = np.nan
            rms_xs_std[i, j] = np.nan
            
        w = np.sum(1 / rms_xs_std[i,j])
        noise[i,j] = 1 / np.sqrt(w)
        
        chi3 = np.sum((xs[i,j] / rms_xs_std[i,j]) ** 3)

        chi2[i, j] = np.sign(chi3) * (np.sum((xs[i,j] / rms_xs_std[i,j]) ** 2) - n_k) / np.sqrt(2 * n_k)
        if abs(chi2[i,j]) < 5.0 and not np.isnan(chi2[i,j]) and i != j: # and i != 7 and j != 7 and i!=1 and j!=1 and i!=10 and j!=10: # and (i == 4 or j == 4):
            xs_sum += xs[i,j] / rms_xs_std[i,j] ** 2
            rms_xs_sum += rms_xs[i,j] / rms_xs_std[i,j][:,None] ** 2
            xs_div += 1 / rms_xs_std[i,j] ** 2
            n_sum += 1


print(n_sum)
xs_mean = xs_sum / xs_div
rms_xs_mean = rms_xs_sum / xs_div[:,None]

xs_sigma = 1 / np.sqrt(xs_div)

rms_xs_sigma = np.std(rms_xs_mean, 1)


print(xs_mean)
print(xs_sigma)
print(rms_xs_sigma)

plt.figure()
vmax = 15

plt.imshow(chi2, interpolation='none', vmin=-vmax, vmax=vmax, extent=(0.5, n_feed + 0.5, n_feed + 0.5, 0.5))
new_tick_locations = np.array(range(n_feed)) + 1
plt.xticks(new_tick_locations)
plt.yticks(new_tick_locations)
plt.xlabel('Feed')
plt.ylabel('Feed')
cbar = plt.colorbar()
cbar.set_label(r'$|\chi^2| \times$ sign($\chi^3$)')
plt.savefig('xs_grid_par_co7_full.png', bbox_inches='tight')

plt.figure()
vmax = 2e2

plt.imshow(noise, interpolation='none', vmin=vmax/2, vmax=vmax, extent=(0.5, n_feed + 0.5, n_feed + 0.5, 0.5))
new_tick_locations = np.array(range(n_feed)) + 1
plt.xticks(new_tick_locations)
plt.yticks(new_tick_locations)
plt.xlabel('Feed')
plt.ylabel('Feed')
cbar = plt.colorbar()
cbar.set_label(r'noise level [$\mu$K]')
plt.savefig('xs_grid_noise.png', bbox_inches='tight')

k_th = np.load('comap2ps/k.npy')
ps_th = np.load('comap2ps/ps.npy')
ps_th_nobeam = np.load('comap2ps/psn.npy')

print(ps_th)

transfer = scipy.interpolate.interp1d(k_th, ps_th / ps_th_nobeam)

print(transfer(k))

ps_copps = 8.746e3 * ps_th / ps_th_nobeam
ps_copps_nobeam = 8.7e3

lim = np.mean(np.abs(xs_mean[4:-2] * k[4:-2])) * 8
#lim = np.mean(np.abs(xs_mean[4:-2])) * 8

fig = plt.figure()


ax1 = fig.add_subplot(211)
ax1.errorbar(k, k * xs_mean / transfer(k), k * xs_sigma / transfer(k), fmt='o', label=r'$k\tilde{C}_{data}(k)$')
#ax1.errorbar(k, k * xs_mean / transfer(k), k * rms_xs_sigma / transfer(k), fmt='o', label=r'$k\tilde{C}_{data}(k)$')
#ax1.errorbar(k, xs_mean / transfer(k), xs_sigma / transfer(k), fmt='o', label=r'$\tilde{C}_{data}(k)$')
ax1.plot(k, 0 * xs_mean, 'k', alpha=0.4)
#ax1.plot(k, , 'r--', label=r'$100 * P_{Theory}(k)$')
ax1.plot(k_th, k_th * ps_th_nobeam * 10, 'r--', label=r'$10 \times kP_{Theory}(k)$')
#ax1.plot(k_th, k_th * ps_copps * 5, 'g', label=r'$5 \times kP_{COPPS}$ (shot)')
#ax1.plot(k_th, k_th * ps_copps_nobeam * 5, 'g--', label=r'$5 \times kP_{COPPS}$ (shot, no beam)')
#ax1.plot(k_th, k_th * ps_copps * 5, 'g', label=r'$5 \times kP_{COPPS}$ (shot)')
#ax1.plot(k_th, k_th * ps_copps_nobeam * 5, 'g--', label=r'$5 \times kP_{COPPS}$ (shot)')
#ax1.plot(k_th, ps_copps_nobeam * 5 + 0 * k_th, 'g--', label=r'$5 \times P_{COPPS}$ (shot)')

#ax1.plot(k_th, k_th * ps_th_nobeam * 100, 'r', label=r'$100 \times P_{Theory}(k)$')

#ax1.plot(k_th, k_th * ps_th * 100, 'r--', label=r'$100 \times P_{Theory}(k)$ (with beam)')
#ax1.set_ylabel(r'$\tilde{C}(k)$ [$\mu$K${}^2$ Mpc${}^3$]')
ax1.set_ylabel(r'$k\tilde{C}(k)$ [$\mu$K${}^2$ Mpc${}^2$]')
ax1.set_ylim(-lim, lim)              # ax1.set_ylim(0, 0.1)
#ax1.set_ylim(5e2, 2e5)              # ax1.set_ylim(0, 0.1)
#ax1.set_yscale('log')
ax1.set_xscale('log')
ax1.grid()
plt.legend()

ax2 = fig.add_subplot(212)
#ax2.errorbar(k, xs_mean / xs_sigma, xs_sigma / xs_sigma, fmt='o', label=r'$\tilde{C}_{data}(k)$')
ax2.errorbar(k, xs_mean / rms_xs_sigma, rms_xs_sigma / rms_xs_sigma, fmt='o', label=r'$\tilde{C}_{data}(k)$')
ax2.plot(k, 0 * xs_mean, 'k', alpha=0.4)
ax2.set_ylabel(r'$\tilde{C}(k) / \sigma_\tilde{C}$')
ax2.set_xlabel(r'$k$ [Mpc${}^{-1}$]')
ax2.set_ylim(-12, 12)
ax2.set_xscale('log')
ax2.grid()
plt.tight_layout()
plt.legend()

#plt.savefig('xs_mean_co2_full_proper_error.png', bbox_inches='tight')
plt.savefig('xs_mean_par_co7_full.png', bbox_inches='tight')
#plt.savefig('xs_mean_good_proper_error.png', bbox_inches='tight')


plt.show()
