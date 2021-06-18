import numpy as np
import matplotlib.pyplot as plt
import numpy.fft as fft
from mpi4py import MPI

import tools
import master


def p(k):
    k0 = 1e0
    return 0.0 * k + 1.0 #1.0 / (k + 0.001) ** (-0.5) * (k > 0) #4.0 * (k < 5) + 0.2 * (k > 5) - 0.19 * (k > 20) #  1 + 1.0/(k + 0.1) ** 9

n = 30
# p = 10
n_x = n
n_y = n
n_z = n  # n

n_k = 15

x = np.linspace(0, 5, n_x)
y = np.linspace(0, 5, n_y)
z = np.linspace(0, 5, n_z)
dx = x[1] - x[0]
dy = y[1] - y[0]
dz = z[1] - z[0]

signal_map = tools.create_map_3d(p, x, y, z)

x_ind, y_ind, z_ind = np.indices(signal_map.shape)

r = np.hypot(x[x_ind] - 2, y[y_ind] - 2, z[z_ind] - 2)

rms_map = (r / np.max(r.flatten()) + 0.05) * np.std(signal_map.flatten()) ** 2.5 / 5.0

datamap = signal_map + rms_map * np.random.randn(*rms_map.shape)

w = np.mean(rms_map.flatten() ** 2) / rms_map ** 2

kmax = np.sqrt(
    np.max(np.abs(fft.fftfreq(n_x - 1, dx))) ** 2
    + np.max(np.abs(fft.fftfreq(n_y - 1, dy))) ** 2
    + np.max(np.abs(fft.fftfreq(n_z - 1, dz))) ** 2
)

k_bin_edges = np.linspace(0 - 1e-5, kmax + 1e-5, n_k)

# lim = 100

# plt.figure()
# plt.imshow(rms_map, interpolation='none', vmin=0, vmax=lim)
# plt.colorbar()
# plt.figure()
# im = plt.imshow(signal_map, interpolation='none', vmin=-lim, vmax=lim)
# plt.colorbar()
# plt.figure()
# plt.imshow(datamap, interpolation='none', vmin=-lim, vmax=lim)
# plt.colorbar()
print('before')
M_inv, k_full = master.mode_coupling_matrix_3d(w, k_bin_edges=k_bin_edges, dx=dx, dy=dy, dz=dz, insert_edge_bins=False)
print('between')
ps_signal, k, _ = tools.compute_power_spec3d(signal_map, k_bin_edges, dx=dx, dy=dy, dz=dz)
ps_with_noise, _, _ = tools.compute_power_spec3d(datamap, k_bin_edges, dx=dx, dy=dy, dz=dz)
ps_with_weight, _, _ = tools.compute_power_spec3d(datamap * w, k_full, dx=dx, dy=dy, dz=dz)
print('after')

n_sim = 1000
rms_ps = np.zeros((len(k_full) - 1, n_sim))
rms_ps_noweight = np.zeros((len(k), n_sim))
for i in range(n_sim):
    rms_ps[:, i] = tools.compute_power_spec3d(
        rms_map * np.random.randn(*rms_map.shape) * w,
        k_full, dx=dx, dy=dy, dz=dz
    )[0]
    rms_ps_noweight[:, i] = tools.compute_power_spec3d(
        rms_map * np.random.randn(*rms_map.shape),
        k_bin_edges, dx=dx, dy=dy, dz=dz
    )[0]
rms_ps_mean = rms_ps.mean(1)
rms_ps_mean_noweight = rms_ps_noweight.mean(1)

ps_estimated = np.dot(M_inv, ps_with_weight - rms_ps_mean)#[1:-1]

ps_signal_mean = np.zeros((len(k), n_sim))
data_ps = np.zeros((len(k), n_sim))
data_ps_simple = np.zeros((len(k), n_sim))
for i in range(n_sim):
    signal_map = tools.create_map_3d(p, x, y, z)
    datamap = signal_map + rms_map * np.random.randn(*rms_map.shape)
    ps_with_weight, _, _ = tools.compute_power_spec3d(datamap * w, k_full, dx=dx, dy=dy, dz=dz)
    data_ps[:, i] = np.dot(M_inv, ps_with_weight - rms_ps_mean)#[1:-1]
    ps_signal_mean[:, i], _, _ = tools.compute_power_spec3d(signal_map, k_bin_edges, dx=dx, dy=dy, dz=dz)
    data_ps_simple[:, i], _, _ = tools.compute_power_spec3d(datamap, k_bin_edges, dx=dx, dy=dy, dz=dz)
data_ps_mean = np.mean(data_ps, axis=1)
data_ps_std = np.std(data_ps, axis=1)

ps_simple_subtracted = np.mean(data_ps_simple, axis=1) - rms_ps_mean_noweight

plt.figure()
# plt.loglog(k, ps_estimated, 'b', label='master')
# plt.loglog(k, ps_with_noise, 'g', label='signal + noise')
plt.loglog(k, ps_simple_subtracted, 'g', label='signal + noise - <noise>')

plt.loglog(k, data_ps_mean + data_ps_std, 'k--', label='master + std')
plt.loglog(k, data_ps_mean, 'k', label='master mean')
plt.loglog(k, data_ps_mean - data_ps_std, 'k--', label='master - std')
plt.loglog(k, p(k), 'r--', label='analytic')
plt.legend()
plt.axis([k[0], k[-1], 2e-2, 1e2])
plt.show()