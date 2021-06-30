import numpy as np
import matplotlib.pyplot as plt
import numpy.fft as fft
from mpi4py import MPI
import time

import tools
import master

comm = MPI.COMM_WORLD
n_list = [10, 15, 20, 25, 30, 35, 40]
time_list = []
# ps_list = []
for n in n_list:

    # n = 50
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
    signal_map = np.random.randn(n_x - 1, n_y - 1, n_z - 1)

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
    before = time.time()
    M_inv, k_full = master.mode_coupling_matrix_3d_mpi(w, k_bin_edges=k_bin_edges, mpi_comm=comm, dx=dx, dy=dy, dz=dz, insert_edge_bins=False)
    after = time.time()
    time_list.append(after - before)
    print(time_list[-1], n)

if comm.rank == 0:
    n_list = np.array(n_list)
    time_list = np.array(time_list)
    # ps_list = np.array(ps_list)
    # np.save('n_list', n_list)
    # np.save('time_list', time_list)
    model_scale = time_list[2] * (n_list / n_list[2]) ** 6
    # print(time_list[2], n_list)
    # print(model_scale)
    # print((n_list / time_list[2]), (25.0 / 5.0) ** 6)
    plt.plot(n_list, time_list, label=r'calclate $M^{-1}$')
    plt.plot(n_list, model_scale, label=r'$a * n^6$')
    plt.legend()
    plt.xlabel('n')
    plt.ylabel('time [s]')
    plt.savefig('time_M_mpi.pdf')

    plt.show()
