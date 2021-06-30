import numpy as np
import numpy.fft as fft
from mpi4py import MPI

import tools


def mode_coupling_matrix_3d(w, k_bin_edges, dx=1.0, dy=1.0, dz=1.0, insert_edge_bins=False):
    w_k = fft.fftn(w)
    n_x, n_y, n_z = w.shape

    k_x_values = fft.fftfreq(n_x, dx)
    k_y_values = fft.fftfreq(n_y, dy)
    k_z_values = fft.fftfreq(n_z, dz)

    if insert_edge_bins:
        n_k = len(k_bin_edges) + 1
        kmax = np.sqrt(np.max(np.abs(k_x_values)) ** 2 + np.max(np.abs(k_y_values)) ** 2 +
                       np.max(np.abs(k_z_values)) ** 2)
        full_k_bins = np.zeros(n_k + 1)
        full_k_bins[0] = -1e-5
        full_k_bins[-1] = kmax + 1e-5
        full_k_bins[1:-1] = k_bin_edges
    else:
        n_k = len(k_bin_edges) - 1
        full_k_bins = k_bin_edges
    
    M = np.zeros((n_k, n_k))
    n_count = np.zeros((n_k, n_k))

    k_2_grid = np.sqrt(sum(ki ** 2 for ki in np.meshgrid(k_x_values, k_y_values, k_z_values, indexing='ij')))
    for i in range(n_x):
        for j in range(n_y):
            for k in range(n_z):
                k_1 = np.sqrt(k_x_values[i] ** 2 + k_y_values[j] ** 2 + k_z_values[k] ** 2)

                ind1 = np.digitize([k_1], full_k_bins) - 1
                # print(k_1, np.max(np.abs(full_k_bins)), i, j, k)
                n_count[ind1, :] += 1

                l_indices = (i - np.arange(n_x)) % n_x
                m_indices = (j - np.arange(n_y)) % n_y
                n_indices = (k - np.arange(n_z)) % n_z
                grid = np.array([ki for ki in np.meshgrid(l_indices, m_indices, n_indices, indexing='ij')])

                k_1_grid = k_1 + 0.0 * np.zeros_like(k_2_grid)
                k_tilde_grid = [grid[0], grid[1], grid[2]]
 
                weights = (np.abs(w_k[k_tilde_grid]) ** 2 / (n_x ** 2 * n_y ** 2 * n_z ** 2)).flatten()
                M += np.histogram2d(k_1_grid.flatten(), k_2_grid.flatten(), bins=full_k_bins, weights=weights)[0]
    M /= n_count
    print(n_count[:, 0])
    return np.linalg.inv(M), full_k_bins


def mode_coupling_matrix_3d_mpi(w, k_bin_edges, mpi_comm, dx=1.0, dy=1.0, dz=1.0, insert_edge_bins=False):
    w_k = fft.fftn(w)
    n_x, n_y, n_z = w.shape

    k_x_values = fft.fftfreq(n_x, dx)
    k_y_values = fft.fftfreq(n_y, dy)
    k_z_values = fft.fftfreq(n_z, dz)

    if insert_edge_bins:
        n_k = len(k_bin_edges) + 1
        kmax = np.sqrt(np.max(np.abs(k_x_values)) ** 2 + np.max(np.abs(k_y_values)) ** 2 +
                       np.max(np.abs(k_z_values)) ** 2)
        full_k_bins = np.zeros(n_k + 1)
        full_k_bins[0] = -1e-5
        full_k_bins[-1] = kmax + 1e-5
        full_k_bins[1:-1] = k_bin_edges
    else:
        n_k = len(k_bin_edges) - 1
        full_k_bins = k_bin_edges

    M = np.zeros((n_k, n_k))
    n_count = np.zeros((n_k, n_k))

    k_2_grid = np.sqrt(sum(ki ** 2 for ki in np.meshgrid(k_x_values, k_y_values, k_z_values, indexing='ij')))

    my_indices = tools.distribute_indices(n_x * n_y * n_z, mpi_comm.size, mpi_comm.rank)
    for p in my_indices:
        i = p // (n_y * n_z)  ## python 2 syntax ? 
        j = (p - i * (n_y * n_z)) // n_z
        k = p % n_z
        k_1 = np.sqrt(k_x_values[i] ** 2 + k_y_values[j] ** 2 + k_z_values[k] ** 2)

        ind1 = np.digitize([k_1], full_k_bins) - 1
        n_count[ind1, :] += 1

        l_indices = (i - np.arange(n_x)) % n_x
        m_indices = (j - np.arange(n_y)) % n_y
        n_indices = (k - np.arange(n_z)) % n_z
        grid = np.array([ki for ki in np.meshgrid(l_indices, m_indices, n_indices, indexing='ij')])

        k_1_grid = k_1 + 0.0 * np.zeros_like(k_2_grid)
        k_tilde_grid = [grid[0], grid[1], grid[2]]

        weights = (np.abs(w_k[k_tilde_grid]) ** 2 / (n_x ** 2 * n_y ** 2 * n_z ** 2)).flatten()
        M += np.histogram2d(k_1_grid.flatten(), k_2_grid.flatten(), bins=full_k_bins, weights=weights)[0]
    M_sum = mpi_comm.reduce(M, op=MPI.SUM, root=0)
    n_count_sum = mpi_comm.reduce(n_count, op=MPI.SUM, root=0)
    if mpi_comm.rank == 0:
        M = M_sum / n_count_sum
    mpi_comm.Bcast([M, MPI.DOUBLE])
    # print "\n\n\n\n\n", mpi_comm.rank, np.linalg.inv(M)
    return np.linalg.inv(M), full_k_bins


def mode_coupling_matrix_2d(w, k_bin_edges, dx=1.0, dy=1.0, insert_edge_bins=False):
    w_k = fft.fftn(w)
    n_x, n_y = w.shape

    k_x_values = fft.fftfreq(n_x, dx)
    k_y_values = fft.fftfreq(n_y, dy)

    if insert_edge_bins:
        n_k = len(k_bin_edges) + 1
        kmax = np.sqrt(np.max(np.abs(k_x_values)) ** 2 + np.max(np.abs(k_y_values)) ** 2)
        full_k_bins = np.zeros(n_k + 1)
        full_k_bins[0] = -1e-5
        full_k_bins[-1] = kmax + 1e-5
        full_k_bins[1:-1] = k_bin_edges
    else:
        n_k = len(k_bin_edges) - 1
        full_k_bins = k_bin_edges

    M = np.zeros((n_k, n_k))
    n_count = np.zeros((n_k, n_k))

    k_2_grid = np.sqrt(sum(ki ** 2 for ki in np.meshgrid(k_x_values, k_y_values, indexing='ij')))
    for i in range(n_x):
        for j in range(n_y):
            k_1 = np.sqrt(k_x_values[i] ** 2 + k_y_values[j] ** 2)

            ind1 = np.digitize([k_1], full_k_bins) - 1
            n_count[ind1, :] += 1

            l_indices = (i - np.arange(n_x)) % n_x
            m_indices = (j - np.arange(n_y)) % n_y
            grid = np.array([ki for ki in np.meshgrid(l_indices, m_indices, indexing='ij')])

            k_1_grid = k_1 + 0.0 * np.zeros_like(k_2_grid)
            k_tilde_grid = [grid[0], grid[1]]

            weights = (np.abs(w_k[k_tilde_grid]) ** 2 / (n_x ** 2 * n_y ** 2)).flatten()
            M += np.histogram2d(k_1_grid.flatten(), k_2_grid.flatten(), bins=full_k_bins, weights=weights)[0]
    M /= n_count
    print(n_count[:, 0])
    return np.linalg.inv(M), full_k_bins

