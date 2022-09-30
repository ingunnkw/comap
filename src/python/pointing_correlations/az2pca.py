import numpy as np
import h5py
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

import os
import sys
from sklearn.decomposition import PCA
from tqdm import tqdm

# import multiprocessing as mp
import warnings
from mpi4py import MPI
import time

warnings.filterwarnings(
    "ignore", category=RuntimeWarning
)  # ignore warnings caused by weights cut-off

field = "co2"
path = f"/mn/stornext/d22/cmbco/comap/protodir/level2/Ka/{field}/"
allfiles = os.listdir(path)
allfiles = [file for file in allfiles if f"{field}" in file]
allids = np.array([int(file[4:-3]) for file in allfiles])
# sort = np.argsort(allids)
# allids = allids[sort]
# allfiles = np.array(allfiles)[sort]
Nfiles = len(allfiles)

savedir = f"/mn/stornext/d22/cmbco/comap/nils/COMAP_general/data/level3/{field}/"

# def ddaz_binned(allnames):
def ddaz_binned(idx):
    # pid = mp.current_process()._identity[0]
    # ids = []
    # with ddaz_binned.lock:
    #    print(
    #     "Progress: ~",
    #     ddaz_binned.iterr.value,
    #     "/",
    #     Nfiles,
    #     "|",
    #     round(ddaz_binned.iterr.value / Nfiles * 100, 4),
    #     "%",
    # )
    # ddaz_binned.iterr.value += 1
    t0 = time.time()
    ncomps = 10
    nbins = 300
    nsb = 4
    nchannels = 64
    nfreqs = nsb * nchannels
    nfeeds = 20

    compsdata = np.zeros((nfeeds, ncomps, nbins))
    ampldata = np.zeros((nfeeds, ncomps, nfreqs))

    compsdata_daz = np.zeros((nfeeds, ncomps, nbins))
    ampldata_daz = np.zeros((nfeeds, ncomps, nfreqs))

    compsdata_ddaz = np.zeros((nfeeds, ncomps, nbins))
    ampldata_ddaz = np.zeros((nfeeds, ncomps, nfreqs))

    name = allfiles[idx]
    id = int(name[4:-3])
    # print("Defining variables", time.time() - t0, "sec", idx)
    t0 = time.time()
    try:
        with h5py.File(path + name, "r") as infile:
            # allids.append(id)
            tod = infile["tod"][()]
            freq = infile["nu"][0, ...]
            rms = infile["sigma0"][()]
            pixels = infile["pixels"][()]
            point_tel = infile["point_tel"][()]
            point_cel = infile["point_cel"][()]
            Time = infile["time"][()]
            feature = infile["feature"][()]

    except:
        print("Skipping:", name)
        return None

    # print("Reading data", time.time() - t0, "sec", idx)
    t0 = time.time()

    freq[0, :] = freq[0, ::-1]
    freq[2, :] = freq[2, ::-1]

    tod[:, 0, :, :] = tod[:, 0, ::-1, :]
    tod[:, 2, :, :] = tod[:, 2, ::-1, :]

    _temp = tod.copy()
    tod = np.zeros((nfeeds, nsb, nchannels, tod.shape[-1]))
    tod[pixels - 1, ...] = _temp
    tod = np.where(np.isnan(tod) == False, tod, 0)

    rms[:, 0, :] = rms[:, 0, ::-1]
    rms[:, 2, :] = rms[:, 2, ::-1]
    _temp = rms.copy()
    rms = np.zeros((nfeeds, *rms.shape[1:]))
    rms[pixels - 1, ...] = _temp

    _temp = point_tel.copy()
    point_tel = np.zeros((nfeeds, *point_tel.shape[1:]))
    point_tel[pixels - 1, ...] = _temp

    _temp = point_cel.copy()
    point_cel = np.zeros((nfeeds, *point_cel.shape[1:]))
    point_cel[pixels - 1, ...] = _temp

    az = point_tel[:, :, 0]
    el = point_tel[:, :, 1]

    ra = point_cel[:, :, 0]
    dec = point_cel[:, :, 1]

    # print("Formatting data", time.time() - t0, "sec", idx)
    t0 = time.time()

    average_mjd = np.mean(Time)
    average_el = np.mean(el)
    average_ra = np.mean(ra)
    average_dec = np.mean(dec)

    Nsamp = ra[0, :].shape[0]

    t0 = time.time()
    freq_idx = np.arange(nfreqs, dtype=int)
    # freq_idx = freq_idx[None, :, None] * np.ones((nfeeds, nfreqs, Nsamp))
    # freq_idx = freq_idx.reshape(nfeeds, nfreqs * Nsamp)

    daz = np.gradient(az, (Time[1] - Time[0]) * 3600 * 24, axis=-1)
    ddaz = np.gradient(daz, (Time[1] - Time[0]) * 3600 * 24, axis=-1)

    # fig, ax = plt.subplots(2, 1)

    # ax[0].plot((az - az[5]).T)
    # ax[1].plot((el - el[5]).T)

    # fig, ax = plt.subplots(2, 1)

    # ax[0].plot((az).T)
    # ax[1].plot((el).T)

    # fig, ax = plt.subplots()

    # ax.plot((daz - daz[5]).T)
    # plt.show()
    # sys.exit()
    bins = np.linspace(az[0].min(), az[0].max(), nbins + 1)
    resolution_az = (az[0].max() - az[0].min()) / (nbins)
    az_idx = np.floor((az[0] - az[0].min()) / resolution_az)
    az_idx = np.clip(az_idx, 0, nbins - 1)
    # az_idx = az_idx[None, None, :] * np.ones((nfeeds, nfreqs, Nsamp))
    # az_idx = az_idx.reshape(nfeeds, nfreqs * Nsamp)
    az_bin_idx = np.ones((nfeeds, 1, 1), dtype=int) * (
        freq_idx[None, :, None] * nbins + az_idx[None, None, :]
    ).astype(int)
    az_bin_idx = az_bin_idx.reshape(nfeeds, nfreqs * Nsamp)
    # az_bin_idx = (freq_idx * nbins + az_idx).astype(int)

    dazbins = np.linspace(-0.5, 0.5, nbins + 1)
    resolution_daz = (dazbins.max() - dazbins.min()) / (nbins)
    daz_idx = np.floor((daz[0] - dazbins.min()) / resolution_daz)
    daz_idx = np.clip(daz_idx, 0, nbins - 1)
    # daz_idx = daz_idx[None, None, :] * np.ones((nfeeds, nfreqs, Nsamp))
    # daz_idx = daz_idx.reshape(nfeeds, nfreqs * Nsamp)
    daz_bin_idx = np.ones((nfeeds, 1, 1), dtype=int) * (
        freq_idx[None, :, None] * nbins + daz_idx[None, None, :]
    ).astype(int)
    daz_bin_idx = daz_bin_idx.reshape(nfeeds, nfreqs * Nsamp)

    ddazbins = np.linspace(-10.0, 10.0, nbins + 1)
    resolution_ddaz = (ddazbins.max() - ddazbins.min()) / (nbins)
    ddaz_idx = np.floor((ddaz[0] - ddazbins.min()) / resolution_ddaz)
    ddaz_idx = np.clip(ddaz_idx, 0, nbins - 1)
    # ddaz_idx = ddaz_idx[None, None, :] * np.ones((nfeeds, nfreqs, Nsamp))
    # ddaz_idx = ddaz_idx.reshape(nfeeds, nfreqs * Nsamp)
    ddaz_bin_idx = np.ones((nfeeds, 1, 1), dtype=int) * (
        freq_idx[None, :, None] * nbins + ddaz_idx[None, None, :]
    ).astype(int)
    ddaz_bin_idx = ddaz_bin_idx.reshape(nfeeds, nfreqs * Nsamp)

    for feed in range(nfeeds):
        # allaz = []
        # alldaz = []
        # allddaz = []

        # tb = time.time()

        # for i in range(4):
        #     for j in range(64):

        #         # t0 = time.time()

        #         histsum, bins = np.histogram(
        #             az[feed], bins=nbins, weights=(tod[feed, i, j, :] / rms[feed, i, j])
        #         )
        #         # print(bins[1] - bins[0], resolution_az)
        #         # print(bins[0], bins[-1], az[0].min(), az[0].max())
        #         nhit = np.histogram(az[feed], bins=nbins)[0]
        #         normhist = histsum / nhit * np.sqrt(nhit)
        #         allaz.append(normhist)
        # # print("Az binning time", time.time() - t0, "sec", idx)
        #         # t0 = time.time()
        #         # # all_hits_az.append(nhit)

        #         daz = np.gradient(az[feed], (Time[1] - Time[0]) * 3600 * 24)
        #         # daz = (az[feed, 1:] - az[feed, :-1]) / (Time[1] - Time[0]) * 3600 * 24
        #         bins = np.linspace(-1.0, 1.0, nbins + 1)
        #         bins[0] = -1e6  # Infinite bin edges
        #         bins[-1] = 1e6
        #         histsum, dazbins = np.histogram(
        #             daz, bins=bins, weights=(tod[feed, i, j, :] / rms[feed, i, j])
        #         )
        #         nhit = np.histogram(daz, bins=bins)[0]
        #         normhist = histsum / nhit * np.sqrt(nhit)
        #         alldaz.append(normhist)
        #         # all_hits_daz.append(nhit)
        # print("DAz binning time", time.time() - t0, "sec", idx)
        #         # t0 = time.time()

        #         ddaz = np.gradient(daz, (Time[1] - Time[0]) * 3600 * 24)
        #         # ddaz = (az[feed, 2:] - 2 * az[feed, 1:-1] + az[feed, 0:-2]) / ((Time[1] - Time[0]) * 3600 * 24) ** 2
        #         bins = np.linspace(-10, 10, nbins + 1)
        #         bins[0] = -1e6
        #         bins[-1] = 1e6
        #         histsum, ddazbins = np.histogram(
        #             ddaz, bins=bins, weights=(tod[feed, i, j, :] / rms[feed, i, j])
        #         )
        #         nhit = np.histogram(ddaz, bins=bins)[0]
        #         normhist = histsum / nhit * np.sqrt(nhit)
        #         allddaz.append(normhist)
        #         # all_hits_ddaz.append(nhit)
        # print("DDAz binning time", time.time() - t0, "sec", idx)
        #         # t0 = time.time()

        # allaz = np.array(allaz)
        # alldaz = np.array(alldaz)
        # allddaz = np.array(allddaz)
        # print("Binning time", time.time() - tb, "sec", idx)
        t0 = time.time()

        # az_array = az[feed, None, :] * np.ones((nfreqs, Nsamp))

        # az_idx = np.digitize(
        #     az_array,
        #     # _bins[:-1],
        #     np.linspace(az[feed].min(), az[feed].max(), nbins + 1)[:-1],
        #     right=False,
        # )
        # az_idx -= 1  # Shifting to base 0 indexing

        # az_idx = az_idx.flatten()
        # freq_idx = freq_idx.flatten()
        # bin_idx = (freq_idx * nbins + az_idx).astype(int)

        hitsum = np.bincount(
            az_bin_idx[feed, ...],
            # bin_idx,
            weights=(tod[feed, :, :, :] / rms[feed, :, :, None]).flatten(),
            minlength=nbins * nfreqs,
        )

        # hits = np.bincount(bin_idx, minlength=nbins * nfreqs)
        hits = np.bincount(az_bin_idx[feed, ...], minlength=nbins * nfreqs)
        allaz = hitsum / hits * np.sqrt(hits)
        allaz = allaz.reshape(nfreqs, nbins)

        hitsum_daz = np.bincount(
            daz_bin_idx[feed, ...],
            # bin_idx,
            weights=(tod[feed, :, :, :] / rms[feed, :, :, None]).flatten(),
            minlength=nbins * nfreqs,
        )

        # hits = np.bincount(bin_idx, minlength=nbins * nfreqs)
        hits_daz = np.bincount(daz_bin_idx[feed, ...], minlength=nbins * nfreqs)
        alldaz = hitsum_daz / hits_daz * np.sqrt(hits_daz)
        alldaz = alldaz.reshape(nfreqs, nbins)

        hitsum_ddaz = np.bincount(
            ddaz_bin_idx[feed, ...],
            # bin_idx,
            weights=(tod[feed, :, :, :] / rms[feed, :, :, None]).flatten(),
            minlength=nbins * nfreqs,
        )

        # hits = np.bincount(bin_idx, minlength=nbins * nfreqs)
        hits_ddaz = np.bincount(ddaz_bin_idx[feed, ...], minlength=nbins * nfreqs)
        allddaz = hitsum_ddaz / hits_ddaz * np.sqrt(hits_ddaz)
        allddaz = allddaz.reshape(nfreqs, nbins)
        # print("Test Binning time", time.time() - t0, "sec", idx)
        # # print(az_idx)
        # # print(az_idx.shape, tod[0, 0, 0, :].shape)

        # sys.exit()

        t0 = time.time()

        # fig, ax = plt.subplots(3, 1, figsize=(12, 12))
        # im1 = ax[0].imshow(
        #     allaz.T,
        #     cmap="RdBu_r",
        #     vmin=-0.8 * np.nanmax(np.abs(allaz)),
        #     vmax=0.8 * np.nanmax(np.abs(allaz)),
        # )
        # divider = make_axes_locatable(ax[0])
        # cax = divider.append_axes("right", size="5%", pad=0.05)
        # fig.colorbar(im1, cax=cax, orientation="vertical")

        # im2 = ax[1].imshow(
        #     test_normhist.T,
        #     cmap="RdBu_r",
        #     vmin=-0.8 * np.nanmax(np.abs(allaz)),
        #     vmax=0.8 * np.nanmax(np.abs(allaz)),
        # )

        # divider = make_axes_locatable(ax[1])
        # cax = divider.append_axes("right", size="5%", pad=0.05)
        # fig.colorbar(im2, cax=cax, orientation="vertical")

        # im3 = ax[2].imshow(
        #     (test_normhist - allaz).T,
        #     cmap="RdBu_r",
        #     vmin=-np.nanmax(np.abs(test_normhist - allaz)),
        #     vmax=np.nanmax(np.abs(test_normhist - allaz)),
        # )

        # divider = make_axes_locatable(ax[2])
        # cax = divider.append_axes("right", size="5%", pad=0.05)
        # fig.colorbar(im3, cax=cax, orientation="vertical")
        # print(np.nanmax(test_normhist - allaz))
        # print(np.nanmin(test_normhist - allaz))

        # ############################################################################

        # fig, ax = plt.subplots(3, 1, figsize=(12, 12))
        # im1 = ax[0].imshow(
        #     alldaz.T,
        #     cmap="RdBu_r",
        #     vmin=-0.8 * np.nanmax(np.abs(alldaz)),
        #     vmax=0.8 * np.nanmax(np.abs(alldaz)),
        # )
        # divider = make_axes_locatable(ax[0])
        # cax = divider.append_axes("right", size="5%", pad=0.05)
        # fig.colorbar(im1, cax=cax, orientation="vertical")

        # im2 = ax[1].imshow(
        #     test_normhist_daz.T,
        #     cmap="RdBu_r",
        #     vmin=-0.8 * np.nanmax(np.abs(alldaz)),
        #     vmax=0.8 * np.nanmax(np.abs(alldaz)),
        # )

        # divider = make_axes_locatable(ax[1])
        # cax = divider.append_axes("right", size="5%", pad=0.05)
        # fig.colorbar(im2, cax=cax, orientation="vertical")

        # im3 = ax[2].imshow(
        #     (test_normhist_daz - alldaz).T,
        #     cmap="RdBu_r",
        #     vmin=-np.nanmax(np.abs(test_normhist_daz - alldaz)),
        #     vmax=np.nanmax(np.abs(test_normhist_daz - alldaz)),
        # )

        # divider = make_axes_locatable(ax[2])
        # cax = divider.append_axes("right", size="5%", pad=0.05)
        # fig.colorbar(im3, cax=cax, orientation="vertical")
        # print(np.nanmax(test_normhist_daz - alldaz))
        # print(np.nanmin(test_normhist_daz - alldaz))

        # ############################################################################

        # fig, ax = plt.subplots(3, 1, figsize=(12, 12))
        # im1 = ax[0].imshow(
        #     allddaz.T,
        #     cmap="RdBu_r",
        #     vmin=-0.8 * np.nanmax(np.abs(allddaz)),
        #     vmax=0.8 * np.nanmax(np.abs(allddaz)),
        # )
        # divider = make_axes_locatable(ax[0])
        # cax = divider.append_axes("right", size="5%", pad=0.05)
        # fig.colorbar(im1, cax=cax, orientation="vertical")

        # im2 = ax[1].imshow(
        #     test_normhist_ddaz.T,
        #     cmap="RdBu_r",
        #     vmin=-0.8 * np.nanmax(np.abs(allddaz)),
        #     vmax=0.8 * np.nanmax(np.abs(allddaz)),
        # )

        # divider = make_axes_locatable(ax[1])
        # cax = divider.append_axes("right", size="5%", pad=0.05)
        # fig.colorbar(im2, cax=cax, orientation="vertical")

        # im3 = ax[2].imshow(
        #     (test_normhist_ddaz - allddaz).T,
        #     cmap="RdBu_r",
        #     vmin=-np.nanmax(np.abs(test_normhist_ddaz - allddaz)),
        #     vmax=np.nanmax(np.abs(test_normhist_ddaz - allddaz)),
        # )

        # divider = make_axes_locatable(ax[2])
        # cax = divider.append_axes("right", size="5%", pad=0.05)
        # fig.colorbar(im3, cax=cax, orientation="vertical")
        # print(np.nanmax(test_normhist_ddaz - allddaz))
        # print(np.nanmin(test_normhist_ddaz - allddaz))
        # plt.show()
        # sys.exit()
        pcadata = np.where(np.isnan(allaz) == False, allaz, 0)

        pca = PCA(n_components=ncomps)
        pca.fit(pcadata)
        comps = pca.components_
        ampl = np.sum(pcadata[:, None, :] * comps[None, :, :], axis=-1).T

        pcadata_daz = np.where(np.isnan(alldaz) == False, alldaz, 0)

        pca_daz = PCA(n_components=ncomps)
        pca_daz.fit(pcadata_daz)
        comps_daz = pca_daz.components_
        ampl_daz = np.sum(pcadata_daz[:, None, :] * comps_daz[None, :, :], axis=-1).T

        pcadata_ddaz = np.where(np.isnan(alldaz) == False, alldaz, 0)

        pca_ddaz = PCA(n_components=ncomps)
        pca_ddaz.fit(pcadata_ddaz)
        comps_ddaz = pca_ddaz.components_
        ampl_ddaz = np.sum(pcadata_ddaz[:, None, :] * comps_ddaz[None, :, :], axis=-1).T

        # print("PCA time", time.time() - t0, "sec", idx)
        t0 = time.time()
        compsdata[feed, ...] = comps
        ampldata[feed, ...] = ampl

        compsdata_daz[feed, ...] = comps_daz
        ampldata_daz[feed, ...] = ampl_daz

        compsdata_ddaz[feed, ...] = comps_ddaz
        ampldata_ddaz[feed, ...] = ampl_ddaz

    result = {
        "comps_az": compsdata,
        "ampl_az": ampldata,
        "comps_daz": compsdata_daz,
        "ampl_daz": ampldata_daz,
        "comps_ddaz": compsdata_ddaz,
        "ampl_ddaz": ampldata_ddaz,
        # "bindata_az": bindata_az,
        # "bindata_daz": bindata_daz,
        # "bindata_ddaz": bindata_ddaz,
        # "hitdata_az": hitdata_az,
        # "hitdata_daz": hitdata_daz,
        # "hitdata_ddaz": hitdata_ddaz,
        "scanid": f"{id:09}",
        "mean_mjd": average_mjd,
        "mean_el": average_el,
        "mean_ra": average_ra,
        "mean_dec": average_dec,
        "azbins": bins,
        "dazbins": dazbins,
        "ddazbins": ddazbins,
        "freq": freq,
    }
    # print("Making output", time.time() - t0, "sec", idx)
    t0 = time.time()
    return result, idx


# def ddaz_binned_init(q, lock, iterr):

#     ddaz_binned.q = q
#     ddaz_binned.lock = lock
#     ddaz_binned.iterr = iterr


def save_result(result, outfile):
    id = result["scanid"]
    if id not in outfile.keys():
        outfile.create_group(id)

    to_be_saved = [
        "comps_az",
        "ampl_az",
        "comps_daz",
        "ampl_daz",
        "comps_ddaz",
        "ampl_ddaz",
        "azbins",
        "mean_mjd",
    ]

    for key in to_be_saved:
        if key not in outfile[id].keys():
            outfile.create_dataset(id + "/" + key, data=result[key])
        else:
            outfile[id + "/" + key][...] = result[key]


def run_l3gen():
    comm = MPI.COMM_WORLD
    Nproc = comm.Get_size()
    rank = comm.Get_rank()
    status = MPI.Status()
    print("hei", rank)
    tags = ["READY", "DONE", "EXIT", "START"]

    if rank == 0:
        idx = 0
        closed_workers = 0
        while closed_workers < Nproc - 1:
            data = comm.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status)
            source = status.Get_source()
            tag = status.Get_tag()

            if tag == tags.index("READY"):
                if idx < Nfiles:
                    comm.send(idx, dest=source, tag=tags.index("START"))
                    idx += 1
                else:
                    comm.send(None, dest=source, tag=tags.index("EXIT"))

            elif tag == tags.index("DONE"):
                result, index = data
                t0 = time.time()
                with h5py.File(
                    savedir + f"edge_effect_data_v4_{field}.h5", "a"
                ) as outfile:
                    # outfile = h5py.File(savedir + f"edge_effect_data_{field}.h5", "a")
                    save_result(result, outfile)

                    # print("RANK:", source, "DONE with:", index)
                print("Rank", rank, "Saving time", time.time() - t0, "sec")

                if index == 0:
                    with h5py.File(
                        savedir + f"edge_effect_data_v4_{field}.h5", "a"
                    ) as outfile:

                        outfile.create_group("parameters")
                        to_be_saved = ["dazbins", "ddazbins", "freq"]

                        for key in to_be_saved:
                            outfile.create_dataset(
                                "parameters/" + key, data=result[key]
                            )

            elif tag == tags.index("EXIT"):
                closed_workers += 1

    else:

        while True:

            comm.send(None, dest=0, tag=tags.index("READY"))
            task = comm.recv(source=0, tag=MPI.ANY_SOURCE, status=status)
            tag = status.Get_tag()
            if tag == tags.index("START"):
                t0 = time.time()
                result, idx = ddaz_binned(task)
                print("Rank", rank, "Computation time", time.time() - t0, "sec")
                comm.send((result, idx), dest=0, tag=tags.index("DONE"))

            elif tag == tags.index("EXIT"):
                print("RANK:", rank, "Finished")
                break

        comm.send(None, dest=0, tag=tags.index("EXIT"))

    MPI.Finalize()


run_l3gen()
# m = mp.Manager()  # Multiprocess manager used to manage Queue and Lock.
# q = m.Queue()
# lock = m.Lock()
# iterr = mp.Value("i", 0)  # Initializing shared iterator value

# NP = 50

# with mp.Pool(NP, ddaz_binned_init, [q, lock, iterr]) as pool:
#     result = pool.map(ddaz_binned, range(len(allfiles)))

# allresults = {}
# print("Saving to file:")

# with h5py.File(savedir + f"edge_effect_data_S2_{field}_v2.h5", "w") as outfile:
#     for i, element in enumerate(tqdm(result)):
#         # for i, element in enumerate(result):
#         if result:
#             for key, value in element.items():
#                 if i == 0:
#                     if key not in ["freq"]:
#                         allresults[key] = value[None, ...]

#                 else:
#                     if key not in ["freq"]:
#                         allresults[key] = np.concatenate(
#                             (allresults[key], value[None, ...])
#                         )
#         else:
#             continue
#     outfile.create_dataset("scanid", data=allids)
#     for key, value in allresults.items():
#         if key not in ["freq", "scanid"]:
#             if key in "scanid":
#                 continue
#             else:
#                 outfile.create_dataset(key, data=value)
#         else:
#             outfile.create_dataset(key, data=value[0, ...])


# with h5py.File("edge_effect_data_v3.h5", "w") as outfile:
#     outfile.create_dataset("comps_az", data = compsdata)
#     outfile.create_dataset("ampl_az", data = ampldata)
#     outfile.create_dataset("comps_daz", data = compsdata_daz)
#     outfile.create_dataset("ampl_daz", data = ampldata_daz)
#     outfile.create_dataset("comps_ddaz", data = compsdata_ddaz)
#     outfile.create_dataset("ampl_ddaz", data = ampldata_ddaz)
#     outfile.create_dataset("azbins", data = bins)
#     outfile.create_dataset("dazbins", data = dazbins)
#     outfile.create_dataset("ddazbins", data = ddazbins)
#     outfile.create_dataset("freq", data = freq)
#     outfile.create_dataset("scanid", data = allids)
