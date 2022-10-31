import argparse
from map_object import COmap
from scipy import linalg
from scipy import sparse
from scipy.stats import skewtest, skew
import numpy as np
from typing import List, Tuple
import re
import warnings
from tqdm import tqdm

import matplotlib.pyplot as plt

# Ignore RuntimeWarning
warnings.filterwarnings("ignore", category=RuntimeWarning)


class PCA_SubTractor:
    """Class for computing PCA components of COMAP maps"""

    def __init__(self, map: COmap, ncomps: int, clean: bool = True, verbose: bool = False):
        """Initializing class instance

        Args:
            map (Map): COMAP map object to compute PCA compoents of.
            ncomp (int): Number of PCA components to compute/subtract.
            verbose (bool, optional): Boolean specifying whether to run in verbose mode.
            clean (bool, optional): Boolean specifying whether to clean subtract PCA modes.
        """

        self.map = map
        self.ncomps = ncomps        
        self.verbose = verbose
        self.clean = clean

        # List of keys to perform PCA on (only per feed hence remove "map" and "rms")
        self.keys_to_pca = [
            key for key in map.keys if ("map" in key) and ("coadd" not in key)
        ]

    def get_svd_basis(
        self, data: np.ndarray
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Method for finding the frequency-transverse PCA basis for input data.

        Args:
            data (np.ndarray): Data from which to compute PCA basis.
                               Must be of shape (sidebands, frequencies, RA, Dec)

        Returns:
            Tuple[np.ndarray, np.ndarray, np.ndarray]: Tuple containing frequency eigen-vector,
                                                       angular eigen-vector and singular values
                                                       making up PCA decomposition.
        """

        # Get length of array data's dimensions
        nfeed, nsb, nfreq, nra, ndec = data.shape

        # Total number of PCA components
        ncomp = nsb * nfreq

        # Define empty buffers
        all_freqvec = np.zeros((nfeed, ncomp, nsb * nfreq))
        all_angvec = np.zeros((nfeed, ncomp, nra * ndec))
        all_singular_values = np.zeros((nfeed, ncomp))

        # Perform PCA decomposition of data per feed
        for feed in tqdm(range(nfeed)):
            feed_data = data[feed, :]
            feed_data = feed_data.reshape(nsb * nfreq, nra * ndec)
            freq_vec, singular_values, ang_vec = linalg.svd(
                feed_data, full_matrices=False
            )

            # feed_data = feed_data.reshape(nsb * nfreq, nra * ndec)
            # freq_vec, singular_values, ang_vec = sparse.linalg.svds(
            #     feed_data, k=self.ncomps
            # )

            # Fill buffers
            all_freqvec[feed, :, :] = freq_vec.T
            all_angvec[feed, :, :] = ang_vec
            all_singular_values[feed, :] = singular_values

        # Return eigenvectors and singular values of PCA
        return (
            all_freqvec.reshape(nfeed, ncomp, nsb, nfreq),
            all_angvec.reshape(nfeed, ncomp, nra, ndec),
            all_singular_values,
        )

    def coadd_feeds(self):
        # To be implemented
        return NotImplemented

    def normalize_data(self, key: str, norm: str) -> np.ndarray:
        """_summary_

        Args:
            key (str): key to map dataset to perform PCA on
            norm (str): String that specifies what RMS normalization to apply to map prior to SVD. Defaults to "rms".
                        Can take values;
                         - 'approx' - for normalizing map by first PCA mode of RMS-map
                         - 'rms' - for normalizing map by RMS-map
                         - 'var' - for normalizing map by variance-map

        Raises:
            ValueError: norm must be either of "approx", "rms", "var".

        Returns:
            np.ndarray: _description_
        """

        # Making sure normalization is valid
        if norm not in ["approx", "rms", "var"]:
            message = 'Make sure that normalization argument norm is either of the three "approx", "rms", "var".'
            raise ValueError(message)

        # Make key for rms dataset that corresponds to map dataset
        rms_key = re.sub(r"map", "rms", key)

        if norm == "rms":
            # rms normalized PCA
            self.norm_exponent = 1
        elif norm == "var":
            # variance normalized PCA
            self.norm_exponent = 2
        else:
            # rms approximation normalized PCA
            return NotImplemented

        norm_data = self.map[key] / self.map[rms_key] ** self.norm_exponent

        # Remove NaN values from indata
        norm_data = np.where(np.isfinite(norm_data), norm_data, 0)

        return norm_data

    def compute_pca(self, norm: str = "rms"):
        """Method to compute PCA of map datasets

        Args:
            norm (str, optional): String that specifies what RMS normalization to apply to map prior to SVD. Defaults to "rms".
                        Can take values;
                         - 'approx' - for normalizing map by first PCA mode of RMS-map
                         - 'rms' - for normalizing map by RMS-map
                         - 'var' - for normalizing map by variance-map
        """
        if self.verbose:
            print(f"Computing {norm} normalized PCA of:")

        # Compute PCA of all feed-map datasets
        for key in self.keys_to_pca:
            if self.verbose:
                print(" " * 4 + "Dataset: " + f"{key}")

            # Normalize data
            indata = self.normalize_data(key, norm)

            # Compute SVD basis of indata
            freqvec, angvec, singular_values = self.get_svd_basis(indata)

            # Save computed PCA components
            self.map[key + "_pca_freqvec"] = freqvec
            self.map[key + "_pca_angvec"] = angvec
            self.map[key + "_pca_sing_val"] = singular_values
            
            if self.clean:
                fig, ax = plt.subplots(figsize = (10, 8))
                rms_key = re.sub(r"map", "rms", key)

                def gaussian(x):
                    return 1 / np.sqrt(2 * np.pi) * np.exp(- 0.5 * x ** 2)
    
    
                x = np.linspace(-5, 5, 100)
                g = gaussian(x)
                unit_norm_data = np.random.normal(0, 1, int(np.sum(np.isfinite(self.map[key] / self.map[rms_key]))))
                ax.hist((self.map[key] / self.map[rms_key]).flatten(), bins = 100, histtype = "step", density = True, label = "Before subr")
                print(skewtest((self.map[key] / self.map[rms_key]).flatten(), nan_policy="omit"))
                print(skew((self.map[key] / self.map[rms_key]).flatten(), nan_policy="omit"))
                # Clean data
                map_reconstruction = self.reconstruct_modes(key)
                self.map[key] -= map_reconstruction
                print(skewtest((self.map[key] / self.map[rms_key]).flatten(), nan_policy="omit"))
                print(skew((self.map[key] / self.map[rms_key]).flatten(), nan_policy="omit"))
                print(skewtest(unit_norm_data))
                print(skew(unit_norm_data))
                ax.hist((self.map[key] / self.map[rms_key]).flatten(), bins = 100, histtype = "step", density = True, label = "After subr", linestyle = "dashed")
                ax.plot(x, g, label = "N(0, 1)", linestyle = "dashed", c = "r")
                ax.legend()
                ax.set_yscale("log")
                plt.show()
                break

        # Assigning parameter specifying that map object is PCA subtracted
        # and what norm was used
        self.map["is_pca_subtr"] = True
        self.map["pca_norm"] = norm

    def reconstruct_modes(self, key: str):
        """Method to compute PCA reconstructions of computed PCA modes for given dataset key.

        Args:
            key (str): Dataset key for which to clean out self.ncomps PCA modes.
        """
        
        # Getting SVD basis
        freqvec = self.map[key + "_pca_freqvec"][:, :self.ncomps, ...]
        angvec = self.map[key + "_pca_angvec"][:, :self.ncomps, ...]
        singular_values = self.map[key + "_pca_sing_val"][:, :self.ncomps]

        # Input map from key
        inmap = self.map[key]
        
        # Make key for rms dataset that corresponds to map dataset
        rms_key = re.sub(r"map", "rms", key)

        # Input rms map from key
        inrms = self.map[rms_key]

        # Perform outer product from basis vecotrs
        map_reconstruction = angvec[:, :, :, :, None, None] 
        map_reconstruction = map_reconstruction * freqvec[:, :, None, None, :, :]
        map_reconstruction = map_reconstruction * singular_values[:, :, None, None, None, None]
        map_reconstruction = np.sum(map_reconstruction, axis = 1)
        
        # Reshape to get original dimensionality
        map_reconstruction = map_reconstruction.reshape(inmap.shape) 
        
        # Get back original units by undoing normalization
        map_reconstruction *= inrms ** self.norm_exponent

        return map_reconstruction

if __name__ == "__main__":
    mappath = "/mn/stornext/d22/cmbco/comap/protodir/maps/"
    mapname = "co7_map_S2.h5"
    mymap = COmap(mappath + mapname)

    pca_sub = PCA_SubTractor(mymap, 10, True, True)

    pca_sub.compute_pca()
