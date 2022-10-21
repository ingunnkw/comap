import argparse
from map_object import COmap
from scipy import linalg
import numpy as np


class PCA_SubTractor:
    """Class for computing PCA components of COMAP maps"""

    def __init__(self, map: COmap):
        """Initializing class instance

        Args:
            map (Map): COMAP map object to compute PCA compoents of.
        """
        self.map = map

    def get_svd_basis(self, data: np.ndarray):  # , x, y, freq):

        nsb, nfreq, nra, ndec = data.shape

        freq_vec, singular_values, ang_vec = linalg.svd(
            data.reshape(nsb * nfreq, nra * ndec), full_matrices=False
        )

        ncomp = singular_values.size

        return (
            freq_vec.T.reshape(ncomp, nsb, nfreq),
            ang_vec.reshape(ncomp, nra, ndec),
            singular_values,
        )


if __name__ == "__main__":
    mappath = "/mn/stornext/d22/cmbco/comap/protodir/maps/"
    mapname = "co7_map_S2.h5"
    mymap = COmap(mappath + mapname)
    mymap.read_map()
    print(type(mymap.data["map"]))
