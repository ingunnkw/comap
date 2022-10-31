import argparse
from typing import Dict, Any
import h5py
import numpy as np
import numpy.typing as ntyping
from dataclasses import dataclass, field


@dataclass
class COmap:
    """COMAP map data class"""

    path: str

    def read_map(self):
        """Function for reading map data from file and fill data dictionary of Map class"""

        # Empty data dict
        self._data = {}

        # Open and read file
        with h5py.File(self.path, "r") as infile:
            for key, value in infile.items():
                if isinstance(value, h5py._hl.group.Group):
                    # If value is a group we want to copy the data in that group
                    if key == "multisplits":
                        # For all parent splits
                        for split_key, split_group in value.items():
                            # For all datasets in parent split
                            for data_key, data_value in split_group.items():
                                # Path to dataset
                                complete_key = f"{key}/{split_key}/{data_key}"
                                self._data[complete_key] = data_value[()]
                    else:
                        # TODO: fill in if new groups are implemented in map file later
                        pass
                else:
                    # Copy dataset data to data dictionary
                    self._data[key] = value[()]

        self.keys = self._data.keys()

    def write_map(self, outpath):
        return NotImplemented

    def __getitem__(self, key: str) -> ntyping.ArrayLike:
        """Method for indexing map data as dictionary

        Args:
            key (str): Dataset key, corresponds to HDF5 map data keys

        Returns:
            dataset (ntyping.ArrayLike): Dataset from HDF5 map file
        """

        return self._data[key]

    def __setitem__(self, key: str, value: ntyping.ArrayLike):
        """Method for saving value corresponding to key

        Args:
            key (str): Key to new dataset
            value (ntyping.ArrayLike): New dataset
        """
        # Set new item
        self._data[key] = value
        # Get new keys
        self.keys = self._data.keys()
