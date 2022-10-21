import argparse
import h5py
import numpy as np


class COmap:
    """COMAP map data class"""

    def __init__(self, path: str):
        """Initializing map object

        Args:
            path (str): Full path to map file from which to initialize map object
        """

        self.path = path
        self.read_map()
        self.keys = self.data.keys()

    def read_map(self):
        """Function for reading map data from file and fill data dictionary of Map class"""

        # Empty data dict
        self.data = {}

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
                                self.data[complete_key] = data_value[()]
                    else:
                        # TODO: fill in if new groups are implemented in map file later
                        pass
                else:
                    # Copy dataset data to data dictionary
                    self.data[key] = value[()]
