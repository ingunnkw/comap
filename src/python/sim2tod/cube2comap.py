import numpy as np
import h5py
import matplotlib.pyplot as plt
import matplotlib
import copy
import WCS
import time
import shutil
from tqdm import trange
import sys
import argparse
import os
import re
import ctypes
from scipy import interpolate


class MapMakerLight():
    def __init__(self):
        """
        Initializing MapMakerLight class and defining class attributes.
        """
        self.nside  = 120           # Number of pixels along RA/Dec
        self.dpix   = 2.0 / 60.0    # Pixel resolution in degrees (2' = 2/60 deg)
        self.nbin   = self.nside ** 2   # Total number of pixels in the image
        self.template_path = "/mn/stornext/d16/cmbco/comap/nils/COMAP_general/data/maps/templates/" # Path to map template
        self.input()                

    def input(self):
        """
        Function parsing the command line input.
        """
        parser = argparse.ArgumentParser()

        parser.add_argument("-i", "--infile", type = str, default = None,
                            help = """Full path to input simulation cube file.""")

        parser.add_argument("-o", "--outfile", type = str, default = None,
                            help = """Full path to output map file.""")

        parser.add_argument("-m", "--inmap", type = str, default = None,
                            help = """Full path to input map file.""")

        parser.add_argument("-f", "--field", type = str, default = None,
                            help = """coX field, where X = 2, 6 or 7.""")

        parser.add_argument("-n", "--norm", type = float, default = 1.0,
                            help = """Normalize simulation cube by input value.""")
        
        parser.add_argument("-r", "--rms", action = "store_false",
                            help = """Set simulation map's rms to one.""")

        args = parser.parse_args()

        if args.infile == None:
            message = """No infile file given, please provide an input file!"""
            raise NameError(message)
       
        if args.outfile == None:
            message = """No output file given, please provide an output file!"""
            raise NameError(message)
       
        if args.field == None:
            message = """No science field chosen. Please provide co2, co6 or co7"""
            raise NameError(message)
       
        else:
            self.infile      = args.infile
            self.outfile     = args.outfile
            self.inmap       = args.inmap
            self.field       = args.field
            self.norm        = args.norm
            self.copy_rms    = args.rms

    def run(self):
        """
        Function executing the projection from TOD to map.
        """
    
        print("Processing simulation cube: ")
        self.load_cube()                  
        self.make_map_cube()
        self.write_map()
        
        print("Cube processed!")

    def load_cube(self):
        """
        Read the simulated datacube into memory.
        """
        if ".npy" in self.infile[-4:]
            cube = np.load(self.infile)
            cubeshape = cube.shape

            cube *= 1e-6 * self.norm    # Normalization of cube by input value
            print("Maximum of cube", np.nanmax(cube), " Kelvin")
            cube = cube.reshape(cubeshape[0], cubeshape[1], 4, 1024)  # Flatten the x/y dims, and split the frequency (depth) dim in 4 sidebands.
            cube = cube.reshape(cubeshape[0], cubeshape[1], 4, 64, 16)
            cube = np.mean(cube, axis = 4)                            # Averaging over 16 frequency channels
            cube = cube.transpose(2, 3, 0, 1)

            self.cube = cube

        else:
            with h5py.File(self.infile, "r") as infile:
                self.cube = infile["simulation"][()]
                self.cube_x = infile["x"][()]   # x bin edges
                self.cube_y = infile["y"][()]   # y bin edges

                self.cube_x = 0.5 * (self.cube_x[:-1] + self.cube_x[1:]) # converting to bin centers
                self.cube_y = 0.5 * (self.cube_y[:-1] + self.cube_y[1:])
    
    def interp_cube(self):

        X, Y = np.meshgrid(self.cube_x, self.cube_y)
        
        self.cube_interp = interpolate.interp2d(self.cube_x, self.cube_y, )


    def make_map_cube(self):
        """
        Function mapping the simulated cube to the pixel regime of a map file. 
        """
        if self.inmap != None:
            inmap   = h5py.File(self.inmap, "r")
            
            try:
                inhits    = np.array(inmap["nhit_coadd"])[()]
            except KeyError:
                inhits    = np.array(inmap["nhit_beam"])[()]

            self.nhit   = inhits.copy()
            if self.copy_rms:     
                try:
                    inrms   = np.array(inmap["rms_coadd"])[()]
                except KeyError:
                    inrms   = np.array(inmap["rms_beam"])[()]
                self.rms    = inrms.copy()
            else:
                rms         = np.ones_like(self.nhit)
                rms         = np.where(self.nhit > 0, rms, 0)
                self.rms    = rms.transpose(0, 1, 3, 2) 
    
            self.cube   = self.cube.transpose(0, 1, 3, 2) 
            self.cube   = np.where(self.nhit > 0, self.cube, 0)
            self.map    = self.cube  
            
            inmap.close()
        else:
            self.cube   = self.cube.transpose(0, 1, 3, 2) 
            self.map    = self.cube.copy()  
            self.nhit    = np.ones_like(self.map)
            self.rms    = np.ones_like(self.map)
            

    def copy_mapfile(self):
        """
        Function copying copying the template file to be filled
        with map and hits coadded over feeds.
        """
        template_file = "/mn/stornext/d16/cmbco/comap/nils/COMAP_general/data/maps/templates/"
        if self.field == "co2":
            template_file = template_file + "co2_template_map.h5"
        elif self.field == "co6":
            template_file = template_file + "co6_template_map.h5"
        elif self.field == "co7":
            template_file = template_file + "co7_template_map.h5"

        print("Filling template with cube:", template_file)
        shutil.copyfile(template_file, self.outfile)


    def write_map(self):

        self.copy_mapfile()

        with h5py.File(self.outfile, "r+") as outfile:  # Write new sim-data to file.
            try:
                map_coadd   = outfile["map_coadd"] 
                nhit_coadd  = outfile["nhit_coadd"] 
                rms_coadd   = outfile["rms_coadd"] 
                
            except KeyError:
                map_coadd   = outfile["map_beam"] 
                nhit_coadd  = outfile["nhit_beam"] 
                rms_coadd   = outfile["rms_beam"] 
            
            map_coadd[...]  = self.map
            nhit_coadd[...] = self.nhit
            rms_coadd[...]  = self.rms 

        outfile.close()

if __name__ == "__main__":
    maker = MapMakerLight()
    maker.run()
