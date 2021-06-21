import numpy as np 
from scipy.sparse import csc_matrix, csr_matrix, identity, diags
import sparse_dot_mkl
import scipy.sparse.linalg as linalg
from scipy import signal
import matplotlib.pyplot as plt 
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable
import multiprocessing as multiproc
import time
import sys
import h5py
import WCS
import copy
import os
import re
import argparse
import warnings
import random

warnings.filterwarnings("ignore", category=RuntimeWarning) #ignore warnings caused by weights cut-off

class Destriper():
    def __init__(self, eps = 0, param_file = None, obsID_map = False):
        """Initializing needed variables for making an instance of the Destriper class

        Parameters
        ----------
        eps : int, optional
            Small number to add to diagonal of matrizes to enable their 'inversion'. 
            Should be as close to 0 as possible. By default 0
        param_file : string or None, optional
            COmplete path to parameter file to read. 
            If None, it reads the parameter file from command-line user input.
            By default None
        obsID_map : bool, optional
            Whether to perform mapmaking per obsID, by default False
        """
        self.obsID_map  = obsID_map 
        self.eps        = eps
        self.param_file = param_file 
        
        if param_file != None:
            self.read_paramfile()   # Reading in parameters from file.
        else:
            self.input()            # Reading user input.
            self.read_paramfile()   # Reading in parameters from file.
        print(self.scheme)
        if self.scheme not in ["destriper", "weighted", "avg", "baselines"]:
            print("Please provide one of the allowed mapmaker schemes: 'destriper', 'weighted', 'avg' or 'baselines'.")
            sys.exit()
        
        self.counter = 0                # Defining and initializing
        self.counter2 = 0               # counters.
        
        """Defining needed constants:"""
        
        self.Nside = 120                # Number of pixels along RA/Dec
        self.Npix   = self.Nside ** 2   # Total number of pixels in the image
        
        self.dpix   = 2.0 / 60.0        # Pixel resolution in degrees (2' = 2/60 deg)
        
        self.cube_filename = "/mn/stornext/d16/cmbco/comap/protodir/cube_real.npy"        
        
        self.get_xy()       # Precompute RA/Dec pixel grid

    def input(self):
        """
        Function parsing the command line input.
        """
        parser = argparse.ArgumentParser()
        parser.add_argument("-p", "--param", type = str, default = None,
                            help = """Full path and name to input parameter file.""")
        
        args = parser.parse_args()

        if args.param == None:
            message = """No input parameter file given, please provide an input parameter file!"""
            raise NameError(message)
        else:
            self.param_file     = args.param
        
    def read_paramfile(self):
        """Function reading the parameter file provided by the command line
           argument or passed to class in __init__.
        """
        param_file  = open(self.param_file, "r")    
        params      = param_file.read()
        
        runlist_path = re.search(r"\nRUNLIST\s*=\s*'(\/.*?)'", params)  # Regex pattern to search for runlist path in parameter file.
        self.runlist_path = str(runlist_path.group(1))                  # Extracting path

        runlist_file = open(self.runlist_path, "r")         
        runlist = runlist_file.read()

        patch_name = re.search(r"\s([a-zA-Z0-9]+)\s", runlist)  # Regex pattern to extract patch name of observations.                                                            
        self.patch_name = str(patch_name.group(1))              

        infile_path = re.search(r"\nLEVEL2_DIR\s*=\s*'(\/.*?)'", params)            # Regex pattern to search for level 2 input file path.
        self.infile_path = str(infile_path.group(1)) + "/" + self.patch_name + "/"  # Generating correct form of infile path                
    
        outfile_path = re.search(r"\nMAP_DIR\s*=\s*'(\/.*?)'", params) # Regex pattern to search for directory where to put the maps.
        self.outfile_path = str(outfile_path.group(1))                          
        
        mapname = re.search(r"\nMAP_NAME\s*=\s*'([0-9A-Za-z\_]*)'", params)   # Defining regex pattern to search for name of output map file.
        self.map_name = str(mapname.group(1))                    

        
        scanIDs = re.findall(r"\s\d{8}\s", runlist)         # Regex pattern to find all scanIDs in runlist
        self.scanIDs = [num.strip() for num in scanIDs]
        self.nscanIDs = len(self.scanIDs)                   # Number of scanIDs in runlist
        
        obsIDs = re.findall(r"\s\d{6}\s", runlist)         # Regex pattern to find all obsIDs in runlist
        self.obsIDs = [num.strip() for num in obsIDs]
        self.nobsIDs = len(self.obsIDs)                    # Number of obsIDs in runlist
        
        patch_def_path = re.search(r"\nPATCH_DEFINITION_FILE\s*=\s*'(\/.*?)'", params)  # Regex pattern to search for patch definition file.
        self.patch_def_path = str(patch_def_path.group(1))

        patch_def_file = open(self.patch_def_path, "r")             # Opening patch definition file
        patch_def = patch_def_file.read()

        fieldcent   = re.search(rf"{self.patch_name}\s*([0-9.]+)\s*([0-9.]+)", patch_def) # Regex pattern to search for patch center
        self.fieldcent = np.array([float(fieldcent.group(1)), float(fieldcent.group(2))])               

        highpass_nu = re.search(r"\nNUCUT_HIGH\s*=\s*([0-9.]+)", params)    # Regex pattern to search for reading cut-off nu for hp-filter of TOD.
        self.highpass_nu = float(highpass_nu.group(1))                          
        
        basefreq = re.search(r"\nBASELINE_NU\s*=\s*([0-9.]+)", params)    # Regex pattern to search for destriper baseline frequency.
        self.basefreq = float(basefreq.group(1))                          
        self.baseline_time = 1 / self.basefreq                            # Baseline length in time.
        
        scheme = re.search(r"\nSCHEME\s*=\s*'(\w+)'", params)  # Regex pattern to search for mapmaking scheme.
        self.scheme = str(scheme.group(1))
        
        if self.scheme == "baselines":
            self.baseline_only = True
        else:
            self.baseline_only = False

        Nproc = re.search(r"\nN_FREQ_PROCESS\s*=\s*(\d+)", params)  # Regex pattern to search for number of processes to run in parallel.
        self.Nproc = int(Nproc.group(1))
        
        perform_split = re.search(r"\nUSE_ACCEPT\s*=\s*(.true.|.false.)", params)  # Regex pattern to search for whether to define split batches. 
        self.perform_split = perform_split.group(1)
        
        if self.perform_split == ".true.":
            self.perform_split = True 
        else:
            self.perform_split = False 
        
        if self.perform_split: 
            # Read in filename for split data file
            accpt_data_path = re.search(r"\nACCEPT_DATA_FOLDER\s*=\s*'(\/.*?)'", params)   # Defining regex pattern to search for path to accept data file.
            accpt_data_path = str(accpt_data_path.group(1))                                

            accpt_id = re.search(r"\nACCEPT_DATA_ID_STRING\s*=\s*'([0-9A-Za-z\_]*)'", params)   # Defining regex pattern to search for accept data ID.
            accpt_id = str(accpt_id.group(1))                                

            split_id = re.search(r"\nJK_DATA_STRING\s*=\s*'(\/.*?)'", params)   # Defining regex pattern to search for name of split file.
            if split_id == None:
                split_id = ""
            else:
                split_id = str(split_id.group(1))                                

            self.acceptfile_name = accpt_data_path + "jk_data_" +  accpt_id + split_id + "_" + self.patch_name + ".h5"  # Complete path to accept file

            # jk_list file from parameter file
            split_def = re.search(r"\nJK_DEF_FILE\s*=\s*'(\/.*?)'", params)   # Defining regex pattern to search for complete path to split definition file.
            self.split_def = str(split_def.group(1))                                


            idx_feed = np.arange(19, dtype = int)   # All feed indices
            idx_sb   = np.arange(4,  dtype = int)   # All sideband indices
            idx_freq = np.arange(64, dtype = int)   # All frequency channel indices
 
            feeds, sidebands, frequencies = np.meshgrid(idx_feed, idx_sb, idx_freq, indexing = "ij")

            self.all_idx = np.array([feeds.flatten(), sidebands.flatten(), frequencies.flatten()])  # Array to get all feed, sideband and frequency indices, i.e. "inverse flattening".

        else:
            sidebands = np.arange(4, dtype = int)   # All sideband indices
            freqs     = np.arange(64, dtype = int)  # All frequency channel indices

            sidebands, freqs = np.meshgrid(sidebands, freqs, indexing = "ij")

            self.all_idx = np.array([sidebands.flatten(), freqs.flatten()]) # Array to get all feed, sideband and frequency indices, i.e. "inverse flattening".


        verbose = re.search(r"VERBOSE_PRINT\s*=\s*(.true.|.false.)", params)  # Regex pattern to search for whether print debug information. 
        self.verbose = perform_split.group(1)
        
        runlist_file.close()    
        param_file.close()



        if self.verbose:
            print("Patch def:", self.patch_def_path)
            print("Patch", self.patch_name)
            print("Field center", self.fieldcent)
            print("Runlist:", self.runlist_path)
            print("Infile path:", self.infile_path)
            print("Outfile path:", self.outfile_path)
            print("scan IDs:", self.scanIDs)
            print("obs IDs:", self.obsIDs)
            print("Number of scans:", self.nscanIDs)
            print("Number of obs:", self.nobsIDs)
            print("Map output name:", self.map_name)
            print("Mapmaker scheme:", self.scheme)
            print("Fit baseline only:", self.baseline_only)

            print("Baseline freq:", self.basefreq)
            print("Highpass cut:", self.highpass_nu)
            print("Number of frequency loop processes:", self.Nproc)

        if self.perform_split:
            """If to perform splits the batches for the destriper are defined
            """
            if self.verbose:
                print("Perform split:", self.perform_split)
                print("Accept data_folder: ", accpt_data_path)
                print("Accept ID: ", accpt_id)
                print("Split ID: ", split_id)
                print("Split data file:", self.acceptfile_name)
                print("Split def file:", self.split_def)
            
            self.read_split_def()
            self.read_split_data()
            self.read_feeds_alive()
            
            
    def run(self, feed = 1, sb = 1, freq = 1, freq_idx = None):
        """Function used to run the destriper for a single feed, sideband and frequency.

        Parameters
        ----------
        feed : int, optional
            Index of feed, by default 1
        sb : int, optional
            Index of sideband, by default 1
        freq : int, optional
            Index of frequency channel, by default 1
        freq_idx : int or None, optional
            Index of flattened [20 (feeds), 4 (sb), 64 (freqs)] array., by default None
        """
        if self.verbose:
            t0 = time.time()
        
        """Defining the current feed, sb and freq indices."""
        if freq_idx == None:
            self.feed = feed
            self.sb   = sb 
            self.freq = freq
        else: 
            self.freq_idx = freq_idx
            if self.perform_split:
                feed, sb, freq = self.all_idx[:, freq_idx]
                self.feed, self.sb, self.freq = feed, sb, freq
            else:    
                sb, freq = self.all_idx[:, freq_idx]
                self.sb, self.freq = sb, freq

        if self.perform_split:
            self.get_split_data()   # Import data for current split batch.
            self.initialize_P_and_F()
        else:
            """ Get tod, sigma and mask data from buffer"""
            self.tod    = self.tod_buffer.reshape(self.Nsamp, 4, 64)[:, sb, freq] 
            self.sigma0 = self.sigma0_buffer.reshape(self.Nscans, 4, 64)[:, sb, freq]
            self.mask   = self.mask_buffer.reshape(self.Nscans, 4, 64)[:, sb, freq] 

        self.get_P()        # Get pointing matrix

        self.sigma0_inv = 1 / self.sigma0
        self.sigma0_inv[self.mask == 0] = 0
        self.get_Cn_inv()                   # Get inverse white noise co-variance matrix
       
        self.get_PCP_inv()      # Compute (P^T C^-1 P)^-1 matrix
        
    def read_feeds_alive(self):
        """Reading in all the feeds alive for each scanID in the runlist provided.
        """
        feeds_alive = {}    # Empty dictionary to save filename and feeds that are alive.
        tod_lens  = []      # Possible values of level 2 file feed list, i.e. infile["pixels"].

        for scan in self.split_scans:
            
            filename = f"{self.patch_name}_0{scan:08}.h5"
            infile = h5py.File(self.infile_path + filename, "r")
            
            feeds  = infile["pixels"][()]
            feeds_alive[filename] = feeds - 1       # Translating base 1 to base 0 indexing
            
            tod_shape  = infile["tod"].shape
            Nfeed, Nsb, Nfreq, Nsamp = tod_shape
            Nfeed  -= 1
            tod_lens.append(Nsamp)                  # Saving length of scan

            infile.close()

        self.feeds_alive = feeds_alive

        tod_lens = np.array(tod_lens)           # Array of TOD scan lengths
         
        """ Computing the length of each baselines """

        infile = h5py.File(self.infile_path + filename, "r")
        tod_time   = infile["time"][:2] * 3600 * 24
        infile.close()
       
        dt    = tod_time[1] - tod_time[0]

        Nperbaseline = int(round(self.baseline_time / dt))      # Number of samples per baseline

        Nbaseline = np.floor(tod_lens / Nperbaseline).astype(int)   # Number of baselines

        excess = tod_lens - Nperbaseline * Nbaseline            # Excess baseline lengths

        Nbaseline_tot = np.sum(Nbaseline) + np.sum(excess != 0) # Total number of baselines

        name_buffer   = []      # Buffer for file names
        Nperbaselines = [0]     # Array of samples per baseline

        for i, scan in enumerate(self.split_scans):
            filename = f"{self.patch_name}_0{scan:08}.h5"
            Nbaseline_in_scan = 0
            for j in range(Nbaseline[i]):
                Nperbaselines.append(Nperbaseline)
                name_buffer.append(filename)

            if excess[i] > 0:
                Nperbaselines.append(excess[i])
                name_buffer.append(filename)

        self.N_buffer = np.array(Nperbaselines[1:])
        self.name_buffer = np.array(name_buffer)        

    def read_split_def(self):
        """Function to read split definition file.
        """

        split_names = []
        n_coadd = n_feedmap = n_ctrl = n_test = n_split = 0

        with open(self.split_def, "r") as split_def_file:
            n_split = int(split_def_file.readline().split()[0]) - 1
            split_def_file.readline()
            for line in split_def_file:
                line_element = line.split()
                split_names.append(line_element[0])
                number   = int(line_element[1])
                
                if number == 0:
                    n_coadd += 1
                elif number == 2:
                    n_test += 1
                elif number == 3:
                    n_ctrl += 1
                else:
                    n_feedmap += 1

        self.n_coadd     = n_coadd      # Number of feed coadded maps
        self.n_feedmap   = n_feedmap    # Number of feed maps
        self.n_test      = n_test       # Number of test variables
        self.n_ctrl      = n_ctrl       # Number of control variables
        self.n_split     = n_split      # Number of total splits
        self.split_names = split_names  # List of split names

        split_def_file.close()

    def read_split_data(self):
        """Function reading split data file and thereby defining the split batches for destriper.
        """
        with h5py.File(self.acceptfile_name, "r") as split_file:
            split_list  = split_file["jk_list"][()]
            split_scans = split_file["scan_list"][()] 
        split_file.close()

        self.split_scans = split_scans              # Scan ID list from split data file
        self.all_scans, nfeed, nsb = split_list.shape   # number of scans, feeds and sidebands from split data file
        self.splits = np.zeros((self.n_split, self.all_scans, nfeed, nsb))  # Generating and filling self.split with bits borresponding to binary index of batch.
       
        for i in range(self.n_split):
            x, y, z = np.where(split_list > 2 ** i)
            self.splits[i, x, y, z] = 1
            split_list[x, y, z] -= 2 ** i

        #self.n_split -= self.n_test + self.n_ctrl  # Updating number of splits

        self.Nbatch = 2 ** self.n_split # Number of total possible independent split batches
        
        """Randomly shuffling input scans, theoretically for better cross-linking."""
        tod_lens  = []
        names     = []
        Nscans    = 0
        t = time.time()
        file_list = os.listdir(self.infile_path)
        
        for scan in self.split_scans:
            Nscans += 1
            infile = h5py.File(self.infile_path + f"{self.patch_name}_0{scan:08}.h5", "r")
            tod_shape  = infile["tod"].shape
            Nfeed, Nsb, Nfreq, Nsamp = tod_shape
            Nfeed  -= 1
            tod_lens.append(Nsamp)
            names.append(f"{self.patch_name}_0{scan:08}.h5")
            infile.close()

        names = np.array(names)
        tod_lens = np.array(tod_lens)
        
        scan_idx = np.arange(Nscans, dtype = int)
        np.random.shuffle(scan_idx)
        
        self.names = names[scan_idx]
        self.tod_lens = tod_lens[scan_idx]
        self.splits = self.splits[:, scan_idx, :, :]
        self.split_scans = self.split_scans[scan_idx]
        self.Nscans = Nscans

        if self.verbose:
            print("Defining batches:")

        """Translating to binary to decimal index, defining the batches."""
        exp_of_two = 2 ** np.arange(self.n_split)
        self.batch_def = np.sum(self.splits * exp_of_two[:, np.newaxis, np.newaxis, np.newaxis], axis = 0)

    def get_split_data(self):
        """Function reading in the and formatting all the data for a split batch at current feed, sb and freq.
        """
        currentNames, feed, sb, freq = self.currentNames, self.feed, self.sb, self.freq
        N_in_batch = len(currentNames)  # Number of scans in current split batch
               
        tod_lens  = []  # Lengths of all the scans included in batch

        """Reading the number of elements in each scan to be included in batch."""
        for i in range(N_in_batch):
            filename = currentNames[i]

            infile = h5py.File(self.infile_path + filename, "r")
            tod_shape  = infile["tod"].shape
            Nfeed, Nsb, Nfreq, Nsamp = tod_shape
            tod_lens.append(Nsamp)
            infile.close()

        self.Nfeed = Nfeed
        self.Nsb   = Nsb
        self.Nfreq = Nfreq


        tod_lens = np.array(tod_lens)
        self.tod_lens = tod_lens
        tod_cumlen = np.zeros(N_in_batch + 1).astype(int)   # Cumulative length of the scans
        tod_cumlen[1:] = np.cumsum(tod_lens).astype(int)
        Nsamp_tot = np.sum(tod_lens)                        # Total number of time samples in batch
    
        """ Computing the length of each baseline """
        infile = h5py.File(self.infile_path + currentNames[0], "r")
        tod_time   = infile["time"][:2] * 3600 * 24                 # MJD array in seconds
        infile.close()
       
        dt    = tod_time[1] - tod_time[0]                           # Temporal resolution.

        Nperbaseline = int(round(self.baseline_time / dt))          # Number of time steps per baseline.

        Nbaseline = np.floor(tod_lens / Nperbaseline).astype(int)   # Number of baselines.

        excess = tod_lens - Nperbaseline * Nbaseline                # Excess baseline lengths.

        Nbaseline_tot = np.sum(Nbaseline) + np.sum(excess != 0)     # Total number of baselines, both of full and partial length. 

        """Filling list with number of elements in each baseline."""
        Nperbaselines = [0]
        scan_per_baseline = []
        for i in range(N_in_batch):
            Nbaseline_in_scan = 0
            for j in range(Nbaseline[i]):
                Nperbaselines.append(Nperbaseline)
                scan_per_baseline.append(currentNames[i])
                Nbaseline_in_scan += 1
            if excess[i] > 0:
                Nperbaselines.append(excess[i])
                scan_per_baseline.append(currentNames[i])
                Nbaseline_in_scan += 1

        """Reading in all needed data for destriping the split batch"""
        self.ra  = np.zeros(Nsamp_tot, dtype = np.float32)              # Array of RA coordinate
        self.dec = np.zeros(Nsamp_tot, dtype = np.float32)              # Array of Dec coordinate
        self.tod = np.zeros(Nsamp_tot, dtype = np.float32)              # Time-ordered data
        
        self.sigma0 = np.zeros(N_in_batch, dtype = np.float32)          # White noise level
        self.mask   = np.zeros(N_in_batch, dtype = np.uint8)            # Mask
        
        self.start_stop = np.zeros((2, N_in_batch), dtype = int)        # Cumulative start and stop indices for each scan.
        
        for i in range(N_in_batch):
            name   = currentNames[i]        # Filenames of current batch
            
            alive = self.feeds_alive[name]          # Number of feeds that are alive

            feed_idx  = np.where(alive == feed)[0]  # Generating correct feed index

            if feed_idx.size == 1:  
                start = tod_cumlen[i]
                end   = tod_cumlen[(i + 1)]
                self.start_stop[0, i] = start
                self.start_stop[1, i] = end
               
                infile = h5py.File(self.infile_path + name, "r")

                self.mask[i]        = infile["freqmask"][feed_idx, sb, freq]
                self.sigma0[i]      = infile["sigma0"][feed_idx, sb, freq]
                self.tod[start:end] = infile["tod"][feed_idx, sb, freq, :] 
                self.ra[start:end]  = infile["point_cel"][feed_idx, :, 0] 
                self.dec[start:end] = infile["point_cel"][feed_idx, :, 1] 
                infile.close()
            else:
                continue
        
        self.dt = dt
        self.Nbaseline    = Nbaseline_tot
        self.Nperbaselines = np.array(Nperbaselines)
        self.Nsamp        = Nsamp_tot
        self.Nscans       = N_in_batch
        self.tod_cumlen   = tod_cumlen
        self.scan_per_baseline = np.array(scan_per_baseline) # Names of scans for each baseline

    def get_data(self):
        """Function reading all data from runlist to buffer array. 
          All feeds are listed up sequentially to form a 
          Nfeed times Nsamples long time-stream.  
        """
        tod_lens  = []  # Lengths of all the scans included
        names     = []  # Names of all input level 2 files
        Nscans    = 0   # Number of scans

        """Reading number of time samples, feeds, sb and freqs per scan"""

        files = os.listdir(self.infile_path)    
        random.shuffle(files)                   # Shuffling scans

     
        for filename in files:
            if np.any([(name in filename and len(filename) < 17) for name in self.scanIDs]):
                Nscans += 1
                infile = h5py.File(self.infile_path + filename, "r")
                tod_shape  = infile["tod"].shape
                Nfeed, Nsb, Nfreq, Nsamp = tod_shape
                Nfeed  -= 1
                for i in range(Nfeed):
                    tod_lens.append(Nsamp)
                names.append(filename)
                infile.close()

        self.Nfeed = Nfeed  # Number of feeds
        self.Nsb = Nsb      # Number of sidebands
        self.Nfreq = Nfreq  # Number of frequency channels
        self.names = names

        infile = h5py.File(self.infile_path + names[1], "r")    
        freq       = infile["nu"][0, ...]                        # Frequency array
        freq[0, :] = freq[0, ::-1]
        freq[2, :] = freq[2, ::-1]   
        self.Freq  = freq
        infile.close()
        
        if self.verbose:
            print("Number of scans to load:", Nscans)
        

        tod_lens = np.array(tod_lens)
        self.tod_lens = tod_lens
        tod_cumlen = np.zeros(Nscans * Nfeed + 1).astype(int)
        tod_cumlen[1:] = np.cumsum(tod_lens).astype(int)        # Cumulative length of elements per scan
        Nsamp_tot = np.sum(tod_lens)
    
        self.tod_buffer = np.zeros((Nsamp_tot, Nsb, Nfreq), dtype = np.float32) # TOD buffer array
       
        self.ra = np.zeros(Nsamp_tot, dtype = np.float32)     # RA coordinate buffer array
        self.dec = np.zeros(Nsamp_tot, dtype = np.float32)    # Dec coordinate buffer array

        self.sigma0_buffer = np.zeros((Nscans, Nfeed, Nsb, Nfreq), dtype = np.float32)  # White noise level buffer array
        self.mask_buffer = np.zeros((Nscans, Nfeed, Nsb, Nfreq), dtype = np.uint8)      # Mask buffer array
        
        Nbaseline_tot = 0       # Total number of baselines
        Nperbaselines = [0]     # Number of elements per baseline
        Nperscan      = [0]     # Number of elements per scan


        self.start_stop = np.zeros((2, Nscans), dtype = int)    # Cumulative start and stop indices per scan
        
        for i in range(Nscans):
            infile = h5py.File(self.infile_path + names[i], "r")
    
            if self.verbose:
                print("Loading scan: ", i, ", ", names[i])
    
            freqmask          = infile["freqmask"][:-1, ...]

            freqmask[:, 0, :] = freqmask[:, 0, ::-1]
            freqmask[:, 2, :] = freqmask[:, 2, ::-1] 
            
            tod                = infile["tod"][:-1, ...]
            
            tod[:, 0, :, :]    = tod[:, 0, ::-1, :] 
            tod[:, 2, :, :]    = tod[:, 2, ::-1, :] 

            tod_time   = infile["time"][()] * 3600 * 24
            
            sigma0 = infile["sigma0"][:-1, ...]
            
            sigma0[:, 0, :]    = sigma0[:, 0, ::-1] 
            sigma0[:, 2, :]    = sigma0[:, 2, ::-1] 

            ra        = infile["point_cel"][:-1, :, 0] 
            dec       = infile["point_cel"][:-1, :, 1] 
            
            Nsamp = tod.shape[-1] 
            Nperscan.append(Nsamp)
            dt    = tod_time[1] - tod_time[0]
                        
            Nperbaseline = int(round(self.baseline_time / dt))
            Nbaseline = int(np.floor(Nsamp / Nperbaseline))
            excess = Nsamp - Nperbaseline * Nbaseline

            for j in range(Nfeed):
                Nbaseline_tot += Nbaseline
                
                for k in range(Nbaseline):
                    Nperbaselines.append(Nperbaseline)
            
                if excess > 0:
                    Nbaseline_tot += 1
                    Nperbaselines.append(excess)
                
                self.sigma0_buffer[i, j, ...] = sigma0[j, ...]
                self.mask_buffer[i, j, ...]   = freqmask[j, ...]
            
            start = tod_cumlen[i * Nfeed]
            end   = tod_cumlen[(i + 1) * Nfeed]
            self.start_stop[0, i] = start
            self.start_stop[1, i] = end

            tod = tod.transpose(0, 3, 1, 2)
            tod = tod.reshape(tod.shape[0] * tod.shape[1], Nsb, Nfreq)
            self.tod_buffer[start:end, ...]  = tod
            self.ra[start:end]   = ra.flatten()
            self.dec[start:end]  = dec.flatten()
            
            infile.close()

        self.dt = dt
        self.sigma0_buffer       = self.sigma0_buffer.reshape(Nscans * Nfeed, Nsb, Nfreq)
        self.mask_buffer         = self.mask_buffer.reshape(Nscans * Nfeed, Nsb, Nfreq)
               
        self.Nbaseline    = Nbaseline_tot
        self.Nperbaselines = np.array(Nperbaselines)
        self.Nsamp        = Nsamp_tot
        self.Nscans       = Nscans * Nfeed
        self.tod_cumlen   = tod_cumlen

        self.tod_buffer          = self.tod_buffer.reshape(   self.Nsamp,  Nsb * Nfreq) 
        self.sigma0_buffer       = self.sigma0_buffer.reshape(self.Nscans, Nsb * Nfreq)
        self.mask_buffer         = self.mask_buffer.reshape(  self.Nscans, Nsb * Nfreq)

    def initialize_P_and_F(self):
        """Initializing pixel indices for pointing matrix P 
           and baseline template matrix F.
        """
        self.get_px_index() # Get pixel index of pointing information.

        self.get_F()        # Get baseline template matrix
     
    def get_xy(self):
        """Function defining the pixel grid of the target field.
        """
        Nside, dpix, fieldcent = self.Nside, self.dpix, self.fieldcent
        
        x = np.zeros(Nside)
        y = np.zeros(Nside)
        dx = dpix / np.abs(np.cos(np.radians(fieldcent[1])))
        dy = dpix 
        
        """Min values in RA/Dec. directions"""
        if Nside % 2 == 0:
            x_min = fieldcent[0] - dx * Nside / 2.0 
            y_min = fieldcent[1] - dy * Nside / 2.0  
            
        else: 
            x_min = fieldcent[0] - dx * Nside / 2.0 - dx / 2.0
            y_min = fieldcent[1] - dy * Nside / 2.0  - dy / 2.0
            
        """Defining pixel centers"""

        x[0] = x_min + dx / 2
        y[0] = y_min + dy / 2
        
        for i in range(1, Nside):
            x[i] = x[i - 1] + dx
            y[i] = y[i - 1] + dy
        
        self.x, self.y = x, y       # Pixel RA/Dec coordinates
        self.dx, self.dy = dx, dy   # Pixel resolution

        
    def get_px_index(self):
        """Funtion defining the pixel index for a given RA/Dec telescope pointing coordinate.
        """
        Nside, dpix, fieldcent, ra, dec, dx, dy = self.Nside, self.dpix, self.fieldcent, self.ra, self.dec, self.dx, self.dy
        
        ra_min, dec_min = self.x[0], self.y[0]

        self.px = np.round((ra - ra_min) / dx) * Nside + np.round((dec - dec_min) / dy)
        self.px = self.px.astype(int)

    def get_P(self):
        """Function that sets up pointing matrix from telescope pointing information.
           Each time sample that its pixel i at time j is assigned P_ij = 1. 
           If the sample is masked it is assigned P_ij = 0.
           P has shape Nsamp vs Nsamp
        """
        Nsamp, Npix, Nscans, cumlen, mask = self.Nsamp, self.Npix, self.Nscans, self.tod_cumlen, self.mask

        tt = time.time()

        hits = np.ones(Nsamp)       # The 1s inside P
    
        for i in range(Nscans):
            start = cumlen[i]
            end = cumlen[i + 1]
            hits[start:end] = mask[i]   # If masked P_ij = 0

        rows = np.arange(0, Nsamp, 1)   # Row indices for sparse elements corresponding to all time samples
        cols = self.px                  # Col indices for sparse elements corresponding to pixel indices
        
        """If pixels outside of grid they are masked"""

        if np.any(cols < 0):
            hits[cols < 0] = 0    
            cols[cols < 0] = 0   

        if np.any(cols >= Npix):
            hits[cols >= Npix] = 0    
            cols[cols >= Npix] = 0

        self.P = csc_matrix((hits, (rows, cols)), shape = (Nsamp, Npix), dtype = np.uint8)  # Setting up sparse matrices for pointing
        self.PT = csc_matrix(self.P.T, dtype = np.uint8)
        
    def get_F(self):
        """Function setting up the baseline templates F used for destriping.
           F has shape Nbaselines vs Nsamp.
        """
        Nsamp, Nbaseline, Nperbaselines = self.Nsamp, self.Nbaseline, self.Nperbaselines

        Nperbaselines_cum = np.zeros(Nbaseline + 1)
        Nperbaselines_cum = np.cumsum(Nperbaselines) # Cumulative number of elements per baseline
        
        ones = np.ones(Nsamp)
        rows = np.arange(0, Nsamp, 1)
        cols = np.zeros(Nsamp)

        for i in range(Nbaseline):
            start = Nperbaselines_cum[i]
            end = Nperbaselines_cum[i + 1]
            cols[start:end] = np.tile(i, Nperbaselines[i + 1])  # Computing column index of each baseline
        
        self.F = csc_matrix((ones, (rows, cols)), shape = (Nsamp, Nbaseline), dtype = np.uint8) # Setting up sparse matrices for template matrix
        self.FT = csc_matrix(self.F.T, dtype = np.uint8)
    
    def get_Cn_inv(self):
        """Function that sets up inverse white noise co-variance matrix quantifying.
           The co-variance matrix has shape Nsamp vs Nsamp.
        """
        Nsamp, Nscans, cumlen = self.Nsamp, self.Nscans, self.tod_cumlen
        C_n_inv = np.zeros(Nsamp)

        for i in range(Nscans):
            start = cumlen[i]
            end = cumlen[i + 1]
            C_n_inv[start:end] = self.sigma0_inv[i] ** 2

        self.C_n_inv = diags(C_n_inv)   # Setting up diagonal matrix with inverse white noise variance from each scan.
    
    def get_PCP_inv(self):
        """Function computing (P^T C_n^-1 P)^-1, and setting up corresponding sparse matrix.
           Matrix has shape Npixel vs Npixel.
        """        
        PCP_inv = self.PT.dot(self.C_n_inv)        
        PCP_inv = PCP_inv.dot(self.P) + diags(self.eps * np.ones(self.Npix))
           
        self.PCP_inv = diags(1 / PCP_inv.diagonal(), format = "csc", dtype = np.float32)
   
    def get_FT_C_P_PCP(self):
        """Function computing matrix product F^T C_n^-1 P (P^T C_n^-1 P)^-1.
           Matrix has shape Nbaseline vs Npixel.
        """
        FT_C = self.FT.dot(self.C_n_inv)        
        FT_C_P = FT_C.dot(self.P)
        
        self.FT_C_P_PCP = FT_C_P.dot(self.PCP_inv)
        
    def get_PT_C(self):
        """Function computing matrix product P^T C_n^-1.
           Matrix has shape Npixel vs Nsamp.
        """
        self.PT_C = self.PT.dot(self.C_n_inv)
        
    def get_FT_C(self):
        """Computing matrix product F^T C_n^-1.
           Matrix has shape Nbaseline vs Nsamp.
        """
        self.FT_C = self.FT.dot(self.C_n_inv)        
    
    def Ma(self, a):
        """Function computing left-hand-side of destriping mapmaker equation.

        Parameters
        ----------
        a : np.ndarray
            Baseline amplitude vector

        Returns
        -------
        np.ndarray
            Result of matrix product: (F^T C_n^-1 F - F^T C_n^-1 P (P^T C_n^-1 P)^-1 P^T C_n^-1 F) a
        """
        temp0 = self.FT_C_F.dot(a)
        temp1 = self.PT_C_F.dot(a)
        temp2 = self.FT_C_P_PCP.dot(temp1)
        
        self.counter += 1   # Counting number of iterations performed in CG solver.
        
        if np.any(np.isnan(a)) or np.any(np.isinf(a)):
            print("NaN or Inf in template vector a!", " All NaN or Inf:", np.all(np.isnan(a)), np.all(np.isinf(a)), "Freq:", self.sb, self.freq, np.any(np.isnan(temp0)), np.any(np.isnan(temp1)), np.any(np.isnan(temp2)))
            sys.exit()

        return temp0 - temp2
        
    def d(self, x):        
        """Function computing right-hand-side of destriping mapmaker equation.

        Parameters
        ----------
        x : np.ndarray
            Baseline tod array

        Returns
        -------
        np.ndarray
            Result of matrix product: (F^T C_n^-1  - F^T C_n^-1 P (P^T C_n^-1 P)^-1 P^T C_n^-1) d
        """
        temp0 = self.FT_C.dot(x)
        temp1 = self.PT_C.dot(x)
        temp2 = self.FT_C_P_PCP.dot(temp1)
        
        if np.any(np.isnan(x)) or np.any(np.isinf(x)):
            print("NaN or Inf in x!" " All NaN or Inf:", np.all(np.isnan(x)), np.all(np.isinf(x)), "Freq:", self.sb, self.freq, np.any(np.isnan(temp0)), np.any(np.isnan(temp1)), np.any(np.isnan(temp2)))
            sys.exit()

        return temp0 - temp2

    def get_preconditioner(self):
        """Function computing the inverse diagonal of matrix 
           M = (F^T C_n^-1 F - F^T C_n^-1 P (P^T C_n^-1 P)^-1 P^T C_n^-1 F)
           for faster convergence of CG solver.
        """
        Nbaseline = self.Nbaseline
        precon = np.zeros((Nbaseline))
        precon += self.FT_C_F.diagonal()
        
        temp = csr_matrix(self.FT_C_P_PCP)
        
        for i in range(Nbaseline):
            diag_elem = temp.getrow(i).dot(self.PT_C_F.getcol(i)).data
           
            if diag_elem.size > 0:
                precon[i] += diag_elem[0]
        
        precon[precon != 0] = 1 / precon[precon != 0]

        self.preconditioner = diags(precon)

    def get_baselines(self):
        """Function computing the baseline fit only for the time-stream by solving 
           the destriper equation Ma = x, i.e.
           (F^T C_n^-1 F - F^T C_n^-1 P (P^T C_n^-1 P)^-1 P^T C_n^-1 F) a = (F^T C_n^-1  - F^T C_n^-1 P (P^T C_n^-1 P)^-1 P^T C_n^-1) x,
           where a is the baseline amplitude vector and x is the TOD.
        """
        
        self.get_FT_C()         # Compute needed matrix products
        self.get_FT_C_P_PCP()
        self.get_PT_C()
        
        self.FT_C_F = self.FT_C.dot(self.F)
        self.PT_C_F = self.PT_C.dot(self.F)

        Ma = linalg.LinearOperator((self.Nbaseline, self.Nbaseline) , matvec = self.Ma) # Defining linear operator from self.Ma since full matrix M is not representable
        d  = self.d(self.tod)

        self.get_preconditioner() # Compute preconditioner matrix M
        
        self.a, info = linalg.cg(Ma, d, M = self.preconditioner) # Computing best fit baseline amplitudes with conjugate gradient method.
        
        if self.verbose:
            print("CG final count: ", self.counter)
        
        self.counter = 0    # Resetting CG iteration counters.
        self.counter2 = 0

    def get_destriped_map(self):
        """Function computing the complete destriped map for single frequency slice.
        """
        self.get_baselines()    # Compute baseline TOD

        m = self.PCP_inv.dot(self.PT).dot(self.C_n_inv).dot(self.tod - self.F.dot(self.a))  # Computing map from mapmaker equation.
        self.m = m.reshape(self.Nside, self.Nside)  # Reshapeing pixel vector into Nside by Nside image.

    def get_bin_averaged_map(self):
        """Function computing bin averaged map for a single frequency slice. 
           [NB! Needs updating!]
        """
        cut = self.highpass_nu      # Whether or not to hiqhpass filter TOD prior
        if self.highpass_nu > 0:    # to mapmaking
            highpass = True        
        else:
            highpass = False
            
        if highpass:
            self.tod = np.fft.rfft(self.tod)                # Highpass filtering TOD
            freqs   = np.fft.rfftfreq(self.Nsamp, self.dt)
            cut_idx = np.argmin(np.abs(freqs - cut))
            self.tod[0:cut_idx] = 0
            self.tod     = np.fft.irfft(self.tod[1:])
            
        M = self.PT.dot(self.P) + diags(self.eps * np.ones(self.Npix))  # Mapmaker matrix on lhs of mapmaker equation.
        x = self.PT.dot(self.tod)                                       # Rhs of mapmaker equation.
        
        self.m = linalg.spsolve(M, b).reshape(self.Nside, self.Nside)   # Solving mapmaker equation as a linear equation. 
                                                                        # [Needs to be updated to CG solver for better efficiency].
    
    def get_noise_weighted_map(self):
        """Function computing the noise averaged binned map 
           for a single frequency slice.
           [NB! Needs updating!]
        """
        
        cut = self.highpass_nu  # Whether or not to highpass filter TOD prior
        if self.highpass_nu> 0: # to mapmaking.
            highpass = True
        else:
            highpass = False
        
        
        if highpass:
            self.tod = np.fft.rfft(self.tod)                    # Highpass filtering TOD
            freqs   = np.fft.rfftfreq(self.Nsamp, self.dt)
            cut_idx = np.argmin(np.abs(freqs - cut))
            self.tod[0:cut_idx] = 0
            self.tod     = np.fft.irfft(self.tod, n = self.Nsamp)
        
        self.m = self.PCP_inv.dot(self.PT).dot(self.C_n_inv).dot(self.tod).reshape(self.Nside, self.Nside) # Solving mapmaker equation.
    
    def make_map(self):
        """Function computing complete map for a given 
           frequency slice according to chosen mapmaker scheme.
        """
        if self.scheme == "destriper":
            self.get_destriped_map()
        elif self.scheme == "weighted":
            self.get_noise_weighted_map()
        else:
            self.get_bin_averaged_map()
    
    def make_baseline_only(self):
        """Function computing TOD baseline fit only
        """
        self.get_baselines()

    def get_hits(self):
        """Function computing hit map for a given
           frequency slice, i.e. the diagonal of the 
           vector P^T P.
        """
        self.hits = self.PT.dot(self.P).diagonal().reshape(self.Nside, self.Nside)
        
    def get_rms(self):
        """Function computing the rms map for a given 
           frequency slice, i.e the diagonal of 
           sqrt((P^T C_n^-1 P)^-1).
        """
        self.rms = np.sqrt(self.PCP_inv.diagonal().reshape(self.Nside, self.Nside))
        
    def write_map(self, full_map, full_hits, full_rms):
        """Function writing full map, hit map and rms noise map
           to HDF5 file.

        Parameters
        ----------
        full_map : np.ndarray, dtype = np.float32
            Complete sky map in units Kelvin for all frequencies.
        full_hits : np.ndarray, dtype = np.int32
            Complete hit map for all frequencies.
        full_rms : np.ndarray, dtype = np.float32
            Complete rms noise map for all frequencies.
        """
        Nside, dpix, fieldcent, ra, dec, freq = self.Nside, self.dpix, self.fieldcent, self.ra, self.dec, self.Freq
        
        outfile = self.outfile_path + self.patch_name + "_" + self.map_name + ".h5"
        

        x = np.zeros(Nside)                           # Defining map grid
        y = np.zeros(Nside)
        dx = dpix / np.cos(np.radians(fieldcent[1]))
        dy = dpix 
        
        if Nside % 2 == 0:
            x_min = fieldcent[0] - dx * Nside / 2 
            y_min = fieldcent[1] - dy * Nside / 2  
            
        else: 
            x_min = fieldcent[0] - dx * Nside / 2 - dx / 2
            y_min = fieldcent[1] - dy * Nside / 2  - dy / 2
            
        x[0] = x_min + dx / 2
        y[0] = y_min + dy / 2
        
        for i in range(1, Nside):
            x[i] = x[i - 1] + dx
            y[i] = y[i - 1] + dy
        
        # Masking possible NaNs to 0, should only corresponging to hits = 0.
        full_map = np.where(np.isnan(full_map) == False, full_map, 0)
        full_hits = np.where(np.isnan(full_hits) == False, full_hits, 0)
        full_rms = np.where(np.isnan(full_rms) == False, full_rms, 0)
        
        # Writing data to file.
        with h5py.File(outfile, "w") as outfile:   
            outfile.create_dataset("map_coadd",    data = full_map, dtype = "float32")
            outfile.create_dataset("nhit_coadd",   data = full_hits, dtype = "int32")
            outfile.create_dataset("rms_coadd",    data = full_rms, dtype = "float32")
            outfile.create_dataset("x",            data = x)
            outfile.create_dataset("y",            data = y)
            outfile.create_dataset("n_x",          data = Nside)
            outfile.create_dataset("n_y",          data = Nside)
            outfile.create_dataset("patch_center", data = fieldcent)
            outfile.create_dataset("freq",         data = freq)
            
        outfile.close()

    def save_baseline_tod(self, baseline_buffer):
        """Function saving the TOD baseline fit for all
           frequencies and feeds to HDF5 file. The outfile will
           be named as the infile, just with '_temp.h5' ending
           added.

        Parameters
        ----------
        baseline_buffer : np.ndarray, dtype = np.float32
            Array containing baselines for all feeds, sidebands 
            and frequency channels.
        """
        tod_lens = self.tod_lens
        tod_lens = tod_lens[::self.Nfeed]

        outfile_path = self.infile_path + "baselines/"
        
        if self.verbose:
            print("Saveing baselines to:", outfile_path)
    
        if not os.path.exists(outfile_path):    # If outfile path is not defined; make it.
            os.mkdir(outfile_path)
        
        for i in range(len(self.names)):            # Formatting baselines back into (feed, sb, freq, time) 
                                                    # array for each scan and subsequently saving it to HDF5 file.
            start, stop = self.start_stop[:, i]

            baseline      = baseline_buffer[start:stop, :]

            shape = baseline.shape
            baseline = baseline.reshape(shape[0], 4, 64)
            baseline = baseline.reshape(self.Nfeed, tod_lens[i], 4, 64)
            baseline = baseline.transpose(0, 2, 3, 1)
            baseline[:, 0, :, :] = baseline[:, 0, ::-1, :]              # Flipping two of the sidebands to match COMAP level 2 format.
            baseline[:, 2, :, :] = baseline[:, 2, ::-1, :]

            baseline = baseline.astype(np.float32)
            
            new_name = self.names[i].split(".")
            new_name = new_name[0] + "_temp." + new_name[1]
            
            infile = h5py.File(outfile_path + new_name, "w")
            infile.create_dataset("tod_baseline", data = baseline, dtype = "float32")
            infile.close()

    def save_baseline_tod_per_freq(self, sb, freq, baseline_buffer):
        """Function saving TOD baseline fit per sideband 
           and frequency channel to HDF5 file. 

        Parameters
        ----------
        sb : [type]
            [description]
        freq : [type]
            [description]
        baseline_buffer : [type]
            [description]
        """
        tod_lens = self.tod_lens
        tod_lens = tod_lens[::self.Nfeed]

        outfile_path = self.infile_path + "baselines/"
        
        if self.verbose:
            print("Saveing baselines to:", outfile_path)
        
        if not os.path.exists(outfile_path):    # If outfile path is not defined; make it.
            os.mkdir(outfile_path)
        
        for i in range(len(self.names)):                    # Formatting baselines back into (feed, time)
            start, stop = self.start_stop[:, i]             # array for each scan and subsequently saving it to HDF5 file.
            baseline      = baseline_buffer[start:stop]

            baseline = baseline.reshape(self.Nfeed, tod_lens[i])    
         
            baseline = baseline.astype(np.float32)
            
            new_name = self.names[i].split(".")
            new_name = new_name[0] + "_temp." + new_name[1]
            
            outfile = h5py.File(outfile_path + new_name, "a")
            if "tod_baseline" not in outfile.keys():
                outfile.create_dataset("tod_baseline", (18, 4, 64, baseline.shape[1]), dtype = "float32")
            
            data = outfile["tod_baseline"]
            if sb == 0 or sb == 2:
                data[:-1, sb, -(freq + 1), :] = baseline    # Scipping blind feed #19 (in base 0 indexing)
            else:
                data[:-1, sb, freq, :] = baseline           # Scipping blind feed #19 (in base 0 indexing)

            outfile.close()

    def save_baselines_from_buffer(self):
        """Function saving baseline amplitudes
           and number of samples per baseline to
           HDF5 file for compact storage.
        """
        names = self.name_buffer

        #outfile_path = self.infile_path + "baselines/"
        outfile_path = "/mn/stornext/d16/cmbco/comap/nils/COMAP_general/data/level2/Ka/wo_sim/highpass/002Hz/default/large_dataset/masked/baselines/dummy/"

        #outfile_path = self.infile_path + "splittest/"
        #outfile_path = self.infile_path + "splittest2/"
        #outfile_path = self.infile_path + "feed_separated_all_in_one/"
        
        #outfile_path = "/mn/stornext/d16/cmbco/comap/nils/COMAP_general/data/level2/Ka/wo_sim/highpass/002Hz/default/large_dataset/masked/baselines/splittest2/"
        #outfile_path = "/mn/stornext/d16/cmbco/comap/nils/COMAP_general/data/level2/Ka/wo_sim/highpass/002Hz/default/large_dataset/masked/baselines/feed_separated_all_in_one/"
        #outfile_path = "/mn/stornext/d16/cmbco/comap/nils/COMAP_general/data/level2/Ka/wo_sim/highpass/002Hz/default/large_dataset/masked/baselines/fullfield/"
        #outfile_path = "/mn/stornext/d16/cmbco/comap/nils/COMAP_general/data/level2/Ka/wo_sim/highpass/002Hz/default/large_dataset/masked/baselines/fullfield2/"
        #outfile_path = "/mn/stornext/d16/cmbco/comap/nils/COMAP_general/data/level2/Ka/wo_sim/highpass/002Hz/default/large_dataset/masked/baselines/fullfield3/"
        #outfile_path = "/mn/stornext/d16/cmbco/comap/nils/COMAP_general/data/level2/Ka/wo_sim/highpass/002Hz/default/large_dataset/masked/baselines/fullfield4/"
        #outfile_path = "/mn/stornext/d16/cmbco/comap/nils/COMAP_general/data/level2/Ka/wo_sim/highpass/002Hz/default/large_dataset/masked/baselines/wopreconditioner/"
        #outfile_path = "/mn/stornext/d16/cmbco/comap/nils/COMAP_general/data/level2/Ka/sim2/default/large_dataset/masked/co6/splittest3/"

        if self.verbose:
            print("Saveing baselines to:", outfile_path)
        
        if not os.path.exists(outfile_path):            # If outfile path is not defined; make it.
            os.mkdir(outfile_path)

        for i in range(len(self.names)):            # Looping over infile names and saving baselines.
            currentName = self.names[i]
            alive = self.feeds_alive[currentName]
            N_alive = alive.shape[0]                # Number of feeds alive

            amplitude = self.a_buffer[..., self.name_buffer == currentName].astype(np.float32)  # Get amplitudes from current scan
            Nperbaseline = self.N_buffer[self.name_buffer == currentName].astype(np.int32)      # Get number of samples per baselines from current scan

            new_name = currentName.split(".")
            new_name = new_name[0] + "_temp." + new_name[1]

            outfile = h5py.File(outfile_path + new_name, "a")

            if "amplitudes" not in outfile.keys() and "Nperbaseline" not in outfile.keys():
                outfile.create_dataset("amplitudes", data = np.zeros((N_alive, 4, 64, amplitude.shape[-1]), dtype = np.float32), dtype = "float32")
                outfile.create_dataset("Nperbaseline", data = np.zeros((Nperbaseline.shape[0],), dtype = np.int32), dtype = "int32")
            
            data_amplitudes = outfile["amplitudes"]
            data_Nperbaseline = outfile["Nperbaseline"]
            
            for feed in range(20):
                feed_idx = np.where(alive == feed)[0]
                if feed_idx.size == 1:
                    data_amplitudes[feed_idx[0], ...] = amplitude[feed_idx[0], ...]
                    data_Nperbaseline[:] = Nperbaseline
            outfile.close()

    
if __name__ =="__main__":
    #datapath    = "/mn/stornext/d16/cmbco/comap/nils/COMAP_general/data/level2/Ka/sim/dynamicTsys/co6/"
    #paramfile = "/mn/stornext/d16/cmbco/comap/nils/COMAP_general/src/sim/Parameterfiles_and_runlists/param_destriper_test_co6.txt"
    #paramfile = "/mn/stornext/d16/cmbco/comap/nils/COMAP_general/src/sim/Parameterfiles_and_runlists/param_destriper_test_wosim_co6.txt"
    
    #paramfile = "/mn/stornext/d16/cmbco/comap/nils/COMAP_general/src/sim/Parameterfiles_and_runlists/param_noise_weighted_test_co6.txt"
    #paramfile = "/mn/stornext/d16/cmbco/comap/nils/COMAP_general/src/sim/Parameterfiles_and_runlists/param_noise_weighted_test_wosim_co6.txt"

    #paramfile0 = "/mn/stornext/d16/cmbco/comap/nils/COMAP_general/src/sim/Parameterfiles_and_runlists/param_noise_weighted_python_wohighpass_co6.txt"
    #paramfile1 = "/mn/stornext/d16/cmbco/comap/nils/COMAP_general/src/sim/Parameterfiles_and_runlists/param_noise_weighted_wohighpass_wosim_co6.txt"
    
    # -------------------------------------
    # Multiple different baseline lengths:
    # -------------------------------------
    
    #paramfile0 = "/mn/stornext/d16/cmbco/comap/nils/COMAP_general/src/sim/Parameterfiles_and_runlists/param_destriper_001s_co6.txt"
    #paramfile1 = "/mn/stornext/d16/cmbco/comap/nils/COMAP_general/src/sim/Parameterfiles_and_runlists/param_destriper_1s_co6.txt"
    #paramfile2 = "/mn/stornext/d16/cmbco/comap/nils/COMAP_general/src/sim/Parameterfiles_and_runlists/param_destriper_10s_co6.txt"
    #paramfile3 = "/mn/stornext/d16/cmbco/comap/nils/COMAP_general/src/sim/Parameterfiles_and_runlists/param_destriper_25s_co6.txt"
    #paramfile4 = "/mn/stornext/d16/cmbco/comap/nils/COMAP_general/src/sim/Parameterfiles_and_runlists/param_destriper_50s_co6.txt"
    #paramfile5 = "/mn/stornext/d16/cmbco/comap/nils/COMAP_general/src/sim/Parameterfiles_and_runlists/param_destriper_100s_co6.txt"
    
    
    #paramfile01 = "/mn/stornext/d16/cmbco/comap/nils/COMAP_general/src/sim/Parameterfiles_and_runlists/param_destriper_001s_wosim_co6.txt"
    #paramfile1 = "/mn/stornext/d16/cmbco/comap/nils/COMAP_general/src/sim/Parameterfiles_and_runlists/param_destriper_1s_wosim_co6.txt"
    #paramfile2 = "/mn/stornext/d16/cmbco/comap/nils/COMAP_general/src/sim/Parameterfiles_and_runlists/param_destriper_10s_wosim_co6.txt"
    #paramfile3 = "/mn/stornext/d16/cmbco/comap/nils/COMAP_general/src/sim/Parameterfiles_and_runlists/param_destriper_25s_wosim_co6.txt"
    #paramfile4 = "/mn/stornext/d16/cmbco/comap/nils/COMAP_general/src/sim/Parameterfiles_and_runlists/param_destriper_50s_wosim_co6.txt"
    #paramfile5 = "/mn/stornext/d16/cmbco/comap/nils/COMAP_general/src/sim/Parameterfiles_and_runlists/param_destriper_100s_wosim_co6.txt"
    
    # -------------------------------------
    # Coadded obsID destriper maps:
    # -------------------------------------
    

    #paramfile0 = "/mn/stornext/d16/cmbco/comap/nils/COMAP_general/src/sim/Parameterfiles_and_runlists/param_destriper_obsID_10s_co6.txt"
    #paramfile1 = "/mn/stornext/d16/cmbco/comap/nils/COMAP_general/src/sim/Parameterfiles_and_runlists/param_destriper_obsID_10s_wosim_co6.txt"


    #paramfiles = [paramfile0, paramfile1, paramfile2, paramfile3, paramfile4]
    #paramfiles = [paramfile0, paramfile01]
    
    #paramfiles = [paramfile0, paramfile1]  
    #paramfiles  = [paramfile1]
   
    #paramfiles = ["/mn/stornext/d16/cmbco/comap/nils/COMAP_general/src/sim/Parameterfiles_and_runlists/param_destriper_10s_wosim_large_dataset_co6.txt"]
    #paramfiles = ["/mn/stornext/d16/cmbco/comap/nils/COMAP_general/src/sim/Parameterfiles_and_runlists/param_destriper_10s_large_dataset_co6.txt"]
    #paramfiles = ["/mn/stornext/d16/cmbco/comap/nils/COMAP_general/src/sim/Parameterfiles_and_runlists/param_liss_CES_mix_freqmask_co6.txt"]
    paramfiles  = ["/mn/stornext/d16/cmbco/comap/nils/COMAP_general/src/sim/Parameterfiles_and_runlists/param_save_baseline_test_co6.txt"]
    
    eps      = 0
    #freq_idx = [113]
    freq_idx = range(4 * 64)
    N_proc = 48
    
    for pfile in paramfiles:
        t = time.time()
        destr = Destriper(param_file = pfile)
        print("Loading data and initializing pointing:")
        t0 = time.time()
        destr.get_data()
        destr.initialize_P_and_F()
        print("Loading time:", time.time() - t0, "sec")
        
        t0 = time.time()
        print("Looping over frequencies:")

        def dummy(idx):
            print("\n", "Processing frequency number:", idx, "\n")
            t = time.time()
            destr.run(freq_idx = idx)

            destr.make_map()
            print("\n", "Making map: ", time.time() - t, "sec \n")

            destr.get_rms()
            print("\n", "Making rms map: ", time.time() - t, "sec \n")

            destr.get_hits()
            print("\n", "Making hit map: ", time.time() - t, "sec \n")

            maps = np.array([destr.m, destr.rms, destr.hits])
            maps = np.where(np.isnan(maps) == False, maps, 0)
            return np.array([destr.m, destr.rms, destr.hits])
    
        with multiproc.Pool(processes = N_proc) as pool:
            full_map = pool.map(dummy, freq_idx)
        pool.close()
        pool.join()
        print("Finished frequency loop:", time.time() - t0, "sec")
    
        print("Formating output:")
        
        full_map = np.array(full_map)

        full_rms = full_map[:, 1, :, :]
        full_hits = full_map[:, 2, :, :]
        full_map = full_map[:, 0, :, :]

        full_map = full_map.reshape(4, 64, 120, 120)
        full_hits = full_hits.reshape(4, 64, 120, 120)
        full_rms = full_rms.reshape(4, 64, 120, 120)
        
        
        full_map = full_map.transpose(0, 1, 3, 2)
        full_hits = full_hits.transpose(0, 1, 3, 2)
        full_rms = full_rms.transpose(0, 1, 3, 2)
        
        print("Writing to file:")
        destr.write_map(full_map, full_hits, full_rms)    
    

    