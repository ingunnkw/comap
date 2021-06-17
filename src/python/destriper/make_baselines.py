import numpy as np 
from scipy.sparse import csc_matrix, identity, diags
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
from destr import Destriper 
from tqdm import trange
import ctypes 


t = time.time()
destr = Destriper()                 # Defining destriper instance.

if destr.perform_split:             # Generating indices to run over in parallel pool.
    freq_idx = range(19 * 4 * 64)
else:
    freq_idx = range(4 * 64)

t0 = time.time()                    # Initializing timer.

if not destr.perform_split:
    print("Loading data and initializing pointing:")
    destr.get_data()
    destr.initialize_P_and_F()
    print("Loading time:", time.time() - t0, "sec")

t0 = time.time()
print("Looping over frequencies:")

def dummy(idx):
    """Dummy function to call destriper class in parallel.

    Parameters
    ----------
    idx : int
        Index characterizing which feed, sideband and frequency
        channel to run the destriper for.

    Returns
    -------
    None
        Dummy return.
    """

    if destr.perform_split:
        """Run split batch destriper"""               

        feed, sb, freq = destr.all_idx[:, idx]          # Get feed, sideband and frequency channel index from idx.
        feed, sb, freq = int(feed), int(sb), int(freq)
        
        destr.freq_idx = idx                            # Divide and define split batches
        current_batch_def = destr.batch_def[:, feed, sb]
        destr.unique_batches, destr.indices = np.unique(destr.batch_def[:, feed, sb], return_inverse = True)
        destr.N_batch_per_freq = destr.unique_batches.shape[0]
        
        destr.batch_buffer = [[] for i in range(destr.N_batch_per_freq)]
        
        for i in range(destr.N_batch_per_freq):
        
            destr.currentNames = destr.names[current_batch_def == destr.unique_batches[i]]
                        
            destr.run(feed, sb, freq)

            destr.make_baseline_only()
            
            n0 = 20                         # Defining constants for baseline buffer arrays
            n1 = 4
            n2 = 64
            n3 = destr.N_buffer.shape[0]

            prod1 = n2 * n3
            prod2 = prod1 * n1

            for j in range(len(destr.currentNames)):
                alive = destr.feeds_alive[destr.currentNames[j]]
                N_alive = alive.shape[0]                            # Number of feeds alive

                feed_idx = np.where(alive == feed)[0]               # Index corresponding to alive feed

                if feed_idx.size == 1:
                    """If current feed is alive data is saved to buffer"""

                    idx1 = np.where(destr.name_buffer == destr.currentNames[j])[0]          # Indices of curret baseline data
                    idx2 = np.where(destr.scan_per_baseline == destr.currentNames[j])[0]    # Indices of amplitudes from current scan

                    for k in range(len(idx1)):
                        idx = int(prod2 * feed_idx + prod1 * sb + n3 * freq + idx1[k])  # Computing 4D flattend buffer index
                        
                        dummy.buffer[idx] = destr.a[idx2[k]]                            # Save amplitudes a from current batch to 
                                                                                        # shared memory buffer
    else:        
        """Run all in one destriper"""
        destr.run(freq_idx = idx)  

        destr.make_baseline_only()
    
        initem = [destr.sb, destr.freq, destr.baseline_tod]     # Putting output product into process queue
        dummy.q.put(initem)
    
        if not dummy.q.empty():
            """Save baselines in queue as long as it is not empty."""
            with dummy.lock:    # Lock saver so only one parallel process opens the outfiles at a time
                                # to prevent outfile corruption.

                outitem = dummy.q.get() # Get baselines from process queue.
                destr.save_baseline_tod_per_freq(outitem[0], outitem[1], outitem[2])

    with dummy.lock:
        if destr.verbose:
            print("Frequency loop progress: ",  dummy.iterr.value, "/", len(freq_idx), "|", round(dummy.iterr.value / len(freq_idx) * 100, 4), "%")
        dummy.iterr.value += 1

    return None # Returning dummy argument

def dummy_init(q, lock, iterr, buffer):
    """Function initializing the parallel process pool.

    Parameters
    ----------
    q : multiprocessing.Manager().Queue()
        Parallel process queue. Used to share 
        data between processes within the parallel
        process pool.
    lock : multiprocessing.Manager().Lock()
        Lock object used to prevent multiple processes
        from running a given code segment at a time.
    iterr : multiprocessing.Value, "i"
        Process iterator value shared between parallel
        processes within pool. Used to keep track of 
        the parallel destriper progress.
    buffer : multiprocessing.RawArray, ctypes.c_float
        Shared (flat) buffer array to save baseline 
        amplitudes from each parallel process.
    """

    dummy.q = q
    dummy.lock = lock
    dummy.iterr = iterr
    dummy.buffer = buffer


m = multiproc.Manager()         # Multiprocess manager used to manage Queue and Lock.
q = m.Queue()
lock = m.Lock()
iterr = multiproc.Value("i", 0)     # Initializing shared iterator value
buffer = multiproc.RawArray(ctypes.c_float, (20 * 4 * 64 * destr.N_buffer.shape[0]))    

with multiproc.Pool(destr.Nproc, dummy_init, [q, lock, iterr, buffer]) as pool:
    baselines = pool.map(dummy, freq_idx)

pool.close()
pool.join()     # Wait for all processes to finish

print("Finished frequency loop:", time.time() - t0, "sec")

"""Loading shared buffer amplitudes into class object for saving"""
destr.a_buffer = np.frombuffer(buffer, dtype = np.float32).reshape(20, 4, 64, destr.N_buffer.shape[0])

destr.save_baselines_from_buffer()

print("Time: ", time.time() - t0, "sec")
