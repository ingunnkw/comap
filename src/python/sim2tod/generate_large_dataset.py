import numpy as np
import os
import time
import re
from tqdm import trange
import sys

param_path = "/mn/stornext/d16/cmbco/comap/nils/COMAP_general/src/sim/Parameterfiles_and_runlists/large_dataset/masked/end2end/second_run/cube7/"
params_and_runlists = os.listdir(param_path)
param_files = [param for param in params_and_runlists if ("param" in param)]# and ("11803" in param)]
print(param_files)
print(len(param_files))
#sys.exit()
python_path = "/mn/stornext/d16/cmbco/comap/nils/COMAP_general/src/sim/sim2TOD_v2.py"
l2gen_path = "/mn/stornext/d16/cmbco/comap/nils/comap/dummy/src/f90/l2gen/l2gen"
t0 = time.time()

for i in trange(len(param_files)):
    param = param_path + param_files[i]
    
    param_file  = open(param, "r")
    params      = param_file.read()

    tod_out_path = re.search(r"\nTOD_OUT_DIR\s*=\s*'(\/.*?)'", params)   # Defining regex pattern to search for level1 file with added simulation path.
    tod_out_path = str(tod_out_path.group(1))                            # Extracting path

    runlist_path = re.search(r"\nRUNLIST\s*=\s*'(\/.*?)'", params)  # Defining regex pattern to search for runlist path in parameter file.
    runlist_path = str(runlist_path.group(1))                      # Extracting path

    runlist_file = open(runlist_path, "r")         # Opening 
    runlist = runlist_file.read()

    tod_in = re.search(r"(\/.*?\.\w+)", runlist)
    L1_file_out = str(tod_in.group(1))
    
    param_file.close()
    runlist_file.close()
    
    run_sim2tod = f"python {python_path} -p {param} -n 220"
    run_l2gen   = f"export OMP_NUM_THREADS=32;time mpirun -n 1 {l2gen_path} {param}"
    run_delete  = f"rm {tod_out_path + L1_file_out}"
    
    #print(run_sim2tod)
    #print("\n")        
    os.system(run_sim2tod)
    
    #print(run_l2gen)
    #print("\n")
    os.system(run_l2gen)
    
    #print(run_delete)
    #print("\n")
    os.system(run_delete)
print("Total time: ", time.time() - t0, "sec")