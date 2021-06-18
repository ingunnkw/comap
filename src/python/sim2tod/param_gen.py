import numpy as np 
import re 
import time 


runlist_in_file = "/mn/stornext/d16/cmbco/comap/nils/COMAP_general/src/sim/Parameterfiles_and_runlists/runlist_good.txt"
runlist_out_path = "/mn/stornext/d16/cmbco/comap/nils/COMAP_general/src/sim/Parameterfiles_and_runlists/large_dataset/masked/end2end/second_run/cube7/"

param_in_path = "/mn/stornext/d16/cmbco/comap/nils/COMAP_general/src/sim/Parameterfiles_and_runlists/"
param_in_name = param_in_path + "param_6obsIDs_mix_sim2_co6.txt"

param_out_path = "/mn/stornext/d16/cmbco/comap/nils/COMAP_general/src/sim/Parameterfiles_and_runlists/large_dataset/masked/end2end/second_run/cube7/"
param_out_raw_name = "param_e2e_cube7_"
param_out_list = []

tod_out_path = "/mn/stornext/d16/cmbco/comap/nils/COMAP_general/data/level1/end2end/cube7"

runlist_out_path = param_out_path 
runlist_out_raw_name = "runlist_"
runlist_out_list = []

runlist_file = open(runlist_in_file, "r")         # Opening runlist file
runlist = runlist_file.read()
tod_in_list = re.findall(r"\/.*?\.\w+", runlist)    # Regex pattern to extract all L1 files to open from runlist.

obsIDs = re.findall(r"\s\d{6}\s", runlist)         # Regex pattern to find all scanIDs in runlist
obsIDs = [num.strip() for num in obsIDs]
nobsIDs = len(obsIDs)                  # Number of scanIDs in runlist
runlist_file.close()

for ID in obsIDs:
    param_out_name = param_out_raw_name + ID + "_co6.txt"
    param_out_name = param_out_path + param_out_name
    param_out_list.append(param_out_name)

#print(param_out_list)
#print(len(param_out_list))

    
for ID in range(nobsIDs):
    with open(param_in_name, "r") as infile:
        print(ID)
        with open(param_out_list[ID], "w") as outfile:
            for line in infile:
                if "RUNLIST" in line:
                    runlist_out_name = runlist_out_path + f"runlist_{obsIDs[ID]}_co6.txt"
                    #print(f"RUNLIST                      = '{runlist_out_name}'")
                    outfile.write(f"RUNLIST                      = '{runlist_out_name}'\n")
                elif "TOD_OUT_DIR" in line:
                    tod_out_name = "TOD_OUT_DIR                 =   '/mn/stornext/d16/cmbco/comap/nils/COMAP_general/data/level1/end2end/cube7'"
                    outfile.write(tod_out_name + "\n")
                    #print(tod_out_name)
                elif "LEVEL1_DIR" in line:
                    L1_out_name = "LEVEL1_DIR                   = '/mn/stornext/d16/cmbco/comap/nils/COMAP_general/data/level1/end2end/cube7'"
                    outfile.write(L1_out_name + "\n")
                elif "LEVEL2_DIR" in line:
                    L2_out_name = "LEVEL2_DIR                   = '/mn/stornext/d16/cmbco/comap/nils/COMAP_general/data/level2/Ka/end2end/cube7'"
                    outfile.write(L2_out_name + "\n")
                elif "MAP_DIR" in line:
                    map_out_name = "MAP_DIR                   = '/mn/stornext/d16/cmbco/comap/nils/COMAP_general/data/maps/end2end/cube7/'"
                    outfile.write(map_out_name + "\n")
                elif "OBSID_MAPS" in line:
                    map_out_name = "OBSID_MAPS                   = .false."
                    outfile.write(map_out_name + "\n")
                elif "DATACUBE" in line and "OUT" not in line:
                    map_out_name = "DATACUBE                    =   '/mn/stornext/d16/cmbco/comap/nils/COMAP_general/cubes/cube_gauss_e2e.npy'"
                    outfile.write(map_out_name + "\n")
                elif "DATACUBE_OUT" in line:
                    map_out_name = "DATACUBE_OUT                =   '/mn/stornext/d16/cmbco/comap/nils/COMAP_general/data/maps/end2end/cubes/"
                    outfile.write(map_out_name + "\n")
                
                else:
                    #print(line)
                    outfile.write(line)

        outfile.close()
    infile.close()

for ID in obsIDs:
    runlist_out_name = runlist_out_raw_name + ID + "_co6.txt"
    runlist_out_name = runlist_out_path + runlist_out_name
    runlist_out_list.append(runlist_out_name)

#print(runlist_out_list)
#print(len(runlist_out_list))

    
for i, ID in enumerate(obsIDs):
    with open(runlist_in_file, "r") as infile:
        first_line = infile.readline()
        second_line = infile.readline()
        #print(i)
        with open(runlist_out_list[i], "w") as outfile:
            outfile.write(first_line)
            outfile.write("co6   1\n")
            for line in infile:
                  
                if ID in line[:11]:
                    if "2020-03" in line:
                        new_line = line.split("/")
                        print(new_line)
                        new_line = new_line[0] + "/" + new_line[1] + "/" + new_line[2] 
                        print(new_line)
                        outfile.write(new_line)
                    else:
                        outfile.write(line)
                else:
                    continue
        outfile.close()
    infile.close()