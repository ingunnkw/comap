from __future__ import print_function
import pickle
import errno
import h5py
import numpy as np
#import matplotlib.pyplot as plt
import glob
import os
import sys

featname = ["f0", "optical_pointing", "radio_pointing", "ground_scan",
            "circular_scan", "constant_elevation_scan", "ambient_load",
            "stationary", "sky_dip"]


class MetaData: 
    ### class containing all interesting info about a file. 
    def __init__(self, obs_id, field, scan_mode, 
                 scan_mode_bit, time_range,
                 scan_ranges, mean_az, mean_el,
                 el_std, file_path):
        self.obs_id = obs_id
        self.field = field
        self.scan_mode = scan_mode
        self.scan_mode_bit = scan_mode_bit
        self.time_range = time_range
        self.scan_ranges = scan_ranges
        self.mean_az = mean_az
        self.mean_el = mean_el
        self.el_std = el_std
        self.file_path = file_path


def find_scan_ranges(mjd, status, params):
    ## find the ranges of each scan and return list with start and end times
    scan_status = 1
    reaquire_status = 2
    stationary = 0
    starts = []
    ends = []
    #print(status[:])
    was_on_scan = False
    ### if all of one file is stationary, treat as single scan. Typical for ambient load etc. 
    if all(s == 0 for s in status): 
        starts.append(mjd[0])
        ends.append(mjd[-1])
        return np.array((starts, ends)).transpose()
    
    dt = mjd[1] - mjd[0]
    i = 0
    for time in mjd:
        if status[i] == scan_status: 
            is_on_scan = True
        elif status[i] == reaquire_status:
            is_on_scan = False
        elif status[i] == stationary:
            is_on_scan = False
        else:
            is_on_scan = False
            print('Unknown status, time:', time, ' status:', status[i])
            sys.exit()
        if (not was_on_scan) and is_on_scan:
            starts.append(time + params['OFFSET_START']*dt)

        if was_on_scan and (not is_on_scan):
            ends.append(time - params['OFFSET_END']*dt)
        i += 1
        was_on_scan = is_on_scan

    if len(starts) == len(ends) + 1:
        ends.append(mjd[-1]) 
    
    return np.array((starts, ends)).transpose()


def find_mean_values(scan_ranges, t_point, az, el):
    ## Find useful stuf for each scan
    n_scans = len(scan_ranges)
    mean_az = np.zeros(n_scans)
    mean_el = np.zeros(n_scans)
    el_std = np.zeros(n_scans)
    i = 0
    for start, end in scan_ranges:
        start_ind = np.argmax(t_point > start)
        if t_point[-1] > end:
            end_ind = np.argmax(t_point > end) - 1
        else:
            end_ind = len(t_point)
        mean_az[i] = np.mean(az[start_ind:end_ind])
        mean_el[i] = np.mean(el[start_ind:end_ind])
        el_std[i] = np.std(el[start_ind:end_ind])
        i += 1
    return mean_az, mean_el, el_std


def find_file_dict(foldername, params, mem={}, bad=[]):
    ## Makes list of files with all relevant info
    file_dict = {}

    # os.chdir(foldername)
    for file in glob.glob(foldername + "/**/*.h5", recursive=True):
        if file in bad:
            pass
        else:
            try:
                with h5py.File(file, mode="r") as fd:
                    # Get the attributes
                    try:
                        att = fd['comap'].attrs
                        hk = fd['hk']
                        obs_id = att['obsId']
                        new_format = True
                    except KeyError:
                        new_format = False
                        print('Wrong format:', file)
                        obs_id = -1
                        bad.append(file)
                    is_late_enough = obs_id >= params['EARLIEST_OBSID']
                    if (str(obs_id) in mem.keys()) and is_late_enough:
                        file_dict[str(obs_id)] = mem[str(obs_id)]
                    elif new_format and is_late_enough:
                        # print('File: '+file)
                        file_path = file.replace(foldername, '')
                        field = att['source'].decode('utf-8')  # 'jupiter' # ###########
                        print('File: '+file, field)
                        features = att['features']  # 2 #
                        is_one_feat = True  # testing for multiple features
                        scan_mode = None
                        for bit, name in enumerate(featname):
                            if (features & 1<<bit):
                                #print("(f%d) %s" % (bit, featname[bit]))
                                scan_mode = featname[bit]
                                ###### These features are treated as separate fields
                                if scan_mode == "ambient_load":
                                    field = scan_mode
                                    break
                                if scan_mode == "ground_scan":
                                    field = scan_mode
                                    break
                                if scan_mode == "stationary":
                                    field = scan_mode
                                    break
                                if scan_mode == "sky_dip":  ## seems like sky dip is associated with other field
                                    field = scan_mode
                                    break
                                scan_mode_bit = bit
                                
                                ########## Test if multiple features ###########
                                if not is_one_feat:
                                    print('Warning, multiple features')
                                    for bit, name in enumerate(featname):
                                        if (features & 1<<bit):
                                            print("(f%d) %s" % (bit, featname[bit]))
                                is_one_feat = False
                                #break  # We can in principle have multiple features, for now we just choose first
                        if scan_mode is None:
                            scan_mode = "other"
                            field = scan_mode

                        scan_mode_bit = features
                        time = fd['time']
                        time_range = (time[0], time[-1])
                        t_status = fd['hk/time_track']
                        status = fd['hk/lissajous_status']
                        el = fd['point_tel'][0, :, 1]
                        az = fd['point_tel'][0, :, 0]
                        
                        scan_ranges = find_scan_ranges(t_status, status, params)
                        t_point = fd['time_point']
                        #plt.plot(t_point, el)
                        #plt.show()
                        if len(scan_ranges) > 0:
                            mean_az, mean_el, el_std = find_mean_values(
                                scan_ranges, t_point, az, el)
                            metadata = MetaData(
                                obs_id, field, scan_mode,
                                scan_mode_bit, time_range,
                                scan_ranges, mean_az,
                                mean_el, el_std, file_path)
                            file_dict[str(obs_id)] = metadata
                            mem[str(obs_id)] = metadata
            except OSError:
                print('\nUnable to read file:')
                print(file, '\n')
                pass
    return file_dict, mem, bad


# hat tip: https://stackoverflow.com/a/19201448/5238625
def save_obj(obj, name, folder):
    with open(folder + name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)


def load_obj(name, folder):
    with open(folder + name + '.pkl', 'rb') as f:
        return pickle.load(f)


# From Tony Li
def ensure_dir_exists(path):
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise
# param_file = sys.argv[1]
try:
    param_file = sys.argv[1]
except IndexError:
    print('You need to provide param file as command-line argument')
    sys.exit()

params = {}
with open(param_file) as f:
    fr = f.readlines()

    fr = [f[:-1] for f in fr]

    frs = [f.split(" = ") for f in fr]

    for stuff in frs:
        try:
            i, j = stuff
            params[str(i).strip()] = eval(j)
        except (ValueError, SyntaxError) as ex:
            pass

# print(params)

runlist_name = params['RUNLIST_FILE_PATH']
foldername = params['LEVEL1_DIR']
aux_data_path = params['AUX_SAVED_DATA']
memory_file = 'mem'
bad_file = 'bad'

ensure_dir_exists(aux_data_path)
#/home/havard/Documents/COMAP/scan_detect'  # ['lvl1', 'ex']#'2018/10/'

try:
    mem = load_obj(memory_file, aux_data_path)
except FileNotFoundError:
    mem = {}


try:
    bad = load_obj(bad_file, aux_data_path)
except FileNotFoundError:
    bad = []

file_dict, mem, bad = find_file_dict(foldername, params, mem=mem, bad=bad)

save_obj(mem, memory_file, aux_data_path)

save_obj(bad, bad_file, aux_data_path)
#obs_ids = list()#[int(f.obs_id for f in file_list]

# def argsort(seq):  # sort files by obsid
#     return sorted(range(len(seq)), key=seq.__getitem__)

# # indices = argsort(obs_id_list)
sorted_obs_ids = sorted([int(obs_id) for obs_id in file_dict.keys()])

file_list = [file_dict[str(obsid)] for obsid in sorted_obs_ids]

n_files = len(file_list)

print('Total number of files in runlist: ', n_files)

# print(file_list[0].time_range[1])

field_names = ('jupiter', 'TauA', 'shela', 'hetdex', 'patch1',
               'patch2', 'co1', 'co2', 'co3', 'mars', 'fg1',
               'fg2', 'fg3', 'fg4', 'fg5', 'ambient_load',
               'ground_scan', 'stationary', 'sky_dip', 'other')



def write_runlist(file_list, runlist_name):
    ## Writes runlist based on the file list
    out_file = open(runlist_name, "w")
    n_files_used = 0
    field_list = []
    for field in field_names:
        files_in_field = []
        for current_file in file_list:
            if current_file.field == field:
                files_in_field.append(current_file)

        n_files_in_field = len(files_in_field)

        if n_files_in_field > 0:
            n_files_used += len(files_in_field)
            field_list.append(files_in_field)

    n_fields = len(field_list)
    out_file.write("%d \n" % n_fields)
    print('Number of fields observed: ', n_fields)
    for field in field_list:
        out_file.write("%s   %d \n" % (field[0].field, len(field)))
        for current_file in field:
            n_scans = len(current_file.scan_ranges)
            out_file.write("  %06i  %17.10f %17.10f %02i %i %s \n" %
                           (current_file.obs_id, current_file.time_range[0],
                            current_file.time_range[-1], n_scans, current_file.scan_mode_bit,
                            current_file.file_path))
            i = 0
            for scan_range in current_file.scan_ranges:
                out_file.write("     %06i_%02i %17.10f %17.10f %10.6f %10.6f %10.6f 0 0 0 0   \n" %
                               (current_file.obs_id, i + 1,
                                scan_range[0],
                                scan_range[1],
                                current_file.mean_az[i],
                                current_file.mean_el[i],
                                current_file.el_std[i]))
                i += 1
    out_file.close()
    print('Number of files used: ', n_files_used)

write_runlist(file_list, runlist_name)
