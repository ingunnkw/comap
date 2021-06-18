import numpy as np
import h5py
import tools

class MapCosmo():
    def __init__(self, mappath, feed=None, jk=None, split=None):
        self.feed = feed
        self.interpret_mapname(mappath)
        with h5py.File(mappath, mode="r") as my_file:
            self.x = np.array(my_file['x'][:])
            self.y = np.array(my_file['y'][:])
            if jk is not None:
                if feed is not None:
                    self.map = np.array(my_file['jackknives/map_' + jk][split,feed-1])
                    self.rms = np.array(my_file['jackknives/rms_' + jk][split,feed-1])
                    
                else:
                    self.map = np.array(my_file['jackknives/map_' + jk][split])
                    self.rms = np.array(my_file['jackknives/rms_' + jk][split])
            else:
                if feed is not None:
                    self.map = np.array(my_file['map'][feed-1])
                    self.rms = np.array(my_file['rms'][feed-1])
                    
                else:
                    try: 
                        self.map = np.array(my_file['map_coadd'][:])
                        self.rms = np.array(my_file['rms_coadd'][:])
                    except:
                        self.map = np.array(my_file['map_beam'][:])
                        self.rms = np.array(my_file['rms_beam'][:])
                
        h = 0.7
        deg2mpc = 76.22 / h  # at redshift 2.9
        dz2mpc = 699.62 / h # redshift 2.4 to 3.4
        K2muK = 1e6
        z_mid = 2.9
        dnu = 32.2e-3  # GHz
        nu_rest = 115  # GHz
        dz = (1 + z_mid) ** 2 * dnu / nu_rest  # conversion 
        n_f = 256  # 64 * 4
        redshift = np.linspace(z_mid - n_f/2*dz, z_mid + n_f/2*dz, n_f + 1)


        self.map = self.map.transpose(3, 2, 0, 1) * K2muK
        self.rms = self.rms.transpose(3, 2, 0, 1) * K2muK
                
        sh = self.map.shape
        self.map = self.map.reshape((sh[0], sh[1], sh[2] * sh[3])) 
        self.rms = self.rms.reshape((sh[0], sh[1], sh[2] * sh[3]))
        self.mask = np.zeros_like(self.rms)
        self.mask[(self.rms != 0.0)] = 1.0
        where = (self.mask == 1.0)
        self.w = np.zeros_like(self.rms)
        self.w[where] = 1 / self.rms[where] ** 2
        meandec = np.mean(self.y)
        self.x = self.x * deg2mpc * np.cos(meandec * np.pi / 180)
        self.y = self.y * deg2mpc 
        self.z = tools.edge2cent(redshift * dz2mpc)
        
        self.nz = len(self.z)
        self.nx = len(self.x)
        self.ny = len(self.y)
        
        self.dx = self.x[1] - self.x[0]
        self.dy = self.y[1] - self.y[0]
        self.dz = self.z[1] - self.z[0]
        
        self.voxel_volume = self.dx * self.dy * self.dz  # voxel volume in (Mpc/h)^3

    def interpret_mapname(self, mappath):
        self.mappath = mappath
        mapname = mappath.rpartition('/')[-1]
        mapname = ''.join(mapname.rpartition('.')[:-2])

        parts = mapname.rpartition('_')
        try:
            self.field = parts[0]
            self.map_string = ''.join(parts[2:])
            if not self.field == '':
                self.save_string = '_' + self.field + '_' + self.map_string
            else:
                self.save_string = '_' + self.map_string
        except:
            print('Unable to find field or map_string')
            self.field = ''
            self.map_string = ''
            self.save_string = ''
        
        if self.feed is not None:
            self.save_string = self.save_string + '_%02i' % self.feed
