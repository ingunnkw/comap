import matplotlib.pyplot as plt
import numpy as np

import healpy as hp
from reproject import reproject_from_healpix, reproject_to_healpix
from astropy import wcs
from astropy.io import fits
from scipy.interpolate import interp1d

import tqdm

def convolve(pixind):
    theta, phi = hp.pix2ang(nside, pixind)
    r = hp.rotator.Rotator(rot=(0, theta, phi), deg=False, eulertype='X')
    beam2 = r.rotate_map_pixel(beam)
    return sum(beam2*pickup)/sum(beam2)

data = np.loadtxt('spherical_grid_fwd_hemisphere.grd', skiprows=28)

ktype = 1
nset, icomp, ncomp, igrid = 1, 2, 2, 5

xs, ys, xe, ye = -0.9e2, -0.9e2, 0.9e2, 0.9e2


nx, ny, klimit = 2001, 2001, 0 

dx = (xe - xs)/(nx-1)
dy = (ye - ys)/(ny-1)
xcen, ycen = 0,0
i_s=1
i_n=nx
F1_r, F1_i, F2_r, F2_i = data.T

Gmap = F1_r**2 + F1_i**2 + F2_r**2 + F2_i**2
x = np.arange(xs, xe+dx, dx)
y = np.arange(ys, ye+dy, dy)
X, Y = np.meshgrid(x,y)

# reproject Gmap to a healpix grid.


#reproject_to_healpix takes as input_data a type containing the ndarray and a
# WCS object, which I will need to create myself.

# The default projection is zenithal equidistant, or 'ARC'
# X = theta\cos \phi and Y = \theta \sin \phi

nside = 64

w = wcs.WCS(naxis=2)
w.wcs.crpix = [1000.5,1000.5]
w.wcs.cdelt = np.array([dx, dy])
w.wcs.crval = [0, 90]
w.wcs.ctype = ['RA---ARC', 'DEC--ARC']
w.wcs.cunit = ['deg', 'deg']
header = w.to_header()

ax = plt.subplot(1,1,1, projection=wcs.WCS(header))
ax.imshow(np.log10(Gmap.reshape(X.shape)))
ax.grid()

Gmap = Gmap.reshape(X.shape)

rbins = np.linspace(0, 90, 26)
r = (rbins[1:] + rbins[:-1])/2
f = np.zeros(len(rbins)-1)
R = np.sqrt(X**2+Y**2)
for i in range(len(rbins)-1):
    inds = (R > rbins[i]) & (R <= rbins[i+1])
    f[i] = Gmap[inds].mean()
plt.figure()
plt.plot(r, f)
plt.xlabel(r'$\theta$')
plt.ylabel(r'$b(\theta)$')
plt.yscale('log')

beam, footprint = reproject_to_healpix((Gmap, wcs.WCS(header)), 'C', nside=nside)

beam[~np.isfinite(beam)] = 0
beam[beam == hp.UNSEEN] = 0
hp.orthview(np.log10(beam))
#plt.show()

d1 = np.loadtxt('horizon_hwt.txt')
az, el = d1.T

pickup = np.zeros(hp.nside2npix(nside))
lon, lat = hp.pix2ang(nside, np.arange(len(pickup)), lonlat=True)
f = interp1d(az, el, fill_value=(el[0], el[-1]), bounds_error=False,
        kind='nearest')
pickup[lat < f(lon)] = 1
hp.orthview(pickup)

i = 0

conv = np.zeros_like(beam)
plt.close('all')
from multiprocessing import Pool
import os
os.environ['OMP_NUM_THREADS'] = '1'
pool = Pool(processes=8)
x = [pool.apply_async(convolve, args=[i]) for i in range(hp.nside2npix(nside))]
for i in tqdm.tqdm(range(len(x))):
    conv[i] = x[i].get()

#for pixind in tqdm.tqdm(np.arange(hp.nside2npix(nside))):
#    conv[pixind] = convolve(pixind)
hp.orthview(conv)
#plt.show()

hp.write_map('convolution_test.fits', conv/conv.max(), overwrite=True)
hp.mollview(pickup, min=0, max=1, title='Ground profile')
hp.graticule()
plt.savefig('pickup.png')
hp.mollview(beam/beam.max(), min=0.01, max=1, norm='log',  title='Beam')
hp.graticule()
plt.savefig('beam_log.png')
hp.mollview(beam/beam.max(), norm='hist',  title='Beam')
hp.graticule()
plt.savefig('beam_hist.png')
hp.mollview(conv/conv.max(), min=0.01, max=1, norm='log',  title='Convolved map')
hp.graticule()
plt.savefig('conv_log.png')
hp.mollview(conv/conv.max(), norm='hist',  title='Convolved map')
hp.graticule()
plt.savefig('conv_hist.png')

plt.figure()
hp.mollview(pickup, min=0, max=1, title='Ground profile', sub=131)
hp.mollview(beam/beam.max(), norm='log', min=0.01, max=1, title='Beam', sub=132)
hp.mollview(conv/conv.max(), norm='log', min=0.01, max=1, title='Convolved map', sub=133)
plt.savefig('all.png')

#plt.show()
