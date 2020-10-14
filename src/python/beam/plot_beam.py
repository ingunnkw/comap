import matplotlib.pyplot as plt
import numpy as np
import healpy as hp



conv = hp.read_map('convolution_test.fits')

pix = np.arange(len(conv))
theta, phi = hp.pix2ang(hp.npix2nside(len(conv)), pix)
conv[theta > (60*np.pi/180)] = hp.UNSEEN

hp.orthview(conv, norm='hist', rot=(0,90,0), half_sky=True, title='Beam pickup')
hp.graticule()
plt.savefig('beam_hist.png', bbox_inches='tight')

hp.orthview(conv, norm='log', rot=(0,90,0), half_sky=True, title='Beam pickup')
hp.graticule()
plt.savefig('beam_log.png', bbox_inches='tight')

hp.orthview(conv, rot=(0,90,0), half_sky=True, title='Beam pickup', min=0,
        max=0.01)
hp.graticule()
plt.savefig('beam_lin.png', bbox_inches='tight')


az = phi*180/np.pi
el = (np.pi/2-theta)*180/np.pi

inds = ((az >=308) & (az <= 320) & (el >= 30) & (el <= 55))

l1 = np.linspace(308, 320)
l2 = np.linspace(30, 55)
l3 = np.ones_like(l1)*30
l4 = np.ones_like(l1)*55
l5 = np.ones_like(l1)*308
l6 = np.ones_like(l1)*320

conv = hp.read_map('convolution_test.fits')
plt.close('all')
hp.orthview(conv, norm='hist', rot=(0,90,0), half_sky=True, title='Beam pickup',
        max=0.01)
#hp.mollview(conv, norm='hist')
hp.graticule()
hp.projplot(l1, l3, lonlat=True, color='r')
hp.projplot(l1, l4, lonlat=True, color='r')
hp.projplot(l5, l2, lonlat=True, color='r')
hp.projplot(l6, l2, lonlat=True, color='r')

bins = np.linspace(308, 320, 11)
prof = np.zeros(len(bins) - 1)

for i in range(1, len(bins)):
    inds = ((az>bins[i-1]) & (az <= bins[i]) & (el >= 30) & (el <= 55))
    prof[i-1] = conv[inds].mean()
plt.figure()
az = (bins[:-1] + bins[1:])/2
print(bins)
print(az)
plt.plot(az, prof, 'o-', label='Output profile')

az2, el2 = np.loadtxt('horizon_hwt.txt').T
inds = np.where((az2 >= min(bins)) & (az2 <= max(bins)))
az2, el2 = az2[inds], el2[inds]
el2 = (el2 - el2.mean())*prof.std()/el2.std() + prof.mean()
plt.plot(az2, el2, label='Normalized ground profile model')
plt.xlabel('az (deg)')
plt.legend(loc='best')
plt.savefig('profile_model.png', bbox_inches='tight')
plt.show()
