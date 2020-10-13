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

hp.orthview(conv, norm='log', rot=(0,90,0), half_sky=True, title='Beam pickup',
        min=0.005, max=0.01)
hp.graticule()
plt.savefig('beam_log.png', bbox_inches='tight')
