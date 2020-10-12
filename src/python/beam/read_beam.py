# Using instructions found at
# https://caltech.app.box.com/s/5tef79ghrwymck2tynhgpijyjy0tshnr

import matplotlib.pyplot as plt
import numpy as np

# data is gridded such that
# X = -theta*cos(phi)
# Y = theta*sin(phi)
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

Gmap0 = F1_r**2 + F1_i**2 + F2_r**2 + F2_i**2
Gmap1 = 10*np.log10(Gmap0)
x = np.arange(xs, xe+dx, dx)
y = np.arange(ys, ye+dy, dy)

X, Y = np.meshgrid(x,y)

plt.figure()
plt.pcolormesh(X,Y,Gmap0.reshape(X.shape)/Gmap0.max())
plt.gca().set_aspect('equal')
plt.xlabel('az')
plt.ylabel('el')
plt.colorbar()
plt.savefig('f1.png', bbox_inches='tight')

plt.figure() 
plt.pcolormesh(X,Y,Gmap1.reshape(X.shape), vmin=-12, vmax=68, cmap='winter')
plt.gca().set_aspect('equal')
plt.xlabel('az')
plt.ylabel('el')
plt.colorbar(label='dB')
plt.savefig('f2.png', bbox_inches='tight', dpi=200)
plt.show()
