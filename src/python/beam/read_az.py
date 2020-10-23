import numpy as np
import matplotlib.pyplot as plt

d1 = np.loadtxt('horizon_hwt.txt')
az, el = d1.T

plt.figure(figsize=(12, 4))
plt.plot(az, el)
plt.xlabel(r'az (deg)')
plt.ylabel(r'el (deg)')
plt.xlim([0, 360])
plt.savefig('el_profile.png', bbox_inches='tight')
plt.show()
