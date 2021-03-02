import numpy as np
from scipy.io import mmwrite

d = 3
h = [0.1,0.1,0.1] # exact values are not relevant but pick the same here
folder = 'mats/'
for n in [50,63,80,101,127,160,202,254]:
	gridm = n
	gridn = n
	coordinates = np.zeros((d,n**d))
	for i in range(n**3):
	    ix = h[0] * (i % gridm); iy = h[1] * ((i % (gridm * gridn)) // gridm); iz = h[2] * (i // (gridm * gridn))
	    coordinates[0,i] = ix
	    coordinates[1,i] = iy
	    coordinates[2,i] = iz
	mmwrite(folder+'neglapl_coord3d_{}_{}_{}'.format(n,n,n),coordinates)