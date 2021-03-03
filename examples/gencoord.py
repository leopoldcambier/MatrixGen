import numpy as np
from scipy.io import mmwrite

d = 2
assert(d == 2 or d == 3)
h = [0.1,0.1,0.1] # exact values are not relevant but pick the same here
folder = 'mats/'
for n in [100,200,400,800,1600,3200,6400,12800,25600]:
	gridm = n
	gridn = n
	coordinates = np.zeros((d,n**d))
	for i in range(n**d):
	    ix = h[0] * (i % gridm); iy = h[1] * ((i % (gridm * gridn)) // gridm); iz = h[2] * (i // (gridm * gridn))
	    coordinates[0,i] = ix
	    coordinates[1,i] = iy
	    if d == 3:
	    	coordinates[2,i] = iz
	if d == 3:
		mmwrite(folder+'neglapl_coord3d_{}_{}_{}'.format(n,n,n),coordinates)
		print('written to ' + 'neglapl_coord3d_{}_{}_{}'.format(n,n,n))
	else:
		mmwrite(folder+'neglapl_coord2d_{}_{}'.format(n,n),coordinates)
		print('written to ' + 'neglapl_coord2d_{}_{}'.format(n,n))		
