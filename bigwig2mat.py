
# coding: utf-8

import sys
import h5py
import numpy as np
import pyBigWig as pbw
import re

if __name__ == '__main__':
        dataname = re.split('.big[Ww]ig',sys.argv[1])[0]
	print 'Opening: ' + dataname + '.bigWig'
	
	f = pbw.open(sys.argv[1])

	chrmlentot = np.array([249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895,51304566])
	res = 200.
	chrmlentotbin = np.cumsum(np.ceil(chrmlentot/res))
	chrmlentotbin = np.insert(chrmlentotbin, 0, 0)
	S = int(np.sum(np.ceil(chrmlentot/res)))


	# In[ ]:

	vect = np.zeros(S)
	for chrm in xrange(1,23):
		  interval = f.intervals("chr{0}".format(chrm))
		  print chrm
		  for i in xrange(len(interval)):
		      chrm = int(chrm)
		      pos1 = int(chrmlentotbin[chrm-1]+ np.floor(int(interval[i][0])/res))
		      pos2 = int(chrmlentotbin[chrm-1]+ np.floor(int(interval[i][1])/res))
		      vect[range(pos1,pos2)] = vect[range(pos1,pos2)] + float(interval[i][2])


	# In[ ]:

	f = h5py.File(dataname + '.hdf5', "w")
	f['data'] = vect
	f.close()

