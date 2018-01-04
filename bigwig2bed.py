
# coding: utf-8

# In[1]:
import sys
import h5py
import numpy as np
import pyBigWig as pbw
import re

if __name__ == '__main__':

	dataname = re.split('.big[Ww]ig',sys.argv[1])[0]

	print 'Opening: ' + dataname + '.bigWig'
	print 'Saving: ' + dataname + '.bed'
	if int(sys.argv[2])==1:
		  print 'Saving: ' + dataname + '.hdf5'

	#BigWig file
	f = pbw.open(sys.argv[1])

	fbed = open(dataname + '.bed','w')

	chrmlentot = np.array([249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895,51304566])
	res = 200
	chrmlentotbin = np.floor(chrmlentot/res)
	chrmlentotbin = np.insert(chrmlentotbin, 0, 0)
	chrmlentotbinCS = np.cumsum(np.floor(chrmlentot/res))
	chrmlentotbinCS = np.insert(chrmlentotbinCS, 0, 0)
	S = np.sum(np.floor(chrmlentot/res))

	vect = np.zeros((S,1))
	for chrm in xrange(1,23):
		interval = f.intervals("chr{0}".format(chrm))
		print chrm
		for i in xrange(len(interval)):
			chrm = int(chrm)
			pos1 = int(chrmlentotbinCS[chrm-1]+ np.floor(int(interval[i][0])/res))
			pos2 = int(chrmlentotbinCS[chrm-1]+ np.floor(int(interval[i][1])/res))
			vect[range(pos1,pos2)] = vect[range(pos1,pos2)] + float(interval[i][2])

		for i in xrange(int(chrmlentotbin[chrm])):
			idx = int(i + chrmlentotbinCS[chrm-1])
			if vect[idx]>0:
				fbed.write('chr{0}\t{1}\t{2}\t{3}\n'.format(chrm,i*res,(i+1)*res,vect[idx]))

	f.close()
	fbed.close()
	#Save HDF5 file if asked
	if int(sys.argv[2])==1:
		fh5 = h5py.File(dataname + '.hdf5', "w")
		fh5['data'] = vect
		fh5.close()
		

