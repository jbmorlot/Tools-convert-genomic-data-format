
import sys
import os
import re
import numpy as np
import h5py

#sys.path.append('~/Documents/Tools/')

def bed2mat(input_file):
	dataname = re.split('.bed',input_file)[0]

	#BED file
	fbed = open(input_file,'r')
	chrmlentot=np.array([ 
249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,59128983,
	63025520,48129895,51304566])
	res = 200.		#Le point est important pour la dicvision qui suit!!!
	chrmlentotbin = np.cumsum(np.ceil(chrmlentot/res))
	chrmlentotbin = np.insert(chrmlentotbin, 0, 0)
	S = np.sum(np.ceil(chrmlentot/res)).astype(np.int32)
	
	vect = np.zeros((S,1))
	k=0
	for line in fbed:
			if re.search(r'chr(\d+)', line):
				#print (line.split('\n')[0]).split('\t')
				values = (line.split('\n')[0]).split('\t')
		
				chrm = int((values[0]).split('chr')[1])
				pos1 = int(chrmlentotbin[chrm-1]+ np.floor(int(values[1])/res))
				pos2 = int(chrmlentotbin[chrm-1]+ np.floor(int(values[2])/res))
				if pos1 != pos2:
									vect[range(pos1,pos2)] = vect[range(pos1,pos2)] + 1
				else:
									vect[range(pos1,pos1+1)] = vect[range(pos1,pos1+1)] + 1

	fbed.close()

	fh5 = h5py.File(dataname + '.hdf5', "w")
	fh5['data'] = vect
	fh5.close()

# DESCRIPTION OF macs2 bdgpeakcall
#NAME
       #bdgpeakcall - Naive call peaks from a single bedGraph track for scores

#SYNOPSIS
       #bdgpeakcall <-i bedGraph> [-c CUTOFF] [-l MIN] [-g MAX] [-o PREFIX]

#DESCRIPTION
       #Call  peaks  from  MACS  pvalue  or  qscore score bedGraph output, with
       #customized settings. Output encodePeak  format  peaks,  combining  peak
       #boundaries, peak summits.

#OPTIONS
       #--version
              #show program's version number and exit

       #-h, --help
              #Show this help message and exit.

       #-i IFILE, --ifile=IFILE
              #MACS pvalue score bedGraph

       #-c CUTOFF, --cutoff=CUTOFF
              #Cutoff  depends on which method you used for score track. If the
              #file contains pvalue scores from MACS2,  score  5  means  pvalue
              #1e-5. DEFAULT: 5

       #-l MINLEN, --min-length=MINLEN
              #minimum  length  of peak, better to set it as d value.  DEFAULT:
              #200

       #-g MAXGAP, --max-gap=MAXGAP
              #maximum gap between significant points in a peak, better to  set
              #it as tag size. DEFAULT: 30

       #-o OPREFIX, --o-prefix=OPREFIX
              #output file prefix, DEFAULT: peak

if __name__ == '__main__':
	
	cutoff = 5

	input_file = sys.argv[1]
	folder,input_file_rel = os.path.split(input_file)
	output_file = folder + '/PEAKS_' + input_file_rel
	
	print '\nComputing MACS...'
	os.system('macs2 bdgpeakcall -c {0}'.format(cutoff) +  ' -i ' + input_file +  ' -o ' + output_file)
	
	print '\nBed2Mat ...'
	bed2mat(output_file)
