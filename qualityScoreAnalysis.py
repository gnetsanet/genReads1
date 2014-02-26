import os
import sys
import numpy as np
import copy
import subprocess
import cPickle as pickle
import matplotlib.pyplot as mpl

from misc_mayo	import *

N_QSCORES = 41
RQ = N_QSCORES+1

OFFSET = 33

def main():
	PICKLE_DIR = '/Users/zach/Sites/SVN/QC_pickles/'

	print 'reading input pickles...'
	allRuns = {}
	for n in [m for m in os.listdir(PICKLE_DIR) if m[-2:] == '.p']:
		print n[:-2]
		allRuns[n[:-2]] = pickle.load(open(PICKLE_DIR+n,'rb'))
		#break
	print 'pickles read.\n'

	allMeans = {}
	allP10   = {}
	allP90   = {}
	instrCount = {}
	for k in allRuns.keys():
		for n in allRuns[k]:
			if 'platform' in n:
				plt = n['platform']
			else:
				plt = 'unknown'
			if plt in instrCount:
				instrCount[plt] += 1
			else:
				instrCount[plt] = 1

			if 'read1' in n and 'read2' in n:# and plt == 'Illumina Hi Seq - 2000':
				for r in ['read1','read2']:
					theD = n[r]['Per base sequence quality']
					pos  = theD[0]
					mean = theD[1]
					med  = theD[2]
					p25  = theD[3]
					p75  = theD[4]
					p10  = theD[5]
					p90  = theD[6]
					for i in xrange(len(pos)):
						if '-' in pos[i]:
							splt = pos[i].split('-')
							pmin = int(splt[0])
							pmax = int(splt[1])
						else:
							pmin = int(pos[i])
							pmax = pmin
						#print (pmin,pmax)
						mv   = int(mean[i]+0.5)
						p10v = p10[i]
						p90v = p90[i]
						for j in range(pmin,pmax+1):
							if j in allMeans:
								allMeans[j].append(mv)
								allP10[j].append(p10v)
								allP90[j].append(p90v)
							else:
								allMeans[j] = [mv]
								allP10[j] = [p10v]
								allP90[j] = [p90v]
	for k in instrCount.keys():
		print k,instrCount[k]
	for k in sorted(allMeans.keys()):
		mv  = float(sum(allMeans[k]))/len(allMeans[k])
		m10 = float(sum(allP10[k]))/len(allP10[k])
		m90 = float(sum(allP90[k]))/len(allP90[k])
		#print k, '{:.3f}'.format(mv), '({:.3f}, {:.3f})'.format(m10,m90)
		allMeans[k] = mv
		allP10[k]   = m10
		allP90[k]   = m90


	# compute prior distributions for each cycle
	allPMFs = []
	for i in xrange(1,max(allMeans.keys())+1):
		std = (allP90[i]-allP10[i])/2.6
		mu  = allMeans[i]
		var = std*std
		myPMF = [0. for n in xrange(RQ)]
		for j in xrange(RQ):
			myPMF[j] = np.exp(-(((j-mu)**2)/(2*var)))/(np.sqrt(2*np.pi*var))
		mySum = sum(myPMF)
		for j in xrange(RQ):
			myPMF[j] /= mySum
		allPMFs.append(myPMF)

	mpl.figure(0)
	mpl.plot(range(RQ),allPMFs[0])
	mpl.xlim([0,N_QSCORES])
	#mpl.show()

	FQREADS = 's_2_1_sequence.fastq'
	f = open(FQREADS,'r')
	nReads = int(subprocess.check_output("wc -l "+FQREADS, shell=True).split(' ')[1])/4

	print 'generating q-score transition matrix... (using '+str(nReads)+' reads)'
	totalQ = np.ones([RQ,RQ])
	rRead  = 0
	while True:
		data1 = f.readline()
		data2 = f.readline()
		data3 = f.readline()
		data4 = f.readline()
		if not all([data1,data2,data3,data4]):
			break

		qscore = [ord(n)-OFFSET for n in data4[:-1]]
		for i in range(1,len(qscore)):
			totalQ[qscore[i-1],qscore[i]] += 1
		rRead += 1
		if rRead%10000 == 0:
			print rRead
		#break
	print totalQ
	f.close()

	probQ  = [[0. for m in xrange(RQ)] for n in xrange(RQ)]
	for i in xrange(RQ):
		rowSum = float(np.sum(totalQ[i,:]))
		for j in xrange(RQ):
			probQ[i][j] = totalQ[i][j]/rowSum

	mpl.figure(1)
	Z = np.array(probQ)
	X, Y = np.meshgrid( range(0,len(Z[0])+1), range(0,len(Z)+1) )
	mpl.pcolormesh(X,Y,Z[::-1],vmin=0.,vmax=0.25)
	mpl.axis([0,len(Z[0]),0,len(Z)])
	mpl.yticks(range(0,len(Z)),range(len(Z)-1,-1,-1))
	mpl.xticks(range(0,len(Z[0])),range(0,len(Z[0])))
	mpl.xlabel('next Q')
	mpl.ylabel('current Q')
	mpl.colorbar()
	mpl.show()

	Qscores = range(RQ)
	pickle.dump([allPMFs,probQ,Qscores,OFFSET],open('qScoreStuff.p','wb'))

if __name__ == '__main__':
	main()