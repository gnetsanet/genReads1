#!/usr/bin/env python
# encoding: utf-8

""" **************************************************

qScoreModel.py

- create q-score model given input FastQ data, maximum q-score, and ASCII q-score offQ

Written by:		Zach Stephens
Date:			March 17, 2014
Contact:		zstephe2@illinois.edu


Usage:

python qScoreModel.py in.fq maxQscore asciiOffset out.p

Example:

python qScoreModel.py input.fq 41 33 outModel.p

************************************************** """

import os
import sys
import numpy as np
import copy
import subprocess
import cPickle as pickle

PLOT_STUFF = 1

if PLOT_STUFF:
	import matplotlib.pyplot as mpl

MAX_READLEN = 1000
MAX_READS   = 10000000

INIT_SMOOTH = 1000.

def main():

	if len(sys.argv) != 5:
		print '\nUsage:\n\npython qScoreModel.py in.fq maxQscore asciiOffset out.p\n\nExample:\n\npython qScoreModel.py input.fq 41 33 outModel.p\n'
		exit(1)

	infq = sys.argv[1]
	maxQ = int(sys.argv[2])
	RQ   = maxQ+1
	offQ = int(sys.argv[3])
	outM = sys.argv[4]

	print 'opening file...'
	f = open(infq,'r')

	#nReads = int(subprocess.check_output("wc -l "+infq, shell=True).split(' ')[1])/4
	#print 'generating q-score transition matrix... (using '+str(min([MAX_READS,nReads]))+'/'+str(nReads)+' reads)'

	totalQ = np.ones([RQ,RQ])
	rRead  = 0
	actual_readlen = 0
	while True:
		data1 = f.readline()
		data2 = f.readline()
		data3 = f.readline()
		data4 = f.readline()
		if not all([data1,data2,data3,data4]):
			break

		qscore = [ord(n)-offQ for n in data4[:-1]]
		if actual_readlen == 0:
			actual_readlen = len(qscore)
			priorQ = np.ones([actual_readlen,RQ])

		for i in range(len(qscore)):
			if i == 0:
				priorQ[i][qscore[i]] += 1
			else:
				totalQ[qscore[i-1],qscore[i]] += 1
				priorQ[i][qscore[i]] += 1
		rRead += 1
		if rRead%10000 == 0:
			print rRead
			#break
		if rRead >= MAX_READS:
			break
	#print totalQ
	f.close()

	probQ  = [[0. for m in xrange(RQ)] for n in xrange(RQ)]
	for i in xrange(RQ):
		rowSum = float(np.sum(totalQ[i,:]))
		for j in xrange(RQ):
			probQ[i][j] = totalQ[i][j]/rowSum

	initQ  = [[INIT_SMOOTH for m in xrange(RQ)] for n in xrange(actual_readlen)]
	for i in xrange(actual_readlen):
		rowSum = float(np.sum(priorQ[i,:]))+INIT_SMOOTH*RQ
		for j in xrange(RQ):
			initQ[i][j] = (priorQ[i][j]+INIT_SMOOTH)/rowSum

	if PLOT_STUFF:
		mpl.figure(1)
		Z = np.array(probQ)
		X, Y = np.meshgrid( range(0,len(Z[0])+1), range(0,len(Z)+1) )
		mpl.pcolormesh(X,Y,Z[::-1],vmin=0.,vmax=0.25)
		mpl.axis([0,len(Z[0]),0,len(Z)])
		mpl.yticks(range(0,len(Z)),range(len(Z)-1,-1,-1))
		mpl.xticks(range(0,len(Z[0])),range(0,len(Z[0])))
		mpl.xlabel('next Q')
		mpl.ylabel('current Q')
		mpl.title('Q-Score Transition Probabilities')
		mpl.colorbar()

		mpl.figure(2)
		Z = np.array(initQ).T
		X, Y = np.meshgrid( range(0,len(Z[0])+1), range(0,len(Z)+1) )
		mpl.pcolormesh(X,Y,Z,vmin=0.,vmax=0.25)
		mpl.axis([0,len(Z[0]),0,len(Z)])
		mpl.yticks(range(0,len(Z),10),range(0,len(Z),10))
		mpl.xticks(range(0,len(Z[0]),10),range(0,len(Z[0]),10))
		mpl.xlabel('Read Position')
		mpl.ylabel('Quality Score')
		mpl.title('Q-Score Prior Probabilities')
		mpl.colorbar()

		mpl.show()

	Qscores = range(RQ)
	allPMFs = initQ
	pickle.dump([allPMFs,probQ,Qscores,offQ],open(outM,'wb'))

if __name__ == '__main__':
	main()