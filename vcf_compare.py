#!/usr/bin/env python
# encoding: utf-8

""" **************************************************

vcf_compare.py

- compare vcf file produced by workflow to golden vcf produced by simulator

Written by:		Zach Stephens
Date:			February 13, 2014
Contact:		zstephe2@illinois.edu


Usage:

python vcf_compare.py ref.fa golden.vcf workflow.vcf workflow.bam

************************************************** """

import sys
import os
import copy
import time
import bisect
import re
import numpy as np
#import matplotlib.pyplot as mpl

from misc	import histoxy, needleman_wunsch

SAMTOOLS_EXEC = '/Users/zach/Desktop/bioinformatics_stuff/samtools'

EV_BPRANGE = 20		# how far to either side of a particular variant location do we want to check for equivalents?
RR_BPRANGE = 70		# how far to either side of a particular variant location do we want to check for region repeats?

RR_THRESH  = RR_BPRANGE/4

def main():

	if len(sys.argv) != 5 and len(sys.argv) != 7:
		print '\nUsage:\n\npython vcf_compare.py ref.fa golden.vcf workflow.vcf workflow.sam\n'
		print 'OR:'
		print '\nUsage:\n\npython vcf_compare.py ref.fa golden.vcf workflow.vcf workflow.sam targetedRegions.bed minRegionLen\n'
		exit(1)

	GOLDEN_VCF   = sys.argv[2]
	WORKFLOW_VCF = sys.argv[3]
	WORKFLOW_BAM = sys.argv[4]
	REFERENCE    = sys.argv[1]

	ref = []
	f = open(REFERENCE,'r')
	nLines = 0
	prevR = None
	prevP = None
	ref_inds = []
	print 'indexing reference fasta...'
	while 1:
		nLines += 1
		data = f.readline()
		if not data:
			ref_inds.append( (prevR, prevP, f.tell()-len(data)) )
			break
		if data[0] == '>':
			if prevP != None:
				ref_inds.append( (prevR, prevP, f.tell()-len(data)) )
			prevP = f.tell()
			prevR = data[1:-1]

	for n_RI in ref_inds:
		refName = n_RI[0]
		# assumes fasta file is sane and has '\n' characters within long sequences
		f.seek(n_RI[1])
		print 'reading '+refName+'...',
		myDat  = f.read(n_RI[2]-n_RI[1]).split('\n')
		myLen  = sum([len(m) for m in myDat])
		if sys.version_info >= (2,7):
			print '{:,} bp\n'.format(myLen)
		else:
			print '{0:} bp\n'.format(myLen)
		inWidth = len(myDat[0])
		if len(myDat[-1]) == 0:	# if last line is empty, remove it.
			del myDat[-1]
		if inWidth*(len(myDat)-1)+len(myDat[-1]) != myLen:
			print 'fasta column-width not consistent.'
			print myLen, inWidth*(len(myDat)-1)+len(myDat[-1])
			for i in xrange(len(myDat)):
				if len(myDat[i]) != inWidth:
					print i, len(myDat[i]), inWidth
			exit(1)

	myDat = bytearray(''.join(myDat)).upper()
	myLen = len(myDat)


	targRegionsFl = []
	if len(sys.argv) == 7:
		BEDFILE  = sys.argv[5]
		bedfile = open(BEDFILE,'r')
		minRegionLen = int(sys.argv[6])
		for line in bedfile:
			splt = line.split('\t')
			if splt[0] == refName:
				targRegionsFl.extend((int(splt[1]),int(splt[2])))
	else:
		BEDFILE  = None
		targRegionsFl = [-1,myLen+1]


	correctVariants = []
	correctHashed   = {}
	correctCov      = {}
	correctReads    = {}
	correctTargLen  = {}
	nBelowMinRLen   = 0
	for line in open(GOLDEN_VCF,'r'):
		if line[0] != '#':
			splt = line.split('\t')
			var  = (int(splt[1]),splt[3],splt[4])
			targInd = bisect.bisect(targRegionsFl,var[0])

			if targInd%2 == 1:
				targLen = targRegionsFl[targInd]-targRegionsFl[targInd-1]
				if (BEDFILE != None and targLen >= minRegionLen) or BEDFILE == None:
					correctVariants.append(var)
					correctHashed[var] = 1
					if 'DP=' in splt[7]:
						correctCov[var] = int(re.findall(r"DP=[0-999999]",splt[7])[0][3:])
					if 'READS=' in splt[7]:
						rstrings = re.findall(r"(READS=.*?)(?=;)",splt[7]) + re.findall(r"(READS=.*?)(?=\n)",splt[7])
						correctReads[var] = ''.join([m[6:] for m in rstrings]).split(',')
					correctTargLen[var] = targLen
				else:
					nBelowMinRLen += 1

	#print correctVariants

	workflowVariants = []
	for line in open(WORKFLOW_VCF,'r'):
		if line[0] != '#':
			splt = line.split('\t')
			var  = (int(splt[1]),splt[3],splt[4])
			targInd = bisect.bisect(targRegionsFl,var[0])

			if targInd%2 == 1:
				targLen = targRegionsFl[targInd]-targRegionsFl[targInd-1]
				if (BEDFILE != None and targLen >= minRegionLen) or BEDFILE == None:
					qual = float(splt[5])
					cov  = int(re.findall(r"DP=[0-999999]",splt[7])[0][3:])
					workflowVariants.append([var,[qual,cov,targLen]])

	nPerfect = 0
	FPvariants = []
	for [var,extraInfo] in workflowVariants:
		if var in correctHashed:
			nPerfect += 1
			correctHashed[var] = 2
		else:
			FPvariants.append([var,extraInfo])

	
	notFound = [n for n in sorted(correctHashed.keys()) if correctHashed[n] == 1]


	totalVariants = nPerfect + len(notFound)


	# let's check for equivalent variants
	delList_i = []
	delList_j = []
	regionsToCheck = []
	for i in xrange(len(FPvariants)):
		pos = FPvariants[i][0][0]
		regionsToCheck.append((max([pos-EV_BPRANGE-1,0]),min([pos+EV_BPRANGE,len(myDat)-1])))

	for n in regionsToCheck:
		refSection = myDat[n[0]:n[1]]

		fpWithin = []
		for i in xrange(len(FPvariants)):
			m  = FPvariants[i][0]
			if (m[0] > n[0] and m[0] < n[1]):
				fpWithin.append((m,i))
		fpWithin = sorted(fpWithin)
		adj = 0
		altSection = copy.deepcopy(refSection)
		for (m,i) in fpWithin:
			lr = len(m[1])
			la = len(m[2])
			dpos = m[0]-n[0]+adj
			altSection = altSection[:dpos-1] + m[2] + altSection[dpos-1+lr:]
			adj += la-lr

		nfWithin = []
		for j in xrange(len(notFound)):
			m = notFound[j]
			if (m[0] > n[0] and m[0] < n[1]):
				nfWithin.append((m,j))
		nfWithin = sorted(nfWithin)
		adj = 0
		altSection2 = copy.deepcopy(refSection)
		for (m,j) in nfWithin:
			lr = len(m[1])
			la = len(m[2])
			dpos = m[0]-n[0]+adj
			altSection2 = altSection2[:dpos-1] + m[2] + altSection2[dpos-1+lr:]
			adj += la-lr

		if altSection == altSection2:
			for (m,i) in fpWithin:
				if i not in delList_i:
					delList_i.append(i)
			for (m,j) in nfWithin:
				if j not in delList_j:
					delList_j.append(j)

	nEquiv = 0
	for i in sorted(list(set(delList_i)),reverse=True):
		del FPvariants[i]
	for j in sorted(list(set(delList_j)),reverse=True):
		del notFound[j]
		nEquiv += 1
	nPerfect += nEquiv


	print '\n**********************************\n'
	if BEDFILE != None:
		print 'ONLY CONSIDERING VARIANTS FOUND WITHIN TARGETED REGIONS\n\n'
	print 'Total Golden Variants:',totalVariants,'\n'
	print 'Perfect Matches:',nPerfect,'({0:.2f}%)'.format(100.*float(nPerfect)/totalVariants)
	print 'FP variants:   ',len(FPvariants),'({0:.2f}%)'.format(100.*float(len(FPvariants))/totalVariants)
	print 'FN variants:   ',len(notFound),'({0:.2f}%)'.format(100.*float(len(notFound))/totalVariants)
	print '\nNumber of equivalent variants denoted differently between the two vcfs:',nEquiv
	if BEDFILE != None:
		print '\nNumber of golden variants located in targeted regions that were too small to be sampled from:',nBelowMinRLen
	print '\n**********************************\n'


	# check if repeats are responsible for any pair of FN / FP variants
	potential_pairs = []
	delList_i = []
	delList_j = []
	#for i in xrange(len(FPvariants)):
	#	pos = FPvariants[i][0][0]
	#	roi = (max([pos-RR_BPRANGE-1,0]),min([pos+RR_BPRANGE,len(myDat)-1]))
	#	refSection1 = myDat[roi[0]:roi[1]]
	#
	#	for j in xrange(len(notFound)):
	#		if (FPvariants[i][0][1],FPvariants[i][0][2]) == (notFound[j][1],notFound[j][2]):
	#			pos = notFound[j][0]
	#			roi = (max([pos-RR_BPRANGE-1,0]),min([pos+RR_BPRANGE,len(myDat)-1]))
	#			refSection2 = myDat[roi[0]:roi[1]]
	#
	#			(score,a1,a2) = needleman_wunsch(1, 1, str(refSection1), str(refSection2))
	#			if score <= RR_THRESH:
	#				potential_pairs.append((i,j,score))
	#				delList_i.append(i)
	#				delList_j.append(j)

	print 'Pairs of variants that originated from within one instance of a repetitive region, but were called by the workflow in another because the reads from the original region were confidently mapped elsewhere:',len(potential_pairs),'\n'
	for n in potential_pairs:
		print 'golden:\t\t',notFound[n[1]]
		print 'workflow:\t',FPvariants[n[0]]
		print 'NW_score:\t',n[2]
	print '\n**********************************\n'

	for i in sorted(list(set(delList_i)),reverse=True):
		del FPvariants[i]
	for j in sorted(list(set(delList_j)),reverse=True):
		del notFound[j]

	


	""" ///////////////////////////////////////////////// """

	for n in FPvariants:
		print n
	print '\n**********************************\n'

	for n in notFound:
		print n,
		#if n in correctReads:
		#	for roi in correctReads[n]:
		#		print roi
		#	die
		if n in correctCov:
			print 'DP =',correctCov[n],
			if BEDFILE != None:
				print 'TL =',correctTargLen[n]
			else:
				print ''
		else:
			print ''


if __name__ == '__main__':
	main()