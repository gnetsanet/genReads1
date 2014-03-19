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
import re
import numpy as np
#import matplotlib.pyplot as mpl

from misc	import histoxy

#GOLDEN_VCF   = 'outData/testData1_results/testDataset1_golden.vcf'
#WORKFLOW_VCF = 'outData/testData1_results/ResultOf_testDataset1_chr1.realrecal.output.bam.chr1.raw.vcf'
#WORKFLOW_BAM = 'outData/testData1_results/testDataset1.wdups.sorted.bam'
#REFERENCE    = 'various_data/biocluster_UCSC_hg19_chr1.fa'

MPILEUP_EXEC = '/Users/zach/Desktop/bioinformatics_stuff/samtools mpileup'

BPRANGE = 20

def main():

	if len(sys.argv) != 5:
		print '\nUsage:\n\npython vcf_compare.py ref.fa golden.vcf workflow.vcf workflow.bam\n'
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
		print '{:,} bp\n'.format(myLen)
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

	correctVariants = []
	correctHashed   = {}
	correctCov      = {}
	for line in open(GOLDEN_VCF,'r'):
		if line[0] != '#':
			splt = line.split('\t')
			var  = (int(splt[1]),splt[3],splt[4])
			correctVariants.append(var)
			correctHashed[var] = 1
			if 'DP=' in splt[7]:
				correctCov[var] = int(splt[7][3:])
	#print correctVariants
	totalVariants = len(correctVariants)

	workflowVariants = []
	for line in open(WORKFLOW_VCF,'r'):
		if line[0] != '#':
			splt = line.split('\t')
			qual = float(splt[5])
			cov  = int(re.findall(r"DP=[0-999999]",splt[7])[0][3:])
			workflowVariants.append([(int(splt[1]),splt[3],splt[4]),[qual,cov]])

	nPerfect = 0
	FPvariants = []
	for [var,extraInfo] in workflowVariants:
		if var in correctHashed:
			nPerfect += 1
			correctHashed[var] = 2
		else:
			FPvariants.append([var,extraInfo])


	notFound = [n for n in sorted(correctHashed.keys()) if correctHashed[n] == 1]


	delList_i = []
	delList_j = []
	regionsToCheck = []
	for i in xrange(len(FPvariants)):
		pos = FPvariants[i][0][0]
		regionsToCheck.append((max([pos-BPRANGE-1,0]),min([pos+BPRANGE,len(myDat)-1])))

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

		#print fpWithin, nfWithin
		#print altSection
		#print altSection2
		if altSection == altSection2:
			for (m,i) in fpWithin:
				if i not in delList_i:
					delList_i.append(i)
			for (m,j) in nfWithin:
				if j not in delList_j:
					delList_j.append(j)

	nEquiv = 0
	for i in sorted(delList_i,reverse=True):
		del FPvariants[i]
	for j in sorted(delList_j,reverse=True):
		del notFound[j]
		nEquiv += 1
	nPerfect += nEquiv

	print '\n**********************************\n'
	print 'Total Golden Variants:',totalVariants,'\n'
	print 'Perfect Matches:',nPerfect,'({:.2f}%)'.format(100.*float(nPerfect)/totalVariants)
	print 'FP variants:   ',len(FPvariants),'({:.2f}%)'.format(100.*float(len(FPvariants))/totalVariants)
	print 'FN variants:   ',len(notFound),'({:.2f}%)'.format(100.*float(len(notFound))/totalVariants)
	print '\nNumber of equivalent variants denoted differently:',nEquiv
	print '\n**********************************\n'
	
	for n in FPvariants:
		print n
	print '\n**********************************\n'
	of = open('FN_pos.txt','w')

	nfsnps = 0
	for n in notFound:
		if len(n[1]) == len(n[2]):
			nfsnps += len(n[1])
		print n,
		if n in correctCov:
			print ', DP =',correctCov[n]
		else:
			print ''
		of.write(refName+'\t'+str(n[0])+'\n')
	print 'rawr:',nfsnps
	of.close()


	""" ///////////////////////////////////////////////// """


	cmd = MPILEUP_EXEC+' -l FN_pos.txt '+WORKFLOW_BAM+' 1>pileOut.txt'
	#os.system(cmd)

	#f = open('pileOut.txt','r')
	#for line in f:
	#	print line
	#f.close()


if __name__ == '__main__':
	main()