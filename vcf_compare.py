#!/usr/bin/env python
# encoding: utf-8

""" **************************************************

vcf_compare.py

************************************************** """

import sys
import os
import copy
import time
import re
import numpy as np
import matplotlib.pyplot as mpl

from misc	import histoxy

GOLDEN_VCF   = 'outData/testData1_results/testDataset1_golden.vcf'
WORKFLOW_VCF = 'outData/testData1_results/ResultOf_testDataset1_chr1.realrecal.output.bam.chr1.raw.vcf'
WORKFLOW_BAM = 'outData/testData1_results/testDataset1.wdups.sorted.bam'
REFERENCE    = 'biocluster_refs/biocluster_UCSC_hg19_chr1.fa'

MPILEUP_EXEC = '/Users/zach/Desktop/test_data/samtools mpileup'

BPRANGE = 20

def main():

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
	print '\n**********************************\n'

	notFound = [n for n in sorted(correctHashed.keys()) if correctHashed[n] == 1]

	#print len(notFound),'+',len(FPvariants),'=',len(notFound)+len(FPvariants)
	#print nPerfect

	#XLIM = [800000,1000000]
	#for n in notFound:
	#	if n[0] > XLIM[1]:
	#		break
	#	elif n[0] >= XLIM[0] and n[0] < XLIM[1]:
	#		print n
	#print '\n---------------\n'
	#for (n,ei) in FPvariants:
	#	if n[0] > XLIM[1]:
	#		break
	#	elif n[0] >= XLIM[0] and n[0] < XLIM[1]:
	#		print n

	# test to see if any of the FP variants are equivalent to a FN
	nEquiv = 0
	delList_i = []
	delList_j = []
	for i in xrange(len(FPvariants)):
		#print FPvariants[i]
		pos = FPvariants[i][0][0]
		lr  = len(FPvariants[i][0][1])
		refSection = myDat[pos-BPRANGE-1:pos+BPRANGE]
		altSection = refSection[:BPRANGE] + FPvariants[i][0][2] + refSection[BPRANGE+lr:]
		altSection = str(altSection)
		#print altSection

		for j in xrange(len(notFound)):
			pos2 = notFound[j][0]
			lr2  = len(notFound[j][1])
			if abs(pos-pos2) <= 2*BPRANGE:
				dpos = pos-pos2
				#print 'rawr',dpos
				#print notFound[j]
				#print refSection[:BPRANGE+dpos], BPRANGE+dpos
				#print refSection[BPRANGE+dpos:BPRANGE+dpos+lr2], BPRANGE+dpos+lr2
				#print refSection[BPRANGE+dpos+lr2:]
				altSection2 = refSection[:BPRANGE-dpos] + notFound[j][2] + refSection[BPRANGE-dpos+lr2:]
				altSection2 = str(altSection2)
				#print refSection
				#print altSection
				#print altSection2
				if altSection2 == altSection:
					#print 'we have a winner, jerry!'
					#print FPvariants[i], notFound[j]
					#print refSection
					#print altSection
					delList_i.append(i)
					delList_j.append(j)
					nEquiv += 1
					#break
	nPerfect += nEquiv

	for i in sorted(delList_i,reverse=True):
		del FPvariants[i]
	for j in sorted(delList_j,reverse=True):
		del notFound[j]

	print '\n**********************************\n'
	print 'Total Golden Variants:',totalVariants,'\n'
	print 'Perfect Matches:',nPerfect,'({:.2f}%)'.format(100.*float(nPerfect)/totalVariants)
	print 'FP variants:   ',len(FPvariants),'({:.2f}%)'.format(100.*float(len(FPvariants))/totalVariants)
	print 'FN variants:   ',len(notFound),'({:.2f}%)'.format(100.*float(len(notFound))/totalVariants)
	print '\nNumber of equivalent variants denoted differently:',nEquiv,'\n'

	for n in FPvariants:
		print n
	print '\n**********************************\n'
	of = open('FN_pos.txt','w')
	for n in notFound:
		print n,
		if n in correctCov:
			print ', DP =',correctCov[n]
		else:
			print ''
		of.write(refName+'\t'+str(n[0])+'\n')
	of.close()


	""" ///////////////////////////////////////////////// """


	cmd = MPILEUP_EXEC+' -l FN_pos.txt '+WORKFLOW_BAM+' 1>pileOut.txt'
	#os.system(cmd)

	f = open('pileOut.txt','r')
	for line in f:
		print line
	f.close()

	#mpl.figure(0)
	#mpl.subplot(211)
	#ind = []
	#err = []
	#for n in notFound:
	#	ind.extend([n[0]-1,n[0],n[0]+1])
	#	err.extend([0,1,0])
	#mpl.plot(ind,err)
	#mpl.xlim(XLIM)
	#mpl.title('False negatives')

	#mpl.subplot(212)
	#ind = []
	#err = []
	#for n in FPvariants:
	#	ind.extend([n[0][0]-1,n[0][0],n[0][0]+1])
	#	err.extend([0,1,0])
	#mpl.plot(ind,err)
	#mpl.xlim(XLIM)
	#mpl.title('False positives')

	#mpl.show()

if __name__ == '__main__':
	main()