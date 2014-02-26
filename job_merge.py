#!/usr/bin/env python
# encoding: utf-8

""" **************************************************

job_merge.py outData/myData

- merge all fq, sam, vcf files with a specific prefix


Written by:		Zach Stephens
Date:			February 18, 2014
Contact:		zstephe2@illinois.edu

************************************************** """

import sys
import os
import time
import re
import glob

BUFFER_SIZE = 512000000

def writeInChunks(of,f):
	outData = f.read(BUFFER_SIZE)
	while outData != '':
		of.write(outData)
		outData = f.read(BUFFER_SIZE)

def main():
	dn = sys.argv[1]

	fq1List = glob.glob(dn+'_job*read1.fq')
	fq2List = glob.glob(dn+'_job*read2.fq')
	samList = glob.glob(dn+'_job*.sam')
	vcfList = glob.glob(dn+'_job*.vcf')

	tt = time.time()
	print 'merging files...'

	of = open(dn+'_read1.fq','wb')
	for fn in fq1List:
		f = open(fn,'r')
		writeInChunks(of,f)
		f.close()
	of.close()

	of = open(dn+'_read2.fq','wb')
	for fn in fq2List:
		f = open(fn,'r')
		writeInChunks(of,f)
		f.close()
	of.close()

	if len(samList) > 0:
		of = open(dn+'_golden.sam','wb')
	for fn in samList:
		f = open(fn,'r')
		if fn == samList[0]:
			writeInChunks(of,f)
		else:
			fch = ''
			fchPrev = ''
			while True:
				(ft,fch) = (f.tell(),f.read(1))
				if fchPrev == '\n' and fch != '@':
					break
				fchPrev = fch
			f.seek(ft)
			writeInChunks(of,f)
		f.close()
	if len(samList) > 0:
		of.close()

	variantCov = {}
	for fn in vcfList:
		f = open(fn,'r')
		fch = ''
		fchPrev = ''
		while True:
			(ft,fch) = (f.tell(),f.read(1))
			if fchPrev == '\n' and fch != '#':
				break
			fchPrev = fch
		f.seek(0)
		header = f.read(ft)
		vcfDat = f.read()
		f.close()

		for n in vcfDat.split('\n'):
			splt = n.split('\t')
			if len(splt) == 8:
				var  = (splt[0],splt[1],splt[2],splt[3],splt[4],splt[5],splt[6])
				if var in variantCov:
					variantCov[var] += int(splt[7][3:])
				else:
					variantCov[var] = int(splt[7][3:])

	if len(vcfList) > 0:
		allVariants = [(n[0],n[1],n[2],n[3],n[4],n[5],n[6],'DP='+str(variantCov[n])) for n in variantCov.keys()]

		allVariants = [x for (y,x) in sorted(zip([int(n[1]) for n in allVariants],allVariants))]

		of = open(dn+'_golden.vcf','wb')
		of.write(header)
		for n in allVariants:
			for m in n:
				if m == n[-1]:
					of.write(m+'\n')
				else:
					of.write(m+'\t')
		of.close()

	if len(samList) > 0:
		print '\nWarning: If your simulated SAM files contained reads for'
		print '         multiple reference sequences then the final merged'
		print '         file will require sorting. Otherwise some downstream'
		print '         applications might be unhappy.\n'

	print time.time()-tt,'(sec)\n\n'

if __name__ == '__main__':
	main()