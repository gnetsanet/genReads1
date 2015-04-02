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

	# sort by job order so we put the files back together in a sensible order
	fq1List = sorted([(int(re.findall(r"_job\d+",n)[0][4:]), n) for n in fq1List])
	fq1List = [n[1] for n in fq1List]
	fq2List = sorted([(int(re.findall(r"_job\d+",n)[0][4:]), n) for n in fq2List])
	fq2List = [n[1] for n in fq2List]
	samList = sorted([(int(re.findall(r"_job\d+",n)[0][4:]), n) for n in samList])
	samList = [n[1] for n in samList]
	vcfList = sorted([(int(re.findall(r"_job\d+",n)[0][4:]), n) for n in vcfList])
	vcfList = [n[1] for n in vcfList]

	tt = time.time()

	#of = open(dn+'_read1.fq','wb')
	#for fn in fq1List:
	#	f = open(fn,'r')
	#	writeInChunks(of,f)
	#	f.close()
	#of.close()
	print 'merging FASTQ files...'
	os.system('cat '+' '.join(fq1List)+' > '+dn+'_read1.fq')

	#of = open(dn+'_read2.fq','wb')
	#for fn in fq2List:
	#	f = open(fn,'r')
	#	writeInChunks(of,f)
	#	f.close()
	#of.close()
	os.system('cat '+' '.join(fq2List)+' > '+dn+'_read2.fq')

	#if len(samList) > 0:
	#	of = open(dn+'_golden.sam','wb')
	#for fn in samList:
	#	f = open(fn,'r')
	#	if fn == samList[0]:
	#		writeInChunks(of,f)
	#	else:
	#		fch = ''
	#		fchPrev = ''
	#		while True:
	#			(ft,fch) = (f.tell(),f.read(1))
	#			if fchPrev == '\n' and fch != '@':
	#				break
	#			fchPrev = fch
	#		f.seek(ft)
	#		writeInChunks(of,f)
	#	f.close()
	#if len(samList) > 0:
	#	of.close()
	if len(samList) > 0:
		print 'merging SAM files...'
		os.system('grep -E "(^@)" '+samList[0]+' > '+dn+'_golden.sam')
		os.system('grep -vhE "(^@|^$)" '+' '.join(samList)+' >> '+dn+'_golden.sam')

	if len(vcfList) > 0:
		print 'merging VCF files...'
	variantInf = {}
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
				if var in variantInf:
					variantInf[var] += splt[7]+';'
				else:
					variantInf[var] = splt[7]+';'

	for n in variantInf.keys():
		targeted = ''
		coverage = ''
		readsCov = ''
		alleleFr = ''
		soi = variantInf[n]

		if 'DP=' in soi:
			coverage = 'DP='+str(sum([int(m[3:]) for m in re.findall(r"(DP=.*?)(?=;)",soi)]))

		if 'TARGETED=1' in soi:
			targeted = 'TARGETED=1;'

		if 'AF=' in soi:
			alleleFr = 'AF='+re.findall(r"(AF=.*?)(?=;)",soi)[0]

		rstrs = [m[6:] for m in re.findall(r"(READS=.*?)(?=;)",soi)]
		if rstrs == ['']:
			readsCov == ''
		else:
			readsCov = 'READS=' + ','.join(rstrs)

		variantInf[n] = coverage+';'+targeted+alleleFr+readsCov

	if len(vcfList) > 0:
		allVariants = [(n[0],int(n[1]),n[2],n[3],n[4],n[5],n[6],variantInf[n]) for n in variantInf.keys()]

		#allVariants = [x for (y,x) in sorted(zip([int(n[1]) for n in allVariants],allVariants))]
		allVariants = sorted(allVariants)

		of = open(dn+'_golden.vcf','wb')
		of.write(header)
		for n in allVariants:
			for m in n:
				if m == n[-1]:
					of.write(m+'\n')
				else:
					of.write(str(m)+'\t')
		of.close()

	if len(samList) > 0:
		print '\nWarning: If your simulated SAM files contained reads for'
		print '         multiple reference sequences then the final merged'
		print '         file will require sorting. Otherwise some downstream'
		print '         applications might be unhappy.\n'

	print time.time()-tt,'(sec)\n\n'

if __name__ == '__main__':
	main()