#!/usr/bin/env python
# encoding: utf-8

""" **************************************************

vcf_compare.py

- compare vcf file produced by workflow to golden vcf produced by simulator

Written by:		Zach Stephens
Date:			January 20, 2015
Contact:		zstephe2@illinois.edu

************************************************** """

import sys
import os
import copy
import time
import bisect
import re
import numpy as np
import optparse

#
#	mappability track functions
#
EIGHTBIT     = {}
REV_EIGHTBIT = {}
for i in xrange(256):
	key = bytearray('\0\0\0\0\0\0\0\0')
	if i&128: key[0] = '\1'
	if i&64: key[1] = '\1'
	if i&32: key[2] = '\1'
	if i&16: key[3] = '\1'
	if i&8: key[4] = '\1'
	if i&4: key[5] = '\1'
	if i&2: key[6] = '\1'
	if i&1: key[7] = '\1'
	EIGHTBIT[str(key)]   = chr(i)
	REV_EIGHTBIT[chr(i)] = key

def compressTrack(t):
	len_t = len(t)
	padded = len_t-(len_t%8)+8
	while len(t) != padded:
		t.append('\0')
	outStr = ''
	for i in xrange(0,padded,8):
		outStr += EIGHTBIT[str(t[i:i+8])]
	return (len_t,outStr.encode("zlib"))

def decompressTrack(ct):
	outBA = bytearray()
	dc = ct[1].decode("zlib")
	for i in xrange(len(dc)):
		outBA.extend(REV_EIGHTBIT[dc[i]])
	return outBA[:ct[0]]
#
#
#

EV_BPRANGE = 50		# how far to either side of a particular variant location do we want to check for equivalents?

DEFAULT_QUAL = -666	# if we can't find a qual score, use this instead so we know it's missing

MAX_VAL = 9999999999999	# an unreasonably large value that no reference fasta could concievably be longer than

DESC   = """%prog: vcf comparison script."""
VERS   = 0.1

PARSER = optparse.OptionParser('python %prog [options] -r <ref.fa> -g <golden.vcf> -w <workflow.vcf>',description=DESC,version="%prog v"+str(VERS))

PARSER.add_option('-r', help='* Reference Fasta',   dest='REFF', action='store', metavar='<ref.fa>')
PARSER.add_option('-g', help='* Golden VCF',        dest='GVCF', action='store', metavar='<golden.vcf>')
PARSER.add_option('-w', help='* Workflow VCF',      dest='WVCF', action='store', metavar='<workflow.vcf>')
PARSER.add_option('-m', help='Mappability Track', dest='MTRK', action='store', metavar='<tracks.dat>')
PARSER.add_option('-t', help='Targetted Regions', dest='TREG', action='store', metavar='<regions.bed>')
PARSER.add_option('-T', help='Min Region Len',    dest='MTRL', action='store', metavar='<int>')

PARSER.add_option('--no-plot', help="No plotting [%default]", dest='NO_PLOT', default=False, action='store_true')
PARSER.add_option('--fast', help="No equivalent variant detection [%default]", dest='FAST', default=False, action='store_true')

(OPTS,ARGS) = PARSER.parse_args()

REFERENCE    = OPTS.REFF
GOLDEN_VCF   = OPTS.GVCF
WORKFLOW_VCF = OPTS.WVCF
MAPTRACK     = OPTS.MTRK
BEDFILE      = OPTS.TREG
NO_PLOT      = OPTS.NO_PLOT
FAST         = OPTS.FAST
if OPTS.MTRL != None:
	MINREGIONLEN = int(OPTS.MTRL)
else:
	MINREGIONLEN = None

if REFERENCE == None:
	print 'Error: No reference provided.'
	exit(1)
if GOLDEN_VCF == None:
	print 'Error: No golden VCF provided.'
	exit(1)
if WORKFLOW_VCF == None:
	print 'Error: No workflow VCF provided.'
	exit(1)
if (BEDFILE != None and MINREGIONLEN == None) or (BEDFILE == None and MINREGIONLEN != None):
	print 'Error: Both -t and -T must be specified'
	exit(1)

if NO_PLOT == False:
	import matplotlib.pyplot as mpl
	from matplotlib_venn import venn2, venn3

AF_STEPS = 20
AF_KEYS  = np.linspace(0.0,1.0,AF_STEPS+1)

def quantize_AF(af):
	if af >= 1.0:
		return AF_STEPS
	elif af <= 0.0:
		return 0
	else:
		return int(af*AF_STEPS)

def parseLine(splt):

	cov  = None
	af   = 1.0
	qual = DEFAULT_QUAL
	alt_alleles = []

	#	any alt alleles?
	alt_split = splt[4].split(',')
	if len(alt_split) > 1:
		alt_alleles = alt_split

	#	do we have the necessary fields to compute DP/AF stats?
	if len(splt) >= 8:

		#	get variant coverage
		#
		if 'DP=' in splt[7]:
			cov = int(re.findall(r"DP=[0-9]+",splt[7])[0][3:])
		# check for different formatting
		elif len(splt) >= 9:
			splt[8] += ':'
			if ':DP:' in splt[8]:
				dpInd = splt[8].split(':').index('DP')
				cov   = int(splt[9].split(':')[dpInd])

		#	get variant AF (for heterozygous genotypes)
		#
		af = 1.0
		if 'AF=' in splt[7]:
			af  = float(re.findall(r"AF=.*?(?=;)",splt[7])[0][3:])
		# check for different formatting
		elif len(splt) >= 9:
			splt[8] += ':'
			if ':AF:' in splt[8]:
				afInd = splt[8].split(':').index('AF')
				af    = float(splt[9].split(':')[afInd])

		#	get variant call quality
		#
		if splt[5] != '.':
			qual = float(splt[5])

	return (cov, af, qual, alt_alleles)

def condenseAlts(listIn,altsList,FNorFP):
	to_condense   = {}
	ext_info_dict = {}
	for i in xrange(len(listIn)):
		if FNorFP: var = listIn[i]
		else: [var,extra_info] = listIn[i]
		if var in altsList:
			concat = (var[0],var[1],','.join([n[2] for n in altsList[var]]))
			if not FNorFP: ext_info_dict[concat] = extra_info
			if concat not in to_condense:
				to_condense[concat] = []
			to_condense[concat].append(i)
	delList = [j for i in to_condense.values() for j in i]
	outList = [listIn[i] for i in xrange(len(listIn)) if i not in delList]
	for n in to_condense.keys():
		if FNorFP: outList.append(n)
		else: outList.append([n,ext_info_dict[n]])
	return outList


def main():

	ref = []
	f = open(REFERENCE,'r')
	nLines = 0
	prevR = None
	prevP = None
	ref_inds = []
	sys.stdout.write('\nindexing reference fasta... ')
	sys.stdout.flush()
	tt = time.time()
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
	print '{0:.3f} (sec)'.format(time.time()-tt)

	ztV = 0
	znP = 0
	zfP = 0
	znF = 0
	znE = 0
	if BEDFILE != None:
		zbM = 0

	mappability_vs_FN = {0:0, 1:0}	# [0] = # of FNs that were in mappable regions, [1] = # of FNs that were in unmappable regions
	coverage_vs_FN    = {}			# [C] = # of FNs that were covered by C reads
	alleleBal_vs_FN   = {}			# [AF] = # of FNs that were heterozygous genotypes with allele freq AF (quantized to multiples of 1/AF_STEPS)
	for n in AF_KEYS:
		alleleBal_vs_FN[n] = 0

	#
	#	read in mappability track
	#
	mappability_tracks = {}
	if MAPTRACK != None:
		print 'reading mappability tracks...'
		mapf = open(MAPTRACK,"rb")
		[headerLen] = np.fromfile(mapf,'<i4',1)
		header = mapf.read(headerLen)
		names  = []
		lens   = []
		bytes  = []
		for line in header.split('\n')[:-1]:
			splt = line.split('\t')
			names.append(splt[0])
			lens.append(int(splt[1]))
			bytes.append(int(splt[2]))

		for i in xrange(len(names)):
			mappability_tracks[names[i]] = decompressTrack((lens[i],mapf.read(bytes[i])))

	#
	#	data for plotting FN analysis
	#
	set1 = []
	set2 = []
	set3 = []
	varAdj = 0

	#
	#
	#	For each sequence in reference fasta...
	#
	#
	for n_RI in ref_inds:
		refName = n_RI[0]

		#
		#	Read reference sequence we're currently interested in
		#
		#	assumes fasta file is sane and has '\n' characters within long sequences
		if FAST == False:
			f.seek(n_RI[1])
			print 'reading '+refName+'...',
			myDat  = f.read(n_RI[2]-n_RI[1]).split('\n')
			myLen  = sum([len(m) for m in myDat])
			if sys.version_info >= (2,7):
				print '{:,} bp'.format(myLen)
			else:
				print '{0:} bp'.format(myLen)
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

		#
		#	Parse relevant targeted regions
		#
		targRegionsFl = []
		if BEDFILE != None:
			bedfile = open(BEDFILE,'r')
			for line in bedfile:
				splt = line.split('\t')
				if splt[0] == refName:
					targRegionsFl.extend((int(splt[1]),int(splt[2])))
			bedfile.close()
		else:
			targRegionsFl = [-1,MAX_VAL+1]

		#
		#	It begins...
		#
		sys.stdout.write('comparing variation in '+refName+'... ')
		sys.stdout.flush()
		tt = time.time()

		#
		#	Parse relevant golden variants
		#
		correctHashed   = {}
		correct_alts    = {}
		correctCov      = {}
		correctAF       = {}
		correctQual     = {}
		correctTargLen  = {}
		nBelowMinRLen   = 0
		for line in open(GOLDEN_VCF,'r'):
			if line[0] != '#':
				splt = line.split('\t')
				if splt[0] == refName:
					var  = (int(splt[1]),splt[3],splt[4])
					targInd = bisect.bisect(targRegionsFl,var[0])

					if targInd%2 == 1:
						targLen = targRegionsFl[targInd]-targRegionsFl[targInd-1]
						if (BEDFILE != None and targLen >= MINREGIONLEN) or BEDFILE == None:
							
							(cov, af, qual, aa) = parseLine(splt)

							if len(aa):
								allVars = [(var[0],var[1],n) for n in aa]
								for i in xrange(len(allVars)):
									correctHashed[allVars[i]] = 1
									correct_alts[allVars[i]]  = allVars
							else:
								correctHashed[var] = 1

							if cov != None:
								correctCov[var]     = cov
							correctAF[var]      = af
							correctQual[var]    = qual
							correctTargLen[var] = targLen
						else:
							nBelowMinRLen += 1

		#
		#	Parse relevant workflow variants
		#
		workflowVariants = []
		workflow_alts    = {}
		for line in open(WORKFLOW_VCF,'r'):
			if line[0] != '#':
				splt = line.split('\t')
				if splt[0] == refName:
					var  = (int(splt[1]),splt[3],splt[4])
					targInd = bisect.bisect(targRegionsFl,var[0])

					if targInd%2 == 1:
						targLen = targRegionsFl[targInd]-targRegionsFl[targInd-1]
						if (BEDFILE != None and targLen >= MINREGIONLEN) or BEDFILE == None:
							
							(cov, af, qual, aa) = parseLine(splt)

							if len(aa):
								allVars = [(var[0],var[1],n) for n in aa]
								for i in xrange(len(allVars)):
									workflowVariants.append([allVars[i],[cov,af,qual,targLen]])
									workflow_alts[allVars[i]] = allVars
							else:
								workflowVariants.append([var,[cov,af,qual,targLen]])

		#
		#	Deduce which variants are FP / FN
		#
		nPerfect = 0
		FPvariants = []
		alts_to_ignore = []
		for [var,extraInfo] in workflowVariants:
			if var in correctHashed:
				nPerfect += 1
				correctHashed[var] = 2
				if var in correct_alts:
					for v2 in correct_alts[var]:
						correctHashed[v2] = 2
				if var in workflow_alts:
					alts_to_ignore.extend(workflow_alts[var])
			else:
				FPvariants.append([var,extraInfo])

		#	remove any trace of workflow variants who were not found, but whose alternate was
		for i in xrange(len(FPvariants)-1,-1,-1):
			if FPvariants[i][0] in alts_to_ignore:
				del FPvariants[i]
		
		notFound = [n for n in sorted(correctHashed.keys()) if correctHashed[n] == 1]
		#print len(notFound), len(FPvariants)

		#
		#	condense all variants who have alternate alleles and were *not* found to have perfect matches
		#	into a single variant again. These will not be included in the candidates for equivalency checking. Sorry!
		#
		notFound   = condenseAlts(notFound,correct_alts,True)
		FPvariants = condenseAlts(FPvariants,workflow_alts,False)

		#
		#
		totalVariants = nPerfect + len(notFound)
		if totalVariants == 0:
			zfP += len(FPvariants)
			print '{0:.3f} (sec)'.format(time.time()-tt)
			continue

		#
		#	let's check for equivalent variants
		#
		if FAST == False:
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

		#
		#	Tally up errors and whatnot
		#
		ztV += totalVariants
		znP += nPerfect
		zfP += len(FPvariants)
		znF += len(notFound)
		if FAST == False:
			znE += nEquiv
		if BEDFILE != None:
			zbM += nBelowMinRLen

		#
		#	try to identify a reason for FN variants:
		#
		if len(correctCov):
			covKeys = [n for n in correctCov.values() if n != None]
			avg_dp = np.mean(covKeys)
			std_dp = np.std(covKeys)

			DP_THRESH = avg_dp - 2 * std_dp		# below this is unusually low
			AF_THRESH = 0.7						# below this is a het variant with potentially low allele balance

			venn_data = [[0,0,0] for n in notFound]		# [i] = (unmappable, low cov, low het)

			for i in xrange(len(notFound)):
				var = notFound[i]

				#	mappability?
				if MAPTRACK != None:
					if mappability_tracks[refName][var[0]]:
						mappability_vs_FN[1] += 1
						venn_data[i][0] = 1
					else:
						mappability_vs_FN[0] += 1

				#	coverage?
				if var in correctCov:
					c = correctCov[var]
					if c != None:
						if c not in coverage_vs_FN:
							coverage_vs_FN[c] = 0
						coverage_vs_FN[c] += 1
						if c < DP_THRESH:
							venn_data[i][1] = 1

				#	heterozygous genotype messing things up?
				if var in correctAF:
					a = AF_KEYS[quantize_AF(correctAF[var])]
					if a not in alleleBal_vs_FN:
						alleleBal_vs_FN[a] = 0
					alleleBal_vs_FN[a] += 1
					if a < AF_THRESH:
						venn_data[i][2] = 1

			for i in xrange(len(notFound)):
				if venn_data[i][0]: set1.append(i+varAdj)
				if venn_data[i][1]: set2.append(i+varAdj)
				if venn_data[i][2]: set3.append(i+varAdj)
			varAdj += len(notFound)

		print '{0:.3f} (sec)'.format(time.time()-tt)

	#
	#	plot some FN stuff
	#
	if NO_PLOT == False:
		nDetected = len(set(set1+set2+set3))
		set1 = set(set1)
		set2 = set(set2)
		set3 = set(set3)

		mpl.figure(0)
		tstr1 = 'False Negative Variants (Missed Detections)'
		tstr2 = str(nDetected)+' / '+str(znF)+' FN variants categorized'
		if MAPTRACK != None:
			v = venn3([set1, set2, set3], ('Unmappable', 'Low Coverage', 'Heterozygous'))
		else:
			v = venn2([set2, set3], ('Low Coverage', 'Heterozygous'))
		mpl.figtext(0.5,0.95,tstr1,fontdict={'size':14,'weight':'bold'},horizontalalignment='center')
		mpl.figtext(0.5,0.03,tstr2,fontdict={'size':14,'weight':'bold'},horizontalalignment='center')

		mpl.show()


	#
	#	spit out results to console
	#
	print '\n**********************************\n'
	if BEDFILE != None:
		print 'ONLY CONSIDERING VARIANTS FOUND WITHIN TARGETED REGIONS\n\n'
	print 'Total Golden Variants:',ztV,'\n'
	if ztV > 0:
		print 'Perfect Matches:',znP,'({0:.2f}%)'.format(100.*float(znP)/ztV)
		print 'FP variants:   ',zfP,'({0:.2f}%)'.format(100.*float(zfP)/ztV)
		print 'FN variants:   ',znF,'({0:.2f}%)'.format(100.*float(znF)/ztV)
	if FAST == False:
		print '\nNumber of equivalent variants denoted differently between the two vcfs:',znE
	if BEDFILE != None:
		print '\nNumber of golden variants located in targeted regions that were too small to be sampled from:',zbM
	if FAST:
		print "\nWarning! Running with '--fast' means that identical variants denoted differently between the two vcfs will not be detected! The values above may be lower than the true accuracy."
	print '\n**********************************\n'

	#if MAPTRACK != None:
	#	print 'mappability:'
	#	print mappability_vs_FN
	#
	#print 'coverage:'
	#for k in sorted(coverage_vs_FN.keys()):
	#	print k,':',coverage_vs_FN[k]
	#
	#print 'allele balance:'
	#for k in sorted(alleleBal_vs_FN.keys()):
	#	print k,':',alleleBal_vs_FN[k]


if __name__ == '__main__':
	main()