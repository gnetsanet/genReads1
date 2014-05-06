#!/usr/bin/env python
# encoding: utf-8

""" **************************************************

python alignAccuracy2.py alignment.sam blackList.p

************************************************** """

import sys
import bisect
import cPickle as pickle
import numpy as np

def main():

	f = open(sys.argv[1],'r')

	[badRegions] = pickle.load(open(sys.argv[2],'rb'))

	print badRegions
	brtemp = []
	totalBadNucleotides = 0
	for n in badRegions:
		totalBadNucleotides += n[1]-n[0]
		brtemp.append(n[0])
		brtemp.append(n[1])
	badRegions = brtemp

	nMatches = 0
	nTotal   = 0
	nUmapped = 0
	nInBadRegion = 0
	nBadAndUnmap = 0

	for line in f:
		if line[0] != '@':
			splt = line.split('\t')
			samFlags    = np.binary_repr(int(splt[1]),16)
			isUnmapped  = int(samFlags[-3])
			isFirstRead = int(samFlags[-7])
			cigar       = splt[5]
			nTotal += 1

			if isUnmapped:
				nUmapped += 1

			#print splt
			#print isFirstRead
			rpos   = int(splt[3])
			rnSplt = splt[0].split('-')
			if len(rnSplt) == 4:
				correctPos = [int(rnSplt[1]),int(rnSplt[2])]
			elif len(rnSplt) == 3:
				isFirstRead = 1
				correctPos  = [int(rnSplt[1])]
			if isFirstRead:
				if correctPos[0] == rpos:
					nMatches += 1
				else:
					dif  = rpos-correctPos[0]
					difS = str(dif)+'S'
					if dif > 0 and difS == cigar[:len(difS)]:
						nMatches += 1
					else:
						if bisect.bisect(badRegions,correctPos[0])%2 == 1 or bisect.bisect(badRegions,rpos)%2 == 1:
							nInBadRegion += 1
							if isUnmapped:
								nBadAndUnmap += 1
						print 'correct:',correctPos[0],'/ workflow:',rpos,'\t',cigar
						#pass
			else:
				if correctPos[1] == rpos:
					nMatches += 1
				else:
					dif  = rpos-correctPos[1]
					difS = str(dif)+'S'
					if dif > 0 and difS == cigar[:len(difS)]:
						nMatches += 1
					else:
						if bisect.bisect(badRegions,correctPos[1])%2 == 1 or bisect.bisect(badRegions,rpos)%2 == 1:
							nInBadRegion += 1
							if isUnmapped:
								nBadAndUnmap += 1
						print 'correct:',correctPos[0],'/ workflow:',rpos,'\t',cigar
						#pass
			#break

	f.close()

	print '\ndataset:\t',sys.argv[1]
	print '\nblacklist:\t',sys.argv[2]
	print '\nnon-unique nucl:',totalBadNucleotides
	print '\ntotal reads:\t',nTotal
	print 'total correct:\t',nMatches,'({0:.3f}%)'.format(100.*float(nMatches)/nTotal)
	print 'total wrong:\t',nTotal-nMatches,'({0:.3f}%)'.format(100.*float(nTotal-nMatches)/nTotal)
	print 'total unmapped:\t',nUmapped,'({0:.3f}%)'.format(100.*float(nUmapped)/nTotal)
	print 'mapped wrong:\t',nTotal-nMatches-nUmapped,'({0:.3f}%)'.format(100.*float(nTotal-nMatches-nUmapped)/nTotal)
	print '\nerroneously mapped reads in non-unique regions:\t',nInBadRegion,'({0:.3f}% of all wrong reads)'.format(100.*nInBadRegion/float(nTotal-nMatches)),'({0:.3f}% of all unmapped reads)'.format(100.*float(nBadAndUnmap)/nUmapped)
	print ''

if __name__ == '__main__':
	main()