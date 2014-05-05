#!/usr/bin/env python
# encoding: utf-8

""" **************************************************


************************************************** """

import sys
import os
import re
import time
#import random
#import cPickle as pickle
#import numpy as np
#import matplotlib.pyplot as mpl

MIN_SLEN = 50
MAX_SLEN = 50000

# probability of non-unique L-mers in a random (1/4 for each nucleotide) sequence of length G
def idealFreq(L,G):
	return 1.-((1.-(4**-L))**(G-L))

def randSeq(L):
	nucl = ['A','C','G','T']
	os   = [nucl[random.randint(0,3)] for n in xrange(L)]
	return ''.join(os)

def printBasesNicely(bp):
	if bp < 1000:
		return str(bp)+' bp'
	elif bp < 1000000:
		return '{0:.2f} kbp'.format(float(bp)/1000.)
	elif bp < 1000000000:
		return '{0:.2f} Mbp'.format(float(bp)/1000000.)
	elif bp < 1000000000000:
		return '{0:.2f} Gbp'.format(float(bp)/1000000000.)
	elif bp < 1000000000000000:
		return '{0:.2f} Tbp'.format(float(bp)/1000000000000.)
	elif bp < 1000000000000000000:
		return '{0:.2f} Pbp'.format(float(bp)/1000000000000000.)
	return str(bp)+' bp'

def str2re(s):
	splt = s.split('.')
	os = re.escape(splt[0])
	for i in xrange(1,len(splt)):
		os = os+r"."+re.escape(splt[i])
	return os

def areYouIn(sl,ss):
	for i in xrange(len(sl)-len(ss)+1):
		wrong = 0
		for j in xrange(len(ss)):
			if ss[j] != sl[i+j] and sl[i+j] != '.':
				wrong = 1
				break
		if not wrong:
			return True
	return False


def main():

	if len(sys.argv) != 6:
		print '\npython repFinding3.py ref.fa out_prefix editDistance jobID nJobs\n'
		print '\nsingle-job example:\n'
		print '\npython repFinding3.py ref.fa myOut 1 1 1\n'
		print '\nmulti-job example:\n'
		print '\npython repFinding3.py ref.fa myOut 1 1 3'
		print 'python repFinding3.py ref.fa myOut 1 2 3'
		print 'python repFinding3.py ref.fa myOut 1 3 3\n'
		exit(1)
	else:
		REF_FILE = sys.argv[1]
		OUT_FILE = sys.argv[2]
		EDIT_DIS = int(sys.argv[3])
		JOB_ID   = int(sys.argv[4])-1
		JOB_TOT  = int(sys.argv[5])

	f = open(REF_FILE,'r')
	fr = f.read()
	f.close()

	of = open(OUT_FILE+'_e'+str(EDIT_DIS)+'_loc_job'+str(JOB_ID+1)+'of'+str(JOB_TOT)+'.txt','w')

	refs = []
	for ref in fr.split('>'):
		if len(ref) > 0:
			splt = ref.split('\n')
			refs.append([splt[0],''.join(splt[1:]).upper()])
	del fr

	for [refName, ref] in refs:

		print 'ref:',refName
		REF_LEN = len(ref)#-ref.count('N')
		print printBasesNicely(REF_LEN)
		#ref = randSeq(REF_LEN)

		howManyNonUnique = {}	# [kmer_len]
		maxSoFar = 0
		for i in xrange(JOB_ID,REF_LEN-1,JOB_TOT):
			stillGoing = True
			kl = MIN_SLEN
			while True:
				if kl > MAX_SLEN or i+kl >= REF_LEN:
					break
				soi = ref[i:i+kl]
				if 'N' in soi:
					break

				if re.subn(soi,'',ref,2)[1] > 1:
					nonUnique = 1
				else:
					nonUnique = 0
					allPatterns = [soi]
					for j in xrange(EDIT_DIS):
						temp = []
						for n in allPatterns:
							for k in xrange(len(n)):
								if k > 0:
									snpPat = str(n[:k])+'.'+str(n[k+1:])
									if re.subn(str2re(snpPat),'',ref,2)[1] > 0 and not areYouIn(snpPat,soi):
										#print snpPat
										nonUnique = 1
										break
									insPat = str(n[:k])+'.'+str(n[k:])
									if j < EDIT_DIS-1:
										temp.append(snpPat)
										temp.append(insPat)
									if re.subn(str2re(insPat),'',ref,2)[1] > 0 and not areYouIn(insPat,soi):
										#print insPat
										nonUnique = 1
										break
							if nonUnique:
								break
						if nonUnique:
							break
						allPatterns = list(set(temp))

				if nonUnique:
					if kl in howManyNonUnique:
						howManyNonUnique[kl] += 1
					else:
						howManyNonUnique[kl] = 1
					if kl > maxSoFar:
						maxSoFar = kl
					of.write(str(i)+'\t'+str(kl)+'\n')
					kl += 1
				else:
					break

			print i, [kl,maxSoFar]

		print howManyNonUnique
		of.close()

		of = open(OUT_FILE+'_e'+str(EDIT_DIS)+'_count_job'+str(JOB_ID+1)+'of'+str(JOB_TOT)+'.txt','w')
		for k in sorted(howManyNonUnique.keys()):
			of.write(str(k)+'\t'+str(howManyNonUnique[k])+'\n')
		of.close()

		# we're only gonna work with first ref in fasta at the moment...
		break

if __name__ == '__main__':
	main()