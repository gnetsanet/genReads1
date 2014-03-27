#!/usr/bin/env python
# encoding: utf-8

""" *******************************************************************

misc.py

******************************************************************* """

import os
import re
import copy
import random
import numpy as np
import bisect

NUCL_NUM = {'a':0,'A':0,'c':1,'C':1,'g':2,'G':2,'t':3,'T':3}
NUM_NUCL = ['A','C','G','T']
NUM_NUCL_LOWER = ['a','c','g','t']
TO_UPPER       = {'a':'A','A':'A','c':'C','C':'C','g':'G','G':'G','t':'T','T':'T'}
TO_UPPER_COMP  = {'a':'T','A':'T','c':'G','C':'G','g':'C','G':'C','t':'A','T':'A'}

NUCL    = 4
ASCII_N = 78
ASCII_n = 110

for k in TO_UPPER.keys():
	TO_UPPER[ord(k)] = TO_UPPER[k]
	TO_UPPER_COMP[ord(k)] = TO_UPPER_COMP[k]
for k in NUCL_NUM.keys():
	NUCL_NUM[ord(k)] = NUCL_NUM[k]

def randEvent(cumP):
	r = random.random()
	return bisect.bisect(cumP,r)-1

# slow, don't use me
def randEvent_old(cumP):
	r = random.random()
	for i in range(len(cumP)-1):
		if r >= cumP[i] and r < cumP[i+1]:
			return i
	return len(cumP)-1

def randSequence(l,upper):
	outSeq = ''
	if upper:
		for i in xrange(l):
			outSeq += NUM_NUCL[random.randint(0,NUCL-1)]
	else:
		for i in xrange(l):
			outSeq += NUM_NUCL_LOWER[random.randint(0,NUCL-1)]
	return outSeq

def posmax(seq):
	m = seq[0]
	index = 0
	for i,x in enumerate(seq):
		if x > m:
			m = x
			index = i
	return index

def histoxy(data,bins=10,myMin=None,myMax=None):
	epsilon = 0.0001
	if myMax == None:
		ma = max(data)
	else:
		ma = myMax
	if myMin == None:
		mi = min(data)
	else:
		mi = myMin
	step = float(ma-mi)/bins
	x = np.zeros(bins*2)
	x[::2] = np.linspace(mi,ma-step,bins)+epsilon
	x[1::2] = np.linspace(mi+step,ma,bins)-epsilon
	y = np.zeros(bins*2)
	for value in data:
		if value == ma:
			y[-1] += 1
			y[-2] += 1
		elif value > ma or value < mi:
			pass
		else:
			ind = 2*int((value-mi)/step)
			y[ind] += 1
			y[ind+1] += 1
	return x, y

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

def printSizeNicely(b):
	if b < 1000:
		return str(b)+' B'
	elif b < 1000000:
		return '{0:.2f} kB'.format(float(b)/1000.)
	elif b < 1000000000:
		return '{0:.2f} MB'.format(float(b)/1000000.)
	elif b < 1000000000000:
		return '{0:.2f} GB'.format(float(b)/1000000000.)
	elif b < 1000000000000000:
		return '{0:.2f} TB'.format(float(b)/1000000000000.)
	elif b < 1000000000000000000:
		return '{0:.2f} PB'.format(float(b)/1000000000000000.)
	return str(b)+' B'

def needleman_wunsch(mismatch, gap, sequence1, sequence2):
	"""
	Calculate the minimum penalty alignment.
	"""

	# Get lengths of strings.
	n = len(sequence1)
	m = len(sequence2)

	# Make two-dimensional list for subproblem solutions.
	subproblems = [[0 for x in range(m+1)] for x in range(n+1)]

	# Fill in zeros on both dimensions with gap penalties.
	for i in range(n+1):
		subproblems[i][0] = i * gap

	for j in range(m+1):
		subproblems[0][j] = j * gap

	# Calculate subproblem solutions.
	for i in range(1, n+1):
		for j in range(1, m+1):
			case1 = subproblems[i-1][j-1]

			if sequence1[i-1] != sequence2[j-1]:
				case1 += mismatch

			case2 = subproblems[i-1][j] + gap
			case3 = subproblems[i][j-1] + gap

			subproblems[i][j] = min([case1, case2, case3])

	penalty = subproblems[n][m]

	# Backtrace to reconstruct optimal alignment.
	alignment1 = ""
	alignment2 = ""

	i = n
	j = m
	while i > 0 or j > 0:
		pos = subproblems[i][j]
		case1_match = subproblems[i-1][j-1]
		case1_mismatch = case1_match + mismatch
		case2 = subproblems[i-1][j] + gap
		case3 = subproblems[i][j-1] + gap

		if i > 0 and pos == case1_match:
			alignment1 = sequence1[i-1] + alignment1
			alignment2 = sequence2[j-1] + alignment2
			i -= 1
			j -= 1
		elif i > 0 and pos == case1_mismatch:
			alignment1 = sequence1[i-1] + alignment1
			alignment2 = sequence2[j-1] + alignment2
			i -= 1
			j -= 1
		elif i > 0 and pos == case2:
			alignment1 = sequence1[i-1] + alignment1
			alignment2 = '-' + alignment2
			i -= 1
		elif j > 0 and pos == case3:
			alignment1 = '-' + alignment1
			alignment2 = sequence2[j-1] + alignment2
			j -= 1

	return (penalty, alignment1, alignment2)
