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
		return '{:.2f} kbp'.format(float(bp)/1000.)
	elif bp < 1000000000:
		return '{:.2f} Mbp'.format(float(bp)/1000000.)
	elif bp < 1000000000000:
		return '{:.2f} Gbp'.format(float(bp)/1000000000.)
	elif bp < 1000000000000000:
		return '{:.2f} Tbp'.format(float(bp)/1000000000000.)
	elif bp < 1000000000000000000:
		return '{:.2f} Pbp'.format(float(bp)/1000000000000000.)
	return str(bp)+' bp'
