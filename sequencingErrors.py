import random
from misc import randEvent, NUCL_NUM, NUM_NUCL

class SequencingErrors:
	def __init__(self, READLEN, cpSSE, cpSIE, SIE_RATE, SIE_INS_FREQ):

		self.RL    = READLEN
		self.cpSSE = cpSSE
		self.cpSIE = cpSIE
		self.sie_R = SIE_RATE
		self.sie_I = SIE_INS_FREQ

		self.numSSE_total = 0
		self.numSIE_total = 0

		self.prevSIE    = {1:-1, 2:-1}
		self.sie_errors = []
		self.rpos       = {1:0, 2:0}
		self.myHitInds  = []

	def newReads(self, r1pos, r2pos=None):
		self.prevSIE    = {1:-1, 2:-1}
		self.sie_errors = []
		self.rpos       = {1:r1pos, 2:r2pos}
		self.myHitInds  = []

	def introduceError(self, rid, r, j):
		if (random.random() < self.sie_R) and (j > self.prevSIE[rid]+1) and (j != 0) and (j != self.RL-1):
			if random.random() < self.sie_I:
				# insert INS error
				#
				#	the way this works is that we insert a het-del with AF=0% that gets "corrected" back to a "ref" allele that contains the SIE error.
				#
				sie = ['D',1,self.rpos[rid]+j,self.rpos[rid]+j+1,NUM_NUCL[randEvent(self.cpSIE)]]
			else:
				# insert DEL error
				sie = ['I',1,self.rpos[rid]+j,self.rpos[rid]+j+2,'X']
			self.myHitInds.append(sie)
			self.sie_errors.append(sie)

			self.numSIE_total += 1
			self.prevSIE[rid]  = j
		else:
			r[j] = NUM_NUCL[randEvent(self.cpSSE[NUCL_NUM[r[j]]])]
			self.numSSE_total += 1

