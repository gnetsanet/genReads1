#!/usr/bin/env python
# encoding: utf-8
""" ////////////////////////////////////////////////////////////////////////////////
   ///                                                                          ///
  ///       genReads.py                                                        ///
 ///                                                                          ///
///////     Variant and read simulator for benchmarking NGS workflows        //////
   ///                                                                         ///
  ///       Written by:     Zach Stephens                                     ///
 ///        For:            DEPEND Research Group                            ///
///////     Date:           January 8, 2013                                 ///////
   ///      Contact:        zstephe2@illinois.edu                              ///
  ///                                                                         ///       
 ///                                                                         ///
/////////////////////////////////////////////////////////////////////////////// """

import os
import sys
import copy
import random
import re
import time
import bisect
import cPickle as pickle
import numpy as np
#import matplotlib.pyplot as mpl
import optparse

# absolute path to this script
SIM_PATH         = '/'.join(os.path.realpath(__file__).split('/')[:-1])+'/'

sys.path.append(SIM_PATH)
from misc	import *
from sequencingErrors	import SequencingErrors


"""//////////////////////////////////////////////////
////////////    PARSE INPUT ARGUMENTS    ////////////
//////////////////////////////////////////////////"""


DESC   = """%prog: Variant and read simulator for benchmarking NGS workflows."""
VERS   = 0.2

DEFAULT_COV = 10.							# default average coverage value
DEFAULT_RNG = None							# default RNG seed
DEFAULT_RLN = 100							# default read-length
DEFAULT_FLN = 250							# default mean fragment-length
DEFAULT_FSD = 10							# default fragment-length std
DEFAULT_WIN = 1000							# default sliding window size (for sampling reads / computing GC%)
DEFAULT_SER = 0.01							# default average sequencing error rate
DEFAULT_VRA = 0.00034						# default average rate of variant occurences
DEFAULT_HVR = 0.2							# default heterozygous variant rate
DEFAULT_ABM = 0.5							# default allele balance mean
DEFAULT_ABS = 0.1							# default allele balance std
DEFAULT_MBD = 95							# default minimum bed-file region size to consider for targeted sequencing
DEFAULT_BCV = 0								# default coverage in non-targeted regions
DEFAULT_QSM = 'models/qModel_AlvaroPhix.p'	# default quality-score model
DEFAULT_SSM = 'models/SSE_model1_pIRS.p'	# default sequencing-substituion-error model
DEFAULT_MUM = 'models/mutModel_default.p'	# default mutation model
DEFAULT_SIE = 0.01							# default sequencing-indel-error frequency
DEFAULT_SII = 0.4							# default frequency of insertion errors

PARSER = optparse.OptionParser('python %prog [options] -r <ref.fa> -o <out.prefix>',description=DESC,version="%prog v"+str(VERS))

PARSER.add_option('-p',    help='Generate paired-end reads [%default]',               dest='BOOL_PE', default=False,       action='store_true')
PARSER.add_option('-b',    help='Input bed file containing regions to sample from',   dest='BED',     default=None,        action='store',      metavar='<targets.bed>')
PARSER.add_option('-a',    help='Allele balance (% ref reads) mean [%default]',       dest='ABM',     default=DEFAULT_ABM, action='store',      metavar='<float>')
PARSER.add_option('-A',    help='Allele balance std [%default]',                      dest='ABS',     default=DEFAULT_ABS, action='store',      metavar='<float>')
PARSER.add_option('-B',    help='Minimum bed file region length [%default]',          dest='MBD',     default=DEFAULT_MBD, action='store',      metavar='<int>')
PARSER.add_option('-c',    help='Average coverage [%default]',                        dest='COV',     default=DEFAULT_COV, action='store',      metavar='<float>')
PARSER.add_option('-C',    help='Avg coverage in non-targeted regions [%default]',    dest='BCV',     default=DEFAULT_BCV, action='store',      metavar='<float>')
PARSER.add_option('-f',    help='Average fragment length [%default]',                 dest='FLEN',    default=DEFAULT_FLN, action='store',      metavar='<int>')
PARSER.add_option('-F',    help='Fragment length std [%default]',                     dest='FSTD',    default=DEFAULT_FSD, action='store',      metavar='<int>')
PARSER.add_option('-i',    help='P(SIE error | sequencing error occurs) [%default]',  dest='SIER',    default=DEFAULT_SIE, action='store',      metavar='<float>')
PARSER.add_option('-I',    help='P(insertion error | SIE error occurs) [%default]',   dest='SIEI',    default=DEFAULT_SII, action='store',      metavar='<float>')
PARSER.add_option('-l',    help='Read length [%default]',                             dest='RLEN',    default=DEFAULT_RLN, action='store',      metavar='<int>')
PARSER.add_option('-m',    help='Mutation model [%default]',                          dest='MUM',     default=DEFAULT_MUM, action='store',      metavar='<model.p>')
PARSER.add_option('-o',    help='Output filename prefix',                             dest='OUT',                          action='store',      metavar='<out.prefix>')
PARSER.add_option('-q',    help='Quality score model [%default]',                     dest='QSM',     default=DEFAULT_QSM, action='store',      metavar='<model.p>')
PARSER.add_option('-r',    help='Reference fasta',                                    dest='REF',                          action='store',      metavar='<ref.fa>')
PARSER.add_option('-R',    help='RNG seed value [%default]',                          dest='RNG',     default=DEFAULT_RNG, action='store',      metavar='<int>')
PARSER.add_option('-s',    help='SSE model [%default]',                               dest='SSM',     default=DEFAULT_SSM, action='store',      metavar='<model.p>')
PARSER.add_option('-S',    help='Average sequencing error rate [%default]',           dest='SER',     default=DEFAULT_SER, action='store',      metavar='<float>')
PARSER.add_option('-v',    help='Input VCF file',                                     dest='VCF',     default=None,        action='store',      metavar='<str>')
PARSER.add_option('-V',    help='Average variant occurence rate [%default]',          dest='VRA',     default=DEFAULT_VRA, action='store',      metavar='<float>')
PARSER.add_option('-w',    help='Window size for computing GC% [%default]',           dest='WIN',     default=DEFAULT_WIN, action='store',      metavar='<int>')
PARSER.add_option('-z',    help='Ratio of heterozygous variants [%default]',          dest='HVR',     default=DEFAULT_HVR, action='store',      metavar='<float>')
PARSER.add_option('--FA',  help='Output modified reference as fasta file [%default]', dest='BOOL_FA', default=False,       action='store_true')
PARSER.add_option('--SAM', help='Output correct alignment as SAM file [%default]',    dest='BOOL_SA', default=False,       action='store_true')
PARSER.add_option('--VCF', help='Output introduced variants as VCF file [%default]',  dest='BOOL_VC', default=False,       action='store_true')
PARSER.add_option('--TXT', help='Output simulation info as .txt file [%default]',     dest='BOOL_RI', default=False,       action='store_true')
PARSER.add_option('--job', help='Job IDs for generating reads in parallel',           dest='MULTI',                        action='store',      metavar='<int> <int>', nargs=2)

PARSER.add_option('--no-snp',     help="Don't introduce SNPs [%default]",                           dest='NO_SNP', default=False, action='store_true')
PARSER.add_option('--no-ind',     help="Don't introduce small indels [%default]",                   dest='NO_IND', default=False, action='store_true')
PARSER.add_option('--no-svs',     help="Don't introduce large structural variants [%default]",      dest='NO_SVS', default=False, action='store_true')
PARSER.add_option('--no-qscores', help="Don't compute quality scores (use dummy value) [%default]", dest='NO_QSC', default=False, action='store_true')
PARSER.add_option('--no-gcbias',  help="Don't account for GC% bias [%default]",                     dest='NO_GCB', default=False, action='store_true')

(OPTS,ARGS) = PARSER.parse_args()
#print OPTS

if len(sys.argv) <= 1:
	PARSER.print_help()
	exit(1)

REFERENCE        = OPTS.REF
INPUT_BED        = OPTS.BED
INPUT_VCF        = OPTS.VCF

PAIRED_END       = OPTS.BOOL_PE
SAVE_SAM         = OPTS.BOOL_SA
SAVE_VCF         = OPTS.BOOL_VC
SAVE_RUNINFO     = OPTS.BOOL_RI
SAVE_CORRECT_REF = OPTS.BOOL_FA

NATURAL_SNPS       = not(OPTS.NO_SNP)
NATURAL_INDELS     = not(OPTS.NO_IND)
NATURAL_SVS        = not(OPTS.NO_SVS)
SEQUENCING_QSCORES = not(OPTS.NO_QSC)
GC_BIAS            = not(OPTS.NO_GCB)

if OPTS.MULTI == None:
	MULTI_JOB  = False
	JOB_ID     = 1
	JOB_TOT    = 1
else:
	MULTI_JOB  = True
	JOB_ID     = int(OPTS.MULTI[0])
	JOB_TOT    = int(OPTS.MULTI[1])
	if JOB_ID > JOB_TOT:
		print 'Error: Invalid job ID'
		exit(1)

READLEN        = int(OPTS.RLEN)
FRAGMENT_SIZE  = int(OPTS.FLEN)
FRAGMENT_STD   = int(OPTS.FSTD)
COV_WINDOW     = int(OPTS.WIN)
MIN_PROBE_LEN  = int(OPTS.MBD)

# if a sequencing error occurs, what are the odds it's an indel?
SIE_RATE = float(OPTS.SIER)
# if a sequencing indel error occurs, what are the odds it's an insertion as opposed to a deletion?
SIE_INS_FREQ = float(OPTS.SIEI)
# if a sequencing insertion error occurs, what's the probability of it being an A, C, G, T...
SIE_NUCL = [0.25, 0.25, 0.25, 0.25]

AVG_SER        = float(OPTS.SER)
if AVG_SER != 0.0:
	SEQUENCING_SNPS    = 1
	if SIE_RATE != 0.0:
		SEQUENCING_INDELS = 1
	else:
		SEQUENCING_INDELS = 0
else:
	SEQUENCING_SNPS    = 0
	SEQUENCING_INDELS  = 0
	SEQUENCING_QSCORES = 0

HET_VAR        = float(OPTS.HVR)
AB_MEAN        = float(OPTS.ABM)
AB_STD         = float(OPTS.ABS)

# consider fragment sizes within +- 3 stds
if PAIRED_END:
	MIN_SAMPLE_SIZE = FRAGMENT_SIZE+3*FRAGMENT_STD
else:
	MIN_SAMPLE_SIZE = READLEN*2

if FRAGMENT_STD == 0:
	FRAGLEN         = [FRAGMENT_SIZE]
	FRAGPROB        = [1.]
else:
	FRAGLEN         = range(FRAGMENT_SIZE-3*FRAGMENT_STD,FRAGMENT_SIZE+3*FRAGMENT_STD+1)
	FRAGPROB        = [np.exp(-(((n-float(FRAGMENT_SIZE))**2)/(2*(FRAGMENT_STD**2)))) for n in FRAGLEN]

GC_COVERAGE   = [n/100. for n in range(101)]
if GC_BIAS:
	RELATIVE_GC_COVERAGE_BIAS = [np.exp(-(((n-0.5)**2)/(2*(0.15**2)))) for n in GC_COVERAGE]
else:
	RELATIVE_GC_COVERAGE_BIAS = [1. for n in GC_COVERAGE]

if AB_STD > 0.0:
	AB      = np.linspace(AB_MEAN-3*AB_STD,AB_MEAN+3*AB_STD,20).tolist()
	AB_PROB = [np.exp(-(((n-AB_MEAN)**2)/(2*(AB_STD**2)))) for n in AB]
	AB_PROB = [n/sum(AB_PROB) for n in AB_PROB]
else:
	AB      = [AB_MEAN]
	AB_PROB = [1.0]

if OPTS.RNG == None:
	RNG_SEED   = random.randint(1,99999999)
else:
	RNG_SEED   = int(OPTS.RNG)
AVG_COVERAGE   = float(OPTS.COV)
AVG_VAR_FREQ   = float(OPTS.VRA)
BED_COVERAGE   = float(OPTS.BCV)

if MULTI_JOB:
	OUTFILE_NAME = OPTS.OUT+'_job'+str(JOB_ID)
else:
	OUTFILE_NAME = OPTS.OUT

if OPTS.MUM == DEFAULT_MUM:
	MUT_MODEL = SIM_PATH+DEFAULT_MUM
else:
	MUT_MODEL = OPTS.MUM


"""//////////////////////////////////////////////////
////////////    OTHER MISC PARAMETERS    ////////////
//////////////////////////////////////////////////"""

# how much random leeway do we have with SNP_FREQ and INDEL_FREQ?
#	e.g. if RNG_MULT = 0.1, during experiments the actual probability
#   of a SNP occuring could be between 0.9*SNP_FREQ and 1.1*SNP_FREQ
RNG_MULT   = .05
# like RNG_MULT, but with coverage
if GC_BIAS:
	COV_MULT   = .1
else:
	COV_MULT   = .01

# dummy quality score if no model is used
DUMMY_QSCORE = 30


"""////////////////////////////////////////////////////////
////////////    LOAD SEQUENCING ERROR MODEL    ////////////
////////////////////////////////////////////////////////"""

## how are the nucleotides misrepresented if a sub sequencing error occurs?
#SEQ_ERR    = [[0.,     0.4918, 0.3377, 0.1705 ],
#			  [0.5238,     0., 0.2661, 0.2101 ],
#			  [0.3754, 0.2355,     0., 0.3890 ],
#			  [0.2505, 0.2552, 0.4942, 0.     ]]
#
## probability of sequencing insertion errors (length 1,2,3..)
#SEQ_INS    = [0., 0., 0.]
#
## probability of sequencing deletion errors (length 1,2,3..)
#SEQ_DEL    = [0., 0., 0.]

if OPTS.SSM == DEFAULT_SSM:
	SSE_MODEL = SIM_PATH+DEFAULT_SSM
else:
	SSE_MODEL = OPTS.SSM


"""////////////////////////////////////////////////
////////////    QUALITY-SCORE MODEL    ////////////
////////////////////////////////////////////////"""

if OPTS.QSM == DEFAULT_QSM:
	QSCORE_MODEL = SIM_PATH+DEFAULT_QSM
else:
	QSCORE_MODEL = OPTS.QSM


"""////////////////////////////////////////////////
////////////      SANITY-CHECKING      ////////////
////////////////////////////////////////////////"""

if OPTS.REF == None:
	print 'Error: No reference provided.'
	exit(1)

if OPTS.OUT == None:
	print 'Error: No outfile prefix provided.'
	exit(1)

if AVG_COVERAGE <= 0.:
	print 'Error: Average coverage must be greater than 0.'
	exit(1)
if AVG_VAR_FREQ < 0. or AVG_VAR_FREQ >= 1.:
	print 'Error: Average variant frequency must be between [0,1)'
	exit(1)
if AVG_SER < 0. or AVG_SER >= 1.:
	print 'Error: Average sequencing error rate (SER) rate must be between [0,1)'
	exit(1)
if SIE_RATE < 0. or SIE_RATE >= 1.:
	print 'Error: Sequencing indel error (SIE) frequency must be between [0,1)'
	exit(1)
if SIE_INS_FREQ < 0. or SIE_INS_FREQ >= 1.:
	print 'Error: Sequencing insertion error frequency must be between [0,1)'
	exit(1)
if HET_VAR < 0. or HET_VAR >= 1.:
	print 'Error: Heterozygous variant ratio must be between [0,1)'
	exit(1)
if READLEN <= 0:
	print 'Error: Read length must be greater than 0.'
	exit(1)
if FRAGMENT_SIZE <= 0:
	print 'Error: Mean fragment length must be greater than 0.'
	exit(1)
if FRAGMENT_STD < 0:
	print 'Error: Fragment length standard deviation must be non-negative.'
	exit(1)
if COV_WINDOW <= 0:
	print 'Error: Sliding sampling window must be greater than 0.'
	exit(1)

if AB_MEAN < 0. or AB_MEAN >= 1.:
	print 'Error: AB_MEAN must be between [0,1)'
	exit(1)
if AB_STD < 0.:
	print 'Error: AB_STD must be non-negative.'
	exit(1)
if AB_MEAN - 3*AB_STD < 0.:
	print 'Error: AB_MEAN - 3*AB_STD must be greater than 0.'
	exit(1)
if AB_MEAN + 3*AB_STD >= 1.:
	print 'Error: AB_MEAN + 3*AB_STD must be less than 1.'
	exit(1)

if PAIRED_END:
	if 2*READLEN >= FRAGMENT_SIZE-3*FRAGMENT_STD:
		print 'Error: Fragment length must be at least 2*READLEN + 3*FRAGMENT_STD'
		exit(1)
	if COV_WINDOW <= FRAGMENT_SIZE+3*FRAGMENT_STD:
		print 'Error: Window size must be at least FRAGMENT_LEN + 3*FRAGMENT_STD'
		exit(1)
else:
	if READLEN >= COV_WINDOW:
		print 'Error: Window size must be at least READLEN'

if not os.path.isfile(REFERENCE):
	print 'Error: Could not open reference file.'
	exit(1)

if not os.path.isfile(QSCORE_MODEL):
	print 'Error: Could not open quality-score model.'
	exit(1)
else:
	maxQScore = pickle.load(open(QSCORE_MODEL,'rb'))[2][-1]
	maxPerror = 10.**(-maxQScore/10.)
	if (AVG_SER > 1.0 or AVG_SER <= maxPerror) and AVG_SER != 0.:
		print 'Error: Average SSE rate must be between: ({0:.6f}, 1.0], or == 0'.format(maxPerror)
		exit(1)

if BED_COVERAGE < 0.:
	print 'Error: Coverage in non-targeted regions must be non-negative.'
	exit(1)

if not os.path.isfile(SSE_MODEL):
	print 'Error: Could not open SSE model.'
	exit(1)

if not os.path.isfile(MUT_MODEL):
	print 'Error: Could not open mutation model.'
	exit(1)


"""////////////////////////////////////////////////
////////////    LOAD MUTATION MODEL    ////////////
////////////////////////////////////////////////"""

## name of model
#MUT_MNAME    = 'DEFAULT'
## probability of SNP occuring at any nucleotide
#SNP_FREQ   = .003
## probability of indels of length 1,2,3,4,... occuring at any nucleotide
#INDEL_MUT  = [.00015,.00002,.00001,.000005,.000002]
## if an indel occurs, what are the odds it's an insertion as opposed to a deletion?
#INS_FREQ   = 0.25
## how do the nucleotides mutate into eachother if a SNP occurs?
#SNP_MUT    = [[0.,   0.15,  0.70,  0.15],
#		      [0.15, 0.,    0.15,  0.70],
#		      [0.70, 0.15,  0.,    0.15],
#		      [0.15, 0.70,  0.15,  0.  ]]

[MUT_MNAME, SNP_FREQ, INDEL_MUT, INS_FREQ, SNP_MUT] = pickle.load(open(MUT_MODEL,'rb'))

# max indel size
MAX_INDEL  = len(INDEL_MUT)
# probability of any indel occuring
INDEL_FREQ = sum(INDEL_MUT)

if SNP_FREQ <= 0.0 or SNP_FREQ >= 1.0:
	NATURAL_SNPS = False
if INDEL_FREQ <= 0.0 or INDEL_FREQ >= 1.0:
	NATURAL_INDELS = False
	INDEL_MUT  = [0.1]		# junk placeholder values that won't be used, but needed to prevent div/0 errors
	INDEL_FREQ = 0.1		#


"""////////////////////////////////////////////////
////////////           MAIN()          ////////////
////////////////////////////////////////////////"""

def main():

	startTime = time.time()
	startDate = time.asctime()

	if RNG_SEED != None:
		random.seed(RNG_SEED)

	# create cumulative probability list for indel lengths
	cpIndel = [0.]+[n/INDEL_FREQ for n in INDEL_MUT][:-1]
	for i in range(1,len(INDEL_MUT)):
		cpIndel[i] = cpIndel[i]+cpIndel[i-1]
	#print cpIndel
	
	# create cumulative probability lists for SNPs
	cpSNP = copy.deepcopy(SNP_MUT)
	for i in range(len(cpSNP)):
		cpSNP[i] = [0]+cpSNP[i][:-1]
		for j in range(1,len(cpSNP[i])):
			cpSNP[i][j] = cpSNP[i][j]+cpSNP[i][j-1]
	#print cpSNP

	# create cumulative probability lists for substitution sequencing errors
	[SEQ_ERR] = pickle.load(open(SSE_MODEL,'rb'))
	cpSSE = copy.deepcopy(SEQ_ERR)
	for i in range(len(cpSSE)):
		cpSSE[i] = [0]+cpSSE[i][:-1]
		for j in range(1,len(cpSSE[i])):
			cpSSE[i][j] = cpSSE[i][j]+cpSSE[i][j-1]
	#print cpSSE

	# create cumulative probability list for allele balance
	cpAB = np.cumsum(AB_PROB).tolist()[:-1]
	cpAB.insert(0,0.)

	# create cumulative probability list for sequencing insertion error nucleotide
	cpSIE = np.cumsum(SIE_NUCL).tolist()[:-1]
	cpSIE.insert(0,0.)

	# load empirical quality score distributions
	[allPMFs,probQ,Qscores,qOffset] = pickle.load(open(QSCORE_MODEL,'rb'))
	if len(allPMFs)-1 < READLEN:
		# we need to extend our qscore model to facilitate longer reads
		nToAdd = READLEN - len(allPMFs) + 1
		nAdded = 0
		while True:
			#print len(allPMFs)
			adj = 0
			for i in xrange(1,len(allPMFs)-1):
				toIns = [(allPMFs[i+adj][j]+allPMFs[i+adj+1][j])/2. for j in xrange(len(allPMFs[0]))]
				toIns = [n/sum(toIns) for n in toIns]
				allPMFs.insert(i+adj,copy.deepcopy(toIns))
				nAdded += 1
				adj += 1
				if nAdded == nToAdd:
					break
			if nAdded == nToAdd:
				break

	if SEQUENCING_QSCORES:
		learnedCycles = len(allPMFs)
		nQscores = len(allPMFs[0])

		targetQscore = int(-10.*np.log10(AVG_SER))
		print '\nDesired P(error):\t\t','{0:.6f}'.format(AVG_SER),'( phred:',targetQscore,')'
		qScoreShift = 0
		prevShift   = 0
		pprevShift  = 0
		printFirst  = 1
		while True:
			initQpmf = copy.deepcopy(allPMFs[0])
			#print initQpmf
			if qScoreShift > 0:
				for k in xrange(qScoreShift):
					del initQpmf[-1]
					initQpmf.insert(0,0.)
			elif qScoreShift < 0:
				for k in xrange(-qScoreShift):
					del initQpmf[0]
					initQpmf.append(0.)
			#print initQpmf
			siqp = sum(initQpmf)
			initQpmf = [n/siqp for n in initQpmf]

			initQProb = np.cumsum(initQpmf).tolist()[:-1]
			initQProb.insert(0,0.)
			allQscoreCumProb = [0 for n in xrange(learnedCycles*nQscores)]

			for i in xrange(learnedCycles):		# for each cycle 'i'...
				for j in xrange(nQscores):		# for each previous qscore 'j'...
					prod = [allPMFs[i][k]*probQ[j][k] for k in xrange(nQscores)]
					#prod = [1.*probQ[j][k] for k in xrange(nQscores)]
					if qScoreShift > 0:
						for k in xrange(qScoreShift):
							del prod[-1]
							prod.insert(0,0.)
					elif qScoreShift < 0:
						for k in xrange(-qScoreShift):
							del prod[0]
							prod.append(0.)
					sprd = sum(prod)
					prod = [n/sprd for n in prod]
					cumProb = np.cumsum(prod).tolist()[:-1]
					cumProb.insert(0,0.)
					allQscoreCumProb[i*nQscores+j] = cumProb

			pError = {}
			for i in xrange(nQscores):
				ind = bytearray([i+qOffset])[0]
				pError[ind] = 10.**(-i/10.)
			# scale error rates such that the average SSE rate is as desired
			# this is done by sampling a number of reads and estimating the current average rate
			sseScaleSample = 2000000/READLEN
			avgError = []
			for i in xrange(sseScaleSample):
				q1 = bytearray()
				q1.append(randEvent(initQProb)+qOffset)
				for j in xrange(1,READLEN):
					q1.append(randEvent(allQscoreCumProb[j*nQscores+q1[j-1]])+qOffset)
				#print [n-qOffset for n in q1]
				avgError.append(sum([pError[n] for n in q1])/READLEN)
			modelError  = (sum(avgError)/len(avgError))
			qScoreShift += int(-10.*np.log10(AVG_SER)) - int(-10.*np.log10(modelError))
			sseScalar = AVG_SER/modelError
			if printFirst:
				printFirst = 0
				print 'Avg P(error) of Input Model:\t','{0:.6f}'.format(modelError),'( phred:',int(-10.*np.log10(modelError)),')'

			if qScoreShift == prevShift or (qScoreShift == pprevShift and prevShift > qScoreShift):
				print 'Calibrated P(error):\t\t','{0:.6f}'.format(modelError),'( phred:',int(-10.*np.log10(modelError)),')'
				break

			pprevShift = prevShift
			prevShift  = qScoreShift

	else:
		plusMinusSSE = 0.5
		positionSSErate = np.linspace((1-plusMinusSSE)*AVG_SER,(1+plusMinusSSE)*AVG_SER,READLEN).tolist()

	# create cumulative probability list for fragment length
	cpFrag = [float(n)/sum(FRAGPROB) for n in FRAGPROB]
	cpFrag = np.cumsum(cpFrag).tolist()[:-1]
	cpFrag.insert(0,0.)

	# initialize class used to introduce sequencing errors
	seqErrObj = SequencingErrors(READLEN, cpSSE, cpSIE, SIE_RATE, SIE_INS_FREQ)


	"""//////////////////////////////////////////////////////
	////////////      INDEX REFERENCE FASTA      ////////////
	//////////////////////////////////////////////////////"""

	ref = []
	REFFILE = open(REFERENCE,'r')
	nLines = 0
	prevR = None
	prevP = None
	ref_inds = []
	sys.stdout.write('\nindexing reference fasta...')
	sys.stdout.flush()
	tt = time.time()
	while 1:
		nLines += 1
		data = REFFILE.readline()
		if not data:
			ref_inds.append( (prevR, prevP, REFFILE.tell()-len(data)) )
			break
		if data[0] == '>':
			if prevP != None:
				ref_inds.append( (prevR, prevP, REFFILE.tell()-len(data)) )

			prevP = REFFILE.tell()
			prevR = data[1:-1]
	print '{0:.3f} (sec)\n'.format(time.time()-tt)


	"""//////////////////////////////////////////////////////////
	////////////   ITERATE THROUGH EACH SEQ IN REF   ////////////
	//////////////////////////////////////////////////////////"""

	# initialize output files
	sys.stdout.write('initializing output files... ')
	sys.stdout.flush()
	tt = time.time()
	OUTFQ1 = open(OUTFILE_NAME+'_read1.fq','wb')
	if PAIRED_END:
		OUTFQ2 = open(OUTFILE_NAME+'_read2.fq','wb')
	if SAVE_SAM:
		OUTSAM = open(OUTFILE_NAME+'_golden.sam','wb')
		OUTSAM.write('@HD\tVN:1.0\tSO:unsorted\n')
		for n_RI in ref_inds:
			refName = n_RI[0]
			REFFILE.seek(n_RI[1])
			myDat   = REFFILE.read(n_RI[2]-n_RI[1]).split('\n')
			myLen   = sum([len(m) for m in myDat])
			OUTSAM.write('@SQ\tSN:'+refName+'\tLN:'+str(myLen)+'\n')
	if SAVE_VCF:
		OUTVCF = open(OUTFILE_NAME+'_golden.vcf','wb')
		OUTVCF.write('##fileformat=VCFv4.1\n')
		OUTVCF.write('##reference='+REFERENCE+'\n')
		OUTVCF.write('##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">\n')
		OUTVCF.write('##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">\n')
		OUTVCF.write('##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">\n')
		OUTVCF.write('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n')
		OUTVCF.write('##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">\n')
		OUTVCF.write('##INFO=<ID=TRANSREF,Number=1,Type=String,Description="Name of reference sequence where inserted bases of translocated sequence originate">\n')
		OUTVCF.write('##INFO=<ID=TRANSPOS,Number=.,Type=Integer,Description="Position within TRANSREF where the inserted bases of translocated sequence originate">\n')
		OUTVCF.write('##ALT=<ID=DEL,Description="Deletion">\n')
		OUTVCF.write('##ALT=<ID=DUP,Description="Duplication">\n')
		OUTVCF.write('##ALT=<ID=INS,Description="Insertion of novel sequence">\n')
		OUTVCF.write('##ALT=<ID=INV,Description="Inversion">\n')
		OUTVCF.write('##ALT=<ID=CNV,Description="Copy number variable region">\n')
		OUTVCF.write('##ALT=<ID=TRANS,Description="Translocation">\n')
		OUTVCF.write('##ALT=<ID=INV-TRANS,Description="Inverted translocation">\n')
		OUTVCF.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
	print '{0:.3f} (sec)\n'.format(time.time()-tt)


	bigReadNameOffset    = 0
	totalSNPs            = 0
	totalInds            = 0
	totalSTVs            = 0
	sequencesSampledFrom = 0
	totalBPtargeted      = 0
	totalBedTargetedBP   = 0
	readsFromThisJob     = 0
	# for each sequence in reference fasta file...
	for n_RI in ref_inds:
		refName = n_RI[0]


		"""//////////////////////////////////////////////////////
		////////////     READ REFERENCE SEQUENCE     ////////////
		//////////////////////////////////////////////////////"""

		# assumes fasta file is sane and has '\n' characters within long sequences
		REFFILE.seek(n_RI[1])
		sys.stdout.write('\n********************************\n\nreading '+refName+'... ')
		sys.stdout.flush()
		tt = time.time()
		myDat  = REFFILE.read(n_RI[2]-n_RI[1]).split('\n')
		myLens = [len(m) for m in myDat]
		myLen  = sum(myLens)
		print '{0:.3f} (sec)'.format(time.time()-tt)
		if sys.version_info >= (2,7):
			print '{:,} bp\n'.format(myLen)
		else:
			print '{0:} bp\n'.format(myLen)
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


		"""////////////////////////////////////////////////////////
		////////////    PRE-COMPUTE INDEL LOCATIONS    ////////////
		////////////////////////////////////////////////////////"""

		if NATURAL_INDELS or NATURAL_SNPS or NATURAL_SVS:
			sys.stdout.write('introducing variants... ')
			sys.stdout.flush()
			indelStart = time.time()

		#keepOriginal = bytearray(''.join(myDat))
		originalLen  = myLen
		if RNG_SEED != None:
			random.seed(RNG_SEED)

		cumLens = np.cumsum(myLens)
		lenAdj  = 0
		del myLens

		# normalize variant frequencies
		var_mult = AVG_VAR_FREQ/(INDEL_FREQ+SNP_FREQ)

		# introduce simulated indels
		indelsToAttempt = []
		nIndels = 0
		if NATURAL_INDELS:
			rIF = ((random.random()-.5)*(2*RNG_MULT) + 1.)*INDEL_FREQ*var_mult
			nIndels = int(rIF*myLen+.5)
			for i in xrange(nIndels):
				p1 = random.randint(0,len(myDat)-1)
				p2 = random.randint(0,len(myDat[p1])-1)
				idlen   = randEvent(cpIndel)+1	# rawr
				if random.random() < INS_FREQ:
					idT = 'I'
					# insertion:
					#
					#  ...p2-1, p2, ][, p2+1, p2+2, ...
				else:
					idT = 'D'
					# deletion:
					#
					#  ...p2-1, p2, [], p2+idlen+1, p2+idlen+2, ...
				indelsToAttempt.append((p1,p2,idlen,idT))


		SVsToAttempt = []
		# introduce simulated structural variants
		if NATURAL_SVS:

			"""
			SVsToAttempt.append(([5000,1000],'BD'))
			SVsToAttempt.append(([20000,''.join(['AG']*250)],'BI'))			
			SVsToAttempt.append(([40000,'ACGTACGTTTGA',3],'R'))
			SVsToAttempt.append(([60000,1000,80000],'T'))
			SVsToAttempt.append(([90000,50],'IV'))
			SVsToAttempt.append(([120000,200,110000],'IT'))
			SVsToAttempt.append(([130000,50],'BD'))
			"""
			pass


		# introduce variants from provided input VCF file
		input_snps    = []
		input_snps_AF = []
		input_inds    = {}
		input_inds_AF = {}
		myDatCat      = bytearray(''.join(myDat))
		if INPUT_VCF != None:
			invcf = open(INPUT_VCF,'r')
			for line in invcf:
				if line[0] != '#':
					splt = line.split('\t')
					if splt[0] == refName:
						pos = int(splt[1])
						rnt = splt[3].upper()
						ant = splt[4].upper()
						inf = ';'+splt[7]+';'
						# structural variants aren't supported yet, sorry!
						if ('[' in ant) or (']' in ant) or (':' in ant) or ('SVTYPE' in splt[7]):
							print "skipping variant [invalid symbols]:",line
							continue
						# if multiple alternate alleles, lets just use first one...
						if ',' in ant:
							ant = ant.split(',')[0]
						# lets use the AF field, if present, to force a variant to be heterozygous
						myAF = 1.
						if ';AF=' in inf:
							myAF = float(re.findall(r";AF=.*?(?=;)",inf)[0][3:].split(',')[0])
							if myAF < 0.01 or myAF > 0.99:
								myAF = 1.
						# snps
						if len(rnt) == len(ant):
							for i in xrange(len(rnt)):
								if chr(myDatCat[pos+i-1]).upper() != rnt[i]:
									print "skipping variant [!=REF]:",line
									continue
								input_snps.append((pos+i,rnt[i],ant[i]))
								input_snps_AF.append(myAF)
						# insertion
						elif len(rnt) == 1 and len(ant) > 1:
							if chr(myDatCat[pos+1]).upper() != rnt:
								print "skipping variant [!=REF]:",line
								continue
							v = ((pos,ant[1:]),'BI')
							SVsToAttempt.append(v)
							input_inds[(pos,rnt,ant)] = 1
							input_inds_AF[tuple(v)] = myAF
						# deletion
						elif len(rnt) > 1 and len(ant) == 1:
							if myDatCat[pos+1:pos+len(rnt)-1].upper() != rnt:
								print myDatCat[pos+1:pos+len(rnt)+1].upper(), rnt
								print "skipping variant [!=REF]:",line
								exit(1)
								continue
							v = ((pos,len(rnt)-1),'BD')
							SVsToAttempt.append(v)
							input_inds[(pos,rnt,ant)] = 1
							input_inds_AF[tuple(v)] = myAF
						# otherwise I have no idea, and don't feel like figuring it out now.
						else:
							print "skipping variant [confused]:",line
							continue
			invcf.close()
		del myDatCat
		nIndels += len(input_inds)


		"""//////////////////////////////////////////////////////////////////
		////////////    TRANSLATE COMPLEX SVS INTO BIG INDELS    ////////////
		//////////////////////////////////////////////////////////////////"""

		# idT:
		#
		#	'I'  : insertion
		#	'D'  : deletion
		#	'BI' : big insertion
		#	'BD' : big deletion
		#	'R'  : repeat
		#	'T'  : translocation
		#	'IV' : inversion
		#	'IT' : inverted translocation

		SVregions = []	# disallow small "normal" indels from happening within regions affected by SVs
		nSVs = -len(input_inds)
		indelsToAttempt_AF = {}		# key: indel variant tuple, value: allele freq
		for i in xrange(len(SVsToAttempt)):
			n = SVsToAttempt[i]

			# some elaborate code goes here to prevent overlapping events
			# TBD

			# check if event is in bounds
			# TBD

			# check if event involves any Ns
			# TBD

			nSVs += 1
			idT = n[1]
			Svv = n[0]

			p1 = Svv[0]/inWidth
			p2 = Svv[0]%inWidth

			if idT == 'T' or idT == 'IV' or idT == 'IT':
				hungryHungry = 0	# feed me nucleotides!
				hp1 = p1
				hp2 = p2
				seq = bytearray()
				while hungryHungry < Svv[1]:
					if len(myDat[hp1][hp2:]) <= Svv[1]-hungryHungry:
						hungryHungry += len(myDat[hp1][hp2:])
						seq.extend(myDat[hp1][hp2:])
						hp1 += 1
						hp2 = 0
					else:
						seq.extend(myDat[hp1][hp2:hp2+Svv[1]-hungryHungry])
						break
				seq = str(seq)

			if idT == 'BD':
				v = (p1,p2,Svv[1],'BD',i,'')
				indelsToAttempt.append(v)
				SVregions.append((Svv[0],Svv[0]+Svv[1]))
				# if event happens to be heterozygous, lets pass that information along in another list
				if tuple(n) in input_inds_AF:
					indelsToAttempt_AF[tuple(v)] = input_inds_AF[tuple(n)]

			elif idT == 'BI':
				v = (p1,p2,len(Svv[1]),'BI',i,Svv[1])
				indelsToAttempt.append(v)
				SVregions.append((Svv[0],Svv[0]+len(Svv[1])))
				if tuple(n) in input_inds_AF:
					indelsToAttempt_AF[tuple(v)] = input_inds_AF[tuple(n)]

			elif idT == 'R':
				indelsToAttempt.append((p1,p2,len(Svv[1])*Svv[2],'BI',i,''.join([Svv[1]]*Svv[2])))
				SVregions.append((Svv[0],Svv[0]+len(Svv[1])*Svv[2]))

			elif idT == 'T':
				indelsToAttempt.append((p1,p2,Svv[1],'BD',i,''))
				p3 = Svv[2]/inWidth
				p4 = Svv[2]%inWidth
				indelsToAttempt.append((p3,p4,Svv[1],'BI',i,seq))
				SVregions.append((Svv[0],Svv[0]+Svv[1]))
				SVregions.append((Svv[2],Svv[2]+Svv[1]))

			elif idT == 'IV':
				newSeq = ''.join([TO_UPPER_COMP[n] for n in seq[::-1]])
				indelsToAttempt.append((p1,p2,Svv[1],'BI',i,newSeq))
				indelsToAttempt.append((p1,p2,Svv[1],'BD',i,'**inv**'))
				SVregions.append((Svv[0],Svv[0]+Svv[1]))

			elif idT == 'IT':
				newSeq = ''.join([TO_UPPER_COMP[n] for n in seq[::-1]])
				indelsToAttempt.append((p1,p2,Svv[1],'BD',i,''))
				p3 = Svv[2]/inWidth
				p4 = Svv[2]%inWidth
				indelsToAttempt.append((p3,p4,Svv[1],'BI',i,newSeq))
				SVregions.append((Svv[0],Svv[0]+Svv[1]))
				SVregions.append((Svv[2],Svv[2]+Svv[1]))
		del input_inds_AF


		"""/////////////////////////////////////////////////////////////
		////////////    INTRODUCE INSERTIONS & DELETIONS    ////////////
		/////////////////////////////////////////////////////////////"""

		indelsToAttempt = sorted(indelsToAttempt)
		for i in xrange(1,len(indelsToAttempt)):
			if indelsToAttempt[i][:3] == indelsToAttempt[i-1][:3] and indelsToAttempt[i][3] == 'BI' and indelsToAttempt[i-1][3] == 'BD':
				(indelsToAttempt[i],indelsToAttempt[i-1]) = (indelsToAttempt[i-1],indelsToAttempt[i])

		prevP1       = -1
		prevP2       = -1
		prevP2_adj   = 0
		svInds       = {}
		indelList    = []
		afterThis    = [-1]
		addThis      = [0]
		indelList_AF = {}
		for n in indelsToAttempt:

			idT   = n[3]
			svInd = -1
			if idT == 'I' or idT == 'D':
				(p1,p2,idlen,idT) = n
			else:
				(p1,p2,idlen,idT,svInd,insSeq) = n

			if p1 != prevP1:
				prevP2 = -1
				prevP2_adj = 0
			if idT != 'BD' and idT != 'BI':
				if len(myDat[p1]) <= p2 or len(myDat[p1]) <= MAX_INDEL:
					continue
				if myDat[p1][p2] == 'N':
					continue
				if p2-prevP2 <= abs(MAX_INDEL):	# disallow overlapping indels
					continue

			if p1 == 0:
				ind = p2 + lenAdj
			else:
				ind = cumLens[p1-1]+p2+lenAdj

			if idT == 'I':
				#seq  = randSequence(idlen,myDat[p1][p2] in NUM_NUCL)
				seq  = randSequence(idlen,1)	# always insert upper-case
				if p2 == 0:
					myDat[p1] = seq + myDat[p1]
				else:
					myDat[p1] = myDat[p1][:p2] + seq + myDat[p1][p2:]
				lenAdj += idlen
				indBp = ind+idlen+1

			elif idT == 'BI':
				if len(insSeq) == 0:
					seq = randSequence(idlen,1)
				else:
					seq = insSeq
				if p2 == 0:
					myDat[p1] = seq + myDat[p1]
				else:
					myDat[p1] = myDat[p1][:p2+prevP2_adj] + seq + myDat[p1][p2+prevP2_adj:]
				idT = 'I'
				lenAdj += idlen
				prevP2_adj += idlen
				indBp = ind+idlen+1

			elif idT == 'D':
				if p2+idlen > len(myDat[p1]):
					# if we need to delete more bases than we have in this substring, let's cheat and shift it back..
					delBp = p2+idlen-len(myDat[p1])
					p2  -= delBp
					ind -= delBp
					if p2-prevP2 <= abs(MAX_INDEL):	# disallow overlapping indels
						continue

				seq = myDat[p1][p2:p2+idlen]
				myDat[p1] = myDat[p1][:p2]+myDat[p1][p2+idlen:]
				#myLens[p1] -= idlen
				lenAdj -= idlen
				indBp = ind+1

			elif idT == 'BD':

				hungryHungry = 0	# feed me nucleotides!
				hp1 = p1
				hp2 = p2+prevP2_adj
				if insSeq == '**inv**':
					hp2 += idlen
				seq = bytearray()
				while hungryHungry < idlen:
					#print '\n', hp1, hp2, len(myDat), len(myDat[hp1]), hungryHungry,'/',idlen
					if len(myDat[hp1][hp2:]) <= idlen-hungryHungry:
						hungryHungry += len(myDat[hp1][hp2:])
						seq.extend(myDat[hp1][hp2:])
						myDat[hp1] = myDat[hp1][:hp2]
						hp1 += 1
						hp2 = 0
					else:
						seq.extend(myDat[hp1][hp2:hp2+idlen-hungryHungry])
						myDat[hp1] = myDat[hp1][:hp2]+myDat[hp1][hp2+idlen-hungryHungry:]
						break

				idT = 'D'
				lenAdj -= idlen
				prevP2_adj -= idlen
				indBp = ind+1

			v = [idT,idlen,ind,indBp,str(seq)]
			indelList.append(v)
			if tuple(n) in indelsToAttempt_AF:
				indelList_AF[tuple(v)] = indelsToAttempt_AF[tuple(n)]
			afterThis.extend([ind,indBp])
			if idT == 'D':
				idlen = -idlen
			addThis.extend([addThis[-1],addThis[-1]-idlen])
			if svInd != -1:
				svInds[len(indelList)-1] = svInd

			prevP1 = p1
			prevP2 = p2+idlen
		del cumLens
		del indelsToAttempt
		del indelsToAttempt_AF


		"""//////////////////////////////////////////////////////
		////////////    CONCAT INTO SINGLE STRING    ////////////
		//////////////////////////////////////////////////////"""

		myDat = bytearray(''.join(myDat))
		myLen = len(myDat)


		"""/////////////////////////////////////////////////
		////////////       INTRODUCE SNPS       ////////////
		/////////////////////////////////////////////////"""

		uppers = [ord(n) for n in NUM_NUCL]
		lowers = [ord(n) for n in NUM_NUCL_LOWER]
		both   = uppers+lowers
		nSNPs  = 0
		snps = {}
		if NATURAL_SNPS:
			rSF = ((random.random()-.5)*(2*RNG_MULT) + 1.)*SNP_FREQ*var_mult
			nSNPs = int(rSF*myLen+.5)
			for i in xrange(nSNPs):
				spos = random.randint(0,myLen-1)
				# have we already inserted a snp here? if so, skip plz
				if spos in snps:
					continue
				prevBase = myDat[spos]
				if prevBase in both:
					newBase = randEvent(cpSNP[NUCL_NUM[prevBase]])
					if prevBase in uppers:
						newBase = NUM_NUCL[newBase]
					elif prevBase in lowers:
						newBase = NUM_NUCL_LOWER[newBase]
					else:
						continue
					snps[spos] = (prevBase,newBase,spos)
					myDat[spos] = newBase
					#print chr(prevBase),'->',newBase,chr(myDat[spos]), spos

		inSNP_AF = {}	# sloppy temp dictionary
		for i in xrange(len(input_snps)):
			n = input_snps[i]
			spos = n[0]-1
			sind = bisect.bisect(afterThis,spos)
			if sind == len(afterThis):
				sind -= 1
			spos -= addThis[sind]

			#print n, spos, chr(myDat[spos])
			if spos not in snps and chr(myDat[spos]).upper() == n[1]:
				snps[spos] = (myDat[spos],n[2],spos)
				myDat[spos] = n[2]
				nSNPs += 1
				inSNP_AF[spos] = input_snps_AF[i]
			else:
				# the snp we're trying to insert was affected by a random variant, skip it, sorry.
				pass
		del input_snps_AF

		nSNPs   = len(snps)			# I'm reusing this variable name? -.-
		nIndels = len(indelList)
		# initialize variant coverage dictionaries
		snpKeys = sorted(snps.keys())
		snpCoverage = {}
		snpReads    = {}
		snpTargeted = {}
		snpAB       = {}
		for k in snpKeys:
			snpCoverage[k] = 0
			snpReads[k]    = []
			if k in inSNP_AF:
				if inSNP_AF[k] < 1.:
					snpAB[k] = inSNP_AF[k]
			else:
				if random.random() < HET_VAR:
					snpAB[k] = AB[randEvent(cpAB)]
		indelCoverage = [0 for n in indelList]
		indelReads    = [[] for n in indelList]
		indelTargeted = {}
		indelAB       = {}
		#print indelList
		#print indelList_AF
		for i in xrange(len(indelList)):
			thij = tuple(indelList[i])
			if thij in indelList_AF:
				if indelList_AF[thij] < 1.:
					indelAB[thij] = indelList_AF[thij]
			else:
				if random.random() < HET_VAR:
					indelAB[thij] = AB[randEvent(cpAB)]
		print '\n'


		"""/////////////////////////////////////////////////////
		////////////    MORE MISCELLANEOUS STUFF    ////////////
		/////////////////////////////////////////////////////"""

		# print number of variants introduced and how long it took to do so
		if NATURAL_INDELS or NATURAL_SNPS or NATURAL_SVS:
			print '{0:.3f} (sec)'.format(time.time()-indelStart)
			print nSNPs, 'SNPs', '('+str(len(snpAB))+' hets)'
			print nIndels, 'indels', '('+str(len(indelAB))+' hets)'
			print nSVs, 'SVs'
			print ''

		# save new reference to file if desired
		if SAVE_CORRECT_REF and JOB_ID == 1:
			colWidth = 50
			of = open(OUTFILE_NAME+'_'+refName+'_correctRef.fa','wb')
			of.write('>'+refName+'\n')
			for i in xrange(0,myLen,colWidth):
				of.write(myDat[i:i+colWidth]+'\n')
			of.close()

		# convert all nucleotides to upper-case
		sys.stdout.write('capitalizing all nucleotides... ')
		sys.stdout.flush()
		tt = time.time()
		myDat = myDat.upper()
		print '{0:.3f} (sec)\n'.format(time.time()-tt)


		"""/////////////////////////////////////////////////////
		////////////   SOME WINDOW PREPROCESSING   /////////////
		/////////////////////////////////////////////////////"""


		# construct regions to sample from (from input bed file, if present)
		nBedTargetedBP = 0
		targRegionsFl  = []
		targetRegionsToSample = []
		if INPUT_BED != None:
			bedfile = open(INPUT_BED,'r')
			for line in bedfile:
				splt = line.split('\t')
				if splt[0] == refName:
					regionLen   = int(splt[2])-int(splt[1])
					if regionLen < MIN_PROBE_LEN:
						continue
					nBedTargetedBP += regionLen
					targRegionsFl.extend((int(splt[1]),int(splt[2])))
					origCoords  = ( max([int(splt[1])-MIN_SAMPLE_SIZE, 0]), min([int(splt[2])+MIN_SAMPLE_SIZE, originalLen-1]) )
					myDatCoords = ( origCoords[0]-addThis[bisect.bisect(afterThis,origCoords[0])-1], origCoords[1]-addThis[bisect.bisect(afterThis,origCoords[1])-1] )
					#print origCoords,'-->',myDatCoords
					addMe = 1
					#
					#
					# WARNING: AT THE MOMENT THIS EXPECTS BED INPUT IS SORTED.
					#
					#
					# detect redundant regions
					for i in xrange(len(targetRegionsToSample)):
						if myDatCoords[0] < targetRegionsToSample[i][1]:
							if myDatCoords[1] <= targetRegionsToSample[i][1]:
								addMe = 0
								break
							else:
								targetRegionsToSample[i] = (targetRegionsToSample[i][0],myDatCoords[1])
								addMe = 0
								break
					if addMe:
						targetRegionsToSample.append(myDatCoords)
			bedfile.close()


		# determine which variants were targeted
		if INPUT_BED != None:
			for i in xrange(len(indelList)):
				origInd = indelList[i][2]+addThis[i*2]
				if bisect.bisect(targRegionsFl,origInd)%2 == 1:
					indelTargeted[i] = 1
			for k in snps.keys():
				spos    = snps[k][2]
				si      = bisect.bisect(afterThis,spos)-1
				origInd = spos+addThis[si]+1
				if bisect.bisect(targRegionsFl,origInd)%2 == 1:
					snpTargeted[k] = 1
		del targRegionsFl


		# compute GC % of each window
		if PAIRED_END:
			windowShift = COV_WINDOW-FRAGMENT_SIZE
		else:
			windowShift = COV_WINDOW-READLEN
		gcp = []
		slidingWindowAfter = 0	# at what index in targetRegionsToSample begin uniform sliding windows?
		if INPUT_BED == None:
			for i in xrange(0,myLen,windowShift):
				(bi,bf) = (i,i+COV_WINDOW)
				targetRegionsToSample.append((bi,bf))
				gcp.append(float(myDat[bi:bf].count('C')+myDat[bi:bf].count('G'))/(bf-bi))

			# account for weird effects around the final window
			del targetRegionsToSample[-1]
			del gcp[-1]
			targetRegionsToSample[-1] = (targetRegionsToSample[-1][0],len(myDat))

			targetCov = [RELATIVE_GC_COVERAGE_BIAS[int(100.*n)] for n in gcp]
			alpha = (len(gcp)*AVG_COVERAGE)/sum(targetCov)
			alpha *= float(windowShift)/COV_WINDOW # account for window overlap
		  	targetCov = [alpha*n for n in targetCov]
		else:
			for i in xrange(len(targetRegionsToSample)):
				(bi,bf) = targetRegionsToSample[i]
				gcp.append(float(myDat[bi:bf].count('C')+myDat[bi:bf].count('G'))/(bf-bi))
			targetCov = [RELATIVE_GC_COVERAGE_BIAS[int(100.*n)] for n in gcp]
			if len(targetCov) > 0 and sum(targetCov) > 0:
				alpha = (len(gcp)*(AVG_COVERAGE-BED_COVERAGE))/sum(targetCov)
				targetCov = [alpha*n for n in targetCov]

			# so we got coverage in desired regions, let's add a little bit of the rest if desired
			if BED_COVERAGE > 0.:
				slidingWindowAfter = len(gcp)
				gcp = []
				for i in xrange(0,myLen,windowShift):
					(bi,bf) = (i,i+COV_WINDOW)
					targetRegionsToSample.append((bi,bf))
					gcp.append(float(myDat[bi:bf].count('C')+myDat[bi:bf].count('G'))/(bf-bi))

				del targetRegionsToSample[-1]
				del gcp[-1]
				targetRegionsToSample[-1] = (targetRegionsToSample[-1][0],len(myDat))

				targetCovR = [RELATIVE_GC_COVERAGE_BIAS[int(100.*n)] for n in gcp]
				alpha = (len(gcp)*BED_COVERAGE)/sum(targetCovR)
				alpha *= float(windowShift)/COV_WINDOW # account for window overlap
				targetCovR = [alpha*n for n in targetCovR]
				targetCov.extend(targetCovR)

		# determine which windows contain Ns and enumerate their non-N regions

		hasN = {}
		for i in xrange(len(targetRegionsToSample)):

			(bi,bf) = targetRegionsToSample[i]

			currentDat = myDat[bi:bf]
			if 'N' in currentDat:
				if currentDat.count('N') != len(currentDat):
					inNind   = 0
					NonNregions = []
					for ci in xrange(1,len(currentDat)):
						if currentDat[ci-1] != ASCII_N and currentDat[ci] == ASCII_N:
							NonNregions.append([inNind,ci-1])
						elif currentDat[ci-1] == ASCII_N and currentDat[ci] != ASCII_N:
							inN = 1
							inNind = ci
					if currentDat[-1] != ASCII_N:
						NonNregions.append([inNind,len(currentDat)-1])
					#print NonNregions
					hasN[i] = NonNregions


		"""//////////////////////////////////////////////////
		////////////       READ SIMULATION       ////////////
		//////////////////////////////////////////////////"""

		if PAIRED_END:
			print 'simulating paired-end reads...'
		else:
			print 'simulating single-end reads...'
		
		# pre-compute number of reads for each window for each job
		jobIterators = [xrange(n,len(targetRegionsToSample),JOB_TOT) for n in xrange(JOB_TOT)]
		jobOffsets   = [0]
		myJobReadsToAdd = []
		myJobRegions    = []
		for i in xrange(len(jobIterators)):
			n = jobIterators[i]
			jobReadCount = 0
			for j in n:

				(bi,bf) = targetRegionsToSample[j]

				currentLen = bf-bi
				if myDat[bi:bf].count('N') == currentLen:
					continue

				# determine regions around this window that we're going to sample reads from
				prevN = j-1 in hasN
				nextN = j+1 in hasN
				regionsToSample = []
				coverageMult    = []
				if j >= slidingWindowAfter:
					if j not in hasN:
						if prevN:
							# merge end of previous window
							if hasN[j-1][-1][1] == COV_WINDOW-1:
								bi += hasN[j-1][-1][0]-windowShift
						if nextN:
							# merge beginning of next window
							if hasN[j+1][0][0] == 0:
								bf += hasN[j+1][0][1]+windowShift-currentLen+1
						regionsToSample.append((bi,bf))
						coverageMult.append(1.)
					else:
						# sample isolated regions
						for region in hasN[j]:
							if region[0] != 0 and region[1] != currentLen-1 and region[1]-region[0] >= MIN_SAMPLE_SIZE:
								regionsToSample.append((bi+region[0],bi+region[1]+1))
								coverageMult.append(1.)
						if prevN:
							# merge end of previous window
							poti = bi-windowShift+hasN[j-1][-1][0]
							potf = bi+hasN[j][0][1]+1
							if hasN[j-1][-1][1] == COV_WINDOW-1 and hasN[j][0][0] == 0 and potf-poti >= MIN_SAMPLE_SIZE:
								regionsToSample.append((poti,potf))
								coverageMult.append(0.5)
						if nextN:
							# merge beginning of next window
							poti = bi+hasN[j][-1][0]
							potf = bi+windowShift+hasN[j+1][0][1]+1
							if hasN[j][-1][1] == COV_WINDOW-1 and hasN[j+1][0][0] == 0 and potf-poti >= MIN_SAMPLE_SIZE:
								regionsToSample.append((poti,potf))
								coverageMult.append(0.5)
				else:
					if j not in hasN:
						if currentLen >= MIN_SAMPLE_SIZE:
							regionsToSample.append((bi,bf))
							coverageMult.append(1.)
					else:
						for region in hasN[j]:
							if region[1]-region[0]-1 >= MIN_SAMPLE_SIZE:
								regionsToSample.append((bi+region[0],bi+region[1]+1))
								coverageMult.append(1.)

				for k in xrange(len(regionsToSample)):
					(bi,bf) = regionsToSample[k]
					cm      = coverageMult[k]

					currentLen = bf-bi
					currentCov = ((random.random()-.5)*(2*COV_MULT) + 1.)*targetCov[j]*cm
					if PAIRED_END:
						readsToAdd = int((currentLen*currentCov)/(2*READLEN)+0.5)
						jobReadCount += readsToAdd*2
					else:
						readsToAdd = int((currentLen*currentCov)/(READLEN)+0.5)
						jobReadCount += readsToAdd
					#print (currentLen,readsToAdd)

					if i == JOB_ID-1:
						myJobReadsToAdd.append(readsToAdd)
						myJobRegions.append(regionsToSample[k])

			jobOffsets.append(jobReadCount+jobOffsets[-1])
		myJobOffset   = jobOffsets[JOB_ID-1]

		tt = time.time()
		if not SEQUENCING_QSCORES:
			if AVG_SER == 0.:
				q1 = bytearray([Qscores[-1]+qOffset]*READLEN)
			else:
				q1 = bytearray([DUMMY_QSCORE+qOffset]*READLEN)
			q2 = q1
		nReads = 0
		PRINT_EVERY_PERCENT = 5
		totalProgress       = PRINT_EVERY_PERCENT
		for i in xrange(len(myJobRegions)):
			#print (i+1),'/',len(myJobRegions)
			while 100.*float(i+1)/len(myJobRegions) >= totalProgress:
				sys.stdout.write(str(totalProgress)+'% ')
				sys.stdout.flush()
				totalProgress += PRINT_EVERY_PERCENT

			(bi,bf)    = myJobRegions[i]
			readsToAdd = myJobReadsToAdd[i]
			currentDat = myDat[bi:bf]
			currentDatRC = [TO_UPPER_COMP[n] for n in currentDat[::-1]]

			if len(indelList) == 0:
				relevantAT = []
				wHitInds    = []
				afwi = 0
			else:
				afwi = bisect.bisect(afterThis,bi-MAX_INDEL)
				afwf = bisect.bisect(afterThis,bf+MAX_INDEL)
				relevantAT = afterThis[afwi:afwf]
				wHitInds    = indelList[afwi/2:min([afwf/2,len(indelList)])]
			if len(snpKeys) == 0:
				relevantSNP = []
				snpi = 0
			else:
				snpi = bisect.bisect(snpKeys,bi-MAX_INDEL)
				snpf = bisect.bisect(snpKeys,bf+MAX_INDEL)
				relevantSNP = snpKeys[snpi:snpf]

			for rta in xrange(readsToAdd):

				hitInds = [n for n in wHitInds]
				if PAIRED_END:
					fl = FRAGLEN[randEvent(cpFrag)]
					n0 = random.randint(0,len(currentDat)-fl-1)
					revoff = len(currentDat)-n0-fl
					r1posMyDat = bi+n0+1
					r2posMyDat = bi+n0+1+fl-READLEN
				else:
					n0 = random.randint(0,len(currentDat)-READLEN-1)
					r1posMyDat = bi+n0+1


				#
				# readData
				#
				r1 = bytearray(currentDat[n0:n0+READLEN])
				if PAIRED_END:
					r2 = bytearray(currentDatRC[revoff:revoff+READLEN])

				# re-init seqErrObj
				if PAIRED_END:
					seqErrObj.newReads(r1posMyDat,r2posMyDat)
				else:
					seqErrObj.newReads(r1posMyDat)

				#
				# quality scores
				#
				if SEQUENCING_QSCORES:
					q1 = bytearray()
					q1.append(randEvent(initQProb)+qOffset)
					if PAIRED_END:
						q2 = bytearray()
						q2.append(randEvent(initQProb)+qOffset)
					for j in xrange(1,READLEN):
						qscore1 = randEvent(allQscoreCumProb[j*nQscores+q1[j-1]])+qOffset
						q1.append(qscore1)
						if PAIRED_END:
							qscore2 = randEvent(allQscoreCumProb[j*nQscores+q2[j-1]])+qOffset
							q2.append(qscore2)

					#
					# introduce sequencing errors based on quality scores
					#
					if SEQUENCING_SNPS:
						for j in xrange(READLEN):
							if random.random() < pError[q1[j]]:
								seqErrObj.introduceError(1, r1, j)
							if PAIRED_END and random.random() < pError[q2[j]]:
								seqErrObj.introduceError(2, r2, j)
				else:
					#
					# if no q-scores, introduce sequencing errors based on desired average frequency (with moderate positional preference towards end of read)
					#
					if SEQUENCING_SNPS:
						for j in xrange(READLEN):
							if random.random() < positionSSErate[j]:
								seqErrObj.introduceError(1, r1, j)
							if PAIRED_END and random.random() < positionSSErate[j]:
								seqErrObj.introduceError(2, r2, j)


				# append sequencing indel errors to list of variation we need to consider for these reads
				hitInds.extend(seqErrObj.myHitInds)

				# oragnize some variable in preparation for constructing cigar strings and shifting read positions to account for variants
				cigarsAndAdjustedPos = []
				if PAIRED_END:
					rtt = [r1posMyDat,r2posMyDat]
				else:
					rtt = [r1posMyDat]

				#if SAVE_VCF:
				(r1i, r1f) = (r1posMyDat, r1posMyDat+READLEN)
				if PAIRED_END:
					(r2i, r2f) = (r2posMyDat, r2posMyDat+READLEN)
				
				#
				# if indel is heterozygous, lets trade it out for the ref allele and remove it from cigar-string computations (with corresponding allele balance probability of course)
				#
				skipTheseInds = []
				adjustedVariants = {tuple(n):0 for n in hitInds}
				for j in xrange(len(hitInds)):
					thij     = tuple(hitInds[j])
					adj      = adjustedVariants[thij]
					indStart = thij[2]
					if (thij in indelAB and random.random() < indelAB[thij]) or (hitInds[j] in seqErrObj.sie_errors):
						#print ''
						#print hitInds[j]
						#print (r1i,r1f),(r2i,r2f),indStart
						if (indStart >= r1i and indStart < r1f-1):
							#print 'read #:',str(nReads+myJobOffset+bigReadNameOffset)+'/1'
							#print r1
							r1Ind = indStart-r1i+1

							# if alt allele is deletion, insert ref bases into read and trim to read length
							if thij[0] == 'D':
								r1 = r1[:r1Ind+adj] + thij[4].upper() + r1[r1Ind+adj:]
								r1 = r1[:READLEN]
								# every variant in the read after this one needs to have its position adjusted, oh joy.
								for j2 in xrange(j+1,len(hitInds)):
									adjustedVariants[tuple(hitInds[j2])] += thij[1]

							# if alt allele is insertion, delete alt bases and fill in with subsequent sequence data.
							# have to make sure we catch any variants that might occur in the "reserve" sequence that we append!
							elif thij[0] == 'I':
								buflen    = min([READLEN,thij[1]])
								extra_adj = max([0,thij[3]+adj-r1f])
								# check if enough room left in window to fill in gap left by removing insertion
								if len(currentDat[n0+READLEN+extra_adj-adj:]) < buflen:
									continue
								else:
									(e1, e2) = (n0+READLEN+extra_adj-adj, n0+READLEN+extra_adj-adj+buflen)
									extra = bytearray(currentDat[e1:e2])
									r1 = r1[:r1Ind+adj] + r1[r1Ind+thij[1]+adj:] + extra
									r1 = r1[:READLEN]
									for j2 in xrange(j+1,len(hitInds)):
										adjustedVariants[tuple(hitInds[j2])] -= thij[1]

						elif PAIRED_END and (indStart >= r2i and indStart < r2f-1):
							#print 'read #:',str(nReads+myJobOffset+bigReadNameOffset)+'/2'
							#print r2
							#print bytearray([TO_UPPER_COMP[n] for n in r2][::-1])
							r2IndComp = indStart-r2i+1

							if thij[0] == 'D':
								r2 = bytearray([TO_UPPER_COMP[n] for n in r2][::-1])
								r2 = r2[:r2IndComp+adj] + thij[4].upper() + r2[r2IndComp+adj:]
								r2 = r2[:READLEN]
								r2 = bytearray([TO_UPPER_COMP[n] for n in r2][::-1])
								for j2 in xrange(j+1,len(hitInds)):
									adjustedVariants[tuple(hitInds[j2])] += thij[1]

							elif thij[0] == 'I':
								buflen    = min([READLEN,thij[1]])
								extra_adj = max([0,thij[3]+adj-r2f])
								if len(currentDatRC[:revoff-extra_adj+adj]) < buflen:
									print '\n\nproblem, captain!\n\n'
								else:
									r2 = bytearray([TO_UPPER_COMP[n] for n in r2][::-1])
									(e1, e2) = (revoff-extra_adj+adj-buflen, revoff-extra_adj+adj)
									extra = bytearray(currentDatRC[e1:e2])
									r2 = r2[:r2IndComp+adj] + r2[r2IndComp+thij[1]+adj:]
									r2 = extra + bytearray([TO_UPPER_COMP[n] for n in r2][::-1])
									r2 = r2[len(r2)-READLEN:]
									for j2 in xrange(j+1,len(hitInds)):
										adjustedVariants[tuple(hitInds[j2])] -= thij[1]

						else:
							continue

						if hitInds[j] not in seqErrObj.sie_errors:
							skipTheseInds.append(hitInds[j])

				for rp in rtt:
					# I want to get nucleotides 'ni' to 'nf' from our modified reference and get the
					# corresponding bases from the original reference, as well as a cigar string that relates them
					ni = rp
					nf = rp+READLEN		# ni + READLEN

					if len(relevantAT) == 0 and len(seqErrObj.sie_errors) == 0:
						if len(indelList) > 0:
							if afwi == len(addThis):
								ni += addThis[-1]
							else:
								ni += addThis[afwi]
						cigarsAndAdjustedPos.extend([ni,str(READLEN)+'M'])
						continue

					ptbi = bisect.bisect(relevantAT,ni)-1
					#if ni in relevantAT:
					#	ptbi -= 1
					ptbf = bisect.bisect(relevantAT,nf)-1
					if nf in relevantAT:
						ptbf -= 1
					ptbi += afwi
					ptbf += afwi
					nir = ni+addThis[ptbi]
					nfr = nf+addThis[ptbf]

					bpi = ptbi/2
					bpf = min([ptbf/2,len(indelList)-1])
					affectedBp = copy.deepcopy(indelList[bpi:bpf+10])	# grab 10 extra variants in case het-insert deleted region fill-in contains variants
					## make sure we catch insertions near end of read
					#if indelList[bpf][2]+1 < nf:
					#	affectedBp.append(indelList[bpf])

					# insert sequencing indel errors to the list of indels that affect the bp in our read
					for n in seqErrObj.sie_errors:
						if n[0] == 'I':
							n = ['D',n[1],n[2],n[2]+1,n[4]]
						elif n[0] == 'D':
							n = ['I',n[1],n[2],n[2]+n[1]+1,n[4]]
						affectedBp.append(n)
					affectedBp.sort(key=lambda x: int(x[2]))

					# take out indels if this read is heterozygous and has ref allele instead
					affectedBp = [n for n in affectedBp if n not in skipTheseInds]
					
					#if len(seqErrObj.sie_errors):
					#	print 'read #:',str(nReads+myJobOffset+bigReadNameOffset)+'/1'
					#	print 'affectedBp:',affectedBp

					# adjust indel offsets of variants following hets that have ref instead
					for i in xrange(len(affectedBp)):
						tabi = tuple(affectedBp[i])
						if tabi in adjustedVariants:
							adj = adjustedVariants[tabi]
							affectedBp[i][2] += adj
							affectedBp[i][3] += adj
					# remove irrelevant variants
					affectedBp = [n for n in affectedBp if n[2] < nf]
					# get rid of deletions that occur at last index of read
					if len(affectedBp) and affectedBp[-1][0] == 'D' and affectedBp[-1][2] == nf-1:
						del affectedBp[-1]


					myCigar = ''
					soFar = ni
					for n in affectedBp:
						chunkSize = n[2]-soFar+1
						if chunkSize > 0:
							soFar += chunkSize
							myCigar = myCigar+str(chunkSize)+'M'
							# spanning a deletion
							if n[0] == 'D':
								myCigar = myCigar+str(n[1])+'D'
							# at least partially spanning the beginning of an insertion
							elif n[0] == 'I':
								chunkSize = min([nf-soFar,n[3]-soFar])
								if chunkSize > 0:
									myCigar = myCigar+str(chunkSize)+'I'
									soFar += chunkSize
						else:
							# at least partially spanning the end of an insertion
							if n[0] == 'I':
								chunkSize = n[3]-soFar
								# is this read fully inside of a large insertion?
								if chunkSize >= READLEN:
									myCigar = myCigar+str(rp-n[2])+'P'+str(READLEN)+'I'
									soFar += READLEN
									nir += chunkSize-n[1]
									break

								if chunkSize > 0:
									myCigar = myCigar+str(chunkSize)+'I'
									soFar += chunkSize
									nir += chunkSize-n[1]

					if soFar < nf:
						myCigar = myCigar+str(nf-soFar)+'M'
					cigarsAndAdjustedPos.extend([nir,myCigar])
				#print cigarsAndAdjustedPos

				if PAIRED_END:
					[r1pos, cigar1, r2pos, cigar2] = cigarsAndAdjustedPos
					r1NameSuffix = str(r1pos)+'-'+str(r2pos)+'-'+str(nReads+myJobOffset+bigReadNameOffset)+'/1'
					r2NameSuffix = str(r1pos)+'-'+str(r2pos)+'-'+str(nReads+myJobOffset+bigReadNameOffset)+'/2'
					r1Name = refName+'-'+r1NameSuffix
					r2Name = refName+'-'+r2NameSuffix
					OUTFQ1.write('@'+r1Name+'\n'+r1+'\n+\n'+q1+'\n')
					OUTFQ2.write('@'+r2Name+'\n'+r2+'\n+\n'+q2+'\n')
				else:
					[r1pos, cigar1] = cigarsAndAdjustedPos
					r1NameSuffix = str(r1pos)+'-'+str(nReads+myJobOffset+bigReadNameOffset)+'/1'
					r1Name = refName+'-'+r1NameSuffix
					OUTFQ1.write('@'+r1Name+'\n'+r1+'\n+\n'+q1+'\n')


				#
				#	Notate SNP coverage
				#
				for sk in relevantSNP:
					if (sk >= r1i and sk < r1f-1):

						# if this snp is heterozygous, lets trade it out for the ref allele with it's corresponding allele balance probability
						if sk in snpAB and random.random() < snpAB[sk]:
							r1[sk-r1i+1] = chr(snps[sk][0]).upper()
						# if not traded for ref allele lets count that we covered it.
						else:
							snpReads[sk].append(r1NameSuffix)
							snpCoverage[sk] += 1
							
					elif PAIRED_END and (sk >= r2i and sk < r2f-1):

						if sk in snpAB and random.random() < snpAB[sk]:
							r2[(READLEN-1)-(sk-r2i+1)] = TO_UPPER_COMP[chr(snps[sk][0])]
						else:
							snpReads[sk].append(r2NameSuffix)
							snpCoverage[sk] += 1

				#
				#	Notate Indel coverage
				#
				for j in xrange(len(hitInds)-1,-1,-1):	# remove SIE errors
					if hitInds[j] in seqErrObj.sie_errors:
						del hitInds[j]
				for j in xrange(len(hitInds)):
					if hitInds[j] not in skipTheseInds:
						indStart = hitInds[j][2]
						if (indStart >= r1i and indStart < r1f-1):
							indelCoverage[j+afwi/2] += 1
							indelReads[j+afwi/2].append(r1NameSuffix)
						elif PAIRED_END and (indStart >= r2i and indStart < r2f-1):
							indelCoverage[j+afwi/2] += 1
							indelReads[j+afwi/2].append(r2NameSuffix)

				#
				#	Write SAM file entry
				#
				if SAVE_SAM:

					if PAIRED_END:

						samFlag1 = str(int('001100011',2))
						pos1     = str(r1pos)	# position of first matching base
						mapQ1    = str(70)		# set to max value
						mateRef  = '='			# read and mate have same reference
						matePos1 = str(r2pos)	# position of mates right-most matching base
						tLen1    = str(fl)
						seq1     = str(r1)
						qstr1    = str(q1)
						aln1     = ''			# no special alignment info..

						outData  = [r1Name[:-2],samFlag1,refName,pos1,mapQ1,cigar1,mateRef,matePos1,tLen1,seq1,qstr1,aln1]
						outData  = [n for n in outData if n != '']
						OUTSAM.write('\t'.join(outData))
						OUTSAM.write('\n')

						samFlag2 = str(int('010010011',2))
						pos2     = str(r2pos)
						mapQ2    = str(70)
						mateRef  = '='
						matePos2 = str(r1pos)
						tLen2    = str(-fl)
						seq2     = str(bytearray([TO_UPPER_COMP[n] for n in r2[::-1]]))
						qstr2    = str(q2)
						aln2     = ''

						outData  = [r1Name[:-2],samFlag2,refName,pos2,mapQ2,cigar2,mateRef,matePos2,tLen2,seq2,qstr2,aln2]
						outData  = [n for n in outData if n != '']
						OUTSAM.write('\t'.join(outData))
						OUTSAM.write('\n')

					else:

						samFlag1 = str(int('000000000',2))
						pos1     = str(r1pos)
						mapQ1    = str(70)
						mateRef  = '*'
						matePos1 = str(0)
						tLen1    = str(0)
						seq1     = str(r1)
						qstr1    = str(q1)
						aln1     = ''

						outData  = [r1Name[:-2],samFlag1,refName,pos1,mapQ1,cigar1,mateRef,matePos1,tLen1,seq1,qstr1,aln1]
						outData  = [n for n in outData if n != '']
						OUTSAM.write('\t'.join(outData))
						OUTSAM.write('\n')

				if PAIRED_END:
					nReads += 2
				else:
					nReads += 1

		print '\n',nReads,'(reads)','{0:.3f} (sec),'.format(time.time()-tt),int((nReads*READLEN)/(time.time()-tt)),'(bp/sec)'
		if nReads > 0:
			print 'SSE rate:',float(seqErrObj.numSSE_total)/(nReads*READLEN)
			print 'SIE rate:',float(seqErrObj.numSIE_total)/(nReads*READLEN)


		"""//////////////////////////////////////////////////
		////////////      WRITE VCF OUTFILE      ////////////
		//////////////////////////////////////////////////"""


		if SAVE_VCF:

			allVariants = []
			snpsAccountedFor = {}
			SVsAccountedFor  = {}
			for k in snps.keys():
				snpsAccountedFor[k] = 0

			SVids = {'BD':0, 'BI':0, 'R':0, 'T':0, 'IV':0, 'IT':0}
			for i in xrange(len(indelList)):
				n = indelList[i]
				origInd = n[2]+addThis[i*2]
				myInd   = n[2]-1
				origBase = chr(myDat[myInd])
				if myInd in snps:
					origBase = chr(snps[myInd][0]).upper()
					snpsAccountedFor[myInd]
				if n[0] == 'I':
					for j in xrange(myInd+1,myInd+n[1]+1):
						if j in snps:
							snpsAccountedFor[j] = 1

				if i in svInds:
					svi = svInds[i]
					Svv = SVsToAttempt[svi][0]
					idT = SVsToAttempt[svi][1]
					svidNum = SVids[idT]
					outDat = ()

					if svi not in SVsAccountedFor:
						SVsAccountedFor[svi] = 1

						VarTargeted = ''
						if i in indelTargeted:
							VarTargeted = ';TARGETED=1'
						VarAF = ''
						if tuple(indelList[i]) in indelAB:
							VarAF = ';AF={0:.3f}'.format(indelAB[tuple(indelList[i])])
						VarReads = ''
						if len(indelReads[i]) > 0:
							VarReads = ';READS='+','.join(indelReads[i])

						if idT == 'BI':
							#print 'BigInsert', Svv
							newStr = str(myDat[myInd:myInd+len(Svv[1])+1])
							outDat = (refName,str(Svv[0]),'.',origBase,newStr.upper(),'.','PASS','SVTYPE=INS;END='+str(Svv[0])+';SVLEN='+str(len(Svv[1])))
							if (int(outDat[1]),outDat[3],outDat[4]) in input_inds:
								outDat = (refName,str(Svv[0]),'.',origBase,newStr.upper(),'.','PASS','DP='+str(indelCoverage[i])+VarTargeted+VarAF+VarReads)
						elif idT == 'BD':
							#print 'BigDeletion', Svv
							outDat = (refName,str(Svv[0]),'.',str(origBase+n[4]).upper(),chr(myDat[myInd]),'.','PASS','SVTYPE=DEL;END='+str(Svv[0]+Svv[1])+';SVLEN='+str(-Svv[1]))
							if (int(outDat[1]),outDat[3],outDat[4]) in input_inds:
								outDat = (refName,str(Svv[0]),'.',str(origBase+n[4]).upper(),chr(myDat[myInd]),'.','PASS','DP='+str(indelCoverage[i])+VarTargeted+VarAF+VarReads)
						elif idT == 'R':
							#print 'Repeat', Svv
							newStr = str(myDat[myInd:myInd+len(Svv[1])*Svv[2]+1])
							outDat = (refName,str(Svv[0]),'.',origBase,newStr.upper(),'.','PASS','SVTYPE=DUP;END='+str(Svv[0]+len(Svv[1])*Svv[2])+';SVLEN='+str(len(Svv[1])*Svv[2]))
					elif SVsAccountedFor[svi] == 1:
						SVsAccountedFor[svi] = 2
						if idT == 'IV':
							#print 'Inversion', Svv
							event = 'INV'+str(svidNum)
							outDat = (refName,str(Svv[0]),event,origBase,'<INV>','.','PASS','SVTYPE=INV;END='+str(Svv[0]+Svv[1])+';SVLEN='+str(Svv[1]))
						elif idT == 'IT':
							#print 'InvertedTranslocation', Svv
							event = 'INV-TRANS'+str(svidNum)
							origBase = chr(myDat[Svv[2]])
							outDat = (refName,str(Svv[2]),event,origBase,'<INV-TRANS>','.','PASS','SVTYPE=INV-TRANS;END='+str(Svv[2])+';SVLEN='+str(Svv[1])+';TRANSREF='+refName+';TRANSPOS='+str(Svv[0]))
						elif idT == 'T':
							#print 'Translocation', Svv
							event = 'TRANS'+str(svidNum)
							origBase = chr(myDat[Svv[2]])
							outDat = (refName,str(Svv[2]),event,origBase,'<TRANS>','.','PASS','SVTYPE=TRANS;END='+str(Svv[2])+';SVLEN='+str(Svv[1])+';TRANSREF='+refName+';TRANSPOS='+str(Svv[0]))
						SVids[idT] += 1
				else:
					VarTargeted = ''
					if i in indelTargeted:
						VarTargeted = ';TARGETED=1'
					VarAF = ''
					if tuple(indelList[i]) in indelAB:
						VarAF = ';AF={0:.3f}'.format(indelAB[tuple(indelList[i])])
					VarReads = ''
					if len(indelReads[i]) > 0:
						VarReads = ';READS='+','.join(indelReads[i])

					if n[0] == 'D':
						outDat = (refName,str(origInd),'.',str(origBase+n[4]).upper(),chr(myDat[myInd]),'.','PASS','DP='+str(indelCoverage[i])+VarTargeted+VarAF+VarReads)

					elif n[0] == 'I':
						newStr = str(myDat[myInd:myInd+n[1]+1])
						outDat = (refName,str(origInd),'.',origBase,newStr.upper(),'.','PASS','DP='+str(indelCoverage[i])+VarTargeted+VarAF+VarReads)

				if outDat != ():
					allVariants.append(outDat)

			for k in snpKeys:
				if snpsAccountedFor[k] == 0:
					spos    = snps[k][2]
					si      = bisect.bisect(afterThis,spos)-1
					origInd = spos+addThis[si]+1

					VarTargeted = ''
					if k in snpTargeted:
						VarTargeted = ';TARGETED=1'
					VarAF = ''
					if k in snpAB:
						VarAF = ';AF={0:.3f}'.format(snpAB[k])
					VarReads = ''
					if len(snpReads[k]) > 0:
						VarReads = ';READS='+','.join(snpReads[k])

					outDat = (refName,str(origInd),'.',chr(snps[k][0]).upper(),chr(myDat[spos]),'.','PASS','DP='+str(snpCoverage[k])+VarTargeted+VarAF+VarReads)
					#print outDat
					allVariants.append(outDat)

			allVariants = [x for (y,x) in sorted(zip([int(n[1]) for n in allVariants],allVariants))]
			for n in allVariants:
				for m in n:
					if m == n[-1]:
						OUTVCF.write(m+'\n')
					else:
						OUTVCF.write(m+'\t')

		# update read name offset for the next reference
		if nReads > 0:
			bigReadNameOffset    += jobOffsets[-1]
			totalSNPs            += nSNPs
			totalInds            += nIndels
			totalSTVs            += nSVs
			sequencesSampledFrom += 1
			totalBedTargetedBP   += nBedTargetedBP
			readsFromThisJob     += nReads
			if INPUT_BED == None:
				totalBPtargeted     += len(myDat)-myDat.count('N')
			else:
				for i in xrange(len(targetRegionsToSample)):
					if i < slidingWindowAfter:
						region = targetRegionsToSample[i]
						totalBPtargeted += region[1]-region[0]

	finalTime = time.time()-startTime

	# close files
	REFFILE.close()
	OUTFQ1.close()
	if PAIRED_END:
		OUTFQ2.close()
	if SAVE_SAM:
		OUTSAM.close()
	if SAVE_VCF:
		OUTVCF.close()


	"""//////////////////////////////////////////////////
	////////////    WRITE RUNINFO OUTFILE    ////////////
	//////////////////////////////////////////////////"""

	if SAVE_RUNINFO:
		rfOut = open(OUTFILE_NAME+'_runInfo.txt','w')

		refSize = float(os.path.getsize(REFERENCE))
		fq1Size = float(os.path.getsize(OUTFILE_NAME+'_read1.fq'))
		if PAIRED_END:
			fq2Size = float(os.path.getsize(OUTFILE_NAME+'_read2.fq'))

		rfOut.write('Reference:\t\t'+REFERENCE+' ('+printSizeNicely(refSize)+')\n')
		if INPUT_BED == None:
			rfOut.write('# Sequences:\t'+str(len(ref_inds))+' ('+printBasesNicely(totalBPtargeted)+' in total)\n')
		else:
			rfOut.write('# Sequences:\t'+str(len(ref_inds))+' ('+str(sequencesSampledFrom)+' sampled from, '+printBasesNicely(totalBPtargeted)+' in total)\n')
		rfOut.write('RunDate:\t\t'+startDate+'\n')
		rfOut.write('Command:\t\t'+' '.join(sys.argv)+'\n')

		rfOut.write('\n\n********* FILES GENERATED *********\n\n')
		rfOut.write('Read files:\t\t'+OUTFILE_NAME+'_read1.fq ('+printSizeNicely(fq1Size)+')\n\t\t\t\t')
		if PAIRED_END:
			rfOut.write(OUTFILE_NAME+'_read2.fq ('+printSizeNicely(fq2Size)+')\n')
		if SAVE_CORRECT_REF:
			crSize = float(os.path.getsize(OUTFILE_NAME+'_correctRef.fa'))
			rfOut.write('\nRef + Variants:\t'+OUTFILE_NAME+'_correctRef.fa ('+printSizeNicely(crSize)+')\n')
		if SAVE_SAM:
			samSize = float(os.path.getsize(OUTFILE_NAME+'_golden.sam'))
			rfOut.write('\nGolden SAM:\t\t'+OUTFILE_NAME+'_golden.sam ('+printSizeNicely(samSize)+')\n')
		if SAVE_VCF:
			vcfSize = float(os.path.getsize(OUTFILE_NAME+'_golden.vcf'))
			rfOut.write('\nGolden VCF:\t\t'+OUTFILE_NAME+'_golden.vcf ('+printSizeNicely(vcfSize)+')\n')

		rfOut.write('\n\n********* PARAMETERS *********\n\n')
		rfOut.write('ReadLen:\t\t'+str(READLEN)+'\n')
		if PAIRED_END:
			rfOut.write('MeanFragLen:\t'+str(FRAGMENT_SIZE)+'\n')
			rfOut.write('FragLen Std:\t'+str(FRAGMENT_STD)+'\n')
		if INPUT_BED == None:
			rfOut.write('Avg Coverage:\t'+str(AVG_COVERAGE)+'x\n')
		else:
			rfOut.write('Avg Coverage:\t'+str(AVG_COVERAGE)+'x in targeted regions, '+str(BED_COVERAGE)+'x elsewhere\n')
		rfOut.write('Variant Freq:\t'+str(AVG_VAR_FREQ)+'\n')

		if SEQUENCING_SNPS:
			rfOut.write('SSE rate:\t\t'+str(AVG_SER)+'\n')
		else:
			rfOut.write('SSE rate:\t\t0.0\n')
		if SEQUENCING_INDELS:
			rfOut.write('SIE rate:\t\t'+str(SIE_RATE)+' of all sequencing errors\n')
		else:
			rfOut.write('SIE rate:\t\t0.0\n')
		if SEQUENCING_QSCORES:
			rfOut.write('QScore Model:\t'+QSCORE_MODEL+' : ['+str(Qscores[0])+','+str(Qscores[-1])+'] + '+str(qOffset)+'\n')
		else:
			rfOut.write('Qscore Model:\tNone (dummy value: '+str(DUMMY_QSCORE)+' + '+str(qOffset)+' )\n')
		rfOut.write('Mutation Model:\t'+MUT_MNAME+'\n')
		if INPUT_BED == None:
			rfOut.write('WindowSize:\t\t'+str(COV_WINDOW)+'\n')
		else:
			rfOut.write('Windowing:\t\t'+INPUT_BED+' ('+printBasesNicely(totalBedTargetedBP)+' targeted)\n')
		rfOut.write('RNG_SEED:\t\t'+str(RNG_SEED)+'\n')

		rfOut.write('\n\n********* STATS *********\n\n')
		if MULTI_JOB:
			rfOut.write('Total Reads:\n')
			rfOut.write('\t- This Job:\t\t\t'+str(readsFromThisJob)+' ('+printBasesNicely(readsFromThisJob*READLEN)+')\n')
			rfOut.write('\t- Among All Jobs:\t'+str(bigReadNameOffset)+' ('+printBasesNicely(bigReadNameOffset*READLEN)+')\n\n')
		else:
			rfOut.write('Total Reads:\t'+str(bigReadNameOffset)+' ('+printBasesNicely(bigReadNameOffset*READLEN)+')\n')
		rfOut.write('Total Runtime:\t'+str(int(finalTime))+' sec\n\n')
		#rfOut.write('Rate:\t\t\t'+printBasesNicely(int(float(bigReadNameOffset*READLEN)/finalTime+0.5))+'/sec\n\n')
		rfOut.write('Variants Introduced:\n')
		rfOut.write('\t- '+str(totalSNPs)+' SNPs, ['+str(len(snpAB))+' hets]\n')
		rfOut.write('\t- '+str(totalInds)+' small indels (length 1-'+str(MAX_INDEL)+') ['+str(len(indelAB))+' hets]\n')
		rfOut.write('\t- '+str(totalSTVs)+' SVs\n')

		rfOut.close()

if __name__ == '__main__':
	main()