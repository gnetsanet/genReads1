#!/usr/bin/env python
# encoding: utf-8

import cPickle as pickle

def main():

	refName = 'chr21'
	LEN     = 280

	f = open('repeat_locations_'+refName+'.txt','r')

	badRegions = []

	for line in f:
		splt = line.split('\t')
		if splt[0] == '>'+str(LEN):
			#print splt
			for i in xrange(1,len(splt),2):
				n1 = int(splt[i])
				n2 = int(splt[i+1])
				if len(badRegions) == 0:
					badRegions.append((n1,n2))
				else:
					# ASSUMES REGIONS LIST IS SORTED
					if n1 < badRegions[-1][1]:
						if n2 > badRegions[-1][1]:
							badRegions[-1] = (badRegions[-1][0],n2)
					else:
						badRegions.append((n1,n2))
	
	pickle.dump([badRegions],open('badRegions_'+refName+'_l'+str(LEN)+'.p','wb'))

	f.close()

if __name__ == '__main__':
	main()