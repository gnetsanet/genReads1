#!/usr/bin/env python
# encoding: utf-8

""" **************************************************


************************************************** """

import sys
import os

def main():

	if len(sys.argv) != 2:
		print '\npython repMerge.py out_prefix\n'
		exit(1)
	else:
		DIR = sys.argv[1]
		op  = DIR.split('/')[-1]
		dr  = '/'.join(DIR.split('/')[:-1])+'/'
		print dr
		if dr == '/':
			dr = './'

	listing = os.listdir(dr)

	locs = [n for n in listing if 'loc' in n and op+'_' == n[:len(op)+1]]
	coun = [n for n in listing if 'count' in n and op+'_' == n[:len(op)+1]]
	#print locs
	#print coun

	cmd = 'cat '
	for n in locs:
		cmd = cmd + dr + n + ' '
	cmd = cmd + '> ' + dr + op + '_mergedLoc.txt'

	print '\n'+cmd+'\n'
	os.system(cmd)

	cmd = 'cat '
	for n in coun:
		cmd = cmd + dr + n + ' '
	cmd = cmd + '> ' + dr + op + '_mergedCount.txt'

	print '\n'+cmd+'\n'
	os.system(cmd)

if __name__ == '__main__':
	main()