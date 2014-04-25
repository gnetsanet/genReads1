#!/usr/bin/env python
# encoding: utf-8

""" **************************************************



************************************************** """

import shlex
from subprocess import Popen, PIPE
import time

def main():

	c1 = 'python genReads.py '
	c2 = ' -r testing/toy_chr21.fa -o testing/test_out --VCF --SAM --TXT'
	
	tt = time.time()
	print '1.) testing basic functionality...',
	cmd = c1+'-p -l 100'+c2
	process = Popen(shlex.split(cmd), stdout=PIPE, stderr=PIPE)
	process.communicate()
	exit_code = process.wait()
	print time.time()-tt,'['+str(exit_code)+']'

	tt = time.time()
	print '2.) long reads...',
	cmd = c1+'-p -l 1000 -f 2500 -F 50 -w 5000'+c2
	process = Popen(shlex.split(cmd), stdout=PIPE)
	process.communicate()
	exit_code1 = process.wait()
	cmd = c1+'-l 1000 -f 2500 -F 50 -w 5000'+c2
	process = Popen(shlex.split(cmd), stdout=PIPE)
	process.communicate()
	exit_code2 = process.wait()
	print time.time()-tt,'['+str(exit_code1)+']','['+str(exit_code2)+']'

	tt = time.time()
	print '3.) input vcf...',
	cmd = c1+'-p -v testing/input.vcf'+c2
	process = Popen(shlex.split(cmd), stdout=PIPE, stderr=PIPE)
	process.communicate()
	exit_code = process.wait()
	print time.time()-tt,'['+str(exit_code)+']'

	tt = time.time()
	print '4.) input target regions...',
	cmd = c1+'-p -b testing/bedTestRegions_toychr21.bed -c 100 -C 1'+c2
	process = Popen(shlex.split(cmd), stdout=PIPE, stderr=PIPE)
	process.communicate()
	exit_code = process.wait()
	print time.time()-tt,'['+str(exit_code)+']'

	tt = time.time()
	print '5.) multiple sequences in reference fasta...',
	cmd = c1+'-p -r testing/micro_chr21_double.fa -o testing/test_out --VCF --SAM --TXT'
	process = Popen(shlex.split(cmd), stdout=PIPE, stderr=PIPE)
	process.communicate()
	exit_code = process.wait()
	print time.time()-tt,'['+str(exit_code)+']'

	tt = time.time()
	print '6.) disabled errors, no qscores...',
	cmd = c1+'-V 0 -F 0 -S 0 --no-qscores'+c2
	process = Popen(shlex.split(cmd), stdout=PIPE, stderr=PIPE)
	process.communicate()
	exit_code = process.wait()
	print time.time()-tt,'['+str(exit_code)+']'

	tt = time.time()
	print '7.) parallel jobs...',
	cmd = c1+'--job 1 2 -R 12345'+c2
	process = Popen(shlex.split(cmd), stdout=PIPE, stderr=PIPE)
	process.communicate()
	exit_code1 = process.wait()
	cmd = c1+'--job 2 2 -R 12345'+c2
	process = Popen(shlex.split(cmd), stdout=PIPE, stderr=PIPE)
	process.communicate()
	exit_code2 = process.wait()
	cmd = 'python job_merge.py testing/test_out'
	process = Popen(shlex.split(cmd), stdout=PIPE, stderr=PIPE)
	process.communicate()
	exit_code3 = process.wait()
	print time.time()-tt,'['+str(exit_code1)+']','['+str(exit_code2)+']','['+str(exit_code3)+']'

	tt = time.time()
	print '8.) the kitchen sink...',
	cmd = c1+'-r testing/micro_chr21_double.fa -p -l 250 -f 600 -F 10 -w 2000 -c 20 -V 0.002 --no-ind -v testing/input2.vcf -R 12345 -o testing/test_out --VCF --SAM --TXT'
	process = Popen(shlex.split(cmd), stdout=PIPE, stderr=PIPE)
	process.communicate()
	exit_code = process.wait()
	print time.time()-tt,'['+str(exit_code)+']'

if __name__ == '__main__':
	main()