#!/usr/bin/env python
# encoding: utf-8
"""
check_homopolymers.py

Count the number of homopolymers in the problematic middle region of the sequence.

Created by Måns Magnusson on 2012-11-30.
Copyright (c) 2012 __MyCompanyName__. All rights reserved.
"""

from __future__ import division
import glob
import sys
import os
import numpy as np
import matplotlib.pyplot as plt



def main():
	dir_path		= sys.argv[1]
	homopolymers	= {'CCCCC':0,'CCCC':0,'CCC':0,'CC':0,'C':0,'no':0}
	nr_of_ind		= 0
	nr_of_seq		= 0
	max_seq			= 0
	sequences		= {}
	for infile in glob.glob(os.path.join(dir_path,'*.*')):
		nr_of_ind	+=1
		ind_seq		= 0
		for line in open(infile,'r'):
			if len(line)>0 and line[0]!='>':
				ind_seq		+= 1
				nr_of_seq	+= 1
				if 'CCCCC' in line[160:170]:
					homopolymers['CCCCC'] += 1
				elif 'CCCC' in line[160:170]:
					homopolymers['CCCC'] += 1
				elif 'CCC' in line[160:170]:
					homopolymers['CCC'] += 1
				elif 'CC' in line[160:170]:
					homopolymers['CC'] += 1
				elif 'C' in line[160:170]:
					homopolymers['C'] += 1
				else:
					homopolymers['no'] += 1
		if ind_seq > max_seq:
			max_seq = ind_seq
		if ind_seq in sequences:
			sequences[ind_seq] += 1
		else:
			sequences[ind_seq] = 1
	print sequences
	median		= 0
	median_seq	= 0
	for key, value in sequences.items():
		if value > median_seq:
			median_seq	= value
			median		= key
	for key,value in homopolymers.items():
		print key, value/nr_of_seq
	print 'nr_of_seq:',nr_of_seq
	print 'nr_of_ind:',nr_of_ind
	print 'Medium number of sequences:',nr_of_seq/nr_of_ind
	print 'Median number of sequences:',median
	print 'Maximum number of sequences for one individual',max_seq


if __name__ == '__main__':
	main()

