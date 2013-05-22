#!/usr/bin/env python
# encoding: utf-8
"""
parse_individuals.py

Created by Måns Magnusson on 2013-05-12.
Copyright (c) 2013 __MyCompanyName__. All rights reserved.
"""

import sys
import os

def parse_individuals(ref_file):
	"""Create a dictionary with the found sequences on the form {<seq_id>:[sequence_1, sequence_2}"""
	individuals = {}
	number_of_individuals = 0
	with open(ref_file, 'r') as f:
		seq_id = ""
		sequence = ""
		beginning = True
		for line in f:
			line = line.rstrip()
			if len(line)>0:
				if beginning:
					seq_id = line[1:-2]
					beginning = False
				else:
					if line[0] == ">":
						if seq_id in individuals:
							individuals[seq_id].append(sequence)
						else:
							individuals[seq_id] = [sequence]
							number_of_individuals += 1
						seq_id = line[1:-2]
						sequence = ""
					else: 
						sequence += line
		if seq_id in individuals:
			individuals[seq_id].append(sequence)
		else:
			individuals[seq_id] = [sequence]
	print 'Number of individuals = ', number_of_individuals
	return individuals


def main():
	pass


if __name__ == '__main__':
	main()

