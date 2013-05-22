#!/usr/bin/env python
# encoding: utf-8
"""
add_tags.py

Add the complete tag sequences to all references.

Created by Måns Magnusson on 2012-03-02.
Copyright (c) 2012 __MyCompanyName__. All rights reserved.
"""

import sys
import os

def get_seq(infile):
	"""Put all sequences in a dictionary on the form:{<seqId>:<sequence>}
	Input: A FASTA file with the sequences for an individual
	Output: A dictionary on the form {<seqId>:<sequence>}"""

	individual_dict={}
	seq_id=""
	sequence=""
	j=0
	for line in open(infile,'r'):
		line=line.rstrip()
		if len(line)>0:
			if line[0] == ">":
				if j != 0:
					individual_dict[seq_id]=sequence
					seq_id=line[1:]
					sequence = ''
				elif j == 0:
					seq_id=line[1:]
					j+=1
			elif line[0] != ">": 
				sequence+=line
	individual_dict[seq_id]=sequence
	return individual_dict

def add_tags(sequences,front_tag,tail_tag,out_file):
	"""Adds sequences to fron and back of all sequences in the dictionary and prints them to the outfile."""
	new_sequences	= {}
	f				= open(out_file,'w')
	for ids, seq in sequences.items():
		new_seq 			= front_tag + seq + tail_tag
		print >>f,'>'+ids
		print >>f,new_seq

def main():
	ref_file		= sys.argv[1]
	ref_file_out	= sys.argv[2]
	front_tag		= 'CGCCCGCTGCGCTCAC'
	tail_tag		= 'CTGTGGGGACGGGGGATC'
		
	ref_sequences	= get_seq(ref_file)
	add_tags(ref_sequences,front_tag,tail_tag,ref_file_out)
	


if __name__ == '__main__':
	main()

