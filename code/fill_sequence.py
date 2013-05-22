#!/usr/bin/env python
# encoding: utf-8
"""
fill_sequence.py

Created by Måns Magnusson on 2013-04-25.
Copyright (c) 2013 __MyCompanyName__. All rights reserved.
"""

import sys
import os
import argparse

def parse_references(ref_file):
	"""Create a dictionary with the reference sequences on the form {<seq_id>:sequence}"""
	references = {}
	with open(ref_file, 'r') as f:
		seq_id = ""
		sequence = ""
		beginning = True
		for line in f:
			line = line.rstrip()
			if len(line)>0:
				if beginning:
					seq_id = line[1:]
					beginning = False
				else:
					if line[0] == ">":
						references[seq_id] = sequence
						seq_id = line[1:]
						sequence = ""
					else: 
						sequence += line
		references[seq_id] = sequence
	return references

def check_similar(sequence, references):
	"""Check which of the references that are most similar to the sequence."""
	best_score = 0
	best_reference = ''
	for seq_id in references:
		score = 0
		for position in range(len(sequence)):
			if sequence[position] == references[seq_id][position]:
				score += 1
		if score > best_score:
			best_score = score
			best_reference = seq_id
	print best_score, best_reference
	return best_reference

def fill_sequence(sequence, closest_reference):
	"""Fill the gaps with the corresponding nucleotide from the closest reference sequence."""
	for position in range(len(sequence)):
		if sequence[position] == '-':
			sequence = sequence[:position] + closest_reference[position] + sequence[position+1:]
	return sequence

def fill_sequences(fasta_file, references, out_file):
	"""Read the sequences from the fasta file and fill the gaps with corresponding nucleotide from the most similar reference sequence."""
	out_handle = open(out_file, 'w')
	with open(fasta_file, 'r') as f:
		seq_id = ""
		sequence = ""
		beginning = True
		for line in f:
			line = line.rstrip()
			if len(line)>0:
				if beginning:
					seq_id = line[1:]
					beginning = False
				else:
					if line[0] == ">":
						closest_reference = check_similar(sequence, references)
						filled_sequence = fill_sequence(sequence, references[closest_reference])
						out_handle.write('>'+seq_id + '\n')
						out_handle.write(filled_sequence + '\n')
						seq_id = line[1:]
						sequence = ""
					else: 
						sequence += line
		closest_reference = check_similar(sequence, references)
		filled_sequence = fill_sequence(sequence, references[closest_reference])
		out_handle.write('>' + seq_id + '\n')
		out_handle.write(filled_sequence + '\n')
	out_handle.close()

def main():
	parser = argparse.ArgumentParser(description="Fill the indels with the nucleotide from the most similar reference sequence.")
	parser.add_argument('fasta_file', type=str, nargs=1, help='Specify the the path to a fasta file containing sequences or a directory with fasta files.')
	parser.add_argument('references', type=str, nargs=1, help='Specify the the path to a fasta file containing reference sequences.')
	parser.add_argument('outfile', type=str, nargs=1, help='Specify the path to a fastafile where we write the results.')
	args = parser.parse_args()
	references = parse_references(args.references[0])
	fill_sequences(args.fasta_file[0], references, args.outfile[0])

if __name__ == '__main__':
	main()

