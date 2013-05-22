#!/usr/bin/env python
# encoding: utf-8
"""
fasta_parser.py

Class for make .fasta parsers. Creates sequence objects out of each line on the form {<seq_id>: sequence}

Created by Måns Magnusson on 2013-04-22.
Copyright (c) 2013 __MyCompanyName__. All rights reserved.
"""

import sys
import os
import glob
import argparse
from DLA_Genotyper.sequences import Sequences

class Fasta_Parser(object):
	"""Parse a fasta file"""
	def __init__(self, fasta_file, variable_positions_treshold=0.3):
		super(Fasta_Parser, self).__init__()
		self.ind_id = os.path.basename(fasta_file)[:-3]
		self.sequences = Sequences(variable_positions_treshold=variable_positions_treshold)
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
							self.sequences.add_sequence(seq_id, sequence)
							seq_id = line[1:]
							sequence = ""
						else: 
							sequence += line
			self.sequences.add_sequence(seq_id, sequence)
	
	def find_alleles(self):
		"""Find the alleles of this individual"""
		self.sequences.find_alleles()
	
	def count_frequencies(self):
		"""Count the frequencies for each position and print to stdout"""
		self.sequences.count_frequencies()
	
	def get_sequences(self):
		"""Return the sequence obejct"""
		return self.sequences
	
	def write_alleles(self, out_file):
		"""Print the alleles to a fasta file"""
		self.sequences.print_fasta(out_file, self.ind_id)

def write_alleles(file_handle, ind_id, sequence_1, sequence_2):
	"""Write the alleles to a file in fasta format"""
	file_handle.write('>'+ ind_id + '_1\n')
	file_handle.write(sequence_1+'\n')
	file_handle.write('>'+ ind_id + '_2\n')
	file_handle.write(sequence_2+'\n')
	
def main():
	parser = argparse.ArgumentParser(description="Put the fasta files in a dictionar")
	parser.add_argument('fasta_file', type=str, nargs=1, help='Specify the the path to a fasta file containing sequences or a directory with fasta files.')
	parser.add_argument('-write_alleles', '--write_alleles', type=str, nargs=1, help='Specify the path to a fastafile where we write the results.')
	parser.add_argument('-find_alleles', '--find_alleles', action='store_true', help='Find the alleles for the individual/individuals')	
	parser.add_argument('-count_frequencies', '--count_frequencies', action='store_true', help='Count the frequencies for the positions in the give fasta file')	
	args = parser.parse_args()
	path_indata = args.fasta_file[0]
	if args.write_alleles:
		file_handle = open(args.write_alleles[0], 'a')
	if os.path.isdir(path_indata):
		for infile in glob.glob(os.path.join(path_indata, '*.fa')):
			if os.path.getsize(infile) > 0:
				my_fasta_sequences = Fasta_Parser(infile)
				print 'individual:', my_fasta_sequences.ind_id
				if args.find_alleles:
					my_fasta_sequences.find_alleles()
					if args.write_alleles:
						write_alleles(file_handle, my_fasta_sequences.ind_id, my_fasta_sequences.sequences.allele_1, my_fasta_sequences.sequences.allele_2)
				elif args.count_frequencies:
					my_fasta_sequences.count_frequencies()
	else:
		infile = path_indata
		if args.find_alleles:
			my_fasta_sequences = Fasta_Parser(infile)
			my_fasta_sequences.find_alleles()
			if args.write_alleles:
				write_alleles(file_handle, my_fasta_sequences.ind_id, my_fasta_sequences.sequences.allele_1, my_fasta_sequences.sequences.allele_2)
		elif args.count_frequencies:
			my_fasta_sequences = Fasta_Parser(infile, variable_positions_treshold=0.00001)
			my_fasta_sequences.count_frequencies()

if __name__ == '__main__':
	main()

