#!/usr/bin/env python
# encoding: utf-8
"""
check_homozygotes.py

Count some metrics for the called alleles.

Created by Måns Magnusson on 2013-05-11.
Copyright (c) 2013 __MyCompanyName__. All rights reserved.
"""

from __future__ import division
import sys
import os
import argparse
from parse_individuals import parse_individuals

def compare_sequences(seq_1, seq_2):
	"""Look if both sequences are the same, '-' counts as a wild card. Sequences needs to be of same length."""
	if len(seq_1) != len(seq_2):
		return 100
	difference = 0
	for position in range(len(seq_1)):
		if seq_1[position] != '-' and seq_2[position] != '-':
			if seq_1[position] != seq_2[position]:
				difference += 1
	return difference

def parse_allele_list(alle_list):
	"""Parse a file with <individual>\t<allele_1>\t(<allele_2>)"""
	individuals = {}
	with open(alle_list, 'r') as f:
		for line in f:
			line = line.rstrip().split()
			individuals[line[0]] = line[1:]
	return individuals

def main():
	parser = argparse.ArgumentParser(description="Print the information for each individual to a specified outfile.")
	parser.add_argument('fasta_file', type=str, nargs=1, help='Specify the the path to a fasta file containing the found sequences.')
	parser.add_argument('-allele_list', '--allele_list', type=str, nargs=1, help='Specify the the path to a fasta file containing the found sequences.')
	args = parser.parse_args()
	# Get all the alleles in a dictionary, {<ind_id>:[sequence_1, sequence_2]}:
	found_alleles = parse_individuals(args.fasta_file[0])
	homozygotes = 0
	for ind, sequences in found_alleles.items():
		if compare_sequences(sequences[0], sequences[1]) == 0:
			homozygotes += 1
	print 'Number of homozygotes:', homozygotes 
	print 'Fraction of homozygotes:', homozygotes/len(found_alleles)
	
	
	if args.allele_list:
		individual_alleles = parse_allele_list(args.allele_list[0])
		new_homozygotes = {}
		for ind, alleles in individual_alleles.items():
			if len(alleles) == 1:
				allele = alleles[0]
				if 'New' in allele:
					if allele in new_homozygotes:
						new_homozygotes[allele].append(ind)
					else:
						new_homozygotes[allele] = [ind]
		for allele, individuals in new_homozygotes.items():
			print allele, individuals
		print 'Number of new alleles found in homozygotes:', len(new_homozygotes)
		
		different_positions = {}
		for ind, alleles in individual_alleles.items():
			if len(alleles) == 2:
				allele_1 = False 
				allele_2 = False
				if 'New' in alleles[0]:
					allele_1 = True
				if 'New' in alleles[1]:
					allele_2 = True
				if allele_1 != allele_2:		
					difference = compare_sequences(found_alleles[ind][0], found_alleles[ind][1])
					if difference in different_positions:
						different_positions[difference].append(ind)
					else:
						different_positions[difference] = [ind]
		for difference in sorted(different_positions.keys()):
			print difference, len(different_positions[difference])


if __name__ == '__main__':
	main()

