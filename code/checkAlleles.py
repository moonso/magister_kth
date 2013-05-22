#!/usr/bin/env python
# encoding: utf-8
"""
checkAlleles.py

Take the alleles found in the individuals and create two dictionarys one on the form {<sequence> : alleleID} and one with {<individualNr> : [alleleId:s]}, also print to a file on the form:

nr	ID	Species	Type(homo/hetero)	Allele1		Allele2


This script alter the alleles in way in the way that we copy and paste the problematic poly-c region in the exon.

Created by Måns Magnusson on 2011-04-13.
Copyright (c) 2011 __MyCompanyName__. All rights reserved.
"""

import sys
import os
import argparse
import operator
from parse_individuals import parse_individuals

def parse_references(ref_file):
	"""Create a dictionary with the reference sequences on the form {<sequence>:seq_id}"""
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
						references[sequence] = seq_id
						seq_id = line[1:]
						sequence = ""
					else: 
						sequence += line
		references[sequence] = seq_id
	return references


def main():
	parser = argparse.ArgumentParser(description="Print the information for each individual to a specified outfile.")
	parser.add_argument('fasta_file', type=str, nargs=1, help='Specify the the path to a fasta file containing the found sequences.')
	parser.add_argument('references', type=str, nargs=1, help='Specify the the path to a fasta file containing reference sequences.')
	parser.add_argument('out_file', type=str, nargs=1, help='Specify the path to a fastafile where we write the results.')
	args = parser.parse_args()
	# Get all the alleles in a dictionary, {<ind_id>:[sequence_1, sequence_2]}:
	found_alleles = parse_individuals(args.fasta_file[0])
	out_file = args.out_file[0]
	#This is the known alleles from litterature, {<sequence>: seq_id}:
	references = parse_references(args.references[0])
	out_handle = open(out_file, 'w')
	new = 1 # Keep track of the number of new alleles
	ref_count = {} # Keep track of how many times the sequences are seen.
	print 'Number of reference squences:', len(references)
	
	with open(out_file, 'w') as f:
		for individual, sequence_list in found_alleles.items():
			print_line = individual+'\t' # Build the outdata for this individual
			if sequence_list[0] == sequence_list[1]:
				sequence_list = sequence_list[:1] # If homozygote, just use one allele
			for sequence in sequence_list:
			# Check if the sequence has been observed before:
				if sequence in references:
					seq_id = references[sequence]
				else:
					seq_id = 'New_' + str(new)
					new += 1
					# If we find a new sequence we add them to the references
					references[sequence] = seq_id
				# Add the allele count for this sequence:
				if seq_id in ref_count:
					ref_count[seq_id] += 1
				else:
					ref_count[seq_id] = 1
				print_line += seq_id + '\t'
			f.write(print_line + '\n')
	
	# Print all sequences seen and the number of times they where observed:		
	with open('all_sequences.txt', 'w') as f:
		for ref_name in ref_count:
			f.write(ref_name + '\t' + str(ref_count[ref_name]) + '\n')
			
	allele_counter = 0 # How many alleles where seen in the data set
	max_count = 0 # How many times did we see the most abundant allele
	max_allele = '' # The allele that was seen most of the times
	for seq_id, count in ref_count.items():
		allele_counter += count
		if count > max_count:
			max_count = count
			max_allele = seq_id
	# Produce a list of sorted tuples:
	sorted_alleles = sorted(ref_count.iteritems(), key=operator.itemgetter(1))
	count_dict = {}
	new_alleles = 0
	old_alleles = 0
	for i in sorted_alleles:
		alle_name = i[0]
		allele_count = i[1]
		if 'New' in alle_name:
			new_alleles += 1
		else:
			old_alleles += 1
		if allele_count in count_dict:
			count_dict[allele_count] += 1
		else:
			count_dict[allele_count] = 1
	for count, number_of_alleles in count_dict.items():
		print count, number_of_alleles
	print 'Number of New alleles:', new_alleles
	print 'Number of old alleles:', old_alleles
	print 'Number of seen unique alleles:', len(ref_count)
	print 'Number of alleles seen in our data:', allele_counter
	print 'Number of individuals with the max allele:', max_count
	print  'Max allele:', max_allele
	# print >>f, "nr".center(6), "Id".center(8), "Zygote".center(8), "Allele1".center(12), "Allele2".center(12), "Method".center(8), "NrOfSeq".center(10)
	# #
	# #
	# #
	# individuals = read_ids(individual_file)
	# old_alleles = make_seq_dict(known_alleles)#Seq is key and id is value.
	# new_nr = 1
	# ind_dict = make_ind_dict(found_alleles)
	# #
	# #
	# #
	# for ind_id, alleles in ind_dict.items():
	# 	for allele in alleles:
	# 		if allele not in old_alleles:
	# 			allele_name 		= 'New'+str(new_nr)
	# 			old_alleles[allele] = allele_name
	# 			new_nr 				+= 1
	# 			g					= open(known_alleles,'a')
	# 			g.write('>'+allele_name+"\n")
	# 			g.write(allele+"\n")
	# 			g.close() 
	# 	getID			  = '../results/2012-03-04/sequences_aligned_noref/' + ind_id + '.fa'
	# 	getIndDict		  = get_seq(getID) 
	# 	sampleType 		  = getType(ind_id,individuals)
	# 	if len(alleles) == 1:
	# 		zygosity = 'Homo'
	# 		allele1 = old_alleles[alleles[0]]
	# 		allele2 = '-'
	# 	else:
	# 		zygosity = 'Hetero'
	# 		allele1 = old_alleles[alleles[0]]
	# 		allele2 = old_alleles[alleles[1]]
	# 	
	# 	print >>f, ind_id.center(6), individuals[ind_id].center(8), zygosity.center(8), allele1.center(12), allele2.center(12), sampleType.center(8), str(len(getIndDict)).center(10)
	# f.close()


if __name__ == '__main__':
	main()

