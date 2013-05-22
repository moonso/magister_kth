#!/usr/bin/env python
# encoding: utf-8
"""
Sequences.py

This is a class with the information about an individual.

So I would like to re write the find allele section in such a way that we divide the sequences into TWO groups.
The alleles are then built in a way such that in every variable position each of the groups have to choose between one of the two candidates and nothing else.
This will result in that the groups will have the maximum humming distance between eachother.

We should divide the sequences into two groups based on the position where the distribution is most even!


Questions:

1. What happends if both groups have the same most abundant nucleotide in one position?
Answer:
We can look at the distribution of the other nucletide.

2. What happends if both groups have the same distribution of both nucleotides in one position? 
Answer:
Maybe we can discard this individual as impossible 

Not one single variable nucleotide had this problem, that means if we divide the sequences into two groups based on what type they have in the first variable position we could just make the consensus sequence out of them to get the true alleles!


Created by Måns Magnusson on 2012-03-05.
Copyright (c) 2012 __MyCompanyName__. All rights reserved.
"""

from __future__ import division
import sys
import os
import Sequences
import numpy as np

class Individual:
	"""Sequences"""
	def __init__(self, individual_id, sequences, treshold = 0.35):
		"""Put all sequences in a dictionary on the form:{<seqId>:<sequence>}"""
		self.individual				= individual_id# This is the numerical id of this individual
		self.heterozygote			= True
		self.sequences				= sequences# Sequences object with the sequences for the ind. like {<seqId>:<sequence>}
		# self.old_alleles			= old_alleles# Dict. with the known sequences like {<seqId>:<sequence>}
		self.var_pos_treshold		= treshold
		self.allele_candidates		= [] # Holds two Sequence objects, one for each set of sequences.
		self.allele_1				= ''
		self.allele_2				= ''
		self.treshold				= treshold
		self.allele_candidates		= []
		self.sequences.make_frequencies()
	
	def get_pos_candidates(self, frequencies):
		"""Returns <first> and <second> which are tuples on the form ('A',0.87)"""
		base_one	= ''
		base_two	= ''
		freq_one	= 0.0
		freq_two	= 0.0
		i			= 0
		for frequence in frequencies:
			if frequence > freq_one:
				freq_two	= freq_one
				base_two	= base_one
				freq_one	= frequence
				if i == 0:
					base_one = 'A'
				elif i == 1:
					base_one = 'C'
				elif i == 2:
					base_one = 'G'
				elif i == 3:
					base_one = 'T'
			elif frequence > freq_two:
				freq_two	= frequence
			 	if i == 0:
			 		base_two = 'A'
			 	elif i == 1:
			 		base_two = 'C'
			 	elif i == 2:
			 		base_two = 'G'
			 	elif i == 3:
			 		base_two = 'T'
			i+=1
		return [(base_one,freq_one),(base_two,freq_two)]
	
	def check_variable_positions(self):
		"""Check if there are any contradictions of the variable positions when dividing the sequences into two groups."""
		self.sequences.find_variable_positions(self.treshold)
		sequence_groups	= self.sequences.sort_sequences()
		for group in sequence_groups:
			self.allele_candidates.append(Sequences.Sequences(group,'dict'))
		for position in self.sequences.variable_positions:
			group_one	= self.get_pos_candidates(self.allele_candidates[0].frequence_dict[position])
			group_two	= self.get_pos_candidates(self.allele_candidates[1].frequence_dict[position])
			if group_one[0][0] == group_two[0][0]:
				print self.individual, position
				print group_one
				print group_two
			# print position, self.sequences.variable_positions[position], self.allele_candidates[0].frequence_dict[position], self.allele_candidates[1].frequence_dict[position], self.sequences.frequence_dict[position]
	
	def find_best_ref(self,sequence_object_one, sequence_object_two):
		"""Finds and returns the sequence from sequence_object_one that is most like the consensus in sequence_object_two."""
		best_score		= 0
		best_sequence	= ''
		for sequence in sequence_object_one.seq_dict.values():
			score	= sequence_object_two.hamming_distance(sequence_object_two.consensus,sequence)
			if score > best_score:
				best_score		= score
				best_sequence	= sequence
		return best_sequence
	
	def fill_consensus(self, ref, sequence_object):
		"""Fills the gaps in the consensensus sequence with the corresponding base in the reference sequence."""
		i = 0
		new_consensus	= sequence_object.consensus
		for position in sequence_object.consensus:
			if position == '-':
		 		new_consensus = new_consensus[:i] + ref[i] + new_consensus[i+1:]
		 	i += 1
		sequence_object.consensus	= new_consensus
	
	def find_alleles(self):
		"""Find the alleles for the individual."""
		self.sequences.find_variable_positions(self.var_pos_treshold)
		sequence_groups	= self.sequences.sort_sequences()
		test_file		= open('my_test.fa','w')
		for group in sequence_groups:
			for key, value in group.items():
				test_file.write('>'+key+'\n')
				test_file.write(value+'\n')
		test_file.close()
		# for group in sequence_groups:
		# 	self.allele_candidates.append(Sequences.Sequences(group,'dict'))
		# for candidate in self.allele_candidates:
		# 	candidate.make_consensus()
		# 	best_ref	= self.find_best_ref(self.old_alleles, candidate)
		# 	self.fill_consensus(best_ref, candidate)
		# self.allele_1	= self.allele_candidates[0].consensus
		# self.allele_2	= self.allele_candidates[1].consensus
	
	def get_zygosity(self):
		"""Get the zygosity for the individual."""
		self.sequences.find_variable_positions(self.var_pos_treshold)
		self.nr_of_var_positions = len(self.sequences.variable_positions)
		if len(self.sequences.variable_positions) > 1:
			self.heterozygote		= True
		else:
			self.heterozygote		= False
	
	def print_sequences(self, mode='normal'):
		"""Prints the sequences to strdout"""
		if mode == 'normal':
			for ids,seq in self.sequences.seq_dict.items():
				print ids
				print seq
				print len(seq)
				print ''
		elif mode == 'old_alleles':
			for ids,seq in self.old_alleles.seq_dict.items():
				print ids
				print seq
				print len(seq)
				print ''
	
	def print_alleles(self, file_handle=False):
		"""Print the alleles to a file in the FASTA format."""
		if file_handle:
			file_handle.write('>'+self.individual+'\t'+'1'+'\n')
			file_handle.write(self.allele_1+'\n')
			file_handle.write('\n')
			file_handle.write('>'+self.individual+'\t'+'2'+'\n')
			file_handle.write(self.allele_2+'\n')
			file_handle.write('\n')
		else:
			print '>'+self.individual+'\t'+'1'+'\n'
			print self.allele_1+'\n'
			print '>'+self.individual+'\t'+'2'+'\n'
			print self.allele_2+'\n'	
	

def printAlleles(path_outdata, samples):
	"""Print the alleles to a file in the FASTA format."""
	f = open(path_outdata,'w')
	for sample in samples:
		i = 1
		for allele, sequence in sample.alleles.consensus_sequences.items():
			f.write('>'+sample.individual+' '+str(i)+'\n')
			f.write(sequence[16:-18]+'\n')
			i	+= 1
	f.close()


def main():
	for tres in tresholds:
		heterozygotes		= 0
		homozygotes			= 0
		mean_nr_of_var_pos	= 0
		for sample in samples:
			sample.get_zygosity(tres)
			mean_nr_of_var_pos	+= sample.nr_of_var_positions
			if sample.heterozygote:
				heterozygotes	+= 1
			else:
				homozygotes		+= 1
		print 'Treshold:', tres
		print 'Heterozygotes:', heterozygotes
		print 'Homozygotes:	', homozygotes
		print 'Mean number of variable positions:', mean_nr_of_var_pos/len(samples)
		

if __name__ == '__main__':
	main()

