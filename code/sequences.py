#!/usr/bin/env python
# encoding: utf-8
"""
sequence.py

This is a class that will have info and methods for dealing with sequences.

Created by Måns Magnusson on 2012-03-10.
Copyright (c) 2012 __MyCompanyName__. All rights reserved.
"""

from __future__ import division
import sys
import os
import glob

class Sequences:
	"""Alleles holds information about the allele candidates found in a set of sequences."""
	def __init__(self, variable_positions_treshold = 0.2):
		self.seq_dict = {}
		self.variable_positions_treshold = variable_positions_treshold
		self.nr_of_seq = len(self.seq_dict)
		self.frequence_dict = {}
		self.variable_positions	= {}
		self.sequences_sorted = {}
		self.consensus = '' # The consensus sequence of the sequences in this object
		self.allele_1 = ''
		self.allele_2 = ''
	
	def add_sequence(self, sequence_id, sequence):
		"""Add a sequence to the sequence dict."""
		self.seq_dict[sequence_id] = sequence
		self.nr_of_seq = len(self.seq_dict)
	
	def make_own_frequencies(self):
		"""Make the frequence dict for all sequences"""
		self.frequence_dict = self.make_frequencies(self.seq_dict)
	
	def make_frequencies(self, sequence_dict):
		"""Make a dictionary with the frequencies of each position"""
		frequence_dict = {}
		pos_dict = {}
		nr_of_seq = len(sequence_dict)
		for ids, seq in sequence_dict.items():
			i = 1
			for nucleotide in seq:
				if i in pos_dict:
					pos_dict[i][nucleotide] += 1
				else:
					if nucleotide in ['A','a']:
						pos_dict[i] = {'A':1,'C':0,'G':0,'T':0,'-':0,'N':0}
					elif nucleotide in ['C','c']:
						pos_dict[i] = {'A':0,'C':1,'G':0,'T':0,'-':0,'N':0}
					elif nucleotide in ['G','g']:
						pos_dict[i] = {'A':0,'C':0,'G':1,'T':0,'-':0,'N':0}
					elif nucleotide in ['T','t']:
						pos_dict[i] = {'A':0,'C':0,'G':0,'T':1,'-':0,'N':0}
					elif nucleotide in ['-','N','n']:
						pos_dict[i] = {'A':0,'C':0,'G':0,'T':0,'-':1,'N':0}
				i += 1
		for pos,values in pos_dict.items():
			# To avoid division with zero:
			if nr_of_seq == 0:
				nr_of_seq	= 1
			A	= values['A']/nr_of_seq
			C	= values['C']/nr_of_seq
			G	= values['G']/nr_of_seq
			T	= values['T']/nr_of_seq
			N	= (values['N'] + values['-'])/nr_of_seq
			frequence_dict[pos]=[A,C,G,T,N]
		return frequence_dict
	
	def count_frequencies(self):
		"""Count the frequency for each position and print to stdout."""
		self.make_own_frequencies()
		self.find_variable_positions()
		for position in sorted(self.variable_positions.keys()):
			print position
		print 'Number of variable positions:', len(self.variable_positions)
	
	def make_seq_dict(self, seq_ids):
		"""Make a sequence dict from the sequence ids in a list."""
		seq_dict	= {}
		for sequence_id in seq_ids:
			seq_dict[sequence_id]	= self.seq_dict[sequence_id]
		return seq_dict
	
	def find_variable_positions(self):
		"""Picks out the variable positions, based on the frequence of the positions in one individual."""
		def get_variants(frequencies):
			"""Returns a string with the nucloetides, eg. 'AC'"""
			variants	= ''
			# Pick out the two most abundant nucleotides
			sorted_frequencies = sorted(frequencies[:4])[2:]
			if frequencies[0] in sorted_frequencies:
				if len(variants) < 2:
					variants	+= 'A'
			if frequencies[1] in sorted_frequencies:
				if len(variants) < 2:
					variants	+= 'C'
			if frequencies[2] in sorted_frequencies:
				if len(variants) < 2:
					variants	+= 'G'
			if frequencies[3] in sorted_frequencies:
				if len(variants) < 2:
					variants	+= 'T'
			return variants
		
		for pos in self.frequence_dict:
		#We look at the second max so we can see if there is distribution
		#between more than one nucleotide.
			if pos > 20 and pos < 280:# Avoid looking at the tags. Front tag is 16bp and end 17bp
			# Only look at the distribution between the nucleotides, not indels:
				second_max	= sorted(self.frequence_dict[pos][:4])[2]
				if second_max > self.variable_positions_treshold:
					variants	= get_variants(self.frequence_dict[pos])
					self.variable_positions[pos] = variants
					#Variants is a string with the two bases.
	
	def sort_sequences(self, position):
		"""Divide the sequences into groups based on the position given. Returns a list with two sequence dictionarys."""
		groups = {}
		for sequence_id in self.seq_dict:
			nucleotide = self.seq_dict[sequence_id][position-1]
			if nucleotide in self.variable_positions[position]:# Check if the nucleotide is one of the two variants
				if nucleotide in groups:
					groups[nucleotide].append(sequence_id)
				else:
					groups[nucleotide] = [sequence_id]
		return groups
	
	def make_consensus(self):
		"""Make the consensus sequence of self.sequences. If there are less than 10% of information for a position we let that one be unknown(insert -)."""
		if len(self.frequence_dict) == 0:
			self.make_frequencies()
		for pos, frequencies in self.frequence_dict.items():
			position	= 0
			highest		= 0
			i			= 0
			for frequency in frequencies:
				if frequency > highest:
					highest		= frequency
					position	= i
				i += 1
			if position == 0:
				self.consensus += 'A'
			elif position == 1:
				self.consensus += 'C'
			elif position == 2:
				self.consensus += 'G'
			elif position == 3:
				self.consensus += 'T'
			elif position in [4,5]:
				self.consensus += '-'
	
	def print_sequences(self):
		"""Prints the sequences to strdout"""
		for ids,seq in self.seq_dict.items():
			print ids
			print seq
			print len(seq)
			print ''
	
	def find_alleles(self):
		"""Find the alleles for these sequences."""
		
		def get_highest_frequency(frequencies):
			"""Returns the nucleotide with the highest frequency."""
			highest = 0
			highest_freq = 0
			for i in range(len(frequencies)):
				if frequencies[i] > highest_freq:
					highest = i
					highest_freq = frequencies[i]
			if highest == 0:
				return 'A'
			elif highest == 1:
				return 'C'
			elif highest == 2:
				return 'G'
			elif highest == 3:
				return 'T'
			elif highest == 4:
				return '-'
		
		self.make_own_frequencies()
		self.find_variable_positions()
		print sorted(self.variable_positions.keys())
		previous_position = 0
		# We look at every position in the sequence and add them to both alleles if not variable
		for position in range(1,len(self.frequence_dict)+1):
			# If the position is variable we sort the sequences based on the previous variable position
			if position in self.variable_positions:
				# If this is the first variable position we add one of each variation to the alleles
				if previous_position != 0:
					freq_dicts = []
					# Sort the sequences based on their variation in the previous variable position
					sequence_groups = self.sort_sequences(previous_position)
					for group in sequence_groups:
						# Turn the dictionary into a sequence dict
						sequence_dict = self.make_seq_dict(sequence_groups[group])
						# Make new frequences of the new sequence dict
						freq_dicts.append(self.make_frequencies(sequence_dict))
					# If the previous variable position for this group is the same as in allele 1 then add the new nucleotide for this group to allele 1, else to allele 2.
					# Find which of the groups that have the strongest vote for its nucleotide:
					if sorted(freq_dicts[0][position])[-1] > sorted(freq_dicts[1][position])[-1]:
						group_1_nucleotide = get_highest_frequency(freq_dicts[0][position])
						if group_1_nucleotide == self.variable_positions[position][0]:
							group_2_nucleotide = self.variable_positions[position][1]
						else:
							group_2_nucleotide = self.variable_positions[position][0]
					else:
						group_2_nucleotide = get_highest_frequency(freq_dicts[1][position])
						if group_2_nucleotide == self.variable_positions[position][0]:
							group_1_nucleotide = self.variable_positions[position][1]
						else:
							group_1_nucleotide = self.variable_positions[position][0]
					if get_highest_frequency(freq_dicts[0][previous_position]) == self.allele_1[previous_position-1]:
						self.allele_1 += group_1_nucleotide
						self.allele_2 += group_2_nucleotide
					else:
						self.allele_1 += group_2_nucleotide
						self.allele_2 += group_1_nucleotide

					previous_position = position
				else:
					# If this is the first variable position, add one variant to each.
					self.allele_1 += self.variable_positions[position][0]
					self.allele_2 += self.variable_positions[position][1]
					previous_position = position
			else:
				# If not variable position we add the nucleotide with the highest frequency:
				nucleotide = get_highest_frequency(self.frequence_dict[position])
				self.allele_1 += nucleotide
				self.allele_2 += nucleotide
	
	def print_fasta(self, fasta_file, ind_id):
		"""Print the sequences to the end of a fasta file"""
		f = open(fasta_file, 'a')
		f.write('>'+ ind_id + '_1\n')
		f.write(self.allele_1+'\n')
		f.write('>'+ ind_id + '_2\n')
		f.write(self.allele_2+'\n')


def main():
	pass


if __name__ == '__main__':
	main()

