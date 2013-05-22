#!/usr/bin/env python
# encoding: utf-8
"""
localParAlign.py

This is the first attempt to write my own local pairwise aligner that will be easy to tune. For example so we can adjust the score to conserved regions in the reference sequences.

Edit 111102:

This is the version that is beeing used!!! The algorithm is not tweeked but i've made some small changes to opening of files and stuff.

Created by Måns Magnusson on 2011-06-22.
Copyright (c) 2011 __MyCompanyName__. All rights reserved.
"""

import sys
import os
import numpy as np
import glob

def getId(infile):
	"""Returns a string that is the id of a individual.
	Input: String with path to a file. 
	Output: The individuals number"""
	
	newId=""
	i=-1
	dot=0
	while infile[i]!="/":
		if dot >= 1:
			newId+=infile[i]
			i-=1
		else:
			if infile[i]==".":
				dot += 1 
				i-=1
			else:
				i-=1
	return newId[::-1]

def getSeq(infile):
	"""Put all sequences in a dictionary on the form:{<seqId>:<sequence>}
	Input: A FASTA file with the sequences for an individual
	Output: A dictionary on the form {<seqId>:<sequence>}"""

	individualDict={}
	f=open(infile, 'r')
	seqId=""
	sequence=""
	j=0
	for line in f:
		line=line.rstrip()
		if len(line)>0:
			if line[0] == ">":
				if j != 0:
					individualDict[seqId]=sequence
					seqId=line[1:]
					sequence = ''
				elif j == 0:
					seqId=line[1:]
					j+=1
			elif line[0] != ">": 
				sequence+=line
	individualDict[seqId]=sequence
	f.close()
	return individualDict


def makeMatrix(dicti,length):
	"""Turns the table into a matrix"""
	myMatrix = np.zeros(shape=(6,length))
	f = open(dicti,'r')
	for line in f:
		line	= line.rstrip()
		words	= line.split()
		pos 	= int(words[0])
		j		= 0
		for i in words:
			if j != 0:
				if int(i) != 0:
					myMatrix[j-1][pos-1] = int(i)
			j+=1
	f.close()
	return myMatrix

def checkScore(refPos,typ,scoreMatrix):
	"""Returns the score of the specific match, det är i denna funktionen som jag kan börja ändra saker. Måste gå på postionsnummer och därmed ha en lista med scores för varje position. Undrar om man kan ha en tabell vid sidan av som scriptet kan titta i eller om jag måste ha den med i scriptet? Bästa måste helt enkelt vara ett dictionary med en lista för varje position. Skulle vara snyggt att ha ett script som läser in alla tillgängliga referenssekvenser och bygger sedan ett dictionary"""
	#Dictionaryt behöver endast innehålla 270 poster eftersom referensen är av denna längd, insertioner kommer ju att ske senare
	if typ == 'A':
		return scoreMatrix[0][refPos]
	elif typ == 'C':
		return scoreMatrix[1][refPos]
	elif typ == 'G':
		return scoreMatrix[2][refPos]
	elif typ == 'T':
		return scoreMatrix[3][refPos]
	elif typ == 'gapSeq':
		return scoreMatrix[4][refPos]
	elif typ == 'gapRef':
		return scoreMatrix[5][refPos]
	else:
		return 0	

def computeFMatrix(ids,ref,seq2,scoreMatrix,indId,f,g):
	"""Computes the scoring matrix. ref : reference, seq2 : other sequence"""
	regularGap	=	0
	rows 		=	len(ref) + 1 
	cols		=	len(seq2) + 1
	#In tracematrix 0=diag, 1=left, 2=up
	traceMatrix	=	np.zeros(shape=(rows,cols))
	#Create the score matrix.
	fMatrix		=	np.zeros(shape=(rows,cols))
	for i in range(0,rows):
		fMatrix[i][0] = i * regularGap
	for j in range(0,cols):
		fMatrix[0][j] = j * regularGap
	for i in range(1, rows):
		refPos		=	i-1
		seqPos		=	1
		for j in seq2:
			if seqPos < cols-1:
				match	= fMatrix[i - 1][seqPos - 1] + checkScore(refPos,j,scoreMatrix)
				delete	= fMatrix[i - 1][seqPos] + checkScore(refPos,'gapSeq',scoreMatrix)
				insert	= fMatrix[i][seqPos - 1] + checkScore(refPos,'gapRef',scoreMatrix)
				if max(match, delete, insert) == match:
					fMatrix[i][seqPos]		= match
					traceMatrix[i][seqPos]	= 0
				elif max(match, delete, insert) == delete:
					fMatrix[i][seqPos]		= delete
					traceMatrix[i][seqPos]	= 2
				elif max(match, delete, insert) == insert:
					fMatrix[i][seqPos]		= insert
					traceMatrix[i][seqPos]	= 1
			seqPos	+= 1
	makeSequences(ids,ref,seq2,traceMatrix,indId,f,g)

def findIndels(refSeq):
	"""Returns a list with the indel positions of the reference sequence."""
	indels=[]
	for i in range(len(refSeq)):
		if refSeq[i]=='-':
			indels.append(i)
	return indels

def removeIndels(refSeq,seq):
	"""Removes the nucleotides in the indelpositions of seqDict"""
	indels 		=	findIndels(refSeq)
	for indel in indels:
		seq		=	seq[:indel]+'H'+seq[(indel+1):]
	seq		=seq.replace('H','')
	return seq

def makeSequences(ids,ref,seq2,traceMatrix,indId,f,g):
	"""returns the result"""
	refAl=''
	seq2Al=''
	i=len(ref)
	j=len(seq2)
	while i > 0 and j > 0:
		if traceMatrix[i][j] == 0:
			refAl = ref[i-1] + refAl
			seq2Al = seq2[j-1] + seq2Al
			i-=1
			j-=1
		elif traceMatrix[i][j] == 1:
			refAl = '-' + refAl
			seq2Al = seq2[j-1] + seq2Al
			j-=1
	 	elif traceMatrix[i][j] == 2:
			refAl = ref[i-1] + refAl
			seq2Al = '-' + seq2Al
			i-=1
	while i > 0:
		refAl = ref[i-1] + refAl
		seq2Al = '-' + seq2Al
		i-=1
	while j > 0:
		refAl = '-' + refAl
		seq2Al = seq2[j-1] + seq2Al
		j-=1
	print >> f, '>ref'
	print >> f, refAl
	print >> f, '>'+ids
	print >> f, seq2Al
	seq2Al	= removeIndels(refAl,seq2Al)
	print >> g, '>'+ids
	print >> g, seq2Al

def main():
	path_indata	=	'../data/raw_data/'	
	ref			=	'CGCCCGCTGCGCTCACCTCGCCGCTGCACCGTGAAGCTCTCAATCACCCCGTAGTTGTGTCTGCAGTAGGTGTCCACCGTTGCCCGCTCCTGCTCCAAGATCTCCTTCTGCCCGTTCCAGGACTCAGCGACGGGCCGCCCGAGCTCCGTGACCGCCCGGTACTCCCCCACGTCGCTGTCGAAGCGCACGAACTCCTCCCGGTTATGGATGTATCTTTCCACGAACCGCACCCGCTCCGTCCCGTTGGTGAAATAGCACTCGGACTTTGCCACCTCCAAGAAATGTGCTGTGGGGACGGGGGATC'
	scoreTable		= 'scoreTable.txt'
	scoreMatrix		= makeMatrix(scoreTable,len(ref))
	for infile in glob.glob(os.path.join(path_indata,'*.*')):
		if os.path.getsize(infile) > 0:
			indId	= getId(infile)
			print indId
			f	= open('../results/2012-03-03/sequences_aligned/'+indId+'.fa','a')
			g	= open('../results/2012-03-03/test/sequences_aligned_noref/' +indId+ '.fa', 'a')
			sequences 		= getSeq(infile)
			print '# of seq:',len(sequences)
			for ids, seq in sequences.items():
				computeFMatrix(ids,ref,seq,scoreMatrix,indId,f,g)
			f.close()
			g.close()
	print scoreMatrix
	


if __name__ == '__main__':
	main()

