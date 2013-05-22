#!/usr/bin/env python
# encoding: utf-8
"""
plot_het_hom.py

Plot how the relations between hetero and homozygotes varies with treshold for variable positions.

Created by Måns Magnusson on 2012-06-14.
Copyright (c) 2012 __MyCompanyName__. All rights reserved.
"""

import sys
import os
import numpy as np
import matplotlib.pyplot as plt


def main():
	N		= 10
	hetero	= (2698,2612,2583,2569,2551,2513,2454,2319,1924,0)
	homo	= (758,844,873,887,905,943,1002,1137,1532,3456)
	ind		= np.arange(N)
	width	= 0.35
	p1		= plt.bar(ind, hetero, width, color='r')
	p2 		= plt.bar(ind, homo, width, color='b', bottom=hetero)
	
	plt.ylabel('Individuals')
	plt.xlabel('Treshold')
	plt.title('Relation Hetero- Homozygotes')
	plt.xticks(ind+width/2., ('0.05', '0.10', '0.15', '0.20', '0.25', '0.30', '0.35', '0.40', '0.45', '0.50') )
	plt.yticks(np.arange(0,3500,200))
	plt.legend( (p1[0], p2[0]), ('Hetero', 'Homo'))

	plt.savefig('../doc/pictures/het_hom.png')


if __name__ == '__main__':
	main()

