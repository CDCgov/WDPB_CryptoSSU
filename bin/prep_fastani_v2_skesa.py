#!/usr/bin/env python3

import os
import sys
import argparse
import numpy as np
import subprocess as sp
from Bio import SeqIO

def fastANIFormat(contigs,outDir):

	out_filename =  os.path.join(outDir,'fastANI_genomeList_skesa.txt')
	contigsDir = os.listdir(contigs)	
	with open(out_filename,"w+") as output_file:
		for file in contigsDir:
			if file.endswith('.fasta.gz'):
				output_file.write(os.path.abspath(file)+"\n")




def main():

	parser = argparse.ArgumentParser()
	parser.add_argument("-assemblies",required=True, help="REQUIRED: directory location of all assemblies")
	parser.add_argument("-output",required=True, help="REQUIRED: Output directory for fastANI query and reference genomes")


	args = parser.parse_args()
	contigs = args.assemblies
	outDir = args.output
	fastANIFormat(contigs,outDir)

if __name__ == '__main__':
	main()
	

