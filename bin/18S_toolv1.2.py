#!/usr/bin/env python

import argparse
import os
import subprocess
import pandas as pd
import shutil
import numpy as np


#Tool details
'''
This script takes the raw blast results to process and filter the best hits for 18s RNA classification of Cryptosporidium
usage: python <script> --localdir <rawblastresultsfolder> --resultsdir <finalresultfolder>
'''
algoname = "Crypto_18S_RNA_typing"
version=1.1

#Set colors for cmdline messages/warnings
CRED = '\033[91m' + '\nError:'
CYEL = '\033[93m' + '\nWarning:'
CGRE = '\033[92m' + '\nStatus:'
CEND = '\033[0m'

class myargumentparser(argparse.ArgumentParser):
	# override to split space based cmdline args
	def convert_arg_line_to_args(self, arg_line):
		return arg_line.split()

def cmdline_args():
	parser = myargumentparser(fromfile_prefix_chars='@')
	#parser.add_argument('--reference_folder',default='',help="Blast Reference Database",type=str) # where reference database is located
	#parser.add_argument('--query',default='', help='filename for the assembled query genome', type=str) # where fasta assemblies are located
	parser.add_argument('--inputdir',default='',help='blast files directory',type=str) #takes in the blast input files 
	parser.add_argument('--resultsdir',default='',help='results directory',type=str) #this is where intermediate files and results are saved
	parser.add_argument('--localdir', help='scratch location for intermediate files', type=str) # just a place holder
	args = parser.parse_args()
	return args

def dirs(args):
	folder= (args.resultsdir)
	for f in folder:
		os.path.join(os.getcwd(),f)
		os.makedirs(f, exist_ok=True)

def subdirs(args):
	mode=0o777
	tempsub = args.localdir+"/sorted_blastresults"
	os.makedirs(tempsub,mode,exist_ok=True)
	# os.makedirs(logsub,mode,exist_ok=True)

def blast_output(args):
	fpath=os.path.join(args.localdir, "sorted_blastresults")
	for blast in os.listdir(args.inputdir):
		if blast.endswith(".blast"):
			genename=os.path.basename(blast)
			genomename=genename.split(".")[0]
			blastresult_path=os.path.join(args.inputdir, blast)
			gene={}			
			blastresultsDF = pd.read_csv(blastresult_path,sep="\t", header=None)			
			blastresultsDF.columns=['qseqid','sseqid','pident','length','qcovs','mismatch','gapopen','qstart','qend','sstart','send','qlen','slen','bitscore','qseq']
			# print(blastresultsDF)
			with open(blastresult_path) as blastresult:
				fname=os.path.join(fpath,"{}.blast".format(genomename))
				for line in blastresult:
					# print(line)
					try:
						line = line.split( )
						qseqid=line[0]
						sseqid=line[1]
						pident=float(line[2])
						length=(line[3])
						qcovs=float(line[4])
						mismatch=(line[5])
						gapopen=(line[6])
						qstart=int(line[7])
						qend=int(line[8])
						sstart=int(line[9])
						send=int(line[10])
						qlen=int(line[11])
						slen=int(line[12])
						# qseq=(line[13])
						bitscore=(line[13])
						qseq=(line[14])
						# sec_qlen=len(qseq)
						# print(qseqid)
					# print(sec_qlen)
						if qstart < qend:
							aligned_query_length = max(qend - qstart + 1, 0)  # Length of aligned region on query
						else:
							aligned_query_length = max(abs(qstart - qend)+ 1, 0) #indicates the alignment is in reverse direction
						aligned_subject_length = max(send - sstart + 1, 0)
						if aligned_subject_length != 0:
							# print(aligned_subject_length)
							query_coverage = (aligned_query_length / aligned_subject_length) * 100
						else:
							if aligned_query_length == 0:
								query_coverage = 0  # If both are zero, coverage is zero
							else:
								query_coverage = 100  # Assume full coverage if aligned_subject_length is zero but aligned_query_length is not
						if (pident> 98) & (query_coverage > 75)  :
							if qseqid in gene:
								gene[qseqid].append(sseqid)
							else:
								gene[qseqid]=[sseqid]
						#print(gene)
					except IOError as err:
						print(err)
			with open(fname,"a") as ofile:
				for key in gene:
					if len(gene[key]) > 1:
						for i in range(0,len(gene[key])):
							# print(i)
							tempResult = blastresultsDF.loc[blastresultsDF['sseqid']==gene[key][i]].values[0]
							finalResults = tempResult.tolist()
							str1 = '\t'.join(str(e) for e in finalResults)
							# print(str1)
							ofile.write(str1+"\n")
					else:
						#print(gene[key][0])
						tempResult = blastresultsDF.loc[blastresultsDF['sseqid']==gene[key][0]].values[0]
						str1 = '\t'.join(str(e) for e in tempResult)
						ofile.write(str1+"\n")
					
				ofile.close				

def write_csvs(args):
	blastpath=os.path.join(args.localdir, "sorted_blastresults")
	csv_dir = os.path.join(blastpath, "blast_csv")
	os.makedirs(csv_dir, exist_ok=True)
	for file in os.listdir(blastpath):
		if file.endswith(".blast"):
			filename=file.split(".")[0]
			blast_file_path = os.path.join(blastpath, file)
			if os.path.getsize(blast_file_path) == 0:
				print(f"Skipping empty file: {file}")
				continue
			try:
				blast_data = pd.read_csv(blast_file_path, sep="\t", header=None)
			except pd.errors.EmptyDataError:
				print(f"File is empty or not found: {file}")
				continue
			try:
				blast_data.columns=["query_genome","db_bestmatch","pident","alignment_length","coverage","mismatch","gapopen","query_start","query_end","subject_start","subject_end","query_length","subject_length","bitscore","qseq"]
				data_filtered=blast_data[blast_data.bitscore == blast_data.bitscore.max()]
				filtered_uniq = data_filtered.drop_duplicates()
				blast_csvfmt=os.path.join(csv_dir, f"{filename}.csv")
				filtered_uniq.to_csv(blast_csvfmt)
			except Exception as e:
				print(f"Error processing file {file}: {e}")

def filter_besthit(args):
	i = 0
	csvpath=os.path.join(args.localdir,"sorted_blastresults","blast_csv")
	csvfiles=[file for file in os.listdir(csvpath) if file.endswith(".csv")]
	if not csvfiles:
		df_empty = pd.DataFrame({"NCE": ["Sample did not meet the BLAST parameters, check input sequence or BLAST results"]})
		df_empty.to_csv(os.path.join(args.resultsdir, "NO_18SResults.csv"), index=False)
	else:
		for file in csvfiles:
			filename=os.path.basename(file)
			fname=filename.split(".")[0]
			data=pd.read_csv(os.path.join(csvpath, file), sep=",", index_col=False)
			data = data.drop(["Unnamed: 0","mismatch","gapopen","query_length","subject_length"],axis=1)
			data["rank"]=data[["query_genome","pident","coverage","bitscore"]].apply(tuple,axis=1).rank(method='dense',ascending=False).astype(int)
			#print("HERE")
			#print(data["db_bestmatch"][0])
			#species=data["db_bestmatch"][0].split('_',2, expand=True)
			species=data["db_bestmatch"][0].split('_',2)
			data["species"] = species[0]+"."+species[1]
			data["sample_name"]=fname
			data_uniq = data.drop_duplicates()
			NCE_1 = ["Alignment length is <700bp - do manual check for Identity and coverage" if x < 700 else "NA" for x in data_uniq["alignment_length"]]
			NCE_2 = ["Sample has multiple hits - Review manually" if len(data_uniq) > 1 else "NA"] * len(data_uniq)
			NCE_combined = [f"{n1};{n2}" for n1, n2 in zip(NCE_1, NCE_2)]
			data_uniq["NCE"]= NCE_combined
			finaldata = data_uniq.drop(columns=["rank"])
			# print(finaldata)
			# Reorder the columns
			column_order = ["sample_name","query_genome", "db_bestmatch", "species", "pident", "alignment_length", "coverage", "query_start", "query_end", "subject_start", "subject_end", "bitscore", "NCE"]
			finaldata = finaldata[column_order]
			output_file = os.path.join(args.resultsdir, f"{fname}.18SResults.csv")
			if i == 0:
				finaldata.to_csv(output_file,index=False)
			else:
				finaldata.to_csv(output_file,mode='a',index=False, header = False)
			i += 1

def main():
	args = cmdline_args()
	dirs(args)
	subdirs(args)
	blast_output(args)
	write_csvs(args)	
	filter_besthit(args)
	#print(CGRE + "Job has completed successfully" + CEND)

if __name__=='__main__':
	main()
