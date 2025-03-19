#!/usr/bin/env python3

#dev:Anusha Ginni,qxu0@cdc.gov, 
#edit: Shatavia Morrison
#1/08/2025
#3/19/2025
#Func:scripts extracts the required elements from the quast report, edit as per your needs
#USAGE: python quast_extract.py <path to report.tsv from quast module> <extractreport.tsv>

# "Assembly": "Sample Name",
#                 "# contigs": "Number of Contigs_Skesa",
#                 "Total length": "Total Length_Skesa",
#                 "L50": "L50_Skesa",
#                 "N50": "N50_Skesa",
#                 "GC (%)":"GC% Content_Skesa",
#                 "# N's per 100 kbp":"# N's per 100 kbp_Skesa"

import csv
import pandas as pd
import os
import argparse

def parse_args(args=None):
    Description = "scripts extracts the required elements from the quast report, edit as per your needs."
    Epilog = "Example usage: python quast_extract.py -infile <path to report.tsv from quast module> -outfile <extractreport.tsv"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("inreport", help="Input quast report file.")
    parser.add_argument("outreport", help="Output extracted quast report file.")

    return parser.parse_args(args)


def extract_quast_report(infile, outfile):
    columns_to_extract = {
    "# contigs": "Number of Contigs_Skesa",
    "Total length": "Total Length_Skesa",
    "L50": "L50_Skesa",
    "N50": "N50_Skesa",
    "GC (%)": "GC Content_Skesa",
    "# N's per 100 kbp": "# N's per 100 kbp_Skesa"
    }
    df = pd.read_csv(infile, sep="\t", index_col=0)
    extracted_data = df.loc[columns_to_extract.keys()].transpose()
    extracted_data.rename(columns=columns_to_extract, inplace=True)
    extracted_data.index.name = "Sample Name"
    extracted_data.index = extracted_data.index.str.split('.').str[0]
    extracted_data.reset_index(inplace=True)
    extracted_data = extracted_data.astype(str)
    extracted_data.to_csv(outfile, index=False)

def main(args=None):
    args = parse_args(args)
    extract_quast_report(args.inreport, args.outreport)


if __name__ == "__main__":
    main()
