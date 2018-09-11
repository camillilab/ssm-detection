#! /usr/bin/env python

"""
Combines a series of fasta files into one fasta record

Usage: merge_fasta.py <files, comma-separated> -o <filename.fa>
"""

from Bio import SeqIO
import argparse, os

parser = argparse.ArgumentParser()

# required arguments
parser.add_argument('files', type=str, help="Comma-separated list of fasta files to combine, no spaces")

# optional arguments
parser.add_argument('-o', '--output', help='Output filename')

args = parser.parse_args()

# read files
files = args.files.split(',')
records = []

for file in files:
    with open(file, 'r') as handle:
        for record in SeqIO.parse(handle, "fasta"):
            records.append(record)

# write out as single record
handle = os.path.join(os.getcwd(), args.output)
SeqIO.write(records, handle, "fasta")
print("Merged {0} into single fasta record {1}".format(', '.join(files), handle))

