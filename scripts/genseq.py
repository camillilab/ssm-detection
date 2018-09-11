#!/usr/bin/env python

"""
Generates a random nucleotide sequence in FASTA format

Arguments: length, GC content

Requirements: biopython
"""

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from Bio import SeqIO
from Bio.SeqUtils import GC
import random
import argparse
import os

parser = argparse.ArgumentParser()

# required arguments
parser.add_argument('length', type=int, help="Length of generated sequence")

# optional arguments
parser.add_argument('-c', '--gc', type=float, help='GC content. Default = 0.5. Keep to three decimals', default=0.5)
parser.add_argument('-o', '--output', help='Output filename')

args = parser.parse_args()

# Generate a random nucleotide sequence of desired length and GC content and wrap string in BioSeq object
seq = ''
seed = 'GC'*int((1000*args.gc)) + 'AT'*int((1000*(1-args.gc)))
for i in range(0, args.length):
    seq += random.choice(seed)
seq = Seq(seq, alphabet=generic_dna)

# Wrap into a SeqRecord object
record = SeqRecord(seq,
                   id='XX_ART_GC'+str(100*args.gc)+'_LEN'+str(args.length),
                   description='Artificial sequence of length {0} and GC% {1:.2f}'.format(args.length, GC(seq)))

print('Sequence of length {0} and GC% {1:.2f} generated.'.format(args.length, GC(seq)))

# Save to output file
handle = os.path.join(os.getcwd(), args.output)
SeqIO.write(record, handle, "fasta")
print("Wrote sequence to {0}".format(handle))
