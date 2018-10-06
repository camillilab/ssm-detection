#!/usr/bin/env python3

"""
Modifies an input sequence at a given nt position
For our purposes, inserts a sequence
"""

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO
import random
import argparse
import os

parser = argparse.ArgumentParser()

# required arguments
parser.add_argument('seq', help='Path to FASTA file')

# optional arguments
parser.add_argument('-p', '--position', help='Position to insert modified sequence. Default is random.', default=-1, type=int)
parser.add_argument('-s', '--seed', help='Seed sequence to insert', default='AGCT', type=str)
parser.add_argument('-c', '--copies', help='Number of copies of seed to insert', default=15, type=int)
parser.add_argument('-o', '--output', help='Output filename')

args = parser.parse_args()

# Read the reference sequence as FASTA - should be single sequence, so use SeqIO.read
record = SeqIO.read(args.seq, "fasta")
seq = record.seq

# generate the insert sequence
insert_seq = Seq(args.seed * args.copies)

# insert the sequence
loc = args.position
if loc == -1:
    loc = random.randint(0, len(seq))

new_seq = seq[0:loc-1] + insert_seq + seq[loc-1:len(seq)]
print("Inserted {0}x{1} at position {2}".format(args.seed, args.copies, loc))

# Wrap into a SeqRecord object
new_record = SeqRecord(new_seq, id=record.id + '_' + args.seed + 'x' + str(args.copies),
                       description=record.description + ' with {1} copies of {0} at {2}'.format(
                           args.seed, args.copies, loc))
f1 = SeqFeature(FeatureLocation(loc, loc + len(args.seed)*args.copies), type='domain')
new_record.features.append(f1)

# write the seq to file
base_name = args.seq.split('.')[0]
handle = os.path.join(os.getcwd(), args.output)
SeqIO.write(new_record, handle, "fasta")

print("Wrote inserted sequence to {0}".format(handle))
