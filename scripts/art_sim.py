#! /usr/bin/env python3

"""
Python script for the successive automation of ART-MountRainier NGS read simulation
Had to switch to ART because GemSim is hella slow...can't beat C++ in speed I guess
Metagenomics is simulated by combining fastq files at the end at define proportions
"""

import argparse
import os
import subprocess
from numpy import arange
import shutil


# construct and subprocess art commands
def art_paired(prof, r, nr, rl, ml, s, out):
    # create the command
    subprocess.run(['./art_illumina', '-ss', prof, '-i', r, '-p', '-na', '-c', nr, '-l', rl, '-m', ml, '-s',
                    s, '-o', out], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, encoding='utf-8')
    return


parser = argparse.ArgumentParser()

# optional arguments
parser.add_argument('-p', '--profile', help='Simulation profile. ART can pick one by default')
parser.add_argument('-n', '--numreads', help='Number of total reads to simulate for any given abundance')
parser.add_argument('-R', '--refpath', help="Path to reference genomes")
parser.add_argument('-pe', help='Switch for paired-end simulation', action='store_true')
parser.add_argument('-s', help="Stddev for fragment length")
parser.add_argument('-l', '--readlength', help="Length of read to simulate")
parser.add_argument('-i', '--insertsize', help="Mean insert length. Default is 3x read length", default=-1)
parser.add_argument('-a', '--start', help="Starting percent of reads", type=int)
parser.add_argument('-b', '--end', help="End percent of reads", type=int)
parser.add_argument('-c', '--step', help="Step percent size", type=int)

# parser.add_argument('-o', help='Output fastq file prefix')


args = parser.parse_args()

# We want to simulate a few controls - all normal, all insert, all deleted, and no insert - as a baseline
# Then, we want to successively simulate various abundance profiles by combining reads
# To keep the analysis simple, let's just test various levels of insertion OR deletion relative to control
# At first, let's just test 25%, 50%, and 75% mixes
step_size = args.step / 100
start = args.start / 100
end = args.end / 100

# Control simulation

# 2: 100% standard insert genome
print("simulating insert genome...")
ref = os.path.join(args.refpath, 'art_ref_ssr.fasta')
art_paired(prof=args.profile,
           r=ref,
           nr=args.numreads,
           rl=args.readlength,
           ml=args.insertsize,
           s=args.s,
           out='art_ref_ssr')


# Now create various combined fastq sequencing files based off step size

# insertion profiles
for i in arange(start, end + step_size, step_size):
    i = round(i, 2)  # fix i for floating errors
    print("Now doing insert {0}".format(i))
    ref_reads = str(int(int(args.numreads) * (1 - i)))
    ins_reads = str(int(int(args.numreads) * i))
    ref = os.path.join(args.refpath, 'art_ref_ssr.fasta')
    ins_ref = os.path.join(args.refpath, 'art_ref_ssr_ins.fasta')

    art_paired(prof=args.profile,
               r=ref,
               nr=ref_reads,
               rl=args.readlength,
               ml=args.insertsize,
               s=args.s,
               out='art_refi')
    art_paired(prof=args.profile,
               r=ins_ref,
               nr=ins_reads,
               rl=args.readlength,
               ml=args.insertsize,
               s=args.s,
               out='art_insi')
    # concatenate files using shutil since these are likely large files
    handle = os.path.join(os.getcwd(), 'art_ins_{0}_1.fq'.format(int(i*100)))
    with open(handle, 'wb') as wfd:
        for f in ['art_refi1.fq', 'art_insi1.fq']:
            with open(f, 'rb') as fd:
                shutil.copyfileobj(fd, wfd)
    handle = os.path.join(os.getcwd(), 'art_ins_{0}_2.fq'.format(int(i * 100)))
    with open(handle, 'wb') as wfd:
        for f in ['art_refi2.fq', 'art_insi2.fq']:
            with open(f, 'rb') as fd:
                shutil.copyfileobj(fd, wfd)
    # remove temp files
    os.remove('art_refi1.fq')
    os.remove('art_insi1.fq')
    os.remove('art_refi2.fq')
    os.remove('art_insi2.fq')

# deletion profiles
for i in arange(start, end + step_size, step_size):
    i = round(i, 2)
    ref_reads = str(int(int(args.numreads) * (1 - i)))
    del_reads = str(int(int(args.numreads) * i))
    print("Now doing del {0}".format(i))
    ref = os.path.join(args.refpath, 'art_ref_ssr.fasta')
    del_ref = os.path.join(args.refpath, 'art_ref_ssr_del.fasta')

    art_paired(prof=args.profile,
               r=ref,
               nr=ref_reads,
               rl=args.readlength,
               ml=args.insertsize,
               s=args.s,
               out='art_refd')
    art_paired(prof=args.profile,
               r=del_ref,
               nr=del_reads,
               rl=args.readlength,
               ml=args.insertsize,
               s=args.s,
               out='art_deld')
    # concatenate files using shutil since these are likely large files
    handle = os.path.join(os.getcwd(), 'art_del_{0}_1.fq'.format(int(i*100)))
    with open(handle, 'wb') as wfd:
        for f in ['art_refd1.fq', 'art_deld1.fq']:
            with open(f, 'rb') as fd:
                shutil.copyfileobj(fd, wfd)
    handle = os.path.join(os.getcwd(), 'art_del_{0}_2.fq'.format(int(i * 100)))
    with open(handle, 'wb') as wfd:
        for f in ['art_refd2.fq', 'art_deld2.fq']:
            with open(f, 'rb') as fd:
                shutil.copyfileobj(fd, wfd)
    # remove temp files
    os.remove('art_refd1.fq')
    os.remove('art_deld1.fq')
    os.remove('art_refd2.fq')
    os.remove('art_deld2.fq')


# Then, bowtie2 each of these back to the inserted reference genome (#2)
# Build the index for the initial genome and ref inserted genome


# Finally, BAMify the SAM files

