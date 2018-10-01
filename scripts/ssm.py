#! usr/env/python

"""
Detection of potential sites of slipped-strand mispairing (SSM) from sequencing sets
Takes in a sequencing scaffold, detects areas vulnerable to SSM, and then uses SAM/BAM data to
determine is SSM happens at a meaningful frequency
"""

import jakelib
from collections import defaultdict
import argparse
import os


# ARGPARSE
parser = argparse.ArgumentParser()

# REQUIRED ARGUMENTS
parser.add_argument('reference', help='Reference sequence. FASTA format. May be a series of contigs in a scaffold.')
parser.add_argument('mapping', help='BAM file of reads on reference sequence.')

# TODO Add these args in
"""
# OPTIONAL ARGUMENTS
parser.add_argument('-p', '--processes', help='Number of processes to use, should be equal to CPU cores', default=1)
parser.add_argument('-s', '--seed_length', help='Minimum/Maximum pattern seed lengths (min, max)', default=[1, 5],
                    nargs=2)
parser.add_argument('-r', '--threshold',
                    help='Repeat threshold. Must be this many iterations of a seed pattern to be a repeat [10]',
                    default=10)
"""
args = parser.parse_args()


# INPUT VARS
reference_file = args.reference
bam_file = args.mapping
#min_seed_repeat_length = args.seed_length[0]
#max_seed_repeat_length = args.seed_length[1]


# DEBUG
# reference_file = os.path.join(os.getcwd(), '.fa')
# bam_file = os.path.join(os.getcwd(), 'BAM.bam')
num_processes = 8
num_threads = 8
min_seed_repeat_length = 2
max_seed_repeat_length = 7
repeat_threshold = 5
coverage_range = 50  # looks this many bases around repeat region to determine average depth

# PREPROCESSING

# load reference sequence
reference = jakelib.SeqFile(reference_file, file_type='fasta')
print("Loaded {0} sequences of {1} nucleotides from {2}".format(len(reference.records),
                                                                reference.seq_size(),
                                                                reference.sequence_file))

# PATTERN DETECTION

# Generate repeat seed sequences...
print("Generating patterns seeds of length {0} to {1}...".format(min_seed_repeat_length, max_seed_repeat_length))
repeat_seeds = jakelib.generate_repeat_seeds(min_seed_repeat_length, max_seed_repeat_length, procs=num_processes)

# For each sequence, look for patterns
repeat_dict = dict()
name_pos_dict = defaultdict(list)

for record in reference.records:
    print("Detecting patterns in reference sequence {0}...".format(record.name))
    sequence = record.seq

    for pattern in repeat_seeds:
        matches = jakelib.scan_for_repeats(sequence, pattern, num_repeat_threshold=repeat_threshold)
        for pos, repeat in matches:
            # organize data into a dict - {(name, pos) : (repeat, length)}
            repeat_dict[(record.name, pos)] = (pattern, int(len(repeat) / len(pattern)))

            # keep a record of which repeats had what position
            name_pos_dict[record.name].append(pos)
print("{0} patterns discovered.".format(len(repeat_dict)))

# Save patterns to file
with open('patterns.txt', 'w') as f:
    f.write('chr_name\tgpos\tpattern\tnum\n')
    for key in repeat_dict:
        f.write('{0}\t{1}\t{2}\t{3}\n'.format(key[0], key[1], repeat_dict[key][0], repeat_dict[key][1]))


# READ ASSESSMENT

# load alignment file
print("Reading, sorting, and indexing alignment file {0}...".format(bam_file))
bam = jakelib.AlignmentFile(bam_file, chr_name='MAIN', use_threads=num_threads)
print("Done.")
# print("{0} reads.".format(bam.num_reads))

ssm_data = dict()
ssm_positive_nodes = list()

print("Finding SSM sites...")
for record in reference.records:

    name = record.name
    seq = record.seq

    outfile = os.path.join(os.getcwd(), name + '_sorted.bam')
    outfile_index = os.path.join(os.getcwd(), name + '_sorted.bam.bai')

    # gross hack for one ref TODO clean this up
    if len(reference.records) == 1:
        print("Only one record detected, renaming MAIN bam and index to {0}...".format(name))
        os.rename('MAIN_sorted.bam', outfile)
        os.rename('MAIN_sorted.bam.bai', outfile_index)
        new_bam = bam
        new_bam.give_name(name)
        new_bam.sorted_alignment_file = outfile
        new_dict = {k: repeat_dict[(name, k)] for k in name_pos_dict[name]}

    else:
        print("Processing {0}...".format(name))

        # Split the main alignment file by read group
        bam.dump(outfile, name=name)

        # Create a new AlignmentFile instance using the read subset, sort and index to grab reads around SSM sites
        # maybe have some skip or name assignment if only one sequence?
        new_bam = jakelib.AlignmentFile(outfile, chr_name=name, use_threads=num_threads)
        new_bam.give_name(name)

        # Grab a subset of the repeat dict involving this read group
        # Can't think of a better way to do this...maybe have a different dict with the repeat positions detected?
        new_dict = {k: repeat_dict[(name, k)] for k in name_pos_dict[name]}

    # Do the magic
    print("Launching {0} processes...".format(num_processes))
    data = jakelib.find_ssm(new_dict, new_bam, seq, procs=num_processes)
    print("{0} SSM site(s) found in {1}".format(len(data), name))

    if len(data) > 0:
        ssm_positive_nodes.append((name, seq))
    ssm_data.update(data)

print("{0} total SSM sites found over {1} sequence(s)".format(len(ssm_data), len(reference.records)))

# add contingency if 0 SSM reads were found?

# Disable temporarily as none of the artificial genomes are going to match....

# print("Detecting most likely reference sequence....")
# accession_num, query_shift_dict = jakelib.detect_reference_sequence_blast('blastn', 'ref_prok_rep_genomes', ssm_positive_nodes)

accession_num = 'NULL'

# now calculate coverage +- some bp around the repeat area
# probably easiest to report the average depth surrounding the area using samtools depth
# is it best to use samtools depth, and compare to the "read agreement"? unsure
# mpileup may have useful information related to phred qualities
# probably better to implement this in the above loop to avoid re-sorting bam files
print("Detecting coverage for SSM sites...")
coverage_data = dict()
for ssm in ssm_data:

    chr_name = ssm[0]
    rpos = ssm[1]
    gpos = rpos  # + query_shift_dict[chr_name]
    seed = ssm[2]
    num_repeats = ssm[3]
    code = ssm[4]
    read_depth = ssm_data[ssm]
    file = os.path.join(os.getcwd(), chr_name + '_sorted.bam')

    bam = jakelib.AlignmentFile(file, chr_name=chr_name, sort=False, index=False)
    average_depth = bam.average_depth(name=chr_name,
                                      start=rpos-coverage_range,
                                      end=rpos+(len(seed) * num_repeats)+coverage_range)

    # try to calculate reads in this region...need to adjust start/end based on read length
    # divide by two if paired end?
    # TODO Make read counting adaptive...samtools stats? Use GATK 25% rule?
    # TODO this runs a samtools stats for every SSM! Move out of the loop
    # read_length = bam.get_read_length()
    reads_in_region = len(bam.fetch(name=chr_name, start=rpos-100, end=rpos-10)) / 2

    coverage_data[chr_name, rpos, gpos, seed, num_repeats, code] = (read_depth, average_depth, reads_in_region)

# TODO add coverage threshold
print('chr_name\trpos\tgpos\trepeat\tlength\tcode\treadcov\tavgcov\ttotalreads')
for cov in coverage_data:
    chr_name = cov[0]
    rpos = cov[1]
    gpos = cov[2]
    seed = cov[3]
    num_repeats = cov[4]
    code = cov[5]
    read_depth = coverage_data[cov][0]
    average_depth = coverage_data[cov][1]
    reads_in_region = coverage_data[cov][2]
    print(chr_name, rpos, gpos, seed, num_repeats, code, read_depth, average_depth, reads_in_region, sep='\t')


# attempt to remove all files?
try:
    os.remove('MAIN_sorted.bam')
    os.remove('MAIN_sorted.bam.bai')
except FileNotFoundError:
    pass

try:
    for record in reference.records:
        name = record.name
        # clean bam files
        # TODO this doesn't work if there was only one record. Change behavior
        os.remove(os.path.join(os.getcwd(), name+'.bam'))
        os.remove(os.path.join(os.getcwd(), name + '_sorted.bam'))
        os.remove(os.path.join(os.getcwd(), name + '_sorted.bam.bai'))
except FileNotFoundError:
    pass

print("Done!")

# TODO Theoretically I could use some Illumina substitution error rate as a possible measure for background noise.






