#! usr/bin/env python

"""
Script for loading the output of ssm.py and annotating the genomic region in which the SSM event occurs, and
what the predicted SSM event does to the coding region if is in one.

Input: comma-separated gb files, ssm.py output
Output: annotated ssm.py

"""

import argparse
import os
import csv
from Bio import SeqIO

# ARGPARSE
parser = argparse.ArgumentParser()

# REQUIRED ARGUMENTS
parser.add_argument('ssm_output')
parser.add_argument('gbk_files', help='comma-separated list of gbk files to use for annotating. No spaces!')

args = parser.parse_args()
outfile = 'ssm_output_annotated.txt'

# load annotation data from gbk files into a dictionary
print(os.getcwd())
gbk_files = [os.path.join(os.getcwd(), k) for k in args.gbk_files.split(',')]
annotations = dict()
seqs = dict()

for gbk_file in gbk_files:
    gbk_cds = list()
    with open(gbk_file, 'r') as handle:
        record = SeqIO.read(handle, format='gb')
        elements = record.features
        accession_num = record.name
        seqs[accession_num] = record.seq
        for element in elements:
            if element.type == "CDS":
                gbk_cds.append(element)
    annotations[accession_num] = gbk_cds

# load the ssm.py output
with open(outfile, 'w') as write_handle:
    writer = csv.writer(write_handle, delimiter='\t')

    writer.writerow(('chr_name', 'rpos', 'gpos', 'repeat', 'repeat_length', 'indel code',
                     'ssm_reads', 'totalreads', 'percent SSM', 'locus_tag', 'descriptions', 'effect'))

    with open(args.ssm_output, 'r') as read_handle:
        reader = csv.DictReader(read_handle, delimiter='\t')

        for row in reader:
            name = row['chr_name']
            ssm_pos = int(row['gpos'])
            ssm_repeat = row['repeat']
            ssm_num_repeats = int(row['repeat_length'])
            ssm_code = row['indel code']

            # first, let's see where the repeat lies in the chromosome
            chr_annot = annotations[name]

            # get the range that the repeat stretches into
            repeat_start = ssm_pos
            repeat_length = len(ssm_repeat) * ssm_num_repeats
            repeat_end = ssm_pos + repeat_length

            # now look through the annotation data and try to see if someone peeks in there
            ssm_locus_tag = 'intergenic'
            ssm_description = ''
            ssm_effect = ''
            for cds in chr_annot:
                cds_loc_data = cds.location

                cds_start = cds_loc_data.start
                cds_end = cds_loc_data.end
                cds_strand = cds_loc_data.strand

                if (repeat_start >= cds_start) and (repeat_start <= cds_end):
                    print("CDS found!")
                    try:
                        ssm_locus_tag = cds.qualifiers['locus_tag']
                        ssm_description = cds.qualifiers['product']
                    except KeyError:
                        pass

                    # now that we know the gene, what effect does an insertion or deletion have on the gene?
                    # first, grab the sequence
                    cds_sequence = seqs[name][cds_start:cds_end]
                    cds_ssr_sequence = cds_sequence

                    # depending on the code, let's modify this sequence accordingly. Use the leftmost pos?
                    # TODO is this splicing correct?
                    if ssm_code == "I":
                        cds_ssr_sequence = cds_sequence[cds_start:repeat_start] + ssm_repeat + cds_sequence[repeat_start:cds_end]
                    if ssm_code == "D":
                        cds_ssr_sequence = cds_sequence[cds_start:repeat_start] + cds_sequence[repeat_start + len(ssm_repeat):cds_end]

                    # if we are on the negative strand, let's take the reverse complement to get the coding sequence
                    if cds_strand == -1:
                        cds_sequence = cds_sequence.reverse_complement()
                        cds_ssr_sequence = cds_ssr_sequence.reverse_complement()

                    # finally, let's get the translations
                    cds_translation = cds_sequence.translate()
                    cds_ssr_translation = cds_ssr_sequence.translate()
                    print("Woah")










