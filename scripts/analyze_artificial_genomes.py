#! usr/bin/env python

"""
Script for the alignment and SSM analysis of every artificial genome
Cause I'm too lazy to type each one in
Author: JB
Camilli Lab
09/09/2018
"""

import subprocess
import os

"""
Artificial genomes should have been outputted in the following format:

No extra insertions/deletions: art_ref_ssr1.fq, art_ref_ssr2.fq
100% insertion: art_ssr_ins1.fq, art_ssr_ins2.fq
100% deletion: art_ssr_del1.fq, art_ssr_del2.fq

And the various grades of insertion or deletion:
For example, in steps of 10%
art_ins_10_1.fq, art_ins_10_2.fq
art_ins_20_1.fq, art_ins_20_2.fq

And so on.
So I want to bowtie2 each of these to the reference genome art_ref_ssr.fasta
Convert the sam to bam using samtools
And then sun ssm.py to generate SSM reports
"""


# uses bowtie2-build to make a reference index
def bowtie2_build(ref, ind):
    subprocess.run(['bowtie2-build', ref, ind], stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    return


# align using bowtie2, capturing alignment output into log
# note bowtie2 uses stderr for output, oddly enough
def bowtie2(ind, fq1, fq2, sam, log):
    subprocess.run(['bowtie2', '-x', ind, '-1', fq1, '-2', fq2, '-S', sam, '-p', '8', '2>', log], check=True)
    return


# use samtools to convert sam to bam
def bamify(sam, bam):
    subprocess.run(['samtools', 'view', '-b', sam, '-o', bam, '-@', '8'], check=True)
    os.remove(sam)  # delete old sam file
    return


# run smm.py, capturing output into log
def ssm(ref, mapping, log):
    with open(log, 'w') as log_file:
        subprocess.run(['python3', '-u', 'ssm.py', ref, mapping], stdout=log_file, stderr=subprocess.STDOUT, check=True)
    return


# main script
def main():

    # alright, define paths and files
    print("Making directories...")

    ref = '/Volumes/camillilab/ART/art_refs/art_ref_ssr.fasta'
    control_fqs = ['/Volumes/camillilab/ART/fq/art_ref_ssr1.fq', '/Volumes/camillilab/ART/fq/art_ref_ssr2.fq']
    bt2_index_name = 'SSM'

    # path to fqs
    ssm_insertion_fq_path_base = '/Volumes/camillilab/ART/fq'
    ssm_deletion_fq_path_base = 'Volumes/camillilab/ART/fq'

    # path to save bowtie2 map logs
    bt2_log_path_base = os.path.join(os.getcwd(), 'bt2_log')
    if not os.path.exists(bt2_log_path_base):
        os.mkdir(bt2_log_path_base)

    # path to save ssm logs
    ssm_log_path_base = os.path.join(os.getcwd(), 'ssm_log')
    if not os.path.exists(ssm_log_path_base):
        os.mkdir(ssm_log_path_base)

    # path to save bam files
    bam_path_base = os.path.join(os.getcwd(), 'bam')
    if not os.path.exists(bam_path_base):
        os.mkdir(bam_path_base)

    # now for the simulation!
    print("Generating bowtie2 index...")
    bowtie2_build(ref, bt2_index_name)

    print("Now processing control reference genome (the patten exists, but should be no SSM detected)...")
    bt2_ref_log = os.path.join(bt2_log_path_base, 'ctrl.log')

    print("Aligning using bowtie2...")
    bowtie2(bt2_index_name, control_fqs[0], control_fqs[1], 'ctrl.sam', bt2_ref_log)

    print("Converting to BAM...")
    bamify('ctrl.sam', os.path.join(bam_path_base, 'art_ref_ssr.bam'))

    print("Mapping using SSM...")
    ssm(ref, os.path.join(bam_path_base, 'art_ref_ssr.bam'), os.path.join(ssm_log_path_base, 'art_ref_ssr.log'))

    # Now lets loop this for errythang
    for i in range(10, 110, 10):
        for j in ['ins', 'del']:

            print("Processing type {1} at percent {0}%...".format(j, i))

            # define fq files
            fq1_name = 'art_{1}_{0}_1.fq'.format(i, j)
            fq2_name = 'art_{1}_{0}_2.fq'.format(i, j)
            fq_path_base = '/Volumes/camillilab/ART/fq'

            fq1 = os.path.join(fq_path_base, fq1_name)
            fq2 = os.path.join(fq_path_base, fq2_name)


            # define output files
            bt2_log_name = 'bt2_{0}_{1}.log'.format(j, i)
            bt2_log = os.path.join(bt2_log_path_base, bt2_log_name)

            ssm_log_name = 'ssm_{0}_{1}.log'.format(j, i)
            ssm_log = os.path.join(ssm_log_path_base, ssm_log_name)

            bam_name = 'art_{1}_{0}.bam'.format(i, j)
            bam_file = os.path.join(bam_path_base, bam_name)

            # bowtie2
            bowtie2(bt2_index_name, fq1, fq2, 'temp.sam', bt2_log)

            # samtools
            bamify('temp.sam', bam_file)

            # ssm
            ssm(ref, bam_file, ssm_log)

    print("Done!")


if __name__ == '__main__':
    main()
