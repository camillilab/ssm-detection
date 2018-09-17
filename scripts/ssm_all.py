#! usr/env/python

"""
Automation of ssm.py for a folder of BAM files
"""

import argparse
import glob
import os
import subprocess


# run smm.py, capturing output into log
def ssm(ref, mapping, log):
    with open(log, 'w') as f:
        subprocess.run(['python3', '-u', 'ssm.py', ref, mapping], stdout=f, stderr=subprocess.STDOUT, check=True)
    return


ssm_path = os.path.join(os.getcwd(), 'ssm_logs/')
os.mkdir(ssm_path)

parser = argparse.ArgumentParser()

parser.add_argument('bampath')
parser.add_argument('ref')

args = parser.parse_args()

bam_path = args.bampath
ref = args.ref

# get list of all bam files in the path
bam_files = [f for f in glob.glob(bam_path+'*.bam')]

i = 0
for bam_file in bam_files:

    print("Now processing {0}".format(bam_file))
    i += 1
    # define log file
    log_file = os.path.join(ssm_path, 'ssm_log_{0}'.format(i))

    # run ssm.py
    ssm(ref, bam_file, log_file)

print("Done!")
