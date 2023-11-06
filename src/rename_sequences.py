#!/usr/bin/env python
# 03 mai2023
# Usage rename_sequences.py (i)   map_file.tsv \
# 			                (ii)  sequences_to_rename \
#                           (iii) out_folder \
#                           (iv)  extension
# e.g. rename_sequences.py \
#           metadata/map_file2.tsv \
#           ./data/mini_reads \
#           ./data/mini_reads/renamed/ \
#           fastq.gz
#
#  map_file example #1  : we want to rename 2090-HCLON-0001 in CFI2-2-U18F-04
#  CFI2-2-U18F-04	2090-HCLON-0001
#  CFI2-3-U03G-05	2090-HCLON-0002
#  CFI2-3-U15G-03	2090-HCLON-0003

#  map file example #2  : we want to rename ugly sequence NS...66b_R1 in 66b_R1
#  66b_R1	NS.2066.002.IDT_i7_21---IDT_i5_21.RP-E1-66b_R1
#  66b_R2	NS.2066.002.IDT_i7_21---IDT_i5_21.RP-E1-66b_R2
#  67b_R1	NS.2066.002.IDT_i7_33---IDT_i5_33.RP-E1-67b_R1
#  67b_R2	NS.2066.002.IDT_i7_33---IDT_i5_33.RP-E1-67b_R2

import sys
import csv
import shutil
import os
from pathlib import Path

def get_files(start_dir, file_type):
    """ Given a directory, make a list of files
    This function will find all files recursively.
    """
    files = []
    for filename in Path(start_dir).glob('**/*.' + file_type):
        files.append(str(filename))
    return files

# Create a folder for renamed assemblies
cwd = os.getcwd()
if not os.path.exists(cwd + sys.argv[3]):
    os.mkdir(cwd + sys.argv[3])

# Create a dict where the key is 2nd column of the .tsv map file
D = {}
with open(sys.argv[1]) as f:
    rd = csv.reader(f, delimiter = "\t")
    next(rd, None)  # skip the headers
    
    for row in rd:
        if ".fastq.gz" in row[1]:
            row[1] = row[1].replace(".fastq.gz", "")
        D[row[1]] = row[0]

# Get a list of large assemblies fasta files
fsa_files = get_files (sys.argv[2], sys.argv[4])


# For all files copy (and rename) when a corresponding key is found
for f in fsa_files:
    fname = (os.path.basename(f).split('.' + sys.argv[4])[0])
    if fname in D:
        print(fname, 'has been renamed in', D[fname])
        #shutil.copyfile(f, sys.argv[3] + D[fname] + '.' + sys.argv[4])
        os.symlink(f, sys.argv[3] + D[fname] + '.' + sys.argv[4])
    else:
        print(fname, 'could not have been renamed, no matching key!')

