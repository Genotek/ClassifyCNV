#!/usr/bin/env python3
# coding=utf-8

import argparse
import sys

"""
The script parses a ClinGen tsv file that contains dosage sensitive regions
The output is two bed files, one for haploinsufficient regions and one for triplosensitive regions
"""
parser = argparse.ArgumentParser()
parser.add_argument('--infile', help='Input file in tsv format, downloaded from ClinGen;')
args = parser.parse_args()

if args.infile.endswith('tsv'):
    outfile_hi = args.infile.replace('tsv', 'HI.bed')
    outfile_ts = args.infile.replace('tsv', 'TS.bed')
else:
    sys.exit('The file is not in the tsv format\n')

infile = open(args.infile, 'r')
hi_out = open(outfile_hi, 'w')
ts_out = open(outfile_ts, 'w')

# go through the file and parse positions, haploinsufficiency and triplosensitivity scores
for line in infile:
    if not line.startswith('#'):
        # parse out chr, start and end
        fields = line.strip().split('\t')
        if fields[3].startswith('chr'):
            (chr, pos) = fields[3].strip().split(':')
            (start, end) = pos.strip().split('-')
            # if the score field is not empty, print to file
            if fields[4]:
                hi_out.write('{}\t{}\t{}\t{}\n'.format(chr, start.strip(), end.strip(), fields[4]))
            if fields[12]:
                ts_out.write('{}\t{}\t{}\t{}\n'.format(chr, start.strip(), end.strip(), fields[12]))

infile.close()
hi_out.close()
ts_out.close()
