#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 20 11:40:15 2020

@author: nilesh-tawari <tawari.nilesh@gmail.com>
"""
import os
import os.path
import ntpath
import argparse


'''
USAGE: python -i BVEN_C01.0.AOI_01.faa -p BSUB -s 105
'''

# Parameters
parser = argparse.ArgumentParser(description = 'Rename contigs in fasta file')
parser.add_argument('-i', '--inp_file', help='Input file')
parser.add_argument('-p', '--prefix', help='prefix for contigs e.g. BSUB')
parser.add_argument('-s', '--strain', required=False, help='strain for file name e.g. 64') 

args = parser.parse_args()
inp_file = os.path.abspath(args.inp_file)
prefix = args.prefix
strain = args.strain

###############################################################################
def path_leaf(path):
    """
    split the path in tail and head
    str -> str, str
    """
    head, tail = ntpath.split(path)
    return head, tail or ntpath.basename(head)

out_dir, filename_in = path_leaf(inp_file)
file_prefix = os.path.splitext(filename_in)[0]
file_suffix = os.path.splitext(filename_in)[1]
filename_out = file_prefix + "_renamed" + file_suffix
out_file = os.path.join(out_dir, filename_out)
#if strain != 'a':
#   out_file = os.path.join(out_dir, strain + '_renamed_' + filename_in)

###############################################################################

i = 1
with open(out_file, 'w') as out_fasta:
    with open(inp_file, 'r') as fasta:
        for line in fasta:
            if line.startswith('>'):
                line = '>{0}_C{1:02d} {2}'.format(prefix, i, line[1:])
                i = i + 1
            out_fasta.writelines(line)
