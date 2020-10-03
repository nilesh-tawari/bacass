#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Sat Oct  3 09:47:24 2020

@author: nilesh-tawari <tawari.nilesh@gmail.com>
"""
import os
import os.path
import ntpath
import argparse


'''
USAGE: python 

perl -lane '$c=$1 if /^LOCUS\s+(\S+)/;s/^ACCESSION\s+/ACCESSION   $c/; s/^VERSION/VERSION     $c\.1/;s/Elanco\://g;s/^DEFINITION  (\W+)/$1/;if (/(DEFINITION\s+\S+ \S+)/){$_=$1};s/Unclassified./Bacteria; Actinobacteria; Streptomycetales; Streptomycetaceae;\n            Streptomyces./;print'  /home/NILESH_RAMESH.TAWARI/Projects/17_bacterial_anno/rereun_samples/iseq_samples/03_annotation/results/105/105_prokka_rrna_rfam_nfast/BSUB105.gbk > /home/NILESH_RAMESH.TAWARI/Projects/17_bacterial_anno/rereun_samples/iseq_samples/03_annotation/results/105/105_prokka_rrna_rfam_nfast/BSUB105.selector.final.gbk ; done


'''

# Parameters
parser = argparse.ArgumentParser(description = 'Convert Prokka GBK to selector compatible format')
parser.add_argument('-i', '--inp_file', help='Input file')


args = parser.parse_args()
inp_file = os.path.abspath(args.inp_file)



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
out_file = os.path.join(out_dir, file_prefix + '.selector.gbk')
###############################################################################

out = open(out_file, 'w')
# read gbk file
with open(inp_file, 'r') as input_gbk:
    for line in input_gbk:
        if line.startswith('LOCUS'):
            tag = line.split()[1]
        if line.startswith('DEFINITION'):
            line = line.split()[0] + '  ' + line.split()[1] + ' ' + line.split()[2] + '\n'
        if line.startswith('ACCESSION'):
            line = line.strip() + '   ' + tag + '\n'
        if line.startswith('VERSION'):    
            line = line.strip() + '     ' + tag + '.1' + '\n'
        if line.startswith('            Unclassified.'):    
            line = '            Bacteria; Actinobacteria; Streptomycetales; Streptomycetaceae;' + '\n' + '            Streptomyces.\n'
        out.writelines(line)
        

out.close()
