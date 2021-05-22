#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 11 16:04:35 2021

@author: nilesh-tawari <tawari.nilesh@gmail.com>
"""
import os
import os.path
import ntpath
import argparse
import pandas as pd

'''
USAGE: python 
'''

# Parameters
parser = argparse.ArgumentParser(description = ' ')
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
out_file = os.path.join(out_dir,  'misc_feature.tab.txt')

###############################################################################

# read big csv
df = pd.read_csv(inp_file, sep='\t', comment = '#', chunksize=1000, \
                low_memory=False, iterator = True, names=['chromosome', 'pro', 'feature', 'start', 'end', '.', 'strand', '2', 'locus_tag' ])
df = pd.concat(list(df), ignore_index=True)

df = df.loc[df['feature'] == 'prophage_region'] 
df['FeatureType'] = 'misc_feature'
df['strand'] = 1
df['inference'] = 'COORDINATES:profile:PhiSpy:3.7.8'
df['note'] = 'prophage'
df['locus_tag'] = df['locus_tag'].str.split("=", expand=True)[1]


if not os.path.exists(out_file): 
    df.to_csv(out_file, sep='\t', index=False, columns=['chromosome', 'FeatureType', 'start', 'end', 'strand', 'inference', 'note', 'locus_tag'])
    print('Completed writing report !!!')
else:
    print('Report exists !!!')
