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
out_file = os.path.join(out_dir, file_prefix + '.coord')
out_file_1 = os.path.join(out_dir, file_prefix + '.crd')
###############################################################################

# read big csv
df = pd.read_csv(inp_file, sep='\t', comment = '#', chunksize=1000, \
                low_memory=False, iterator = True, names=['chr', 'pro', 'feature', 'start', 'end', '.', 'strand', '2', 'info' ])
df = pd.concat(list(df), ignore_index=True)

df = df.loc[df['feature'] == 'CDS'] 

df['chr'] = df['chr'].str.strip('gnl|Elanco|')

df['new_start'] = df.apply(lambda row: row['end'] if row['strand'] == '-' else row['start'],  axis=1)
df['new_end'] = df.apply(lambda row: row['start'] if row['strand'] == '-' else row['end'],  axis=1)
df['new_start'] = df['new_start'].astype(int)
df['new_end'] = df['new_end'].astype(int)



if not os.path.exists(out_file): 
    df.to_csv(out_file, sep='\t', index=False, columns=['chr', 'new_start', 'new_end'], header=False)
    df.to_csv(out_file_1, sep='\t', index=False, columns=['chr', 'new_start', 'new_end', 'chr'], header=False)
    print('Completed writing report !!!')
else:
    print('Report exists !!!')
