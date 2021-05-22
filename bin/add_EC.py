#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 17 14:20:55 2020

@author: nilesh-tawari <tawari.nilesh@gmail.com>
"""
import os
import os.path
import ntpath
import argparse
import pandas as pd

'''
USAGE: python add_EC.py -e sequenceECs.txt -c cds.tab
'''

# Parameters
parser = argparse.ArgumentParser(description = 'Add EC numbers to CDS tab file')
parser.add_argument('-e', '--ec_file', help='EC file')
parser.add_argument('-c', '--cds_file', help='CDS tab file')

args = parser.parse_args()
inp_file = os.path.abspath(args.ec_file)
cds_file = os.path.abspath(args.cds_file)

###############################################################################
def path_leaf(path):
    """
    split the path in tail and head
    str -> str, str
    """
    head, tail = ntpath.split(path)
    return head, tail or ntpath.basename(head)

out_dir, filename_in = path_leaf(cds_file)
file_prefix = os.path.splitext(filename_in)[0] 
out_file = os.path.join(out_dir, file_prefix + '.ecFixed.txt')
###############################################################################

i = 0
ec_dir = {}
# read big csv
with open(inp_file, 'r') as seq_ec:
    for line in seq_ec:
        if line.startswith('>'):
            pid = line.split()[0].strip()[1:]
            ec_dir[pid] = [[], []]
        if line.startswith('#'):
            other_ec = line.split()[0].strip()
            if other_ec not in ec_dir[pid][1]:
                ec_dir[pid][1].append(other_ec)
        if not line.startswith('#') and not line.startswith('>') and line != '\n':
            ec = line.split()[0].strip()
            if ec not in ec_dir[pid][0]:
                ec_dir[pid][0].append(ec)

df_ec = pd.DataFrame.from_dict(ec_dir,orient='index')
df_ec.columns = ['EC', 'OtherEC']
#df1 = pd.concat({k: pd.Series(v) for k, v in oderded_dir.items()})

df_ec['EC'] = df_ec['EC'].apply(lambda x: ', '.join(x))
df_ec['OtherEC'] = df_ec['OtherEC'].apply(lambda x: ', '.join(x))
df_ec = pd.DataFrame(df_ec, columns=['EC'])
df_ec.replace(r'', '-', inplace = True)
df_ec.fillna(value='-', inplace = True)
#df = pd.DataFrame(list(ec_dir.items()))

# read CDS
df = pd.read_csv(cds_file, sep='\t', comment = '#', chunksize=1000, \
                low_memory=False, iterator = True)
df = pd.concat(list(df), ignore_index=True)

df_merge = pd.merge(df, df_ec, left_on='locus_tag', right_index=True, how='left')
df_merge['EC'].fillna(value='-', inplace = True)
def new_ec(row):
    val = row['EC_number']
    if row['EC_number'] == '-' and row['EC'] != '-' and row['product'] != 'hypothetical protein':
        val = row['EC'].split(',')[0]
    elif '-' in row['EC_number'] and '-' not in row['EC'] and row['product'] != 'hypothetical protein':
        val = row['EC'].split(',')[0]
    return val

df_merge['EC_new'] = df_merge.apply(new_ec, axis=1)

df_merge['EC_number'] = df_merge['EC_new']
del df_merge['EC_new']
del df_merge['EC']
if not os.path.exists(out_file): 
    df_merge.to_csv(out_file, sep='\t', index=False)
    print("Completed writing: {}".format(out_file))
else:
    print('Report exists {}'.format(out_file))