#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 11 16:25:25 2020

@author: nilesh-tawari <tawari.nilesh@gmail.com>
"""
import os
import os.path
import sys
import ntpath
import argparse
import pandas as pd
from Bio import SeqIO
try:
    # Python 2
    xrange
except NameError:
    # Python 3, xrange is now named range
    xrange = range

'''
USAGE: python summarize_abricate.py -d 1016_abricate -g 1016.gbk -o . -p 1016
'''
# Parameters
parser = argparse.ArgumentParser(description = 'Create aggregated excel report for ABRicate by searching current directory')
parser.add_argument('-d', '--abricate_results', help='Full path for ABRicate results directory')
parser.add_argument('-g', '--prokka_gbk',  help='Full path for Prokka GBK file')
parser.add_argument('-o', '--out_dir', help='Path of output directory')
parser.add_argument('-p', '--out_prefix',  help='Output filename prefix')

args = parser.parse_args()
path = os.path.abspath(args.abricate_results)
prokka_gbk = args.prokka_gbk
out_dir = os.path.abspath(args.out_dir)
out_prefix = args.out_prefix

###########################################################################
def path_leaf(path):
    """
    split the path in tail and head
    str -> str, str
    """
    head, tail = ntpath.split(path)
    return head, tail or ntpath.basename(head)

files = []
dfs = {}

file_ext = ['.abr.summ.txt', '.abr.plasmidfinder.out', '.abr.argannot.out', '.abr.ncbi.out', \
         '.abr.resfinder.out', '.abr.vfdb.out']
 
for dirpath, dirnames, filenames in os.walk(path):
    filenames = [file for file in filenames if any(t in file for t in file_ext)]
    for filename in filenames:
        if filename not in files:
            files.append(filename)
        else:
            print('ERROR: Duplicate filename {0}'.format(filename))
            sys.exit()

for file in files:
    full_path = os.path.join(dirpath, file)
    print('CURATED: {0}'.format(full_path))
    df = pd.read_csv(full_path, sep='\t', chunksize=1000, \
                low_memory=False, iterator = True)
    df = pd.concat(list(df), ignore_index=True)
    if 'COVERAGE_MAP' in df.columns:
        df['COVERAGE_MAP'] = "'" +  df['COVERAGE_MAP'] + "'"
    if not file.endswith('.abr.summ.txt'):
        for gb_record in SeqIO.parse(open(prokka_gbk,"r"), "genbank") :
            # now do something with the record
            #print("Name %s, %i features" % (gb_record.name, len(gb_record.features)))
            for i in range(len(df)):
                if gb_record.name == df.iloc[i].SEQUENCE:
                    start, end = [df.iloc[i].START, df.iloc[i].END]
                    desired = set(xrange(int(start),int(end),1))
                    for f in gb_record.features:
                        span = set(xrange(f.location._start.position, f.location._end.position))
                        if span & desired:
                            matched_feature = f
                            s = matched_feature.location.start.position + 1
                            e = matched_feature.location.end.position
                            st = matched_feature.location.strand
                            if st == 1:
                                st = '+'
                            else:
                                st = '-'
                            df.at[i, 'PROKKA FEATURE'] = '{0}:{1}-{2}({3})'.format(gb_record.name,s,e,st)
                            for key, value in matched_feature.qualifiers.items():
                                if key == 'locus_tag':
                                    df.at[i, 'LOCUS_TAG'] = ', '.join(value)
                                if key == 'inference':
                                    df.at[i, 'INFERENCE'] = ', '.join(value)
                                if key == 'product':
                                    df.at[i, 'PROKKA PRODUCT'] = ', '.join(value)
    dfs[file] = df        
    
# writing excel 
def set_format(df, worksheet1):
    '''
    set column width in excel sheet based on len(column)
    df ->
    '''
    for i, col in enumerate(df.columns):
        column_len = df[col].astype(str).str.len().max()
        column_len = max(column_len, len(col)) + 2
        if column_len > 20:
            column_len = len(col) +2
        worksheet1.set_column(i,i,column_len)

out_file_name = out_prefix + '.ABRicate.report.xlsx' 
out_file = os.path.join(out_dir, out_file_name)


sort_map = {day: next((idx for idx, val in enumerate(file_ext) if val in day),
                 len(files)) for day in files}
files = sorted(files, key=sort_map.__getitem__)

if not os.path.exists(out_file): 
    writer = pd.ExcelWriter(out_file, engine = 'xlsxwriter')
    for df_name in files:
        #print(df_name)
        df_write = dfs[df_name] 
        df_write.to_excel(writer, df_name[:30], index = False, startrow = 0)#, float_format ="%.2g")
        workbook = writer.book
        format_highlight = workbook.add_format({'bg_color': '#FFC7CE', 'font_color': '#9C0006'})
        worksheet1 = writer.sheets[df_name[:30]]
        worksheet1.freeze_panes(1,1)
        # Highlight the > 90% values in Red
        if not df_name.endswith('.abr.summ.txt'):
            worksheet1.freeze_panes(1,6)
            worksheet1.conditional_format('J2:K{0}'.format(len(df_write)+1), {'type': 'cell',
                                                   'criteria' : '>=',
                                                   'value': 90,
                                                   'format': format_highlight})
        
        set_format(df_write, worksheet1)
    writer.save()
    print('Completed writing report {0} !!!'.format(out_file))
else:
    print('Report exists {0} !!!'.format(out_file))
 
        
        