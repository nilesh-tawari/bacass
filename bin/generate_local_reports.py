#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 15 17:55:42 2020

@author: nilesh-tawari <tawari.nilesh@gmail.com>
"""
import os
import os.path as op
import io
import string
import json
import ntpath
import argparse

'''
USAGE: python  generate_local_reports.py -r PATH_FOR_BAGEL4RESULTS
'''

# Parameters
parser = argparse.ArgumentParser(description = 'Create standalone HTML report files from output of Bagel4')
parser.add_argument('-r', '--bagel4results', help='Full path for Bagel4 results directory')


args = parser.parse_args()
bagel4results = os.path.abspath(args.bagel4results)

###############################################################################

def path_leaf(path):
    """
    split the path in tail and head
    str -> str, str
    """
    head, tail = ntpath.split(path)
    return head, tail or ntpath.basename(head)

inp_file = os.path.join(bagel4results, '00.OverviewGeneTables.json')
out_idx_file = os.path.join(bagel4results, 'index.html')
templates_dir = op.abspath(op.join(op.dirname(__file__), 'templates'))
###############################################################################
vals = {}

# create index.html from OverviewGeneTables.json file
ovhtml_tmpl = string.Template(open(op.join(templates_dir, "OverviewTable_template.html")).read())

with open(inp_file) as OverviewGeneTables:
  OverviewGeneTables = json.load(OverviewGeneTables)

if 'ResultsTable' in OverviewGeneTables:  
    all_aoi = []
    #replace href links as individual html result files
    for i in range(len(OverviewGeneTables['ResultsTable'])):
        aoi_id = OverviewGeneTables['ResultsTable'][i]['AOI'].split('><b>')[1].split('</b></a>')[0]
        entry = '<a href={0}.html target=_blank><b>{1}</b></a>'.format(aoi_id, aoi_id)
        OverviewGeneTables['ResultsTable'][i]['AOI'] = entry
        all_aoi.append(aoi_id)
        vals[aoi_id] = aoi_id
    vals['OverviewGeneTables'] = OverviewGeneTables
else:
    vals['OverviewGeneTables'] = OverviewGeneTables

# substitute vals in index template
tmpl = ovhtml_tmpl.substitute(**vals)
with io.open(out_idx_file, "w", encoding='utf8') as fh:
    fh.write(tmpl)
    
if 'ResultsTable' in vals['OverviewGeneTables']:
    # create individual GeneTableSingleGraph files
    gthtml_tmpl = string.Template(open(op.join(templates_dir, "GeneTableSingleGraph_template.html")).read())
    
    for aoi in all_aoi:
        print('Generating report for {0}'.format(aoi))
        sp_vals = {'d1': '$(document)', 'd2': '$( window )', 'aoi_id': aoi, 'blast_hit': ''}
        out_aoi_file = os.path.join(bagel4results, '{0}.html'.format(aoi))
        #BVEN24_C08.11.AOI_01.GeneTable.json
        with open(os.path.join(bagel4results, aoi + '.GeneTable.json')) as GeneTable:
            sp_vals['GeneTable'] = json.load(GeneTable)
        for i in range(len(sp_vals['GeneTable']['AOIs'])):
            results_dir, filename_in = path_leaf(sp_vals['GeneTable']['AOIs'][i]['filename'])
            sp_vals['GeneTable']['AOIs'][i]['filename'] = filename_in # replace full path with filename 
            
        #BVEN24_C08.11.AOI_01.prom_term.json
        with open(os.path.join(bagel4results, aoi + '.prom_term.json')) as prom_term:
            sp_vals['dataPT'] = json.load(prom_term)
        try:
            #BVEN24_C04.8.AOI_01.orf00022.blast_1.json
            blast_1 = [file for file in os.listdir(bagel4results) if '.blast_1.json' in file and aoi in file]
            blast_1_file = open(os.path.join(bagel4results, blast_1[0]))
            sp_vals['blast_hit'] = json.load(blast_1_file)['blasthits'][0]['htmlrecord']
        except:
            pass
        try:
            #BVEN24_C04.8.AOI_01.orf00022.blast_2.json
            blast_2 = [file for file in os.listdir(bagel4results) if '.blast_2.json' in file and aoi in file]
            blast_2_file = open(os.path.join(bagel4results, blast_2[0]))
            sp_vals['blast_hit'] = json.load(blast_2_file)['blasthits'][0]['htmlrecord']
        except:
            pass
        
        # substitute vals in individual template
        tmpl = gthtml_tmpl.substitute(**sp_vals)
        with io.open(out_aoi_file, "w", encoding='utf8') as fh:
            fh.write(tmpl)
        
print("Completed creating HTML reports in {0} directory!!!".format(bagel4results))