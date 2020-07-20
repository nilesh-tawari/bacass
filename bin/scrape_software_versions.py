#!/usr/bin/env python
from __future__ import print_function
from collections import OrderedDict
import re

regexes = {
    'nf-core/bacass': ['v_pipeline.txt', r"(\S+)"],
    'Nextflow': ['v_nextflow.txt', r"(\S+)"],
    'FastQC': ['v_fastqc.txt', r"FastQC v(\S+)"],
    'MultiQC': ['v_multiqc.txt', r"multiqc, version (\S+)"],
    'Skewer': ['v_skewer.txt', r"version: (\S+)"],
    'Kraken2': ['v_kraken2.txt', r"([\d\.]+)-beta"],
    'Quast': ['v_quast.txt', r"QUAST v(\S+)"],
    'Prokka': ['v_prokka.txt', r"prokka (\S+)"],
    'Bandage': ['v_bandage.txt', r"Version: (\S+)"]}
'''
,
    'Nanopolish': ['v_nanopolish.txt', r"nanopolish version (\S+)"],
    'Miniasm': ['v_miniasm.txt', r"(\S+)"],
    'Racon': ['v_racon.txt', r"v(\S+)"],
    'Canu': ['v_canu.txt', r"snapshot (\S+)"],
    'Samtools': ['v_samtools.txt', r"samtools (\S+)"],
    'Minimap2': ['v_minimap2.txt', r"(\S+)"],
    'NanoPlot': ['v_nanoplot.txt', r"NanoPlot (\S+)"],
}
'''
results = OrderedDict()
results['nf-core/bacass'] = '<span style="color:#999999;\">N/A</span>'
results['Nextflow'] = '<span style="color:#999999;\">N/A</span>'
results['FastQC'] = '<span style="color:#999999;\">N/A</span>'
results['MultiQC'] = '<span style="color:#999999;\">N/A</span>'
results['Skewer'] = '<span style="color:#999999;\">N/A</span>'
results['Kraken2'] = '<span style="color:#999999;\">N/A</span>'
results['Quast'] = '<span style="color:#999999;\">N/A</span>'
results['Prokka'] = '<span style="color:#999999;\">N/A</span>'
results['Bandage'] = '<span style="color:#999999;\">N/A</span>'
#results['Porechop'] = '<span style="color:#999999;\">N/A</span>'
#results['Nanopolish'] = '<span style="color:#999999;\">N/A</span>'
#results['Miniasm'] = '<span style="color:#999999;\">N/A</span>'
#results['Racon'] = '<span style="color:#999999;\">N/A</span>'
#results['Canu'] = '<span style="color:#999999;\">N/A</span>'
#results['DFAST'] = '<span style="color:#999999;\">N/A</span>'
#results['Samtools'] = '<span style="color:#999999;\">N/A</span>'
#results['Minimap2'] = '<span style="color:#999999;\">N/A</span>'
#results['NanoPlot'] = '<span style="color:#999999;\">N/A</span>'


# Search each file using its regex
for k, v in regexes.items():
    with open(v[0]) as x:
        versions = x.read()
        match = re.search(v[1], versions)
        if match:
            results[k] = "v{}".format(match.group(1))

# Remove software set to false in results
for k in results:
    if not results[k]:
        del(results[k])

# Dump to YAML
print ('''
id: 'software_versions'
section_name: 'nf-core/bacass Software Versions'
section_href: 'https://github.com/nf-core/bacass'
plot_type: 'html'
description: 'are collected at run time from the software output.'
data: |
    <dl class="dl-horizontal">
''')
for k,v in results.items():
    print("        <dt>{}</dt><dd><samp>{}</samp></dd>".format(k,v))
print ("    </dl>")

# Write out regexes as csv file:
with open('software_versions.csv', 'w') as f:
    for k,v in results.items():
        f.write("{}\t{}\n".format(k,v))
