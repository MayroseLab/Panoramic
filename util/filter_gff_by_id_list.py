"""
Takes a gff3 file and a file with
list of IDs (one per line).
Output is a gff3 file containing
only features with ID or Parent
from the IDs list.
"""

from __future__ import print_function
import sys

in_gff = sys.argv[1]
in_list = sys.argv[2]
out_gff = sys.argv[3]

with open(in_list) as f:
  ids_list = set([l.strip() for l in f.readlines()])

with open(in_gff) as f, open(out_gff, 'w') as fo:
  for line in f:
    line = line.strip()
    if line.startswith('#'):
      print(line, file=fo)
      continue
    if line.endswith(';'):
      line = line[:-1]
    attr = {k.split('=')[0]: k.split('=')[1] for k in line.split('\t')[8].split(';')}
    if attr['ID'] in ids_list or ('Parent' in attr and attr['Parent'] in ids_list):
      print(line, file=fo)
