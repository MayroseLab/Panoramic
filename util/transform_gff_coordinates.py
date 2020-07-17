"""
Takes a gff file and a bed file which allows
conversion of sequence names and coordinates
"""
from __future__ import print_function
import sys

in_gff = sys.argv[1]
in_bed = sys.argv[2]
out_gff = sys.argv[3]

# parse bed file
chunks = {}
with open(in_bed) as f:
  for line in f:
    fields = line.strip().split('\t')
    chunks[fields[3]] = [fields[0], int(fields[1])]

with open(in_gff) as f, open(out_gff,'w') as fo:
  for line in f:
    line = line.strip()
    if line.startswith('#'):
      print(line,file=fo)
      continue
    fields = line.split('\t')
    if len(fields) != 9:
      print(line,file=fo)
      continue
    if fields[2] == "contig":
      continue
    rec_name = fields[0]
    if rec_name not in chunks:
      print(line,file=fo)
      continue
    real_rec_name = chunks[rec_name][0]
    real_rec_start = chunks[rec_name][1]
    fields[0] = real_rec_name
    fields[3] = str(int(fields[3]) + real_rec_start)
    fields[4] = str(int(fields[4]) + real_rec_start)
    print('\t'.join(fields), file=fo)
