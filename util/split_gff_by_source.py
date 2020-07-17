"""
Takes a GFF3 file and a list of
sources. Creates new GFF3 files, each
with features of one source.
Sources not in the list are ignored.
"""

from __future__ import print_function
import sys
import os

in_gff = sys.argv[1]
sources_list = sys.argv[2]
out_dir = sys.argv[3]

sources = sources_list.split(',')
gff_base, gff_ext = os.path.splitext(os.path.basename(in_gff))
file_handles = {}
for src in sources:
  file_path = os.path.join(out_dir, gff_base + ".%s%s" %(src, gff_ext)) 
  file_handle = open(file_path,'w')
  file_handles[src] = file_handle

def write_to_all(x, file_handles):
  for fh in file_handles:
    print(x, file=file_handles[fh])

with open(in_gff) as f:
  for line in f:
    line = line.strip()
    if line.startswith('#'):
      write_to_all(line, file_handles)
    else:
      line_source = line.split('\t')[1]
      if line_source in file_handles:
        print(line, file=file_handles[line_source])

for fh in file_handles.values():
  fh.close()
