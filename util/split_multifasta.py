"""
Split a multi-record fasta into
multiple single-record fasta files.
if create_dirs is on, put each file
in a new dir.
"""

import sys
import os
from Bio import SeqIO

in_fasta = sys.argv[1]
out_dir = sys.argv[2]
create_dirs = False
if len(sys.argv) > 3:
  create_dirs = True
   
for rec in SeqIO.parse(in_fasta, 'fasta'):
  rec_id = rec.id
  if create_dirs:
    rec_dir = os.path.join(out_dir, rec_id)
    if not os.path.isdir(rec_dir):	# if dir doesn't exist
      os.mkdir(rec_dir)
    rec_fasta = os.path.join(rec_dir, rec_id+'.fa')
  else:
    rec_fasta = os.path.join(out_dir, rec_id+'.fa')
  with open(rec_fasta, 'w') as fo:
    print(rec.format('fasta'), file=fo)
