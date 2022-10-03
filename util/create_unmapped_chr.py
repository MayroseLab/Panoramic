"""
Concat unmapped scaffolds in a
fasta file to create an "unmapped chromosome".
If a corresponding gff3 file if also given,
convert relevant coordinates.
"""

import sys
import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import argparse
import re

def create_um_chr(genome, pattern, out_fasta, gap_size, um_chr_name):
  """
  Concatenate unmapped scaffolds.
  Return a dict with scaffold
  possitions on unmapped chr.
  """

  regex = re.compile(pattern)
  um_scaffolds = []

  if out_fasta == '-':
    fo = sys.stdout
  else:
    fo = open(out_fasta,'w')
  for rec in SeqIO.parse(genome, 'fasta'):
    if not regex.match(rec.id):
      print(rec.format('fasta'), file=fo)
    else:
      um_scaffolds.append(rec)
  un_chr_seq = ('N'*gap_size).join([str(u.seq) for u in um_scaffolds])
  um_chr_rec = SeqRecord(id=um_chr_name, name='', description='', seq=Seq(un_chr_seq))
  print(um_chr_rec.format('fasta'), file=fo)
  if out_fasta != '-':
    fo.close()

  um_scaffold_coords = {}
  i = 1
  for u in um_scaffolds:
    um_scaffold_coords[u.id] = i
    i += len(u.seq) + gap_size
  return um_scaffold_coords

def convert_coords(gff3, gff3_out, um_chr_name, um_scaffolds):
  """
  Convert coordinates of features
  on unmapped scaffolds to unmapped
  chromosome coordinates.
  """
  with open(gff3) as f, open(gff3_out,'w') as fo:
    for line in f:
      line = line.strip()
      if line.startswith('#'):
        print(line, file=fo)
        continue
      fields = line.split('\t')
      if fields[0] not in um_scaffolds:
        print(line, file=fo)
        continue
      # if feature on unmapped scaffold
      fields[3] = str(int(fields[3]) + um_scaffolds[fields[0]] - 1)
      fields[4] = str(int(fields[4]) + um_scaffolds[fields[0]] - 1)
      fields[0] = um_chr_name
      print('\t'.join(fields), file=fo)

if __name__ == "__main__":

  parser = argparse.ArgumentParser(description='Create Unmapped Chromosome')
  parser.add_argument('genome', help='Genome fasta file')
  parser.add_argument('pattern', help='Pattern to use for identifying unmapped scaffolds')
  parser.add_argument('out_fasta', help='Output fasta file')
  parser.add_argument('--gap_size', '-n', type=int, default=100, help='Number of Ns to separate between unmapped scaffolds')
  parser.add_argument('--um_chr_name', '-u', default='Chr0', help='Name of unmapped chromosome')
  parser.add_argument('--gff3', '-a', default=None, help='Annotation gff3 file')
  parser.add_argument('--gff3_out', '-o', default=None, help='Output gff3 file')
  args = parser.parse_args()

  um_scaffolds = create_um_chr(args.genome, args.pattern, args.out_fasta, args.gap_size, args.um_chr_name)
  if args.gff3:
    convert_coords(args.gff3, args.gff3_out, args.um_chr_name, um_scaffolds)
