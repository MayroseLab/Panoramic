"""
Takes a directory path with orderings files
created by RaGOO, and produces a bed file with
contig borders. Also requires a fasta index
(faidx) of the contigs and the padding gap used
when RaGOO was run.
"""

from __future__ import print_function
import sys
import os

orderings_dir = sys.argv[1]
faidx = sys.argv[2]
gap = int(sys.argv[3])
out_bed = sys.argv[4]

# get contig lengths from faidx
with open(faidx) as f:
  ctg_len = { line[0]: int(line[1]) for line in [line.strip().split('\t') for line in f]}

# go over all orderings files and create bed
ord_files = [ os.path.join(orderings_dir, f) for f in os.listdir(orderings_dir) if f.endswith('_orderings.txt')]
with open(out_bed,'w') as fo:
  for ordf in ord_files:
    chrom = os.path.basename(ordf).replace('_orderings.txt','') + '_RaGOO'
    with open(ordf) as f:
      start_pos = 0
      for line in f:
        contig = line.strip().split('\t')[0]
        contig_len = ctg_len[contig]
        end_pos = start_pos + contig_len
        print('\t'.join([chrom,str(start_pos),str(end_pos),contig]), file=fo)
        start_pos = end_pos + gap
