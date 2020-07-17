"""
Filter MAKER proteins - only keep proteins with AED <= X
"""

from Bio import SeqIO
import sys

in_fasta = sys.argv[1]
max_aed = float(sys.argv[2])
out_fasta = sys.argv[3]

out_recs = []
for rec in SeqIO.parse(in_fasta,'fasta'):
  aed = float(rec.description.split(' ')[2].split(':')[1])
  if aed <= max_aed:
    out_recs.append(rec)
SeqIO.write(out_recs, out_fasta, "fasta")
