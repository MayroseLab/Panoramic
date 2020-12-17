"""
Filter contigs based on length
"""

from Bio import SeqIO
import sys

in_fasta = sys.argv[1]
min_len = int(sys.argv[2])
out_fasta = sys.argv[3]

out_recs = []
for rec in SeqIO.parse(in_fasta,'fasta'):
  contig_len = len(rec.seq)
  if contig_len > min_len:
    out_recs.append(rec)
SeqIO.write(out_recs, out_fasta, "fasta")
