"""
Filter SPAdes contigs based on length and coverage
"""

from Bio import SeqIO
import sys

in_fasta = sys.argv[1]
min_len = int(sys.argv[2])
min_cov = int(sys.argv[3])
out_fasta = sys.argv[4]

out_recs = []
for rec in SeqIO.parse(in_fasta,'fasta'):
  desc = rec.description.split('_')
  contig_len = int(desc[3])
  contig_cov = float(desc[5])
  if contig_len > min_len and contig_cov > min_cov:
    out_recs.append(rec)
SeqIO.write(out_recs, out_fasta, "fasta")
