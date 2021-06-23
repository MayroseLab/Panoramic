"""
Create the partitions list file
for EVM. Each line has the structure:
chromosome	chromosome_dir	Y	partition_dir (chromosome_dir/chromosome_start-end)
"""

import sys
from Bio import SeqIO

genome_fasta = sys.argv[1]
segment_size = int(sys.argv[2])
overlap_size = int(sys.argv[3])
out_dir = sys.argv[4]
partitions_out = sys.argv[5]

# get chromosome lengths
chr_lens = {}
for rec in SeqIO.parse(genome_fasta, 'fasta'):
  chr_lens[rec.id] = len(rec.seq)

# compute partitions
with open(partitions_out, 'w') as fo:
  for chrom in sorted(chr_lens.keys()):
    chrom_dir = out_dir + '/' + chrom
    partitions = ["%s-%s" %(start, min(start+segment_size-1, chr_lens[chrom])) for start in range(1, chr_lens[chrom], segment_size - overlap_size) if (min(start+segment_size-1, chr_lens[chrom]) - start) > overlap_size]
    for p in partitions:
      partition_dir = chrom_dir + '/' + "%s_%s" %(chrom, p)
      print('\t'.join([chrom, chrom_dir, 'Y', partition_dir]), file=fo)
