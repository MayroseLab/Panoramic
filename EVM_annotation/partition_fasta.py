"""
Given a fasta file and the EVM partitions list,
split the fasta into multiple partitions.
"""

import sys
import os
from Bio import SeqIO

in_fasta = sys.argv[1]
partitions_list = sys.argv[2]

# parse partitions list
partitions_dict = {}
with open(partitions_list) as f:
  for line in f:
    partition_path = line.strip().split('\t')[3]
    partition_base = os.path.basename(partition_path)
    chrom, coords = partition_base.rsplit('_',1)
    start, end = [int(x) for x in coords.split('-')]
    if chrom not in partitions_dict:
      partitions_dict[chrom] = {}
    partitions_dict[chrom][(start, end)] = partition_path

# parse fasta
genome_dict = SeqIO.to_dict(SeqIO.parse(in_fasta, "fasta"))

# create partitions
fasta_basename = os.path.basename(in_fasta)
for chrom in sorted(partitions_dict.keys()):
  for p in partitions_dict[chrom]:
    start, end = p
    partition_rec = genome_dict[chrom][start-1:end]
    partition_rec.description = ''
    # write fasta
    partition_dir = partitions_dict[chrom][p]
    if not os.path.isdir(partition_dir):
      os.makedirs(partition_dir)
    partition_fasta = os.path.join(partition_dir, fasta_basename)
    with open(partition_fasta, 'w') as fo:
      SeqIO.write([partition_rec], fo, "fasta")
