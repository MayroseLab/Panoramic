"""
Create a fasta file with one protein
sequence per gene in the pan genome.
In case a gene represents an OG with
multiple genes, take the reference
sequence if one exists, otherwise
take the longest protein.
"""

from Bio import SeqIO
import sys
import os
from itertools import chain
import csv

OG_seq_dir = sys.argv[1]
in_mapping_file = sys.argv[2]	# OG to gene name table
in_OG_file = sys.argv[3]	# OGS in OrthoFinder format
out_mapping_file = sys.argv[4]	# same as input mapping + name of selected representative

# Create dict {OG1 -> {sample1: [g1], sample2: [g2,g3]}, OG2 -> {sample2: [g4], sample3: [g5]}}
og_dict = {}
with open(in_OG_file) as f:
  reader = csv.DictReader(f, delimiter='\t')
  for row in reader:
    og_dict[row['Orthogroup']] = { sample: set(row[sample].split(', ')) for sample in row.keys() if (sample != 'Orthogroup' and row[sample]) }

# find protein sequences and select one per gene
with open(in_mapping_file) as f, open(out_mapping_file,'w') as fo:
  f.readline()	# skip header
  for line in f:
    og_name, gene_name = line.strip().split('\t')
    print(og_name)
    og_orig_name = og_name.split('.')[0]
    og_file = os.path.join(OG_seq_dir, og_orig_name + '.fa')
    og_records_dict = SeqIO.to_dict(SeqIO.parse(og_file, "fasta"))
    og_genes_list = set(chain(*og_dict[og_name].values()))
    og_records_dict = { g: og_records_dict[g] for g in og_genes_list }
    # if OG has a reference gene
    if not gene_name.startswith('PanGene'):
      rep_name = gene_name
    # if not
    else:
      rep_name = sorted(list(og_records_dict.keys()), key=lambda x: len(og_records_dict[x].seq))[-1]
    rep_record = og_records_dict[rep_name]
    rep_record.id = gene_name
    rep_record.description = ''
    # find the sample from which rep_name originates
    for sample in og_dict[og_name]:
      if rep_name in og_dict[og_name][sample]:
        rep_sample = sample
        break
    # print selected representatives to file
    print("%s\t%s\t%s\t%s" %(og_name, gene_name, rep_sample, rep_name), file=fo)
