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

OG_seq_dir = sys.argv[1]
in_mapping_file = sys.argv[2]	# OG to gene name table
in_OG_file = sys.argv[3]	# OGS in OrthoFinder format
out_fasta = sys.argv[4]

# Create dict {OG1 -> [g1,g2,g3], OG2 -> [g4,g5]}
og_dict = {}
with open(in_OG_file) as f:
  f.readline()	# skip header
  for line in f:
    fields = line.strip().split('\t')
    og_dict[fields[0]] = list(chain(*[ x.split(', ') for x in fields[1:] if x ]))

# find protein sequences and select one per gene
records = []
with open(in_mapping_file) as f:
  f.readline()	# skip header
  for line in f:
    og_name, gene_name = line.strip().split('\t')
    og_orig_name = og_name.split('.')[0]
    og_file = os.path.join(OG_seq_dir, og_orig_name + '.fa')
    og_records_dict = SeqIO.to_dict(SeqIO.parse(og_file, "fasta"))
    og_genes_list = og_dict[og_name]
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
    records.append(rep_record)

# print representative sequences to fasta
SeqIO.write(records, out_fasta, "fasta")
