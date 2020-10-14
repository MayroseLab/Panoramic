"""
Filter a fasta file according to the
features in a gff file. Outputs a new
fasta.
"""

from __future__ import print_function
import gffutils
from Bio import SeqIO
import sys
from os import remove
from time import time
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('in_gff', help='input GFF3 file')
parser.add_argument('in_fasta', help='input FASTA file')
parser.add_argument('out_fasta', help='output filtered FASTA file')
parser.add_argument('feature_type', help='feature type to look at in GFF3 file (e.g. gene or mRNA)')
parser.add_argument('name_attribute', help='attribute in GFF3 file that identifies the record (e.g. ID or transcript_id)')
parser.add_argument('--min_len', default=0, type=int, help='discard fasta records shorter than min_len')
args = parser.parse_args()

in_gff = args.in_gff
in_fasta = args.in_fasta
out_fasta = args.out_fasta
feature_type = args.feature_type
name_attribute = args.name_attribute

# create list of desired names
db_path = "tmp_%s.sqlite3" % time()
try:
  gff_db = gffutils.create_db(in_gff, db_path, force=True, merge_strategy="create_unique")
  gff = gffutils.FeatureDB(db_path)
  names = {feat.attributes[name_attribute][0] for feat in gff.features_of_type(feature_type)}
except ValueError:
  print("WARNING: gff file %s seems to be empty. Therefore, no fasta records were written." % args.in_gff)
  names = []

remove(db_path)

# go over fasta and select records by name
records = []
with open(in_fasta) as f:
  for rec in SeqIO.parse(f, "fasta"):
    if rec.id in names and len(rec.seq) >= args.min_len:
      rec.description = ''
      records.append(rec)

with open(out_fasta, 'w') as fo:
  SeqIO.write(records, fo, "fasta")
