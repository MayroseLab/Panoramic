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

in_gff = sys.argv[1]
in_fasta = sys.argv[2]
out_fasta = sys.argv[3]
feature_type = sys.argv[4]  # e.g. gene or mRNA
name_attribute = sys.argv[5] # e.g. ID or transcript_id

# create list of desired names
db_path = "tmp_%s.sqlite3" % time()
gff_db = gffutils.create_db(in_gff, db_path, force=True, merge_strategy="create_unique")
gff = gffutils.FeatureDB(db_path)

names = {feat.attributes[name_attribute][0] for feat in gff.features_of_type(feature_type)}
remove(db_path)

# go over fasta and select records by name
records = []
with open(in_fasta) as f:
  for rec in SeqIO.parse(f, "fasta"):
    if rec.id in names:
      rec.description = ''
      records.append(rec)

with open(out_fasta, 'w') as fo:
  SeqIO.write(records, fo, "fasta")
