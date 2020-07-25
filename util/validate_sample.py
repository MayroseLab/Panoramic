"""
==========================================
Use this script to validate reference
and HQ sample inputs. Specifically,
checks that protein and transcript
FASTA headers match the ID attribute of
mRNA features in the GFF3 file.
Usage:
python valudate_sample.py <GFF3> <proteins.fasta> <transcripts.fasta>
==========================================
"""

from __future__ import print_function
import argparse

parser = argparse.ArgumentParser(description='Validate sample')
parser.add_argument('gff3', help="Path to GFF3 file")
parser.add_argument('prot', help="Path to proteins FASTA")
parser.add_argument('trans', help="Path to transcripts FASTA")

args = parser.parse_args()

# get mRNA ID list
ids = set()
with open(args.gff3) as f:
  for line in f:
    if line.startswith('#'):
      continue
    fields = line.strip().split('\t')
    if fields[2] != "mRNA":
      continue
    attr = {k.split('=')[0]: k.split('=')[1] for k in fields[8].split(';')}
    if attr['ID'] in ids:
      print("PROBLEM: ID %s encountered is non-unique" % attr['ID'])
    ids.add(attr['ID'])
print("%s unique mRNA IDs were detected" % len(ids))

# check matching with proteins and transcripts
for fasta in [args.prot, args.trans]:
  print("Validating file %s" % fasta)
  names = set()
  with open(fasta) as f:
    for line in f:
      if not line.startswith('>'):
        continue
      name = line.strip()[1:]
      names.add(name)
  missing_from_fasta = ids - names
  missing_from_gff = names - ids
  print("%s IDs found in GFF3 but not in FASTA:" % len(missing_from_fasta))
  print(', '.join(list(missing_from_fasta)))
  print("%s IDs found in FASTA but not in GFF3:" % len(missing_from_gff))
  print(', '.join(list(missing_from_gff)))
