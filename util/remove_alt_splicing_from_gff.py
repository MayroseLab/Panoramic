"""
Takes a gff file and for genes with
multiple mRNAs, only keeps the longest
(largest sum of exon lengths). Prints
out a new gff.
If a proteins fasta is also provided,
removes genes shorter than a given
cutoff
"""

from __future__ import print_function
import gffutils
import sys
from os import remove
from time import time
from Bio import SeqIO

in_gff = sys.argv[1]
out_gff = sys.argv[2]
gene_mrna_out = out_gff + '.gene_to_mRNA'
if len(sys.argv) == 6:
  proteins_fasta = sys.argv[3]
  min_len = int(sys.argv[4])
  name_attribute = sys.argv[5] # mRNA attribute where the protein name is stored (must mastch with fasta)
else:
  proteins_fasta = None

db_path = "tmp_%s.sqlite3" % time()
gff_db = gffutils.create_db(in_gff, db_path, force=True, merge_strategy="create_unique")
gff = gffutils.FeatureDB(db_path)

prot_lens = {}
if proteins_fasta:
  # go over fasta and save protein lengths
  prot_lens = {rec.id: len(rec.seq) for rec in SeqIO.parse(proteins_fasta, "fasta")}

with open(out_gff, 'w') as fo, open(gene_mrna_out,'w') as fo2:
  for feature in gff.all_features():
    if feature.featuretype not in {'gene', 'mRNA', 'exon', 'CDS', 'five_prime_UTR', 'three_prime_UTR'}:
      print(str(feature), file=fo)
      continue
    if feature.featuretype != 'gene':	# mRNA and exon features
      continue
    gene = feature
    print(str(gene), file=fo)
    mrnas = list(gff.children(gene, featuretype='mRNA'))
    longest_transcript = mrnas[0]
    max_len = 0
    for mrna in mrnas:
      exons = list(gff.children(mrna, featuretype='exon'))
      total_exons_len = sum([exon.end - exon.start for exon in exons])
      if total_exons_len > max_len:
        max_len = total_exons_len
        longest_transcript = mrna
    if proteins_fasta:
      transcript_name = longest_transcript[name_attribute][0]
      if transcript_name in prot_lens and prot_lens[transcript_name] < min_len:
        continue
    print(str(longest_transcript), file=fo)
    print("%s\t%s" %(gene['ID'][0], longest_transcript['ID'][0]),file=fo2)
    for x in gff.children(longest_transcript):
      print(str(x), file=fo)
        
remove(db_path)
