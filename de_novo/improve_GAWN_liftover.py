"""
Inputs:
1) GFF3 generated by GAWN
2) blastp format 6 file created by:
   a. Extracting transcript sequences based on same GFF3
   b. Running TransDecoder on the transcripts
   c. Running blastp of TransDecoder .pep against ref proteins (max_target_seqs = 1)
      The qlen and slen columns are also expected (-outfm6 "6 ... qlen slen")
3) TransDecoder .gff3 output from the same run

The script goes over the genes in the GAWN GFF3
and discards genes with no high quality protein product. This includes:
a. best blastp match is to the ref protein with the same name
b. match has % identity > min_identity
c. 1 - max_ratio_diff < slen/qlen < 1 + max_ratio_diff
When there is more than one high quality gene or mRNA
for a given liftover transcript, the best one is taken (based on bitscore)
In addition, the script replaces all CDS features created by GAWN
with new ones based on TransDecoder results
"""

from __future__ import print_function
import gffutils
import sys
from time import time
from os import remove
from copy import deepcopy
import itertools
import numpy as np

in_gff = sys.argv[1]
blastp_res = sys.argv[2]
transdec_gff = sys.argv[3]
min_identity = int(sys.argv[4])
max_ratio_diff = float(sys.argv[5])
out_gff = sys.argv[6]

db_path = "tmp_%s.sqlite3" % time()
gff_db = gffutils.create_db(in_gff, db_path, force=False, merge_strategy="create_unique", verbose=True)
gff = gffutils.FeatureDB(db_path)

def choose_mrna(mrna1, mrna2, min_identity, max_ratio_diff, blast_dict):
  """
  """
  best_mrna = None
  best_mrna_length = 0
  for mrna in [mrna1, mrna2]:
    if not mrna:
      continue
    if mrna not in blast_dict:
      continue
    transcript_name = mrna.split(".mrna")[0]
    # if the best match is a protein from another transcript - disqualify
    if transcript_name != blast_dict[mrna]["sseqid"]:
      continue
    perc_identity = blast_dict[mrna]["pident"]
    # if protein is too different from the ref protein
    if perc_identity < min_identity:
      continue
    qlen = blast_dict[mrna]["qlen"]
    slen = blast_dict[mrna]["slen"]
    # if query and subject lengths are too different
    if not (1 - max_ratio_diff < slen/qlen < 1 + max_ratio_diff):
      continue
    if best_mrna and blast_dict[mrna]["bitscore"] < blast_dict[best_mrna]["bitscore"]:
      continue
    # if we get here, it means this is the best mRNA so far
    best_mrna = mrna
  return best_mrna

def split_array(a, d, start=None, end=None):
  """
  Splits numpy array 'a' into sub-arrays
  every time difference between subsequent
  elements drops below d. If start and/or
  end are given, will only work on this range
  """
  a = a[start:end]
  diff = np.diff(a)
  split = np.where(diff > d)[0] + 1
  return np.split(a,split) 

# Parse blast format 6 TSV
blast_dict = {}
with open(blastp_res) as f:
  for line in f:
    fields = line.strip().split('\t')
    name = fields[0]
    blast_dict[name] = {"sseqid": fields[1], "pident": float(fields[2]), "bitscore": float(fields[11]), "qlen": int(fields[12]), "slen": int(fields[13])}

# go over all genes in the liftover gff and find the best mRNA
# There might be multiple mRNAs per input genes (.mrna1, .mrna2, ...)
# and also  multiple genes per input transcript (.path1, .path2, ...) 
genes = {}
for gene in gff.features_of_type('gene'):
  gene_name = gene['Name'][0]
  # choose best mRNA for gene (if any qualifies)
  best_mrna = None
  best_mrna_length = 0
  for mrna in gff.children(gene, featuretype='mRNA'):
    mrna_name = mrna['ID'][0]
    best_mrna = choose_mrna(best_mrna, mrna_name, min_identity, max_ratio_diff, blast_dict)
  if best_mrna:
    if gene_name not in genes:
      genes[gene_name] = [gene, best_mrna]
    else:
      # we get here if there is more than one gene from the same transcript
      # we choose the gene with the better mRNA
      better_mrna = choose_mrna(genes[gene_name][1], best_mrna, min_identity, max_ratio_diff, blast_dict)
      if better_mrna == best_mrna:
        genes[gene_name] = [gene, best_mrna]

# filter gff (take only genes with best mRNA) and fix CDSs
# start by parsing TransDecoder gff
transdec_db_path = "tmp_%s.sqlite3" % time()
transdec_gff_db = gffutils.create_db(transdec_gff, transdec_db_path, force=False, merge_strategy="create_unique", verbose=True)
transdec_gff = gffutils.FeatureDB(transdec_db_path)
transcripts = {}
for cds in transdec_gff.features_of_type('CDS'):
  transcripts[cds.seqid] = cds
# now go over all genes again. Only take genes with
# their best mRNA. Fix CDS features based on TransDecoder gff
fo = open(out_gff,'w')
for gene in gff.features_of_type('gene'):
  gene_name = gene['Name'][0]
  print(gene['ID'][0])
  if gene_name not in genes:
    continue
  best_gene = genes[gene_name][0]['ID'][0]
  if gene['ID'][0] != best_gene:
    continue
  print(gene, file=fo)
  best_mrna = genes[gene_name][1]
  for mrna in gff.children(gene, featuretype='mRNA'):
    if mrna['ID'][0] != best_mrna:
      continue
    print(mrna, file=fo)
    exons_ranges = []
    for exon in gff.children(mrna, featuretype='exon'):
      print(exon, file=fo)
      exons_ranges.append(range(exon.start, exon.end+1))
    if exons_ranges == []:
      continue
    exons_positions = [list(r) for r in exons_ranges]
    exons_positions = np.array(list(itertools.chain.from_iterable(exons_positions)))
    exons_positions = np.sort(exons_positions)
    transdec_cds = transcripts[mrna['ID'][0]]
    exons_positions_split = split_array(exons_positions,1,transdec_cds.start-1,transdec_cds.end)
    cds_i = 1
    for cds_a in exons_positions_split:
      cds = deepcopy(mrna)
      cds.featuretype = "CDS"
      cds.source = "Panoramic_fix"
      cds.start = cds_a[0]
      cds.end = cds_a[-1]
      cds['ID'][0] = [cds['ID'][0] + '.cds%s' % cds_i]
      cds['Parent'] = [mrna['ID'][0]]
      cds.attributes = {k: cds.attributes[k] for k in cds.attributes if k in {"ID","Name","Parent"}}
      print(cds, file=fo)
      cds_i += 1

fo.close()
remove(db_path)
remove(transdec_db_path)
