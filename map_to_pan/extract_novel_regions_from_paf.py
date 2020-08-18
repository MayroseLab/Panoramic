"""
This script is intended for extracting
novel genomic regions present in the
query genome but not in the target genome.
The expected input is a PAF file generated
using minimap2 -x asm5 -c .
Regions not mapped at all + large insertions
within mapped regions are written as fasta
records. In addition, if a gff file is
provided, genes located on novel regions are
extracted.
"""

from __future__ import print_function
import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from pafpy import PafFile
from intervaltree import IntervalTree, Interval
import re
import gffutils
from os import remove

parser = argparse.ArgumentParser()
parser.add_argument('in_paf', help='Input PAF file')
parser.add_argument('in_fasta', help='Input query genome fasta')
parser.add_argument('min_region', type=int, help='Min size for extracted novel regions')
parser.add_argument('--in_gff', default=None, help='Input gff annotation of query genome')
parser.add_argument('out_fasta', help='Output fasta file')
parser.add_argument('--out_gff', default=None, help='Output gff file')
parser.add_argument('--genome_name', default='', help='Name to include in contig names')
args = parser.parse_args()
if args.genome_name:
  args.genome_name += "_"

# Read input fasta
genome_dict = SeqIO.to_dict(SeqIO.parse(args.in_fasta, "fasta"))
# create per-chromosome intervals
# this is used for finding unmapped regions of the query genome
# initialize by creating full-length chromosome intervals
# and chop whenerver mapped regions are detected while parsing PAF
chrom_intervals_dict = { chrom: IntervalTree([Interval(0,len(genome_dict[chrom].seq)-1)]) for chrom in genome_dict }

# Read PAF and find novel regions
with PafFile(args.in_paf) as paf:
  for record in paf:
    # if entire scaffold is unmapped
    if record.tname == "*":
      continue
    # if at least part of the record is mapped
    # chop mapped region from the chromosome interval
    # fetch CIGAR and use it to extract matched regions
    cigar = str(record.tags["cg"])[5:]
    cigar_tuples = re.findall(r'(\d+)([A-Z]{1})', cigar)
    start = record.qstart
    for ct in cigar_tuples:
      if ct[1] in {"D","N"}:
        continue
      if ct[1] == "M":
        end = start + int(ct[0])
        chrom_intervals_dict[record.qname].chop(start,end)
      start += int(ct[0])

  # now extract all unmapped regions (larger than cutoff)
  out_fasta_records = []
  for chrom in chrom_intervals_dict:
    if chrom not in genome_dict:
      continue
    for iv in chrom_intervals_dict[chrom]:
      iv_len = iv.end - iv.begin
      if iv_len > args.min_region:
        unmapped_seq = genome_dict[chrom].seq[iv.begin:iv.end]
        # if unmapped sequence is entirely made of gaps - skip
        if unmapped_seq.replace('N','') == '':
          continue
        unmapped_name = "%s%s_%s-%s" %(args.genome_name, chrom, iv.begin, iv.begin+len(unmapped_seq))
        unmapped_rec = SeqRecord(unmapped_seq, id=unmapped_name, description='')
        out_fasta_records.append(unmapped_rec)

# Write records to output fasta
SeqIO.write(out_fasta_records, args.out_fasta, "fasta")

# If provided - read gff and output genes found on output records

def interval_contains(interval_tree, start, end):
  res = interval_tree[start:end]
  if not res or len(res) > 1:
    return None
  iv = res.pop()
  if iv.begin <= start and iv.end >= end:
    return iv
  return None

def convert_feature_coords(feature, start, end):
  feature.seqid = "%s%s_%s-%s" %(args.genome_name, feature.seqid,start,end)
  feature.start -= start
  feature.end -= start
  return feature

if args.in_gff:
  # create gff DB
  db_path = "tmp.sqlite3"
  gff_db = gffutils.create_db(args.in_gff, db_path, force=True, merge_strategy="create_unique")
  gff = gffutils.FeatureDB(db_path)
  # extract genes on output features
  with open(args.out_gff,'w') as fo:
    print("##gff-version 3", file=fo)
    for gene in gff.features_of_type("gene"):
      seqid = gene.seqid
      gene_start = gene.start
      gene_end = gene.end
      seqid_intervals = chrom_intervals_dict[seqid]
      gene_iv = interval_contains(seqid_intervals, gene_start, gene_end)
      if gene_iv and gene_iv.end - gene_iv.begin > args.min_region:
        print(convert_feature_coords(gene, gene_iv.begin, gene_iv.end), file=fo)
        for f in gff.children(gene):
          print(convert_feature_coords(f, gene_iv.begin, gene_iv.end), file=fo)

  remove(db_path)
