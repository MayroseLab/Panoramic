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
parser.add_argument('--min_protein', default=0, type=int, help='Min length of novel proteins in output gff (determined as len(CDS)/3)')
parser.add_argument('--remap_genome', default=None, help='Indicates that this is a remap run of provided fasta')
args = parser.parse_args()
if args.genome_name:
  args.genome_name += "__"

def parse_region_name(name):
  # assumes sequence names like "genome1__chr1_1000-2000"
  genome, rest = name.split('__')
  start, end = [int(n) for n in rest.split('_')[-1].split('-')]
  chrom = '_'.join(rest.split('_')[:-1])
  return genome, chrom, start, end

# Read input fasta
if not args.remap_genome:
  genome_dict = SeqIO.to_dict(SeqIO.parse(args.in_fasta, "fasta"))
else:
  genome_dict = SeqIO.to_dict(SeqIO.parse(args.remap_genome, "fasta"))
# create per-chromosome intervals
# this is used for finding unmapped regions of the query genome
# initialize by creating full-length chromosome intervals
# and chop whenerver mapped regions are detected while parsing PAF
if not args.remap_genome:
  chrom_intervals_dict = { chrom: IntervalTree([Interval(0,len(genome_dict[chrom].seq)-1)]) for chrom in genome_dict }
else:
  chrom_intervals_dict = {}
  for rec in SeqIO.parse(args.in_fasta, 'fasta'):
    genome, chrom, start, end = parse_region_name(rec.id)
    if chrom not in chrom_intervals_dict:
      chrom_intervals_dict[chrom] = IntervalTree()
    chrom_intervals_dict[chrom][start:end] = 1

# Read PAF and find novel regions
with PafFile(args.in_paf) as paf:
  for record in paf:
    if args.remap_genome:
      genome, orig_chrom, orig_start, orig_end = parse_region_name(record.qname)
    else:
      orig_chrom = record.qname
      orig_start = 0
      orig_end = len(genome_dict[orig_chrom].seq)-1
      
    # if entire scaffold is unmapped
    if record.tname == "*":
      continue
    # if at least part of the record is mapped
    # chop mapped region from the chromosome interval
    # fetch CIGAR and use it to extract matched regions
    cigar = str(record.tags["cg"])[5:]
    cigar_tuples = re.findall(r'(\d+)([A-Z]{1})', cigar)
    start = record.qstart + orig_start
    for ct in cigar_tuples:
      if ct[1] in {"D","N"}:
        continue
      if ct[1] == "M":
        end = start + int(ct[0])
        chrom_intervals_dict[orig_chrom].chop(start,end)
      start += int(ct[0])

  # go over all unmapped regions, strip N's from the ends and update if needed
  # only keep intervals larger than cutoff
  filter_chrom_intervals_dict = {}
  # and keep fasta records
  out_fasta_records = []
  for chrom in chrom_intervals_dict:
    filter_chrom_intervals_dict[chrom] = IntervalTree()
    for iv in chrom_intervals_dict[chrom]:
      iv_len = iv.end - iv.begin
      if iv_len < args.min_region:
        continue
      # get unmapped sequence
      unmapped_seq = genome_dict[chrom].seq[iv.begin:iv.end]
      # if unmapped sequence is entirely made of gaps - skip
      if str(unmapped_seq).replace('N','') == '':
        continue
      # strip N's from the ends of sequence and update coordinates
      unmapped_seq = unmapped_seq.lstrip('N')
      strip_begin = iv.begin + (iv_len - len(unmapped_seq))
      iv_len = iv.end - strip_begin
      unmapped_seq = unmapped_seq.rstrip('N')
      strip_end = iv.end - (iv_len - len(unmapped_seq))
      strip_iv_len = strip_end - strip_begin
      if strip_iv_len >= args.min_region:
        unmapped_name = "%s%s_%s-%s" %(args.genome_name, chrom, strip_begin, strip_end)
        filter_chrom_intervals_dict[chrom][strip_begin:strip_end] = unmapped_name
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

def convert_feature_coords(feature, iv):
  feature.seqid = iv.data
  feature.start -= iv.begin
  feature.end -= iv.begin
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
      if seqid not in filter_chrom_intervals_dict:
        continue
      seqid_intervals = filter_chrom_intervals_dict[seqid]
      gene_iv = interval_contains(seqid_intervals, gene_start, gene_end)
      if gene_iv:
        # if min prot len was defined, filter short proteins
        if args.min_protein > 0:
          n_cds = 0
          total_cds_len = 0
          for cds in gff.children(gene, featuretype='CDS'):
            n_cds += 1
            total_cds_len += cds.end - cds.start
          if n_cds > 0 and total_cds_len/3 < args.min_protein:
            continue
        print(convert_feature_coords(gene, gene_iv), file=fo)
        for f in gff.children(gene):
          print(convert_feature_coords(f, gene_iv), file=fo)

  remove(db_path)
