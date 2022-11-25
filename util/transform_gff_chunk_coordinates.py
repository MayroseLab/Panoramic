"""
Convert gff from chunks coordinates
to original contig coordinates, based
on a bed file.
* If a gene spans multiple original
contigs or is not contained within
a single contig, it will be discarded.
"""

from __future__ import print_function
import sys
from intervaltree import Interval, IntervalTree
import gffutils

in_gff = sys.argv[1]
in_bed = sys.argv[2]
out_gff = sys.argv[3]

# parse bed file
chunks = {}
with open(in_bed) as f:
  for line in f:
    chunk, start, end, name = line.strip().split('\t')
    start = int(start) + 1
    end = int(end)
    if chunk not in chunks:
      chunks[chunk] = IntervalTree()
    chunks[chunk].add(Interval(start,end,name))

db_path = in_gff + '.db'
try:
  gff_db = gffutils.create_db(in_gff, dbfn=db_path, force=True, merge_strategy='create_unique')
except ValueError:
  # input gff is empty - create empty output gff and exit
  with open(out_gff, 'w') as fo:
    print('##gff-version 3', file=fo)
  sys.exit(0)
  
gff_db = gffutils.FeatureDB(db_path, keep_order=True)

def gff_feature_to_list(feat):
  return str(feat).split('\t')

with open(out_gff,'w') as fo:
  print('##gff-version 3', file=fo)
  for gene in gff_db.features_of_type('gene'):
    res = chunks[gene.seqid][gene.start:gene.end]
    # if gene spans more than one interval or not completely within interval - skip
    if len(res) > 1:
      print(str(gene))
      print(res)
      continue
    interval = res.pop()
    if not (interval.begin <= gene.start and interval.end >= gene.end):
      print(str(gene))
      print(res)
      continue
   
    orig_seqid = interval.data
    
    gene_gff_fields = gff_feature_to_list(gene)
    gene_gff_fields[0] = orig_seqid
    gene_gff_fields[3] = str(gene.start - interval.begin + 1)
    gene_gff_fields[4] = str(gene.end - interval.begin + 1)
    print('\t'.join(gene_gff_fields), file=fo)
    ftype_ord = ['mRNA','CDS','exon']
    for ftype in ftype_ord:
      for feat in gff_db.children(gene, featuretype=ftype):
        gff_fields = gff_feature_to_list(feat)
        gff_fields[0] = orig_seqid
        gff_fields[3] = str(feat.start - interval.begin + 1)
        gff_fields[4] = str(feat.end - interval.begin + 1)
        print('\t'.join(gff_fields), file=fo)
