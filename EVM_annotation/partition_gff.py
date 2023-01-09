"""
Given a gff3 file and the EVM partitions list,
split the gff to multiple partitions.
"""

import sys
import os
import gffutils

in_gff = sys.argv[1]
partitions_list = sys.argv[2]

def write_features_to_gff(features_list, out_gff):
  with open(out_gff,'w') as fo:
    print('##gff-version 3', file=fo)
    for feat in features_list:
      print(str(feat), file=fo)

def convert_coords(feat, partition_start):
  """
  Convert start and end coordinates
  of GFF feature so they refer to
  the partition coordinates (starting
  from 1)
  """
  feat.start = feat.start - partition_start + 1
  feat.end = feat.end - partition_start + 1
  return feat

# parse partitions list
partitions_dict = {}
with open(partitions_list) as f:
  for line in f:
    partition_path = line.strip().split('\t')[3]
    partition_base = os.path.basename(partition_path)
    chrom, coords = partition_base.rsplit('_',1)
    start, end = [int(x) for x in coords.split('-')]
    if chrom not in partitions_dict:
      partitions_dict[chrom] = {}
    partitions_dict[chrom][(start, end)] = partition_path

# build gff DB
gff_db_path = in_gff + '.db'
empty_input = False
try:
  gff_db = gffutils.create_db(in_gff, dbfn=gff_db_path, force=True, merge_strategy='create_unique')
  gff_db = gffutils.FeatureDB(gff_db_path)
except ValueError:
  empty_input = True

# create partitions
gff_basename = os.path.basename(in_gff)
for chrom in sorted(partitions_dict.keys()):
  for p in partitions_dict[chrom]:
    start, end = p
    if empty_input:
      features_in_region = []
      parents_in_region = []
    else:
      features_in_region = list(gff_db.region((chrom,start,end), completely_within=True))
      parents_in_region = set([feat['ID'][0] for feat in features_in_region if 'ID' in feat.attributes])
    # remove features whos parent is not in the list - repeat until no more paretless features are found
    parentless_count = 1
    while parentless_count > 0:
      features_in_region_new = list(filter(lambda feat: ('Parent' not in feat.attributes) or (feat['Parent'][0] in parents_in_region), features_in_region))
      parentless_count = len(features_in_region) - len(features_in_region_new)
      features_in_region = features_in_region_new
      parents_in_region = set([feat['ID'][0] for feat in features_in_region if 'ID' in feat.attributes])
    # convert coordinates
    features_in_region = [convert_coords(feat, start) for feat in features_in_region]
    # write features
    partition_dir = partitions_dict[chrom][p]
    if not os.path.isdir(partition_dir):
      os.makedirs(partition_dir)
    region_gff = os.path.join(partition_dir, gff_basename)
    write_features_to_gff(features_in_region, region_gff)
