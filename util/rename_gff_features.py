"""
Given a gff3 file, renames all
genes (and child features) to:
<PREFIX><i> where <PREFIX> is given
by the user and i is a running index.
"""

import sys
import gffutils

in_gff = sys.argv[1]
out_gff = sys.argv[2]
pref = sys.argv[3]
id_map_out = out_gff + '.id_map.tsv'

db_path = in_gff + '.db'
try:
  gff_db = gffutils.create_db(in_gff, dbfn=db_path, force=True, merge_strategy='create_unique')
except ValueError:
  # input gff is empty - create empty output gff and table, then exit
  with open(out_gff,'w') as go, open(id_map_out,'w') as mo:
    print('##gff-version 3', file=go)
  sys.exit(0)

gff_db = gffutils.FeatureDB(db_path, keep_order=True)

i = 1
with open(out_gff,'w') as go, open(id_map_out,'w') as mo:
  print('##gff-version 3', file=go)
  for gene in gff_db.features_of_type('gene'):
    gene_id = pref + str(i)
    i += 1
    orig_id = gene['ID'][0]
    gene['ID'] = [gene_id]
    print(str(gene), file=go)
    print('%s\t%s' %(orig_id,gene_id), file=mo)
    mrna_i = 1
    for mrna in gff_db.children(gene, featuretype='mRNA'):
      mrna['Parent'] = [gene_id]
      orig_mrna_id = mrna['ID'][0]
      mrna_id = gene_id + '_mRNA%s' % mrna_i
      mrna['ID'] = [mrna_id]
      print(str(mrna), file=go)
      print('%s\t%s' %(orig_mrna_id,mrna_id), file=mo)
      mrna_i += 1
      exon_i = 1
      for exon in gff_db.children(mrna, featuretype='exon'):
        exon['Parent'] = [mrna_id]
        exon_id = mrna_id + '_exon%s' % exon_i
        exon['ID'] = [exon_id]
        print(str(exon), file=go)
        exon_i += 1
      cds_i = 1
      for cds in gff_db.children(mrna, featuretype='CDS'):
        cds['Parent'] = [mrna_id]
        cds_id = mrna_id + '_CDS%s' % cds_i
        cds['ID'] = [cds_id]
        print(str(cds), file=go)
        cds_i += 1
