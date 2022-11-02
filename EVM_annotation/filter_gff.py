import sys
import gffutils

in_gff = sys.argv[1]
max_aed = float(sys.argv[2])
min_prot = int(sys.argv[3])

# create GFF DB
gff_db_path = in_gff + '.db'
gff_db = gffutils.create_db(in_gff, dbfn=gff_db_path, force=True, merge_strategy='create_unique')
gff_db = gffutils.FeatureDB(gff_db_path)

print('##gff-version 3')
for mrna in gff_db.features_of_type('mRNA', order_by='start'):
  # if AED or wAED > max_aed - discard
  if float(mrna['AED'][0]) > max_aed or float(mrna['wAED'][0]) > max_aed:
    continue
  # calculate protein length
  prot_len = sum([cds.end - cds.start + 1 for cds in gff_db.children(mrna, featuretype='CDS')])/3
  # if protein length < min_prot - discard
  if prot_len < min_prot:
    continue
  # if mRNA passed filtration, find parent gene and print all mRNA children
  parent_gene = next(gff_db.parents(mrna, featuretype='gene'))
  print(str(parent_gene))
  print(str(mrna))
  for cf in gff_db.children(mrna):
    print(str(cf))
