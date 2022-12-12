import sys
import gffutils

in_gff = sys.argv[1]
out_db = sys.argv[2]

gffutils.create_db(in_gff, dbfn=out_db, force=True, keep_order=True, merge_strategy='create_unique', sort_attribute_values=True)
