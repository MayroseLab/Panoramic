"""
Parses the output of OrthoFinder
and creates PAV and CNV matrices.
The gene names are selected as follows:
If an OG contains a reference gene, it
is given as the gene name. If multiple
ref genes exist - take the first. If
no ref genes exist, call it PanGeneX,
where X is a running ID.
PAV and CNV are the same except in PAV
matrix, CNs > 1 are converted to 1.
"""

import pandas as pd
import sys

in_orthogroups_tsv = sys.argv[1]
ref_name = sys.argv[2]
out_pav_tsv = sys.argv[3]
out_cnv_tsv = sys.argv[4]
out_names_mapping = sys.argv[5]

def choose_gene_name(row):
  if pd.notna(row[ref_name]):
    return row[ref_name].split(', ')[0]
  else:
    return "PanGene%s" %(row.name - first_pan_index + 1)

df = pd.read_csv(in_orthogroups_tsv, sep='\t', index_col='Orthogroup')
# remove _LQ, _HQ, _REF from genome names
suffixes = ['_LQ', '_HQ', '_REF']
new_cols = []
#for col in df.columns:
#  for suff in suffixes:
#    if col.endswith(suff):
#      new_cols.append(col.replace(suff, ''))
#df.columns = new_cols

# give each orthogroup a gene name
df = df.rename(columns={ref_name+"_REF": ref_name})
df.sort_values(by = ref_name, inplace=True)
df = df.reset_index()
first_pan_index = df.loc[df[ref_name].isna()].iloc[0].name
df['gene'] = df.apply(choose_gene_name, axis=1)
# write out name mapping
df.to_csv(out_names_mapping, sep='\t', columns=['Orthogroup','gene'], index=False)
# then keep manipulating table
df = df.drop('Orthogroup', axis=1)
df = df.set_index('gene')
# convert OG content to PAV (or CNV)
for col in df:
  df[col] = df[col].str.split(', ').str.len()
df = df.fillna(0)
df = df.astype(int)
# print CNV matrix
df.to_csv(out_cnv_tsv, sep='\t')
# convert to PAV and print
df[df > 1] = 1
df.to_csv(out_pav_tsv, sep='\t')
