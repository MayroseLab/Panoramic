"""
Combine results from multiple SGSGeneLoss
.excov files into one PAV matrix:
	sample_1|sample_2|...|sample_n
gene_1	1       |1       |...|0
gene_2	0       |1       |...|1
...
gene_m	1       |0       |...|1

The script takes a list of .excov files
(one per sample) sepearated by spaces.
File names should be:
<sample name>.all.PAV
after the list there are 4 additional parameters:
1. name of the reference genome
2. List of ref gene names (to create ref PAV column)
3. a tsv with two columns: gene name and new name,
   used for modifying gene names
4. output PAV tsv
"""

import pandas as pd
import sys
import os

in_files = sys.argv[1:-4]
ref_name = sys.argv[-4]
ref_genes_list = sys.argv[-3]
name_sub_tsv = sys.argv[-2]
out_tsv = sys.argv[-1]


samples_pav = []
for f in in_files:
  df = pd.read_csv(f, sep='\t')
  pav_s = df.iloc[:,1]
  pav_s.index = df['Gene']
  samples_pav.append(pav_s)

pav_df = pd.concat(samples_pav, axis=1)
ref_genes_df = pd.read_csv(ref_genes_list, names=["Gene"], index_col = 0, usecols = [0])
ref_genes_df[ref_name] = 1
pav_df[ref_name] = pd.concat([pav_df,ref_genes_df], axis=1).iloc[:,-1].fillna(0)
name_sub_df = pd.read_csv(name_sub_tsv, sep="\t", index_col=0, names=['new_name'])
pav_df.index = pav_df.apply(lambda row: name_sub_df.loc[row.name]['new_name'] if row.name in name_sub_df.index else row.name, axis=1)
pav_df.index.rename('gene', inplace=True)
pav_df = pav_df.astype(int)
pav_df.to_csv(out_tsv, sep='\t')
