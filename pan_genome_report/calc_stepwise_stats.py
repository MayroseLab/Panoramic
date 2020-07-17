"""
Produce statistics relevant to the
genome step-wise addition plots.
This requires sampling from the
full list of samples multiple times
to produce stdev information.
"""

import sys
import pandas as pd
from random import sample

in_pav_table = sys.argv[1]
n_iter = int(sys.argv[2])
out_tsv = sys.argv[3]

# read PAV TSV
pav_df = pd.read_csv(in_pav_table, sep='\t', index_col = 0)

# generate random samples orders
all_samples = list(pav_df.columns)
n_samples = len(all_samples)
orders = [ sample(all_samples, n_samples) for i in range(n_iter) ]

# calculate stats per order
def calc_stats(pav_df, samples):
  """
  Given the PAV df and a list
  of samples (column names),
  calculate pan genome stats
  """
  df = pav_df[samples]
  n_samples = len(samples)
  row_sums = df.sum(axis=1)
  tot_genes = (row_sums > 0).sum()
  core_genes = (row_sums == n_samples).sum()
  singletons = (row_sums == 1).sum()
  return tot_genes, core_genes, singletons

rows = []
ord_id = 1
for curr_ord in orders:
  for s in range(1,n_samples+1):
    samples_list = curr_ord[:s]
    tot_genes, core_genes, singletons = calc_stats(pav_df, samples_list)
    row1 = pd.Series([ord_id, s, "Total genes", tot_genes])
    row2 = pd.Series([ord_id, s, "Core genes", core_genes])
    row3 = pd.Series([ord_id, s, "Singletons", singletons])
    rows.extend([row1,row2,row3])
  ord_id += 1

stepwise_df = pd.DataFrame(rows)
stepwise_df.columns = ["ord_id", "n_samples", "stat", "n_genes"]
stepwise_df.to_csv(out_tsv, sep='\t', index=False)
