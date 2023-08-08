"""
Input tsv should contain lines with the following values:
<sample name>	<quast report tsv> <busco short sum>	[ragoo fasta]
"""

from __future__ import division, print_function
import sys
import pandas as pd
from Bio import SeqIO


in_tsv = sys.argv[1]
out_tsv = sys.argv[2]

columns = []
with open(in_tsv) as f:
  for line in f:
    if line == '' or line == '\n' or line.startswith('#'):
      continue
    fields = line.strip().split('\t')
    sample, quast_report, busco_sum, ragoo_fasta, data_stats_file, contamination = fields
    # get data stats
    stats_df = pd.read_csv(data_stats_file, sep='\t', index_col = 0, names=[sample])
    # get quast stats
    quast_df = pd.read_csv(quast_report, sep='\t', index_col = 0, names=[sample])
    stats_df = stats_df.append(quast_df)
    # get % complete BUSCOs
    with open(busco_sum) as f:
      for line in f:
        if line.startswith('\tC:'):
          complete = line.split('[')[0][3:-1]
    s = pd.Series([complete], name="% Complete BUSCOs", index=[sample])
    stats_df = stats_df.append(s)
    # calculate % unmapped from RG assembly (if exists
    if ragoo_fasta:
      total_len = 0
      chr0_len = 0
      for rec in SeqIO.parse(ragoo_fasta, 'fasta'):
        if rec.id == "Chr0_RagTag":
          chr0_len = len(rec.seq)
          total_len += len(rec.seq)
        else:
          total_len += len(rec.seq)
      unmapped_perc = chr0_len/total_len*100
      s = pd.Series([unmapped_perc], name="% unmapped (Chr0)", index=[sample])
      stats_df = stats_df.append(s)
    # QUAST report HTML link
    html_quast_report = "file://" + quast_report.replace('.tsv', '.html')
    s = pd.Series([html_quast_report], name="QUAST report", index=[sample])
    stats_df = stats_df.append(s)
    columns.append(stats_df)
    with open(contamination) as con_file:
      columns.append(con_file.readlines()[0])


if not columns:
  headers = ["Assembly","# contigs (>= 0 bp)","# contigs (>= 1000 bp)","# contigs (>= 5000 bp)",
           "# contigs (>= 10000 bp) # contigs (>= 25000 bp)","# contigs (>= 50000 bp)",
           "Total length (>= 0 bp)","Total length (>= 1000 bp)","Total length (>= 5000 bp)",
           "Total length (>= 10000 bp)","Total length (>= 25000 bp)","Total length (>= 50000 bp)",
           "# contigs","Largest contig","Total length","GC (%%)","N50","N75","L50","L75",
           "# total reads","# left","# right Mapped (%%)","Properly paired (%%)","Avg. coverage depth",
           "Coverage >= 1x (%%)","# N's per 100 kbp","%% Complete BUSCOs","%% unmapped (Chr0)",
           "QUAST report","Read length (bp)", "% Removed contamination"]
  print('\t'.join(headers), file=out_tsv)
else:
  stats_df = pd.concat(columns, axis=1)
  stats_df = stats_df.transpose()
  stats_df.to_csv(out_tsv, sep='\t')

