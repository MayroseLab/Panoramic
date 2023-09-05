import sys
import pandas as pd

data = pd.read_csv(sys.argv[1], sep="\t")
data[(data.annotation_gff == "--") & (data.proteins_fasta == "--")][['sample', 'genome_fasta']].to_csv(sys.argv[2], sep="\t", index=False)
data[~((data.annotation_gff == "--") & (data.proteins_fasta == "--"))].to_csv(sys.argv[3], sep="\t", index=False)

