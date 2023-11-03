"""
The code reads high-quality (HQ) data from a TSV file, filters it based on HQ samples that are fully annotated, and saves the fully annotated HQ samples into a new TSV file.
"""

import sys
import pandas as pd

data = pd.read_csv(sys.argv[1], sep="\t")
data[~((data.annotation_gff == "--") & (data.proteins_fasta == "--"))].to_csv(sys.argv[2], sep="\t", index=False)
