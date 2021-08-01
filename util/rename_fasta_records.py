import sys
from Bio import SeqIO

in_fasta = sys.argv[1]
in_id_table = sys.argv[2]
out_fasta = sys.argv[3]

# parse table
with open(in_id_table) as f:
  id_map_dict = dict([line.strip().split('\t')[:2] for line in f])

# convert record headers
with open(out_fasta,'w') as fo:
  for rec in SeqIO.parse(in_fasta, 'fasta'):
    orig_id = rec.id
    new_id = id_map_dict[orig_id]
    rec.id = new_id
    rec.name = ''
    rec.description = ''
    print(rec.format('fasta'), file=fo)
