"""
Given a multi-fasta file, create
a new fasta where records are
concatenated (with gaps) to create
records (chunks) of size up to
chunk_size.
If a record longer than chunk_size
is encountered - print as is (will
not split records).
Also outputs a BED file with positions
of the original records in the chunks.
"""

from __future__ import print_function
import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

in_fasta = sys.argv[1]
out_fasta = sys.argv[2]
out_bed = out_fasta + '.chunks.bed'
fo = open(out_fasta,'w')
chunk_size = int(sys.argv[3])


# Partition records into chunks with max chunk_size
gap_len = 100
chunks = {}
chunks_remain = {}
i = 1
for rec in SeqIO.parse(in_fasta, 'fasta'):
  rec_seq = rec.seq
  rec_seq_len = len(rec_seq)
  if rec_seq_len >= chunk_size:
    chunks[i] = [rec]
    i += 1
  else:
    found = False
    for s in chunks_remain:
      if chunks_remain[s] >= rec_seq_len + gap_len:
        chunks[s].append(rec)
        chunks_remain[s] -= rec_seq_len + gap_len
        found = True
        break
    if not found:
      chunks[i] = [rec]
      chunks_remain[i] = chunk_size - rec_seq_len
      i += 1
    
# print chunks to fasta and bed
print(i)
print(sorted(chunks.keys()))
with open(out_fasta,'w') as of, open(out_bed,'w') as ob:
  for c in range(1,i):
    records = chunks[c]
    chunk_id = 'chunk%s' % c
    # concatenate sequences
    chunk_seq = Seq(('N' * gap_len).join(str(rec.seq) for rec in records))
    chunk_rec = SeqRecord(chunk_seq, id=chunk_id, name='', description='')
    print(chunk_rec.format('fasta'), file=of)
    # print to bed
    start = 0
    for rec in records:
      orig_id = rec.id
      orig_len = len(rec.seq)
      end = start + orig_len
      print("%s\t%s\t%s\t%s" %(chunk_id, start, end, orig_id), file=ob)
      start = end + gap_len
