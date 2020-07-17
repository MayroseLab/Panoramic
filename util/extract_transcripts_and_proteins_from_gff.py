"""
Given a gff file and the corresponding fasta,
extract, transcript, cDNA and protein sequences
into fasta files (based on exon and CDS features)
"""

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import gffutils
import sys
from time import time
from os import remove

in_gff = sys.argv[1]
in_fasta = sys.argv[2]
out_trans = sys.argv[3]
out_cdna = sys.argv[4]
out_prot = sys.argv[5]

# read genome fasta
genome_dict = SeqIO.to_dict(SeqIO.parse(in_fasta, "fasta"))
# read gff
db_path = "tmp_%s.sqlite3" % time()
gff_db = gffutils.create_db(in_gff, db_path, force=True, merge_strategy="create_unique")
gff = gffutils.FeatureDB(db_path)

# go over mRNA features and create transcripts and cDNA sequences
trans_recs = []
cdna_recs = []
prot_recs = []
for trans in gff.features_of_type("mRNA"):
  name = trans["ID"][0]
  chrom = trans.seqid
  trans_seq = ''
  for exon in gff.children(trans, featuretype='exon', order_by='start'):
    start = exon.start - 1
    end = exon.end
    exon_seq = genome_dict[chrom].seq[start:end]
    trans_seq += exon_seq
  if trans.strand == "-":
    trans_seq = trans_seq.reverse_complement()
  trans_rec = SeqRecord(trans_seq, id=name, description='')
  trans_recs.append(trans_rec)

  cdna_seq = ''
  for cds in gff.children(trans, featuretype='CDS', order_by='start'):
    start = cds.start - 1
    end = cds.end
    cds_seq = genome_dict[chrom].seq[start:end]
    cdna_seq += cds_seq
  if len(cdna_seq) > 0:
    if trans.strand == "-":
      cdna_seq = cdna_seq.reverse_complement()
    cdna_rec = SeqRecord(cdna_seq, id=name, description='')
    cdna_recs.append(cdna_rec)

    prot_seq = cdna_seq.translate()
    prot_rec = SeqRecord(prot_seq, id=name, description='')
    prot_recs.append(prot_rec)

# write to fasta files
SeqIO.write(trans_recs, out_trans, 'fasta')
SeqIO.write(cdna_recs, out_cdna, 'fasta')
SeqIO.write(prot_recs, out_prot, 'fasta')

remove(db_path)
