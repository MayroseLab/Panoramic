"""
This script iteratively maps chromosome-
or contig-level assemblies to a reference
genome, while adding novel sequences to
the reference at each iteration, thereby
creating a pan genome.
The input is a TSV file with header line:
sample	genome_fasta  annotation_gff (optional)
Genes on novel sequences will be extracted from
the gff3 files.
*** The order of genomes in the input TSV may change the results ***
"""

from __future__ import print_function
import sys
import os
import csv
import argparse
from Bio import SeqIO

if __name__ == "__main__":

  # command line args
  parser = argparse.ArgumentParser()
  parser.add_argument('ref_fasta', help='Path to reference genome fasta')
  parser.add_argument('ref_gff', help='Path to reference annotation gff')
  parser.add_argument('ref_proteins', help='Path to reference proteins fasta')
  parser.add_argument('in_tsv', help='Path to TSV file with samples data')
  parser.add_argument('out_dir', help='Path to output directory')
  parser.add_argument('--cpus', default=1, type=int, help='Number of CPUs to use in Minimap')
  parser.add_argument('--min_len', default=1, type=int, help='Minimal sequence length to add to pan')
  parser.add_argument('--min_protein', default=0, type=int, help='Minimal protein length of genes to be included in output gff')
  args = parser.parse_args()

  # make sure required scripts are available
  script_dir = os.path.dirname(sys.argv[0])
  extract_script = os.path.join(script_dir, "extract_novel_regions_from_paf.py")
  filter_fasta_script = os.path.join(os.path.dirname(script_dir),"util/filter_fasta_by_gff.py")
  for script in [extract_script, filter_fasta_script]:  
    if not os.path.isfile(script):
      exit("Required script %s file is missing" % script)

  # re-write ref fasta so record
  # names only include simple IDs
  ref_recs = SeqIO.parse(args.ref_fasta, "fasta")
  simp_records = []
  for rec in ref_recs:
    rec.description = ''
    simp_records.append(rec)
  ref_simp_fasta = os.path.join(args.out_dir, "ref_genome.fasta")
  SeqIO.write(simp_records, ref_simp_fasta, "fasta")

  # initialize
  pan_fasta = ref_simp_fasta
  pan_gff = args.ref_gff
  pan_proteins = args.ref_proteins

  # read TSV
  with open(args.in_tsv) as f:
    reader = csv.DictReader(f, delimiter = '\t')
    for row in reader:
      genome_name = row['sample']
      assembly_fasta = row['genome_fasta']
      genome_gff = None
      proteins_fasta = None
      if 'annotation_gff' in row and row['annotation_gff']:
        genome_gff = row['annotation_gff']
        if 'proteins_fasta' in row and row['proteins_fasta']:
          proteins_fasta = row['proteins_fasta']

      # run minimap2
      out_paf = os.path.join(args.out_dir, "%s_vs_pan.paf" % genome_name)
      os.system("minimap2 -x asm5 -L -t %s %s %s -o %s -c" % (args.cpus, pan_fasta, assembly_fasta, out_paf))

      # extract novel sequences
      # first run
      novel_seq_fasta = os.path.join(args.out_dir, "%s_novel.fasta" % genome_name)
      os.system("python %s %s %s %s %s --genome_name %s" %(extract_script, out_paf, assembly_fasta, args.min_len, novel_seq_fasta, genome_name))
      # validate by remapping extracted sequences
      out_paf_remap = os.path.join(args.out_dir, "%s_vs_pan_remap.paf" % genome_name)
      os.system("minimap2 -x asm5 -L -t %s %s %s -o %s -c" % (args.cpus, pan_fasta, novel_seq_fasta, out_paf_remap))
      novel_seq_fasta_remap = os.path.join(args.out_dir, "%s_novel_remap.fasta" % genome_name)
      # extract from remap (use gff if provided)
      max_dist = int(args.min_len/5)
      if genome_gff:
        novel_gff = os.path.join(args.out_dir, "%s_novel.gff" % genome_name)
        os.system("python %s %s %s %s %s --in_gff %s --out_gff %s --genome_name %s --min_protein %s --remap_genome %s --max_merge_dist %s" %(extract_script, out_paf_remap, novel_seq_fasta, args.min_len, novel_seq_fasta_remap, genome_gff, novel_gff, genome_name, args.min_protein, assembly_fasta, max_dist))
      else:
        os.system("python %s %s %s %s %s --genome_name %s --remap_genome %s --max_merge_dist %s" %(extract_script, out_paf_remap, novel_seq_fasta, args.min_len, novel_seq_fasta_remap, genome_name, assembly_fasta, max_dist))

      # get novel proteins (if available)
      if proteins_fasta:
        novel_proteins = os.path.join(args.out_dir, "%s_novel_proteins.fasta" % genome_name)
        os.system("python %s %s %s %s mRNA ID" %(filter_fasta_script, novel_gff, proteins_fasta, novel_proteins))

      # create new pan genome
      curr_pan = os.path.join(args.out_dir, "current_pan.fasta")
      os.system("cat %s %s > %s" %(pan_fasta, novel_seq_fasta_remap, curr_pan + '.tmp'))
      os.replace(curr_pan + '.tmp', curr_pan)

      # if gff exists, add to pan genome annotation
      if genome_gff:
        curr_gff = os.path.join(args.out_dir, "current_pan.gff")
        os.system("cat %s %s > %s" %(pan_gff, novel_gff, curr_gff + '.tmp'))
        os.replace(curr_gff + '.tmp', curr_gff)
        # if proteins fasta exist, add to pan proteome
        if proteins_fasta:
          curr_prot = os.path.join(args.out_dir, "current_pan_prot.fasta")
          os.system("cat %s %s > %s" %(pan_proteins, novel_proteins, curr_prot + '.tmp'))
          os.replace(curr_prot + '.tmp', curr_prot)

      # update pan
      pan_fasta = curr_pan
      if genome_gff:
        pan_gff = curr_gff
        if proteins_fasta:
          pan_proteins = curr_prot

  # create final pan genome and gff
  final_pan = os.path.join(args.out_dir, "pan_genome.fasta")
  os.link(pan_fasta, final_pan)
  final_gff = os.path.join(args.out_dir, "pan_genes.gff")
  os.link(pan_gff, final_gff)
  final_prot = os.path.join(args.out_dir, "pan_proteins.fasta")
  os.link(pan_proteins, final_prot)
  # create final novel fasta
  all_novel_fasta = os.path.join(args.out_dir, "all_novel.fasta")
  os.system("cat %s/*_novel_remap.fasta > %s" %(args.out_dir, all_novel_fasta))
