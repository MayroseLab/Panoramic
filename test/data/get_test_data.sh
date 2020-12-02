#!/bin/bash

# Reference (S288C) 
wget http://sgd-archive.yeastgenome.org/sequence/S288C_reference/genome_releases/S288C_reference_genome_R64-2-1_20150113.tgz
tar -zxvf S288C_reference_genome_R64-2-1_20150113.tgz
mv S288C_reference_genome_R64-2-1_20150113/{S288C_reference_sequence_R64-2-1_20150113.fsa,saccharomyces_cerevisiae_R64-2-1_20150113.gff,orf_coding_all_R64-2-1_20150113.fasta,orf_trans_all_R64-2-1_20150113.fasta} ./
rm -rf S288C_reference_genome_R64-2-1_20150113.tgz S288C_reference_genome_R64-2-1_20150113/
sed -i -e 's/>.*chromosome=\([^]]*\)]/>chr\1/' -e 's/>ref|NC_001224|.*/>mt/' S288C_reference_sequence_R64-2-1_20150113.fsa
head -23076 saccharomyces_cerevisiae_R64-2-1_20150113.gff > saccharomyces_cerevisiae_R64-2-1_20150113.gff.tmp
mv saccharomyces_cerevisiae_R64-2-1_20150113.gff.tmp saccharomyces_cerevisiae_R64-2-1_20150113.gff
sed -i 's/ .*/_mRNA/' orf_coding_all_R64-2-1_20150113.fasta
sed -i 's/ .*/_mRNA/' orf_trans_all_R64-2-1_20150113.fasta

# RM-11a
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/149/365/GCA_000149365.1_ASM14936v1/GCA_000149365.1_ASM14936v1_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/149/365/GCA_000149365.1_ASM14936v1/GCA_000149365.1_ASM14936v1_genomic.gff.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/149/365/GCA_000149365.1_ASM14936v1/GCA_000149365.1_ASM14936v1_protein.faa.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/149/365/GCA_000149365.1_ASM14936v1/GCA_000149365.1_ASM14936v1_rna_from_genomic.fna.gz
gzip -d *.gz
sed -i 's/>.*locus_tag=\([^]]*\)].*/>rna-\1/' GCA_000149365.1_ASM14936v1_rna_from_genomic.fna
python fix_ids.py > GCA_000149365.1_ASM14936v1_protein.faa.fix
mv GCA_000149365.1_ASM14936v1_protein.faa.fix GCA_000149365.1_ASM14936v1_protein.faa

# T73
wget http://sgd-archive.yeastgenome.org/sequence/strains/T73/T73_WashU_2011_AFDF01000000/T73_WashU_2011_AFDF01000000.fsa.gz
wget http://sgd-archive.yeastgenome.org/sequence/strains/T73/T73_WashU_2011_AFDF01000000/T73_AFDF01000000_pep.fsa.gz
wget http://sgd-archive.yeastgenome.org/sequence/strains/T73/T73_WashU_2011_AFDF01000000/T73_AFDF01000000_cds.fsa.gz
wget http://sgd-archive.yeastgenome.org/sequence/strains/T73/T73_WashU_2011_AFDF01000000/T73_AFDF01000000.gff.gz
gzip -d *.gz
head -13983 T73_AFDF01000000.gff | awk '$3 != "contig"' | awk '$3 == "gene" {split($9,a,";"); split(a[1],b,"="); print $0; print $1"\t"$2"\tmRNA\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\tID="b[2]"_mRNA;Parent="b[2]; print $1"\t"$2"\texon\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\tID="b[2]"_exon;Parent="b[2]"_mRNA"; print $1"\t"$2"\tCDS\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\tID="b[2]"_CDS;Parent="b[2]"_mRNA"}' > T73_AFDF01000000.gff.tmp
mv T73_AFDF01000000.gff.tmp T73_AFDF01000000.gff
sed -i 's/ .*/_mRNA/' T73_AFDF01000000_pep.fsa
sed -i 's/ .*/_mRNA/' T73_AFDF01000000_cds.fsa

# dummy repeats
echo -e ">dummy\nAGAG" > repeats.fa
