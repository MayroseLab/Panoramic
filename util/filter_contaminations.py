import sys
import pandas as pd
from Bio import SeqIO


def find_base_level(report):
    f = pd.read_csv(report, sep='\t', names=['precent', 'contigs_number', 'appears', 'identifier', 'num_identifier', 'class_name'])
    filtered_df = f[f.appears > 0]
    class_lst = ['unclassified', 'root', 'other sequences']
    base_level = filtered_df[~filtered_df['class_name'].str.strip().isin(class_lst)]['num_identifier'].tolist()
    return base_level


def find_contig_names(classification, base_level):
    f = pd.read_csv(classification, sep='\t', names=['classify', 'contig_id', 'class_info', 'len', 'other_info'])
    f['tax_id'] = f['class_info'].apply(lambda s: s.split("(taxid")[-1].strip(")").strip())
    table = f[(f['tax_id'].apply(lambda x: int(x) not in base_level)) | (f.classify == 'U')]
    return table.contig_id.tolist()


def filter_contaminations(contigs, output_fasta, filtered_contigs):
    with open(contigs, "r") as fasta_file, open(output_fasta, "w") as filtered_fasta:
        for record in SeqIO.parse(fasta_file, "fasta"):
            if record.id in filtered_contigs:
                SeqIO.write(record, filtered_fasta, "fasta")


def calc_contamination(contigs, output_fasta, result_dir):
    contigs_sum = sum([len(seq_record) for seq_record in SeqIO.parse(contigs, "fasta")])
    filtered = sum([len(seq_record) for seq_record in SeqIO.parse(output_fasta, "fasta")])
    precent = (1-(filtered / contigs_sum)) * 100
    with open(result_dir, "a") as f:
        f.write(str(round(precent, 4)))


if __name__ == '__main__':
    base_level = find_base_level(sys.argv[1])
    filtered_contigs = find_contig_names(sys.argv[2], base_level)
    filter_contaminations(sys.argv[3], sys.argv[4], filtered_contigs)
    calc_contamination(sys.argv[3], sys.argv[4], sys.argv[5])
