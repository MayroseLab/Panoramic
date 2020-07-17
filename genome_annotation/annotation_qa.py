from __future__ import print_function, division
import argparse
from gff3 import Gff3
import pandas as pd
from interval_map import intervalmap


### FUNCTIONS

def mrna_chrom(gff):
    """
    Takes a gff3 object and returns a dictionary
    of chromosomes with mRNA names as keys.
    """
    res = {}
    mrna_lines = [line for line in gff.lines if line['line_type'] == 'feature' and line['type'] == 'mRNA']
    for mrna_line in mrna_lines:
        gene_name = mrna_line['attributes']['Name']
        mrna_chrom = mrna_line['seqid']
        res[gene_name] = mrna_chrom
    return res

def mrna_exons_count(gff):
    """
    Takes a gff3 object and returns a dictionary
    of exons counts with mRNA names as keys.
    """
    res = {}
    mrna_lines = [line for line in gff.lines if line['line_type'] == 'feature' and line['type'] == 'mRNA']
    for mrna_line in mrna_lines:
        gene_name = mrna_line['attributes']['Name']
        mrna_qi = mrna_line['attributes']['_QI']
        n_exons = mrna_qi.split('|')[6]
        res[gene_name] = n_exons
    return res

def mrna_aed(gff):
    """
    Takes a gff3 object and returns a dictionary
    of AED scores with mRNA names as keys.
    """
    res = {}
    mrna_lines = [line for line in gff.lines if line['line_type'] == 'feature' and line['type'] == 'mRNA']
    for mrna_line in mrna_lines:
        gene_name = mrna_line['attributes']['Name']
        mrna_aed = mrna_line['attributes']['_AED']
        res[gene_name] = mrna_aed
    return res

def mrna_utr(gff):
    """
    Takes a gff3 object and returns a dictionary.
    Keys are mRNA IDs and values are UTR status
    -1 - illegal UTR; 0 - no UTRs; 1 - missing 3' or 5' UTR; 2 - legal UTRs
    """
    res = {}
    translate = {'five_prime_UTR':'5', 'three_prime_UTR':'3', 'CDS': 'C'}
    mrna_lines = [line for line in gff.lines if line['line_type'] == 'feature' and line['type'] == 'mRNA']
    for mrna_line in mrna_lines:
        gene_name = mrna_line['attributes']['Name']
        descendants = gff.descendants(mrna_line)
        descendants_features_order = [f['type'] for f in descendants if f['type'] in ['five_prime_UTR', 'CDS', 'three_prime_UTR']]
        descendants_features_order_simp = ''.join([translate[t] for t in descendants_features_order])
        if '5' not in descendants_features_order_simp and '3' not in descendants_features_order_simp:
            stat = 0
        elif 'C5C' in descendants_features_order_simp or 'C3C' in descendants_features_order_simp:
            stat = -1
        elif '5' not in descendants_features_order_simp or '3' not in descendants_features_order_simp:
            stat = 1
        else:
            stat = 2
        res[gene_name] = stat
    return res

def prot_busco(busco_full):
    """
    Reads a BUSCO full report and returns a dict
    with proteins in which BUSCOs were found. {gene name: BUSCO ID}
    """
    res = {}
    with open(busco_full) as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.split('\t')
            if fields[1] == "Complete" or fields[1] == "Duplicated":
                res[fields[2]] = fields[0]
    return res

def prot_similarity(blast_res):
    """
    Reads a blast tsv result and returns a dict with
    gene names as keys and query coverage as values.
    Assumes query name at column 2 and query coverage at column 10.
    """
    res = {}
    with open(blast_res) as f:
        for line in f:
            fields = line.strip().split('\t')
            res[fields[0]] = fields[9]
    return res

def prot_domains(ips_result):
    """
    Reads the tsv output of an InterProScan run and
    returns a dict with gene name and minimal E-value found
    for any domain for any application for the protein.
    """
    res = {}
    with open(ips_result) as f:
        for line in f:
            fields = line.split('\t')
            gene, e_value = fields[0], fields[8]
            e_value = float(e_value)
            if gene not in res:
                res[gene] = e_value
            else:
                res[gene] = min(res[gene], e_value)
        return res

def repeats_overlap(genes_gff, repeats_gff):
    """
    Parses gene models and repeats GFF file and 
    returns the % of overlap per gene (dict).
    """
    # find repeats - store as interval maps
    repeats_dict = {}
    rep_gff_lines = list(repeats_gff.lines)
    repeat_lines = [line for line in rep_gff_lines if line['line_type'] == 'feature' and line['source'] == "repeatmasker" and line['type'] == 'match']
    for line in repeat_lines:
        chrom, start, end = line['seqid'], line['start'], line['end']
        if chrom not in repeats_dict:
           repeats_dict[chrom] = intervalmap()
        repeats_dict[chrom][start:end] = 'repeat'
    # find mrna overlaps
    res = {}
    genes_gff_lines = list(genes_gff.lines)
    mrna_lines = [line for line in genes_gff_lines if line['line_type'] == 'feature' and line['type'] == 'mRNA']
    for mrna_line in mrna_lines:
        gene_name = mrna_line['attributes']['Name']
        chrom, mrna_start, mrna_end = mrna_line['seqid'], mrna_line['start'], mrna_line['end']
        chrom_repeats_interval_map = repeats_dict[chrom]
        mrna_len = mrna_end - mrna_start
        repeat_overlap_len = len(chrom_repeats_interval_map.slice(mrna_start,mrna_end))
        res[gene_name] = repeat_overlap_len/mrna_len*100
    return res


### MAIN
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('genes_gff', help="Input gff3 annotation file including only genes features")
    parser.add_argument('--busco_result', default=None, help="BUSCO full report output")
    parser.add_argument('--ips_result', help="InterProScan tsv output")
    parser.add_argument('--blast_result', help="Blast search result of proteins vs. DB")
    parser.add_argument('--repeats_gff', help="gff3 file with repeats features only")
    parser.add_argument('out_report', help="Output QA report")
    args = parser.parse_args()

    gff_obj = Gff3(args.genes_gff)
    qa_methods = [("Chromosome",mrna_chrom,(gff_obj,)), ("AED",mrna_aed,(gff_obj,)), ("Exons",mrna_exons_count,(gff_obj,)), ("UTR",mrna_utr,(gff_obj,))]
    if args.repeats_gff:
        rep_gff_obj = Gff3(args.repeats_gff)
        qa_methods.append(('Repeats',repeats_overlap,(gff_obj, rep_gff_obj)))
    if args.busco_result:
        qa_methods.append(('BUSCO',prot_busco,(args.busco_result,)))
    if args.blast_result:
        qa_methods.append(('BLAST',prot_similarity,(args.blast_result,)))
    if args.ips_result:
        qa_methods.append(('IPS',prot_domains,(args.ips_result,)))

    qa_data = [ pd.Series(m[1].__call__(*m[2]), name=m[0]) for m in qa_methods ]
    qa_df = pd.concat(qa_data, axis = 1)
    qa_df.to_csv(args.out_report, sep='\t', index_label="gene", na_rep='NA')
