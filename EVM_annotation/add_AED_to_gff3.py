"""
This script calculates Annotation Edit Distances
(AED) and weighted AED (wAED) for mRNA features
predicted by EVidencdModeler. These scores can
be used as per-gene quality measures, based on
the support of evidence for gene models.
The idea of AED is further explain in
Holt, C., & Yandell, M. (2011).
The output of the script is a GFF3 file, identical
to the input GFF3, except all mRNA features have
two additional attributes: AED (raw distance), and
wAED (considering the weight assigned to each
source of evidence). Use:
python add_AED_to_gff3.py -h
for usage instructions.
Requirements: python3, pandas, gffutils, intervaltree
"""

# IMPORTS
import gffutils
import os
import sys
import pandas as pd
import argparse
from intervaltree import Interval, IntervalTree

# FUNCTIONS

def create_gff_db(gff_path, force=False):
    print('Reading GFF %s' % gff_path)
    db_path = gff_path + '.db'
    if os.path.isfile(db_path) and not force:
        print('Using existing DB %s' % db_path)
    else:
        gff_db = gffutils.create_db(gff_path, dbfn=db_path, force=True, keep_order=True, merge_strategy='create_unique', verbose=True)
    return gffutils.FeatureDB(db_path, keep_order=True)

def gff_to_intervals(gff_path, feature_types=None, merge_overlaps=False):
    """
    Reads a GFF3 file and creates a data structure:
    {chr1: IntervalTree(), chr2: IntervalTree(),...}
    if feature_types (iterable) is given, only
    include features of the given types.
    """
    print('Reading GFF %s' % gff_path)
    res = {}
    with open(gff_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            fields = line.split('\t')
            chrom = fields[0]
            source = fields[1]
            ftype = fields[2]
            start = int(fields[3])
            end = int(fields[4])
            if source not in feature_types or (feature_types is not None and ftype not in feature_types[source]):
                continue
            if start == end:	# not allowing empty intervals
                continue
            if chrom not in res:
                res[chrom] = IntervalTree()
            res[chrom].add(Interval(start, end, source))
    if merge_overlaps:
        for chrom in res:
            res[chrom].merge_overlaps(data_reducer=lambda x, y: x)
    return res

def merge_gff_iv_trees(gff_ds_list):
    """
    Takes a list of data structures created
    by gff_to_intervals and merges them into
    one structure of the same type.
    """
    res = {}
    for ds in gff_ds_list:
        for chrom in ds:
            if chrom not in res:
                res[chrom] = IntervalTree()
            res[chrom].update(ds[chrom])
    return res

def overlap_length(range1, range2):
    """
    Compute the length of overlap
    between two ranges given as (start,end)
    """
    return min(range1[1], range2[1]) - max(range1[0], range2[0])

def main():
    # command line arguments
    DESC = "Add Annotation Edit Distance (AED) to EvidenceModeler GFF3"
    USAGE = "python %s -h (display full usage doc)" % sys.argv[0]
    parser = argparse.ArgumentParser(description=DESC, usage=USAGE)
    parser.add_argument('-g', '--gff', help='EVM output GFF3', required=True)
    parser.add_argument('-w', '--weights', help='EVM weights TSV', required=True)
    parser.add_argument('-l', '--feature_types', help='A TSV file denoting which feature types to use as evidence for each source', required=True)
    parser.add_argument('-t', '--transcript', help='Transcript evidence GFF3', default=None)
    parser.add_argument('-p', '--protein', help='Protein evidence GFF3', default=None)
    parser.add_argument('-a', '--prediction', help='Gene predictions GFF3', default=None)
    parser.add_argument('-f', '--force', help='Force override existng gff DBs', action='store_true', default=False)
    parser.add_argument('-i', '--ignore_ab_initio', help='Ab-initio predictions will not be considered when calculating AED', action='store_true', default=False)
    parser.add_argument('-o', '--out_gff', help='Output GFF with AEDs', required=True)
    args = parser.parse_args()

    # if input GFF is empty, print empty GFF output and exit
    if os.path.getsize(args.gff) == 0:
        print("Input GFF is empty.")
        with open(args.out_gff, 'w') as fo:
            print('##gff-version 3', file=fo)
        sys.exit(0)

    # Read weights TSV
    weights_df = pd.read_csv(args.weights, sep='\t', names=['class','source','weight'], index_col='source')
    # Read feature types TSV
    ftypes_df = pd.read_csv(args.feature_types, sep='\t', names=['source','feature_types'], index_col='source', converters={'feature_types': lambda x: set(x.split(','))})
    assert set(ftypes_df.index) == set(weights_df.index), "Must list the same sources in the weights and in the feature types files"
    weights_df = weights_df.join(ftypes_df)
    if args.ignore_ab_initio:
        weights_df = weights_df.query('`class` != "ABINITIO_PREDICTION"')
    weights_df['relative_weight'] = weights_df['weight'] / sum(weights_df['weight'])

    # Read EVM GFF
    evm_gff = create_gff_db(args.gff, args.force)

    # Read ab-initio and evidence GFFs
    evidence = []
    for gff in [args.transcript, args.protein, args.prediction]:
        if not gff:
            continue
        evidence.append(gff_to_intervals(gff, feature_types=weights_df['feature_types'], merge_overlaps=True))
    all_evidence = merge_gff_iv_trees(evidence)

    # Assign AEDs
    print('Calculating AEDs...')
    with open(args.out_gff,'w') as fo:
        print('##gff-version 3', file=fo)
        for evm_feat in evm_gff.all_features():
            if evm_feat.featuretype != 'mRNA':
                print(str(evm_feat), file=fo)
                continue
            mrna = evm_feat
            chrom = mrna.seqid
            mrna_sup = {}
            mrna_cov = 0
            sup_len = 0
            tot_exon_len = 0
            for exon in evm_gff.children(mrna, featuretype='exon'):
                tot_exon_len += exon.end - exon.start
                if chrom in all_evidence:
                    sup_features = all_evidence[chrom][exon.start:exon.end]
                else:
                    sup_features = []
                # for non-weighted AED
                sup_features_tree = IntervalTree(sup_features)
                sup_features_tree.merge_overlaps()
                for supp_feature in sup_features_tree:
                    mrna_cov += overlap_length((exon.start,exon.end), (supp_feature.begin,supp_feature.end))
                    sup_len += supp_feature.end - supp_feature.begin
                # for weighted AED
                for supp_feature in sup_features:
                    overlap_len = overlap_length((exon.start,exon.end), (supp_feature.begin,supp_feature.end))
                    source = supp_feature.data
                    if source not in mrna_sup:
                        mrna_sup[source] = [0,0]
                    mrna_sup[source][0] += overlap_len
                    mrna_sup[source][1] += supp_feature.end - supp_feature.begin

            # calculate non-weighted AED
            if sup_len == 0:
                mrna_aed = 1
            else:
                nw_sensitivity = mrna_cov / tot_exon_len
                nw_specificity = mrna_cov / sup_len
                nw_congruence = (nw_sensitivity + nw_specificity)/2
                mrna_aed = 1 - nw_congruence
            mrna_aed = ("%.2f" % mrna_aed)
            mrna['AED'] = [mrna_aed]
            # calculate weighted AED
            mrna_waed = 0
            for source in weights_df.index:
                if source in mrna_sup:
                    sensitivity = mrna_sup[source][0] / tot_exon_len
                    specificity = mrna_sup[source][0] / mrna_sup[source][1]
                    congruence = (sensitivity + specificity)/2
                    aed = 1 - congruence
                else:
                   aed = 1
                aed_weighted = aed * weights_df.loc[source]['relative_weight']
                mrna_waed += aed_weighted
            mrna_waed = ("%.2f" % mrna_waed)
            mrna['wAED'] = [mrna_waed]
            
            print(str(mrna), file=fo)

# MAIN
if __name__ == "__main__":
    main()
