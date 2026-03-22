#!/usr/bin/env python
# coding: utf-8

import sys
import os
import argparse
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from contextlib import redirect_stdout


def _get_args():
    parser = argparse.ArgumentParser(
        description='Extract the longest isoforms of protein-coding genes from a NCBI RefSeq euaryotic genome assembly')
    parser.add_argument(
        '-i', '--input', help='GenBank flat file format of the genomic sequence(s) (required)', type=str, required=True)
    parser.add_argument(
        '-o', '--output', help='output fasta-formatted file (optional)', type=str, default="out.faa")
    parser.add_argument(
        '-t', '--table', help='output tab-separated table (optional)', type=str, default="out.tbl")
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()
    return args

def _gb_record_to_proteins(gb_record):
    protein_records = defaultdict(list)
    for feature in gb_record.features:
        if feature.type == 'CDS':
            if 'protein_id' in feature.qualifiers.keys():
                protein_id = feature.qualifiers['protein_id'][0]
            else:
                continue
            if 'locus_tag' in feature.qualifiers.keys():
                gene_id = feature.qualifiers['locus_tag'][0]
            elif 'gene' in feature.qualifiers.keys():
                gene_id = feature.qualifiers['gene'][0]
            else:
                continue
            product_name = feature.qualifiers['product'][0]
            aa_seq = Seq(feature.qualifiers['translation'][0])
            aa_record = SeqRecord(aa_seq, id=protein_id,
                                  name=gene_id, description="{} [{}]".format(product_name, gb_record.annotations['organism']))
            protein_records[gene_id].append(aa_record)
        else:
            continue
    out_records = [max(protein_records[gene_id], key=lambda k: len(
        protein_records[gene_id])) for gene_id in protein_records.keys()]
    return out_records

def _gb_to_proteins(in_gb):
    gb_records = SeqIO.parse(in_gb, 'genbank')
    protein_records = [_gb_record_to_proteins(
        gb_record) for gb_record in gb_records]
    protein_records = [
        protein for proteins in protein_records for protein in proteins]
    return protein_records

def main():
    args = _get_args()
    in_gb = args.input
    out_file = args.output
    out_table = args.table
    protein_records = _gb_to_proteins(in_gb)
    with open(out_table, "w") as output_table:
        with redirect_stdout(output_table):
            for record in protein_records:
                print('{}\t{}'.format(record.name, record.id))
    with open(out_file, "w") as output_handle:
        SeqIO.write(protein_records, output_handle, "fasta")

if __name__ == "__main__":
    main()

