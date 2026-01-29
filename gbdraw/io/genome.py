#!/usr/bin/env python
# coding: utf-8

import os
import sys
import logging
from typing import Generator, List, Dict, Set

from Bio import SeqIO
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord
from BCBio import GFF

from .record_select import parse_record_selector, reverse_records, select_record

logger = logging.getLogger(__name__)


def load_gbks(
    gbk_list: List[str],
    mode: str,
    load_comparison: bool = False,
    record_selectors: list[str] | None = None,
    reverse_flags: list[bool] | None = None,
) -> list[SeqRecord]:
    record_list: list[SeqRecord] = []
    id_list: list[str] = []
    logger.info("INFO: Loading GenBank file(s)...")
    if load_comparison:
        logger.info(
            "INFO: BLAST results were provided. Only the first entry from each GenBank file will be loaded."
        )
    for file_idx, gbk_file in enumerate(gbk_list):
        if not os.path.isfile(gbk_file):
            logger.error(f"ERROR: File does not exist or is not accessible: {gbk_file}")
            sys.exit(1)
        try:
            logger.info("INFO: Loading GenBank file {}".format(gbk_file))
            records_list = list(SeqIO.parse(gbk_file, "genbank"))
            if mode == "linear":
                selector_raw = record_selectors[file_idx] if record_selectors and file_idx < len(record_selectors) else None
                selector = parse_record_selector(selector_raw)
                if selector is not None:
                    records_list = select_record(records_list, selector, log=logger)
                if len(gbk_list) > 1 and load_comparison and len(records_list) > 1:
                    records_list = [records_list[0]]
                reverse_flag = reverse_flags[file_idx] if reverse_flags and file_idx < len(reverse_flags) else False
                records_list = reverse_records(records_list, reverse_flag, log=logger)
                if len(gbk_list) == 1 or not load_comparison:
                    logger.info("INFO: Importing all entries...")
                for record in records_list:
                    if record.id in id_list:
                        logger.warning(
                            f"WARNING: Record {record.id} seems to have been already loaded. Check for duplicates!"
                        )
                    record_list.append(record)
                    id_list.append(record.id)  # type: ignore
            elif mode == "circular":
                logger.info("INFO: Importing all entries...")
                for record in records_list:
                    if record.id in id_list:
                        logger.warning(
                            f"WARNING: Record {record.id} seems to have been already loaded. Check for duplicates!"
                        )
                    if "topology" in record.annotations:
                        topology: str = record.annotations["topology"]  # type: ignore
                        if topology == "linear":
                            logger.warning(
                                f"WARNING: The annotation indicates that record {record.id} is linear. Are you sure you want to visualize it as circular?"
                            )
                        elif topology == "circular":
                            pass
                        else:
                            logger.warning(
                                f"WARNING: Topology information not available for {record.id}."
                            )
                    record_list.append(record)
                    id_list.append(record.id)  # type: ignore
        except ValueError as e:  # Catching common exception when parsing GenBank files
            logger.error(
                f"ERROR: error parsing GenBank file {gbk_file}. It may be corrupt or in the wrong format. Error: {e}"
            )
            sys.exit(1)
        except Exception as e:  # A more generic catch-all for other unexpected issues
            logger.error(
                f"ERROR: an unexpected error occurred while processing {gbk_file}: {e}"
            )
            sys.exit(1)
    logger.info("INFO:              ... finished loading GenBank file(s)")
    logger.info(f"INFO: Number of sequences loaded to gbdraw: {len(id_list)}")
    if len(id_list) < 1:
        logger.error(
            "ERROR: No valid GenBank records were loaded. Please check your input files."
        )
        sys.exit(1)
    return record_list


def merge_gff_fasta_records(
    gff_records: List[SeqRecord], fasta_records: List[SeqRecord]
) -> List[SeqRecord]:
    merged_records: List[SeqRecord] = []
    fasta_dict: Dict[str, SeqRecord] = {record.id: record for record in fasta_records}
    for gff_record in gff_records:
        if gff_record.id in fasta_dict:
            try:
                fasta_record = fasta_dict[gff_record.id]
                # If FASTA record exists, record.seq is FASTA sequence
                gff_record.seq = fasta_record.seq
            except Exception as e:
                logger.error(
                    f"ERROR: Failed to merge sequences for record {gff_record.id} with corresponding FASTA entry: {e}"
                )
                sys.exit(1)
            merged_records.append(gff_record)
        else:
            # If no matching FASTA record, raise error and stop the process
            logger.error(
                f"ERROR: No matching FASTA record found for GFF record {gff_record.id}. Please ensure that all GFF records have corresponding FASTA entries."
            )
            sys.exit(1)
    return merged_records


def scan_features_recursive(
    features: List[SeqFeature], feature_types_to_keep: Set[str]
) -> List[SeqFeature]:
    filtered_list: list[SeqFeature] = []
    for feature in features:
        if feature.type in feature_types_to_keep:
            new_feature = SeqFeature(
                location=feature.location,
                type=feature.type,
                qualifiers=feature.qualifiers,
                id=feature.id,
            )
            filtered_list.append(new_feature)

        if feature.sub_features:
            # Recursively call and extend the list with the results
            filtered_list.extend(
                scan_features_recursive(feature.sub_features, feature_types_to_keep)
            )

    return filtered_list


def filter_features_by_type(record: SeqRecord, feature_types_to_keep: Set[str]) -> SeqRecord:
    new_record = SeqRecord(
        seq=record.seq,
        id=record.id,
        name=record.name,
        description=record.description,
        dbxrefs=record.dbxrefs,
        annotations=record.annotations,
    )

    filtered_features = scan_features_recursive(record.features, feature_types_to_keep)
    new_record.features = filtered_features

    logger.info(
        f"INFO: For record {record.id}, filtered to {len(filtered_features)} features of types: {', '.join(feature_types_to_keep)}."
    )
    return new_record


def load_gff_fasta(
    gff_list: List[str],
    fasta_list: List[str],
    mode: str,
    selected_features_set=None,
    load_comparison: bool = False,
    record_selectors: list[str] | None = None,
    reverse_flags: list[bool] | None = None,
) -> list[SeqRecord]:
    record_list: list[SeqRecord] = []
    id_list: list[str] = []
    logger.info("INFO: Loading GFF3/FASTA file(s)...")
    if load_comparison:
        logger.info(
            "INFO: BLAST results were provided. Only the first entry from each GFF3/FASTA pair will be loaded."
        )

    if len(gff_list) != len(fasta_list):
        logger.error("ERROR: Number of GFF3 files does not match number of FASTA files.")
        sys.exit(1)

    for file_idx, (gff_file, fasta_file) in enumerate(zip(gff_list, fasta_list)):
        if not os.path.isfile(gff_file):
            logger.error(f"ERROR: File does not exist or is not accessible: {gff_file}")
            sys.exit(1)
        if not os.path.isfile(fasta_file):
            logger.error(f"ERROR: File does not exist or is not accessible: {fasta_file}")
            sys.exit(1)

        try:
            logger.info("INFO: Loading GFF3 file {}".format(gff_file))
            gff_records: list[SeqRecord] = [
                filter_features_by_type(record, selected_features_set) for record in GFF.parse(gff_file)
            ]
            logger.info("INFO: Loading FASTA file {}".format(fasta_file))
            fasta_records: list[SeqRecord] = list(SeqIO.parse(fasta_file, "fasta"))
            merged_records = merge_gff_fasta_records(gff_records, fasta_records)

            if mode == "linear":
                selector_raw = record_selectors[file_idx] if record_selectors and file_idx < len(record_selectors) else None
                selector = parse_record_selector(selector_raw)
                if selector is not None:
                    merged_records = select_record(merged_records, selector, log=logger)
                if load_comparison and len(merged_records) > 1:
                    merged_records = [merged_records[0]]
                reverse_flag = reverse_flags[file_idx] if reverse_flags and file_idx < len(reverse_flags) else False
                merged_records = reverse_records(merged_records, reverse_flag, log=logger)
            elif mode == "circular":
                for record in merged_records:
                    if "topology" in record.annotations:
                        topology: str = record.annotations["topology"]  # type: ignore
                        if topology == "linear":
                            logger.warning(
                                f"WARNING: The annotation indicates that record {record.id} is linear. Are you sure you want to visualize it as circular?"
                            )
                        elif topology == "circular":
                            pass
                        else:
                            logger.warning(
                                f"WARNING: Topology information not available for {record.id}."
                            )

            for record in merged_records:
                if record.id in id_list:
                    logger.warning(
                        f"WARNING: Record {record.id} seems to have been already loaded. Check for duplicates!"
                    )
                record_list.append(record)
                id_list.append(record.id)  # type: ignore
        except ValueError as e:
            logger.error(
                f"ERROR: error parsing GFF3/FASTA files ({gff_file}, {fasta_file}). Error: {e}"
            )
            sys.exit(1)
        except Exception as e:
            logger.error(
                f"ERROR: an unexpected error occurred while processing {gff_file} or {fasta_file}: {e}"
            )
            sys.exit(1)

    logger.info("INFO:              ... finished loading GFF3 and FASTA files")
    logger.info(f"INFO: Number of sequences loaded to gbdraw: {len(record_list)}")
    if len(record_list) < 1:
        logger.error(
            "ERROR: No valid records were loaded after merging GFF3 and FASTA files. Please check your input files."
        )
        sys.exit(1)
    return record_list


__all__ = [
    "load_gbks",
    "load_gff_fasta",
    "merge_gff_fasta_records",
    "scan_features_recursive",
    "filter_features_by_type",
]


