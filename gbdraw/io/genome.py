#!/usr/bin/env python
# coding: utf-8

import os
import logging
from typing import List, Dict, Set

from Bio import SeqIO
from Bio.SeqFeature import CompoundLocation, SeqFeature
from Bio.SeqRecord import SeqRecord
from BCBio import GFF

from .record_select import parse_record_selector, reverse_records, select_record
from ..exceptions import InputFileError, ParseError, ValidationError

logger = logging.getLogger(__name__)


def _attach_source_annotations(record: SeqRecord, source_file: str) -> None:
    if getattr(record, "annotations", None) is None:
        record.annotations = {}
    record.annotations["gbdraw_source_file"] = source_file
    record.annotations["gbdraw_source_basename"] = os.path.basename(source_file)
    # Coordinate mapping from current left-to-right index (0-based) to absolute
    # genomic coordinate: coord = base + step * index.
    # Defaults for untouched records: 1..len (step=+1).
    if "gbdraw_coord_base" not in record.annotations:
        record.annotations["gbdraw_coord_base"] = 1
    if "gbdraw_coord_step" not in record.annotations:
        record.annotations["gbdraw_coord_step"] = 1


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
            raise InputFileError(f"File does not exist or is not accessible: {gbk_file}")
        try:
            logger.info("INFO: Loading GenBank file {}".format(gbk_file))
            records_list = list(SeqIO.parse(gbk_file, "genbank"))
            if mode in {"linear", "circular"}:
                selector_raw = record_selectors[file_idx] if record_selectors and file_idx < len(record_selectors) else None
                selector = parse_record_selector(selector_raw)
                if selector is not None:
                    records_list = select_record(records_list, selector, log=logger)
            if mode == "linear":
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
                    _attach_source_annotations(record, gbk_file)
                    record_list.append(record)
                    id_list.append(record.id)  # type: ignore
            elif mode == "circular":
                reverse_flag = reverse_flags[file_idx] if reverse_flags and file_idx < len(reverse_flags) else False
                records_list = reverse_records(records_list, reverse_flag, log=logger)
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
                    _attach_source_annotations(record, gbk_file)
                    record_list.append(record)
                    id_list.append(record.id)  # type: ignore
        except ValueError as e:  # Catching common exception when parsing GenBank files
            if str(e).startswith("Record selector"):
                logger.error(f"ERROR: {e}")
                raise ValidationError(str(e)) from e
            logger.error(
                f"ERROR: error parsing GenBank file {gbk_file}. It may be corrupt or in the wrong format. Error: {e}"
            )
            raise ParseError(
                f"Error parsing GenBank file {gbk_file}. It may be corrupt or in the wrong format."
            ) from e
        except Exception as e:  # A more generic catch-all for other unexpected issues
            logger.error(
                f"ERROR: an unexpected error occurred while processing {gbk_file}: {e}"
            )
            raise ParseError(
                f"Unexpected error while processing {gbk_file}: {e}"
            ) from e
    logger.info("INFO:              ... finished loading GenBank file(s)")
    logger.info(f"INFO: Number of sequences loaded to gbdraw: {len(id_list)}")
    if len(id_list) < 1:
        logger.error(
            "ERROR: No valid GenBank records were loaded. Please check your input files."
        )
        raise ValidationError(
            "No valid GenBank records were loaded. Please check your input files."
        )
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
                raise ParseError(
                    f"Failed to merge sequences for record {gff_record.id} with corresponding FASTA entry: {e}"
                ) from e
            merged_records.append(gff_record)
        else:
            # If no matching FASTA record, raise error and stop the process
            logger.error(
                f"ERROR: No matching FASTA record found for GFF record {gff_record.id}. Please ensure that all GFF records have corresponding FASTA entries."
            )
            raise ValidationError(
                f"No matching FASTA record found for GFF record {gff_record.id}. Please ensure that all GFF records have corresponding FASTA entries."
            )
    return merged_records


def _gff3_feature_id(feature: SeqFeature) -> str | None:
    """Return the canonical GFF3 ID for a parsed feature, if present."""
    raw_id = feature.qualifiers.get("ID")
    if isinstance(raw_id, (list, tuple)):
        raw_id = raw_id[0] if raw_id else None
    feature_id = str(raw_id or "").strip()
    return feature_id or None


def _merge_gff3_qualifiers(features: list[SeqFeature]) -> dict:
    """Merge qualifiers while retaining each multipart CDS phase in part order."""
    merged: dict = {}
    for feature in features:
        for key, raw_values in feature.qualifiers.items():
            values = list(raw_values) if isinstance(raw_values, (list, tuple)) else [raw_values]
            target_values = merged.setdefault(key, [])
            for value in values:
                if key == "phase" or value not in target_values:
                    target_values.append(value)
    return merged


def _normalize_gff3_feature_list(features: list[SeqFeature]) -> tuple[list[SeqFeature], int]:
    """Collapse same-ID GFF3 rows into multipart Biopython features."""
    grouped_features: list[list[SeqFeature]] = []
    group_indexes: dict[tuple[str, str], int] = {}

    for feature in features:
        feature_id = _gff3_feature_id(feature)
        location = getattr(feature, "location", None)
        if feature_id is None or location is None:
            grouped_features.append([feature])
            continue

        key = (str(feature.type), feature_id)
        group_index = group_indexes.get(key)
        if group_index is None:
            group_indexes[key] = len(grouped_features)
            grouped_features.append([feature])
        else:
            grouped_features[group_index].append(feature)

    normalized: list[SeqFeature] = []
    merged_count = 0
    for feature_group in grouped_features:
        feature = feature_group[0]
        if len(feature_group) > 1:
            parts = [
                part
                for grouped_feature in feature_group
                for part in grouped_feature.location.parts
            ]
            feature.location = CompoundLocation(parts, operator="join")
            feature.qualifiers = _merge_gff3_qualifiers(feature_group)
            feature.id = _gff3_feature_id(feature) or feature.id
            feature.sub_features = [
                sub_feature
                for grouped_feature in feature_group
                for sub_feature in (getattr(grouped_feature, "sub_features", None) or [])
            ]
            merged_count += 1

        sub_features = getattr(feature, "sub_features", None) or []
        if sub_features:
            feature.sub_features, child_merge_count = _normalize_gff3_feature_list(sub_features)
            merged_count += child_merge_count
        normalized.append(feature)

    return normalized, merged_count


def _normalize_gff3_multipart_features(record: SeqRecord) -> SeqRecord:
    """Normalize repeated GFF3 IDs into ``CompoundLocation`` features in-place.

    GFF3 uses multiple rows with the same ID to represent a single discontinuous
    feature. BCBio-GFF retains those rows as separate ``SeqFeature`` objects, so
    normalize them before feature filtering and rendering.
    """
    record.features, merged_count = _normalize_gff3_feature_list(record.features)
    if merged_count:
        logger.info(
            "INFO: Normalized %d multipart GFF3 feature(s) for record %s.",
            merged_count,
            record.id,
        )
    return record


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
    keep_all_features: bool = False,
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
        raise ValidationError("Number of GFF3 files does not match number of FASTA files.")

    for file_idx, (gff_file, fasta_file) in enumerate(zip(gff_list, fasta_list)):
        if not os.path.isfile(gff_file):
            logger.error(f"ERROR: File does not exist or is not accessible: {gff_file}")
            raise InputFileError(f"File does not exist or is not accessible: {gff_file}")
        if not os.path.isfile(fasta_file):
            logger.error(f"ERROR: File does not exist or is not accessible: {fasta_file}")
            raise InputFileError(f"File does not exist or is not accessible: {fasta_file}")

        try:
            logger.info("INFO: Loading GFF3 file {}".format(gff_file))
            parsed_gff_records = [
                _normalize_gff3_multipart_features(record) for record in GFF.parse(gff_file)
            ]
            if keep_all_features or selected_features_set is None:
                gff_records = parsed_gff_records
            else:
                feature_types_to_keep = set(selected_features_set)
                gff_records = [
                    filter_features_by_type(record, feature_types_to_keep)
                    for record in parsed_gff_records
                ]
            logger.info("INFO: Loading FASTA file {}".format(fasta_file))
            fasta_records: list[SeqRecord] = list(SeqIO.parse(fasta_file, "fasta"))
            merged_records = merge_gff_fasta_records(gff_records, fasta_records)

            if mode in {"linear", "circular"}:
                selector_raw = record_selectors[file_idx] if record_selectors and file_idx < len(record_selectors) else None
                selector = parse_record_selector(selector_raw)
                if selector is not None:
                    merged_records = select_record(merged_records, selector, log=logger)
            if mode == "linear":
                if load_comparison and len(merged_records) > 1:
                    merged_records = [merged_records[0]]
                reverse_flag = reverse_flags[file_idx] if reverse_flags and file_idx < len(reverse_flags) else False
                merged_records = reverse_records(merged_records, reverse_flag, log=logger)
            elif mode == "circular":
                reverse_flag = reverse_flags[file_idx] if reverse_flags and file_idx < len(reverse_flags) else False
                merged_records = reverse_records(merged_records, reverse_flag, log=logger)
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
                _attach_source_annotations(record, gff_file)
                record_list.append(record)
                id_list.append(record.id)  # type: ignore
        except ValueError as e:
            if str(e).startswith("Record selector"):
                logger.error(f"ERROR: {e}")
                raise ValidationError(str(e)) from e
            logger.error(
                f"ERROR: error parsing GFF3/FASTA files ({gff_file}, {fasta_file}). Error: {e}"
            )
            raise ParseError(
                f"Error parsing GFF3/FASTA files ({gff_file}, {fasta_file})."
            ) from e
        except Exception as e:
            logger.error(
                f"ERROR: an unexpected error occurred while processing {gff_file} or {fasta_file}: {e}"
            )
            raise ParseError(
                f"Unexpected error while processing {gff_file} or {fasta_file}: {e}"
            ) from e

    logger.info("INFO:              ... finished loading GFF3 and FASTA files")
    logger.info(f"INFO: Number of sequences loaded to gbdraw: {len(record_list)}")
    if len(record_list) < 1:
        logger.error(
            "ERROR: No valid records were loaded after merging GFF3 and FASTA files. Please check your input files."
        )
        raise ValidationError(
            "No valid records were loaded after merging GFF3 and FASTA files. Please check your input files."
        )
    return record_list


__all__ = [
    "load_gbks",
    "load_gff_fasta",
    "merge_gff_fasta_records",
    "scan_features_recursive",
    "filter_features_by_type",
]


