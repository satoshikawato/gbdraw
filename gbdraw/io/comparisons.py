#!/usr/bin/env python
# coding: utf-8

import os
import logging
from typing import List

import pandas as pd
from pandas import DataFrame

from ..configurators import BlastMatchConfigurator

logger = logging.getLogger(__name__)


def load_comparisons(
    comparison_files: List[str], blast_config: BlastMatchConfigurator
) -> List[DataFrame]:
    evalue_threshold: float = blast_config.evalue
    bitscore_threshold: float = blast_config.bitscore
    identity_threshold: float = blast_config.identity
    logger.info(
        "INFO: BLAST output visualization settings: e-value threshold: {}; bitscore threshold: {}; identity threshold: {}".format(
            evalue_threshold, bitscore_threshold, identity_threshold
        )
    )
    comparison_list: list[DataFrame] = []
    logger.info("INFO: Loading comparison file(s)...")
    for comparison_file in comparison_files:
        logger.info("INFO: Loading {}".format(comparison_file))
        if not os.path.isfile(comparison_file):
            logger.warning(
                f"WARNING: File does not exist or is not accessible: {comparison_file}"
            )
            continue
        try:
            df: DataFrame = pd.read_csv(
                comparison_file,
                sep="\t",
                comment="#",
                names=(
                    "query",
                    "subject",
                    "identity",
                    "alignment_length",
                    "mismatches",
                    "gap_opens",
                    "qstart",
                    "qend",
                    "sstart",
                    "send",
                    "evalue",
                    "bitscore",
                ),
            )
            df = df[
                (df["evalue"] <= evalue_threshold)
                & (df["bitscore"] >= bitscore_threshold)
                & (df["identity"] >= identity_threshold)
            ]
            comparison_list.append(df)
        except ValueError as e:  # Catching common exception when parsing comparison files
            logger.warning(
                f"WARNING: Error parsing comparison file {comparison_file}. It may be corrupt or in the wrong format. Error: {e}"
            )
        except Exception as e:  # A more generic catch-all for other unexpected issues
            logger.error(
                f"ERROR: An unexpected error occurred while processing {comparison_file}: {e}"
            )
    logger.info("INFO:             ... finished loading comparison file(s)")
    return comparison_list


__all__ = ["load_comparisons"]


