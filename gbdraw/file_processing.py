#!/usr/bin/env python
# coding: utf-8
import os
import sys
import logging
import tomllib
import cairosvg


from pathlib import Path
from typing import Optional, Generator, List, Dict, Set
from importlib import resources
from importlib.abc import Traversable
import pandas as pd
from pandas import DataFrame
from svgwrite import Drawing
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord
from BCBio import GFF


from .object_configurators import BlastMatchConfigurator
logger = logging.getLogger()
handler = logging.StreamHandler(sys.stdout)


def load_gbks(gbk_list: List[str], mode: str, load_comparison=False) -> list[SeqRecord]:
    """
    Loads GenBank records from the given list of file paths based on the specified mode.

    Args:
        gbk_list (List[str]): List of paths to GenBank files.
        mode (str): Mode of operation, either 'linear' or 'circular'.
        load_comparison (bool, optional): Flag to load only the first entry for comparison. Defaults to False.

    Returns:
        list[SeqRecord]: List of SeqRecord objects parsed from the GenBank files.

    The function processes each GenBank file in the provided list. In 'linear' mode, it loads all entries
    from a single file or the first entry from multiple files. In 'circular' mode, it loads all entries,
    issuing warnings for sequences with linear topology.
    """
    """
    Loads GenBank records from provided file paths based on the specified mode (linear or circular).

    This function processes a list of GenBank file paths, loading and appending the corresponding 
    SeqRecord objects to a list. It adapts its behavior based on the specified mode:
    - In 'linear' mode, it imports all entries from a single file or the first entry from multiple files.
    - In 'circular' mode, it imports all entries and logs warnings if the topology of any sequence is 
      linear or unspecified, as circular diagrams may not be suitable.

    Args:
        gbk_list (list): A list of file paths to GenBank files.
        mode (str): The mode of operation, either 'linear' or 'circular', determining how the files
                    are processed and how many entries are loaded.

    Returns:
        list: A list of BioPython SeqRecord objects loaded from the specified GenBank files.

    The function iterates through each file in the provided list, logging the file being processed.
    It utilizes BioPython's SeqIO.parse to load records and conditionally appends them to the output list
    based on the specified mode and the content of each file. Warnings are logged to inform the user of 
    potential issues with sequence topology in relation to the desired diagram type.
    """
    record_list: list[SeqRecord] = []
    id_list: list[str] = []
    logger.info("INFO: Loading GenBank file(s)...")
    if load_comparison:
        logger.info(
            "INFO: BLAST results were provided. Only the first entry from each GenBank file will be loaded.")
    for gbk_file in gbk_list:
        if not os.path.isfile(gbk_file):
            logger.error(
                f"ERROR: File does not exist or is not accessible: {gbk_file}")
            sys.exit(1)
        try:
            logger.info("INFO: Loading GenBank file {}".format(gbk_file))
            records: Generator[SeqRecord, None,
                               None] = SeqIO.parse(gbk_file, 'genbank')
            if mode == "linear":
                if len(gbk_list) == 1:
                    logger.info("INFO: Importing all entries...")
                    for record in records:
                        if record.id in id_list:
                            logger.warning(f"WARNING: Record {record.id} seems to have been already loaded. Check for duplicates!")
                        record_list.append(record)
                        id_list.append(record.id)  # type: ignore
                elif len(gbk_list) > 1 and load_comparison:
                    record: SeqRecord = next(records)
                    if record.id in id_list:
                        logger.warning(f"WARNING: Record {record.id} seems to have been already loaded. Check for duplicates!")
                    record_list.append(record)
                    id_list.append(record.id)  # type: ignore
                else:
                    logger.info("INFO: Importing all entries...")
                    for record in records:
                        if record.id in id_list:
                            logger.warning(f"WARNING: Record {record.id} seems to have been already loaded. Check for duplicates!")
                        record_list.append(record)
                        id_list.append(record.id)           # type: ignore
            elif mode == "circular":
                logger.info("INFO: Importing all entries...")
                for record in records:
                    if record.id in id_list:
                        logger.warning(f"WARNING: Record {record.id} seems to have been already loaded. Check for duplicates!")
                    if "topology" in record.annotations:
                        # type: ignore
                        topology: str = record.annotations["topology"]
                        if topology == "linear":
                            logger.warning(f"WARNING: The annotation indicates that record {record.id} is linear. Are you sure you want to visualize it as circular?")
                        elif topology == "circular":
                            pass
                        else:
                            logger.warning(f"WARNING: Topology information not available for {record.id}.")
                    record_list.append(record)
                    id_list.append(record.id)  # type: ignore
        except ValueError as e:  # Catching common exception when parsing GenBank files
            logger.error(f"ERROR: error parsing GenBank file {gbk_file}. It may be corrupt or in the wrong format. Error: {e}")
            sys.exit(1)
        except Exception as e:  # A more generic catch-all for other unexpected issues
            logger.error(f"ERROR: an unexpected error occurred while processing {gbk_file}: {e}")
            sys.exit(1)
    logger.info("INFO:              ... finished loading GenBank file(s)")
    logger.info(f"INFO: Number of sequences loaded to gbdraw: {len(id_list)}")
    if len(id_list) < 1:
        logger.error("ERROR: No valid GenBank records were loaded. Please check your input files.")
        sys.exit(1)
    return record_list

def merge_gff_fasta_records(gff_records: List[SeqRecord], fasta_records: List[SeqRecord]) -> List[SeqRecord]:
    """
    Merges GFF and FASTA records based on matching IDs.

    Args:
        gff_records (List[SeqRecord]): List of SeqRecord objects from GFF files.
        fasta_records (List[SeqRecord]): List of SeqRecord objects from FASTA files.

    Returns:
        List[SeqRecord]: Merged list of SeqRecord objects with annotations from GFF files.

    This function iterates through the GFF records and searches for matching FASTA records by ID.
    If a match is found, it merges the features from the GFF record into the corresponding FASTA record.
    If no match is found, the GFF record is added to the merged list as is. Warnings are logged for
    any GFF records without matching FASTA entries.
    """
    merged_records: List[SeqRecord] = []
    fasta_dict: Dict[str, SeqRecord] = {record.id: record for record in fasta_records}
    for gff_record in gff_records:
        if gff_record.id in fasta_dict:
            try:
                fasta_record = fasta_dict[gff_record.id]
                # If FASTA record exists, record.seq is FASTA sequence
                gff_record.seq = fasta_record.seq
            except Exception as e:
                logger.error(f"ERROR: Failed to merge sequences for record {gff_record.id} with corresponding FASTA entry: {e}")
                sys.exit(1)
            merged_records.append(gff_record)
        else:
            # If no matching FASTA record, raise error and stop the process
            logger.error(f"ERROR: No matching FASTA record found for GFF record {gff_record.id}. Please ensure that all GFF records have corresponding FASTA entries.")
            sys.exit(1)
    return merged_records

def scan_features_recursive(features: List[SeqFeature], feature_types_to_keep: Set[str]) -> List[SeqFeature]:
    """
    Recursively scans a list of features and their sub-features,
    extracting only those of the specified types.
    """
    filtered_list = []
    for feature in features:
        if feature.type in feature_types_to_keep:
            new_feature = SeqFeature(
                location=feature.location,
                type=feature.type,
                qualifiers=feature.qualifiers,
                id=feature.id
            )
            filtered_list.append(new_feature)

        if feature.sub_features:
            # Recursively call and extend the list with the results
            filtered_list.extend(scan_features_recursive(feature.sub_features, feature_types_to_keep))
            
    return filtered_list

def filter_features_by_type(record: SeqRecord, feature_types_to_keep: Set[str]) -> SeqRecord:
    """
    Creates a new SeqRecord containing only features of specified types.
    Sub-features of the desired types are extracted and promoted to top-level features.
    """
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
    
    logger.info(f"INFO: For record {record.id}, filtered to {len(filtered_features)} features of types: {', '.join(feature_types_to_keep)}.")
    return new_record


def load_gff_fasta(gff_list: List[str], fasta_list: List[str], mode: str, selected_features_set=None, load_comparison=False) -> list[SeqRecord]:
    """
    """
    gff_record_list_tmp : list[SeqRecord] = []
    gff_id_list_tmp: list[str] = []
    fasta_record_list_tmp : list[SeqRecord] = []
    fasta_id_list_tmp: list[str] = []
    record_list: list[SeqRecord] = []
    logger.info("INFO: Loading GFF3 file(s)...")
    if load_comparison:
        logger.info(
            "INFO: BLAST results were provided. Only the first entry from each GFF3 file will be loaded.")
    for gff_file in gff_list:
        if not os.path.isfile(gff_file):
            logger.error(
                f"ERROR: File does not exist or is not accessible: {gff_file}")
            sys.exit(1)
        try:
            logger.info("INFO: Loading GFF3 file {}".format(gff_file))
            records: Generator[SeqRecord, None,
                               None] = GFF.parse(gff_file)
            if mode == "linear":
                if len(gff_list) == 1:
                    logger.info("INFO: Importing all entries...")
                    for record in records:
                        record = filter_features_by_type(record, selected_features_set)
                        if record.id in gff_id_list_tmp:
                            logger.warning(f"WARNING: Record {record.id} seems to have been already loaded. Check for duplicates!")
                        gff_record_list_tmp.append(record)
                        gff_id_list_tmp.append(record.id)  # type: ignore
                elif len(gff_list) > 1 and load_comparison:
                    record: SeqRecord = next(records)
                    record = filter_features_by_type(record, selected_features_set)
                    if record.id in gff_id_list_tmp:
                        logger.warning(f"WARNING: Record {record.id} seems to have been already loaded. Check for duplicates!")
                    gff_record_list_tmp.append(record)
                    gff_id_list_tmp.append(record.id)  # type: ignore
                else:
                    logger.info("INFO: Importing all entries...")
                    for record in records:
                        record = filter_features_by_type(record, selected_features_set)
                        if record.id in gff_id_list_tmp:
                            logger.warning(f"WARNING: Record {record.id} seems to have been already loaded. Check for duplicates!")
                        gff_record_list_tmp.append(record)
                        gff_id_list_tmp.append(record.id)           # type: ignore
            elif mode == "circular":
                logger.info("INFO: Importing all entries...")
                for record in records:
                    record = filter_features_by_type(record, selected_features_set)
                    if record.id in gff_id_list_tmp:
                        logger.warning(f"WARNING: Record {record.id} seems to have been already loaded. Check for duplicates!")
                    if "topology" in record.annotations:
                        # type: ignore
                        topology: str = record.annotations["topology"]
                        if topology == "linear":
                            logger.warning(f"WARNING: The annotation indicates that record {record.id} is linear. Are you sure you want to visualize it as circular?")
                        elif topology == "circular":
                            pass
                        else:
                            logger.warning(f"WARNING: Topology information not available for {record.id}.")
                    gff_record_list_tmp.append(record)
                    gff_id_list_tmp.append(record.id)  # type: ignore
        except ValueError as e:  # Catching common exception when parsing GenBank files
            logger.error(f"ERROR: error parsing GFF3 file {gff_file}. It may be corrupt or in the wrong format. Error: {e}")
            sys.exit(1)
        except Exception as e:  # A more generic catch-all for other unexpected issues
            logger.error(f"ERROR: an unexpected error occurred while processing {gff_file}: {e}")
            sys.exit(1)
    logger.info("INFO:              ... finished loading GFF3 file(s)")
    logger.info(f"INFO: Number of sequences loaded from GFF3 files: {len(gff_id_list_tmp)}")
    logger.info("INFO: Loading FASTA file(s)...")
    for fasta_file in fasta_list:
        if not os.path.isfile(fasta_file):
            logger.error(
                f"ERROR: File does not exist or is not accessible: {fasta_file}")
            sys.exit(1)
            continue
        try:
            logger.info("INFO: Loading FASTA file {}".format(fasta_file))
            records: Generator[SeqRecord, None,
                               None] = SeqIO.parse(fasta_file, 'fasta')
            if mode == "linear":
                if len(fasta_list) == 1:
                    logger.info("INFO: Importing all entries...")
                    for record in records:
                        if record.id in fasta_id_list_tmp:
                            logger.warning(f"WARNING: Record {record.id} seems to have been already loaded. Check for duplicates!")
                        fasta_record_list_tmp.append(record)
                        fasta_id_list_tmp.append(record.id)  # type: ignore
                elif len(fasta_list) > 1 and load_comparison:
                    record: SeqRecord = next(records)
                    if record.id in fasta_id_list_tmp:
                        logger.warning(f"WARNING: Record {record.id} seems to have been already loaded. Check for duplicates!")
                    fasta_record_list_tmp.append(record)
                    fasta_id_list_tmp.append(record.id)  # type: ignore
                else:
                    logger.info("INFO: Importing all entries...")
                    for record in records:
                        if record.id in fasta_id_list_tmp:
                            logger.warning(f"WARNING: Record {record.id} seems to have been already loaded. Check for duplicates!")
                        fasta_record_list_tmp.append(record)
                        fasta_id_list_tmp.append(record.id)           # type: ignore
            elif mode == "circular":
                logger.info("INFO: Importing all entries...")
                for record in records:
                    if record.id in fasta_id_list_tmp:
                        logger.warning(f"WARNING: Record {record.id} seems to have been already loaded. Check for duplicates!")
                    if "topology" in record.annotations:
                        # type: ignore
                        topology: str = record.annotations["topology"]
                        if topology == "linear":
                            logger.warning(f"WARNING: The annotation indicates that record {record.id} is linear. Are you sure you want to visualize it as circular?")
                        elif topology == "circular":
                            pass
                        else:
                            logger.warning(f"WARNING: Topology information not available for {record.id}.")
                    fasta_record_list_tmp.append(record)
                    fasta_id_list_tmp.append(record.id)  # type: ignore
        except ValueError as e:  # Catching common exception when parsing GenBank files
            logger.error(f"ERROR: error parsing FASTA file {fasta_file}. It may be corrupt or in the wrong format. Error: {e}")
            sys.exit(1)
        except Exception as e:  # A more generic catch-all for other unexpected issues
            logger.error(f"ERROR: an unexpected error occurred while processing {fasta_file}: {e}")
            sys.exit(1)
    logger.info("INFO:              ... finished loading FASTA file(s)")
    logger.info(f"INFO: Number of sequences loaded from FASTA files: {len(fasta_id_list_tmp)}")
    record_list = merge_gff_fasta_records(gff_record_list_tmp, fasta_record_list_tmp)
    logger.info(f"INFO: Number of sequences after merging GFF3 and FASTA files: {len(record_list)}")
    logger.info("INFO:              ... finished loading GFF3 and FASTA files")
    logger.info(f"INFO: Number of sequences loaded to gbdraw: {len(record_list)}")
    if len(record_list) < 1:
        logger.error("ERROR: No valid records were loaded after merging GFF3 and FASTA files. Please check your input files.")
        sys.exit(1)
    return record_list
            

def load_comparisons(
        comparison_files: List[str],
        blast_config: BlastMatchConfigurator) -> List[DataFrame]:
    """
    Loads and filters comparison data from tab-separated BLAST output files.

    Args:
        comparison_files (list[str]): A list of file paths to tab-separated BLAST output files.
        blast_config (BlastMatchConfigurator): Configuration object containing evalue, bitscore, and identity thresholds.

    Returns:
        list[DataFrame]: List of DataFrames containing the filtered comparison data.

    This function reads comparison data from each specified file and applies filters based on the thresholds
    in the blast_config. It retains records that meet the specified criteria for evalue, bitscore, and identity.
    """
    evalue_threshold: float = blast_config.evalue
    bitscore_threshold: float = blast_config.bitscore
    identity_threshold: float = blast_config.identity
    logger.info(
        "INFO: BLAST output visualization settings: e-value threshold: {}; bitscore threshold: {}; identity threshold: {}".format(evalue_threshold,bitscore_threshold,identity_threshold))
    comparison_list: list[DataFrame] = []
    logger.info("INFO: Loading comparison file(s)...")
    for comparison_file in comparison_files:
        logger.info("INFO: Loading {}".format(comparison_file))
        if not os.path.isfile(comparison_file):
            logger.warning(f"WARNING: File does not exist or is not accessible: {comparison_file}")
            continue
        try:
            df: DataFrame = pd.read_csv(
                comparison_file,
                sep='\t',
                comment='#',
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
                    "bitscore"))
            df = df[(df['evalue'] <= evalue_threshold) & (df['bitscore'] >=bitscore_threshold) & (df['identity'] >= identity_threshold)]
            comparison_list.append(df)
        except ValueError as e:  # Catching common exception when parsing GenBank files
            logger.warning(f"WARNING: Error parsing comparison file {comparison_file}. It may be corrupt or in the wrong format. Error: {e}")
        except Exception as e:  # A more generic catch-all for other unexpected issues
            logger.error(f"ERROR: An unexpected error occurred while processing {comparison_file}: {e}")
    logger.info("INFO:             ... finished loading comparison file(s)")
    return comparison_list


def load_default_colors(user_defined_default_colors: str, palette: str = "default") -> DataFrame:
    """
    Load the built-in colour palette from *color_palettes.toml* (section
    `[default]` by default) and merge any user-supplied TSV overrides.

    ```
    Parameters
    ----------
    user_defined_default_colors : str
        Path to a TSV file whose rows are ``feature_type<TAB>HEX``.
        If empty, the built-in palette is used as-is.
    palette : str, optional
        Name of the palette inside *color_palettes.toml*.
        Falls back to ``[default]`` if the name is not found.

    Returns
    -------
    pandas.DataFrame
        Two columns: ``feature_type`` and ``color``.
    """
    column_names = ["feature_type", "color"]

    # ── 1) Load TOML
    try:
        toml_path = resources.files("gbdraw.data").joinpath("color_palettes.toml")
        with toml_path.open("rb") as fh:
            palettes_dict = tomllib.load(fh)
    except Exception as exc:
        logger.error(f"ERROR: failed to read colour_palettes.toml – {exc}")
        raise

    if palette not in palettes_dict:
        logger.warning(f"Palette '{palette}' not found; using [default]")
        palette_dict = palettes_dict.get("default", {})
    else:
        palette_dict = palettes_dict[palette]

    default_colors = (
        pd.DataFrame(palette_dict.items(), columns=column_names)
        .set_index("feature_type")
    )

    # ── 2) Apply user TSV overrides
    if user_defined_default_colors:
        try:
            user_df = (
                pd.read_csv(user_defined_default_colors,
                            sep="\t",
                            names=column_names,
                            header=None,
                            dtype=str)
                .set_index("feature_type")
            )
            # Drop rows with missing colour cells
            missing = user_df["color"].isna()
            if missing.any():
                for ft in user_df[missing].index.tolist():
                    logger.warning(
                        f"WARNING: colour missing for feature '{ft}' "
                        f"in '{user_defined_default_colors}' – "
                        "keeping built-in value."
                    )
                user_df = user_df[~missing]

            default_colors.update(user_df)
            logger.info(f"User overrides applied: {user_defined_default_colors}")

        except FileNotFoundError:
            logger.error(
                f"ERROR: override file '{user_defined_default_colors}' not found")
        except Exception as exc:
            logger.error(
                f"ERROR: failed to read '{user_defined_default_colors}' – {exc}")

    # ── 3) Return tidy DataFrame (index reset for downstream code)
    return default_colors.reset_index()


def read_color_table(color_table_file: str) -> Optional[DataFrame]:
    """
    Reads a color table from a TSV file and errors out immediately if any row has
    the wrong number of fields or missing values.

    Expected columns: feature_type, qualifier_key, value, color, caption

    Returns:
        DataFrame with the user-defined color mappings, or None if no file was provided.
    """
    required_cols = ['feature_type', 'qualifier_key', 'value', 'color', 'caption']

    # If user did not supply -t, just skip and return None
    if not color_table_file:
        return None

    try:
        df = pd.read_csv(
            color_table_file,
            sep='\t',
            header=None,
            names=required_cols,
            dtype=str,
            on_bad_lines='error',  # raise on any row with wrong number of fields
            engine='python'        # required for on_bad_lines
        )
    except pd.errors.ParserError as e:
        logger.error(f"ERROR: Malformed line in '{color_table_file}': {e}")
        sys.exit(1)
    except FileNotFoundError as e:
        logger.error(f"ERROR: Color table file not found: {e}")
        sys.exit(1)
    except Exception as e:
        logger.error(f"ERROR: Failed to read '{color_table_file}': {e}")
        sys.exit(1)

    # Check for any rows with missing values and error out if found
    null_rows = df[df.isnull().any(axis=1)]
    if not null_rows.empty:
        for idx, row in null_rows.iterrows():
            missing = [c for c in required_cols if pd.isna(row[c])]
            logger.error(
                f"ERROR: Missing values in '{color_table_file}' at line {idx+1}. "
                f"Missing columns: {missing}. Row data: {row.to_dict()}"
            )
        sys.exit(1)

    return df


def parse_formats(out_formats: str) -> list[str]:
    """
    Parses and validates the output file formats from a comma-separated string.

    Args:
        out_formats (str): Comma-separated string of output file formats.

    Returns:
        list[str]: List of valid output file formats. Defaults to ['png'] if no valid formats are provided.

    The function splits the input string, trims whitespace, and converts each format to lowercase.
    It validates each format against a list of accepted formats and logs warnings for unrecognized formats.
    """
    list_of_formats: list[str] = [format.strip().lower()
                                  for format in out_formats.split(',')]
    accepted_formats: list[str] = ["svg", "png", "eps", "ps", "pdf"]

    # Create a new list to hold only the accepted formats
    accepted_list_of_formats = []

    for format in list_of_formats:
        if format in accepted_formats:
            accepted_list_of_formats.append(format)
        else:
            logger.warning(f"WARNING: Unaccepted/unrecognized output file format: {format}")

    # If no valid formats are found, default to 'png'
    if not accepted_list_of_formats:
        logger.warning("WARNING: No valid output file format was specified; generate a PNG file.")
        accepted_list_of_formats: list[str] = ["png"]
    # Remove duplicates if exist
    accepted_list_of_formats = list(set(accepted_list_of_formats))
    return accepted_list_of_formats


def save_figure(canvas: Drawing, list_of_formats: List[str]) -> None:
    """
    Saves the rendered figure to files in the specified formats.

    Args:
        canvas: Canvas object containing the rendered figure.
        list_of_formats (List[str]): List of file formats to save the figure in.

    This function saves the figure in each of the specified formats. It supports 'svg', 'png', 'eps', 'ps', and 'pdf'.
    It logs information about each generated file. If no valid format is specified, it defaults to 'png'.
    """
    # Extract the base filename
    list_of_formats_copy = list_of_formats.copy()
    outfile_base: str = os.path.splitext(canvas.filename)[0]
    svg_data: str = canvas.tostring()
    outfile_count: int = 0
    # Define accepted formats and the conversion function for each
    accepted_formats = {
        "png": cairosvg.svg2png,
        "eps": cairosvg.svg2eps,
        "ps": cairosvg.svg2ps,
        "pdf": cairosvg.svg2pdf
    }
    # Handling svg format separately
    if "svg" in list_of_formats_copy:
        canvas.save()
        logger.info(f"Generated {canvas.filename}!")
        outfile_count += 1
        # Remove 'svg' to avoid processing it again
        list_of_formats_copy.remove('svg')

    for out_format in list_of_formats_copy:
        if out_format in accepted_formats:
            # Call the appropriate Cairo conversion function
            conversion_func = accepted_formats[out_format]
            if conversion_func:  # Ensure conversion function is not None
                out_filename: str = f'{outfile_base}.{out_format}'
                conversion_func(bytestring=svg_data, write_to=out_filename, dpi=96, scale=3)
                logger.info(f"Generated {out_filename}!")
                outfile_count += 1
            else:
                logger.warning(f"WARNING: Conversion function for {out_format} is not defined.")
        else:
            logger.warning(f"WARNING: Unaccepted/unrecognized output file format: {out_format}")

    # Default to PNG if no valid format is specified
    if outfile_count < 1:
        logger.warning("WARNING: No valid output file format was specified. By default, generating a PNG file.")
        out_filename = f'{outfile_base}.png'
        conversion_func = accepted_formats["png"]
        conversion_func(bytestring=svg_data, write_to=out_filename, dpi=96, scale=3)
        logger.info(f"Generated {out_filename}!")


def load_config_toml(config_directory: str, config_file: str) -> dict:
    """
    Loads configuration settings from a TOML file.

    Returns:
        dict: Dictionary containing the loaded configuration settings.

    This function attempts to load configuration settings from a TOML file located within the package.
    If the file is not found or an error occurs during loading, it logs the error and returns an empty dictionary.
    """
    config_dict = {}  # Initialize to empty dict to handle cases where loading fails.
    absolute_config_path = None  # Initialize outside of try for scope in exception block
    try:
        # Generate the path object for the 'config.toml' file
        config_path: Traversable = resources.files(config_directory).joinpath(config_file)
        # Convert the path to an absolute path
        absolute_config_path = config_path.resolve()  # type: ignore
        # Display or log the absolute path
        logger.info(f"INFO: Loading config file: {absolute_config_path}")
        # Open the file and load the configuration
        with open(absolute_config_path, 'rb') as config_toml:
            config_dict: dict = tomllib.load(config_toml)
    except FileNotFoundError as e:
        logger.error(f"Failed to load configs from {absolute_config_path}: {e}")
    return config_dict



def read_color_table(color_table_file: str) -> Optional[DataFrame]:
    """
    Reads a color table from a TSV file and errors out immediately if any row has
    the wrong number of fields or missing values.

    Expected columns: feature_type, qualifier_key, value, color, caption

    Returns:
        DataFrame with the user-defined color mappings, or None if no file was provided.
    """
    required_cols = ['feature_type', 'qualifier_key', 'value', 'color', 'caption']

    # If user did not supply -t, just skip and return None
    if not color_table_file:
        return None

    try:
        df = pd.read_csv(
            color_table_file,
            sep='\t',
            header=None,
            names=required_cols,
            dtype=str,
            on_bad_lines='error',  # raise on any row with wrong number of fields
            engine='python'        # required for on_bad_lines
        )
    except pd.errors.ParserError as e:
        logger.error(f"ERROR: Malformed line in '{color_table_file}': {e}")
        sys.exit(1)
    except FileNotFoundError as e:
        logger.error(f"ERROR: Color table file not found: {e}")
        sys.exit(1)
    except Exception as e:
        logger.error(f"ERROR: Failed to read '{color_table_file}': {e}")
        sys.exit(1)

    # Check for any rows with missing values and error out if found
    null_rows = df[df.isnull().any(axis=1)]
    if not null_rows.empty:
        for idx, row in null_rows.iterrows():
            missing = [c for c in required_cols if pd.isna(row[c])]
            logger.error(
                f"ERROR: Missing values in '{color_table_file}' at line {idx+1}. "
                f"Missing columns: {missing}. Row data: {row.to_dict()}"
            )
        sys.exit(1)

    return df

def read_qualifier_priority_file(filepath: str) -> Optional[DataFrame]:
    """
    Reads a qualifier priority file (TSV) and returns a DataFrame.
    Errors out immediately if any row has the wrong number of fields or missing values.

    Expected columns: feature_type, priorities (comma-separated string)

    Returns:
        DataFrame with the user-defined qualifier priorities, or None if no file was provided.
    """
    # 
    if not filepath:
        return None

    required_cols = ['feature_type', 'priorities']

    try:
        df = pd.read_csv(
            filepath,
            sep='\t',
            header=None,
            names=required_cols,
            dtype=str,
            on_bad_lines='error',  
            engine='python'
        )
    except pd.errors.ParserError as e:
        logger.error(f"ERROR: Malformed line in qualifier priority file '{filepath}': {e}")
        sys.exit(1)
    except FileNotFoundError as e:
        logger.error(f"ERROR: Qualifier priority file not found: {e}")
        sys.exit(1)
    except Exception as e:
        logger.error(f"ERROR: Failed to read '{filepath}': {e}")
        sys.exit(1)

    null_rows = df[df.isnull().any(axis=1)]
    if not null_rows.empty:
        for idx, row in null_rows.iterrows():
            missing = [c for c in required_cols if pd.isna(row[c])]
            logger.error(
                f"ERROR: Missing values in '{filepath}' at line {idx+1}. "
                f"Missing columns: {missing}. Row data: {row.to_dict()}"
            )
        sys.exit(1)

    logger.info(f"Successfully loaded qualifier priority from {filepath}")
    return df