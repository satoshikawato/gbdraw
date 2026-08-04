#!/usr/bin/env python3
"""Execute canonical CLI documentation recipes in clean directories."""

from __future__ import annotations

import argparse
import gzip
import json
import math
import os
import re
import shlex
import shutil
import struct
import subprocess
import zlib
from collections import Counter
from contextlib import ExitStack
from pathlib import Path
from tempfile import TemporaryDirectory
from xml.etree import ElementTree

if __package__:
    from ._scenario_support import (
        PUBLISHED_IMAGE_ROOT,
        RecipeContractError,
        assert_exact_workdir_files,
        copy_declared_inputs,
        extract_executable_block,
        inspect_standard_svg,
        load_chapter,
        parse_translate_chain,
        publish_output,
        validate_standard_svg,
    )
else:
    from _scenario_support import (
        PUBLISHED_IMAGE_ROOT,
        RecipeContractError,
        assert_exact_workdir_files,
        copy_declared_inputs,
        extract_executable_block,
        inspect_standard_svg,
        load_chapter,
        parse_translate_chain,
        publish_output,
        validate_standard_svg,
    )


IMPLEMENTED_SCENARIOS = (
    "T-CLI-01",
    "T-CLI-02",
    "T-CLI-03",
    "T-CLI-05",
    "T-CLI-06",
    "T-CLI-07",
    "T-CLI-08",
    "T-CLI-09",
    "T-CLI-10",
    "T-CLI-11",
    "H-CLI-01",
    "H-CLI-02",
    "H-CLI-03",
    "H-CLI-04",
    "H-CLI-05",
    "H-CLI-06",
    "H-CLI-07",
    "H-CLI-08",
    "H-CLI-09",
    "H-CLI-10",
    "H-CLI-11",
    "H-CLI-12",
    "H-CLI-13",
)
RUNNER_PATH = "docs/recipes/run_cli_scenarios.py"
_PATH_COORDINATE_RE = re.compile(r"[ML]\s+(-?[0-9.]+),(-?[0-9.]+)")
_BGC_RECORD_IDS = (
    "BGC0000708",
    "BGC0000709",
    "BGC0000711",
    "BGC0000712",
    "BGC0000713",
)
_BGC_ADJACENT_PAIRS = (
    ("BGC0000708", "BGC0000709"),
    ("BGC0000709", "BGC0000711"),
    ("BGC0000711", "BGC0000712"),
    ("BGC0000712", "BGC0000713"),
)
_METAZOAN_CONSERVATION_RINGS = (
    (
        "danio-human.tlosatx.tsv",
        "NC_002333.2",
        "Danio rerio (NC_002333.2)",
        "#4e79a7",
        68,
    ),
    (
        "drosophila-human.tlosatx.tsv",
        "NC_024511.2",
        "Drosophila melanogaster (NC_024511.2)",
        "#f28e2b",
        24,
    ),
    (
        "caenorhabditis-human.tlosatx.tsv",
        "NC_001328.1",
        "Caenorhabditis elegans (NC_001328.1)",
        "#59a14f",
        14,
    ),
)
_HUMAN_MT_CDS_GENES = {
    "ND1",
    "ND2",
    "COX1",
    "COX2",
    "ATP8",
    "ATP6",
    "COX3",
    "ND3",
    "ND4L",
    "ND4",
    "ND5",
    "ND6",
    "CYTB",
}
_HUMAN_MT_CDS_PRODUCTS = {
    "NADH dehydrogenase subunit 1",
    "NADH dehydrogenase subunit 2",
    "cytochrome c oxidase subunit I",
    "cytochrome c oxidase subunit II",
    "ATP synthase F0 subunit 8",
    "ATP synthase F0 subunit 6",
    "cytochrome c oxidase subunit III",
    "NADH dehydrogenase subunit 3",
    "NADH dehydrogenase subunit 4L",
    "NADH dehydrogenase subunit 4",
    "NADH dehydrogenase subunit 5",
    "NADH dehydrogenase subunit 6",
    "cytochrome b",
}

_FEATURE_PRESENTATION_TABLES = {
    "tables/presentation_colors.tsv": (
        "CDS\tgene\t^ND[1-6]$\t#3B82F6\tNADH dehydrogenase\n"
        "CDS\tgene\t^COX[1-3]$\t#EF4444\tCytochrome c oxidase\n"
        "CDS\tgene\t^ATP[68]$\t#F59E0B\tATP synthase\n"
        "CDS\tgene\t^CYTB$\t#8B5CF6\tCytochrome b\n"
        "rRNA\tgene\t^RNR[12]$\t#10B981\tRibosomal RNA\n"
    ),
    "tables/presentation_labels.tsv": (
        "CDS\tgene\t^(ND1|COX2|ATP8|ATP6|COX3|CYTB)$\n"
        "rRNA\tgene\t^RNR[12]$\n"
    ),
    "tables/presentation_label_overrides.tsv": (
        "record_id\tfeature_type\tqualifier\tvalue\tlabel_text\n"
        "NC_012920.1\tCDS\tlabel\t^ND1$\tComplex I (ND1)\n"
        "NC_012920.1\tCDS\tlabel\t^COX2$\tOxidase II\n"
        "NC_012920.1\trRNA\tlabel\t^s-rRNA$\t12S rRNA\n"
        "NC_012920.1\trRNA\tlabel\t^l-rRNA$\t16S rRNA\n"
    ),
}

_TUTORIAL_FEATURE_PRESENTATION_TABLES = {
    **_FEATURE_PRESENTATION_TABLES,
    "tables/mitochondrial_regions.tsv": (
        "set_id\tid\tmark\trecord\tstart\tend\tcoordinate_space\twraps_origin\tlabel\tlane\tstroke\tstroke_width\tline_cap\tlabel_color\tlabel_font_size\tlabel_orientation\tlabel_offset\n"
        "mitochondrial_regions\td_loop\tbracket\tNC_012920.1\t16024\t576\tsource\ttrue\tD-loop\t0\t#202020\t3\ttick\t#202020\t14\ttangent\t7\n"
    ),
}

_GENERATED_TABLES = {
    "H-CLI-02": {
        "tables/records.tsv": (
            "gbk\trecord_label\trecord_subtitle\trecord_id\torder\trow\tcolumn\n"
            "../HmmtDNA.gbk\tHuman mitochondrion\tComplete RefSeq record\tNC_012920.1\t1\t1\t1\n"
            "../NC_002333.2.gb\tZebrafish mitochondrion\tComplete RefSeq record\tNC_002333.2\t2\t1\t2\n"
            "../NC_024511.2.gb\tFruit fly mitochondrion\tComplete RefSeq record\tNC_024511.2\t3\t2\t1\n"
            "../NC_001328.1.gb\tNematode mitochondrion\tComplete RefSeq record\tNC_001328.1\t4\t2\t2\n"
        ),
        "tables/comparisons.tsv": (
            "blast\tquery\tsubject\n"
            "../lambda-de3.losatn.tsv\tNC_001416.1\tNC_042057.1\n"
        ),
        "tables/conservation.tsv": (
            "blast\tlabel\tcolor\n"
            "../danio-human.tlosatx.tsv\tDanio rerio (NC_002333.2)\t#4E79A7\n"
            "../drosophila-human.tlosatx.tsv\tDrosophila melanogaster (NC_024511.2)\t#F28E2B\n"
            "../caenorhabditis-human.tlosatx.tsv\tCaenorhabditis elegans (NC_001328.1)\t#59A14F\n"
        ),
        "tables/tracks.tsv": (
            "id\trenderer\tside\tr\tw\tinner_gap_px\touter_gap_px\tz\tparams\n"
            "ticks\tticks\toutside\t\t\t\t\t\ttick_label_layout=label_out_tick_in\n"
            "features\tfeatures\taxis\t\t\t\t\t\t\n"
            "gc_content\tdinucleotide_content\tinside\t\t0.10\t3\t3\t\tnt=GC,legend_label=GC content\n"
            "gc_skew\tdinucleotide_skew\tinside\t\t0.10\t3\t3\t\tnt=GC,legend_label=GC skew\n"
            "at_skew\tdinucleotide_skew\tinside\t\t0.10\t3\t3\t\tnt=AT,positive_color=#deaf6e,negative_color=#7294e3,legend_label=AT skew\n"
        ),
    },
    "H-CLI-11": _FEATURE_PRESENTATION_TABLES,
    "T-CLI-03": _TUTORIAL_FEATURE_PRESENTATION_TABLES,
}

_SCENARIO_REQUIRED_FIXTURE_FILES = {
    "H-CLI-02": (
        "HmmtDNA.gbk",
        "NC_002333.2.gb",
        "NC_024511.2.gb",
        "NC_001328.1.gb",
        "lambda-de3.losatn.tsv",
        "danio-human.tlosatx.tsv",
        "drosophila-human.tlosatx.tsv",
        "caenorhabditis-human.tlosatx.tsv",
    ),
}

_OUTPUT_INPUT_FILES = {
    ("T-CLI-11", "restored_interactive_figure.svg"): {
        "HmmtDNA.gbk",
        "cds_gene_qualifier_priority.tsv",
    },
    ("T-CLI-03", "mitochondrial_features_baseline.svg"): {
        "HmmtDNA.gbk",
    },
    ("T-CLI-03", "mitochondrial_features_highlighted.svg"): {
        "HmmtDNA.gbk",
        "cds_gene_qualifier_priority.tsv",
    },
    ("T-CLI-05", "quantitative_genome_baseline.svg"): {
        "AP027133.gb",
    },
    ("T-CLI-05", "quantitative_genome_map.svg"): {
        "AP027133.gb",
        "AP027133.DRR394922.depth-1kb.tsv",
    },
    ("H-CLI-02", "record_table.svg"): {
        "HmmtDNA.gbk",
        "NC_002333.2.gb",
        "NC_024511.2.gb",
        "NC_001328.1.gb",
        "cds_gene_qualifier_priority.tsv",
    },
    ("H-CLI-02", "comparison_table.svg"): {
        "NC_001416.gb",
        "NC_042057.1.gb",
        "lambda-de3.losatn.tsv",
    },
    ("H-CLI-02", "conservation_table.svg"): {
        "HmmtDNA.gbk",
        "danio-human.tlosatx.tsv",
        "drosophila-human.tlosatx.tsv",
        "caenorhabditis-human.tlosatx.tsv",
        "cds_gene_qualifier_priority.tsv",
    },
    ("H-CLI-02", "annotation_table.svg"): {
        "NC_001879.gbk",
        "nicotiana-tabacum-regions.tsv",
        "modified_default_colors.tsv",
        "qualifier_priority.tsv",
    },
    ("H-CLI-02", "track_table.svg"): {
        "HmmtDNA.gbk",
        "cds_gene_qualifier_priority.tsv",
    },
    ("H-CLI-09", "cli_quantitative_tracks.svg"): {
        "AP027133.gb",
        "AP027133.DRR394922.depth-1kb.tsv",
    },
    ("H-CLI-10", "cli_annotations_slots.svg"): {
        "NC_001879.gbk",
        "nicotiana-tabacum-regions.tsv",
        "modified_default_colors.tsv",
        "qualifier_priority.tsv",
    },
    ("H-CLI-11", "cli_feature_presentation.svg"): {
        "HmmtDNA.gbk",
        "HmmtDNA_feature_visibility.tsv",
        "cds_gene_qualifier_priority.tsv",
    },
    ("H-CLI-12", "cli_session_roundtrip.svg"): {
        "HmmtDNA.gbk",
        "cds_gene_qualifier_priority.tsv",
    },
    ("H-CLI-13", "cli_export.svg"): {
        "HmmtDNA.gbk",
        "cds_gene_qualifier_priority.tsv",
    },
}


def _parse_commands(recipe: str, *, scenario_id: str) -> list[list[str]]:
    commands = [
        shlex.split(line)
        for line in recipe.replace("\\\n", " ").splitlines()
        if line.strip()
    ]
    if not commands:
        raise RecipeContractError(f"{scenario_id} has no command to execute.")
    for command in commands:
        if command[:1] != ["gbdraw"] or command[1:2] not in (["circular"], ["linear"]):
            raise RecipeContractError(
                f"{scenario_id} must call a documented gbdraw diagram command."
            )
        if any(token in {";", "&&", "||", "|", ">", ">>"} for token in command):
            raise RecipeContractError(f"{scenario_id} must not use a shell pipeline.")
    return commands


def _command_entries(
    command: list[str],
    *,
    used_entries: list[dict[str, object]],
) -> list[dict[str, object]]:
    filenames = set(command)
    return [
        entry
        for entry in used_entries
        if Path(str(entry["relativePath"])).name in filenames
    ]


def _output_entries(
    scenario_id: str,
    output_name: str,
    *,
    command: list[str],
    used_entries: list[dict[str, object]],
) -> list[dict[str, object]]:
    filenames = _OUTPUT_INPUT_FILES.get((scenario_id, output_name))
    if filenames is None:
        return _command_entries(command, used_entries=used_entries)
    return [
        entry
        for entry in used_entries
        if Path(str(entry["relativePath"])).name in filenames
    ]


def _copy_recipe_source(scenario_id: str, recipe: str) -> str:
    source = recipe
    for relative_path in _GENERATED_TABLES.get(scenario_id, {}):
        source = source.replace(relative_path, "")
    required = _SCENARIO_REQUIRED_FIXTURE_FILES.get(scenario_id, ())
    if required:
        source += "\n" + " ".join(required) + "\n"
    return source


def _materialize_generated_tables(scenario_id: str, workdir: Path) -> set[str]:
    generated: set[str] = set()
    for relative_path, source in _GENERATED_TABLES.get(scenario_id, {}).items():
        destination = workdir / relative_path
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_text(source, encoding="utf-8", newline="\n")
        generated.add(relative_path)
    return generated


def _artifact_comparison_payload(scenario_id: str, path: Path) -> bytes:
    payload = path.read_bytes()
    normalized_session_names = {
        "H-CLI-12": {"cli_session.json", "cli_session.json.gz"},
        "T-CLI-11": {"interactive_handoff.gbdraw-session.json.gz"},
    }
    if path.name in normalized_session_names.get(scenario_id, set()):
        session = json.loads(gzip.decompress(payload) if path.suffix == ".gz" else payload)

        def without_file_times(value: object) -> object:
            if isinstance(value, dict):
                return {
                    key: without_file_times(item)
                    for key, item in value.items()
                    if key not in {"createdAt", "lastModified"}
                }
            if isinstance(value, list):
                return [without_file_times(item) for item in value]
            return value

        return json.dumps(
            without_file_times(session),
            sort_keys=True,
            separators=(",", ":"),
        ).encode("utf-8")
    if scenario_id == "H-CLI-13" and path.suffix == ".pdf":
        object_pattern = re.compile(
            rb"(?ms)^(\d+) 0 obj\r?\n(.*?)\r?\nendobj$",
        )
        objects = object_pattern.findall(payload)
        length_objects = {
            match.group(1)
            for _, body in objects
            for match in re.finditer(rb"/Length (\d+) 0 R", body)
        }
        normalized_objects: list[bytes] = []
        for object_number, body in objects:
            if object_number in length_objects or b"/Type /XRef" in body:
                continue
            stream = re.fullmatch(
                rb"(.*?)\r?\nstream\r?\n(.*?)\r?\nendstream",
                body,
                flags=re.DOTALL,
            )
            if stream is not None:
                try:
                    decoded = zlib.decompress(stream.group(2))
                except zlib.error:
                    decoded = stream.group(2)
                body = stream.group(1) + b"\nstream\n" + decoded + b"\nendstream"
            body = re.sub(
                rb"/CreationDate \(D:\d{14}Z\)",
                b"/CreationDate (normalized)",
                body,
            )
            normalized_objects.append(
                object_number + b" 0 obj\n" + body + b"\nendobj"
            )
        return b"\n".join(normalized_objects)
    if scenario_id == "H-CLI-13" and path.suffix in {".eps", ".ps"}:
        return re.sub(
            rb"^%%CreationDate:.*$",
            b"%%CreationDate: normalized",
            payload,
            flags=re.MULTILINE,
        )
    return payload


def _assert_gff_id_mismatch_is_rejected(
    *,
    command: list[str],
    workdir: Path,
    environment: dict[str, str],
) -> None:
    fasta_name = "NC_001416.fna"
    if fasta_name not in command:
        raise RecipeContractError("H-CLI-01 is missing its GFF3/FASTA command.")
    source = (workdir / fasta_name).read_text(encoding="utf-8")
    mismatched_name = "mismatched_sequence_id.fna"
    (workdir / mismatched_name).write_text(
        source.replace(">NC_001416.1", ">MISMATCHED_ID", 1),
        encoding="utf-8",
    )
    failure_command = [
        mismatched_name if token == fasta_name else token
        for token in command
    ]
    output_index = failure_command.index("-o") + 1
    failure_command[output_index] = "must_not_exist"
    result = subprocess.run(
        failure_command,
        cwd=workdir,
        env=environment,
        capture_output=True,
        text=True,
        timeout=120,
        check=False,
    )
    evidence = result.stdout + result.stderr
    expected = "No matching FASTA record found for GFF record NC_001416.1"
    if result.returncode == 0 or expected not in evidence:
        raise RecipeContractError(
            "H-CLI-01 did not reject a mismatched GFF3/FASTA sequence ID "
            f"with the documented message.\nstdout:\n{result.stdout}\nstderr:\n{result.stderr}"
        )
    if (workdir / "must_not_exist.svg").exists():
        raise RecipeContractError("H-CLI-01 wrote output after an input-ID failure.")


def _assert_multi_record_circular_layout(
    output_path: Path,
    *,
    used_entries: list[dict[str, object]],
) -> None:
    expected_ids = [
        "NC_012920.1",
        "NC_002333.2",
        "NC_024511.2",
        "NC_001328.1",
    ]
    expected_lengths = {
        "NC_012920.1": 16_569,
        "NC_002333.2": 16_596,
        "NC_024511.2": 19_524,
        "NC_001328.1": 13_794,
    }
    expected_feature_counts = {
        "NC_012920.1": 37,
        "NC_002333.2": 37,
        "NC_024511.2": 37,
        "NC_001328.1": 36,
    }
    source_records = {
        str(record["id"]): (entry, record)
        for entry in used_entries
        if entry.get("inputType") == "genbank"
        for record in entry.get("records", ())
    }
    if set(source_records) != set(expected_ids):
        raise RecipeContractError("H-CLI-03 has the wrong source records.")
    for record_id in expected_ids:
        entry, record = source_records[record_id]
        if (
            entry.get("role") != "raw"
            or entry.get("derivation") is not None
            or record.get("topology") != "circular"
            or record.get("length") != expected_lengths[record_id]
            or not str(record.get("description", "")).endswith("complete genome")
        ):
            raise RecipeContractError(
                "H-CLI-03 requires complete, circular, unmodified source records."
            )

    root = ElementTree.parse(output_path).getroot()
    record_groups = [
        element
        for element in root
        if element.tag.rsplit("}", 1)[-1] == "g"
        and element.attrib.get("id", "").startswith("record_")
        and "data-gbdraw-record-id" in element.attrib
    ]
    record_groups.sort(key=lambda element: int(element.attrib["data-gbdraw-record-index"]))
    if [
        element.attrib["data-gbdraw-record-id"] for element in record_groups
    ] != expected_ids:
        raise RecipeContractError("H-CLI-03 output has the wrong record order.")

    positions: list[tuple[float, float]] = []
    expected_slots = ["ticks", "features", "gc_content", "gc_skew"]
    feature_counts: dict[str, int] = {}
    for group in record_groups:
        translation = parse_translate_chain(group.attrib.get("transform", ""))
        if translation is None:
            raise RecipeContractError("H-CLI-03 record placement metadata is missing.")
        positions.append(translation)
        slots = [
            child.attrib["data-gbdraw-slot-id"]
            for child in group
            if "data-gbdraw-slot-id" in child.attrib
        ]
        if slots != expected_slots:
            raise RecipeContractError("H-CLI-03 output has the wrong track order.")
        axes = [
            child
            for child in group
            if child.attrib.get("id", "").startswith("Axis_")
        ]
        if len(axes) != 1:
            raise RecipeContractError("H-CLI-03 output lacks one axis per record.")
        record_id = group.attrib["data-gbdraw-record-id"]
        feature_counts[record_id] = len(
            {
                element.attrib["data-gbdraw-feature-id"]
                for element in group.iter()
                if "data-gbdraw-feature-id" in element.attrib
            }
        )
    if feature_counts != expected_feature_counts:
        raise RecipeContractError("H-CLI-03 output has the wrong feature counts.")
    _assert_multi_record_cds_gene_labels(
        record_groups,
        source_records=source_records,
    )
    if not _positions_form_documented_2x2_grid(positions):
        raise RecipeContractError("H-CLI-03 records do not form the documented 2x2 grid.")


def _positions_form_documented_2x2_grid(
    positions: list[tuple[float, float]],
    *,
    row_origin_tolerance: float = 25.0,
    minimum_separation: float = 100.0,
) -> bool:
    if len(positions) != 4:
        return False
    top_row = positions[:2]
    bottom_row = positions[2:]
    rows = (top_row, bottom_row)
    if any(abs(left[1] - right[1]) > row_origin_tolerance for left, right in rows):
        return False
    if min(position[1] for position in bottom_row) - max(
        position[1] for position in top_row
    ) < minimum_separation:
        return False
    return all(
        right[0] - left[0] >= minimum_separation
        for left, right in rows
    )


def _assert_multi_record_cds_gene_labels(
    record_groups: list[ElementTree.Element],
    *,
    source_records: dict[str, tuple[dict[str, object], dict[str, object]]],
) -> None:
    from Bio import SeqIO

    fixture_root = (
        Path(__file__).resolve().parents[2] / "gbdraw" / "web" / "tutorial-data"
    )
    for group in record_groups:
        record_id = group.attrib["data-gbdraw-record-id"]
        entry, _record = source_records[record_id]
        source = SeqIO.read(fixture_root / str(entry["relativePath"]), "genbank")
        genes = {
            feature.qualifiers["gene"][0]
            for feature in source.features
            if feature.type == "CDS" and feature.qualifiers.get("gene")
        }
        products = {
            feature.qualifiers["product"][0]
            for feature in source.features
            if feature.type == "CDS" and feature.qualifiers.get("product")
        }
        labels = {
            "".join(element.itertext()).strip()
            for element in group.iter()
            if element.tag.rsplit("}", 1)[-1] == "text"
        }
        if not genes or not genes <= labels:
            raise RecipeContractError(
                f"Multi-record CLI output is missing CDS gene labels for {record_id}."
            )
        product_labels = (products - genes) & labels
        if product_labels:
            raise RecipeContractError(
                "Multi-record CLI output used CDS product text as labels for "
                f"{record_id}: {sorted(product_labels)}"
            )


def _assert_linear_regions_orientation_layout(
    chapter: dict[str, object], output_path: Path
) -> None:
    evidence = inspect_standard_svg(chapter, output_path=output_path)
    expected_ids = ["NC_001416.1", "BGC0000708", "BGC0000713"]
    if evidence.record_ids != frozenset(expected_ids):
        raise RecipeContractError("H-CLI-04 output has the wrong records.")

    root = ElementTree.parse(output_path).getroot()
    record_groups = [
        element
        for element in root
        if element.tag.rsplit("}", 1)[-1] == "g"
        and element.attrib.get("id", "").startswith("record_group_")
        and "data-gbdraw-record-id" in element.attrib
    ]
    record_groups.sort(key=lambda element: int(element.attrib["data-gbdraw-record-index"]))
    if [
        element.attrib["data-gbdraw-record-id"] for element in record_groups
    ] != expected_ids:
        raise RecipeContractError("H-CLI-04 record order metadata is wrong.")

    positions: list[tuple[float, float]] = []
    feature_counts: dict[str, int] = {}
    for group in record_groups:
        translation = parse_translate_chain(group.attrib.get("transform", ""))
        if translation is None:
            raise RecipeContractError("H-CLI-04 record placement metadata is missing.")
        positions.append(translation)
        record_id = group.attrib["data-gbdraw-record-id"]
        feature_counts[record_id] = len(
            {
                element.attrib["data-gbdraw-feature-id"]
                for element in group.iter()
                if "data-gbdraw-feature-id" in element.attrib
            }
        )

    if not (
        positions[0][1] < positions[1][1]
        and positions[1][1] == positions[2][1]
        and positions[1][0] < positions[2][0]
    ):
        raise RecipeContractError("H-CLI-04 row topology is wrong.")
    if feature_counts != {"NC_001416.1": 43, "BGC0000708": 30, "BGC0000713": 26}:
        raise RecipeContractError("H-CLI-04 cropped/full-record feature counts are wrong.")

    if "5001-35500" not in evidence.text_nodes:
        raise RecipeContractError("H-CLI-04 Lambda coordinate range is wrong.")

    expected_ticks = [
        "5 kbp",
        "10 kbp",
        "15 kbp",
        "20 kbp",
        "25 kbp",
        "30 kbp",
    ]
    for record_id, group in zip(expected_ids, record_groups, strict=True):
        if record_id not in {"NC_001416.1", "BGC0000713"}:
            continue
        ticks = [
            ("".join(element.itertext()).strip(), float(element.attrib["x"]))
            for element in group.iter()
            if element.tag.rsplit("}", 1)[-1] == "text"
            and "".join(element.itertext()).strip().endswith("kbp")
            and "x" in element.attrib
        ]
        if [label for label, _x in ticks] != expected_ticks:
            raise RecipeContractError(f"H-CLI-04 {record_id} ruler labels are wrong.")
        if any(
            left_x >= right_x
            for (_, left_x), (_, right_x) in zip(ticks, ticks[1:])
        ):
            raise RecipeContractError(
                f"H-CLI-04 {record_id} ruler does not ascend from left to right."
            )

    bgc_0713_group = record_groups[2]
    gene_x = {
        "".join(element.itertext()).strip(): float(element.attrib["x"])
        for element in bgc_0713_group.iter()
        if element.tag.rsplit("}", 1)[-1] == "text"
        and "".join(element.itertext()).strip() in {"racG", "racP"}
        and "x" in element.attrib
    }
    if set(gene_x) != {"racG", "racP"} or gene_x["racG"] >= gene_x["racP"]:
        raise RecipeContractError("H-CLI-04 reverse-complement feature order is wrong.")
    required_text = {
        "Lambda selected region",
        "NC_001416.1 positions 5,001–35,500",
        "Lividomycin cluster",
        "Complete BGC0000708 record",
        "Ribostamycin cluster",
        "Complete BGC0000713 reverse complement",
        "40,579 bp",
        "31,892 bp",
    }
    if not required_text <= evidence.text_nodes:
        raise RecipeContractError("H-CLI-04 record labels or subtitles are missing.")


def _match_elements(root: ElementTree.Element) -> list[ElementTree.Element]:
    return [
        element
        for element in root.iter()
        if "data-gbdraw-match-id" in element.attrib
    ]


def _assert_complete_metazoan_mtdna_conservation(
    output_path: Path,
    *,
    title: str = "Complete metazoan mitochondrial TLOSATX evidence",
) -> None:
    root = ElementTree.parse(output_path).getroot()
    record_ids = {
        element.attrib["data-gbdraw-record-id"]
        for element in root.iter()
        if "data-gbdraw-record-id" in element.attrib
    }
    feature_ids = {
        element.attrib["data-gbdraw-feature-id"]
        for element in root.iter()
        if "data-gbdraw-feature-id" in element.attrib
    }
    if record_ids != {"NC_012920.1"} or len(feature_ids) != 37:
        raise RecipeContractError(
            "Circular mtDNA conservation must display the complete human record."
        )

    expected_by_query = {
        query_id: (label, color, count)
        for _, query_id, label, color, count in _METAZOAN_CONSERVATION_RINGS
    }
    matches = _match_elements(root)
    if len(matches) != sum(item[4] for item in _METAZOAN_CONSERVATION_RINGS):
        raise RecipeContractError("Circular mtDNA conservation hit count changed.")
    counts: Counter[str] = Counter()
    subject_intervals: list[tuple[int, int]] = []
    for element in matches:
        attributes = element.attrib
        query_id = attributes.get("data-query-record-id", "")
        expected = expected_by_query.get(query_id)
        if expected is None:
            raise RecipeContractError("Circular mtDNA conservation query changed.")
        label, color, _ = expected
        if (
            attributes.get("data-subject-record-id") != "NC_012920.1"
            or attributes.get("data-reference-record-id") != "NC_012920.1"
            or attributes.get("data-reference-side") != "subject"
            or attributes.get("data-match-kind") != "homology"
            or attributes.get("data-track-label") != label
            or attributes.get("data-track-color", "").lower() != color
        ):
            raise RecipeContractError(
                "Circular mtDNA conservation endpoint or ring metadata changed."
            )
        counts[query_id] += 1
        start = int(float(attributes["data-sstart"]))
        end = int(float(attributes["data-send"]))
        subject_intervals.append((min(start, end), max(start, end)))
    if counts != Counter(
        {query_id: count for _, query_id, _, _, count in _METAZOAN_CONSERVATION_RINGS}
    ):
        raise RecipeContractError("Circular mtDNA conservation ring counts changed.")

    covered = 0
    current_start = current_end = None
    for start, end in sorted(subject_intervals):
        if current_start is None:
            current_start, current_end = start, end
        elif start > current_end + 1:
            covered += current_end - current_start + 1
            current_start, current_end = start, end
        else:
            current_end = max(current_end, end)
    if current_start is not None:
        covered += current_end - current_start + 1
    if covered != 9_813:
        raise RecipeContractError("Circular mtDNA conservation coverage changed.")

    slots = [
        (
            element.attrib["data-gbdraw-slot-id"],
            element.attrib["data-gbdraw-slot-renderer"],
        )
        for element in root.iter()
        if "data-gbdraw-slot-id" in element.attrib
    ]
    if slots != [
        ("features", "features"),
        ("conservation_1", "sequence_conservation"),
        ("conservation_2", "sequence_conservation"),
        ("conservation_3", "sequence_conservation"),
        ("ticks", "ticks"),
        ("gc_content", "dinucleotide_content"),
        ("gc_skew", "dinucleotide_skew"),
    ]:
        raise RecipeContractError("Circular mtDNA conservation track order changed.")

    text_nodes = _text_nodes(root)
    required_text = _HUMAN_MT_CDS_GENES | {
        item[2] for item in _METAZOAN_CONSERVATION_RINGS
    } | {
        "NC_012920.1",
        "16,569 bp",
        title,
    }
    if not required_text <= text_nodes:
        raise RecipeContractError(
            "Human mtDNA labels or conservation-ring legend text are missing."
        )
    if _HUMAN_MT_CDS_PRODUCTS & text_nodes:
        raise RecipeContractError("Human mtDNA CDS product labels leaked into output.")


def _assert_precomputed_comparison(
    chapter: dict[str, object], output_path: Path
) -> None:
    evidence = inspect_standard_svg(chapter, output_path=output_path)
    root = ElementTree.parse(output_path).getroot()
    expected_pair = ("NC_001416.1", "NC_042057.1")
    if output_path.name == "linear_precomputed_comparison.svg":
        matches = _match_elements(root)
        alignment_lengths = {
            int(float(element.attrib["data-alignment-length"])) for element in matches
        }
        expected_lengths = {21_232, 6_412, 5_205, 1_914, 1_620, 254}
        if len(matches) != 6 or alignment_lengths != expected_lengths:
            raise RecipeContractError("H-CLI-05 retained the wrong LOSATN HSP set.")
        if evidence.record_ids != frozenset(expected_pair):
            raise RecipeContractError("H-CLI-05 Linear records are wrong.")
        if len(evidence.feature_ids) != 130:
            raise RecipeContractError("H-CLI-05 Linear feature count is wrong.")
        expected_match_ids = {
            f"comparison1_match{index}" for index in range(1, 7)
        }
        if {
            element.attrib["data-gbdraw-match-id"] for element in matches
        } != expected_match_ids:
            raise RecipeContractError("H-CLI-05 Linear match IDs are wrong.")
        for element in matches:
            attributes = element.attrib
            if (
                attributes.get("data-query-record-id"),
                attributes.get("data-subject-record-id"),
            ) != expected_pair or attributes.get("data-match-kind") != "pairwise":
                raise RecipeContractError(
                    "H-CLI-05 Linear query/subject mapping is wrong."
                )
            if attributes.get("data-pairwise-match-style") != "curve":
                raise RecipeContractError("H-CLI-05 Linear links are not curves.")
        return

    if output_path.name != "circular_conservation_ring.svg":
        raise RecipeContractError(f"Unexpected H-CLI-05 output: {output_path.name}")
    _assert_complete_metazoan_mtdna_conservation(output_path)


def _assert_tutorial_losatn_comparison(
    chapter: dict[str, object], output_path: Path
) -> None:
    evidence = inspect_standard_svg(chapter, output_path=output_path)
    root = ElementTree.parse(output_path).getroot()
    matches = _match_elements(root)
    expected_pair = ("NC_001416.1", "NC_042057.1")
    expected_lengths = {21_232, 6_412, 5_205, 1_914, 1_620, 254}
    if evidence.record_ids != frozenset(expected_pair):
        raise RecipeContractError("T-CLI-07 rendered the wrong complete records.")
    if len(matches) != 6 or {
        int(float(element.attrib["data-alignment-length"])) for element in matches
    } != expected_lengths:
        raise RecipeContractError("T-CLI-07 retained the wrong LOSATN evidence.")
    for element in matches:
        attributes = element.attrib
        if (
            attributes.get("data-query-record-id"),
            attributes.get("data-subject-record-id"),
        ) != expected_pair or attributes.get("data-match-kind") != "pairwise":
            raise RecipeContractError("T-CLI-07 changed its comparison endpoints.")
        if attributes.get("data-pairwise-match-style") != "ribbon":
            raise RecipeContractError("T-CLI-07 must retain the GUI ribbon style.")


def _assert_hepatoplasmataceae_collinear(
    chapter: dict[str, object], output_path: Path
) -> None:
    evidence = inspect_standard_svg(chapter, output_path=output_path)
    expected_records = {
        "AP027078.1",
        "AP027131.1",
        "AP027133.1",
        "AP027132.1",
        "NZ_CP006932.1",
    }
    if evidence.record_ids != expected_records or len(evidence.feature_ids) != 2_994:
        raise RecipeContractError("T-CLI-10 rendered the wrong complete genomes.")
    root = ElementTree.parse(output_path).getroot()
    matches = [
        element
        for element in root.iter()
        if element.attrib.get("data-match-kind") == "collinear"
    ]
    expected_pairs = {
        ("AP027078.1", "AP027131.1"),
        ("AP027131.1", "AP027133.1"),
        ("AP027133.1", "AP027132.1"),
        ("AP027132.1", "NZ_CP006932.1"),
    }
    observed_pairs = {
        (
            element.attrib.get("data-query-record-id"),
            element.attrib.get("data-subject-record-id"),
        )
        for element in matches
    }
    if len(matches) != 500 or observed_pairs != expected_pairs:
        raise RecipeContractError("T-CLI-10 changed its adjacent Collinear evidence.")
    required_text = {
        "LOSATP Collinear blocks across Hepatoplasmataceae",
        "GC content",
        "GC skew (+)",
        "GC skew (-)",
    }
    if not required_text <= evidence.text_nodes:
        raise RecipeContractError("T-CLI-10 lost its title or quantitative tracks.")


def _assert_bgc_svg(
    chapter: dict[str, object], output_path: Path
) -> tuple[ElementTree.Element, list[ElementTree.Element]]:
    evidence = inspect_standard_svg(chapter, output_path=output_path)
    if evidence.record_ids != frozenset(_BGC_RECORD_IDS):
        raise RecipeContractError(f"{chapter['id']} output has the wrong BGC records.")
    if len(evidence.feature_ids) != 155:
        raise RecipeContractError(f"{chapter['id']} expected 155 CDS features.")
    root = ElementTree.parse(output_path).getroot()
    return root, _match_elements(root)


def _assert_pinned_losat(command: list[str], *, scenario_id: str) -> None:
    mode_index = command.index("--protein_blastp_mode") + 1
    threads_index = command.index("--losatp_threads") + 1
    expected_mode = {
        "T-CLI-08": "orthogroup",
        "T-CLI-10": "collinear",
        "H-CLI-06": "pairwise",
        "H-CLI-07": "orthogroup",
        "H-CLI-08": "collinear",
    }[scenario_id]
    expected_threads = "32" if scenario_id == "T-CLI-10" else "1"
    if (
        command[mode_index] != expected_mode
        or command[threads_index] != expected_threads
    ):
        raise RecipeContractError(
            f"{scenario_id} must run its documented LOSATP mode and thread count."
        )
    if "--losatp_bin" in command or "--ncbi_blastp_bin" in command:
        raise RecipeContractError(
            f"{scenario_id} must exercise automatic bundled-runtime selection."
        )

    from gbdraw.analysis.protein_colinearity import (  # noqa: PLC2701
        _resolve_protein_blastp_runtime,
    )

    with ExitStack() as stack:
        runtime = _resolve_protein_blastp_runtime("losat", None, stack)
        if runtime.kind != "losat" or runtime.source != "bundled":
            raise RecipeContractError(
                f"{scenario_id} did not resolve the bundled LOSAT runtime."
            )
        result = subprocess.run(
            [runtime.executable, "--version"],
            capture_output=True,
            text=True,
            timeout=15,
            check=False,
        )
    if result.returncode != 0 or result.stdout.strip() != "losat 0.1.0":
        raise RecipeContractError(
            f"{scenario_id} expected bundled LOSAT 0.1.0, got "
            f"{result.stdout.strip() or result.stderr.strip()!r}."
        )


def _assert_pairwise_protein_search(
    chapter: dict[str, object], *, tsv_path: Path, svg_path: Path
) -> None:
    source = tsv_path.read_text(encoding="utf-8")
    expected_filenames = tuple(
        f"{query}.{subject}.losatp.tsv"
        for query, subject in _BGC_ADJACENT_PAIRS
    )
    section_counts: list[int] = []
    filenames: list[str] = []
    rows: list[list[str]] = []
    for line in source.splitlines():
        entry_match = re.fullmatch(r"# entry ([1-9][0-9]*): (.+)", line)
        if entry_match is not None:
            if int(entry_match.group(1)) != len(filenames) + 1:
                raise RecipeContractError("H-CLI-06 raw entry order is wrong.")
            filenames.append(entry_match.group(2))
            section_counts.append(0)
            continue
        if not line or line.startswith("#"):
            continue
        if not section_counts:
            raise RecipeContractError("H-CLI-06 raw row precedes its entry header.")
        columns = line.split("\t")
        if len(columns) != 12:
            raise RecipeContractError("H-CLI-06 raw TSV is not BLAST outfmt 6.")
        section_counts[-1] += 1
        rows.append(columns)
    if tuple(filenames) != expected_filenames or section_counts != [204, 220, 160, 207]:
        raise RecipeContractError("H-CLI-06 raw pair evidence or row counts changed.")
    if len(rows) != 791 or "h_" in source:
        raise RecipeContractError(
            "H-CLI-06 raw evidence is incomplete or exposes runtime handles."
        )

    root, matches = _assert_bgc_svg(chapter, svg_path)
    del root
    expected_pair_counts = Counter(
        {
            ("BGC0000708", "BGC0000709"): 19,
            ("BGC0000709", "BGC0000711"): 21,
            ("BGC0000711", "BGC0000712"): 19,
            ("BGC0000712", "BGC0000713"): 17,
        }
    )
    pair_counts = Counter(
        (
            element.attrib.get("data-query-record-id"),
            element.attrib.get("data-subject-record-id"),
        )
        for element in matches
    )
    if len(matches) != 76 or pair_counts != expected_pair_counts:
        raise RecipeContractError("H-CLI-06 Pairwise link counts changed.")
    if any(
        element.attrib.get("data-match-kind") != "pairwise"
        or "data-orthogroup-id" in element.attrib
        or "data-collinearity-block-id" in element.attrib
        for element in matches
    ):
        raise RecipeContractError("H-CLI-06 contains grouped or Collinear links.")


def _assert_similarity_groups(
    chapter: dict[str, object], output_path: Path
) -> None:
    root, matches = _assert_bgc_svg(chapter, output_path)
    if chapter["id"] == "T-CLI-08":
        _assert_gallery_bgc_definitions(root, scenario_id="T-CLI-08")
    expected_group_counts = Counter(
        {
            "og_1": 4,
            "og_2": 4,
            "og_3": 4,
            "og_4": 4,
            "og_5": 4,
            "og_6": 5,
            "og_7": 4,
            "og_8": 3,
            "og_9": 4,
            "og_10": 4,
            "og_11": 4,
            "og_12": 4,
            "og_13": 4,
            "og_15": 4,
            "og_16": 4,
            "og_17": 4,
            "og_18": 5,
            "og_19": 3,
            "og_21": 3,
            "og_22": 1,
            "og_23": 1,
        }
    )
    group_counts = Counter(
        element.attrib.get("data-orthogroup-id") for element in matches
    )
    pair_counts = Counter(
        (
            element.attrib.get("data-query-record-id"),
            element.attrib.get("data-subject-record-id"),
        )
        for element in matches
    )
    if len(matches) != 77 or group_counts != expected_group_counts:
        raise RecipeContractError("H-CLI-07 similarity-group IDs or links changed.")
    if pair_counts != Counter(
        {
            ("BGC0000708", "BGC0000709"): 19,
            ("BGC0000709", "BGC0000711"): 21,
            ("BGC0000711", "BGC0000712"): 19,
            ("BGC0000712", "BGC0000713"): 18,
        }
    ) or any(element.attrib.get("data-match-kind") != "orthogroup" for element in matches):
        raise RecipeContractError("H-CLI-07 rendered the wrong group endpoints.")

    expected_labels = {
        "livZ",
        "livY",
        "livW",
        "livV",
        "livO",
        "livX",
        "livA",
        "livG",
        "livH",
        "livI",
        "livK",
        "livB",
        "livL",
        "livD",
        "livF",
        "livP",
        "livN",
        "livQ",
        "livU",
        "livT",
        "livM",
        "livC",
        "livS",
        "livE",
    }
    text_nodes = {
        "".join(element.itertext()).strip()
        for element in root.iter()
        if element.tag.rsplit("}", 1)[-1] == "text"
    }
    if not expected_labels <= text_nodes:
        raise RecipeContractError("H-CLI-07 first-record gene labels changed.")

    og1_matches = [
        element
        for element in matches
        if element.attrib.get("data-orthogroup-id") == "og_1"
    ]
    og1_pairs = {
        (
            element.attrib["data-query-record-id"],
            element.attrib["data-subject-record-id"],
        )
        for element in og1_matches
    }
    if og1_pairs != set(_BGC_ADJACENT_PAIRS):
        raise RecipeContractError("H-CLI-07 alignment group does not span five records.")

    record_x = {}
    for element in root.iter():
        if not element.attrib.get("id", "").startswith("record_group_"):
            continue
        translation = parse_translate_chain(element.attrib.get("transform", ""))
        if translation is not None:
            record_x[element.attrib["data-gbdraw-record-id"]] = translation[0]
    feature_to_record: dict[str, str] = {}
    for element in og1_matches:
        feature_to_record[element.attrib["data-query-feature-svg-id"]] = (
            element.attrib["data-query-record-id"]
        )
        feature_to_record[element.attrib["data-subject-feature-svg-id"]] = (
            element.attrib["data-subject-record-id"]
        )
    feature_elements = {
        element.attrib["id"]: element
        for element in root.iter()
        if element.attrib.get("id") in feature_to_record
    }
    aligned_centers = []
    for feature_id, record_id in feature_to_record.items():
        coordinates = _PATH_COORDINATE_RE.findall(
            feature_elements[feature_id].attrib.get("d", "")
        )
        x_values = [float(x) for x, _y in coordinates]
        if not x_values:
            raise RecipeContractError("H-CLI-07 alignment feature geometry is missing.")
        aligned_centers.append(record_x[record_id] + (min(x_values) + max(x_values)) / 2)
    if max(aligned_centers) - min(aligned_centers) > 1e-6:
        raise RecipeContractError("H-CLI-07 CAG38695.1 group is not aligned.")


def _assert_gallery_bgc_definitions(
    root: ElementTree.Element, *, scenario_id: str
) -> None:
    definition_groups = [
        element
        for element in root.iter()
        if element.attrib.get("data-gbdraw-role") == "record-definition"
    ]
    translations = [
        parse_translate_chain(group.attrib.get("transform", ""))
        for group in definition_groups
    ]
    if len(definition_groups) != 5 or any(item is None for item in translations):
        raise RecipeContractError(f"{scenario_id} definition groups are incomplete.")
    x_positions = [translation[0] for translation in translations if translation]
    if max(x_positions) - min(x_positions) > 1e-6:
        raise RecipeContractError(
            f"{scenario_id} definitions are not locked to one left column."
        )

    expected = {
        "name": (20.0, "bold", "black"),
        "subtitle": (20.0, "normal", "black"),
        "accession": (20.0, "normal", "#7b7c7d"),
        "length": (20.0, "normal", "#7b7c7d"),
    }
    for group in definition_groups:
        lines = {
            element.attrib.get("data-definition-line-kind", ""): (
                float(element.attrib.get("font-size", "0")),
                element.attrib.get("font-weight", ""),
                element.attrib.get("fill", "").lower(),
            )
            for element in group.iter()
            if element.attrib.get("data-definition-line-kind")
        }
        anchors = {
            element.attrib.get("text-anchor", "")
            for element in group.iter()
            if element.attrib.get("data-definition-line-kind")
        }
        if lines != expected or anchors != {"start"}:
            raise RecipeContractError(
                f"{scenario_id} definition typography differs from the Gallery."
            )


def _assert_collinear_blocks(
    chapter: dict[str, object], output_path: Path
) -> None:
    _root, matches = _assert_bgc_svg(chapter, output_path)
    expected = {
        "block_0001": ("BGC0000708", "BGC0000709", "13", "plus", "#b9c2d8"),
        "block_0002": ("BGC0000708", "BGC0000709", "3", "minus", "#f3aeaf"),
        "block_0003": ("BGC0000709", "BGC0000711", "21", "plus", "#8b9cc1"),
        "block_0004": ("BGC0000711", "BGC0000712", "2", "plus", "#e1e4ed"),
        "block_0005": ("BGC0000711", "BGC0000712", "15", "plus", "#bdc6db"),
        "block_0006": ("BGC0000712", "BGC0000713", "13", "minus", "#ef9d9e"),
        "block_0007": ("BGC0000712", "BGC0000713", "2", "minus", "#f9cecf"),
    }
    observed = {
        element.attrib["data-collinearity-block-id"]: (
            element.attrib["data-query-record-id"],
            element.attrib["data-subject-record-id"],
            element.attrib["data-collinearity-anchor-count"],
            element.attrib["data-collinearity-orientation"],
            element.attrib["fill"],
        )
        for element in matches
    }
    if len(matches) != 7 or observed != expected:
        raise RecipeContractError("H-CLI-08 Collinear block semantics changed.")
    for element in matches:
        attributes = element.attrib
        anchor_count = int(attributes["data-collinearity-anchor-count"])
        query_units = attributes["data-query-unit-id"].split(";")
        subject_units = attributes["data-subject-unit-id"].split(";")
        if len(query_units) != anchor_count or len(subject_units) != anchor_count:
            raise RecipeContractError("H-CLI-08 block anchor membership is wrong.")
        if (
            attributes.get("data-match-kind") != "collinear"
            or attributes.get("data-collinearity-block-kind") != "cluster"
            or attributes.get("data-collinearity-anchor-index") != "1"
            or attributes.get("data-collinearity-color-mode")
            != "orientation_identity"
        ):
            raise RecipeContractError("H-CLI-08 block metadata is incomplete.")


def _semantic_slots(root: ElementTree.Element) -> list[tuple[str, str]]:
    return [
        (
            element.attrib["data-gbdraw-slot-id"],
            element.attrib["data-gbdraw-slot-renderer"],
        )
        for element in root.iter()
        if "data-gbdraw-slot-id" in element.attrib
        and not element.attrib["data-gbdraw-slot-id"].startswith("__")
    ]


def _text_nodes(root: ElementTree.Element) -> set[str]:
    return {
        "".join(element.itertext()).strip()
        for element in root.iter()
        if element.tag.rsplit("}", 1)[-1] == "text"
    }


def _assert_baseline_svg(
    chapter: dict[str, object],
    output_path: Path,
    *,
    record_id: str,
    feature_count: int,
) -> None:
    evidence = inspect_standard_svg(chapter, output_path=output_path)
    if evidence.record_ids != {record_id} or len(evidence.feature_ids) != feature_count:
        raise RecipeContractError(
            f"{chapter['id']} baseline does not contain the complete expected record."
        )
    if _match_elements(ElementTree.parse(output_path).getroot()):
        raise RecipeContractError(f"{chapter['id']} baseline unexpectedly contains matches.")


def _assert_generated_tables_are_documented(
    chapter: dict[str, object], scenario_id: str
) -> None:
    tables = _GENERATED_TABLES.get(scenario_id, {})
    if not tables:
        return
    destination = Path(__file__).resolve().parents[2] / str(chapter["destination"])
    source = destination.read_text(encoding="utf-8")
    for table_source in tables.values():
        documented = f"```tsv\n{table_source.rstrip()}\n```"
        if documented not in source:
            raise RecipeContractError(
                f"{scenario_id} generated table content differs from its guide."
            )


def _assert_table_inputs(workdir: Path) -> None:
    from gbdraw.io.cli_tables import (
        read_circular_track_table,
        read_comparisons_table,
        read_conservation_table,
        read_records_table,
    )

    records = read_records_table(str(workdir / "tables/records.tsv"))
    expected_ids = [
        "NC_012920.1",
        "NC_002333.2",
        "NC_024511.2",
        "NC_001328.1",
    ]
    if (
        records.record_ids != expected_ids
        or records.multi_record_positions() != ["#1@1", "#2@1", "#3@2", "#4@2"]
        or [row.column for row in records.rows] != [1, 2, 1, 2]
    ):
        raise RecipeContractError("H-CLI-02 records table semantics changed.")
    if [Path(path).name for path in records.dependency_paths] != [
        "HmmtDNA.gbk",
        "NC_002333.2.gb",
        "NC_024511.2.gb",
        "NC_001328.1.gb",
    ] or any(Path(path).parent != workdir for path in records.dependency_paths):
        raise RecipeContractError(
            "H-CLI-02 records table paths did not resolve from the table directory."
        )

    comparisons = read_comparisons_table(str(workdir / "tables/comparisons.tsv"))
    if len(comparisons.rows) != 1 or (
        comparisons.rows[0].query,
        comparisons.rows[0].subject,
        Path(comparisons.rows[0].blast).name,
    ) != ("NC_001416.1", "NC_042057.1", "lambda-de3.losatn.tsv"):
        raise RecipeContractError("H-CLI-02 comparisons table semantics changed.")
    conservation = read_conservation_table(str(workdir / "tables/conservation.tsv"))
    if (
        conservation.labels
        != [item[2] for item in _METAZOAN_CONSERVATION_RINGS]
        or conservation.colors
        != [item[3].upper() for item in _METAZOAN_CONSERVATION_RINGS]
        or [Path(path).name for path in conservation.dependency_paths]
        != [item[0] for item in _METAZOAN_CONSERVATION_RINGS]
    ):
        raise RecipeContractError("H-CLI-02 conservation table semantics changed.")
    if any(
        Path(path).parent != workdir
        for path in (*comparisons.dependency_paths, *conservation.dependency_paths)
    ):
        raise RecipeContractError(
            "H-CLI-02 BLAST paths did not resolve from the table directory."
        )

    tracks = read_circular_track_table(str(workdir / "tables/tracks.tsv"))
    if tracks.axis_index != 1 or [
        spec.split(":", 1)[0] for spec in tracks.slot_specs
    ] != ["ticks", "features", "gc_content", "gc_skew", "at_skew"]:
        raise RecipeContractError("H-CLI-02 track table semantics changed.")


def _assert_hcli02_output(
    output_path: Path,
    *,
    used_entries: list[dict[str, object]],
) -> None:
    root = ElementTree.parse(output_path).getroot()
    if output_path.name == "record_table.svg":
        _assert_multi_record_circular_layout(output_path, used_entries=used_entries)
        return

    if output_path.name == "comparison_table.svg":
        matches = _match_elements(root)
        lengths = {
            int(float(element.attrib["data-alignment-length"]))
            for element in matches
        }
        if len(matches) != 6 or lengths != {21_232, 6_412, 5_205, 1_914, 1_620, 254}:
            raise RecipeContractError("H-CLI-02 comparison-table HSPs changed.")
        for element in matches:
            attributes = element.attrib
            if (
                attributes.get("data-query-record-id"),
                attributes.get("data-subject-record-id"),
            ) != ("NC_001416.1", "NC_042057.1"):
                raise RecipeContractError("H-CLI-02 comparison endpoints changed.")
            if (
                attributes.get("data-match-kind") != "pairwise"
                or attributes.get("data-pairwise-match-style") != "curve"
            ):
                raise RecipeContractError("H-CLI-02 Linear link semantics changed.")
        return

    if output_path.name == "conservation_table.svg":
        _assert_complete_metazoan_mtdna_conservation(output_path)
        return

    if output_path.name == "annotation_table.svg":
        _assert_annotation_slots(output_path, command=None, workdir=None)
        return

    if output_path.name == "track_table.svg":
        if _semantic_slots(root) != [
            ("ticks", "ticks"),
            ("features", "features"),
            ("gc_content", "dinucleotide_content"),
            ("gc_skew", "dinucleotide_skew"),
            ("at_skew", "dinucleotide_skew"),
        ]:
            raise RecipeContractError("H-CLI-02 track-table render order changed.")
        required = {"GC content", "GC skew (+)", "GC skew (-)", "AT skew (+)", "AT skew (-)"}
        if not required <= _text_nodes(root):
            raise RecipeContractError("H-CLI-02 track-table legend changed.")
        return
    raise RecipeContractError(f"Unexpected H-CLI-02 output: {output_path.name}")


def _slot_specs_from_command(command: list[str]) -> list[str]:
    return [
        command[index + 1]
        for index, token in enumerate(command[:-1])
        if token == "--circular_track_slot"
    ]


def _assert_quantitative_tracks(
    output_path: Path,
    *,
    command: list[str],
    workdir: Path,
    scenario_id: str = "H-CLI-09",
) -> None:
    from gbdraw.tracks.circular import parse_circular_track_slots

    root = ElementTree.parse(output_path).getroot()
    expected_slots = [
        ("ticks", "ticks"),
        ("features", "features"),
        ("depth_1", "depth"),
        ("gc_content", "dinucleotide_content"),
        ("gc_skew", "dinucleotide_skew"),
    ]
    include_at_skew = any(spec.startswith("at_skew:") for spec in _slot_specs_from_command(command))
    if include_at_skew:
        expected_slots.append(("at_skew", "dinucleotide_skew"))
    if _semantic_slots(root) != expected_slots:
        raise RecipeContractError(f"{scenario_id} quantitative track order changed.")
    slots = parse_circular_track_slots(_slot_specs_from_command(command))
    if (
        [slot.id for slot in slots] != [item[0] for item in expected_slots]
        or [slot.side for slot in slots]
        != ["outside", "overlay", "inside", "inside", "inside"]
        + (["inside"] if include_at_skew else [])
        or slots[2].params.get("track_index") != "0"
        or slots[4].params.get("nt") != "GC"
        or (include_at_skew and slots[5].params.get("nt") != "AT")
    ):
        raise RecipeContractError(f"{scenario_id} slot ownership changed.")

    rows = [
        line.split("\t")
        for line in (workdir / "AP027133.DRR394922.depth-1kb.tsv")
        .read_text(encoding="utf-8")
        .splitlines()
        if line.strip() and not line.startswith("reference_name\t")
    ]
    positions = [int(row[1]) for row in rows]
    values = [float(row[2]) for row in rows]
    if (
        len(rows) != 607
        or {row[0] for row in rows} != {"AP027133.1"}
        or positions != list(range(1, 606_002, 1000))
        or min(values) != 12.446
        or max(values) != 74.546
    ):
        raise RecipeContractError(f"{scenario_id} depth fixture semantics changed.")

    paths_by_slot: dict[str, list[ElementTree.Element]] = {}
    for slot in root.iter():
        slot_id = slot.attrib.get("data-gbdraw-slot-id")
        if slot_id in {"depth_1", "gc_content", "gc_skew", "at_skew"}:
            paths_by_slot[slot_id] = [
                child
                for child in slot.iter()
                if child.tag.rsplit("}", 1)[-1] == "path"
            ]
    depth_paths = paths_by_slot.get("depth_1", [])
    if (
        len(depth_paths) != 1
        or depth_paths[0].attrib.get("fill") != "#2563EB"
        or depth_paths[0].attrib.get("d", "").count("L") != 1214
    ):
        raise RecipeContractError(f"{scenario_id} depth series geometry changed.")
    if [path.attrib.get("d", "").count("L") for path in paths_by_slot.get("gc_skew", [])[1:]] != [607, 607]:
        raise RecipeContractError(f"{scenario_id} GC-skew series changed.")
    if include_at_skew and [path.attrib.get("d", "").count("L") for path in paths_by_slot.get("at_skew", [])[1:]] != [607, 607]:
        raise RecipeContractError(f"{scenario_id} AT-skew series changed.")
    required_text = {
        "0x",
        "20x",
        "40x",
        "60x",
        "80x",
        "10%",
        "20%",
        "30%",
        "40%",
        "50%",
        "55%",
        "DRR394922 mean depth",
        "GC content (%)",
        "GC skew (+)",
        "GC skew (-)",
    }
    if include_at_skew:
        required_text.update({"AT skew (+)", "AT skew (-)"})
    if not required_text <= _text_nodes(root):
        raise RecipeContractError(f"{scenario_id} quantitative axes or legend changed.")
    required_args = {
        "--depth_window": "1",
        "--depth_step": "1000",
        "--depth_min": "0",
        "--depth_max": "80",
        "--gc_content_min_percent": "10",
        "--gc_content_max_percent": "55",
    }
    if "--no_depth_log_scale" not in command or any(
        command[command.index(option) + 1] != value
        for option, value in required_args.items()
    ):
        raise RecipeContractError(f"{scenario_id} numeric scale options changed.")


def _assert_annotation_slots(
    output_path: Path,
    *,
    command: list[str] | None,
    workdir: Path | None,
) -> None:
    root = ElementTree.parse(output_path).getroot()
    if _semantic_slots(root) != [
        ("ticks", "ticks"),
        ("features", "features"),
        ("plastome_regions", "annotations"),
        ("gc_content", "dinucleotide_content"),
        ("gc_skew", "dinucleotide_skew"),
    ]:
        raise RecipeContractError("H-CLI-10 annotation slot order changed.")
    annotations = [
        element
        for element in root.iter()
        if "data-gbdraw-annotation-id" in element.attrib
    ]
    if [element.attrib["data-gbdraw-annotation-id"] for element in annotations] != [
        "lsc",
        "irb",
        "ssc",
        "ira",
    ] or any(
        element.attrib.get("data-gbdraw-annotation-set-id") != "plastome_regions"
        or element.attrib.get("data-gbdraw-annotation-track-id") != "plastome_regions"
        or element.attrib.get("data-gbdraw-annotation-mark") != "bracket"
        for element in annotations
    ):
        raise RecipeContractError("H-CLI-10 annotation identities changed.")
    if len([element for element in root.iter() if element.attrib.get("id") == "Axis"]) != 1:
        raise RecipeContractError("H-CLI-10 feature-axis ownership changed.")
    if not {"LSC", "IRb", "SSC", "IRa"} <= _text_nodes(root):
        raise RecipeContractError("H-CLI-10 region labels changed.")
    if command is None or workdir is None:
        return

    from gbdraw.annotations import read_annotation_table
    from gbdraw.tracks.circular import parse_circular_track_slots

    slots = parse_circular_track_slots(_slot_specs_from_command(command))
    region_slot = slots[2]
    if (
        [slot.id for slot in slots]
        != ["ticks", "features", "plastome_regions", "gc_content", "gc_skew"]
        or [slot.side for slot in slots]
        != ["outside", "overlay", "inside", "inside", "inside"]
        or slots[1].params.get("lane_direction") != "split"
        or region_slot.renderer != "annotations"
        or region_slot.params.get("set_id") != "plastome_regions"
        or region_slot.width is None
        or (region_slot.width.value, region_slot.width.unit) != (30.0, "px")
    ):
        raise RecipeContractError("H-CLI-10 parsed slot ownership changed.")
    sets = read_annotation_table(
        str(workdir / "nicotiana-tabacum-regions.tsv"), mode="circular"
    )
    if len(sets) != 1 or sets[0].id != "plastome_regions":
        raise RecipeContractError("H-CLI-10 annotation set changed.")
    rows = sets[0].annotations
    if (
        [row.id for row in rows] != ["lsc", "irb", "ssc", "ira"]
        or [row.lane for row in rows] != [0, 0, 0, 0]
        or [(row.target.start, row.target.end) for row in rows]
        != [(1, 86_686), (86_687, 112_029), (112_030, 130_600), (130_601, 155_943)]
    ):
        raise RecipeContractError("H-CLI-10 annotation lanes or ranges changed.")


def _assert_gallery_chloroplast(
    chapter: dict[str, object],
    output_path: Path,
    *,
    command: list[str],
) -> None:
    def option_value(option: str) -> str:
        try:
            return command[command.index(option) + 1]
        except (ValueError, IndexError) as exc:
            raise RecipeContractError(
                f"T-CLI-06 is missing the documented {option} value."
            ) from exc

    expected_values = {
        "--gbk": "NC_001879.gbk",
        "-t": "chloroplast_specific_table.tsv",
        "-k": "CDS,rRNA,tRNA,tmRNA,ncRNA,misc_RNA,rep_origin",
        "--species": "<i>Nicotiana tabacum</i>",
        "--track_type": "tuckin",
        "--labels": "both",
        "--label_placement": "radial",
        "--outer_label_x_radius_offset": "0.9",
        "--outer_label_y_radius_offset": "0.9",
        "--inner_label_x_radius_offset": "0.975",
        "--inner_label_y_radius_offset": "0.975",
        "--qualifier_priority": "qualifier_priority.tsv",
        "--annotation_table": "nicotiana-tabacum-regions.tsv",
        "--block_stroke_color": "black",
        "--block_stroke_width": "1",
        "--line_stroke_width": "2",
        "--axis_stroke_width": "3",
        "--definition_font_size": "28",
        "--legend": "upper_left",
    }
    if any(option_value(option) != value for option, value in expected_values.items()):
        raise RecipeContractError("T-CLI-06 Gallery presentation options changed.")
    if not {"--separate_strands", "--gc", "--no-skew"} <= set(command):
        raise RecipeContractError("T-CLI-06 strand, GC, or skew settings changed.")
    slot_specs = [
        command[index + 1]
        for index, token in enumerate(command)
        if token == "--circular_track_slot"
    ]
    if slot_specs != [
        "features:features@side=overlay,lane_direction=split",
        "plastome_regions:annotations@set_id=plastome_regions,side=inside,r=0.65,w=20px,inner_gap_px=1,outer_gap_px=1,show_labels=true,padding_px=1,overflow=compress",
        "gc_content:dinucleotide_content@side=inside,r=0.56,w=0.08,nt=GC,legend_label=GC content",
    ]:
        raise RecipeContractError("T-CLI-06 custom Circular slot contract changed.")

    evidence = inspect_standard_svg(chapter, output_path=output_path)
    if evidence.record_ids != {"NC_001879.2"} or len(evidence.feature_ids) != 147:
        raise RecipeContractError(
            "T-CLI-06 must render the complete Gallery feature selection."
        )
    if evidence.slot_renderers != {
        "features",
        "annotations",
        "dinucleotide_content",
    }:
        raise RecipeContractError("T-CLI-06 SVG contains an unexpected track.")
    required_text = {
        "Nicotiana tabacum",
        "NC_001879.2",
        "155,943 bp",
        "LSC",
        "IRb",
        "SSC",
        "IRa",
        "GC content",
        "matK",
        "psaA",
        "photosystem I",
        "photosystem II",
        "RNA polymerase",
        "rep_origin",
    }
    if not required_text <= evidence.text_nodes:
        raise RecipeContractError("T-CLI-06 SVG is missing Gallery labels.")
    if {"GC skew (+)", "GC skew (-)", "AT skew (+)", "AT skew (-)"} & (
        evidence.text_nodes
    ):
        raise RecipeContractError("T-CLI-06 must not draw a skew track.")

    root = ElementTree.parse(output_path).getroot()
    slots = [
        (
            element.attrib["data-gbdraw-slot-id"],
            element.attrib["data-gbdraw-slot-renderer"],
        )
        for element in root.iter()
        if "data-gbdraw-slot-id" in element.attrib
    ]
    if slots != [
        ("features", "features"),
        ("plastome_regions", "annotations"),
        ("gc_content", "dinucleotide_content"),
    ]:
        raise RecipeContractError("T-CLI-06 SVG slot order changed.")
    annotations = {
        (
            element.attrib["data-gbdraw-annotation-id"],
            element.attrib.get("data-gbdraw-annotation-set-id"),
            element.attrib.get("data-gbdraw-annotation-track-id"),
            element.attrib.get("data-gbdraw-record-id"),
            element.attrib.get("data-gbdraw-annotation-label"),
        )
        for element in root.iter()
        if "data-gbdraw-annotation-id" in element.attrib
    }
    if annotations != {
        (key, "plastome_regions", "plastome_regions", "NC_001879.2", label)
        for key, label in (
            ("lsc", "LSC"),
            ("irb", "IRb"),
            ("ssc", "SSC"),
            ("ira", "IRa"),
        )
    }:
        raise RecipeContractError("T-CLI-06 annotation identities changed.")
    fills = {element.attrib.get("fill") for element in root.iter()}
    if not {"#00662c", "#328925", "#bd1220", "#e95d0f", "#ffec00"} <= fills:
        raise RecipeContractError("T-CLI-06 lost functional chloroplast colors.")

    python_output = (
        Path(__file__).resolve().parents[1]
        / "images"
        / "t-py-02"
        / "python_annotated_chloroplast.svg"
    )
    if not python_output.is_file():
        raise RecipeContractError("T-CLI-06 parity target T-PY-02 is missing.")
    cli_xml = ElementTree.tostring(root, encoding="utf-8")
    python_xml = ElementTree.tostring(
        ElementTree.parse(python_output).getroot(),
        encoding="utf-8",
    )
    if cli_xml != python_xml:
        raise RecipeContractError(
            "T-CLI-06 and T-PY-02 no longer render the same SVG tree."
        )


def _assert_feature_presentation(
    output_path: Path,
    *,
    command: list[str],
    workdir: Path,
    scenario_id: str = "H-CLI-11",
) -> None:
    import math

    from Bio import SeqIO

    from gbdraw.features.ids import compute_feature_hash

    root = ElementTree.parse(output_path).getroot()
    record = SeqIO.read(workdir / "HmmtDNA.gbk", "genbank")
    rendered = {
        element.attrib["data-gbdraw-feature-id"]: element
        for element in root.iter()
        if "data-gbdraw-feature-id" in element.attrib
        and element.attrib.get("data-gbdraw-auto-feature-underlay") != "true"
    }
    underlays = [
        element
        for element in root.iter()
        if element.attrib.get("data-gbdraw-auto-feature-underlay") == "true"
    ]
    source_by_gene = {
        feature.qualifiers.get("gene", [""])[0]: feature
        for feature in record.features
        if feature.qualifiers.get("gene")
    }
    nd1_id = compute_feature_hash(source_by_gene["ND1"], record_id=record.id)
    cox1_id = compute_feature_hash(source_by_gene["COX1"], record_id=record.id)
    atp6_id = compute_feature_hash(source_by_gene["ATP6"], record_id=record.id)
    rnr1_id = compute_feature_hash(source_by_gene["RNR1"], record_id=record.id)
    first_trna = next(feature for feature in record.features if feature.type == "tRNA")
    trna_id = compute_feature_hash(first_trna, record_id=record.id)
    if scenario_id == "T-CLI-03":
        annotation_elements = [
            element
            for element in root.iter()
            if element.attrib.get("data-gbdraw-annotation-id") == "d_loop"
        ]
        if (
            len(rendered) != 37
            or cox1_id not in rendered
            or nd1_id not in rendered
            or atp6_id not in rendered
            or underlays
            or not annotation_elements
            or any(
                element.attrib.get("data-gbdraw-annotation-set-id")
                != "mitochondrial_regions"
                or element.attrib.get("data-gbdraw-annotation-mark") != "bracket"
                or element.attrib.get("data-gbdraw-annotation-label") != "D-loop"
                for element in annotation_elements
            )
        ):
            raise RecipeContractError(
                "T-CLI-03 must retain every CDS and render D-loop as a region bracket."
            )
    elif (
        len({*rendered, *(element.attrib["data-gbdraw-feature-id"] for element in underlays)})
        != 37
        or cox1_id in rendered
        or nd1_id not in rendered
        or atp6_id not in rendered
        or len(underlays) != 2
        or {element.attrib["data-gbdraw-feature-id"] for element in underlays}
        != {compute_feature_hash(next(feature for feature in record.features if feature.type == "D-loop"), record_id=record.id)}
    ):
        raise RecipeContractError(f"{scenario_id} visibility or underlay semantics changed.")
    if (
        rendered[nd1_id].attrib.get("d", "").count("L") != 4
        or rendered[rnr1_id].attrib.get("d", "").count("L") != 1
        or rendered[trna_id].attrib.get("d", "").count("L") != 2
    ):
        raise RecipeContractError(f"{scenario_id} arrow or rectangle shapes changed.")
    fills = {element.attrib.get("fill") for element in rendered.values()}
    if not {"#3B82F6", "#EF4444", "#F59E0B", "#8B5CF6", "#10B981"} <= fills:
        raise RecipeContractError(f"{scenario_id} specific feature colors changed.")
    required_text = {
        "Complex I (ND1)",
        "Oxidase II",
        "12S rRNA",
        "16S rRNA",
        "ATP8",
        "ATP6",
        "COX3",
        "CYTB",
        "NADH dehydrogenase",
        "Cytochrome c oxidase",
        "ATP synthase",
        "Cytochrome b",
        "Ribosomal RNA",
    }
    text_nodes = _text_nodes(root)
    if not required_text <= text_nodes or {"s-rRNA", "l-rRNA", "COX1"} & text_nodes:
        raise RecipeContractError(f"{scenario_id} label filtering or overrides changed.")
    foreground = list(rendered.values())
    if any(
        element.attrib.get("stroke") != "#1F2937"
        or element.attrib.get("stroke-width") != "1.5"
        for element in foreground
    ):
        raise RecipeContractError(f"{scenario_id} feature block strokes changed.")
    axes = [element for element in root.iter() if element.attrib.get("id") == "Axis"]
    axis_shapes = [
        element
        for element in axes[0].iter()
        if element is not axes[0] and element.attrib.get("stroke") is not None
    ] if len(axes) == 1 else []
    if len(axes) != 1 or not axis_shapes or (
        axis_shapes[0].attrib.get("stroke"), axis_shapes[0].attrib.get("stroke-width")
    ) != ("#374151", "4.0"):
        raise RecipeContractError(f"{scenario_id} axis stroke changed.")
    radii = set()
    for element in foreground:
        match = re.search(
            r"M\s*([-+0-9.eE]+)[, ]+([-+0-9.eE]+)", element.attrib.get("d", "")
        )
        if match is not None:
            radii.add(
                round(math.hypot(float(match.group(1)), float(match.group(2))))
            )
    if "--resolve_overlaps" not in command or len(radii) < 4:
        raise RecipeContractError(f"{scenario_id} overlap lanes changed.")


def _assert_session_overwrite_refused(
    command: list[str], *, workdir: Path, environment: dict[str, str]
) -> bytes:
    replay = [token for token in command if token != "--overwrite"]
    svg_path = workdir / "cli_session_roundtrip.svg"
    before = svg_path.read_bytes()
    result = subprocess.run(
        replay,
        cwd=workdir,
        env=environment,
        capture_output=True,
        text=True,
        timeout=120,
        check=False,
    )
    evidence = result.stdout + result.stderr
    if result.returncode == 0 or "already exist" not in evidence or "overwrite" not in evidence:
        raise RecipeContractError("H-CLI-12 did not refuse an existing output.")
    if svg_path.read_bytes() != before or (workdir / "cli_session.json.gz").exists():
        raise RecipeContractError("H-CLI-12 changed files after overwrite refusal.")
    return before


def _assert_session_roundtrip(
    workdir: Path,
    *,
    first_svg: bytes,
    replay_command: list[str],
    environment: dict[str, str],
) -> None:
    plain = json.loads((workdir / "cli_session.json").read_text(encoding="utf-8"))
    with gzip.open(workdir / "cli_session.json.gz", "rt", encoding="utf-8") as handle:
        compressed = json.load(handle)
    for payload in (plain, compressed):
        if (
            payload.get("format") != "gbdraw-session"
            or payload.get("version") != 40
            or payload.get("renderRequest", {}).get("schema") != 5
            or payload.get("renderRequest", {}).get("mode") != "circular"
        ):
            raise RecipeContractError("H-CLI-12 session schema changed.")
    ignored_session_metadata = {"createdAt", "runMetadata"}
    normalized_plain = {
        key: value for key, value in plain.items() if key not in ignored_session_metadata
    }
    normalized_compressed = {
        key: value
        for key, value in compressed.items()
        if key not in ignored_session_metadata
    }
    if normalized_plain != normalized_compressed:
        changed = sorted(
            key
            for key in normalized_plain.keys() | normalized_compressed.keys()
            if normalized_plain.get(key) != normalized_compressed.get(key)
        )
        raise RecipeContractError(
            f"H-CLI-12 plain and gzip session semantics differ: {changed}."
        )
    run_metadata = compressed.get("runMetadata", {})
    geometry = run_metadata.get("trackSlotGeometry", {})
    records = geometry.get("records", [])
    if (
        geometry.get("mode") != "circular"
        or len(records) != 1
        or records[0].get("recordId") != "NC_012920.1"
        or records[0].get("resultName") != "cli_session_roundtrip"
    ):
        raise RecipeContractError("H-CLI-12 replay geometry metadata changed.")
    if (workdir / "cli_session_roundtrip.svg").read_bytes() != first_svg:
        raise RecipeContractError("H-CLI-12 replayed SVG differs from the first render.")

    incompatible_path = workdir / "incompatible_session.json"
    incompatible = dict(plain)
    incompatible["version"] = 38
    incompatible_path.write_text(json.dumps(incompatible), encoding="utf-8")
    failure_command = list(replay_command)
    failure_command[failure_command.index("--session") + 1] = incompatible_path.name
    failure_command[failure_command.index("--session_output") + 1] = "must_not_exist.json.gz"
    failure_command[failure_command.index("-o") + 1] = "must_not_exist"
    result = subprocess.run(
        failure_command,
        cwd=workdir,
        env=environment,
        capture_output=True,
        text=True,
        timeout=120,
        check=False,
    )
    evidence = result.stdout + result.stderr
    incompatible_path.unlink()
    if result.returncode == 0 or "Unsupported session version: 38" not in evidence:
        raise RecipeContractError("H-CLI-12 accepted an incompatible session version.")
    if any(workdir.glob("must_not_exist*")):
        raise RecipeContractError("H-CLI-12 wrote output for an incompatible session.")


def _assert_export_set(workdir: Path) -> None:
    static_path = workdir / "cli_export.svg"
    interactive_path = workdir / "cli_export.interactive.svg"
    interactive_root = ElementTree.parse(interactive_path).getroot()
    if interactive_root.attrib.get("data-gbdraw-interactive-svg") != "true":
        raise RecipeContractError("H-CLI-13 interactive SVG marker is missing.")
    scripts = [
        element
        for element in interactive_root.iter()
        if element.tag.rsplit("}", 1)[-1] == "script"
    ]
    metadata = [
        element
        for element in interactive_root.iter()
        if element.attrib.get("id") == "gbdraw-interactive-feature-metadata"
    ]
    interactive_source = interactive_path.read_text(encoding="utf-8")
    if (
        len(scripts) != 1
        or len(metadata) != 1
        or metadata[0].attrib.get("data-schema") != "3"
        or any(
            token not in interactive_source
            for token in (
                "gbdraw-feature-search-controls",
                "zoomViewBy",
                "resetView",
                "Original",
                "Search",
            )
        )
    ):
        raise RecipeContractError("H-CLI-13 interactive controls or metadata changed.")
    static_root = ElementTree.parse(static_path).getroot()
    if any(
        element.tag.rsplit("}", 1)[-1] == "script" for element in static_root.iter()
    ):
        raise RecipeContractError("H-CLI-13 static SVG contains a script.")

    png = (workdir / "cli_export.png").read_bytes()
    if png[:8] != b"\x89PNG\r\n\x1a\n" or len(png) < 100_000:
        raise RecipeContractError("H-CLI-13 PNG signature or size changed.")
    width, height = struct.unpack(">II", png[16:24])
    try:
        svg_width = float(static_root.attrib["width"].removesuffix("px"))
        svg_height = float(static_root.attrib["height"].removesuffix("px"))
        view_box = tuple(float(value) for value in static_root.attrib["viewBox"].split())
    except (KeyError, ValueError) as exc:
        raise RecipeContractError("H-CLI-13 SVG canvas geometry is invalid.") from exc
    if (
        len(view_box) != 4
        or not math.isclose(view_box[2], svg_width)
        or not math.isclose(view_box[3], svg_height)
    ):
        raise RecipeContractError("H-CLI-13 SVG size and viewBox differ.")
    expected_png_dimensions = (int(svg_width), int(svg_height))
    if (
        (width, height) != expected_png_dimensions
        or min(expected_png_dimensions) <= 0
        or not png.endswith(b"IEND\xaeB`\x82")
    ):
        raise RecipeContractError(
            "H-CLI-13 PNG dimensions do not match the SVG master or its trailer changed."
        )

    pdf = (workdir / "cli_export.pdf").read_bytes()
    match = list(re.finditer(rb"startxref\s+(\d+)", pdf))
    if (
        not pdf.startswith(b"%PDF-")
        or not pdf.rstrip().endswith(b"%%EOF")
        or not match
    ):
        raise RecipeContractError("H-CLI-13 PDF signature or trailer changed.")
    xref_offset = int(match[-1].group(1))
    if not 0 < xref_offset < len(pdf) or b"/Type /XRef" not in pdf[xref_offset : xref_offset + 256]:
        raise RecipeContractError("H-CLI-13 PDF cross-reference stream is invalid.")

    for filename in ("cli_export.eps", "cli_export.ps"):
        postscript = (workdir / filename).read_bytes()
        if (
            not postscript.startswith(b"%!PS-Adobe-3.0")
            or b"%%Trailer" not in postscript
            or not postscript.rstrip().endswith(b"%%EOF")
            or len(postscript) < 50_000
        ):
            raise RecipeContractError(f"H-CLI-13 {filename} is not valid DSC PostScript.")


def _assert_tutorial_interactive_handoff(workdir: Path) -> None:
    initial_path = workdir / "interactive_human_mitochondrion.svg"
    restored_path = workdir / "restored_interactive_figure.svg"
    if initial_path.read_bytes() != restored_path.read_bytes():
        raise RecipeContractError("T-CLI-11 session replay changed the static SVG.")

    interactive_path = workdir / "interactive_human_mitochondrion.interactive.svg"
    interactive_root = ElementTree.parse(interactive_path).getroot()
    interactive_source = interactive_path.read_text(encoding="utf-8")
    feature_ids = {
        element.attrib["data-gbdraw-feature-id"]
        for element in interactive_root.iter()
        if "data-gbdraw-feature-id" in element.attrib
    }
    metadata = [
        element
        for element in interactive_root.iter()
        if element.attrib.get("id") == "gbdraw-interactive-feature-metadata"
    ]
    if (
        interactive_root.attrib.get("data-gbdraw-interactive-svg") != "true"
        or len(feature_ids) != 37
        or len(metadata) != 1
        or metadata[0].attrib.get("data-schema") != "3"
        or "COX1" not in interactive_source
        or any(
            token not in interactive_source
            for token in (
                "gbdraw-feature-search-controls",
                "zoomViewBy",
                "resetView",
                "Search",
            )
        )
    ):
        raise RecipeContractError(
            "T-CLI-11 interactive controls or feature metadata changed."
        )

    with gzip.open(
        workdir / "interactive_handoff.gbdraw-session.json.gz",
        "rt",
        encoding="utf-8",
    ) as handle:
        session = json.load(handle)
    resources = session.get("resources", {})
    if (
        session.get("format") != "gbdraw-session"
        or session.get("version") != 40
        or session.get("renderRequest", {}).get("schema") != 5
        or session.get("renderRequest", {}).get("mode") != "circular"
        or len(resources) < 2
        or not any(
            str(resource.get("name", "")).endswith("HmmtDNA.gbk")
            for resource in resources.values()
            if isinstance(resource, dict)
        )
    ):
        raise RecipeContractError("T-CLI-11 session payload changed.")

def run_scenario(
    scenario_id: str,
    *,
    output_root: Path = PUBLISHED_IMAGE_ROOT,
    check: bool = False,
) -> tuple[Path, ...]:
    chapter = load_chapter(
        scenario_id,
        expected_kind="cli-recipe",
        runner_path=RUNNER_PATH,
    )
    recipe = extract_executable_block(chapter, language="bash")
    _assert_generated_tables_are_documented(chapter, scenario_id)
    commands = _parse_commands(recipe, scenario_id=scenario_id)
    expected_mode = chapter["settings"].get("mode")
    if expected_mode is not None and any(
        command[:2] != ["gbdraw", expected_mode] for command in commands
    ):
        raise RecipeContractError(
            f"{scenario_id} must call gbdraw {chapter['settings']['mode']}."
        )
    expected_outputs = chapter["execution"]["expected_outputs"]
    if scenario_id == "H-CLI-06":
        if len(commands) != 1 or expected_outputs != [
            "cli_losatp_pairwise.tsv",
            "cli_losatp_pairwise.svg",
        ]:
            raise RecipeContractError(
                "H-CLI-06 must produce its raw TSV and SVG in one search run."
            )
        command_jobs = [(commands[0], tuple(expected_outputs))]
    elif scenario_id == "H-CLI-12":
        if len(commands) != 2 or expected_outputs != [
            "cli_session.json",
            "cli_session.json.gz",
            "cli_session_roundtrip.svg",
        ]:
            raise RecipeContractError(
                "H-CLI-12 must save plain JSON, replay it to gzip, and replace one SVG."
            )
        command_jobs = [
            (commands[0], ("cli_session.json", "cli_session_roundtrip.svg")),
            (commands[1], ("cli_session.json.gz", "cli_session_roundtrip.svg")),
        ]
    elif scenario_id == "T-CLI-11":
        if len(commands) != 2 or expected_outputs != [
            "interactive_human_mitochondrion.svg",
            "interactive_human_mitochondrion.interactive.svg",
            "interactive_handoff.gbdraw-session.json.gz",
            "restored_interactive_figure.svg",
        ]:
            raise RecipeContractError(
                "T-CLI-11 must export SVG, Interactive SVG, session, and replayed SVG."
            )
        command_jobs = [
            (commands[0], tuple(expected_outputs[:3])),
            (commands[1], (expected_outputs[3],)),
        ]
    elif scenario_id == "H-CLI-13":
        if len(commands) != 1 or expected_outputs != [
            "cli_export.svg",
            "cli_export.interactive.svg",
            "cli_export.png",
            "cli_export.pdf",
            "cli_export.eps",
            "cli_export.ps",
        ]:
            raise RecipeContractError(
                "H-CLI-13 must create the six documented export files in one render."
            )
        command_jobs = [(commands[0], tuple(expected_outputs))]
    elif len(commands) != len(expected_outputs):
        raise RecipeContractError(
            f"{scenario_id} must have one command per declared output."
        )
    else:
        command_jobs = [
            (command, (output_name,))
            for command, output_name in zip(commands, expected_outputs, strict=True)
        ]
    executable = shutil.which("gbdraw")
    if executable is None:
        raise RecipeContractError("The gbdraw executable is not available on PATH.")
    for command in commands:
        command[0] = executable

    with TemporaryDirectory(prefix=f"gbdraw-{scenario_id.lower()}-") as temp_name:
        workdir = Path(temp_name)
        copied_inputs, used_entries = copy_declared_inputs(
            chapter,
            recipe_source=_copy_recipe_source(scenario_id, recipe),
            workdir=workdir,
        )
        copied_inputs.update(_materialize_generated_tables(scenario_id, workdir))
        environment = os.environ.copy()
        environment.update({"LC_ALL": "C.UTF-8", "PYTHONHASHSEED": "0", "TZ": "UTC"})
        session_first_svg: bytes | None = None
        for job_index, (command, output_names) in enumerate(command_jobs):
            if scenario_id == "H-CLI-12" and job_index == 1:
                session_first_svg = _assert_session_overwrite_refused(
                    command,
                    workdir=workdir,
                    environment=environment,
                )
            result = subprocess.run(
                command,
                cwd=workdir,
                env=environment,
                capture_output=True,
                text=True,
                timeout=300 if scenario_id in {"T-CLI-08", "T-CLI-10", "H-CLI-06", "H-CLI-07", "H-CLI-08"} else 120,
                check=False,
            )
            if result.returncode != 0:
                raise RecipeContractError(
                    f"{scenario_id} exited {result.returncode}.\n"
                    f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}"
                )
            for output_name in output_names:
                generated_path = workdir / output_name
                if not generated_path.is_file():
                    raise RecipeContractError(
                        f"{scenario_id} did not create {output_name}."
                    )
                if generated_path.suffix != ".svg" or generated_path.name.endswith(
                    ".interactive.svg"
                ):
                    continue
                if scenario_id == "H-CLI-04":
                    _assert_linear_regions_orientation_layout(chapter, generated_path)
                elif scenario_id in {"T-CLI-06", "T-CLI-10"}:
                    inspect_standard_svg(chapter, output_path=generated_path)
                else:
                    validate_standard_svg(
                        chapter,
                        output_path=generated_path,
                        used_entries=_output_entries(
                            scenario_id,
                            output_name,
                            command=command,
                            used_entries=used_entries,
                        ),
                    )
                if scenario_id == "H-CLI-03":
                    _assert_multi_record_circular_layout(
                        generated_path,
                        used_entries=_output_entries(
                            scenario_id,
                            output_name,
                            command=command,
                            used_entries=used_entries,
                        ),
                    )
                elif scenario_id == "T-CLI-03":
                    if output_name == "mitochondrial_features_baseline.svg":
                        _assert_baseline_svg(
                            chapter,
                            generated_path,
                            record_id="NC_012920.1",
                            feature_count=37,
                        )
                    else:
                        _assert_feature_presentation(
                            generated_path,
                            command=command,
                            workdir=workdir,
                            scenario_id=scenario_id,
                        )
                elif scenario_id == "T-CLI-05":
                    if output_name == "quantitative_genome_baseline.svg":
                        _assert_baseline_svg(
                            chapter,
                            generated_path,
                            record_id="AP027133.1",
                            feature_count=576,
                        )
                    else:
                        _assert_quantitative_tracks(
                            generated_path,
                            command=command,
                            workdir=workdir,
                            scenario_id=scenario_id,
                        )
                elif scenario_id == "T-CLI-06":
                    _assert_gallery_chloroplast(
                        chapter,
                        generated_path,
                        command=command,
                    )
                elif scenario_id == "T-CLI-07":
                    _assert_tutorial_losatn_comparison(chapter, generated_path)
                elif scenario_id == "T-CLI-09":
                    if "--track_type" not in command or command[
                        command.index("--track_type") + 1
                    ] != "middle":
                        raise RecipeContractError(
                            "T-CLI-09 must use the Middle Circular track preset."
                        )
                    if "--plot_title_position" not in command or command[
                        command.index("--plot_title_position") + 1
                    ] != "bottom":
                        raise RecipeContractError(
                            "T-CLI-09 must place its plot title at the bottom."
                        )
                    _assert_complete_metazoan_mtdna_conservation(
                        generated_path,
                        title="Precomputed TLOSATX rings around Homo sapiens mtDNA",
                    )
                elif scenario_id == "T-CLI-10":
                    _assert_hepatoplasmataceae_collinear(chapter, generated_path)
                elif scenario_id == "H-CLI-05":
                    _assert_precomputed_comparison(chapter, generated_path)
                elif scenario_id == "H-CLI-02":
                    _assert_hcli02_output(
                        generated_path,
                        used_entries=_output_entries(
                            scenario_id,
                            output_name,
                            command=command,
                            used_entries=used_entries,
                        ),
                    )
                elif scenario_id == "H-CLI-09":
                    _assert_quantitative_tracks(
                        generated_path,
                        command=command,
                        workdir=workdir,
                        scenario_id=scenario_id,
                    )
                elif scenario_id == "H-CLI-10":
                    _assert_annotation_slots(
                        generated_path,
                        command=command,
                        workdir=workdir,
                    )
                elif scenario_id == "H-CLI-11":
                    _assert_feature_presentation(
                        generated_path,
                        command=command,
                        workdir=workdir,
                        scenario_id=scenario_id,
                    )
            if scenario_id in {"T-CLI-08", "T-CLI-10", "H-CLI-06", "H-CLI-07", "H-CLI-08"}:
                _assert_pinned_losat(command, scenario_id=scenario_id)
            if scenario_id == "H-CLI-06":
                _assert_pairwise_protein_search(
                    chapter,
                    tsv_path=workdir / "cli_losatp_pairwise.tsv",
                    svg_path=workdir / "cli_losatp_pairwise.svg",
                )
            elif scenario_id == "H-CLI-07":
                _assert_similarity_groups(chapter, workdir / output_names[0])
            elif scenario_id == "T-CLI-08":
                _assert_similarity_groups(chapter, workdir / output_names[0])
            elif scenario_id == "H-CLI-08":
                _assert_collinear_blocks(chapter, workdir / output_names[0])
        if scenario_id == "H-CLI-02":
            _assert_table_inputs(workdir)
        elif scenario_id == "H-CLI-12":
            if session_first_svg is None:
                raise RecipeContractError("H-CLI-12 did not run its replay probe.")
            _assert_session_roundtrip(
                workdir,
                first_svg=session_first_svg,
                replay_command=commands[1],
                environment=environment,
            )
        elif scenario_id == "H-CLI-13":
            _assert_export_set(workdir)
        elif scenario_id == "T-CLI-11":
            _assert_tutorial_interactive_handoff(workdir)
        assert_exact_workdir_files(
            chapter,
            workdir=workdir,
            copied_inputs=copied_inputs,
        )
        if scenario_id == "H-CLI-01":
            _assert_gff_id_mismatch_is_rejected(
                command=commands[1],
                workdir=workdir,
                environment=environment,
            )
        return tuple(
            publish_output(
                chapter,
                generated_path=workdir / output_name,
                output_root=output_root,
                check=check,
                comparison_payload=lambda path: _artifact_comparison_payload(
                    scenario_id,
                    path,
                ),
            )
            for output_name in expected_outputs
        )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    selection = parser.add_mutually_exclusive_group(required=True)
    selection.add_argument(
        "--all", action="store_true", help="Run all implemented CLI recipes."
    )
    selection.add_argument("--scenario", choices=IMPLEMENTED_SCENARIOS)
    parser.add_argument(
        "--output-root",
        type=Path,
        default=PUBLISHED_IMAGE_ROOT,
        help="Root directory for scenario-owned SVG artifacts.",
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help="Regenerate in a clean directory and fail if the published SVG differs.",
    )
    args = parser.parse_args()

    scenario_ids = IMPLEMENTED_SCENARIOS if args.all else (args.scenario,)
    try:
        for scenario_id in scenario_ids:
            destinations = run_scenario(
                scenario_id,
                output_root=args.output_root.resolve(),
                check=args.check,
            )
            action = "verified" if args.check else "wrote"
            for destination in destinations:
                print(f"{scenario_id}: {action} {destination}")
    except RecipeContractError as exc:
        parser.exit(1, f"recipe contract failed: {exc}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
