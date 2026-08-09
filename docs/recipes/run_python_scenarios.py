#!/usr/bin/env python3
"""Execute the canonical Python onboarding recipes in clean directories."""

from __future__ import annotations

import argparse
import gzip
import io
import json
import os
from contextlib import redirect_stdout
from pathlib import Path
from tempfile import TemporaryDirectory
from xml.etree import ElementTree

from Bio.SeqRecord import SeqRecord
from pandas import DataFrame, read_csv

from gbdraw import (
    CircularOptions,
    CircularTrackOptions,
    DepthTrackOptions,
    Diagram,
    FeatureOptions,
    LabelOptions,
    LinearComparisonOptions,
)
from gbdraw.api import (
    AnnotationOptions,
    CircularDiagramRequest,
    CircularRequestPlan,
    CircularTrackSlot,
    PreparedDiagramRequest,
    RequestRenderResult,
    ScalarSpec,
    SessionDocument,
    load_session_document,
)
from gbdraw.session_io import CURRENT_SESSION_VERSION
from gbdraw.session_request_codec import CANONICAL_REQUEST_SCHEMA

if __package__:
    from ._scenario_support import (
        PUBLISHED_IMAGE_ROOT,
        RecipeContractError,
        assert_gallery_bgc_definitions,
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
        assert_gallery_bgc_definitions,
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
    "T-PY-01",
    "T-PY-02",
    "T-PY-03",
    "T-PY-04",
    "T-PY-05",
    "T-PY-06",
    "T-PY-07",
    "T-PY-08",
    "T-PY-09",
    "T-PY-11",
    "H-PY-01",
    "H-PY-02",
    "H-PY-03",
    "H-PY-04",
    "H-PY-05",
)
RUNNER_PATH = "docs/recipes/run_python_scenarios.py"
def run_scenario(
    scenario_id: str,
    *,
    output_root: Path = PUBLISHED_IMAGE_ROOT,
    check: bool = False,
) -> tuple[Path, ...]:
    chapter = load_chapter(
        scenario_id,
        expected_kind="python-recipe",
        runner_path=RUNNER_PATH,
    )
    recipe = extract_executable_block(chapter, language="python")

    with TemporaryDirectory(prefix=f"gbdraw-{scenario_id.lower()}-") as temp_name:
        workdir = Path(temp_name)
        copied_inputs, used_entries = copy_declared_inputs(
            chapter,
            recipe_source=recipe,
            workdir=workdir,
        )
        namespace: dict[str, object] = {
            "__name__": "__gbdraw_documented_recipe__",
        }
        output = io.StringIO()
        previous_cwd = Path.cwd()
        try:
            os.chdir(workdir)
            with redirect_stdout(output):
                exec(compile(recipe, chapter["destination"], "exec"), namespace)
        finally:
            os.chdir(previous_cwd)

        generated_paths = tuple(
            workdir / name for name in chapter["execution"]["expected_outputs"]
        )
        missing = [path.name for path in generated_paths if not path.is_file()]
        if missing:
            raise RecipeContractError(
                f"{scenario_id} did not create its declared output(s): {missing}."
            )
        assert_exact_workdir_files(
            chapter,
            workdir=workdir,
            copied_inputs=copied_inputs,
        )
        _validate_scenario(
            chapter,
            namespace=namespace,
            stdout=output.getvalue(),
            generated_paths=generated_paths,
            used_entries=used_entries,
        )
        return tuple(
            publish_output(
                chapter,
                generated_path=generated_path,
                output_root=output_root,
                check=check,
            )
            for generated_path in generated_paths
        )


def _validate_scenario(
    chapter: dict[str, object],
    *,
    namespace: dict[str, object],
    stdout: str,
    generated_paths: tuple[Path, ...],
    used_entries: list[dict[str, object]],
) -> None:
    scenario_id = str(chapter["id"])
    outputs = {path.name: path for path in generated_paths}
    if scenario_id == "T-PY-01":
        _validate_first_python_tutorial(
            chapter,
            namespace=namespace,
            stdout=stdout,
            output_path=generated_paths[0],
            used_entries=used_entries,
        )
    elif scenario_id == "T-PY-02":
        _validate_annotated_chloroplast_tutorial(
            chapter,
            namespace=namespace,
            stdout=stdout,
            output_path=generated_paths[0],
        )
    elif scenario_id in {"T-PY-03", "T-PY-04", "T-PY-05", "T-PY-06", "T-PY-07"}:
        _validate_migrated_tutorial(
            chapter,
            namespace=namespace,
            stdout=stdout,
            output_path=generated_paths[0],
            used_entries=used_entries,
        )
    elif scenario_id == "T-PY-08":
        _validate_interactive_handoff_tutorial(
            chapter,
            namespace=namespace,
            stdout=stdout,
            outputs=outputs,
            used_entries=used_entries,
        )
    elif scenario_id in {"T-PY-09", "T-PY-11"}:
        _validate_cli_origin_tutorial(
            chapter,
            namespace=namespace,
            stdout=stdout,
            output_path=generated_paths[0],
            used_entries=used_entries,
        )
    elif scenario_id == "H-PY-01":
        _validate_circular_and_multi_record_howto(
            chapter,
            namespace=namespace,
            stdout=stdout,
            outputs=outputs,
        )
    elif scenario_id == "H-PY-02":
        _validate_linear_comparison_howto(
            chapter,
            namespace=namespace,
            stdout=stdout,
            outputs=outputs,
        )
    elif scenario_id == "H-PY-03":
        _validate_tracks_annotations_howto(
            chapter,
            namespace=namespace,
            stdout=stdout,
            output_path=generated_paths[0],
            used_entries=used_entries,
        )
    elif scenario_id == "H-PY-04":
        _validate_gff_memory_and_bytes_howto(
            chapter,
            namespace=namespace,
            stdout=stdout,
            outputs=outputs,
        )
    elif scenario_id == "H-PY-05":
        _validate_typed_request_session_howto(
            chapter,
            namespace=namespace,
            stdout=stdout,
            outputs=outputs,
        )
    else:  # pragma: no cover - IMPLEMENTED_SCENARIOS owns this dispatch.
        raise RecipeContractError(f"No Python validator for {scenario_id}.")


def _validate_migrated_tutorial(
    chapter: dict[str, object],
    *,
    namespace: dict[str, object],
    stdout: str,
    output_path: Path,
    used_entries: list[dict[str, object]],
) -> None:
    scenario_id = str(chapter["id"])
    expected_type = RequestRenderResult if scenario_id == "T-PY-05" else Diagram
    if not isinstance(namespace.get("diagram"), expected_type):
        raise RecipeContractError(
            f"{scenario_id} must retain its public render result in `diagram`."
        )
    if namespace.get("saved_path") != Path(output_path.name):
        raise RecipeContractError(f"{scenario_id} did not retain its saved path.")
    if f"Saved {output_path.name}" not in stdout:
        raise RecipeContractError(f"{scenario_id} did not report its saved diagram.")
    if scenario_id != "T-PY-07":
        validate_standard_svg(
            chapter,
            output_path=output_path,
            used_entries=used_entries,
        )
    evidence = inspect_standard_svg(chapter, output_path=output_path)
    expected_records = {
        "T-PY-03": frozenset({"NC_001416.1"}),
        "T-PY-04": frozenset({"NC_001416.1", "NC_042057.1"}),
        "T-PY-05": frozenset({
            "BGC0000708",
            "BGC0000709",
            "BGC0000711",
            "BGC0000712",
            "BGC0000713",
        }),
        "T-PY-06": frozenset({"NC_012920.1"}),
        "T-PY-07": frozenset({
            "AP027078.1",
            "AP027131.1",
            "AP027133.1",
            "AP027132.1",
            "NZ_CP006932.1",
        }),
    }[scenario_id]
    if evidence.record_ids != expected_records:
        raise RecipeContractError(f"{scenario_id} rendered the wrong record set.")
    if scenario_id == "T-PY-04":
        _assert_losatn_matches(chapter, output_path=output_path)
    elif scenario_id == "T-PY-05":
        root = ElementTree.parse(output_path).getroot()
        assert_gallery_bgc_definitions(root, scenario_id="T-PY-05")
        matches = [
            element
            for element in root.iter()
            if element.attrib.get("data-match-kind") == "orthogroup"
        ]
        if len(matches) != 77 or not any(
            element.attrib.get("data-orthogroup-id") == "og_1"
            for element in matches
        ):
            raise RecipeContractError(
                "T-PY-05 changed the Similarity-group link set."
            )
    elif scenario_id == "T-PY-06":
        options = namespace.get("options")
        if (
            not isinstance(options, CircularOptions)
            or dict(options.config_overrides or {}).get("canvas.circular.track_type")
            != "middle"
            or options.title.position != "bottom"
        ):
            raise RecipeContractError(
                "T-PY-06 must use the Middle Circular track preset and a bottom plot title."
            )
        root = ElementTree.parse(output_path).getroot()
        matches = [
            element
            for element in root.iter()
            if element.attrib.get("data-match-kind") == "homology"
        ]
        if len(matches) != 106 or {
            element.attrib.get("data-track-label") for element in matches
        } != {
            "Danio rerio (NC_002333.2)",
            "Drosophila melanogaster (NC_024511.2)",
            "Caenorhabditis elegans (NC_001328.1)",
        }:
            raise RecipeContractError("T-PY-06 changed its three comparison rings.")
    elif scenario_id == "T-PY-07":
        if len(evidence.feature_ids) != 2_994:
            raise RecipeContractError(
                "T-PY-07 changed its complete-record feature set."
            )
        root = ElementTree.parse(output_path).getroot()
        matches = [
            element
            for element in root.iter()
            if element.attrib.get("data-match-kind") == "collinear"
        ]
        if len(matches) != 500:
            raise RecipeContractError("T-PY-07 changed its Collinear match set.")


def _validate_interactive_handoff_tutorial(
    chapter: dict[str, object],
    *,
    namespace: dict[str, object],
    stdout: str,
    outputs: dict[str, Path],
    used_entries: list[dict[str, object]],
) -> None:
    request = namespace.get("request")
    result = namespace.get("result")
    restored_request = namespace.get("restored_request")
    restored_result = namespace.get("restored_result")
    if not isinstance(request, CircularDiagramRequest):
        raise RecipeContractError("T-PY-08 did not build a Circular request.")
    if not isinstance(result, RequestRenderResult):
        raise RecipeContractError("T-PY-08 did not return a render result.")
    if not isinstance(restored_request, CircularDiagramRequest):
        raise RecipeContractError("T-PY-08 did not decode its session request.")
    if not isinstance(restored_result, RequestRenderResult):
        raise RecipeContractError("T-PY-08 replay returned the wrong result type.")
    if not isinstance(namespace.get("session_document"), SessionDocument):
        raise RecipeContractError("T-PY-08 did not retain its session document.")
    if not isinstance(namespace.get("loaded_document"), SessionDocument):
        raise RecipeContractError("T-PY-08 did not reload its saved session.")

    initial_path = outputs["interactive_human_mitochondrion.svg"]
    restored_path = outputs["restored_interactive_figure.svg"]
    validate_standard_svg(
        chapter,
        output_path=initial_path,
        used_entries=used_entries,
    )
    validate_standard_svg(
        chapter,
        output_path=restored_path,
        used_entries=used_entries,
    )
    if initial_path.read_bytes() != restored_path.read_bytes():
        raise RecipeContractError("T-PY-08 session replay changed the static SVG.")

    interactive_path = outputs["interactive_human_mitochondrion.interactive.svg"]
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
        or "gbdraw-feature-search-controls" not in interactive_source
    ):
        raise RecipeContractError("T-PY-08 interactive metadata changed.")

    with gzip.open(
        outputs["interactive_handoff.gbdraw-session.json.gz"],
        "rt",
        encoding="utf-8",
    ) as handle:
        session = json.load(handle)
    if (
        session.get("format") != "gbdraw-session"
        or session.get("version") != CURRENT_SESSION_VERSION
        or session.get("createdAt") != "2026-08-04T00:00:00+00:00"
        or session.get("renderRequest", {}).get("schema")
        != CANONICAL_REQUEST_SCHEMA
        or session.get("renderRequest", {}).get("mode") != "circular"
    ):
        raise RecipeContractError("T-PY-08 wrote a non-current session payload.")
    if (
        "Exported the interactive figure and session" not in stdout
        or "Restored restored_interactive_figure.svg" not in stdout
    ):
        raise RecipeContractError("T-PY-08 did not report its handoff outputs.")


def _validate_cli_origin_tutorial(
    chapter: dict[str, object],
    *,
    namespace: dict[str, object],
    stdout: str,
    output_path: Path,
    used_entries: list[dict[str, object]],
) -> None:
    scenario_id = str(chapter["id"])
    if not isinstance(namespace.get("diagram"), Diagram):
        raise RecipeContractError(
            f"{scenario_id} retained the wrong public render result."
        )
    if namespace.get("saved_path") != Path(output_path.name):
        raise RecipeContractError(f"{scenario_id} did not retain its saved path.")
    if f"Saved {output_path.name}" not in stdout:
        raise RecipeContractError(f"{scenario_id} did not report its saved diagram.")
    validate_standard_svg(
        chapter,
        output_path=output_path,
        used_entries=used_entries,
    )

    cli_reference = {
        "T-PY-09": PUBLISHED_IMAGE_ROOT
        / "t-cli-03"
        / "mitochondrial_features_highlighted.svg",
        "T-PY-11": PUBLISHED_IMAGE_ROOT
        / "t-cli-05"
        / "quantitative_genome_map.svg",
    }[scenario_id]
    generated_tree = ElementTree.tostring(ElementTree.parse(output_path).getroot())
    reference_tree = ElementTree.tostring(ElementTree.parse(cli_reference).getroot())
    if generated_tree != reference_tree:
        raise RecipeContractError(
            f"{scenario_id} SVG differs from its CLI project figure."
        )

def _validate_first_python_tutorial(
    chapter: dict[str, object],
    *,
    namespace: dict[str, object],
    stdout: str,
    output_path: Path,
    used_entries: list[dict[str, object]],
) -> None:
    if not isinstance(namespace.get("diagram"), Diagram):
        raise RecipeContractError(
            "T-PY-01 must leave its returned Diagram in `diagram`."
        )
    if namespace.get("saved_path") != Path(output_path.name):
        raise RecipeContractError(
            "T-PY-01 must leave its saved Path in `saved_path`."
        )
    options = namespace.get("options")
    if (
        not isinstance(options, CircularOptions)
        or Path(str(options.labels.qualifier_priority)).name
        != "cds_gene_qualifier_priority.tsv"
        or options.species != "<i>Homo sapiens</i>"
        or options.legend != "right"
        or dict(options.config_overrides or {})
        != {
            "canvas.strandedness": True,
            "canvas.circular.track_type": "middle",
            "labels.circular.scope": "outer",
            "labels.circular.placement": "horizontal",
        }
    ):
        raise RecipeContractError(
            "T-PY-01 no longer matches the shared first Circular figure."
        )
    if f"Saved {output_path.name}" not in stdout:
        raise RecipeContractError("T-PY-01 did not report its visible result.")
    validate_standard_svg(
        chapter,
        output_path=output_path,
        used_entries=used_entries,
    )
    cli_output = (
        Path(__file__).resolve().parents[1]
        / "images"
        / "t-cli-01"
        / "human_mitochondrion.svg"
    )
    if ElementTree.tostring(
        ElementTree.parse(output_path).getroot(), encoding="utf-8"
    ) != ElementTree.tostring(
        ElementTree.parse(cli_output).getroot(), encoding="utf-8"
    ):
        raise RecipeContractError(
            "T-PY-01 and T-CLI-01 no longer render the same SVG tree."
        )


def _validate_annotated_chloroplast_tutorial(
    chapter: dict[str, object],
    *,
    namespace: dict[str, object],
    stdout: str,
    output_path: Path,
) -> None:
    _assert_diagram_boundary(
        chapter,
        namespace=namespace,
        prefix="chloroplast",
        output_path=output_path,
    )
    record = namespace.get("record")
    if not isinstance(record, SeqRecord) or (
        record.id,
        len(record),
        record.annotations.get("topology"),
    ) != ("NC_001879.2", 155_943, "circular"):
        raise RecipeContractError(
            "T-PY-02 must use the complete circular tobacco plastome."
        )

    options = namespace.get("options")
    if not isinstance(options, CircularOptions):
        raise RecipeContractError("T-PY-02 must retain its typed CircularOptions.")
    if tuple(options.features.types or ()) != (
        "CDS",
        "rRNA",
        "tRNA",
        "tmRNA",
        "ncRNA",
        "misc_RNA",
        "rep_origin",
    ) or Path(str(options.features.color_table)).name != (
        "chloroplast_specific_table.tsv"
    ):
        raise RecipeContractError(
            "T-PY-02 lost the Gallery feature selection or functional colors."
        )
    if Path(str(options.labels.qualifier_priority)).name != "qualifier_priority.tsv":
        raise RecipeContractError("T-PY-02 lost its CDS gene-label priority.")
    if not isinstance(options.annotations, AnnotationOptions) or Path(
        str(options.annotations.table_file)
    ).name != "nicotiana-tabacum-regions.tsv":
        raise RecipeContractError("T-PY-02 lost its plastome region table.")
    if options.species != "<i>Nicotiana tabacum</i>" or options.legend != "upper_left":
        raise RecipeContractError(
            "T-PY-02 lost the Gallery definition or legend placement."
        )

    expected_overrides = {
        "canvas.strandedness": True,
        "canvas.circular.track_type": "tuckin",
        "labels.circular.scope": "both",
        "labels.circular.placement": "radial",
        "labels.unified_adjustment.outer_labels.x_radius_offset": 0.9,
        "labels.unified_adjustment.outer_labels.y_radius_offset": 0.9,
        "labels.unified_adjustment.inner_labels.x_radius_offset": 0.975,
        "labels.unified_adjustment.inner_labels.y_radius_offset": 0.975,
        "objects.definition.circular.font_size": 28,
        "objects.definition.circular.interval": 30,
        "objects.features.block_stroke_color": "black",
        "objects.features.block_stroke_width.long": 1,
        "objects.features.line_stroke_width.long": 2,
        "objects.axis.circular.stroke_width.long": 3,
    }
    if dict(options.config_overrides or {}) != expected_overrides:
        raise RecipeContractError("T-PY-02 Gallery presentation settings changed.")

    slots = tuple(options.tracks.slots or ())
    if not all(isinstance(slot, CircularTrackSlot) for slot in slots) or [
        (slot.id, slot.renderer, slot.side) for slot in slots
    ] != [
        ("features", "features", "overlay"),
        ("plastome_regions", "annotations", "inside"),
        ("gc_content", "dinucleotide_content", "inside"),
    ]:
        raise RecipeContractError("T-PY-02 custom Circular slot order changed.")
    feature_slot, annotation_slot, gc_slot = slots
    if feature_slot.params != {"lane_direction": "split"}:
        raise RecipeContractError("T-PY-02 must split the two feature strands.")
    if (
        annotation_slot.radius != ScalarSpec(0.65)
        or annotation_slot.width != ScalarSpec(20, "px")
        or annotation_slot.inner_gap_px != 1
        or annotation_slot.outer_gap_px != 1
        or annotation_slot.params
        != {
            "set_id": "plastome_regions",
            "show_labels": True,
            "padding_px": 1,
            "overflow": "compress",
        }
    ):
        raise RecipeContractError("T-PY-02 plastome region geometry changed.")
    if (
        gc_slot.radius != ScalarSpec(0.56)
        or gc_slot.width != ScalarSpec(0.08)
        or gc_slot.params != {"nt": "GC", "legend_label": "GC content"}
    ):
        raise RecipeContractError("T-PY-02 GC-content geometry changed.")

    _assert_svg_records(
        chapter,
        output_path=output_path,
        records=(("NC_001879.2", 155_943),),
        displayed_feature_count=147,
        required_text=(
            "Nicotiana tabacum",
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
        ),
    )
    evidence = inspect_standard_svg(chapter, output_path=output_path)
    if evidence.slot_renderers != {
        "features",
        "annotations",
        "dinucleotide_content",
    }:
        raise RecipeContractError("T-PY-02 SVG contains an unexpected track.")
    if {"GC skew (+)", "GC skew (-)", "AT skew (+)", "AT skew (-)"} & (
        evidence.text_nodes
    ):
        raise RecipeContractError("T-PY-02 must not draw a skew track.")

    root = ElementTree.parse(output_path).getroot()
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
        raise RecipeContractError(
            f"T-PY-02 plastome annotation identities changed: {sorted(annotations)!r}."
        )
    fills = {element.attrib.get("fill") for element in root.iter()}
    if not {"#00662c", "#328925", "#bd1220", "#e95d0f", "#ffec00"} <= fills:
        raise RecipeContractError("T-PY-02 SVG lost functional chloroplast colors.")
    gallery_source = (
        Path(__file__).resolve().parents[2]
        / "gbdraw"
        / "web"
        / "gallery"
        / "sources"
        / "tobacco-chloroplast.svg"
    )
    if ElementTree.tostring(root, encoding="utf-8") != ElementTree.tostring(
        ElementTree.parse(gallery_source).getroot(),
        encoding="utf-8",
    ):
        raise RecipeContractError(
            "T-PY-02 no longer reproduces the Gallery chloroplast SVG tree."
        )
    if f"Saved {output_path.name}" not in stdout:
        raise RecipeContractError("T-PY-02 did not report its saved diagram.")


def _validate_circular_and_multi_record_howto(
    chapter: dict[str, object],
    *,
    namespace: dict[str, object],
    stdout: str,
    outputs: dict[str, Path],
) -> None:
    circular = outputs["python_circular.svg"]
    multi_record = outputs["python_multi_record.svg"]
    _assert_diagram_boundary(
        chapter,
        namespace=namespace,
        prefix="circular",
        output_path=circular,
    )
    _assert_diagram_boundary(
        chapter,
        namespace=namespace,
        prefix="multi_record",
        output_path=multi_record,
    )
    circular_record = namespace.get("circular_record")
    multi_records = namespace.get("multi_records")
    if not isinstance(circular_record, SeqRecord) or (
        circular_record.id,
        len(circular_record),
    ) != ("NC_012920.1", 16_569):
        raise RecipeContractError("H-PY-01 used the wrong Circular record.")
    if (
        not isinstance(multi_records, list)
        or len(multi_records) != 4
        or not all(isinstance(record, SeqRecord) for record in multi_records)
        or [(record.id, len(record)) for record in multi_records]
        != [
            ("NC_012920.1", 16_569),
            ("NC_002333.2", 16_596),
            ("NC_024511.2", 19_524),
            ("NC_001328.1", 13_794),
        ]
        or any(record.annotations.get("topology") != "circular" for record in multi_records)
    ):
        raise RecipeContractError(
            "H-PY-01 must use four complete circular mitochondrial records."
        )
    labels = namespace.get("labels")
    if not isinstance(labels, LabelOptions) or (
        Path(str(labels.qualifier_priority)).name
        != "cds_gene_qualifier_priority.tsv"
    ):
        raise RecipeContractError("H-PY-01 lost its CDS gene qualifier priority.")
    cds_label_whitelist = labels.whitelist
    if (
        not isinstance(cds_label_whitelist, DataFrame)
        or list(cds_label_whitelist.columns)
        != ["feature_type", "qualifier", "keyword"]
        or cds_label_whitelist.to_dict(orient="records")
        != [{"feature_type": "CDS", "qualifier": "gene", "keyword": ".+"}]
    ):
        raise RecipeContractError("H-PY-01 must label CDS features only.")
    expected_label_overrides = {
        "labels.circular.scope": "outer",
        "labels.circular.placement": "horizontal",
        "labels.font_size.short": 18,
        "labels.font_size.long": 18,
    }
    for option_name in ("circular_options", "multi_record_options"):
        options = namespace.get(option_name)
        if not isinstance(options, CircularOptions) or (
            options.labels is not labels
            or dict(options.config_overrides or {}) != expected_label_overrides
        ):
            raise RecipeContractError(
                f"H-PY-01 {option_name} lost its readable outer CDS labels."
            )
    _assert_svg_records(
        chapter,
        output_path=circular,
        records=(("NC_012920.1", 16_569),),
        displayed_feature_count=37,
    )
    _assert_svg_records(
        chapter,
        output_path=multi_record,
        records=(
            ("NC_012920.1", 16_569),
            ("NC_002333.2", 16_596),
            ("NC_024511.2", 19_524),
            ("NC_001328.1", 13_794),
        ),
        displayed_feature_count=147,
        required_text=(
            "Complete metazoan mitochondrial genomes",
            "Homo sapiens",
            "Danio rerio",
            "Drosophila melanogaster",
            "Caenorhabditis elegans",
        ),
    )
    _assert_cds_gene_labels(
        chapter,
        output_path=circular,
        records=(circular_record,),
        font_size=18,
    )
    _assert_cds_gene_labels(
        chapter,
        output_path=multi_record,
        records=tuple(multi_records),
        font_size=18,
    )
    if "Saved python_circular.svg and python_multi_record.svg" not in stdout:
        raise RecipeContractError("H-PY-01 did not report both saved diagrams.")


def _validate_linear_comparison_howto(
    chapter: dict[str, object],
    *,
    namespace: dict[str, object],
    stdout: str,
    outputs: dict[str, Path],
) -> None:
    comparison = outputs["python_linear_comparison.svg"]
    _assert_diagram_boundary(
        chapter,
        namespace=namespace,
        prefix="comparison",
        output_path=comparison,
    )

    comparison_records = namespace.get("comparison_records")
    if (
        not isinstance(comparison_records, list)
        or [(record.id, len(record)) for record in comparison_records]
        != [("NC_001416.1", 48_502), ("NC_042057.1", 42_925)]
    ):
        raise RecipeContractError(
            "H-PY-02 must compare the complete Lambda and DE3 records."
        )
    comparison_options = namespace.get("comparison_options")
    if not isinstance(comparison_options, LinearComparisonOptions) or (
        comparison_options.blast_files != ("lambda-de3.losatn.tsv",)
    ):
        raise RecipeContractError("H-PY-02 did not use the fixed LOSATN evidence.")

    _assert_svg_records(
        chapter,
        output_path=comparison,
        records=(("NC_001416.1", 48_502), ("NC_042057.1", 42_925)),
        displayed_feature_count=130,
    )
    _assert_linear_rows(chapter, output_path=comparison)
    _assert_losatn_matches(chapter, output_path=comparison)
    if "Saved python_linear_comparison.svg" not in stdout:
        raise RecipeContractError("H-PY-02 did not report its saved diagram.")


def _validate_tracks_annotations_howto(
    chapter: dict[str, object],
    *,
    namespace: dict[str, object],
    stdout: str,
    output_path: Path,
    used_entries: list[dict[str, object]],
) -> None:
    _assert_diagram_boundary(
        chapter,
        namespace=namespace,
        prefix="tracks",
        output_path=output_path,
    )
    records = namespace.get("records")
    if (
        not isinstance(records, list)
        or not all(isinstance(record, SeqRecord) for record in records)
        or [(record.id, len(record)) for record in records]
        != [
            ("AP027133.1", 606_194),
            ("NC_001879.2", 155_943),
            ("NC_012920.1", 16_569),
        ]
        or any(record.annotations.get("topology") != "circular" for record in records)
    ):
        raise RecipeContractError(
            "H-PY-03 must use the three complete circular source records."
        )

    options = namespace.get("options")
    if not isinstance(options, CircularOptions):
        raise RecipeContractError("H-PY-03 must retain its typed CircularOptions.")
    if options.depth_window != 1 or options.depth_step != 1000:
        raise RecipeContractError(
            "H-PY-03 must plot the pre-aggregated means with window 1 and step 1000."
        )
    if not isinstance(options.features, FeatureOptions) or (
        Path(str(options.features.default_colors)).name != "modified_default_colors.tsv"
    ):
        raise RecipeContractError("H-PY-03 lost its documented feature colors.")
    if not isinstance(options.labels, LabelOptions) or (
        Path(str(options.labels.qualifier_priority)).name != "qualifier_priority.tsv"
    ):
        raise RecipeContractError("H-PY-03 lost its gene-first label rule.")
    label_whitelist = options.labels.whitelist
    if not isinstance(label_whitelist, DataFrame) or list(label_whitelist.columns) != [
        "feature_type",
        "qualifier",
        "keyword",
    ]:
        raise RecipeContractError("H-PY-03 lost its typed in-memory label whitelist.")
    if not isinstance(options.annotations, AnnotationOptions) or (
        Path(str(options.annotations.table_file)).name
        != "nicotiana-tabacum-regions.tsv"
    ):
        raise RecipeContractError("H-PY-03 lost its plastome annotation table.")

    depth_tracks = tuple(options.depth_tracks)
    if len(depth_tracks) != 1 or not isinstance(depth_tracks[0], DepthTrackOptions):
        raise RecipeContractError("H-PY-03 must define one logical depth series.")
    depth_sources = tuple(depth_tracks[0].source)
    if (
        len(depth_sources) != 3
        or Path(str(depth_sources[0])).name
        != "AP027133.DRR394922.depth-1kb.tsv"
        or depth_sources[1:] != (None, None)
    ):
        raise RecipeContractError(
            "H-PY-03 depth evidence must be bound only to AP027133.1."
        )

    if not isinstance(options.tracks, CircularTrackOptions):
        raise RecipeContractError("H-PY-03 must use typed Circular track options.")
    slots = tuple(options.tracks.slots or ())
    if not all(isinstance(slot, CircularTrackSlot) for slot in slots) or [
        (slot.id, slot.renderer) for slot in slots
    ] != [
        ("features", "features"),
        ("plastome_regions", "annotations"),
        ("depth", "depth"),
        ("gc_content", "dinucleotide_content"),
    ]:
        raise RecipeContractError("H-PY-03 custom Circular slot order changed.")

    validate_standard_svg(
        chapter,
        output_path=output_path,
        used_entries=used_entries,
    )
    evidence = inspect_standard_svg(chapter, output_path=output_path)
    required_renderers = {
        "features",
        "annotations",
        "depth",
        "dinucleotide_content",
    }
    if not required_renderers <= evidence.slot_renderers:
        raise RecipeContractError("H-PY-03 SVG is missing a documented track renderer.")

    expected_gene_labels = {
        "rpoB",
        "secA",
        "polC",
        "psaA",
        "atpB",
        "rbcL",
        "ndhF",
        "ND1",
        "COX1",
        "ATP6",
        "CYTB",
    }
    if not expected_gene_labels <= evidence.text_nodes:
        raise RecipeContractError("H-PY-03 SVG is missing selected CDS gene labels.")
    product_labels = {
        "DNA-directed RNA polymerase subunit beta",
        "photosystem I P700 chlorophyll a apoprotein A1",
        "cytochrome c oxidase subunit I",
    }
    if product_labels & evidence.text_nodes:
        raise RecipeContractError(
            "H-PY-03 uses product descriptions instead of selected gene labels."
        )
    required_text = {
        "Candidatus Hepatoplasma scabrum",
        "Nicotiana tabacum",
        "Homo sapiens",
        "LSC",
        "IRb",
        "SSC",
        "IRa",
    }
    if not required_text <= evidence.text_nodes:
        raise RecipeContractError("H-PY-03 SVG is missing record or region labels.")

    root = ElementTree.parse(output_path).getroot()
    annotations = {
        (
            element.attrib["data-gbdraw-annotation-id"],
            element.attrib.get("data-gbdraw-record-id"),
            element.attrib.get("data-gbdraw-annotation-label"),
        )
        for element in root.iter()
        if "data-gbdraw-annotation-id" in element.attrib
    }
    if annotations != {
        ("lsc", "NC_001879.2", "LSC"),
        ("irb", "NC_001879.2", "IRb"),
        ("ssc", "NC_001879.2", "SSC"),
        ("ira", "NC_001879.2", "IRa"),
    }:
        raise RecipeContractError("H-PY-03 plastome annotation identities changed.")
    fills = {element.attrib.get("fill") for element in root.iter()}
    if not {"#d3d3d3", "#009e73", "#e69f00", "#2563eb"} <= fills:
        raise RecipeContractError("H-PY-03 SVG lost a documented feature or depth color.")

    depth_path = output_path.parent / "AP027133.DRR394922.depth-1kb.tsv"
    depth = read_csv(depth_path, sep="\t")
    if (
        list(depth.columns) != ["reference_name", "position", "depth"]
        or len(depth) != 607
        or set(depth["reference_name"]) != {"AP027133.1"}
        or depth["position"].tolist() != list(range(1, 606_195, 1000))
        or float(depth["depth"].min()) != 12.446
        or float(depth["depth"].max()) != 74.546
    ):
        raise RecipeContractError("H-PY-03 depth derivative semantics changed.")
    if "Saved python_tracks_annotations.svg" not in stdout:
        raise RecipeContractError("H-PY-03 did not report its saved diagram.")


def _validate_gff_memory_and_bytes_howto(
    chapter: dict[str, object],
    *,
    namespace: dict[str, object],
    stdout: str,
    outputs: dict[str, Path],
) -> None:
    gff_output = outputs["python_gff3.svg"]
    memory_output = outputs["python_memory.svg"]
    _assert_diagram_boundary(
        chapter,
        namespace=namespace,
        prefix="gff",
        output_path=gff_output,
    )
    _assert_diagram_boundary(
        chapter,
        namespace=namespace,
        prefix="memory",
        output_path=memory_output,
    )
    gff_records = namespace.get("gff_records")
    memory_record = namespace.get("memory_record")
    if not isinstance(gff_records, list) or len(gff_records) != 1:
        raise RecipeContractError("H-PY-04 must load one whole GFF3 record.")
    gff_record = gff_records[0]
    if not isinstance(gff_record, SeqRecord) or (
        gff_record.id,
        len(gff_record),
        len(gff_record.features),
    ) != ("NC_001416.1", 48_502, 165):
        raise RecipeContractError("H-PY-04 changed the whole Lambda record.")
    if not isinstance(memory_record, SeqRecord) or (
        memory_record.id,
        len(memory_record),
    ) != ("NC_012920.1", 16_569):
        raise RecipeContractError("H-PY-04 did not draw its in-memory SeqRecord.")
    _assert_svg_records(
        chapter,
        output_path=gff_output,
        records=(("NC_001416.1", 48_502),),
        displayed_feature_count=73,
    )
    _assert_svg_records(
        chapter,
        output_path=memory_output,
        records=(("NC_012920.1", 16_569),),
        displayed_feature_count=37,
    )
    if "Wrote python_gff3.svg and python_memory.svg" not in stdout:
        raise RecipeContractError("H-PY-04 did not report both byte outputs.")


def _validate_typed_request_session_howto(
    chapter: dict[str, object],
    *,
    namespace: dict[str, object],
    stdout: str,
    outputs: dict[str, Path],
) -> None:
    svg_output = outputs["typed_request.svg"]
    session_output = outputs["typed_request.session.json"]
    if not isinstance(namespace.get("request"), CircularDiagramRequest):
        raise RecipeContractError("H-PY-05 did not build a Circular request.")
    plan = namespace.get("plan")
    prepared = namespace.get("prepared")
    result = namespace.get("result")
    if not isinstance(plan, CircularRequestPlan) or plan.mode != "circular":
        raise RecipeContractError("H-PY-05 did not produce a Circular plan.")
    if not isinstance(prepared, PreparedDiagramRequest):
        raise RecipeContractError("H-PY-05 did not build its resolved plan.")
    if not isinstance(result, RequestRenderResult) or result.mode != "circular":
        raise RecipeContractError("H-PY-05 did not return a Circular render result.")
    if tuple(path.name for path in result.output_paths) != (svg_output.name,):
        raise RecipeContractError("H-PY-05 render result has the wrong output path.")
    if not isinstance(namespace.get("session_document"), SessionDocument):
        raise RecipeContractError("H-PY-05 did not retain its saved session document.")
    if not isinstance(namespace.get("loaded_document"), SessionDocument):
        raise RecipeContractError("H-PY-05 did not reload its session document.")
    if not isinstance(namespace.get("replay_request"), CircularDiagramRequest):
        raise RecipeContractError("H-PY-05 session did not decode to a Circular request.")
    replay_result = namespace.get("replay_result")
    if not isinstance(replay_result, RequestRenderResult):
        raise RecipeContractError("H-PY-05 session replay returned the wrong type.")
    rendered_bytes = namespace.get("rendered_svg_bytes")
    replayed_bytes = namespace.get("replayed_svg_bytes")
    if not isinstance(rendered_bytes, bytes) or b"<svg" not in rendered_bytes[:512]:
        raise RecipeContractError("H-PY-05 render did not return SVG bytes.")
    if rendered_bytes != replayed_bytes:
        raise RecipeContractError("H-PY-05 session replay changed the SVG bytes.")
    if rendered_bytes != svg_output.read_bytes():
        raise RecipeContractError("H-PY-05 saved SVG differs from its render result.")
    wrong_mode_error = namespace.get("wrong_mode_error")
    if not isinstance(wrong_mode_error, str) or (
        "Circular request cannot contain Linear track values" not in wrong_mode_error
    ):
        raise RecipeContractError("H-PY-05 did not capture the wrong-mode error.")
    if namespace.get("resources_expired") is not True:
        raise RecipeContractError("H-PY-05 retained expired session resources.")
    if namespace.get("replay_output_removed") is not True:
        raise RecipeContractError("H-PY-05 left its replay scratch output behind.")
    payload = json.loads(session_output.read_text(encoding="utf-8"))
    if (
        payload.get("format") != "gbdraw-session"
        or payload.get("version") != CURRENT_SESSION_VERSION
        or payload.get("createdAt") != "2026-08-03T00:00:00+00:00"
        or payload.get("renderRequest", {}).get("schema")
        != CANONICAL_REQUEST_SCHEMA
        or payload.get("renderRequest", {}).get("mode") != "circular"
    ):
        raise RecipeContractError("H-PY-05 wrote a non-current session payload.")
    if load_session_document(session_output).mode != "circular":
        raise RecipeContractError("H-PY-05 published session does not reload.")
    _assert_svg_records(
        chapter,
        output_path=svg_output,
        records=(("NC_012920.1", 16_569),),
        displayed_feature_count=37,
    )
    if (
        "Rendered typed_request.svg" not in stdout
        or "Replayed the current Circular session" not in stdout
        or "Rejected wrong-mode session" not in stdout
    ):
        raise RecipeContractError("H-PY-05 did not report its session checks.")


def _assert_diagram_boundary(
    chapter: dict[str, object],
    *,
    namespace: dict[str, object],
    prefix: str,
    output_path: Path,
) -> None:
    diagram = namespace.get(f"{prefix}_diagram")
    svg = namespace.get(f"{prefix}_svg")
    payload = namespace.get(f"{prefix}_bytes")
    saved_path = namespace.get(f"{prefix}_path")
    if not isinstance(diagram, Diagram):
        raise RecipeContractError(
            f"{chapter['id']} must retain `{prefix}_diagram` as a Diagram."
        )
    if not isinstance(svg, str) or not svg.startswith("<svg"):
        raise RecipeContractError(f"{chapter['id']} returned invalid SVG text.")
    if not isinstance(payload, bytes) or not payload.startswith(b"<svg"):
        raise RecipeContractError(f"{chapter['id']} returned invalid SVG bytes.")
    if svg.encode("utf-8") != payload or output_path.read_bytes() != payload:
        raise RecipeContractError(
            f"{chapter['id']} file, text, and byte output do not match."
        )
    if saved_path != Path(output_path.name):
        raise RecipeContractError(
            f"{chapter['id']} must retain `{prefix}_path` as the saved Path."
        )


def _assert_svg_records(
    chapter: dict[str, object],
    *,
    output_path: Path,
    records: tuple[tuple[str, int], ...],
    displayed_feature_count: int,
    required_text: tuple[str, ...] = (),
) -> None:
    evidence = inspect_standard_svg(chapter, output_path=output_path)
    expected_ids = frozenset(record_id for record_id, _length in records)
    if evidence.record_ids != expected_ids:
        raise RecipeContractError(
            f"{chapter['id']} SVG record IDs changed: {sorted(evidence.record_ids)}."
        )
    for record_id, length in records:
        if (
            record_id not in evidence.text_content
            or f"{length:,} bp" not in evidence.text_content
        ):
            raise RecipeContractError(
                f"{chapter['id']} SVG is missing {record_id} metadata."
            )
    if len(evidence.feature_ids) != displayed_feature_count:
        raise RecipeContractError(
            f"{chapter['id']} expected {displayed_feature_count} displayed features, "
            f"found {len(evidence.feature_ids)}."
        )
    missing_text = [text for text in required_text if text not in evidence.text_content]
    if missing_text:
        raise RecipeContractError(
            f"{chapter['id']} SVG is missing required labels: {missing_text}."
        )


def _assert_cds_gene_labels(
    chapter: dict[str, object],
    *,
    output_path: Path,
    records: tuple[SeqRecord, ...],
    font_size: float,
) -> None:
    root = ElementTree.parse(output_path).getroot()
    record_containers = {
        element.attrib["data-gbdraw-record-id"]: element
        for element in root
        if element.attrib.get("id", "").startswith("record_")
        and "data-gbdraw-record-id" in element.attrib
    }
    if len(records) == 1:
        record_containers = {records[0].id: root}
    expected_ids = {record.id for record in records}
    if set(record_containers) != expected_ids:
        raise RecipeContractError(
            f"{chapter['id']} could not resolve per-record label containers."
        )

    for record in records:
        cds_features = [feature for feature in record.features if feature.type == "CDS"]
        expected_genes = {
            str(value)
            for feature in cds_features
            for value in feature.qualifiers.get("gene", ())
            if str(value)
        }
        product_labels = {
            str(value)
            for feature in cds_features
            for value in feature.qualifiers.get("product", ())
            if str(value)
        }
        if len(expected_genes) != len(cds_features):
            raise RecipeContractError(
                f"{chapter['id']} source CDS genes are missing or duplicated for {record.id}."
            )

        text_elements = [
            element
            for element in record_containers[record.id].iter()
            if element.tag.rsplit("}", 1)[-1] == "text"
        ]
        text_by_value: dict[str, list[ElementTree.Element]] = {}
        for element in text_elements:
            text = "".join(element.itertext()).strip()
            text_by_value.setdefault(text, []).append(element)
        text_nodes = set(text_by_value)
        missing_genes = sorted(expected_genes - text_nodes)
        if missing_genes:
            raise RecipeContractError(
                f"{chapter['id']} SVG is missing CDS gene labels for {record.id}: "
                f"{missing_genes}."
            )
        visible_products = sorted(product_labels & text_nodes)
        if visible_products:
            raise RecipeContractError(
                f"{chapter['id']} SVG uses CDS product labels for {record.id}: "
                f"{visible_products}."
            )
        wrong_font = sorted(
            gene
            for gene in expected_genes
            if not any(
                float(styled.attrib.get("font-size", "nan")) == font_size
                for element in text_by_value[gene]
                for styled in element.iter()
            )
        )
        if wrong_font:
            raise RecipeContractError(
                f"{chapter['id']} SVG uses the wrong CDS label font for "
                f"{record.id}: {wrong_font}."
            )


def _assert_losatn_matches(
    chapter: dict[str, object], *, output_path: Path
) -> None:
    expected_source = output_path.parent / "lambda-de3.losatn.tsv"
    if not expected_source.is_file():
        raise RecipeContractError("H-PY-02 lost its copied LOSATN evidence.")
    expected_matches = {
        tuple(line.split("\t")[index] for index in (0, 1, 6, 7, 8, 9))
        for line in expected_source.read_text(encoding="utf-8").splitlines()
        if line
    }
    root = ElementTree.parse(output_path).getroot()
    actual_matches = {
        (
            element.attrib["data-query-record-id"],
            element.attrib["data-subject-record-id"],
            element.attrib["data-qstart"],
            element.attrib["data-qend"],
            element.attrib["data-sstart"],
            element.attrib["data-send"],
        )
        for element in root.iter()
        if "data-gbdraw-match-id" in element.attrib
    }
    if actual_matches != expected_matches:
        raise RecipeContractError(
            "H-PY-02 SVG comparison links differ from lambda-de3.losatn.tsv."
        )


def _assert_linear_rows(chapter: dict[str, object], *, output_path: Path) -> None:
    root = ElementTree.parse(output_path).getroot()
    rows: list[tuple[int, float]] = []
    for element in root.iter():
        if (
            element.tag.rsplit("}", 1)[-1] != "g"
            or not element.attrib.get("id", "").startswith("record_group_")
        ):
            continue
        index = int(element.attrib["data-gbdraw-record-index"])
        translation = parse_translate_chain(element.attrib.get("transform", ""))
        if translation is None:
            raise RecipeContractError("H-PY-02 record rows lost their translation.")
        rows.append((index, translation[1]))
    if rows != sorted(rows) or len(rows) != 2 or rows[0][1] >= rows[1][1]:
        raise RecipeContractError("H-PY-02 SVG does not retain two ordered Linear rows.")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    selection = parser.add_mutually_exclusive_group(required=True)
    selection.add_argument(
        "--all", action="store_true", help="Run all implemented Python recipes."
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
            print(
                f"{scenario_id}: {action} "
                + ", ".join(str(destination) for destination in destinations)
            )
    except RecipeContractError as exc:
        parser.exit(1, f"recipe contract failed: {exc}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
