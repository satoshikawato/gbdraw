from __future__ import annotations

import hashlib
import json
import re
import shutil
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path
from typing import Any
from xml.etree import ElementTree


REPO_ROOT = Path(__file__).resolve().parents[2]
SCENARIO_MANIFEST_PATH = REPO_ROOT / "docs" / "scenarios" / "manifest.json"
FIXTURE_ROOT = REPO_ROOT / "gbdraw" / "web" / "tutorial-data"
FIXTURE_MANIFEST_PATH = FIXTURE_ROOT / "manifest.json"
PUBLISHED_IMAGE_ROOT = REPO_ROOT / "docs" / "images"


class RecipeContractError(RuntimeError):
    """A documented recipe no longer satisfies its manifest contract."""


def load_chapter(
    scenario_id: str,
    *,
    expected_kind: str,
    runner_path: str,
) -> dict[str, Any]:
    manifest = json.loads(SCENARIO_MANIFEST_PATH.read_text(encoding="utf-8"))
    chapter = next(
        (item for item in manifest["chapters"] if item["id"] == scenario_id),
        None,
    )
    if chapter is None:
        raise RecipeContractError(f"Unknown scenario ID: {scenario_id}")

    execution = chapter["execution"]
    if execution["kind"] != expected_kind:
        raise RecipeContractError(
            f"{scenario_id} uses {execution['kind']}, not {expected_kind}."
        )
    if execution["path"] != runner_path:
        raise RecipeContractError(
            f"{scenario_id} is owned by {execution['path']}, not {runner_path}."
        )
    expected_outputs = execution["expected_outputs"]
    if not expected_outputs or len(expected_outputs) != len(set(expected_outputs)):
        raise RecipeContractError(
            f"{scenario_id} must declare one or more unique outputs."
        )
    if any(Path(name).name != name for name in expected_outputs):
        raise RecipeContractError(
            f"{scenario_id} outputs must be filename-only paths."
        )
    return chapter


def extract_executable_block(chapter: dict[str, Any], *, language: str) -> str:
    scenario_id = chapter["id"]
    destination = (REPO_ROOT / chapter["destination"]).resolve()
    if not destination.is_relative_to(REPO_ROOT) or not destination.is_file():
        raise RecipeContractError(
            f"{scenario_id} destination is missing: {chapter['destination']}"
        )

    source = destination.read_text(encoding="utf-8")
    start = f"<!-- executable:{scenario_id}:start -->"
    end = f"<!-- executable:{scenario_id}:end -->"
    if source.count(start) != 1 or source.count(end) != 1:
        raise RecipeContractError(
            f"{scenario_id} must own exactly one executable marker pair."
        )
    marked = source.split(start, 1)[1].split(end, 1)[0]
    match = re.fullmatch(
        rf"\s*```{re.escape(language)}\r?\n(?P<code>.*?)\r?\n```\s*",
        marked,
        re.DOTALL,
    )
    if match is None:
        raise RecipeContractError(
            f"{scenario_id} marker must contain one {language} code block."
        )
    return match.group("code").strip() + "\n"


def copy_declared_inputs(
    chapter: dict[str, Any],
    *,
    recipe_source: str,
    workdir: Path,
) -> tuple[set[str], list[dict[str, Any]]]:
    fixture_manifest = json.loads(FIXTURE_MANIFEST_PATH.read_text(encoding="utf-8"))
    declared_file_ids: list[str] = []
    for fixture_id in chapter["fixtures"]:
        try:
            fixture = fixture_manifest["fixtures"][fixture_id]
        except KeyError as exc:
            raise RecipeContractError(
                f"{chapter['id']} references unknown fixture {fixture_id}."
            ) from exc
        declared_file_ids.extend(fixture.get("fileIds", ()))
        declared_file_ids.extend(fixture.get("fileReferences", ()))

    declared_files = {
        Path(fixture_manifest["files"][file_id]["relativePath"]).name: (
            fixture_manifest["files"][file_id]
        )
        for file_id in declared_file_ids
    }
    mentioned_inputs = set(
        re.findall(
            r"(?<![A-Za-z0-9_.-])([A-Za-z0-9_.-]+\.(?:gbk|gb|tsv|gff3|fna|fasta))",
            recipe_source,
        )
    )
    mentioned_inputs.difference_update(chapter["execution"]["expected_outputs"])
    undeclared = mentioned_inputs - set(declared_files)
    if undeclared:
        raise RecipeContractError(
            f"{chapter['id']} uses undeclared fixture files: {sorted(undeclared)}"
        )

    copied: set[str] = set()
    used_entries: list[dict[str, Any]] = []
    for filename in sorted(mentioned_inputs):
        entry = declared_files[filename]
        source = FIXTURE_ROOT / entry["relativePath"]
        if not source.is_file():
            raise RecipeContractError(f"Missing public tutorial fixture: {source}")
        payload = source.read_bytes()
        digest = hashlib.sha256(payload).hexdigest()
        if len(payload) != entry["sizeBytes"] or digest != entry["sha256"]:
            raise RecipeContractError(
                f"Fixture manifest mismatch for {entry['relativePath']}."
            )
        shutil.copyfile(source, workdir / filename)
        copied.add(filename)
        used_entries.append(entry)

    if not used_entries:
        raise RecipeContractError(f"{chapter['id']} does not use a declared fixture.")
    return copied, used_entries


def assert_exact_workdir_files(
    chapter: dict[str, Any],
    *,
    workdir: Path,
    copied_inputs: set[str],
) -> None:
    expected = copied_inputs | set(chapter["execution"]["expected_outputs"])
    actual = {
        path.relative_to(workdir).as_posix()
        for path in workdir.rglob("*")
        if path.is_file()
    }
    if actual != expected:
        raise RecipeContractError(
            f"{chapter['id']} created unexpected files: "
            f"expected {sorted(expected)}, got {sorted(actual)}"
        )


@dataclass(frozen=True)
class StandardSvgEvidence:
    """Semantic evidence extracted from one validated standard SVG."""

    size_bytes: int
    text_content: str
    text_nodes: frozenset[str]
    record_ids: frozenset[str]
    feature_ids: frozenset[str]
    slot_renderers: frozenset[str]


def inspect_standard_svg(
    chapter: dict[str, Any],
    *,
    output_path: Path,
) -> StandardSvgEvidence:
    size_bytes = output_path.stat().st_size
    if size_bytes < 10_000:
        raise RecipeContractError(f"{chapter['id']} produced an implausibly small SVG.")
    try:
        root = ElementTree.parse(output_path).getroot()
    except ElementTree.ParseError as exc:
        raise RecipeContractError(f"{chapter['id']} produced invalid XML.") from exc
    if _local_name(root.tag) != "svg" or "viewBox" not in root.attrib:
        raise RecipeContractError(f"{chapter['id']} output is not a complete SVG.")

    elements = list(root.iter())
    for element in elements:
        if _local_name(element.tag) == "script":
            raise RecipeContractError(f"{chapter['id']} standard SVG contains a script.")
        for raw_name, raw_value in element.attrib.items():
            name = _local_name(raw_name).lower()
            value = raw_value.strip().lower()
            if name.startswith("on"):
                raise RecipeContractError(
                    f"{chapter['id']} standard SVG contains an event handler."
                )
            if name == "href" and (
                value.startswith(("http:", "https:", "javascript:"))
                or "data:text/html" in value
            ):
                raise RecipeContractError(
                    f"{chapter['id']} standard SVG contains an unsafe external link."
                )
    if root.attrib.get("data-gbdraw-interactive-svg") == "true":
        raise RecipeContractError(
            f"{chapter['id']} promised standard SVG but produced interactive SVG."
        )

    text_content = " ".join(
        text.strip() for text in root.itertext() if text.strip()
    )
    text_nodes = frozenset(
        "".join(element.itertext()).strip()
        for element in elements
        if _local_name(element.tag) == "text"
    )
    record_ids = frozenset(
        element.attrib["data-gbdraw-record-id"]
        for element in elements
        if "data-gbdraw-record-id" in element.attrib
    )
    feature_ids = frozenset(
        element.attrib["data-gbdraw-feature-id"]
        for element in elements
        if "data-gbdraw-feature-id" in element.attrib
    )
    slot_renderers = frozenset(
        element.attrib["data-gbdraw-slot-renderer"]
        for element in elements
        if "data-gbdraw-slot-renderer" in element.attrib
    )
    return StandardSvgEvidence(
        size_bytes=size_bytes,
        text_content=text_content,
        text_nodes=text_nodes,
        record_ids=record_ids,
        feature_ids=feature_ids,
        slot_renderers=slot_renderers,
    )


def validate_standard_svg(
    chapter: dict[str, Any],
    *,
    output_path: Path,
    used_entries: list[dict[str, Any]],
) -> None:
    evidence = inspect_standard_svg(chapter, output_path=output_path)
    annotated_records = [
        record
        for entry in used_entries
        if entry.get("inputType") in {"genbank", "gff3"}
        for record in entry.get("records", ())
    ]
    if not annotated_records:
        raise RecipeContractError(
            f"{chapter['id']} validation expects annotated GenBank or GFF3 records."
        )
    expected_record_ids = {record["id"] for record in annotated_records}
    if len(expected_record_ids) != len(annotated_records):
        raise RecipeContractError(f"{chapter['id']} fixture record IDs are not unique.")
    for record in annotated_records:
        if (
            record["id"] not in evidence.text_content
            or f"{record['length']:,} bp" not in evidence.text_content
        ):
            raise RecipeContractError(
                f"{chapter['id']} output is missing the accession or length for "
                f"{record['id']}."
            )
    if evidence.record_ids != expected_record_ids:
        raise RecipeContractError(
            f"{chapter['id']} output record metadata does not match its inputs."
        )

    expected_features = 0
    for record in annotated_records:
        record_feature_count = record.get("displayedFeatureCount")
        if record_feature_count is None:
            record_feature_count = record.get("featureCounts", {}).get("CDS")
        if not isinstance(record_feature_count, int):
            raise RecipeContractError(
                f"{chapter['id']} fixture lacks a displayed CDS count for "
                f"{record['id']}."
            )
        expected_features += record_feature_count
    if len(evidence.feature_ids) != expected_features:
        raise RecipeContractError(
            f"{chapter['id']} expected {expected_features} displayed features, "
            f"found {len(evidence.feature_ids)}."
        )

    if chapter["id"] in {"T-CLI-01", "T-PY-01"}:
        required = {"ticks", "dinucleotide_content", "dinucleotide_skew"}
        if not required <= evidence.slot_renderers:
            raise RecipeContractError(
                f"{chapter['id']} is missing ticks, GC content, or GC skew."
            )
        if chapter["id"] == "T-CLI-01":
            gene_labels = {
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
            product_labels = {
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
            if not gene_labels <= evidence.text_nodes:
                raise RecipeContractError(
                    "T-CLI-01 output is missing one or more CDS gene labels."
                )
            if product_labels & evidence.text_nodes:
                raise RecipeContractError(
                    "T-CLI-01 output uses product descriptions instead of gene labels."
                )
            if "Homo sapiens" not in evidence.text_content:
                raise RecipeContractError(
                    "T-CLI-01 output is missing the requested species label."
                )
    elif chapter["id"] == "T-CLI-02":
        if not {"5 kbp", "45 kbp"} <= evidence.text_nodes:
            raise RecipeContractError("T-CLI-02 output is missing ruler labels.")
        if not {"A", "B", "J", "int"} <= evidence.text_nodes:
            raise RecipeContractError("T-CLI-02 output is missing concise gene labels.")


def publish_output(
    chapter: dict[str, Any],
    *,
    generated_path: Path,
    output_root: Path,
    check: bool,
    comparison_payload: Callable[[Path], bytes] | None = None,
) -> Path:
    destination = output_root / chapter["id"].lower() / generated_path.name
    if check:
        if not destination.is_file():
            raise RecipeContractError(f"Published artifact is missing: {destination}")
        read_payload = comparison_payload or Path.read_bytes
        if read_payload(destination) != read_payload(generated_path):
            display_path = (
                destination.relative_to(REPO_ROOT)
                if destination.is_relative_to(REPO_ROOT)
                else destination
            )
            raise RecipeContractError(
                f"Published artifact is stale: {display_path}"
            )
        return destination

    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(generated_path, destination)
    return destination


def _local_name(name: str) -> str:
    return name.rsplit("}", 1)[-1]
