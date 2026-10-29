"""Harness for the warm-equals-cold render-stage tests (STEP2 design 3.2).

The harness runs the CLI entry (``gbdraw.circular.circular_main`` or
``gbdraw.linear.linear_main``) with ``render_request`` replaced in the CLI
module by a function that keeps the typed request and stops the CLI. It then
renders that request through ``render_request`` inside a prepared-input cache
transaction, as the Web Worker does, with a test resource identity for every
input file named on the command line.

A render's reply is the SVG text of every Result plus the metadata the Web
reply carries (feature catalog, Legend rows, track-slot geometry), built by the
same production functions, the group manifest of each Result (keys, order,
positions, content revisions), and the remaining fields of the typed result
(interactive context, Linear metadata, LOSAT entries, warnings and notices) in
a canonical JSON form. A render that fails is an outcome too: its error type and
message are compared, and only a user-input validation error
(``USER_INPUT_ERRORS``) may be a step's accepted outcome.
"""

from __future__ import annotations

import dataclasses
import enum
import hashlib
import math
from collections.abc import Callable, Iterator, Mapping, Sequence
from contextlib import contextmanager
from dataclasses import dataclass, field
from pathlib import Path
from typing import Literal, Protocol, Union, get_args
from unittest import mock

import numpy as np
import pandas as pd
import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

import gbdraw.circular as circular_cli
import gbdraw.linear as linear_cli
from gbdraw.api import render_request
from gbdraw.api.prepared import PreparedBiologicalInputCache, PreparedResourceIdentity
from gbdraw.api.request_render import (
    CircularBatchRenderResult,
    CurrentRequestArtifacts,
    DiagramRequest,
    RequestRenderResult,
)
from gbdraw.exceptions import LabelsDoNotFitError, ValidationError
from gbdraw.io.cli_tables import (
    TablePathDependency,
    read_circular_track_table,
    read_comparisons_table,
    read_conservation_table,
    read_records_table,
)
from gbdraw.legend.row_facts import collect_legend_row_facts
from svgwrite import Drawing

from gbdraw.render.group_content import group_manifest
from gbdraw.render.track_slot_metadata import (
    build_track_slot_geometry_run_metadata,
    collect_track_slot_geometry_records,
)
from gbdraw.render_stages import RenderStageName
from gbdraw.web_support.feature_catalog import (
    build_feature_catalog,
    build_feature_catalog_item,
)

Mode = Literal["circular", "linear"]
JsonValue = Union[None, bool, int, float, str, list["JsonValue"], dict[str, "JsonValue"]]

STAGE_NAMES: frozenset[str] = frozenset(get_args(RenderStageName))

# The errors that reject a user's input. Some variants of the universe and some
# walk combinations are invalid input on the base branch too: a slot that does
# not fit inside (``ValidationError`` and its ``_InsideFitError``), a
# ``--multi_record_position`` list of the wrong length (``ValidationError``) and
# labels that do not fit (``LabelsDoNotFitError``, #1049). Every other error is
# a defect, even when warm and cold raise it alike.
USER_INPUT_ERRORS: tuple[type[Exception], ...] = (ValidationError, LabelsDoNotFitError)
_CLI_MODULES = {"circular": circular_cli, "linear": linear_cli}
_CLI_MAINS: dict[Mode, Callable[[list[str]], None]] = {
    "circular": circular_cli.circular_main,
    "linear": linear_cli.linear_main,
}


class _RequestCaptured(Exception):
    """Stops the CLI once its typed request is complete."""


@dataclass(frozen=True)
class CapturedRequest:
    """The typed request one CLI command line builds, with its input files."""

    mode: Mode
    request: DiagramRequest
    artifacts: CurrentRequestArtifacts | None
    input_paths: tuple[Path, ...]


@dataclass(frozen=True)
class StageCounts:
    """Per-stage counts of one render, read from its metrics."""

    builds: Mapping[str, int]
    hits: Mapping[str, int]
    uncacheable: Mapping[str, int]
    undeclared_reads: Mapping[str, int]
    build_total: int
    hit_total: int
    uncacheable_total: int
    undeclared_read_total: int


@dataclass(frozen=True)
class RenderOutcome:
    """The reply (or the error) and the metrics of one render."""

    reply: JsonValue
    error: tuple[str, str] | None
    metrics: Mapping[str, int] = field(compare=False)
    user_input_error: bool = field(default=False, compare=False)

    @property
    def stages(self) -> StageCounts:
        return stage_counts(self.metrics)


def capture_cli_request(mode: Mode, argv: Sequence[str]) -> CapturedRequest:
    """Run the CLI entry until it calls ``render_request`` and keep the request."""

    captured: list[tuple[DiagramRequest, CurrentRequestArtifacts | None]] = []

    def keep(
        request: DiagramRequest,
        *,
        artifacts: CurrentRequestArtifacts | None = None,
        include_feature_catalog: bool = False,
    ) -> RequestRenderResult | CircularBatchRenderResult:
        del include_feature_catalog  # the harness always builds the catalog
        captured.append((request, artifacts))
        raise _RequestCaptured

    with mock.patch.object(_CLI_MODULES[mode], "render_request", keep):
        try:
            _CLI_MAINS[mode](list(argv))
        except _RequestCaptured:
            pass
    if len(captured) != 1:
        raise AssertionError(
            f"The {mode} CLI did not call render_request exactly once "
            f"({len(captured)} calls) for {list(argv)!r}."
        )
    request, artifacts = captured[0]
    return CapturedRequest(
        mode=mode,
        request=request,
        artifacts=artifacts,
        input_paths=input_files(argv),
    )


class _NamesFiles(Protocol):
    @property
    def path_dependencies(self) -> tuple[TablePathDependency, ...]: ...


# The table options whose tables name other input files, each with its reader;
# the reader's ``path_dependencies`` lists those files. A Worker stages them, so
# the harness does too: otherwise every stage that reads them is uncacheable.
TABLE_FILE_READERS: Mapping[str, Callable[[str], _NamesFiles]] = {
    "--records_table": read_records_table,
    "--comparisons_table": read_comparisons_table,
    "--conservation_table": read_conservation_table,
    "--circular_track_table": read_circular_track_table,
}
# The table options whose tables hold values only and name no other file.
VALUE_TABLE_OPTIONS = frozenset(
    {
        "--annotation_table",
        "--feature_override_table",
        "--feature_placement_table",
        "--feature_visibility_table",
        "--label_table",
    }
)


def input_files(argv: Sequence[str]) -> tuple[Path, ...]:
    """Every command-line token that names an existing file, and the files a
    table option's table names (``TABLE_FILE_READERS``), once, in order."""

    paths: dict[Path, None] = {}
    for previous, token in zip(("", *argv), argv):
        if token.startswith("-"):
            continue
        path = Path(token)
        if not path.is_file():
            continue
        paths[path.resolve()] = None
        reader = TABLE_FILE_READERS.get(previous)
        for dependency in reader(str(path)).path_dependencies if reader is not None else ():
            if Path(dependency.path).is_file():
                paths[Path(dependency.path).resolve()] = None
    return tuple(paths)


def resource_identities(
    paths: Sequence[Path],
) -> dict[str | Path, PreparedResourceIdentity]:
    """A test resource identity per file: its path names it, its bytes stamp it."""

    identities: dict[str | Path, PreparedResourceIdentity] = {}
    for path in paths:
        content = path.read_bytes()
        identities[path] = PreparedResourceIdentity(
            resource_id="test-resource-" + hashlib.sha256(str(path).encode()).hexdigest()[:16],
            cache_token="test-token-" + hashlib.sha256(content).hexdigest()[:32],
            size=len(content),
        )
    return identities


def render_captured(
    captured: CapturedRequest,
    cache: PreparedBiologicalInputCache | None,
) -> RenderOutcome:
    """Render one captured request; inside a transaction when ``cache`` is given."""

    diagnostics: dict[str, object] = {"timingsMs": {}, "metrics": {}}
    reply: JsonValue = None
    error: tuple[str, str] | None = None
    user_input_error = False
    try:
        with _transaction(cache, captured.input_paths, diagnostics):
            reply = _reply(captured, diagnostics if cache is not None else None)
    except Exception as exc:  # noqa: BLE001 - a failed render is an outcome
        error = (type(exc).__name__, str(exc))
        user_input_error = isinstance(exc, USER_INPUT_ERRORS)
    return RenderOutcome(
        reply=reply,
        error=error,
        metrics=_metrics(diagnostics),
        user_input_error=user_input_error,
    )


@contextmanager
def _transaction(
    cache: PreparedBiologicalInputCache | None,
    paths: Sequence[Path],
    diagnostics: dict[str, object],
) -> Iterator[None]:
    if cache is None:
        yield
        return
    with cache.transaction(
        resource_paths=resource_identities(paths),
        diagnostics=diagnostics,
    ):
        yield


def _reply(
    captured: CapturedRequest,
    diagnostics: dict[str, object] | None,
) -> JsonValue:
    # No label-fit warning search (as the Web): it is bounded by wall-clock
    # time, so the number of verify renders, and the stages they build, would
    # depend on machine load (OV-705).
    rendered = render_request(
        captured.request,
        artifacts=captured.artifacts,
        include_feature_catalog=True,
        label_fit_warn_search_seconds=0.0,
    )
    items: tuple[RequestRenderResult, ...] = (
        rendered.items if isinstance(rendered, CircularBatchRenderResult) else (rendered,)
    )
    results: list[JsonValue] = []
    catalog_items: list[dict[str, object]] = []
    legend_rows: list[JsonValue] = []
    geometry: list[dict[str, object]] = []
    for result_index, item in enumerate(items):
        svg_paths = [path for path in item.output_paths if Path(path).suffix == ".svg"]
        svg_texts = [Path(path).read_text(encoding="utf-8") for path in svg_paths]
        result_name = Path(svg_paths[0]).name if svg_paths else f"result-{result_index}"
        catalog_items.append(
            build_feature_catalog_item(
                svg_texts[0] if svg_texts else "",
                item.interactive_context,
                result_index=result_index,
                result_name=result_name,
                drawing=item.drawing,
                _diagnostics=diagnostics,
            )
        )
        legend_rows.append(
            canonical(
                collect_legend_row_facts(
                    item.drawing,
                    result_index=result_index,
                    result_name=result_name,
                )
            )
        )
        geometry.extend(
            collect_track_slot_geometry_records(
                item.drawing,
                result_index=result_index,
                result_name=result_name,
            )
        )
        results.append(
            {
                "svg": [canonical(text) for text in svg_texts],
                "interactiveContext": canonical(item.interactive_context),
                "linearMetadata": canonical(item.linear_metadata),
                "losatCacheEntries": canonical(item.losat_cache_entries),
                "losatDerivedCacheEntries": canonical(item.losat_derived_cache_entries),
                "proteinIdentityManifest": canonical(item.protein_identity_manifest),
                "annotationWarnings": canonical(item.annotation_warnings),
                "featureIdentityNotices": canonical(item.feature_identity_notices),
                "manifest": _manifest(item.drawing),
            }
        )
    return {
        "results": results,
        "batchLosatCacheEntries": canonical(
            rendered.losat_cache_entries
            if isinstance(rendered, CircularBatchRenderResult)
            else ()
        ),
        "featureCatalog": canonical(build_feature_catalog(catalog_items)),
        "legendRows": legend_rows,
        "geometry": canonical(
            build_track_slot_geometry_run_metadata(mode=captured.mode, records=geometry)
        ),
    }


def _manifest(drawing: Drawing) -> JsonValue:
    """The placed groups of a drawing: keys, order, positions and content revisions (design 2.4)."""

    return [
        {
            "key": canonical((entry.key.scope, entry.key.kind.value, entry.key.slot)),
            "domId": entry.dom_id,
            "position": entry.position,
            "revision": entry.content.revision,
            "definitions": list(entry.definition_ids),
        }
        for entry in group_manifest(drawing)
    ]


def _metrics(diagnostics: Mapping[str, object]) -> dict[str, int]:
    metrics = diagnostics.get("metrics")
    if not isinstance(metrics, Mapping):
        return {}
    return {
        str(name): int(value)
        for name, value in metrics.items()
        if isinstance(value, (int, float)) and not isinstance(value, bool)
    }


def stage_counts(metrics: Mapping[str, int]) -> StageCounts:
    """Per-stage build, hit, uncacheable and undeclared-read counts of a render."""

    def per_stage(kind: str) -> dict[str, int]:
        prefix = f"renderStage{kind}:"
        counts = {
            name[len(prefix):]: value
            for name, value in metrics.items()
            if name.startswith(prefix)
        }
        unknown = set(counts) - STAGE_NAMES
        if unknown:
            raise AssertionError(f"Unknown render stage names in metrics: {sorted(unknown)}")
        return counts

    return StageCounts(
        builds=per_stage("Build"),
        hits=per_stage("Hit"),
        uncacheable=per_stage("Uncacheable"),
        undeclared_reads=per_stage("UndeclaredRead"),
        build_total=metrics["renderStageBuildCount"],
        hit_total=metrics["renderStageHitCount"],
        uncacheable_total=metrics["renderStageUncacheableCount"],
        undeclared_read_total=metrics["renderStageUndeclaredReadCount"],
    )


def canonical(value: object) -> JsonValue:
    """A JSON form of a reply value for equality; it keeps order and types."""

    # Exact JSON types first: most of a reply is strings, numbers, lists and
    # dicts, and the ``isinstance`` chain below costs a few ABC checks each.
    kind = type(value)
    if kind is str or kind is int or kind is bool:
        return value
    if kind is float:
        return {"$float": "nan"} if math.isnan(value) else value
    if kind is list or kind is tuple:
        return [canonical(item) for item in value]
    if kind is dict:
        return {"$mapping": [[canonical(key), canonical(entry)] for key, entry in value.items()]}
    if value is None or isinstance(value, (bool, str)):
        return value
    if isinstance(value, enum.Enum):
        return {"$enum": f"{type(value).__name__}.{value.name}"}
    if isinstance(value, (int, np.integer)):
        return int(value)
    if isinstance(value, (float, np.floating)):
        number = float(value)
        return {"$float": "nan"} if math.isnan(number) else number
    if isinstance(value, (bytes, bytearray)):
        return {"$bytes": hashlib.sha256(bytes(value)).hexdigest()}
    if isinstance(value, Path):
        return {"$path": str(value)}
    if isinstance(value, pd.DataFrame):
        return {
            "$frame": canonical(
                {
                    "columns": [str(column) for column in value.columns],
                    "index": list(value.index),
                    "rows": value.to_numpy().tolist(),
                }
            )
        }
    if isinstance(value, pd.Series):
        return {"$series": canonical({"index": list(value.index), "values": value.tolist()})}
    if isinstance(value, np.ndarray):
        return {"$array": canonical(value.tolist())}
    if isinstance(value, SeqRecord):
        return {
            "$record": canonical(
                {
                    "id": value.id,
                    "sequence": hashlib.sha256(str(value.seq).encode()).hexdigest(),
                    "features": [repr(feature.location) + feature.type for feature in value.features],
                }
            )
        }
    if dataclasses.is_dataclass(value) and not isinstance(value, type):
        return {
            "$dataclass": type(value).__name__,
            "fields": {
                item.name: canonical(getattr(value, item.name))
                for item in dataclasses.fields(value)
            },
        }
    if isinstance(value, Mapping):
        return {
            "$mapping": [[canonical(key), canonical(entry)] for key, entry in value.items()]
        }
    if isinstance(value, (set, frozenset)):
        return {"$set": sorted((canonical(item) for item in value), key=repr)}
    if isinstance(value, Sequence):
        return [canonical(item) for item in value]
    raise TypeError(f"No canonical form for a reply value of type {type(value).__qualname__}.")


def first_difference(left: JsonValue, right: JsonValue, path: str = "$") -> str | None:
    """The path and the two values of the first difference, or None when equal."""

    if type(left) is not type(right):
        return f"{path}: {_short(left)} != {_short(right)}"
    if isinstance(left, dict) and isinstance(right, dict):
        for key in sorted(set(left) | set(right)):
            if key not in left or key not in right:
                return f"{path}.{key}: present on one side only"
            found = first_difference(left[key], right[key], f"{path}.{key}")
            if found:
                return found
        return None
    if isinstance(left, list) and isinstance(right, list):
        for index, (a, b) in enumerate(zip(left, right)):
            found = first_difference(a, b, f"{path}[{index}]")
            if found:
                return found
        if len(left) != len(right):
            return f"{path}: length {len(left)} != {len(right)}"
        return None
    if isinstance(left, str) and isinstance(right, str) and left != right:
        index = next(
            (i for i, (a, b) in enumerate(zip(left, right)) if a != b),
            min(len(left), len(right)),
        )
        return (
            f"{path} at character {index}: "
            f"{left[max(0, index - 80):index + 80]!r} != {right[max(0, index - 80):index + 80]!r}"
        )
    if left != right:
        return f"{path}: {_short(left)} != {_short(right)}"
    return None


def _short(value: JsonValue) -> str:
    text = repr(value)
    return text if len(text) <= 300 else text[:300] + "..."


def assert_same_outcome(warm: RenderOutcome, cold: RenderOutcome, context: str) -> None:
    """Fail with the first difference between a warm and a cold render."""

    if warm.error != cold.error:
        raise AssertionError(f"{context}: warm error {warm.error!r} != cold error {cold.error!r}")
    difference = first_difference(warm.reply, cold.reply)
    if difference is not None:
        raise AssertionError(f"{context}: warm reply differs from cold: {difference}")


class EquivalenceHarness:
    """One shared cache for a chain of renders, each compared with a fresh cache."""

    def __init__(self, mode: Mode, output_directory: Path) -> None:
        self.mode: Mode = mode
        self.output_directory = output_directory
        self.cache = PreparedBiologicalInputCache()
        self.empty = True
        output_directory.mkdir(parents=True, exist_ok=True)

    def argv(self, arguments: Sequence[str]) -> list[str]:
        return [
            *arguments,
            "-o", str(self.output_directory / "equivalence"),
            "-f", "svg",
            "--overwrite",
        ]

    def capture(self, arguments: Sequence[str]) -> CapturedRequest:
        return capture_cli_request(self.mode, self.argv(arguments))

    def warm(self, arguments: Sequence[str]) -> RenderOutcome:
        """Render on the shared cache (the state before is the previous render)."""

        outcome = render_captured(self.capture(arguments), self.cache)
        self.empty = False
        return outcome

    def cold(self, arguments: Sequence[str]) -> RenderOutcome:
        """Render on a new cache."""

        return render_captured(self.capture(arguments), PreparedBiologicalInputCache())

    def uncached(self, arguments: Sequence[str]) -> RenderOutcome:
        """Render without a transaction, as the CLI and the Python API do."""

        return render_captured(self.capture(arguments), None)

    def step(self, arguments: Sequence[str], context: str) -> tuple[RenderOutcome, RenderOutcome]:
        """Render warm then cold and assert they are equal."""

        warm = self.warm(arguments)
        cold = self.cold(arguments)
        assert_same_outcome(warm, cold, context)
        return warm, cold


# --- Option universe (the DEPMAP cases C1-C3, L1-L3 and their variants) -------------------------

REPOSITORY = Path(__file__).resolve().parents[2]
TEST_INPUTS = REPOSITORY / "tests" / "test_inputs"
EXAMPLES = REPOSITORY / "examples"


@dataclass(frozen=True)
class Variant:
    """One value of an option: tokens removed from the base, then tokens added."""

    label: str
    add: tuple[str, ...]
    remove: tuple[tuple[str, ...], ...] = ()


@dataclass(frozen=True)
class OptionSpec:
    """One option of the universe, with the prerequisites both states share."""

    option_id: str
    modes: frozenset[Mode]
    variants: tuple[Variant, ...]
    pre: tuple[str, ...] = ()
    pre_remove: tuple[tuple[str, ...], ...] = ()
    cases: frozenset[str] | None = None

    @property
    def key(self) -> str:
        """Unique within the universe: an option may have one entry per mode."""

        cases = ",".join(sorted(self.cases)) if self.cases else "all"
        return f"{self.option_id}[{''.join(sorted(m[0] for m in self.modes))}:{cases}]"


@dataclass(frozen=True)
class Case:
    case_id: str
    mode: Mode
    base: tuple[str, ...]
    tables: Mapping[str, str]


_CL: frozenset[Mode] = frozenset({"circular", "linear"})
_C: frozenset[Mode] = frozenset({"circular"})
_L: frozenset[Mode] = frozenset({"linear"})


def _variants(*values: Sequence[str]) -> tuple[Variant, ...]:
    return tuple(Variant(" ".join(value) or "(none)", tuple(value)) for value in values)


def _opt(
    option_id: str,
    modes: frozenset[Mode],
    *values: Sequence[str],
    pre: Sequence[str] = (),
    pre_remove: Sequence[Sequence[str]] = (),
    cases: Sequence[str] | None = None,
) -> OptionSpec:
    return OptionSpec(
        option_id=option_id,
        modes=modes,
        variants=_variants(*values),
        pre=tuple(pre),
        pre_remove=tuple(tuple(tokens) for tokens in pre_remove),
        cases=frozenset(cases) if cases is not None else None,
    )


# Small tables, written into a temporary directory by ``write_tables``.
_TABLE_TEXT: Mapping[str, str] = {
    "ann_br_hmmt.tsv": "set_id\tid\tmark\trecord\tstart\tend\tlabel\ns2\tb1\tbracket\tNC_012920.1\t8000\t12000\tB\n",
    "ann_br_mj.tsv": "set_id\tid\tmark\trecord\tstart\tend\tlabel\ns2\tb1\tbracket\tLC738868.1\t80000\t120000\tB\n",
    "ann_hi_hmmt.tsv": "set_id\tid\tmark\trecord\tstart\tend\tlabel\ns1\ta1\thighlight\tNC_012920.1\t1000\t5000\tA\n",
    "ann_hi_mj.tsv": "set_id\tid\tmark\trecord\tstart\tend\tlabel\ns1\ta1\thighlight\tLC738868.1\t10000\t50000\tA\n",
    "col_hmmt.tsv": "CDS\tproduct\t^NADH.*\t#ff0000\tNADH\n",
    "default_colors.tsv": "CDS\t#D0D0D0\n",
    "label_hmmt.tsv": "*\tCDS\tgene\t^ND1$\tNADH one\n*\tCDS\tproduct\t.*subunit 1$\tNADH one\n",
    "label_mj.tsv": "*\tCDS\tproduct\t^ubiquitin$\tUbq\n",
    "ovr_hmmt.tsv": (
        "record\tfeature_selector\tfeature_visibility\tlabel_visibility\tlabel_text\n"
        "NC_012920.1\tproduct=NADH dehydrogenase subunit 1\toff\t\t\n"
        "NC_012920.1\tproduct=NADH dehydrogenase subunit 2\t\t\tRenamed label\n"
    ),
    "ovr_mj.tsv": (
        "record\tfeature_selector\tfeature_visibility\tlabel_visibility\tlabel_text\n"
        "LC738868.1\tproduct=ubiquitin\t\t\tRenamed label\n"
    ),
    "place_hmmt.tsv": "record\tfeature_selector\tplacement\nNC_012920.1\tproduct=NADH dehydrogenase subunit 1\toutward\n",
    "place_hmmt_l.tsv": "record\tfeature_selector\tplacement\nNC_012920.1\tproduct=NADH dehydrogenase subunit 1\tabove\n",
    "place_mj.tsv": "record\tfeature_selector\tplacement\nLC738868.1\tproduct=ubiquitin\toutward\n",
    "place_mj_l.tsv": "record\tfeature_selector\tplacement\nLC738868.1\tproduct=ubiquitin\tabove\n",
    "prio.tsv": "CDS\tproduct,gene\n",
    "prio_gene.tsv": "CDS\tgene,product\n",
    "tracks_hmmt.tsv": (
        "id\trenderer\tside\tr\tw\tparams\nfeatures\tfeatures\taxis\t\t\t\n"
        "gc_content\tdinucleotide_content\tinside\t\t0.1\tnt=GC\n"
        "gc_skew\tdinucleotide_skew\tinside\t\t0.1\tnt=GC\n"
        "a_skew_2\tdinucleotide_skew\tinside\t\t0.1\t"
        "nt=AT,positive_color=#deaf6e,negative_color=#7294e3,legend_label=AT skew\n"
        "ticks\tticks\tinside\t\t\ttick_label_layout=label_in_tick_out\n"
    ),
    "vis_hmmt.tsv": "*\tCDS\tgene\t^ND[1-3]$\toff\n",
    "vis_mj.tsv": "*\tCDS\tproduct\t^hypothetical protein$\toff\n",
    "wl_hmmt.tsv": "CDS\tproduct\t^NADH dehydrogenase subunit [1-3]$\n",
    "wl_mj.tsv": "CDS\tproduct\t^(ubiquitin|tyrosine recombinase)$\n",
    "comp_l3.tsv": f"blast\tquery\tsubject\n{EXAMPLES}/MelaMJNV.PemoMJNVA.tblastx.out\t#2\t#3\n",
}


def write_depth_table(path: Path, records: Sequence[tuple[str, int]], *, offset: int = 60) -> Path:
    """A deterministic depth table: one row per position of each record."""

    with path.open("w", encoding="utf-8") as handle:
        for record_id, length in records:
            for position in range(1, length + 1):
                depth = int(offset + 40 * math.sin(position / 900.0) + (position % 7))
                handle.write(f"{record_id}\t{position}\t{depth}\n")
    return path


def write_tables(directory: Path) -> Path:
    """Write the universe's tables into ``directory`` and return it."""

    directory.mkdir(parents=True, exist_ok=True)
    for name, text in _TABLE_TEXT.items():
        (directory / name).write_text(text, encoding="utf-8")
    write_depth_table(
        directory / "depth_hmmt.tsv",
        (("NC_012920.1", 16569), ("NC_019843.3", 30119), ("NC_001416.1", 48502)),
    )
    return directory


def _hm_tables() -> dict[str, str]:
    return {
        "VIS": "{IN}/vis_hmmt.tsv", "LBL": "{IN}/label_hmmt.tsv", "WL": "{IN}/wl_hmmt.tsv",
        "PRIO": "{IN}/prio.tsv", "OVR": "{IN}/ovr_hmmt.tsv", "COL": "{IN}/col_hmmt.tsv",
        "ANN_HI": "{IN}/ann_hi_hmmt.tsv", "ANN_BR": "{IN}/ann_br_hmmt.tsv",
        "PLACE": "{IN}/place_hmmt.tsv", "PLACE_L": "{IN}/place_hmmt_l.tsv",
        "DEPTH": "{IN}/depth_hmmt.tsv", "BL": "NADH,tRNA",
    }


def _mj_tables() -> dict[str, str]:
    return {
        "VIS": "{IN}/vis_mj.tsv", "LBL": "{IN}/label_mj.tsv", "WL": "{IN}/wl_mj.tsv",
        "PRIO": "{IN}/prio.tsv", "OVR": "{IN}/ovr_mj.tsv", "COL": "{E}/color_table.txt",
        "ANN_HI": "{IN}/ann_hi_mj.tsv", "ANN_BR": "{IN}/ann_br_mj.tsv",
        "PLACE": "{IN}/place_mj.tsv", "PLACE_L": "{IN}/place_mj_l.tsv",
        "DEPTH": "{T}/MjeNMV.DRR271272.depth.tsv",
        "BL": "hypothetical,ubiquitin",
    }


_TITLES = ("--plot_title", "Plot title")
CASES: Mapping[str, Case] = {
    "C1": Case("C1", "circular", (
        "--gbk", "{T}/HmmtDNA.gbk", "--labels", "out", "-l", "right", *_TITLES,
        "--plot_title_position", "top", "--gc", "--skew"), _hm_tables()),
    # C2 resolves overlaps: without it its embedded labels do not fit and the
    # base ends in LabelsDoNotFitError (OV-706).
    "C2": Case("C2", "circular", (
        "--gbk", "{T}/HmmtDNA.gbk", "{T}/MERS-CoV.gbk", "{T}/NC_001416.gb", "--multi_record_canvas",
        "--labels", "out", "-l", "right", *_TITLES, "--plot_title_position", "top",
        "--gc", "--skew", "--resolve_overlaps"), _hm_tables()),
    "C3": Case("C3", "circular", (
        "--gbk", "{E}/MjeNMV.gb", "--labels", "out", "-l", "right", *_TITLES,
        "--plot_title_position", "top", "--gc", "--skew",
        "--depth_track", "{DEPTH}", "--depth_track_label", "Depth",
        "--conservation_blast", "{E}/MjeNMV.MelaMJNV.tblastx.out", "--conservation_labels", "MelaMJNV",
        "--conservation_reference", "query"), _mj_tables()),
    "L1": Case("L1", "linear", (
        "--gbk", "{T}/HmmtDNA.gbk", "--show_labels", "all", "-l", "right", *_TITLES,
        "--gc", "--skew", "--scale_style", "ruler"), _hm_tables()),
    "L2": Case("L2", "linear", (
        "--gbk", "{E}/MjeNMV.gb", "{E}/MelaMJNV.gb", "-b", "{E}/MjeNMV.MelaMJNV.tblastx.out",
        "--show_labels", "all", "-l", "right", *_TITLES, "--gc", "--skew",
        "--scale_style", "ruler",
        "--depth_track", "{T}/MjeNMV.DRR271272.depth.tsv", "{T}/MelaMJNV.DRR394934.depth.tsv",
        "--depth_track_label", "Depth"), _mj_tables()),
    "L3": Case("L3", "linear", (
        "--gbk", "{E}/MjeNMV.gb", "{E}/MelaMJNV.gb", "{E}/PemoMJNVA.gb",
        "-b", "{E}/MjeNMV.MelaMJNV.tblastx.out", "{E}/MelaMJNV.PemoMJNVA.tblastx.out",
        "--show_labels", "all", "-l", "right", *_TITLES, "--gc", "--skew",
        "--scale_style", "ruler"), _mj_tables()),
}

_PERCENT = ("--gc_content_mode", "percent")
_DEPTH_CASES = ("C3", "L2")
_ROW2 = ("--multi_record_position", "#1@1", "--multi_record_position", "#2@1",
         "--multi_record_position", "#3@2")
_ROWPRE = (*_ROW2, "--comparisons_table", "{IN}/comp_l3.tsv")
_ROWRM = (("-b", "{E}/MjeNMV.MelaMJNV.tblastx.out", "{E}/MelaMJNV.PemoMJNVA.tblastx.out"),)

UNIVERSE: tuple[OptionSpec, ...] = (
    # Both modes.
    _opt("palette", _CL, ["-p", "orchid"], ["-p", "soft_pastels"]),
    _opt("table", _CL, ["-t", "{COL}"]),
    _opt("default_colors", _CL, ["-d", "{IN}/default_colors.tsv"]),
    _opt("nt", _CL, ["-n", "AT"]),
    _opt("window", _CL, ["-w", "500"], ["-w", "5000"]),
    _opt("step", _CL, ["-s", "50"], ["-s", "2000"]),
    _opt("features", _CL, ["-k", "CDS"],
         ["-k", "CDS,rRNA,tRNA,tmRNA,ncRNA,misc_RNA,repeat_region,gene"]),
    _opt("feature_shape", _CL, ["--feature_shape", "CDS=rectangle"],
         ["--feature_shape", "CDS=underlay"]),
    _opt("arrow_head_length_ratio", _CL, ["--arrow_head_length_ratio", "0.3"],
         ["--arrow_head_length_ratio", "1.5"]),
    _opt("arrow_shaft_width_ratio", _CL, ["--arrow_shaft_width_ratio", "0.5"]),
    _opt("block_stroke_color", _CL, ["--block_stroke_color", "red"]),
    _opt("block_stroke_width", _CL, ["--block_stroke_width", "0"], ["--block_stroke_width", "8"]),
    _opt("axis_stroke_color", _CL, ["--axis_stroke_color", "red"]),
    _opt("axis_stroke_width", _CL, ["--axis_stroke_width", "1"], ["--axis_stroke_width", "12"]),
    _opt("line_stroke_color", _CL, ["--line_stroke_color", "red"]),
    _opt("line_stroke_width", _CL, ["--line_stroke_width", "1"], ["--line_stroke_width", "12"]),
    _opt("definition_font_size", _CL, ["--definition_font_size", "10"],
         ["--definition_font_size", "40"]),
    _opt("plot_title", _CL, ["--plot_title", "Another title text"]),
    _opt("plot_title_font_size", _CL, ["--plot_title_font_size", "16"],
         ["--plot_title_font_size", "80"]),
    _opt("plot_title_position", _CL, ["--plot_title_position", "bottom"],
         ["--plot_title_position", "none"], cases=["C1", "C2", "C3"]),
    _opt("plot_title_position", _L, ["--plot_title_position", "top"],
         ["--plot_title_position", "center"], cases=["L1", "L2"]),
    _opt("label_font_size", _CL, ["--label_font_size", "8"], ["--label_font_size", "30"]),
    OptionSpec("gc", _CL, (Variant("--no-gc", ("--no-gc",), (("--gc",),)),)),
    OptionSpec("skew", _CL, (Variant("--no-skew", ("--no-skew",), (("--skew",),)),)),
    _opt("gc_content_mode", _CL, ["--gc_content_mode", "percent"]),
    _opt("gc_content_min_percent", _CL, ["--gc_content_min_percent", "20"], pre=_PERCENT),
    _opt("gc_content_max_percent", _CL, ["--gc_content_max_percent", "60"], pre=_PERCENT),
    _opt("gc_content_large_tick_interval", _CL, ["--gc_content_large_tick_interval", "10"],
         pre=_PERCENT),
    _opt("gc_content_small_tick_interval", _CL, ["--gc_content_small_tick_interval", "5"],
         pre=_PERCENT),
    _opt("gc_content_tick_font_size", _CL, ["--gc_content_tick_font_size", "24"], pre=_PERCENT),
    _opt("show_gc_content_axis", _CL, ["--show_gc_content_axis"], ["--hide_gc_content_axis"],
         pre=_PERCENT),
    _opt("show_gc_content_ticks", _CL, ["--show_gc_content_ticks"], ["--hide_gc_content_ticks"],
         pre=_PERCENT),
    _opt("depth_track (add)", _CL, ["--depth_track", "{DEPTH}", "--depth_track_label", "Depth"],
         cases=["C1", "C2", "L1"]),
    _opt("depth_track_label", _CL, ["--depth_track_label", "Cov"], cases=_DEPTH_CASES),
    _opt("depth_track_color", _CL, ["--depth_track_color", "red"], cases=_DEPTH_CASES),
    _opt("depth_track_large_tick_interval", _CL, ["--depth_track_large_tick_interval", "20"],
         cases=_DEPTH_CASES),
    _opt("depth_track_small_tick_interval", _CL, ["--depth_track_small_tick_interval", "10"],
         cases=_DEPTH_CASES),
    _opt("depth_track_tick_font_size", _CL, ["--depth_track_tick_font_size", "24"],
         cases=_DEPTH_CASES),
    _opt("depth_color", _CL, ["--depth_color", "#ff0000"], cases=_DEPTH_CASES),
    _opt("depth_window", _CL, ["--depth_window", "3000"], cases=_DEPTH_CASES),
    _opt("depth_step", _CL, ["--depth_step", "3000"], cases=_DEPTH_CASES),
    _opt("share_depth_axis", _CL, ["--share_depth_axis"], cases=_DEPTH_CASES),
    _opt("depth_min", _CL, ["--depth_min", "100"], cases=_DEPTH_CASES),
    _opt("depth_max", _CL, ["--depth_max", "800"], cases=_DEPTH_CASES),
    _opt("depth_log_scale", _CL, ["--depth_log_scale"], cases=_DEPTH_CASES),
    _opt("show_depth_axis", _CL, ["--show_depth_axis"], ["--hide_depth_axis"], cases=_DEPTH_CASES),
    _opt("show_depth_ticks", _CL, ["--show_depth_ticks"], ["--hide_depth_ticks"], cases=_DEPTH_CASES),
    _opt("depth_large_tick_interval", _CL, ["--depth_large_tick_interval", "500"], cases=_DEPTH_CASES),
    _opt("depth_small_tick_interval", _CL, ["--depth_small_tick_interval", "100"], cases=_DEPTH_CASES),
    _opt("depth_tick_font_size", _CL, ["--depth_tick_font_size", "30"], cases=_DEPTH_CASES),
    _opt("legend", _CL, ["-l", "left"], ["-l", "bottom"], ["-l", "none"], ["-l", "top"]),
    _opt("legend", _C, ["-l", "upper_left"], ["-l", "lower_right"]),
    _opt("legend_box_size", _CL, ["--legend_box_size", "10"], ["--legend_box_size", "50"]),
    _opt("legend_font_size", _CL, ["--legend_font_size", "10"], ["--legend_font_size", "36"]),
    _opt("separate_strands", _CL, ["--separate_strands"]),
    _opt("resolve_overlaps", _CL, ["--resolve_overlaps"], cases=["C1", "C3", "L1", "L2", "L3"]),
    _opt("feature_overlap_tolerance_bp", _CL, ["--feature_overlap_tolerance_bp", "200"],
         pre=["--resolve_overlaps"]),
    _opt("label_rendering", _CL, ["--label_rendering", "embedded_only"],
         ["--label_rendering", "external_only"]),
    _opt("label_whitelist", _CL, ["--label_whitelist", "{WL}"]),
    _opt("label_blacklist", _CL, ["--label_blacklist", "{BL}"]),
    _opt("qualifier_priority", _CL, ["--qualifier_priority", "{PRIO}"],
         ["--qualifier_priority", "{IN}/prio_gene.tsv"]),
    _opt("label_table", _CL, ["--label_table", "{LBL}"]),
    _opt("feature_visibility_table", _CL, ["--feature_visibility_table", "{VIS}"]),
    _opt("feature_override_table", _CL, ["--feature_override_table", "{OVR}"]),
    _opt("feature_placement_table", _C, ["--feature_placement_table", "{PLACE}"],
         pre=["--separate_strands", "--circular_track_slot",
              "features:features@side=overlay,lane_direction=split"]),
    # Linear mode places on the sides above and below.
    _opt("feature_placement_table", _L, ["--feature_placement_table", "{PLACE_L}"]),
    _opt("annotation_table", _CL, ["--annotation_table", "{ANN_HI}"],
         ["--annotation_table", "{ANN_BR}"]),
    _opt("scale_interval", _CL, ["--scale_interval", "4000"], ["--scale_interval", "50000"]),
    _opt("hide_scale", _CL, ["--hide_scale"]),
    _opt("record_topology", _CL, ["--record_topology", "linear"], cases=["C1", "L1"]),
    _opt("display_start_coordinate", _CL, ["--display_start_coordinate", "5000"],
         cases=["C1", "L1"]),
    # Circular only.
    _opt("species", _C, ["--species", "<i>Homo sapiens</i>"]),
    _opt("strain", _C, ["--strain", "K-12"]),
    _opt("species (no plot title)", _C, ["--species", "<i>Homo sapiens</i>"],
         pre=["--plot_title_position", "none"]),
    _opt("strain (no plot title)", _C, ["--strain", "K-12"], pre=["--plot_title_position", "none"]),
    _opt("keep_full_definition_with_plot_title", _C, ["--keep_full_definition_with_plot_title"],
         pre=["--species", "<i>Homo sapiens</i>"]),
    _opt("center_reserved_radius", _C, ["--center_reserved_radius", "0"],
         ["--center_reserved_radius", "400"], ["--center_reserved_radius", "250"]),
    _opt("tick_label_font_size", _C, ["--tick_label_font_size", "8"],
         ["--tick_label_font_size", "30"]),
    _opt("circular_label_spacing", _C, ["--circular_label_spacing", "40"]),
    _opt("labels", _C, ["--labels", "none"], ["--labels", "both"]),
    _opt("label_placement", _C, ["--label_placement", "radial"]),
    *(
        _opt(name, _C, [f"--{name}", "0.8"], [f"--{name}", "1.3"],
             pre=["--labels", "both"] if name.startswith("inner") else ())
        for name in (
            "outer_label_x_radius_offset", "outer_label_y_radius_offset",
            "inner_label_x_radius_offset", "inner_label_y_radius_offset",
        )
    ),
    _opt("feature_width", _C, ["--feature_width", "20"], ["--feature_width", "120"]),
    _opt("track_type", _C, ["--track_type", "middle"], ["--track_type", "spreadout"]),
    _opt("circular_track_order", _C,
         ["--circular_track_order", "features,ticks,gc_skew,gc_content"]),
    _opt("circular_track_slot", _C,
         ["--circular_track_slot", "gc_content:dinucleotide_content@w=60px"]),
    _opt("circular_track_table (adds AT skew)", _C,
         ["--circular_track_table", "{IN}/tracks_hmmt.tsv"], cases=["C1", "C2"]),
    _opt("gc_content_width", _C, ["--gc_content_width", "20"], ["--gc_content_width", "90"]),
    _opt("gc_content_radius", _C, ["--gc_content_radius", "0.5"], ["--gc_content_radius", "0.9"],
         ["--gc_content_radius", "0.7"]),
    _opt("gc_skew_width", _C, ["--gc_skew_width", "20"], ["--gc_skew_width", "90"]),
    _opt("gc_skew_radius", _C, ["--gc_skew_radius", "0.4"], ["--gc_skew_radius", "0.7"]),
    _opt("depth_width", _C, ["--depth_width", "20"], ["--depth_width", "80"], cases=["C3"]),
    _opt("conservation_labels", _C, ["--conservation_labels", "Other"], cases=["C3"]),
    _opt("conservation_colors", _C, ["--conservation_colors", "#ff0000"], cases=["C3"]),
    _opt("conservation_ring_width", _C, ["--conservation_ring_width", "15"],
         ["--conservation_ring_width", "60"], cases=["C3"]),
    _opt("conservation_ring_gap", _C, ["--conservation_ring_gap", "2"],
         ["--conservation_ring_gap", "30"], cases=["C3"]),
    _opt("conservation_reference", _C, ["--conservation_reference", "auto"],
         ["--conservation_reference", "subject"], cases=["C3"]),
    _opt("evalue", _C, ["--evalue", "1e-50"], cases=["C3"]),
    _opt("bitscore", _C, ["--bitscore", "500"], cases=["C3"]),
    _opt("identity", _C, ["--identity", "95"], cases=["C3"]),
    _opt("alignment_length", _C, ["--alignment_length", "300"], cases=["C3"]),
    _opt("multi_record_size_mode", _C, ["--multi_record_size_mode", "linear"],
         ["--multi_record_size_mode", "equal"], cases=["C2"]),
    _opt("multi_record_min_radius_ratio", _C, ["--multi_record_min_radius_ratio", "0.3"],
         cases=["C2"]),
    _opt("multi_record_column_gap_ratio", _C, ["--multi_record_column_gap_ratio", "0.5"],
         cases=["C2"]),
    _opt("multi_record_row_gap_ratio", _C, ["--multi_record_row_gap_ratio", "0.4"], cases=["C2"]),
    _opt("multi_record_position", _C, _ROW2, cases=["C2"]),
    # Linear only.
    _opt("align_center", _L, ["--align_center"], cases=["L2"]),
    _opt("keep_definition_left_aligned", _L, ["--keep_definition_left_aligned"]),
    _opt("evalue", _L, ["--evalue", "1e-50"], cases=["L2"]),
    _opt("bitscore", _L, ["--bitscore", "500"], cases=["L2"]),
    _opt("identity", _L, ["--identity", "95"], cases=["L2"]),
    _opt("alignment_length", _L, ["--alignment_length", "300"], cases=["L2"]),
    _opt("pairwise_match_style", _L, ["--pairwise_match_style", "curve"], cases=["L2"]),
    _opt("definition_line_style", _L,
         ["--definition_line_style", "accession:size=20,weight=bold,color=#cc0000"]),
    _opt("record_label", _L, ["--record_label", "Label A", "--record_label", "Label B"]),
    _opt("record_subtitle", _L, ["--record_subtitle", "Sub A", "--record_subtitle", "Sub B"]),
    _opt("show_replicon", _L, ["--show_replicon"]),
    _opt("hide_accession", _L, ["--hide_accession"]),
    _opt("hide_length", _L, ["--hide_length"]),
    _opt("label_placement", _L, ["--label_placement", "above_feature"]),
    _opt("label_rotation", _L, ["--label_rotation", "30"]),
    _opt("linear_label_spacing", _L, ["--linear_label_spacing", "40"]),
    _opt("show_labels", _L, ["--show_labels", "first"], ["--show_labels", "none"], cases=["L2"]),
    _opt("show_labels", _L, ["--show_labels", "none"], cases=["L1"]),
    _opt("track_layout", _L, ["--track_layout", "above"], ["--track_layout", "below"]),
    _opt("track_axis_gap", _L, ["--track_axis_gap", "60"], pre=["--track_layout", "above"]),
    _opt("linear_track_order", _L, ["--linear_track_order", "features,gc_skew,gc_content"]),
    _opt("linear_track_slot", _L,
         ["--linear_track_slot", "gc_content:dinucleotide_content@h=60px"]),
    _opt("linear_track_axis_index", _L, ["--linear_track_axis_index", "0"],
         pre=["--linear_track_order", "features,gc_content,gc_skew"]),
    _opt("ruler_on_axis", _L, ["--ruler_on_axis"], pre=["--track_layout", "above"]),
    _opt("feature_height", _L, ["--feature_height", "20"], ["--feature_height", "120"]),
    _opt("gc_height", _L, ["--gc_height", "10"], ["--gc_height", "80"]),
    _opt("comparison_height", _L, ["--comparison_height", "20"], ["--comparison_height", "150"],
         cases=["L2"]),
    _opt("depth_track_height", _L, ["--depth_track_height", "20"],
         ["--depth_track_height", "90"], cases=["L2"]),
    _opt("depth_height", _L, ["--depth_height", "20"], ["--depth_height", "90"], cases=["L2"]),
    _opt("scale_style", _L, ["--scale_style", "bar"]),
    _opt("scale_stroke_color", _L, ["--scale_stroke_color", "red"]),
    _opt("scale_stroke_width", _L, ["--scale_stroke_width", "1"], ["--scale_stroke_width", "10"]),
    _opt("scale_font_size", _L, ["--scale_font_size", "10"], ["--scale_font_size", "40"]),
    _opt("ruler_label_font_size", _L, ["--ruler_label_font_size", "10"],
         ["--ruler_label_font_size", "50"]),
    _opt("ruler_label_color", _L, ["--ruler_label_color", "red"]),
    _opt("normalize_length", _L, ["--normalize_length"], cases=["L2"]),
    _opt("region", _L, ["--region", "#1:2000-9000"], cases=["L1"]),
    # Comparison rows use the coordinates of the cropped records, so this pair has no BLAST file.
    _opt("region (per record)", _L, ["--region", "#1:1000-100000", "--region", "#2:5000-120000"],
         pre_remove=[("-b", "{E}/MjeNMV.MelaMJNV.tblastx.out")], cases=["L2"]),
    _opt("record_id", _L, ["--record_id", "#1"], cases=["L1"]),
    _opt("reverse_complement", _L, ["--reverse_complement", "true"], cases=["L1"]),
    _opt("reverse_complement (per file)", _L,
         ["--reverse_complement", "false", "--reverse_complement", "true"], cases=["L2"]),
    _opt("multi_record_position", _L, _ROW2, pre=["--comparisons_table", "{IN}/comp_l3.tsv"],
         pre_remove=_ROWRM, cases=["L3"]),
    _opt("linear_record_gap", _L, ["--linear_record_gap", "80"], pre=_ROWPRE, pre_remove=_ROWRM,
         cases=["L3"]),
    _opt("align_center (rows)", _L, ["--align_center"], pre=_ROWPRE, pre_remove=_ROWRM,
         cases=["L3"]),
    _opt("keep_definition_left_aligned (rows)", _L, ["--keep_definition_left_aligned"],
         pre=_ROWPRE, pre_remove=_ROWRM, cases=["L3"]),
    _opt("record_label (rows)", _L,
         ["--record_label", "Label A", "--record_label", "Label B", "--record_label", "Label C"],
         pre=_ROWPRE, pre_remove=_ROWRM, cases=["L3"]),
    _opt("track_layout (rows)", _L, ["--track_layout", "above"], pre=_ROWPRE, pre_remove=_ROWRM,
         cases=["L3"]),
)

# Variants the fast set leaves out (design 3.2: 16-18 s each as CLI runs).
SLOW_VARIANTS: frozenset[tuple[str, str]] = frozenset({
    ("circular_label_spacing", "--circular_label_spacing 40"),
    ("label_font_size", "--label_font_size 30"),
})


def options_for(case: Case) -> tuple[OptionSpec, ...]:
    """The options of the universe that apply to one case."""

    return tuple(
        option
        for option in UNIVERSE
        if case.mode in option.modes
        and (option.cases is None or case.case_id in option.cases)
    )


def _remove_tokens(arguments: list[str], removals: Sequence[Sequence[str]]) -> list[str]:
    for tokens in removals:
        size = len(tokens)
        for index in range(len(arguments) - size + 1):
            if tuple(arguments[index:index + size]) == tuple(tokens):
                del arguments[index:index + size]
                break
    return arguments


@dataclass(frozen=True)
class Placeholders:
    """Directories the case arguments name."""

    tables: Path

    def expand(self, case: Case, tokens: Sequence[str]) -> list[str]:
        values = {"T": str(TEST_INPUTS), "E": str(EXAMPLES), "IN": str(self.tables)}
        resolved = {
            name: value.format(**values) for name, value in case.tables.items()
        }
        return [token.format(**values, **resolved) if "{" in token else token for token in tokens]


def state_arguments(
    case: Case,
    selection: Mapping[OptionSpec, Variant],
    placeholders: Placeholders,
) -> list[str]:
    """The command line of one state: the case base with each selected variant applied."""

    arguments = list(case.base)
    for option, variant in selection.items():
        arguments = _remove_tokens(arguments, (*option.pre_remove, *variant.remove))
    for option, variant in selection.items():
        arguments.extend(option.pre)
        arguments.extend(variant.add)
    return placeholders.expand(case, arguments)


def base_arguments(case: Case, option: OptionSpec, placeholders: Placeholders) -> list[str]:
    """The base of one option: the case base with the option's prerequisites."""

    arguments = _remove_tokens(list(case.base), option.pre_remove)
    arguments.extend(option.pre)
    return placeholders.expand(case, arguments)


# --- Chains and walks ---------------------------------------------------------------------------


@dataclass(frozen=True)
class StepResult:
    """One state of a chain or a walk: warm and cold outcomes, or why it was skipped."""

    label: str
    arguments: tuple[str, ...]
    warm: RenderOutcome | None
    cold: RenderOutcome | None
    refused: str | None = None


def _try_capture(harness: EquivalenceHarness, arguments: Sequence[str]) -> str | None:
    try:
        harness.capture(arguments)
    except SystemExit as exc:
        return f"argument error (exit {exc.code})"
    except Exception as exc:  # noqa: BLE001 - a refused command line is skipped
        return f"{type(exc).__name__}: {exc}"
    return None


def run_step(harness: EquivalenceHarness, label: str, arguments: Sequence[str]) -> StepResult:
    """Render one state warm (on the shared cache) and cold (on a new cache)."""

    refused = _try_capture(harness, arguments)
    if refused is not None:
        return StepResult(label, tuple(arguments), None, None, refused)
    if harness.empty:
        # The first render of a chain is a cold render: it is not repeated.
        cold = harness.warm(arguments)
        return StepResult(label, tuple(arguments), cold, cold)
    return StepResult(label, tuple(arguments), harness.warm(arguments), harness.cold(arguments))


def matrix_variants(
    case: Case,
    *,
    option_ids: Sequence[str] | None = None,
    first_variant_only: bool = False,
    include_slow: bool = True,
) -> list[tuple[OptionSpec, Variant]]:
    """The variants of the selected options of a case, without ``LABEL_FIT_FAILURES``.

    They come in universe order, or in the order of ``option_ids``: the cache
    keeps only what the previous request used, so a chain names the step a
    later step reuses right before it.
    """

    options = options_for(case)
    if option_ids is not None:
        options = tuple(option for option_id in option_ids for option in options if option.option_id == option_id)
    selected: list[tuple[OptionSpec, Variant]] = []
    for option in options:
        variants = [
            variant
            for variant in option.variants
            if (case.case_id, option.option_id, variant.label) not in LABEL_FIT_FAILURES
            and (include_slow or (option.option_id, variant.label) not in SLOW_VARIANTS)
        ]
        selected.extend((option, variant) for variant in variants[:1 if first_variant_only else None])
    return selected


def matrix_steps(
    case: Case,
    placeholders: Placeholders,
    variants: Sequence[tuple[OptionSpec, Variant]],
) -> list[tuple[str, list[str]]]:
    """The chain of a case: its base, then each variant (a single-option change of the base)."""

    return [
        (f"{case.case_id} base", placeholders.expand(case, case.base)),
        *(
            (f"{case.case_id} {option.key} {variant.label}", state_arguments(case, {option: variant}, placeholders))
            for option, variant in variants
        ),
    ]


def run_chain(
    harness: EquivalenceHarness,
    steps: Sequence[tuple[str, Sequence[str]]],
) -> list[StepResult]:
    return [run_step(harness, label, arguments) for label, arguments in steps]


def random_walk(
    harness: EquivalenceHarness,
    case: Case,
    placeholders: Placeholders,
    *,
    seed: int,
    steps: int,
) -> list[StepResult]:
    """Seeded single-option changes from the case base, with returns to earlier states.

    Each step either sets one option to a variant, resets one option to the
    base, or undoes the previous step (a return to the state before it). A
    walk leaves out the slow and heavy variants (``SLOW_VARIANTS``,
    ``HEAVY_RENDER_SECONDS``, ``LABEL_FIT_FAILURES``): a heavy setting stays in the state and makes
    every later step heavy. It also leaves out the variants that end in a
    user-input error (``WALK_INPUT_ERRORS``), which would fail every later
    step. The matrix chunks cover them.
    """

    def walkable(option: OptionSpec, variant: Variant) -> bool:
        key = (option.option_id, variant.label)
        full = (case.case_id, *key)
        left_out = HEAVY_RENDER_SECONDS.keys() | LABEL_FIT_FAILURES | WALK_INPUT_ERRORS
        return key not in SLOW_VARIANTS and full not in left_out

    rng = np.random.default_rng(seed)
    options = [option for option in options_for(case) if any(walkable(option, v) for v in option.variants)]
    order = {option: index for index, option in enumerate(UNIVERSE)}
    state: dict[OptionSpec, Variant] = {}
    history: list[dict[OptionSpec, Variant]] = []
    results = [run_step(harness, f"{case.case_id} seed {seed} start", placeholders.expand(case, case.base))]
    attempts = 0
    while len(results) <= steps and attempts < steps * 5:
        attempts += 1
        draw = float(rng.random())
        if history and draw < 0.2:
            candidate = history[-1]
            description = "return to the previous state"
        else:
            option = options[int(rng.integers(len(options)))]
            candidate = dict(state)
            if option in state and draw < 0.45:
                del candidate[option]
                description = f"reset {option.key}"
            else:
                choices = [
                    variant
                    for variant in option.variants
                    if variant != state.get(option) and walkable(option, variant)
                ]
                if not choices:
                    continue
                candidate[option] = choices[int(rng.integers(len(choices)))]
                description = f"set {option.key} {candidate[option].label}"
        ordered = dict(sorted(candidate.items(), key=lambda item: order[item[0]]))
        result = run_step(
            harness,
            f"{case.case_id} seed {seed} step {len(results)}: {description}",
            state_arguments(case, ordered, placeholders),
        )
        if result.refused is not None:
            continue
        results.append(result)
        if description == "return to the previous state":
            history.pop()
        else:
            history.append(state)
        state = ordered
    return results


# --- Slow chunks -------------------------------------------------------------------------------

# One cold render per case, in seconds, measured 2026-10-14 on the 32-core
# development machine (load average 2.5-6; a single render ran within 15% of
# that at a load average of 8-10; C2 on its ``--resolve_overlaps`` base, six
# renders at once, T/logs/pr4/c2-variants-2a74.log). A step renders twice:
# warm and cold in a chain, or the variant and the base again in the
# expected-builds test. L2 and L3 doubled after the label fit (#1053) slowed
# Linear renders: L2-full-1of4 took 100 s at a load average of 5 against an
# estimate of 118 s, and timed out (300 s) at a load average of 6-25
# (T/logs/pr4/chunk-L2-1of4-incr-pr4.log).
RENDER_SECONDS: Mapping[str, float] = {"C1": 0.5, "C2": 6.0, "C3": 1.3, "L1": 0.2, "L2": 2.6, "L3": 2.8}
# Cold renders far above their case: the slower label settings. A heavy step
# counts three cold renders; the third is a margin for machine load, which can
# slow a long render by 25-70 % (a warm render is no slower than a cold one).
HEAVY_RENDER_SECONDS: Mapping[tuple[str, str, str], float] = {
    ("C1", "label_font_size", "--label_font_size 30"): 21.0,
    ("C2", "circular_track_order", "--circular_track_order features,ticks,gc_skew,gc_content"): 33.0,
    ("C2", "circular_track_slot", "--circular_track_slot gc_content:dinucleotide_content@w=60px"): 25.0,
    ("C2", "circular_track_table (adds AT skew)", "--circular_track_table {IN}/tracks_hmmt.tsv"): 23.0,
    # C3: half a warm-and-cold step under pytest at a load average of 1.
    ("C3", "labels", "--labels both"): 5.0,
    ("C3", "outer_label_y_radius_offset", "--outer_label_y_radius_offset 0.8"): 7.0,
    ("C3", "inner_label_x_radius_offset", "--inner_label_x_radius_offset 0.8"): 6.0,
    ("C3", "inner_label_x_radius_offset", "--inner_label_x_radius_offset 1.3"): 7.0,
    ("C3", "inner_label_y_radius_offset", "--inner_label_y_radius_offset 0.8"): 5.0,
    ("C3", "inner_label_y_radius_offset", "--inner_label_y_radius_offset 1.3"): 7.0,
}
# C2 variants left out of the universe. Each ends in LabelsDoNotFitError, so the
# test only checks that warm and cold raise the same error; cold 40-252 s each
# (T/logs/pr4/c2-variants-2a74.log). The cause is the slow failing label fit,
# bug OV-701, which dev has too.
LABEL_FIT_FAILURES: frozenset[tuple[str, str, str]] = frozenset({
    ("C2", "label_font_size", "--label_font_size 30"),
    ("C2", "label_rendering", "--label_rendering external_only"),
    ("C2", "circular_label_spacing", "--circular_label_spacing 40"),
    ("C2", "outer_label_x_radius_offset", "--outer_label_x_radius_offset 0.8"),
    ("C2", "outer_label_y_radius_offset", "--outer_label_y_radius_offset 0.8"),
    ("C2", "inner_label_y_radius_offset", "--inner_label_y_radius_offset 1.3"),
})
# C2 variants left out of the random walks. Each ends at once in a user-input
# error (ValidationError, _InsideFitError) on the C2 base; in a walk it stays in
# the state, so every later step fails with it and the walk tests nothing
# (C2 seed 24). The matrix chain keeps them (the error rule). Measured at the
# PR 4 review fixes: T/logs/pr4/fix8-c2-walk-pool.log.
WALK_INPUT_ERRORS: frozenset[tuple[str, str, str]] = frozenset({
    ("C2", "features", "-k CDS,rRNA,tRNA,tmRNA,ncRNA,misc_RNA,repeat_region,gene"),
    ("C2", "definition_font_size", "--definition_font_size 40"),
    ("C2", "separate_strands", "--separate_strands"),
    ("C2", "keep_full_definition_with_plot_title", "--keep_full_definition_with_plot_title"),
    ("C2", "center_reserved_radius", "--center_reserved_radius 400"),
    ("C2", "center_reserved_radius", "--center_reserved_radius 250"),
    ("C2", "feature_width", "--feature_width 120"),
    ("C2", "gc_content_radius", "--gc_content_radius 0.9"),
    ("C2", "gc_content_radius", "--gc_content_radius 0.7"),
    ("C2", "gc_skew_radius", "--gc_skew_radius 0.7"),
})
# The estimated time of one chunk.
CHUNK_BUDGET_SECONDS = 120.0
# The universe tests are slow and run only with ``--render-stage-universe``
# (tests/conftest.py); ``slow`` also keeps them out of the release job.
UNIVERSE_MARKS = (pytest.mark.slow, pytest.mark.render_stage_universe)


@dataclass(frozen=True)
class Chunk:
    """Variants ``start:stop`` of a case, rendered as one test that starts from the case base."""

    chunk_id: str
    case_id: str
    start: int
    stop: int
    seconds: float


def cost_chunks(
    case_id: str,
    prefix: str,
    variants: Sequence[tuple[OptionSpec, Variant]],
) -> list[Chunk]:
    """Contiguous chunks of ``variants`` whose estimated time stays within the budget.

    A chunk costs one base render plus two renders per variant (three for a
    heavy one; ``RENDER_SECONDS``, ``HEAVY_RENDER_SECONDS``). A variant above
    the budget closes its chunk; when that leaves a chunk of heavy variants
    alone, it takes a light neighbour, so each chunk holds a step that can
    render (heavy C2 variants end in a user-input error; ``error_rule_failures``).
    """

    base = RENDER_SECONDS[case_id]
    heavy = [(case_id, option.option_id, variant.label) in HEAVY_RENDER_SECONDS for option, variant in variants]
    costs = [
        3 * HEAVY_RENDER_SECONDS[(case_id, option.option_id, variant.label)] if is_heavy else 2 * base
        for (option, variant), is_heavy in zip(variants, heavy)
    ]
    starts = [0]
    total = base
    for index, cost in enumerate(costs):
        if index > starts[-1] and total + cost > CHUNK_BUDGET_SECONDS:
            starts.append(index)
            total = base
        total += cost
    stops = [*starts[1:], len(variants)]
    for number, (start, stop) in enumerate(zip(starts, stops)):
        if not all(heavy[start:stop]):
            continue
        if number > 0 and start - starts[number - 1] >= 2 and not heavy[start - 1]:
            starts[number] = stops[number - 1] = start - 1
        elif number + 1 < len(starts) and stops[number + 1] - stop >= 2 and not heavy[stop]:
            stops[number] = starts[number + 1] = stop + 1
    return [
        Chunk(f"{prefix}{number + 1}of{len(starts)}", case_id, start, stop, base + sum(costs[start:stop]))
        for number, (start, stop) in enumerate(zip(starts, stops))
    ]


def full_chain_chunks(case_ids: Sequence[str]) -> dict[str, Chunk]:
    """The chunks of the full chain of each case, by chunk id ("C2-full-3of20")."""

    return {
        chunk.chunk_id: chunk
        for case_id in case_ids
        for chunk in cost_chunks(case_id, f"{case_id}-full-", matrix_variants(CASES[case_id]))
    }


def error_rule_failures(
    base: tuple[str, RenderOutcome | None],
    steps: Sequence[tuple[str, RenderOutcome]],
) -> list[str]:
    """The error rule for the renders of one test.

    The base state must render. A step may fail only with a user-input
    validation error (``USER_INPUT_ERRORS``): some variants and walk
    combinations are invalid input on the base branch too. Any other error
    fails. The accepted errors must be fewer than the steps, so a test cannot
    pass with every step failing.
    """

    label, outcome = base
    failures = []
    if outcome is None or outcome.error is not None:
        failures.append(f"{label}: the base state must render, got {outcome and outcome.error!r}")
    accepted = 0
    for label, outcome in steps:
        if outcome.error is None:
            continue
        if outcome.user_input_error:
            accepted += 1
        else:
            failures.append(f"{label}: render failed with {outcome.error!r}")
    if steps and accepted >= len(steps):
        failures.append(f"all {len(steps)} steps failed with a user-input validation error")
    return failures


def _rendered_steps(results: Sequence[StepResult]) -> list[tuple[str, RenderOutcome, RenderOutcome]]:
    """The steps after the base that were not refused: label, warm and cold outcomes."""

    return [
        (result.label, result.warm, result.cold)
        for result in results[1:]
        if result.warm is not None and result.cold is not None
    ]


def equivalence_failures(results: Sequence[StepResult]) -> list[str]:
    """Every step whose warm outcome differs from its cold outcome, and the error rule.

    ``results[0]`` is the base state of the chain or walk. Errors are compared
    by type and message; an error equal on both sides passes only under
    ``error_rule_failures``.
    """

    failures = []
    steps = _rendered_steps(results)
    for label, warm, cold in steps:
        try:
            assert_same_outcome(warm, cold, label)
        except AssertionError as error:
            failures.append(str(error))
    base = (results[0].label, results[0].cold)
    return failures + error_rule_failures(base, [(label, cold) for label, _, cold in steps])


def declaration_failures(results: Sequence[StepResult]) -> list[str]:
    """Every render with an undeclared read or an uncacheable stage, and the error rule."""

    base = (results[0].label, results[0].cold)
    outcomes = [
        (f"{label} ({name})", outcome)
        for label, warm, cold in _rendered_steps(results)
        for name, outcome in (("warm", warm), ("cold", cold))
    ]
    failures = []
    for label, outcome in (base, *outcomes):
        if outcome is None or outcome.error is not None:
            continue
        counts = outcome.stages
        if counts.undeclared_read_total or counts.uncacheable_total:
            failures.append(
                f"{label}: undeclared reads "
                f"{counts.undeclared_read_total} {dict(counts.undeclared_reads)}, "
                f"uncacheable {counts.uncacheable_total} {dict(counts.uncacheable)}"
            )
    return failures + error_rule_failures(base, outcomes)


# --- Synthetic inputs for the identity cases ------------------------------------------------------


def synthetic_record(
    record_id: str,
    *,
    length: int = 3000,
    seed: int = 0,
    gc: float = 0.5,
    genes: Sequence[tuple[str, int]] = (("dnaA", 120), ("gyrB", 900), ("recA", 1800)),
) -> SeqRecord:
    """A small linear DNA record with a few CDS and one tRNA; ``gc`` sets its composition."""

    rng = np.random.default_rng(seed)
    weights = [(1 - gc) / 2, gc / 2, gc / 2, (1 - gc) / 2]
    sequence = "".join(rng.choice(["A", "C", "G", "T"], size=length, p=weights))
    record = SeqRecord(Seq(sequence), id=record_id, name=record_id, description=f"{record_id} synthetic")
    record.annotations.update(topology="linear", molecule_type="DNA")
    features = [
        SeqFeature(
            SimpleLocation(start, start + 300, strand=1 if index % 2 == 0 else -1),
            type="CDS",
            qualifiers={
                "gene": [gene],
                "product": [f"{gene} protein"],
                "locus_tag": [f"{record_id}_{index:04d}"],
                "protein_id": [f"{record_id}-{gene}"],
                "translation": ["M" * 99],
            },
        )
        for index, (gene, start) in enumerate(genes)
    ]
    features.append(
        SeqFeature(
            SimpleLocation(length - 200, length - 120, strand=1),
            type="tRNA",
            qualifiers={"product": ["tRNA-Leu"], "locus_tag": [f"{record_id}_t001"]},
        )
    )
    record.features = features
    return record


def write_genbank(path: Path, records: Sequence[SeqRecord]) -> Path:
    SeqIO.write(list(records), path, "genbank")
    return path


def write_gff_fasta(directory: Path, record: SeqRecord) -> tuple[Path, Path]:
    """The same record as GFF3 plus FASTA (genes, CDS and the tRNA)."""

    gff = directory / f"{record.id}.gff3"
    fasta = directory / f"{record.id}.fasta"
    lines = ["##gff-version 3", f"##sequence-region {record.id} 1 {len(record.seq)}"]
    for index, feature in enumerate(record.features):
        start = int(feature.location.start) + 1
        end = int(feature.location.end)
        strand = "+" if feature.location.strand == 1 else "-"
        tag = feature.qualifiers["locus_tag"][0]
        product = feature.qualifiers["product"][0]
        lines.append(f"{record.id}\ttest\tgene\t{start}\t{end}\t.\t{strand}\t.\tID=gene{index};locus_tag={tag}")
        lines.append(
            f"{record.id}\ttest\t{feature.type}\t{start}\t{end}\t.\t{strand}\t0\t"
            f"ID=f{index};Parent=gene{index};product={product};locus_tag={tag}"
        )
    gff.write_text("\n".join(lines) + "\n", encoding="utf-8")
    fasta.write_text(f">{record.id}\n{record.seq}\n", encoding="utf-8")
    return gff, fasta
