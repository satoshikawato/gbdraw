#!/usr/bin/env python3
from __future__ import annotations

import argparse
import base64
import copy
import csv
import dataclasses
import gc
import gzip
import hashlib
import io
import json
import math
import os
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
from collections.abc import Callable, Mapping
from contextlib import contextmanager
from pathlib import Path, PurePosixPath, PureWindowsPath
from typing import TYPE_CHECKING, Any
from xml.etree import ElementTree as ET


REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from gbdraw.session_io import (  # noqa: E402
    CURRENT_SESSION_VERSION,
    DRAWINGS_SESSION_MIN_VERSION,
    LOSAT_DERIVED_CACHE_SCHEMA,
    PROTEIN_IDENTITY_MANIFEST_SCHEMA,
    SessionFileBinding,
    canonicalize_cli_invocation,
    classify_raw_losat_cache_entry,
    load_session,
    session_mode,
    validate_session,
    write_session_json,
)
from gbdraw.session import _keep_resources  # noqa: E402
from gbdraw.session_drawings import (  # noqa: E402
    _drawing_view,
    replace_session_drawings,
    write_session_drawings,
)
from gbdraw.session_contract import GALLERY_SESSION_DRAWN_DRAWINGS  # noqa: E402
from gbdraw.session_resources import drawing_resource_ids  # noqa: E402
from gbdraw.session_request_codec import CANONICAL_REQUEST_SCHEMA  # noqa: E402
from gbdraw.tracks.circular import (  # noqa: E402
    CircularTrackSlot,
    normalize_circular_track_slots_with_axis,
    parse_circular_track_slots,
)
from gbdraw.tracks.scalars import ScalarSpec  # noqa: E402

if TYPE_CHECKING:
    from tools.prepare_interactive_gallery_assets import GallerySessionExample


GALLERY_ROOT = REPO_ROOT / "gbdraw" / "web" / "gallery"
SESSION_ROOT = GALLERY_ROOT / "sessions"
TEST_INPUT_SESSION_ROOT = REPO_ROOT / "tests" / "test_inputs"
SESSION_PUBLICATION_BRIDGE = REPO_ROOT / "tools" / "publish_gallery_session.mjs"
VIBRIO_SESSION_NAME = "vibrio-harveyi-group-collinear.gbdraw-session.json.gz"
VIBRIO_RAW_ENTRY_COUNT = 12
VIBRIO_GZIP_HARD_LIMIT = 90_000_000
VIBRIO_EXPANDED_HARD_LIMIT = 400_000_000
VIBRIO_GZIP_REGRESSION_CEILING = 95_000_000
VIBRIO_EXPANDED_REGRESSION_CEILING = 420_000_000
GEOMETRY_TOLERANCE_PX = 1e-6
LINEAR_SHARED_SPACING_RENDERERS = frozenset(
    {"features", "dinucleotide_content", "dinucleotide_skew"}
)


def _directed_cross_pairs(
    left: tuple[str, ...],
    right: tuple[str, ...],
) -> set[tuple[str, str]]:
    return {
        pair
        for query in left
        for subject in right
        for pair in ((query, subject), (subject, query))
    }


_VIBRIO_RECORD_KEYS = tuple(f"record-{index}" for index in range(1, 5))
_VIBRIO_ROW_GROUPS = (
    _VIBRIO_RECORD_KEYS[0:2],
    _VIBRIO_RECORD_KEYS[2:4],
)
VIBRIO_EXPECTED_RAW_PAIRS = frozenset(
    {(record_key, record_key) for record_key in _VIBRIO_RECORD_KEYS}
    | {
        pair
        for index in range(len(_VIBRIO_ROW_GROUPS) - 1)
        for pair in _directed_cross_pairs(
            _VIBRIO_ROW_GROUPS[index],
            _VIBRIO_ROW_GROUPS[index + 1],
        )
    }
)

TEST_INPUT_SESSION_FILES = (
    "BGC0000708-BGC0000713.gbdraw-session.json",
)


def _is_test_input_session(session_path: Path) -> bool:
    """A test input keeps the unresolved pipelines its specs run on Generate."""

    return (
        session_path.parent == TEST_INPUT_SESSION_ROOT
        and session_path.name in TEST_INPUT_SESSION_FILES
    )


def _public_gallery_session_files() -> tuple[str, ...]:
    from tools.prepare_interactive_gallery_assets import EXAMPLES

    return tuple(example.session_path.name for example in EXAMPLES)


def _validate_gallery_session_inventory() -> None:
    """Keep the physical, refresh-tool, and Gallery example inventories aligned."""

    from tools.prepare_interactive_gallery_assets import EXAMPLES

    public_files = _public_gallery_session_files()
    configured = set(public_files)
    if len(configured) != len(public_files):
        raise ValueError("EXAMPLES contains duplicate session names")
    physical = {
        path.name
        for pattern in ("*.gbdraw-session.json", "*.gbdraw-session.json.gz")
        for path in SESSION_ROOT.glob(pattern)
    }
    examples = {example.session_path.name for example in EXAMPLES}
    if configured != physical or configured != examples:
        raise ValueError(
            "Gallery session inventory mismatch: "
            f"missing-physical={sorted(configured - physical)}, "
            f"unconfigured-physical={sorted(physical - configured)}, "
            f"missing-examples={sorted(configured - examples)}, "
            f"unconfigured-examples={sorted(examples - configured)}"
        )


def _gallery_mutation_targets(
    session_names: tuple[str, ...] | None,
    *,
    include_assets: bool,
) -> tuple[Path, ...]:
    targets = set(_refresh_session_paths(session_names))
    if include_assets:
        from tools.prepare_interactive_gallery_assets import (
            EXAMPLES,
            EXAMPLE_ROOT,
            GALLERY_ROOT as ASSET_GALLERY_ROOT,
            SOURCE_ROOT,
            THUMBNAIL_ROOT,
        )

        for example in EXAMPLES:
            targets.update(
                {
                    example.session_path,
                    example.source_svg_path,
                    example.gallery_svg_path,
                    example.thumbnail_path,
                }
            )
        targets.add(ASSET_GALLERY_ROOT / "examples.json")
        targets.add(ASSET_GALLERY_ROOT / "artifact-manifest.json")
        for root, pattern in (
            (EXAMPLE_ROOT, "*.svg"),
            (SOURCE_ROOT, "*.svg"),
            (THUMBNAIL_ROOT, "*.webp"),
        ):
            if root.exists():
                targets.update(root.glob(pattern))
    return tuple(sorted(targets, key=lambda path: str(path)))


@contextmanager
def _gallery_file_transaction(paths: tuple[Path, ...]):
    """Restore every Gallery output if a refresh phase raises an exception."""

    with tempfile.TemporaryDirectory(prefix="gbdraw-gallery-backup-") as tmpdir:
        backup_root = Path(tmpdir)
        snapshots: list[tuple[Path, Path | None]] = []
        for index, path in enumerate(paths):
            if path.is_file():
                backup_path = backup_root / f"{index:04d}.bak"
                shutil.copy2(path, backup_path)
                snapshots.append((path, backup_path))
            else:
                snapshots.append((path, None))
        try:
            yield
        except BaseException:
            for path, backup_path in snapshots:
                if backup_path is None:
                    if path.is_file():
                        path.unlink()
                    continue
                path.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(backup_path, path)
            raise


def _session_path(name_or_id: str) -> Path:
    name = name_or_id
    if name.endswith((".gbdraw-session.json", ".gbdraw-session.json.gz")):
        return SESSION_ROOT / name
    compressed_path = SESSION_ROOT / f"{name}.gbdraw-session.json.gz"
    if compressed_path.exists():
        return compressed_path
    return SESSION_ROOT / f"{name}.gbdraw-session.json"


def _refresh_session_paths(session_names: tuple[str, ...] | None) -> tuple[Path, ...]:
    if session_names is not None:
        return tuple(_session_path(name) for name in session_names)
    return (
        *(SESSION_ROOT / name for name in _public_gallery_session_files()),
        *(TEST_INPUT_SESSION_ROOT / name for name in TEST_INPUT_SESSION_FILES),
    )


def _run_publication_bridge(
    command: str,
    *paths: Path,
    env: Mapping[str, str],
) -> None:
    node = shutil.which("node")
    if node is None:
        raise RuntimeError(
            "Gallery publication requires Node.js for the shared Web request contract."
        )
    if not SESSION_PUBLICATION_BRIDGE.is_file():
        raise FileNotFoundError(
            "Missing Gallery publication bridge: "
            f"{SESSION_PUBLICATION_BRIDGE.relative_to(REPO_ROOT)}"
        )
    subprocess.run(
        [node, str(SESSION_PUBLICATION_BRIDGE), command, *(str(path) for path in paths)],
        cwd=REPO_ROOT,
        env=dict(env),
        check=True,
    )


def _geometry_number(value: object) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError("Resolved track geometry has an invalid number") from exc
    if not math.isfinite(number):
        raise ValueError("Resolved track geometry has a non-finite number")
    return number


def _geometry_band(slot: Mapping[str, Any]) -> tuple[float, float]:
    band = slot.get("reserveBand")
    if not isinstance(band, Mapping):
        raise ValueError("Resolved track geometry has no reserveBand")
    top = _geometry_number(band.get("topPx"))
    bottom = _geometry_number(band.get("bottomPx"))
    if bottom + GEOMETRY_TOLERANCE_PX < top:
        raise ValueError("Resolved track geometry has an inverted reserveBand")
    return top, bottom


def _geometry_records(geometry: Mapping[str, Any]) -> list[Mapping[str, Any]]:
    records = geometry.get("records")
    if (
        not isinstance(records, list)
        or not records
        or any(not isinstance(record, Mapping) for record in records)
    ):
        raise ValueError("Resolved track geometry has no valid records")
    return records


def _request_tracks(request: Mapping[str, Any]) -> Mapping[str, Any]:
    options = request.get("diagramOptions")
    tracks = options.get("tracks") if isinstance(options, Mapping) else None
    return tracks if isinstance(tracks, Mapping) else {}


def _request_default_linear_track_spacing(
    request: Mapping[str, Any],
) -> float | None:
    options = request.get("diagramOptions")
    config = options.get("config") if isinstance(options, Mapping) else None
    canvas = config.get("canvas") if isinstance(config, Mapping) else None
    linear = canvas.get("linear") if isinstance(canvas, Mapping) else None
    if not isinstance(linear, Mapping):
        return None
    if linear.get("track_spacing") is not None:
        return max(0.0, _geometry_number(linear["track_spacing"]))
    return 0.0


def _linear_spacing_after(
    slot: Mapping[str, Any],
    *,
    default_track_spacing: float | None,
) -> float:
    actual = max(0.0, _geometry_number(slot.get("spacingAfterPx")))
    if (
        default_track_spacing is not None
        and slot.get("renderer") in LINEAR_SHARED_SPACING_RENDERERS
    ):
        if not math.isclose(
            actual,
            default_track_spacing,
            rel_tol=1e-9,
            abs_tol=GEOMETRY_TOLERANCE_PX,
        ):
            raise ValueError(
                "Resolved Linear default spacing metadata is inconsistent "
                f"at slot {slot.get('slotId')!r}"
            )
        return default_track_spacing
    return actual


def _validate_circular_feature_geometry(
    request: Mapping[str, Any],
    geometry: Mapping[str, Any],
) -> None:
    tracks = _request_tracks(request)
    specs = tracks.get("circularTrackSlots")
    if specs is None:
        return
    decoded_specs = [
        _canonical_circular_track_slot(spec) if isinstance(spec, Mapping) else spec
        for spec in specs
    ]
    expected_slots = normalize_circular_track_slots_with_axis(
        parse_circular_track_slots(decoded_specs),
        tracks.get("circularTrackAxisIndex"),
    )
    expected_lanes = {
        str(slot.id): str(slot.params.get("lane_direction") or "")
        for slot in expected_slots
        if slot.enabled and slot.renderer == "features"
    }

    for record in _geometry_records(geometry):
        axis_radius = _geometry_number(record.get("axisRadiusPx"))
        if axis_radius <= 0:
            raise ValueError("Resolved circular track geometry has a nonpositive Axis")
        slots = record.get("slots")
        if not isinstance(slots, list):
            raise ValueError("Resolved circular track geometry has no slots")
        actual_features = {
            str(slot.get("slotId")): slot
            for slot in slots
            if isinstance(slot, Mapping) and slot.get("renderer") == "features"
        }
        for slot_id, lane in expected_lanes.items():
            actual = actual_features.get(slot_id)
            if actual is None:
                raise ValueError(
                    f"Resolved circular Feature slot {slot_id!r} is missing"
                )

            center = (
                _geometry_number(actual.get("radiusFactor"))
                * axis_radius
            )
            half_width = 0.5 * _geometry_number(actual.get("widthPx"))
            inner = center - half_width
            outer = center + half_width
            valid_band = {
                "split": (
                    inner < axis_radius - GEOMETRY_TOLERANCE_PX
                    and outer > axis_radius + GEOMETRY_TOLERANCE_PX
                ),
                "inside": outer <= axis_radius + GEOMETRY_TOLERANCE_PX,
                "outside": inner >= axis_radius - GEOMETRY_TOLERANCE_PX,
            }.get(lane, False)
            if not valid_band:
                raise ValueError(
                    f"Resolved circular Feature lane geometry for {slot_id!r} "
                    f"is inconsistent with lane_direction={lane}"
                )


def _canonical_circular_track_slot(value: Mapping[str, Any]) -> CircularTrackSlot:
    """Decode the current canonical slot shape for geometry validation."""

    def scalar(raw: object) -> ScalarSpec | None:
        if raw is None:
            return None
        if not isinstance(raw, Mapping):
            raise ValueError("Canonical Circular track scalar must be an object")
        unit = raw["unit"]
        if unit not in {"px", "factor"}:
            raise ValueError(f"Unsupported circular track scalar unit: {unit}")
        return ScalarSpec(value=float(raw["value"]), unit=unit)

    if value.get("kind") != "circularTrackSlot":
        raise ValueError("Canonical Circular track slot has an invalid kind")
    params = value.get("params")
    if not isinstance(params, Mapping):
        raise ValueError("Canonical Circular track slot params must be an object")
    return CircularTrackSlot(
        id=str(value.get("id") or ""),
        renderer=str(value.get("renderer") or ""),
        enabled=bool(value.get("enabled")),
        side=value.get("side"),
        radius=scalar(value.get("radius")),
        width=scalar(value.get("width")),
        z=int(value.get("z", 0)),
        params=dict(params),
        inner_gap_px=value.get("innerGapPx"),
        outer_gap_px=value.get("outerGapPx"),
    )


def _validate_linear_side_adjacency(
    slots: list[Mapping[str, Any]],
    *,
    direction: str,
    feature_band: tuple[float, float] | None,
    feature_spacing: float,
    exact_spacing: bool,
    default_track_spacing: float | None,
) -> None:
    occupied_band = feature_band
    incoming_spacing = max(0.0, feature_spacing)
    for slot in sorted(
        (slot for slot in slots if slot.get("side") == direction),
        key=lambda slot: int(slot.get("slotIndex", -1)),
        reverse=direction == "above",
    ):
        reserve_band = _geometry_band(slot)
        if (
            slot.get("paintBand") is None
            and reserve_band[1] - reserve_band[0] <= GEOMETRY_TOLERANCE_PX
        ):
            continue
        if occupied_band is None:
            occupied_band = reserve_band
            incoming_spacing = _linear_spacing_after(
                slot,
                default_track_spacing=default_track_spacing,
            )
            continue
        gap = (
            occupied_band[0] - reserve_band[1]
            if direction == "above"
            else reserve_band[0] - occupied_band[1]
        )
        if gap < incoming_spacing - GEOMETRY_TOLERANCE_PX:
            raise ValueError(
                "Resolved Linear adjacent reserve bands overlap or violate "
                f"declared spacing at slot {slot.get('slotId')!r}"
            )
        if (
            exact_spacing
            and not math.isclose(
                gap,
                incoming_spacing,
                rel_tol=1e-9,
                abs_tol=GEOMETRY_TOLERANCE_PX,
            )
        ):
            raise ValueError(
                "Resolved Linear default adjacent reserve gap "
                f"is inconsistent at slot {slot.get('slotId')!r}"
            )
        occupied_band = reserve_band
        incoming_spacing = _linear_spacing_after(
            slot,
            default_track_spacing=default_track_spacing,
        )


def _validate_linear_reserve_geometry(
    request: Mapping[str, Any],
    geometry: Mapping[str, Any],
) -> None:
    default_layout = _request_tracks(request).get("linearTrackSlots") is None
    default_track_spacing = (
        _request_default_linear_track_spacing(request)
        if default_layout
        else None
    )
    for record in _geometry_records(geometry):
        slots_value = record.get("slots")
        if not isinstance(slots_value, list):
            raise ValueError("Resolved Linear track geometry has no slots")
        slots = [slot for slot in slots_value if isinstance(slot, Mapping)]
        structural = [
            slot
            for slot in slots
            if slot.get("side") == "overlay"
            and slot.get("renderer") == "features"
            and slot.get("dataAvailable") is not False
            and slot.get("paintBand") is not None
        ]
        if structural:
            feature_bands = [_geometry_band(slot) for slot in structural]
            feature_band = (
                min(top for top, _bottom in feature_bands),
                max(bottom for _top, bottom in feature_bands),
            )
            feature_spacing = max(
                _linear_spacing_after(
                    slot,
                    default_track_spacing=default_track_spacing,
                )
                for slot in structural
            )
        else:
            feature_band = None
            feature_spacing = 0.0
        for direction in ("above", "below"):
            _validate_linear_side_adjacency(
                slots,
                direction=direction,
                feature_band=feature_band,
                feature_spacing=feature_spacing,
                exact_spacing=default_layout,
                default_track_spacing=default_track_spacing,
            )


def _serialized_component_bytes(value: object) -> int:
    if value is None:
        return 0
    encoder = json.JSONEncoder(
        ensure_ascii=False,
        separators=(",", ":"),
    )
    return sum(len(chunk.encode("utf-8")) for chunk in encoder.iterencode(value))


def _session_artifact_measurements(
    session: Mapping[str, Any],
    artifact_path: Path,
) -> dict[str, int]:
    compressed_bytes = artifact_path.stat().st_size
    with artifact_path.open("rb") as session_file:
        is_gzip = session_file.read(2) == b"\x1f\x8b"
    opener = gzip.open if is_gzip else Path.open
    with opener(artifact_path, "rb") as session_file:
        expanded_bytes = sum(
            len(chunk)
            for chunk in iter(lambda: session_file.read(1024 * 1024), b"")
        )

    editor_state = session.get("editorState")
    editor_without_catalog = (
        dict(editor_state) if isinstance(editor_state, Mapping) else {}
    )
    feature_catalog = editor_without_catalog.pop("featureCatalog", None)
    return {
        "sessionGzipBytes": compressed_bytes,
        "sessionExpandedBytes": expanded_bytes,
        "resultsBytes": _serialized_component_bytes(session.get("results")),
        "resourcesBytes": _serialized_component_bytes(
            session.get("resources")
        ),
        "webFilesBytes": _serialized_component_bytes(session.get("webFiles")),
        "featureCatalogBytes": _serialized_component_bytes(feature_catalog),
        "editorStateBytes": _serialized_component_bytes(
            editor_without_catalog
        ),
        "losatCacheBytes": _serialized_component_bytes(
            session.get("losatCache")
        ),
        "losatDerivedCacheBytes": _serialized_component_bytes(
            session.get("losatDerivedCache")
        ),
    }


def _format_session_measurements(measurements: Mapping[str, int]) -> str:
    return ", ".join(
        f"{key}={value}" for key, value in measurements.items()
    )


def _validate_current_session_catalog_structure(
    session_path: Path,
    session: Mapping[str, Any],
) -> None:
    resources = session.get("resources")
    seen_payloads: dict[str, str] = {}
    if isinstance(resources, Mapping):
        for resource_id, resource in resources.items():
            if not isinstance(resource, Mapping):
                continue
            data = resource.get("data")
            if not isinstance(data, str) or not data:
                continue
            encoding = str(resource.get("encoding") or "text")
            hasher = hashlib.sha256()
            hasher.update(encoding.encode("utf-8"))
            hasher.update(b"\0")
            for offset in range(0, len(data), 1024 * 1024):
                hasher.update(
                    data[offset : offset + 1024 * 1024].encode("utf-8")
                )
            digest = hasher.hexdigest()
            previous = seen_payloads.get(digest)
            if previous is not None:
                raise ValueError(
                    f"{session_path.name} duplicates one resource payload as "
                    f"{previous} and {resource_id}"
                )
            seen_payloads[digest] = str(resource_id)

    editor_state = session.get("editorState")
    catalog = (
        editor_state.get("featureCatalog")
        if isinstance(editor_state, Mapping)
        else None
    )
    if not isinstance(catalog, Mapping):
        return
    items = catalog.get("items")
    if not isinstance(items, list):
        return
    rendered_allowed = {
        "svgId",
        "recordKey",
        "biologicalFeatureId",
        "fillColor",
        "drawnSelector",
    }
    duplicated_payload_keys = {
        "qualifiers",
        "nucleotide_sequence",
        "nucleotideSequence",
        "amino_acid_sequence",
        "aminoAcidSequence",
    }

    def duplicated_keys(value: object) -> set[str]:
        if isinstance(value, Mapping):
            found = duplicated_payload_keys & set(value)
            for nested in value.values():
                found.update(duplicated_keys(nested))
            return found
        if isinstance(value, list):
            found: set[str] = set()
            for nested in value:
                found.update(duplicated_keys(nested))
            return found
        return set()

    for item in items:
        if not isinstance(item, Mapping):
            continue
        biological = item.get("biologicalFeatures")
        biological = biological if isinstance(biological, list) else []
        biological_keys = [
            (
                str(feature.get("recordKey") or ""),
                str(feature.get("biologicalFeatureId") or ""),
            )
            for feature in biological
            if isinstance(feature, Mapping)
        ]
        if len(biological_keys) != len(set(biological_keys)):
            raise ValueError(
                f"{session_path.name} duplicates a biological feature payload"
            )
        rendered = item.get("features")
        rendered = rendered if isinstance(rendered, list) else []
        for feature in rendered:
            if not isinstance(feature, Mapping):
                continue
            unexpected = set(feature) - rendered_allowed
            if unexpected:
                raise ValueError(
                    f"{session_path.name} rendered feature duplicates "
                    f"biological fields: {', '.join(sorted(unexpected))}"
                )
        duplicated = set()
        for section in ("features", "orthogroups", "comparisonMatches"):
            duplicated.update(duplicated_keys(item.get(section)))
        if duplicated:
            raise ValueError(
                f"{session_path.name} copies biological qualifiers/sequences "
                f"outside biologicalFeatures: {', '.join(sorted(duplicated))}"
            )


def _referenced_resource_ids(value: object):
    """Yield explicit resource references from a canonical request."""

    if isinstance(value, Mapping):
        resource_id = value.get("resourceId")
        if isinstance(resource_id, str) and resource_id:
            yield resource_id
        for nested in value.values():
            yield from _referenced_resource_ids(nested)
    elif isinstance(value, list):
        for nested in value:
            yield from _referenced_resource_ids(nested)


def _validate_staged_gallery_session(
    session_path: Path,
    session: dict[str, Any],
    *,
    artifact_path: Path | None = None,
    require_track_geometry: bool = False,
) -> None:
    from tools.prepare_interactive_gallery_assets import (
        EXAMPLES,
        _validate_session_interactive_orthogroups,
    )

    validate_session(session)
    # A staged Session is a Session 49 that shows one drawing: read its drawn
    # drawing as the one-drawing view (this checks the Gallery rule).
    session, _drawing_id = _gallery_view(session)
    expected_session_version = CURRENT_SESSION_VERSION
    expected_request_schema = CANONICAL_REQUEST_SCHEMA
    if session.get("version") != expected_session_version:
        raise ValueError(
            f"{session_path.name} has session version {session.get('version')}; "
            f"expected {expected_session_version}"
        )
    request = session.get("renderRequest")
    if (
        not isinstance(request, Mapping)
        or request.get("schema") != expected_request_schema
    ):
        raise ValueError(
            f"{session_path.name} has no canonical schema-"
            f"{expected_request_schema} render request"
        )
    # A Gallery Session stores its protein comparison results, so the first
    # Generate after Load reuses them (perf-014x D-17 / docs-restructure D-20:
    # the regenerated Vibrio Session must store the finished protein comparison
    # results; replaces the empty-derived-cache-only rule of 7aab08f8). A test
    # input keeps its pipeline, which the Generate-pipeline specs run.
    if not _is_test_input_session(session_path) and any(
        isinstance(comparison, Mapping)
        and comparison.get("kind") == "generatedProteinComparison"
        and comparison.get("mode") != "none"
        for comparison in request.get("comparisons") or []
    ):
        raise ValueError(
            f"{session_path.name} must store its resolved protein comparison "
            "results instead of a pipeline that Generate runs again"
        )
    if (
        session_path.name == "lambda_basic_linear.gbdraw-session.json"
        and request.get("comparisons")
    ):
        raise ValueError(
            "lambda_basic_linear must contain zero comparison descriptors"
        )
    run_metadata = session.get("runMetadata")
    geometry = (
        run_metadata.get("trackSlotGeometry")
        if isinstance(run_metadata, Mapping)
        else None
    )
    if require_track_geometry and not isinstance(geometry, Mapping):
        raise ValueError(
            f"{session_path.name} has no resolved track-slot run metadata"
        )
    if isinstance(geometry, Mapping):
        if (
            geometry.get("schema") != 1
            or geometry.get("source") != "resolved"
            or geometry.get("mode") != request.get("mode")
        ):
            raise ValueError(
                f"{session_path.name} has invalid resolved track-slot run metadata"
            )
        if request["mode"] == "circular":
            _validate_circular_feature_geometry(request, geometry)
        elif request["mode"] == "linear":
            _validate_linear_reserve_geometry(request, geometry)
    resources = session.get("resources")
    if not isinstance(resources, Mapping):
        raise ValueError(f"{session_path.name} has no canonical resources")
    _validate_current_session_catalog_structure(session_path, session)

    manifest = session.get("proteinIdentityManifest")
    if (
        not isinstance(manifest, Mapping)
        or manifest.get("schema") != PROTEIN_IDENTITY_MANIFEST_SCHEMA
    ):
        raise ValueError(
            f"{session_path.name} has no protein identity manifest schema "
            f"{PROTEIN_IDENTITY_MANIFEST_SCHEMA}"
        )
    raw_cache = session.get("losatCache")
    raw_entries = raw_cache.get("entries", []) if isinstance(raw_cache, Mapping) else []
    if not isinstance(raw_entries, list) or any(
        classify_raw_losat_cache_entry(entry)
        not in {"protein-current", "nucleotide-current"}
        for entry in raw_entries
    ):
        raise ValueError(f"{session_path.name} contains a non-current raw LOSAT artifact")
    if session_path.name == VIBRIO_SESSION_NAME:
        if artifact_path is None:
            raise ValueError("Vibrio size validation requires the staged session path")
        derived_cache = session.get("losatDerivedCache")
        retained_derived_entries = (
            derived_cache.get("entries", [])
            if isinstance(derived_cache, Mapping)
            else []
        )
        # The stored collinear result replaces derived entries; the raw cache
        # lets a changed setting recompute without LOSATP.
        if retained_derived_entries:
            raise ValueError(
                f"{session_path.name} stores its collinear result in the request, "
                "not as derived LOSATP entries"
            )
        protein_entries = [
            entry
            for entry in raw_entries
            if classify_raw_losat_cache_entry(entry) == "protein-current"
        ]
        actual_pairs = {
            (
                str(entry.get("queryRecordInstanceKey") or ""),
                str(entry.get("subjectRecordInstanceKey") or ""),
            )
            for entry in protein_entries
        }
        if (
            len(protein_entries) != VIBRIO_RAW_ENTRY_COUNT
            or actual_pairs != VIBRIO_EXPECTED_RAW_PAIRS
        ):
            raise ValueError(
                f"{session_path.name} must retain the exact "
                f"{VIBRIO_RAW_ENTRY_COUNT}-pair protein cache"
            )
        measurements = _session_artifact_measurements(
            session,
            artifact_path,
        )
        print(
            f"Gallery session metrics [{session_path.name}]: "
            f"{_format_session_measurements(measurements)}"
        )
        if measurements["sessionGzipBytes"] > VIBRIO_GZIP_HARD_LIMIT:
            raise ValueError(
                f"{session_path.name} gzip size exceeds the "
                f"{VIBRIO_GZIP_HARD_LIMIT}-byte hard limit "
                f"({_format_session_measurements(measurements)})"
            )
        if (
            measurements["sessionExpandedBytes"]
            > VIBRIO_EXPANDED_HARD_LIMIT
        ):
            raise ValueError(
                f"{session_path.name} expanded size exceeds the "
                f"{VIBRIO_EXPANDED_HARD_LIMIT}-byte hard limit "
                f"({_format_session_measurements(measurements)})"
            )
    derived_cache = session.get("losatDerivedCache")
    derived_entries = (
        derived_cache.get("entries", [])
        if isinstance(derived_cache, Mapping)
        else []
    )
    if not isinstance(derived_entries, list) or any(
        not isinstance(entry, Mapping)
        or entry.get("schema") != LOSAT_DERIVED_CACHE_SCHEMA
        or entry.get("idEncoding") != "runtime-handle-v1"
        for entry in derived_entries
    ):
        raise ValueError(
            f"{session_path.name} contains a non-current derived LOSATP artifact"
        )
    legacy_artifacts = session.get("legacyArtifacts")
    candidates = (
        legacy_artifacts.get("proteinRawCandidates")
        if isinstance(legacy_artifacts, Mapping)
        else None
    )
    if isinstance(candidates, Mapping) and candidates.get("entries"):
        raise ValueError(
            f"{session_path.name} retained legacy protein raw candidates after refresh"
        )
    protein_artifacts = {
        key: session.get(key)
        for key in (
            "results",
            "features",
            "editorState",
            "orthogroupState",
            "losatCache",
            "losatDerivedCache",
            "proteinIdentityManifest",
        )
    }
    serialized_protein_artifacts = json.dumps(
        protein_artifacts, ensure_ascii=False
    )
    if "p_r_" in serialized_protein_artifacts:
        raise ValueError(
            f"{session_path.name} contains unresolved legacy protein identifiers"
        )
    if re.search(
        r"@.+\|.+~f_[0-9a-f]{64}",
        serialized_protein_artifacts,
    ):
        raise ValueError(
            f"{session_path.name} contains unsupported long protein transport "
            "identifiers"
        )
    del protein_artifacts, serialized_protein_artifacts

    for resource_id in set(_referenced_resource_ids(request)):
        resource = resources.get(resource_id)
        if not isinstance(resource, Mapping) or not str(resource.get("data") or ""):
            raise ValueError(
                f"{session_path.name} references missing or empty resource {resource_id}"
            )

    results = session.get("results")
    svg_content = next(
        (
            result.get("content")
            for result in results
            if isinstance(result, Mapping)
            and isinstance(result.get("content"), str)
            and "<svg" in result["content"]
        ),
        None,
    ) if isinstance(results, list) else None
    if not isinstance(svg_content, str):
        raise ValueError(f"{session_path.name} has no generated SVG result")
    ET.fromstring(svg_content)

    example = next(
        (item for item in EXAMPLES if item.session_path.name == session_path.name),
        None,
    )
    if example is None:
        return
    _validate_session_interactive_orthogroups(example, session)


def _gallery_id(session_path: Path) -> str:
    """The Gallery id of a Session file: its name without the Session suffix."""

    return session_path.name.split(".gbdraw-session", 1)[0]


def _set_output_prefix(session: dict[str, Any], prefix: str) -> int:
    """Store ``prefix`` in every place a Session keeps its Output Prefix.

    A Generate after Load names its files after this prefix, so the form field
    and the render request must agree with the Gallery id the card command uses.
    Returns the number of fields changed.
    """

    containers: list[object] = []
    # Session 46 keeps the form field in each mode slice (``modes.<mode>``).
    modes = session.get("modes")
    for mode_slice in modes.values() if isinstance(modes, dict) else ():
        config = mode_slice.get("config") if isinstance(mode_slice, dict) else None
        containers.append(config.get("form") if isinstance(config, dict) else None)
    request = session.get("renderRequest")
    containers.append(request.get("output") if isinstance(request, dict) else None)
    changed = 0
    for container in containers:
        if isinstance(container, dict) and container.get("prefix") != prefix:
            container["prefix"] = prefix
            changed += 1
    return changed


_DEFAULT_COLORS_HEADER = "feature_type\tcolor"


def _default_color_rows(text: str) -> dict[str, str]:
    lines = [line.split("\t", 1) for line in text.splitlines() if line.strip()]
    return {key.strip(): color.strip() for key, color in lines if f"{key}\t{color}" != _DEFAULT_COLORS_HEADER}


def _resolve_default_colors(session: dict[str, Any]) -> int:
    """Store the request's Default colors resolved over its palette, as the CLI does.

    The Web writes the rows that differ from the palette; the CLI replay writes
    the palette with those rows applied, and ``finalize`` compares the two
    (OV-375). Returns 1 when the table changed, else 0.
    """

    from gbdraw.io.colors import load_default_colors

    colors = session["renderRequest"]["diagramOptions"].get("colors") or {}
    ref = colors.get("defaultColorsFile") or colors.get("defaultColors")
    if not ref:
        return 0
    resource = session["resources"][ref["resourceId"]]
    text = base64.b64decode(resource["data"], validate=True).decode("utf-8")
    with tempfile.NamedTemporaryFile("w", suffix=".tsv", encoding="utf-8") as handle:
        handle.write(text)
        handle.flush()
        table = load_default_colors(handle.name, colors.get("defaultColorsPalette") or "default")
    resolved = dict(zip(table["feature_type"].astype(str), table["color"].astype(str)))
    if resolved == _default_color_rows(text):
        return 0
    header = f"{_DEFAULT_COLORS_HEADER}\n" if text.startswith(_DEFAULT_COLORS_HEADER) else ""
    payload = (header + "".join(f"{key}\t{color}\n" for key, color in resolved.items())).encode("utf-8")
    resource["data"] = base64.b64encode(payload).decode("ascii")
    resource["size"] = len(payload)
    return 1


def _materialize_one_record_sources(session: dict[str, Any]) -> int:
    """Store a one-record GenBank source of ``all`` as ``exactly_one``, as the CLI does.

    The Web writes ``all`` for a Linear file without selector or region; the CLI
    replay materializes each record it draws as ``exactly_one``, and
    ``finalize`` compares the two (OV-376). It reads the prepared request, which
    the Web reader has promoted: a request schema 1-5 stores no cardinality.
    Returns the number of records changed.
    """

    from io import StringIO

    from Bio import SeqIO

    changed = 0
    for record in session["renderRequest"]["records"]:
        source = record.get("source") or {}
        if record.get("cardinality") != "all" or record.get("selector") or record.get("region") or source.get("kind") != "genbank":
            continue
        payload = base64.b64decode(session["resources"][source["resourceId"]]["data"], validate=True)
        if sum(1 for _ in SeqIO.parse(StringIO(payload.decode("utf-8")), "genbank")) == 1:
            record["cardinality"] = "exactly_one"
            changed += 1
    return changed


def _recorded_input_name(arg: str) -> str:
    """The bare file name of a recorded absolute input path in the repository.

    Gallery Sessions record inputs by file name, as the card commands do. A path
    inside this checkout, or recorded in another checkout of the repository (a
    maintainer's machine) that holds the same file here, becomes its file name.
    Any other absolute path is outside the repository and cannot be published.
    """

    windows_path = re.match(r"^[A-Za-z]:[\\/]", arg) is not None
    pure = PureWindowsPath(arg) if windows_path else PurePosixPath(arg)
    if not pure.is_absolute():
        return arg
    if not windows_path:
        try:
            Path(os.path.realpath(arg)).relative_to(REPO_ROOT)
            return pure.name
        except ValueError:
            pass
    checkout_names = {REPO_ROOT.name, "gbdraw"}
    for index, part in enumerate(pure.parts[:-1]):
        if part not in checkout_names:
            continue
        candidate = REPO_ROOT.joinpath(*pure.parts[index + 1 :])
        try:
            Path(os.path.realpath(candidate)).relative_to(REPO_ROOT)
        except ValueError:
            continue
        if candidate.is_file():
            return pure.name
    raise ValueError(
        f"Recorded input path {arg!r} is outside the repository ({REPO_ROOT}); "
        "Gallery Sessions record inputs of the repository by file name."
    )


def _publishable_recorded_args(
    args: list[str],
    bindings: list[SessionFileBinding],
    *,
    gallery_id: str,
    declared_command: bool,
) -> list[str]:
    """Return ``args`` with the Gallery id as ``-o`` and inputs by file name.

    The declared command of an example is authoritative for its inputs (they
    resolve in the example's directory), so only the stored-Session route
    rewrites input paths.
    """

    input_indexes = set() if declared_command else {binding.argIndex for binding in bindings}
    published: list[str] = []
    for index, arg in enumerate(args):
        if index > 0 and args[index - 1] in {"-o", "--output"}:
            published.append(gallery_id)
        elif arg.startswith("--output="):
            published.append(f"--output={gallery_id}")
        elif index in input_indexes:
            published.append(_recorded_input_name(arg))
        else:
            published.append(arg)
    return published


def _canonicalize_recorded_cli_invocation(
    session_path: Path,
    *,
    mode: str,
    gallery_id: str,
    declared_command: bool = False,
) -> None:
    """Make the recorded ``cliInvocation`` of a staged session publishable.

    The recorded argv is provenance for Gallery commands. Current sessions do not
    replay it, so the refresh moves it to the current flag names here, names the
    output after the Gallery id, and records inputs by file name.
    """

    def canonicalize(session: dict[str, Any]) -> int:
        invocation = session.get("cliInvocation")
        if not isinstance(invocation, dict) or not isinstance(invocation.get("args"), list):
            return 0
        args, bindings = canonicalize_cli_invocation(
            [str(arg) for arg in invocation["args"]],
            [
                SessionFileBinding(**binding)
                for binding in invocation.get("fileBindings") or []
            ],
            mode=mode,
        )
        invocation["args"] = _publishable_recorded_args(
            args,
            bindings,
            gallery_id=gallery_id,
            declared_command=declared_command,
        )
        invocation["fileBindings"] = [dataclasses.asdict(binding) for binding in bindings]
        return 1

    _edit_staged_session(session_path, canonicalize)


def _read_session_file(session_path: Path) -> tuple[dict[str, Any], bool]:
    """Return a Session file's JSON and whether the file is gzip-compressed."""

    payload = session_path.read_bytes()
    compressed = payload[:2] == b"\x1f\x8b"
    return json.loads(gzip.decompress(payload) if compressed else payload), compressed


def _write_session_file(session_path: Path, session: Mapping[str, Any], *, compressed: bool) -> None:
    """Write Session JSON with no time in the gzip header, so equal JSON gives equal bytes."""

    text = json.dumps(session, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    session_path.write_bytes(
        gzip.compress(text, compresslevel=6, mtime=0) if compressed else text
    )


def _gallery_view(document: dict[str, Any]) -> tuple[dict[str, Any], str | None]:
    """The one-drawing view a refresh reads and edits, and the ID of its drawing.

    A Gallery Session shows one drawing: a Session 49's only drawing with a
    committed render (``GALLERY_SESSION_DRAWN_DRAWINGS`` in
    gbdraw/session_contract.py, which ``galleryPublicationView`` reads
    generated). A Session 46-48 is its own view (ID ``None``).
    """

    drawings = document.get("drawings")
    if not isinstance(drawings, list):
        return document, None
    drawn = [drawing for drawing in drawings if drawing.get("renderRequest") is not None]
    if len(drawn) != GALLERY_SESSION_DRAWN_DRAWINGS:
        raise ValueError(
            f"A Gallery Session shows {GALLERY_SESSION_DRAWN_DRAWINGS} drawn drawing; this one has {len(drawn)}."
        )
    return _drawing_view(document, drawn[0]), str(drawn[0]["id"])


def _gallery_document(
    document: dict[str, Any], view: dict[str, Any], drawing_id: str | None
) -> dict[str, Any]:
    """``document`` with the drawing of :func:`_gallery_view` replaced by ``view``.

    A resource that no drawing names any more (a stored comparison an edit
    replaced) is left out, as Session 49 requires.
    """

    if drawing_id is None:
        return view
    others = [_drawing_view(document, drawing) for drawing in document["drawings"] if drawing["id"] != drawing_id]
    _keep_resources(view, set().union(*(drawing_resource_ids(fields) for fields in (view, *others))))
    return replace_session_drawings(document, {drawing_id: view})


def _one_drawing_document(session: Mapping[str, Any]) -> Mapping[str, Any]:
    """A document the Session readers take: a Session 49 drawing's one-drawing
    view (:func:`_gallery_view`) as the Session 49 of that one drawing, with
    the resources it names; any other document as it is."""

    if "drawings" in session or int(session["version"]) < DRAWINGS_SESSION_MIN_VERSION:
        return session
    view = dict(session)
    _keep_resources(view, drawing_resource_ids(view))
    return write_session_drawings([view], version=DRAWINGS_SESSION_MIN_VERSION)


def _edit_staged_session(
    session_path: Path, edit: Callable[[dict[str, Any]], int]
) -> None:
    """Rewrite a staged Session file in its own encoding when ``edit`` changes
    the view of its drawing (:func:`_gallery_view`)."""

    document, compressed = _read_session_file(session_path)
    session, drawing_id = _gallery_view(document)
    if edit(session):
        _write_session_file(
            session_path, _gallery_document(document, session, drawing_id), compressed=compressed
        )


def _keep_created_at(tracked_path: Path, staged_path: Path) -> None:
    """Give a refreshed Session the creation time of the Session it replaces.

    The CLI replay stamps the time it ran; a refresh republishes the same
    Session, so two refreshes of one input give the same bytes.
    """

    created_at = _read_session_file(tracked_path)[0]["createdAt"]

    def edit(session: dict[str, Any]) -> int:
        if session.get("createdAt") == created_at:
            return 0
        session["createdAt"] = created_at
        return 1

    _edit_staged_session(staged_path, edit)


def _with_comparisons(comparisons: list[Any], session: dict[str, Any]) -> int:
    """Give a staged Session's request ``comparisons``; 1 when that changes it."""

    request = session["renderRequest"]
    if request["comparisons"] == comparisons:
        return 0
    request["comparisons"] = comparisons
    return 1


def _cli_env() -> dict[str, str]:
    env = os.environ.copy()
    env["PYTHONPATH"] = (
        str(REPO_ROOT)
        if not env.get("PYTHONPATH")
        else f"{REPO_ROOT}{os.pathsep}{env['PYTHONPATH']}"
    )
    return env


def _declared_command_session(
    command: str,
    session_path: Path,
    *,
    env: Mapping[str, str],
    cwd: Path = REPO_ROOT,
) -> None:
    """Write the CLI Session of a declared Gallery command to ``session_path``.

    The command runs in ``cwd`` (the repository root unless the example names
    the directory of its inputs), so its relative inputs resolve as documented;
    only the output prefix moves next to ``session_path``. The recorded
    invocation keeps the declared arguments.
    """

    argv = shlex.split(command)
    if len(argv) < 3 or argv[0] != "gbdraw" or argv[1] not in {"circular", "linear"}:
        raise ValueError(f"Gallery command must start with 'gbdraw circular|linear': {command}")
    mode, args = argv[1], argv[2:]
    output_indexes = [index + 1 for index, arg in enumerate(args) if arg in {"-o", "--output"}]
    if len(output_indexes) != 1 or output_indexes[0] >= len(args):
        raise ValueError(f"Gallery command must name one output prefix: {command}")
    output_index = output_indexes[0]
    staged_output = str(session_path.parent / args[output_index])
    staged_args = [*args[:output_index], staged_output, *args[output_index + 1 :]]
    subprocess.run(
        [
            sys.executable,
            "-m",
            "gbdraw.cli",
            mode,
            *staged_args,
            "--session_output",
            str(session_path),
        ],
        cwd=cwd,
        env=dict(env),
        check=True,
    )
    def edit(session: dict[str, Any]) -> int:
        invocation = session["cliInvocation"]
        invocation["args"] = [
            args[output_index] if arg == staged_output else arg
            for arg in invocation["args"]
        ]
        return 1

    _edit_staged_session(session_path, edit)


# The request comparisons that hold a computed protein result. Next to a
# pipeline marker of mode "none" they are what the Web wrote after a Generate.
_RESOLVED_PROTEIN_KINDS = frozenset(
    {"precomputedProteinComparison", "orthogroupResult", "collinearityResult"}
)
# The comparisons whose resource is the typed analysis result itself.
_PROTEIN_RESULT_KINDS = frozenset({"orthogroupResult", "collinearityResult"})


def _draft_losat_blastp(session: Mapping[str, Any]) -> dict[str, Any] | None:
    """The protein settings of the Web Linear draft (Session 46 keeps them in the Linear slice)."""

    modes = session.get("modes")
    linear = modes.get("linear") if isinstance(modes, Mapping) else None
    config = linear.get("config") if isinstance(linear, Mapping) else None
    losat = config.get("losat") if isinstance(config, Mapping) else None
    blastp = losat.get("blastp") if isinstance(losat, Mapping) else None
    return blastp if isinstance(blastp, dict) else None


def _declared_protein_mode(session: Mapping[str, Any]) -> str | None:
    """The LOSATP mode of the Session's recorded CLI invocation, if it ran LOSATP."""

    invocation = session.get("cliInvocation")
    if not isinstance(invocation, Mapping) or invocation.get("mode") != "linear":
        return None
    from gbdraw.api.options import LOSATP_MODE_WIRE
    from gbdraw.linear import _get_args

    args = _get_args([str(arg) for arg in invocation.get("args") or ()])
    return LOSATP_MODE_WIRE[args.losatp_mode] if args.losat == "losatp" else None


def _reopen_stored_protein_comparisons(session: dict[str, Any], mode: str) -> int:
    """Turn a request's stored protein results back into the pipeline that computes them.

    The stored comparisons and their resources leave the request, the
    pipeline marker asks for ``mode``, and the Web draft's protein mode
    follows it, so the next resolve computes the results from the raw LOSATP
    cache in the mode the Web also rebuilds in (OV-1031, OV-1048).
    """

    request = session["renderRequest"]
    comparisons = list(request.get("comparisons") or [])
    markers = [
        comparison
        for comparison in comparisons
        if comparison.get("kind") == "generatedProteinComparison"
    ]
    stored = [
        comparison
        for comparison in comparisons
        if comparison.get("kind") in _RESOLVED_PROTEIN_KINDS
    ]
    if len(markers) != 1 or markers[0].get("mode") != "none" or not stored:
        return 0
    request["comparisons"] = [
        comparison
        for comparison in comparisons
        if comparison is not markers[0] and comparison not in stored
    ] + [{**markers[0], "mode": mode}]
    for resource_id in set(_referenced_resource_ids(stored)):
        session["resources"].pop(resource_id, None)
    blastp = _draft_losat_blastp(session)
    if blastp is not None:
        blastp["mode"] = mode
    return 1


def _reopen_for_publication(session_path: Path, session: dict[str, Any]) -> int:
    """Reopen a public Gallery Session's stored protein results in its declared mode.

    The declared mode is the recorded CLI invocation's, else the Web draft's.
    A test input keeps its pipeline, which the Generate-pipeline specs run.
    """

    if _is_test_input_session(session_path):
        return 0
    blastp = _draft_losat_blastp(session) or {}
    mode = _declared_protein_mode(session) or blastp.get("mode")
    if not isinstance(mode, str):
        return 0
    return _reopen_stored_protein_comparisons(session, mode)


def _resolved_protein_comparisons(
    session: Mapping[str, Any],
) -> tuple[list[dict[str, Any]], dict[str, Mapping[str, Any]]] | None:
    """Compute a Session's protein pipeline once from its raw LOSATP cache.

    Returns the codec's encoding of the comparisons every Python Session
    writer stores after that render (``stored_protein_result_request``: the
    form the Web writes after a Generate) with their resource descriptors;
    ``None`` when the request runs no protein pipeline.
    A result that needed a LOSATP search the raw cache does not hold is
    refused: a refresh keeps the Session's search evidence.
    """

    from gbdraw.api.request_render import stored_protein_result_request
    from gbdraw.session import materialize_session, render_session, session_to_request
    from gbdraw.session_request_codec import encode_canonical_request
    from gbdraw.session_resources import SessionResourceTable

    modes = {
        comparison.get("mode")
        for comparison in session["renderRequest"].get("comparisons") or []
        if comparison.get("kind") == "generatedProteinComparison"
    } - {"none"}
    if not modes:
        return None
    if not modes <= {"collinear", "orthogroup"} or len(modes) != 1:
        raise ValueError(
            "Gallery publication stores collinear or similarity-group protein "
            f"results only, not {sorted(modes)}"
        )
    (mode,) = modes
    raw_keys = {
        str(entry.get("key") or "")
        for entry in (session.get("losatCache") or {}).get("entries") or []
    }
    with tempfile.TemporaryDirectory(prefix="gbdraw-gallery-resolve-") as output_dir:
        with materialize_session(_one_drawing_document(session), output_directory=output_dir) as materialized:
            typed = session_to_request(materialized)
            rendered = render_session(materialized)
            metadata = getattr(rendered, "linear_metadata", None)
            searched = {
                str(entry.get("key") or "")
                for entry in getattr(rendered, "losat_cache_entries", ())
            } - raw_keys
            if searched:
                raise ValueError(
                    f"the {mode} result needs {len(searched)} LOSATP search(es) "
                    "that the raw cache does not hold"
                )
            table = SessionResourceTable(session["resources"])
            comparisons = encode_canonical_request(
                stored_protein_result_request(typed, metadata), table=table
            ).payload["comparisons"]
            return comparisons, table.descriptors()


def _store_resolved_protein_comparisons(session_path: Path) -> None:
    """Replace a generated protein pipeline with its stored result.

    A Gallery Session stores its protein comparison results, so the first
    Generate after Load reuses them instead of computing them again (perf-014x
    D-17 / docs-restructure D-20: the regenerated Vibrio Session must store the
    finished protein comparison results; replaces the empty-derived-cache-only
    rule of 7aab08f8). The results are computed from the raw LOSATP cache with
    the current code at every refresh (D-65, OV-1031). The replay and
    ``finalize`` then check that the CLI and the Web rebuild draw and request
    the same.
    """

    document = load_session(session_path)
    session, drawing_id = _gallery_view(document)
    try:
        resolved = _resolved_protein_comparisons(session)
    except ValueError as exc:
        raise ValueError(f"{session_path.name}: {exc}") from exc
    if resolved is None:
        return
    comparisons, descriptors = resolved
    session["renderRequest"]["comparisons"] = comparisons
    for resource_id in _referenced_resource_ids(comparisons):
        session["resources"].setdefault(resource_id, dict(descriptors[resource_id]))
    write_session_json(session_path, _gallery_document(document, session, drawing_id))


def _resource_bytes(session: Mapping[str, Any], resource_id: str) -> bytes:
    resource = session["resources"][resource_id]
    data = str(resource.get("data") or "")
    return base64.b64decode(data) if resource.get("encoding") == "base64" else data.encode("utf-8")


def _canonical_table_rows(payload: bytes) -> list[dict[str, str]]:
    """A canonical TSV's rows by column name, as text: equal tables compare equal however written."""

    text = payload.decode("utf-8")
    rows = list(csv.DictReader(io.StringIO(text), delimiter="\t"))
    return [
        {key: _canonical_cell(value) for key, value in row.items()}
        for row in rows
    ]


def _canonical_cell(value: str | None) -> str:
    """One table cell by value: the Web writes ``true`` where pandas writes ``True``, and numbers vary in spelling."""

    text = "" if value is None else str(value)
    if text.lower() in {"true", "false"}:
        return text.lower()
    try:
        number = float(text)
    except ValueError:
        return text
    return repr(number) if math.isfinite(number) else text


def stored_protein_result_differences(session: dict[str, Any]) -> list[str]:
    """How a Session's stored protein results differ from a rebuild from its raw cache.

    The rebuild runs in the Web draft's protein mode, the mode the Web
    computes in when a setting changes, with the current code. An empty list
    means the stored results are the ones the Web would compute (OV-1031).
    """

    request = session["renderRequest"]
    stored = [
        comparison
        for comparison in request.get("comparisons") or []
        if comparison.get("kind") in _RESOLVED_PROTEIN_KINDS
    ]
    mode = (_draft_losat_blastp(session) or {}).get("mode")
    if not stored or not isinstance(mode, str):
        return []
    rebuild = copy.deepcopy(session)
    if not _reopen_stored_protein_comparisons(rebuild, mode):
        return []
    try:
        resolved = _resolved_protein_comparisons(rebuild)
    except Exception as exc:  # the rebuild itself failing is the difference to report
        return [f"cannot rebuild the {mode} result from the raw cache: {exc}"]
    if resolved is None:
        return []
    rebuilt_comparisons, descriptors = resolved
    rebuilt = [
        comparison
        for comparison in rebuilt_comparisons
        if comparison.get("kind") in _RESOLVED_PROTEIN_KINDS
    ]
    stored_kinds = [comparison["kind"] for comparison in stored]
    rebuilt_kinds = [comparison["kind"] for comparison in rebuilt]
    if stored_kinds != rebuilt_kinds:
        return [f"stores {stored_kinds}; the {mode} rebuild gives {rebuilt_kinds}"]
    rebuilt_session = {"resources": descriptors}
    differences = []
    for index, (old, new) in enumerate(zip(stored, rebuilt)):
        old_bytes = _resource_bytes(session, old["resourceId"])
        new_bytes = _resource_bytes(rebuilt_session, new["resourceId"])
        same = (
            old_bytes == new_bytes
            if old["kind"] in _PROTEIN_RESULT_KINDS
            else _canonical_table_rows(old_bytes) == _canonical_table_rows(new_bytes)
        )
        if not same:
            differences.append(
                f"comparisons[{index}] {old['kind']} differs from the {mode} rebuild"
            )
    return differences


def check_stored_protein_results(session_names: tuple[str, ...] | None = None) -> int:
    """Print every Gallery Session whose stored protein results are stale; return their count."""

    stale = 0
    for session_path in _refresh_session_paths(session_names):
        differences = stored_protein_result_differences(_gallery_view(load_session(session_path))[0])
        for difference in differences:
            print(f"{session_path.name}: {difference}")
        stale += bool(differences)
    return stale


def _assert_declared_figure(declared_path: Path, published_path: Path) -> None:
    """The published Session draws the figure of its declared command.

    Publication writes the Web's configuration overrides in place of the CLI's
    resolved configuration, so the replayed figure must equal the declared one.
    """

    declared, published = (
        [result.get("content") for result in _gallery_view(load_session(path))[0].get("results") or []]
        for path in (declared_path, published_path)
    )
    if not declared or declared != published:
        raise ValueError(
            f"{published_path.name} does not draw the figure of its declared command"
        )


def _refresh_one_session(
    session_path: Path,
    *,
    destination_path: Path | None = None,
    declared_command: bool = False,
) -> None:
    gallery_id = _gallery_id(session_path)
    document = load_session(session_path)
    session, drawing_id = _gallery_view(document)
    mode = session_mode(session)
    if mode not in {"circular", "linear"}:
        raise RuntimeError(f"Could not determine gallery session mode: {session_path}")

    env = _cli_env()
    with tempfile.TemporaryDirectory(prefix="gbdraw-gallery-session-") as tmpdir:
        staging_root = Path(tmpdir)
        prepared_path = staging_root / f"prepared-{session_path.name}"
        replayed_path = staging_root / f"replayed-{session_path.name}"
        finalized_path = staging_root / session_path.name
        publication_input_path = session_path
        # The Output Prefix is part of the request the publication rebuilds and
        # `finalize` keeps, so it is set on the input, before `prepare`.
        if _reopen_for_publication(session_path, session) + _set_output_prefix(
            session, gallery_id
        ) + _resolve_default_colors(session):
            publication_input_path = (
                staging_root / f"canonicalized-{session_path.name}"
            )
            write_session_json(publication_input_path, _gallery_document(document, session, drawing_id))
        del session, document
        gc.collect()
        _run_publication_bridge(
            "prepare",
            publication_input_path,
            prepared_path,
            env=env,
        )
        _canonicalize_recorded_cli_invocation(
            prepared_path,
            mode=mode,
            gallery_id=gallery_id,
            declared_command=declared_command,
        )
        _edit_staged_session(prepared_path, _materialize_one_record_sources)
        if not _is_test_input_session(session_path):
            _store_resolved_protein_comparisons(prepared_path)
        subprocess.run(
            [
                sys.executable,
                "-m",
                "gbdraw.cli",
                mode,
                "--session",
                str(prepared_path),
                "-f",
                "interactive_svg",
                "-o",
                gallery_id,
                "--session_output",
                str(replayed_path),
            ],
            cwd=staging_root,
            env=env,
            check=True,
        )
        if _is_test_input_session(session_path):
            # A CLI re-save stores the protein results it computed (OV-1068);
            # a test input keeps the pipelines its specs run on Generate.
            pipelines = _gallery_view(load_session(prepared_path))[0]["renderRequest"]["comparisons"]
            _edit_staged_session(replayed_path, lambda staged: _with_comparisons(pipelines, staged))
        _run_publication_bridge(
            "finalize",
            prepared_path,
            replayed_path,
            finalized_path,
            env=env,
        )
        finalized = load_session(finalized_path)
        validate_session(finalized)
        del finalized
        gc.collect()
        shutil.move(str(finalized_path), destination_path or session_path)


def _declared_session_example(session_path: Path) -> GallerySessionExample | None:
    """The Gallery example whose declared command builds this Session, if any."""

    from tools.prepare_interactive_gallery_assets import EXAMPLES

    return next(
        (
            example
            for example in EXAMPLES
            if example.session_from_command
            and example.session_path.name == session_path.name
        ),
        None,
    )


def refresh_gallery_sessions(
    session_names: tuple[str, ...] | None = None,
) -> None:
    if session_names is None:
        _validate_gallery_session_inventory()
    session_paths = list(_refresh_session_paths(session_names))
    for session_path in session_paths:
        if not session_path.exists():
            raise FileNotFoundError(
                f"Missing gallery session: {session_path.relative_to(REPO_ROOT)}"
            )

    # Render every requested session successfully before replacing any tracked session.
    # A failed promotion therefore cannot leave a half-refreshed gallery behind.
    with tempfile.TemporaryDirectory(
        prefix=".gbdraw-gallery-refresh-",
        dir=SESSION_ROOT,
    ) as staging_dir:
        staging_root = Path(staging_dir)
        staged_paths: list[tuple[Path, Path]] = []
        for index, session_path in enumerate(session_paths):
            print(f"Refreshing gallery session: {session_path.relative_to(REPO_ROOT)}")
            staged_path = staging_root / f"{index:02d}-{session_path.name}"
            example = _declared_session_example(session_path)
            if example is not None:
                with tempfile.TemporaryDirectory(
                    prefix="gbdraw-gallery-command-"
                ) as command_dir:
                    declared_path = Path(command_dir) / session_path.name
                    _declared_command_session(
                        example.command,
                        declared_path,
                        env=_cli_env(),
                        cwd=example.command_cwd,
                    )
                    _refresh_one_session(
                        declared_path,
                        destination_path=staged_path,
                        declared_command=True,
                    )
                    _assert_declared_figure(declared_path, staged_path)
            else:
                _refresh_one_session(session_path, destination_path=staged_path)
            _keep_created_at(session_path, staged_path)
            staged_session = load_session(staged_path)
            _validate_staged_gallery_session(
                session_path,
                staged_session,
                artifact_path=staged_path,
                require_track_geometry=True,
            )
            del staged_session
            gc.collect()
            staged_paths.append((session_path, staged_path))

        for session_path, staged_path in staged_paths:
            staged_path.replace(session_path)


def prepare_gallery_assets() -> None:
    from tools.prepare_interactive_gallery_assets import (
        prepare_gallery_assets as prepare_assets,
    )

    prepare_assets(refresh_sources=True)
    from tools.gallery_artifact_manifest import write_gallery_artifact_manifest

    write_gallery_artifact_manifest()


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Regenerate bundled session JSON files with the current CLI/session schema. "
            "By default, gallery SVG sources, thumbnails, and examples.json are refreshed afterwards."
        )
    )
    parser.add_argument(
        "--session",
        action="append",
        default=None,
        help=(
            "Refresh one gallery session by file name or gallery id. "
            "May be repeated; defaults to all Gallery and test-input sessions."
        ),
    )
    parser.add_argument(
        "--no-assets",
        action="store_true",
        help="Only refresh session JSON files; skip gallery SVG/thumbnail/examples.json preparation.",
    )
    parser.add_argument(
        "--skip-session-refresh",
        action="store_true",
        help="Only prepare gallery SVG/thumbnail/examples.json assets from the existing session files.",
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help=(
            "Write nothing: rebuild each stored protein comparison result from its "
            "Session's raw LOSATP cache with the current code, and exit 1 when one differs."
        ),
    )
    args = parser.parse_args(argv)

    session_names = tuple(args.session) if args.session else None
    if args.check:
        return 1 if check_stored_protein_results(session_names) else 0
    targets = _gallery_mutation_targets(
        () if args.skip_session_refresh else session_names,
        include_assets=not args.no_assets,
    )
    with _gallery_file_transaction(targets):
        if not args.skip_session_refresh:
            refresh_gallery_sessions(session_names)
        if not args.no_assets:
            prepare_gallery_assets()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
