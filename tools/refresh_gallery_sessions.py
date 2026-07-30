#!/usr/bin/env python3
from __future__ import annotations

import argparse
import copy
import gc
import gzip
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
from collections.abc import Mapping
from contextlib import contextmanager
from pathlib import Path
from typing import Any
from xml.etree import ElementTree as ET


REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from gbdraw.session_io import (  # noqa: E402
    CURRENT_SESSION_VERSION,
    LOSAT_DERIVED_CACHE_SCHEMA,
    PROTEIN_IDENTITY_MANIFEST_SCHEMA,
    classify_raw_losat_cache_entry,
    load_session,
    normalize_current_session_artifacts,
    session_mode,
    validate_session,
    write_session_json,
)
from gbdraw.session_request_codec import CANONICAL_REQUEST_SCHEMA  # noqa: E402
from gbdraw.render.formats import INTERACTIVE_SVG_FORMAT  # noqa: E402
from gbdraw.tracks.circular import (  # noqa: E402
    CircularTrackSlot,
    normalize_circular_track_slots_with_axis,
    parse_circular_track_slots,
)
from gbdraw.tracks.scalars import ScalarSpec  # noqa: E402


GALLERY_ROOT = REPO_ROOT / "gbdraw" / "web" / "gallery"
SESSION_ROOT = GALLERY_ROOT / "sessions"
SESSION_PROMOTER = REPO_ROOT / "tools" / "promote_gallery_session.mjs"
VIBRIO_SESSION_NAME = "vibrio-harveyi-group-collinear.gbdraw-session.json.gz"
VIBRIO_RAW_ENTRY_COUNT = 59
VIBRIO_GZIP_HARD_LIMIT = 90_000_000
VIBRIO_EXPANDED_HARD_LIMIT = 400_000_000
VIBRIO_GZIP_REGRESSION_CEILING = 95_000_000
VIBRIO_EXPANDED_REGRESSION_CEILING = 420_000_000
GEOMETRY_TOLERANCE_PX = 1e-6
LINEAR_SHARED_SPACING_RENDERERS = frozenset(
    {"features", "dinucleotide_content", "dinucleotide_skew"}
)
REFRESHED_GALLERY_ARTIFACT_KEYS = (
    "version",
    "createdAt",
    "results",
    "features",
    "editorState",
    "orthogroupState",
    "runMetadata",
    "losatCache",
    "losatDerivedCache",
    "proteinIdentityManifest",
    "legacyArtifacts",
)
SESSION_ENVELOPE_KEYS = (
    "format",
    "version",
    "createdAt",
    "renderRequest",
    "resources",
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


_VIBRIO_RECORD_KEYS = tuple(f"record-{index}" for index in range(1, 12))
_VIBRIO_ROW_GROUPS = (
    _VIBRIO_RECORD_KEYS[0:3],
    _VIBRIO_RECORD_KEYS[3:5],
    _VIBRIO_RECORD_KEYS[5:7],
    _VIBRIO_RECORD_KEYS[7:9],
    _VIBRIO_RECORD_KEYS[9:11],
)
VIBRIO_EXPECTED_RAW_PAIRS = frozenset(
    {(record_key, record_key) for record_key in _VIBRIO_RECORD_KEYS}
    | {
        pair
        for index in range(len(_VIBRIO_RECORD_KEYS) - 1)
        for pair in (
            (_VIBRIO_RECORD_KEYS[index], _VIBRIO_RECORD_KEYS[index + 1]),
            (_VIBRIO_RECORD_KEYS[index + 1], _VIBRIO_RECORD_KEYS[index]),
        )
    }
    | {
        pair
        for index in range(len(_VIBRIO_ROW_GROUPS) - 1)
        for pair in _directed_cross_pairs(
            _VIBRIO_ROW_GROUPS[index],
            _VIBRIO_ROW_GROUPS[index + 1],
        )
    }
)

GALLERY_SESSION_FILES = (
    "BGC0000708-BGC0000713.gbdraw-session.json",
    "HmmtDNA_basic_circular.gbdraw-session.json",
    "HmmtDNA_ATskew.gbdraw-session.json",
    "tobacco-chloroplast.gbdraw-session.json",
    "Vnig_TUMSAT-TG-2018.gbdraw-session.json.gz",
    "WSSV_genome_comparison.gbdraw-session.json",
    "hepatoplasmataceae_collinear.gbdraw-session.json.gz",
    "vibrio-harveyi-group-collinear.gbdraw-session.json.gz",
    "hepatoplasmataceae_orthogroup.gbdraw-session.json.gz",
    "majanivirus_orthogroup.gbdraw-session.json.gz",
    "lambda_basic_linear.gbdraw-session.json",
)


def _validate_gallery_session_inventory() -> None:
    """Keep the physical, refresh-tool, and Gallery example inventories aligned."""

    from tools.prepare_interactive_gallery_assets import EXAMPLES

    configured = set(GALLERY_SESSION_FILES)
    if len(configured) != len(GALLERY_SESSION_FILES):
        raise ValueError("GALLERY_SESSION_FILES contains duplicate session names")
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
    session_names: tuple[str, ...],
    *,
    include_assets: bool,
) -> tuple[Path, ...]:
    targets = {_session_path(name) for name in session_names}
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


def _session_cli_invocation(session: Mapping[str, Any]) -> Mapping[str, Any] | None:
    cli = session.get("cliInvocation")
    if isinstance(cli, Mapping) and isinstance(cli.get("args"), list):
        return cli
    config = session.get("config")
    if isinstance(config, Mapping):
        cli = config.get("cliInvocation")
        if isinstance(cli, Mapping) and isinstance(cli.get("args"), list):
            return cli
    return None


def _with_interactive_svg_format(args: list[Any]) -> list[str]:
    updated: list[str] = []
    index = 0
    found_format = False
    while index < len(args):
        token = str(args[index])
        if token in {"-f", "--format"}:
            updated.append(token)
            updated.append("interactive_svg")
            index += 2 if index + 1 < len(args) else 1
            found_format = True
            continue
        if token.startswith("--format="):
            updated.append("--format=interactive_svg")
            index += 1
            found_format = True
            continue
        updated.append(token)
        index += 1
    if not found_format:
        updated.extend(["-f", "interactive_svg"])
    return updated


def _preserve_gallery_cli_invocation(
    source_session: Mapping[str, Any],
    refreshed_session: dict[str, Any],
    *,
    mode: str,
) -> bool:
    source_cli = _session_cli_invocation(source_session)
    if source_cli is None:
        return False

    preserved_cli = copy.deepcopy(dict(source_cli))
    preserved_cli["schema"] = 1
    preserved_cli["mode"] = mode
    preserved_cli["args"] = _with_interactive_svg_format(list(source_cli["args"]))
    preserved_cli["renderFormats"] = [INTERACTIVE_SVG_FORMAT]
    preserved_cli.setdefault("fileBindings", [])
    preserved_cli.setdefault("generatedBy", "gbdraw")
    refreshed_session["cliInvocation"] = preserved_cli
    return True


def _promote_gallery_session(
    source_path: Path,
    output_path: Path,
    *,
    env: Mapping[str, str],
) -> None:
    if source_path.is_file():
        source_session = load_session(source_path)
        if (
            source_session.get("renderRequest", {}).get("schema")
            == CANONICAL_REQUEST_SCHEMA
        ):
            write_session_json(output_path, source_session)
            return
    node = shutil.which("node")
    if node is None:
        raise RuntimeError(
            "Gallery session promotion requires Node.js because the browser and "
            "refresh command share the canonical session projection code."
        )
    if not SESSION_PROMOTER.is_file():
        raise FileNotFoundError(
            f"Missing gallery session promoter: {SESSION_PROMOTER.relative_to(REPO_ROOT)}"
        )
    subprocess.run(
        [
            node,
            str(SESSION_PROMOTER),
            str(source_path),
            str(output_path),
        ],
        cwd=REPO_ROOT,
        env=dict(env),
        check=True,
    )
    load_session(output_path)


def _merge_refreshed_gallery_artifacts(
    promoted_session: Mapping[str, Any],
    refreshed_session: Mapping[str, Any],
) -> dict[str, Any]:
    """Merge fresh canonical output and artifacts with promoted gallery metadata."""

    merged = dict(promoted_session)
    refreshed_request = refreshed_session.get("renderRequest")
    if (
        refreshed_session.get("version") == CURRENT_SESSION_VERSION
        and isinstance(refreshed_request, Mapping)
        and refreshed_request.get("schema") == CANONICAL_REQUEST_SCHEMA
    ):
        merged["renderRequest"] = refreshed_request
    for key in REFRESHED_GALLERY_ARTIFACT_KEYS:
        if key in refreshed_session:
            merged[key] = refreshed_session[key]
        elif key in {
            "editorState",
            "runMetadata",
            "losatCache",
            "losatDerivedCache",
            "proteinIdentityManifest",
            "legacyArtifacts",
        }:
            merged.pop(key, None)

    refreshed_resources = refreshed_session.get("resources")
    promoted_resources = promoted_session.get("resources")
    merged["resources"] = {
        **(dict(refreshed_resources) if isinstance(refreshed_resources, Mapping) else {}),
        **(dict(promoted_resources) if isinstance(promoted_resources, Mapping) else {}),
    }
    envelope = {
        key: merged.pop(key)
        for key in SESSION_ENVELOPE_KEYS
        if key in merged
    }
    envelope.update(merged)
    return envelope


def _omit_regenerable_gallery_derived_cache(
    session_path: Path,
    session: dict[str, Any],
) -> None:
    """Keep the oversized Vibrio artifact reconstructible from its raw hits."""

    if session_path.name == VIBRIO_SESSION_NAME:
        session["losatDerivedCache"] = {"entries": []}


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
        if retained_derived_entries:
            raise ValueError(
                f"{session_path.name} must rebuild its derived LOSATP cache "
                "from the retained raw cache"
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

    def referenced_resource_ids(value: object):
        if isinstance(value, Mapping):
            resource_id = value.get("resourceId")
            if isinstance(resource_id, str) and resource_id:
                yield resource_id
            for nested in value.values():
                yield from referenced_resource_ids(nested)
        elif isinstance(value, list):
            for nested in value:
                yield from referenced_resource_ids(nested)

    for resource_id in set(referenced_resource_ids(request)):
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


def _refresh_one_session(
    session_path: Path,
    *,
    destination_path: Path | None = None,
) -> None:
    session = load_session(session_path)
    if session.get("version") == CURRENT_SESSION_VERSION:
        normalize_current_session_artifacts(session)
    mode = session_mode(session)
    if mode not in {"circular", "linear"}:
        raise RuntimeError(f"Could not determine gallery session mode: {session_path}")
    env = os.environ.copy()
    env["PYTHONPATH"] = (
        str(REPO_ROOT)
        if not env.get("PYTHONPATH")
        else f"{REPO_ROOT}{os.pathsep}{env['PYTHONPATH']}"
    )
    source_cli = _session_cli_invocation(session)
    cli_source_session = (
        {"cliInvocation": copy.deepcopy(dict(source_cli))}
        if source_cli is not None
        else {}
    )
    render_request = session.get("renderRequest")
    is_canonical = (
        isinstance(render_request, Mapping)
        and render_request.get("schema") == CANONICAL_REQUEST_SCHEMA
    )
    with tempfile.TemporaryDirectory(prefix="gbdraw-gallery-session-") as tmpdir:
        tmpdir_path = Path(tmpdir)
        source_session = tmpdir_path / f"source-{session_path.name}"
        write_session_json(source_session, session)
        if is_canonical:
            render_session = source_session
            promoted_payload = session
        else:
            render_session = tmpdir_path / "promoted.gbdraw-session.json"
            del session
            gc.collect()
            _promote_gallery_session(source_session, render_session, env=env)
            promoted_payload = load_session(render_session)
        for key in REFRESHED_GALLERY_ARTIFACT_KEYS:
            promoted_payload.pop(key, None)
        if is_canonical:
            del session
        gc.collect()
        refreshed_session = tmpdir_path / session_path.name
        subprocess.run(
            [
                sys.executable,
                "-m",
                "gbdraw.cli",
                mode,
                "--session",
                str(render_session),
                "-f",
                "interactive_svg",
                "-o",
                "out",
                "--session_output",
                str(refreshed_session),
            ],
            cwd=tmpdir_path,
            env=env,
            check=True,
        )
        rendered_payload = load_session(refreshed_session)
        refreshed_payload = _merge_refreshed_gallery_artifacts(
            promoted_payload,
            rendered_payload,
        )
        _omit_regenerable_gallery_derived_cache(
            session_path,
            refreshed_payload,
        )
        del promoted_payload, rendered_payload
        gc.collect()
        _preserve_gallery_cli_invocation(
            cli_source_session,
            refreshed_payload,
            mode=mode,
        )
        write_session_json(refreshed_session, refreshed_payload)
        del refreshed_payload
        gc.collect()
        shutil.move(str(refreshed_session), destination_path or session_path)


def refresh_gallery_sessions(
    session_names: tuple[str, ...] = GALLERY_SESSION_FILES,
) -> None:
    if session_names == GALLERY_SESSION_FILES:
        _validate_gallery_session_inventory()
    session_paths = [_session_path(session_name) for session_name in session_names]
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
        for session_path in session_paths:
            print(f"Refreshing gallery session: {session_path.relative_to(REPO_ROOT)}")
            staged_path = staging_root / session_path.name
            _refresh_one_session(session_path, destination_path=staged_path)
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


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Regenerate web gallery session JSON files with the current CLI/session schema. "
            "By default, gallery SVG sources, thumbnails, and examples.json are refreshed afterwards."
        )
    )
    parser.add_argument(
        "--session",
        action="append",
        default=None,
        help=(
            "Refresh one gallery session by file name or gallery id. "
            "May be repeated; defaults to all gallery sessions."
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
    args = parser.parse_args(argv)

    session_names = tuple(args.session or GALLERY_SESSION_FILES)
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
