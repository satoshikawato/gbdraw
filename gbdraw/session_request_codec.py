"""Internal codec for canonical, CLI-independent render request payloads.

Canonical request schemas are intentionally strict: unknown fields are rejected instead of being
silently discarded.  File content is represented by stable resource IDs; this
module does not own session envelopes, base64 encoding, or temporary files.
"""

from __future__ import annotations

from collections.abc import Mapping as MappingABC, Sequence as SequenceABC
from copy import deepcopy
from dataclasses import asdict, dataclass, fields, is_dataclass
from io import StringIO
import json
import math
from pathlib import Path
import re
import types
from typing import Any, Literal, Mapping, Sequence, Union, get_args, get_origin, get_type_hints

from Bio import SeqIO  # type: ignore[reportMissingImports]
from pandas import DataFrame, read_csv  # type: ignore[reportMissingImports]

from gbdraw.analysis.collinearity import (  # type: ignore[reportMissingImports]
    CollinearityAnchor,
    CollinearityBlock,
    CollinearityParameters,
    CollinearityResult,
    LosslessCollinearityParameters,
)
from gbdraw.analysis.protein_colinearity import (  # type: ignore[reportMissingImports]
    OrthogroupMember,
    OrthogroupNameCandidate,
    OrthogroupResult,
    OrthologEdge,
    OrthologPath,
)
from gbdraw.config.models import GbdrawConfig  # type: ignore[reportMissingImports]
from gbdraw.exceptions import ValidationError
from gbdraw.io.record_select import RecordSelector
from gbdraw.io.regions import RegionSpec
from gbdraw.io.comparisons import COMPARISON_COLUMNS
from gbdraw.linear_comparison import LinearComparison
from gbdraw.tracks import CircularTrackSlot, LinearTrackSlot, ScalarSpec
from gbdraw.annotations import (
    AnnotationOptions,
    AnnotationSet,
    CoordinateSpan,
    FeatureSelector,
    FeatureSpan,
    HatchStyle,
    RegionAnnotation,
    RegionAnnotationStyle,
)

from .api.options import (
    CircularMultiRecordOptions,
    ColorOptions,
    DiagramOptions,
    LinearMultiRecordOptions,
    OutputOptions,
    TrackOptions,
)
from .api.requests import (
    CircularDiagramRequest,
    DiagramRequest,
    GenBankInputSource,
    GffFastaInputSource,
    InMemoryRecordSource,
    LinearDiagramRequest,
    RecordInput,
    RecordPresentation,
    RenderOutputRequest,
)


CANONICAL_REQUEST_SCHEMA = 3
UNKNOWN_FIELD_POLICY = "reject"

_RESOURCE_ID_RE = re.compile(r"^[a-z][a-z0-9]*(?:-[a-z0-9]+)*$")
_TOP_LEVEL_FIELDS = frozenset(
    {"schema", "mode", "records", "diagramOptions", "layout", "comparisons", "output"}
)
_DEFAULT_OPTIONS = DiagramOptions()

_COMPARISON_SOURCE_FIELDS = frozenset(
    {
        "blast_files",
        "linear_comparisons",
        "protein_comparisons",
        "protein_comparison_pairs",
        "orthogroups",
        "collinearity_blocks",
    }
)
_PIPELINE_FIELDS = (
    "protein_blastp_mode",
    "collinearity_params",
    "collinearity_unit_mode",
    "collinearity_anchor_mode",
    "collinearity_search_scope",
    "collinearity_color_mode",
    "losatp_bin",
    "ncbi_blastp_bin",
    "losatp_threads",
    "protein_blastp_max_hits",
    "protein_blastp_candidate_limit",
    "orthogroup_membership_mode",
    "orthogroup_member_max_hits",
    "collinear_max_paralog_links_per_orthogroup",
    "align_orthogroup_feature",
)
_COMPARISON_FIELDS = _COMPARISON_SOURCE_FIELDS | frozenset(_PIPELINE_FIELDS)

_TABLE_FIELDS = frozenset(
    {
        "feature_table",
        "feature_visibility_table",
        "label_whitelist_table",
        "qualifier_priority_table",
        "label_override_table",
        "depth_table",
    }
)
_FILE_FIELDS = frozenset(
    {
        "feature_table_file",
        "feature_visibility_table_file",
        "label_whitelist_file",
        "qualifier_priority_file",
        "label_override_file",
        "depth_file",
    }
)
_TABLE_SEQUENCE_FIELDS = frozenset({"depth_tables", "conservation_dataframes"})
_FILE_SEQUENCE_FIELDS = frozenset({"depth_files", "conservation_blast_files"})
_TABLE_MATRIX_FIELDS = frozenset({"depth_track_tables"})
_FILE_MATRIX_FIELDS = frozenset({"depth_track_files"})

_TYPED_TREE_CLASSES = {
    cls.__name__: cls
    for cls in (
        CollinearityAnchor,
        CollinearityBlock,
        CollinearityResult,
        OrthogroupMember,
        OrthogroupNameCandidate,
        OrthogroupResult,
        OrthologEdge,
        OrthologPath,
    )
}


class CanonicalRequestCodecError(ValidationError):
    """Base error for canonical request conversion."""


class CanonicalRequestEncodingError(CanonicalRequestCodecError):
    """Raised when a typed request cannot be represented canonically."""


class CanonicalRequestDecodingError(CanonicalRequestCodecError):
    """Raised when a canonical payload cannot produce a typed request."""


@dataclass(frozen=True)
class CanonicalRequestResource:
    """One resource referenced by a canonical request payload."""

    resource_id: str
    kind: str
    name: str
    source_path: Path | None = None
    content: bytes | None = None

    def __post_init__(self) -> None:
        if not _RESOURCE_ID_RE.fullmatch(self.resource_id):
            raise CanonicalRequestEncodingError(
                f"Invalid canonical resource ID: {self.resource_id!r}."
            )
        if not isinstance(self.kind, str) or not self.kind:
            raise CanonicalRequestEncodingError("Canonical resource kind is required.")
        if not isinstance(self.name, str) or not self.name:
            raise CanonicalRequestEncodingError("Canonical resource name is required.")
        if (self.source_path is None) == (self.content is None):
            raise CanonicalRequestEncodingError(
                "A canonical resource must have exactly one path or byte content source."
            )


@dataclass(frozen=True)
class EncodedCanonicalRequest:
    """JSON-compatible request payload plus its out-of-band resources."""

    payload: dict[str, Any]
    resources: tuple[CanonicalRequestResource, ...]


class _ResourceBuilder:
    def __init__(self) -> None:
        self._resources: list[CanonicalRequestResource] = []
        self._ids: set[str] = set()

    def add_path(self, resource_id: str, *, kind: str, value: object) -> str:
        if not isinstance(value, (str, Path)) or not str(value).strip():
            raise CanonicalRequestEncodingError(
                f"Resource {resource_id!r} must identify a materialized file."
            )
        path = Path(str(value))
        if not path.is_file():
            raise CanonicalRequestEncodingError(
                f"Canonical request resource is not a file: {path}."
            )
        self._add(
            CanonicalRequestResource(
                resource_id=resource_id,
                kind=kind,
                name=path.name,
                source_path=path,
            )
        )
        return resource_id

    def add_bytes(
        self,
        resource_id: str,
        *,
        kind: str,
        name: str,
        content: bytes,
    ) -> str:
        self._add(
            CanonicalRequestResource(
                resource_id=resource_id,
                kind=kind,
                name=name,
                content=content,
            )
        )
        return resource_id

    def _add(self, resource: CanonicalRequestResource) -> None:
        if resource.resource_id in self._ids:
            raise CanonicalRequestEncodingError(
                f"Duplicate canonical resource ID: {resource.resource_id}."
            )
        self._ids.add(resource.resource_id)
        self._resources.append(resource)

    def result(self) -> tuple[CanonicalRequestResource, ...]:
        return tuple(self._resources)


def encode_canonical_request(request: DiagramRequest) -> EncodedCanonicalRequest:
    """Encode one typed request and normalize all conversion failures."""

    try:
        return _encode_canonical_request(request)
    except CanonicalRequestCodecError:
        raise
    except (TypeError, ValueError, ValidationError) as exc:
        raise CanonicalRequestEncodingError(
            f"Typed request could not be converted: {exc}"
        ) from exc


def _encode_canonical_request(request: DiagramRequest) -> EncodedCanonicalRequest:
    """Encode without embedding files or changing the session version."""

    if not isinstance(request, (CircularDiagramRequest, LinearDiagramRequest)):
        raise CanonicalRequestEncodingError("Unsupported typed diagram request.")
    _validate_dataclass_contract(request.options, path="diagramOptions", error="encode")
    _validate_dataclass_contract(request.output, path="output", error="encode")
    if request.layout is not None:
        _validate_dataclass_contract(request.layout, path="layout", error="encode")

    resources = _ResourceBuilder()
    mode: Literal["circular", "linear"] = (
        "circular" if isinstance(request, CircularDiagramRequest) else "linear"
    )
    payload = {
        "schema": CANONICAL_REQUEST_SCHEMA,
        "mode": mode,
        "records": [
            _encode_record(record, index=index, resources=resources)
            for index, record in enumerate(request.records, start=1)
        ],
        "diagramOptions": _encode_diagram_options(request.options, resources=resources),
        "layout": _encode_layout(request),
        "comparisons": _encode_comparisons(request.options, mode=mode, resources=resources),
        "output": {
            "prefix": request.output.output_prefix,
            "formats": list(request.output.formats),
            "overwrite": request.output.overwrite,
            "interactiveMetadataPolicy": request.output.interactive_metadata_policy,
        },
    }
    return EncodedCanonicalRequest(payload=payload, resources=resources.result())


def decode_canonical_request(
    payload: Mapping[str, Any],
    *,
    resource_paths: Mapping[str, str | Path],
    output_directory: str | Path,
) -> DiagramRequest:
    """Decode canonical schemas 1, 2, or 3 and normalize validation failures."""

    try:
        return _decode_canonical_request(
            payload,
            resource_paths=resource_paths,
            output_directory=output_directory,
        )
    except CanonicalRequestCodecError:
        raise
    except (TypeError, ValueError, ValidationError) as exc:
        raise CanonicalRequestDecodingError(
            f"Canonical request could not be converted: {exc}"
        ) from exc


def _decode_canonical_request(
    payload: Mapping[str, Any],
    *,
    resource_paths: Mapping[str, str | Path],
    output_directory: str | Path,
) -> DiagramRequest:
    """Decode a supported schema using caller-owned materialized resources."""

    top = _object(payload, path="renderRequest", required=_TOP_LEVEL_FIELDS)
    schema = top["schema"]
    if not isinstance(schema, int) or isinstance(schema, bool):
        raise CanonicalRequestDecodingError("renderRequest.schema must be an integer.")
    if schema not in {1, 2, CANONICAL_REQUEST_SCHEMA}:
        raise CanonicalRequestDecodingError(
            f"Unsupported canonical request schema: {schema!r}."
        )
    mode = top["mode"]
    if mode not in {"circular", "linear"}:
        raise CanonicalRequestDecodingError(
            f"Unsupported canonical request mode: {mode!r}."
        )
    if output_directory is None or not isinstance(output_directory, (str, Path)):
        raise CanonicalRequestDecodingError(
            "A replay output directory must be supplied by the caller."
        )
    if not str(output_directory).strip():
        raise CanonicalRequestDecodingError(
            "A replay output directory must be supplied by the caller."
        )

    raw_records = top["records"]
    if not isinstance(raw_records, list) or not raw_records:
        raise CanonicalRequestDecodingError(
            "renderRequest.records must be a non-empty array."
        )
    records = tuple(
        _decode_record(value, index=index, schema=schema, resource_paths=resource_paths)
        for index, value in enumerate(raw_records, start=1)
    )
    options_kwargs = _decode_diagram_options(
        top["diagramOptions"], schema=schema, resource_paths=resource_paths
    )
    comparison_kwargs = _decode_comparisons(
        top["comparisons"], mode=mode, schema=schema, resource_paths=resource_paths
    )
    options = DiagramOptions(**options_kwargs, **comparison_kwargs)
    _validate_dataclass_contract(options, path="diagramOptions", error="decode")
    output = _decode_output(top["output"], output_directory=output_directory)

    if mode == "linear":
        if schema == 1:
            layout_payload = _object(top["layout"], path="renderRequest.layout")
            if layout_payload:
                raise CanonicalRequestDecodingError(
                    "A schema 1 Linear canonical request must have an empty layout object."
                )
            linear_layout = None
        else:
            linear_layout = _decode_linear_layout(top["layout"])
        return LinearDiagramRequest(
            records=records,
            options=options,
            layout=linear_layout,
            output=output,
        )

    layout = _decode_circular_layout(top["layout"])
    return CircularDiagramRequest(
        records=records,
        options=options,
        layout=layout,
        output=output,
    )


def _encode_record(
    record: RecordInput,
    *,
    index: int,
    resources: _ResourceBuilder,
) -> dict[str, Any]:
    source = record.source
    if isinstance(source, GenBankInputSource):
        resource_id = resources.add_path(
            f"record-{index}-genbank", kind="genbank", value=source.path
        )
        source_payload = {"kind": "genbank", "resourceId": resource_id}
    elif isinstance(source, GffFastaInputSource):
        gff_id = resources.add_path(
            f"record-{index}-gff3", kind="gff3", value=source.gff_path
        )
        fasta_id = resources.add_path(
            f"record-{index}-fasta", kind="fasta", value=source.fasta_path
        )
        source_payload = {
            "kind": "gffFasta",
            "gffResourceId": gff_id,
            "fastaResourceId": fasta_id,
        }
    elif isinstance(source, InMemoryRecordSource):
        stream = StringIO()
        try:
            serializable_record = deepcopy(source.record)
            serializable_record.annotations.setdefault("molecule_type", "DNA")
            SeqIO.write((serializable_record,), stream, "genbank")
        except Exception as exc:
            raise CanonicalRequestEncodingError(
                f"Could not serialize in-memory record {index} as GenBank."
            ) from exc
        resource_id = resources.add_bytes(
            f"record-{index}-genbank",
            kind="genbank",
            name=f"record-{index}.gbk",
            content=stream.getvalue().encode("utf-8"),
        )
        source_payload = {"kind": "genbank", "resourceId": resource_id}
    else:  # pragma: no cover - RecordInput validates its source union.
        raise CanonicalRequestEncodingError(f"Unsupported record source at index {index}.")

    presentation = record.presentation
    return {
        "recordKey": record.record_key or f"record-{index}",
        "source": source_payload,
        "selector": _encode_selector(record.selector),
        "region": _encode_region(record.region),
        "presentation": {
            "label": presentation.label,
            "subtitle": presentation.subtitle,
            "reverseComplement": presentation.reverse_complement,
            "gridRow": presentation.grid_row,
            "gridColumn": presentation.grid_column,
        },
    }


def _decode_record(
    value: object,
    *,
    index: int,
    schema: int,
    resource_paths: Mapping[str, str | Path],
) -> RecordInput:
    path = f"renderRequest.records[{index - 1}]"
    required = {"source", "selector", "region", "presentation"}
    if schema >= 2:
        required.add("recordKey")
    item = _object(value, path=path, required=required)
    source_payload = _object(
        item["source"], path=f"{path}.source", required={"kind"}, exact=False
    )
    kind = source_payload["kind"]
    if kind == "genbank":
        _require_exact_fields(
            source_payload, path=f"{path}.source", required={"kind", "resourceId"}
        )
        source = GenBankInputSource(
            _resolve_resource(
                source_payload["resourceId"],
                path=f"{path}.source.resourceId",
                resource_paths=resource_paths,
            )
        )
    elif kind == "gffFasta":
        _require_exact_fields(
            source_payload,
            path=f"{path}.source",
            required={"kind", "gffResourceId", "fastaResourceId"},
        )
        source = GffFastaInputSource(
            _resolve_resource(
                source_payload["gffResourceId"],
                path=f"{path}.source.gffResourceId",
                resource_paths=resource_paths,
            ),
            _resolve_resource(
                source_payload["fastaResourceId"],
                path=f"{path}.source.fastaResourceId",
                resource_paths=resource_paths,
            ),
        )
    else:
        raise CanonicalRequestDecodingError(
            f"Unsupported record source kind at {path}.source: {kind!r}."
        )

    presentation_payload = _object(
        item["presentation"],
        path=f"{path}.presentation",
        required={"label", "subtitle", "reverseComplement", "gridRow", "gridColumn"},
    )
    presentation = RecordPresentation(
        label=_optional_string(presentation_payload["label"], f"{path}.presentation.label"),
        subtitle=_optional_string(
            presentation_payload["subtitle"], f"{path}.presentation.subtitle"
        ),
        reverse_complement=_boolean(
            presentation_payload["reverseComplement"],
            f"{path}.presentation.reverseComplement",
        ),
        grid_row=_optional_integer(
            presentation_payload["gridRow"], f"{path}.presentation.gridRow"
        ),
        grid_column=_optional_integer(
            presentation_payload["gridColumn"], f"{path}.presentation.gridColumn"
        ),
    )
    return RecordInput(
        source=source,
        selector=_decode_selector(item["selector"], path=f"{path}.selector"),
        region=_decode_region(item["region"], path=f"{path}.region"),
        presentation=presentation,
        record_key=(
            _required_string(item["recordKey"], f"{path}.recordKey")
            if schema >= 2
            else f"record-{index}"
        ),
    )


def _encode_selector(selector: RecordSelector | None) -> dict[str, Any] | None:
    if selector is None:
        return None
    if selector.record_id is not None:
        return {"kind": "recordId", "value": selector.record_id}
    return {"kind": "recordIndex", "index": selector.record_index}


def _decode_selector(value: object, *, path: str) -> RecordSelector | None:
    if value is None:
        return None
    selector = _object(value, path=path, required={"kind"}, exact=False)
    kind = selector["kind"]
    if kind == "recordId":
        _require_exact_fields(selector, path=path, required={"kind", "value"})
        record_id = selector["value"]
        if not isinstance(record_id, str) or not record_id:
            raise CanonicalRequestDecodingError(f"{path}.value must be a non-empty string.")
        return RecordSelector(raw=record_id, record_id=record_id, record_index=None)
    if kind == "recordIndex":
        _require_exact_fields(selector, path=path, required={"kind", "index"})
        record_index = _integer(selector["index"], f"{path}.index")
        if record_index < 0:
            raise CanonicalRequestDecodingError(f"{path}.index must be non-negative.")
        return RecordSelector(
            raw=f"#{record_index + 1}", record_id=None, record_index=record_index
        )
    raise CanonicalRequestDecodingError(f"Unsupported selector kind at {path}: {kind!r}.")


def _encode_region(region: RegionSpec | None) -> dict[str, Any] | None:
    if region is None:
        return None
    selector: dict[str, Any] | None
    if region.record_id is not None:
        selector = {"kind": "recordId", "value": region.record_id}
    elif region.record_index is not None:
        selector = {"kind": "recordIndex", "index": region.record_index}
    else:
        selector = None
    return {
        "selector": selector,
        "start": region.start,
        "end": region.end,
        "reverseComplement": region.reverse_complement,
    }


def _decode_region(value: object, *, path: str) -> RegionSpec | None:
    if value is None:
        return None
    region = _object(
        value,
        path=path,
        required={"selector", "start", "end", "reverseComplement"},
    )
    selector = _decode_selector(region["selector"], path=f"{path}.selector")
    start = _integer(region["start"], f"{path}.start")
    end = _integer(region["end"], f"{path}.end")
    reverse = _boolean(region["reverseComplement"], f"{path}.reverseComplement")
    selector_text = selector.raw if selector is not None else ""
    coordinate_text = f"{start}-{end}" + (":rc" if reverse else "")
    raw = f"{selector_text}:{coordinate_text}" if selector_text else coordinate_text
    return RegionSpec(
        raw=raw,
        file_selector=None,
        record_id=selector.record_id if selector else None,
        record_index=selector.record_index if selector else None,
        start=start,
        end=end,
        reverse_complement=reverse,
    )


def _encode_layout(request: DiagramRequest) -> dict[str, Any]:
    if isinstance(request, LinearDiagramRequest):
        if request.layout is None:
            return {}
        return {
            "recordGapPx": request.layout.record_gap_px,
            "multiRecordPositions": (
                list(request.layout.multi_record_positions)
                if request.layout.multi_record_positions is not None
                else None
            ),
        }
    if not isinstance(request, CircularDiagramRequest) or request.layout is None:
        return {}
    layout = request.layout
    return {
        "multiRecordSizeMode": layout.multi_record_size_mode,
        "multiRecordMinRadiusRatio": layout.multi_record_min_radius_ratio,
        "multiRecordColumnGapRatio": layout.multi_record_column_gap_ratio,
        "multiRecordRowGapRatio": layout.multi_record_row_gap_ratio,
        "multiRecordPositions": (
            list(layout.multi_record_positions)
            if layout.multi_record_positions is not None
            else None
        ),
    }


def _decode_circular_layout(value: object) -> CircularMultiRecordOptions | None:
    layout = _object(value, path="renderRequest.layout")
    if not layout:
        return None
    required = {
        "multiRecordSizeMode",
        "multiRecordMinRadiusRatio",
        "multiRecordColumnGapRatio",
        "multiRecordRowGapRatio",
        "multiRecordPositions",
    }
    _require_exact_fields(layout, path="renderRequest.layout", required=required)
    positions = layout["multiRecordPositions"]
    if positions is not None and (
        not isinstance(positions, list)
        or not all(isinstance(item, str) for item in positions)
    ):
        raise CanonicalRequestDecodingError(
            "renderRequest.layout.multiRecordPositions must be a string array or null."
        )
    result = CircularMultiRecordOptions(
        multi_record_size_mode=layout["multiRecordSizeMode"],
        multi_record_min_radius_ratio=layout["multiRecordMinRadiusRatio"],
        multi_record_column_gap_ratio=layout["multiRecordColumnGapRatio"],
        multi_record_row_gap_ratio=layout["multiRecordRowGapRatio"],
        multi_record_positions=tuple(positions) if positions is not None else None,
    )
    _validate_dataclass_contract(result, path="layout", error="decode")
    return result


def _decode_linear_layout(value: object) -> LinearMultiRecordOptions | None:
    layout = _object(value, path="renderRequest.layout")
    if not layout:
        return None
    _require_exact_fields(
        layout,
        path="renderRequest.layout",
        required={"recordGapPx", "multiRecordPositions"},
    )
    positions = layout["multiRecordPositions"]
    if positions is not None and (
        not isinstance(positions, list)
        or not all(isinstance(item, str) for item in positions)
    ):
        raise CanonicalRequestDecodingError(
            "renderRequest.layout.multiRecordPositions must be a string array or null."
        )
    result = LinearMultiRecordOptions(
        record_gap_px=layout["recordGapPx"],
        multi_record_positions=tuple(positions) if positions is not None else None,
    )
    _validate_dataclass_contract(result, path="layout", error="decode")
    return result


def _encode_diagram_options(
    options: DiagramOptions,
    *,
    resources: _ResourceBuilder,
) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for item in fields(DiagramOptions):
        name = item.name
        if name in _COMPARISON_FIELDS:
            continue
        value = getattr(options, name)
        default = getattr(_DEFAULT_OPTIONS, name)
        if _same_default(value, default):
            continue
        result[_camel(name)] = _encode_option_value(name, value, resources=resources)
    return result


def _decode_diagram_options(
    value: object,
    *,
    schema: int,
    resource_paths: Mapping[str, str | Path],
) -> dict[str, Any]:
    payload = _object(value, path="renderRequest.diagramOptions")
    known = {
        _camel(item.name): item.name
        for item in fields(DiagramOptions)
        if item.name not in _COMPARISON_FIELDS
    }
    unknown = set(payload) - set(known)
    if unknown:
        raise CanonicalRequestDecodingError(
            "Unknown field(s) at renderRequest.diagramOptions: "
            + ", ".join(sorted(unknown))
            + "."
        )
    decoded = {
        known[key]: _decode_option_value(
            known[key], raw, resource_paths=resource_paths
        )
        for key, raw in payload.items()
    }
    if schema in {1, 2}:
        feature_shapes = dict(decoded.get("feature_shapes") or {})
        feature_shapes.setdefault("repeat_region", "rectangle")
        decoded["feature_shapes"] = feature_shapes
    return decoded


def _encode_option_value(
    name: str,
    value: object,
    *,
    resources: _ResourceBuilder,
) -> Any:
    if name == "config":
        if isinstance(value, GbdrawConfig):
            config = value
            value = asdict(config)
            value["labels"]["filtering"] = deepcopy(
                config.labels.filtering.as_dict()
            )
        return _json_value(value, path="diagramOptions.config")
    if name == "colors":
        return _encode_colors(value, resources=resources)
    if name == "tracks":
        return _encode_tracks(value)
    if name == "annotations":
        return _encode_annotations(value, resources=resources)
    if name == "output":
        return _encode_assembly_output(value)
    if name in _TABLE_FIELDS:
        return _table_ref(name, value, resources=resources)
    if name in _FILE_FIELDS:
        return _file_ref(name, value, resources=resources)
    if name in _TABLE_SEQUENCE_FIELDS:
        return [
            _table_ref(f"{name}-{index}", item, resources=resources)
            for index, item in enumerate(_sequence(value, name=name), start=1)
        ]
    if name in _FILE_SEQUENCE_FIELDS:
        return [
            _file_ref(f"{name}-{index}", item, resources=resources)
            for index, item in enumerate(_sequence(value, name=name), start=1)
        ]
    if name in _TABLE_MATRIX_FIELDS:
        return _encode_resource_matrix(name, value, table=True, resources=resources)
    if name in _FILE_MATRIX_FIELDS:
        return _encode_resource_matrix(name, value, table=False, resources=resources)
    return _json_value(value, path=f"diagramOptions.{_camel(name)}")


def _decode_option_value(
    name: str,
    value: object,
    *,
    resource_paths: Mapping[str, str | Path],
) -> Any:
    if name == "config":
        return dict(_object(value, path="renderRequest.diagramOptions.config"))
    if name == "colors":
        return _decode_colors(value, resource_paths=resource_paths)
    if name == "tracks":
        return _decode_tracks(value)
    if name == "annotations":
        return _decode_annotations(value, resource_paths=resource_paths)
    if name == "output":
        return _decode_assembly_output(value)
    if name in _TABLE_FIELDS:
        return _decode_table_ref(value, name=name, resource_paths=resource_paths)
    if name in _FILE_FIELDS:
        return str(_decode_file_ref(value, name=name, resource_paths=resource_paths))
    if name in _TABLE_SEQUENCE_FIELDS:
        raw = _array(value, path=f"renderRequest.diagramOptions.{_camel(name)}")
        return tuple(
            _decode_table_ref(item, name=f"{name}-{index}", resource_paths=resource_paths)
            for index, item in enumerate(raw, start=1)
        )
    if name in _FILE_SEQUENCE_FIELDS:
        raw = _array(value, path=f"renderRequest.diagramOptions.{_camel(name)}")
        return tuple(
            str(
                _decode_file_ref(
                    item, name=f"{name}-{index}", resource_paths=resource_paths
                )
            )
            for index, item in enumerate(raw, start=1)
        )
    if name in _TABLE_MATRIX_FIELDS:
        return _decode_resource_matrix(
            value, name=name, table=True, resource_paths=resource_paths
        )
    if name in _FILE_MATRIX_FIELDS:
        return _decode_resource_matrix(
            value, name=name, table=False, resource_paths=resource_paths
        )
    return value


def _encode_colors(value: object, *, resources: _ResourceBuilder) -> dict[str, Any]:
    if not isinstance(value, ColorOptions):
        raise CanonicalRequestEncodingError("diagramOptions.colors must be ColorOptions.")
    return {
        "colorTable": (
            _table_ref("colors-color-table", value.color_table, resources=resources)
            if value.color_table is not None
            else None
        ),
        "colorTableFile": (
            _file_ref("colors-color-table-file", value.color_table_file, resources=resources)
            if value.color_table_file is not None
            else None
        ),
        "defaultColors": (
            _table_ref("colors-default-colors", value.default_colors, resources=resources)
            if value.default_colors is not None
            else None
        ),
        "defaultColorsPalette": value.default_colors_palette,
        "defaultColorsFile": (
            _file_ref(
                "colors-default-colors-file", value.default_colors_file, resources=resources
            )
            if value.default_colors_file is not None
            else None
        ),
    }


def _decode_colors(
    value: object,
    *,
    resource_paths: Mapping[str, str | Path],
) -> ColorOptions:
    path = "renderRequest.diagramOptions.colors"
    payload = _object(
        value,
        path=path,
        required={
            "colorTable",
            "colorTableFile",
            "defaultColors",
            "defaultColorsPalette",
            "defaultColorsFile",
        },
    )
    result = ColorOptions(
        color_table=(
            _decode_table_ref(
                payload["colorTable"],
                name="colors-color-table",
                resource_paths=resource_paths,
            )
            if payload["colorTable"] is not None
            else None
        ),
        color_table_file=(
            str(
                _decode_file_ref(
                    payload["colorTableFile"],
                    name="colors-color-table-file",
                    resource_paths=resource_paths,
                )
            )
            if payload["colorTableFile"] is not None
            else None
        ),
        default_colors=(
            _decode_table_ref(
                payload["defaultColors"],
                name="colors-default-colors",
                resource_paths=resource_paths,
            )
            if payload["defaultColors"] is not None
            else None
        ),
        default_colors_palette=payload["defaultColorsPalette"],
        default_colors_file=(
            str(
                _decode_file_ref(
                    payload["defaultColorsFile"],
                    name="colors-default-colors-file",
                    resource_paths=resource_paths,
                )
            )
            if payload["defaultColorsFile"] is not None
            else None
        ),
    )
    _validate_dataclass_contract(result, path="diagramOptions.colors", error="decode")
    return result


def _encode_tracks(value: object) -> dict[str, Any]:
    if not isinstance(value, TrackOptions):
        raise CanonicalRequestEncodingError("diagramOptions.tracks must be TrackOptions.")
    return {
        "circularTrackSlots": _encode_track_slots(value.circular_track_slots),
        "circularTrackAxisIndex": value.circular_track_axis_index,
        "linearTrackSlots": _encode_track_slots(value.linear_track_slots),
        "linearTrackAxisIndex": value.linear_track_axis_index,
        "centerReservedRadius": value.center_reserved_radius,
    }


def _decode_tracks(value: object) -> TrackOptions:
    path = "renderRequest.diagramOptions.tracks"
    payload = _object(
        value,
        path=path,
        required={
            "circularTrackSlots",
            "circularTrackAxisIndex",
            "linearTrackSlots",
            "linearTrackAxisIndex",
            "centerReservedRadius",
        },
    )
    result = TrackOptions(
        circular_track_slots=_decode_track_slots(
            payload["circularTrackSlots"], mode="circular", path=f"{path}.circularTrackSlots"
        ),
        circular_track_axis_index=payload["circularTrackAxisIndex"],
        linear_track_slots=_decode_track_slots(
            payload["linearTrackSlots"], mode="linear", path=f"{path}.linearTrackSlots"
        ),
        linear_track_axis_index=payload["linearTrackAxisIndex"],
        center_reserved_radius=payload["centerReservedRadius"],
    )
    _validate_dataclass_contract(result, path="diagramOptions.tracks", error="decode")
    return result


def _encode_annotation_style(style: RegionAnnotationStyle) -> dict[str, Any]:
    return {
        "stroke": style.stroke,
        "strokeWidth": style.stroke_width,
        "strokeDasharray": list(style.stroke_dasharray),
        "lineCap": style.line_cap,
        "fill": style.fill,
        "fillOpacity": style.fill_opacity,
        "hatch": (
            {
                "angle": style.hatch.angle,
                "spacing": style.hatch.spacing,
                "color": style.hatch.color,
                "width": style.hatch.width,
                "cross": style.hatch.cross,
            }
            if style.hatch is not None
            else None
        ),
        "labelColor": style.label_color,
        "labelFontSize": style.label_font_size,
        "labelOrientation": style.label_orientation,
        "labelPosition": style.label_position,
        "labelOffset": style.label_offset,
    }


def _decode_annotation_style(value: object, *, path: str) -> RegionAnnotationStyle:
    payload = _object(
        value,
        path=path,
        required={
            "stroke", "strokeWidth", "strokeDasharray", "lineCap", "fill",
            "fillOpacity", "hatch", "labelColor", "labelFontSize",
            "labelOrientation", "labelPosition", "labelOffset",
        },
    )
    hatch_payload = payload["hatch"]
    hatch = None
    if hatch_payload is not None:
        raw_hatch = _object(
            hatch_payload,
            path=f"{path}.hatch",
            required={"angle", "spacing", "color", "width", "cross"},
        )
        hatch = HatchStyle(**raw_hatch)
    dasharray = _array(payload["strokeDasharray"], path=f"{path}.strokeDasharray")
    return RegionAnnotationStyle(
        stroke=payload["stroke"],
        stroke_width=payload["strokeWidth"],
        stroke_dasharray=tuple(dasharray),
        line_cap=payload["lineCap"],
        fill=payload["fill"],
        fill_opacity=payload["fillOpacity"],
        hatch=hatch,
        label_color=payload["labelColor"],
        label_font_size=payload["labelFontSize"],
        label_orientation=payload["labelOrientation"],
        label_position=payload["labelPosition"],
        label_offset=payload["labelOffset"],
    )


def _encode_annotations(value: object, *, resources: _ResourceBuilder) -> dict[str, Any]:
    if not isinstance(value, AnnotationOptions):
        raise CanonicalRequestEncodingError("diagramOptions.annotations must be AnnotationOptions.")
    sets = []
    for annotation_set in value.sets:
        items = []
        for annotation in annotation_set.annotations:
            target = annotation.target
            if isinstance(target, CoordinateSpan):
                target_payload = {
                    "kind": "coordinateSpan",
                    "record": _encode_selector(target.record),
                    "start": target.start,
                    "end": target.end,
                    "coordinateSpace": target.coordinate_space,
                    "wrapsOrigin": target.wraps_origin,
                    "outOfBounds": target.out_of_bounds,
                }
            elif isinstance(target, FeatureSpan):
                target_payload = {
                    "kind": "featureSpan",
                    "record": _encode_selector(target.record),
                    "selectors": [
                        {"key": selector.key, "value": selector.value}
                        for selector in target.selectors
                    ],
                    "envelope": target.envelope,
                    "circularPath": target.circular_path,
                }
            else:  # pragma: no cover - RegionAnnotation validates the union.
                raise CanonicalRequestEncodingError("Unsupported annotation target.")
            items.append(
                {
                    "id": annotation.id,
                    "target": target_payload,
                    "label": annotation.label,
                    "mark": annotation.mark,
                    "lane": annotation.lane,
                    "style": _encode_annotation_style(annotation.style) if annotation.style else None,
                    "legendLabel": annotation.legend_label,
                    "metadata": dict(annotation.metadata),
                }
            )
        sets.append(
            {
                "id": annotation_set.id,
                "annotations": items,
                "defaultStyle": _encode_annotation_style(annotation_set.default_style),
                "legendLabel": annotation_set.legend_label,
            }
        )
    return {
        "sets": sets,
        "table": (
            _table_ref("annotations-table", value.table, resources=resources)
            if value.table is not None else None
        ),
        "tableFile": (
            _file_ref("annotations-table-file", value.table_file, resources=resources)
            if value.table_file is not None else None
        ),
    }


def _decode_annotations(
    value: object,
    *,
    resource_paths: Mapping[str, str | Path],
) -> AnnotationOptions:
    path = "renderRequest.diagramOptions.annotations"
    payload = _object(value, path=path, required={"sets", "table", "tableFile"})
    raw_sets = _array(payload["sets"], path=f"{path}.sets")
    sets: list[AnnotationSet] = []
    for set_index, raw_set in enumerate(raw_sets):
        set_path = f"{path}.sets[{set_index}]"
        set_payload = _object(
            raw_set,
            path=set_path,
            required={"id", "annotations", "defaultStyle", "legendLabel"},
        )
        raw_items = _array(set_payload["annotations"], path=f"{set_path}.annotations")
        annotations: list[RegionAnnotation] = []
        for item_index, raw_item in enumerate(raw_items):
            item_path = f"{set_path}.annotations[{item_index}]"
            item = _object(
                raw_item,
                path=item_path,
                required={"id", "target", "label", "mark", "lane", "style", "legendLabel", "metadata"},
            )
            target_payload = _object(item["target"], path=f"{item_path}.target", required={"kind"}, exact=False)
            kind = target_payload["kind"]
            if kind == "coordinateSpan":
                _require_exact_fields(
                    target_payload,
                    path=f"{item_path}.target",
                    required={"kind", "record", "start", "end", "coordinateSpace", "wrapsOrigin", "outOfBounds"},
                )
                target = CoordinateSpan(
                    record=_decode_selector(target_payload["record"], path=f"{item_path}.target.record"),
                    start=target_payload["start"],
                    end=target_payload["end"],
                    coordinate_space=target_payload["coordinateSpace"],
                    wraps_origin=target_payload["wrapsOrigin"],
                    out_of_bounds=target_payload["outOfBounds"],
                )
            elif kind == "featureSpan":
                _require_exact_fields(
                    target_payload,
                    path=f"{item_path}.target",
                    required={"kind", "record", "selectors", "envelope", "circularPath"},
                )
                raw_selectors = _array(target_payload["selectors"], path=f"{item_path}.target.selectors")
                selectors = tuple(
                    FeatureSelector(**_object(raw, path=f"{item_path}.target.selectors[{index}]", required={"key", "value"}))
                    for index, raw in enumerate(raw_selectors)
                )
                target = FeatureSpan(
                    record=_decode_selector(target_payload["record"], path=f"{item_path}.target.record"),
                    selectors=selectors,
                    envelope=target_payload["envelope"],
                    circular_path=target_payload["circularPath"],
                )
            else:
                raise CanonicalRequestDecodingError(f"Unsupported annotation target kind at {item_path}: {kind!r}.")
            metadata = _object(item["metadata"], path=f"{item_path}.metadata")
            annotations.append(
                RegionAnnotation(
                    id=item["id"], target=target, label=item["label"], mark=item["mark"],
                    lane=item["lane"],
                    style=(_decode_annotation_style(item["style"], path=f"{item_path}.style") if item["style"] is not None else None),
                    legend_label=item["legendLabel"], metadata=metadata,
                )
            )
        sets.append(
            AnnotationSet(
                id=set_payload["id"],
                annotations=tuple(annotations),
                default_style=_decode_annotation_style(set_payload["defaultStyle"], path=f"{set_path}.defaultStyle"),
                legend_label=set_payload["legendLabel"],
            )
        )
    table = (
        _decode_table_ref(payload["table"], name="annotations-table", resource_paths=resource_paths)
        if payload["table"] is not None else None
    )
    table_file = (
        str(_decode_file_ref(payload["tableFile"], name="annotations-table-file", resource_paths=resource_paths))
        if payload["tableFile"] is not None else None
    )
    return AnnotationOptions(sets=tuple(sets), table=table, table_file=table_file)


def _encode_track_slots(value: object) -> list[Any] | None:
    if value is None:
        return None
    result: list[Any] = []
    for slot in _sequence(value, name="track slots"):
        if isinstance(slot, str):
            result.append(slot)
        elif isinstance(slot, CircularTrackSlot):
            result.append(_encode_track_slot(slot, kind="circularTrackSlot"))
        elif isinstance(slot, LinearTrackSlot):
            result.append(_encode_track_slot(slot, kind="linearTrackSlot"))
        else:
            raise CanonicalRequestEncodingError(
                f"Unsupported track slot value: {type(slot).__name__}."
            )
    return result


def _encode_track_slot(
    slot: CircularTrackSlot | LinearTrackSlot,
    *,
    kind: str,
) -> dict[str, Any]:
    result = {"kind": kind}
    for item in fields(slot):
        raw = getattr(slot, item.name)
        if isinstance(raw, ScalarSpec):
            raw = {"value": raw.value, "unit": raw.unit}
        elif item.name == "params" and isinstance(raw, MappingABC):
            raw = dict(raw)
            if isinstance(raw.get("style_override"), RegionAnnotationStyle):
                raw["style_override"] = _encode_annotation_style(raw["style_override"])
        result[_camel(item.name)] = _json_value(raw, path=f"trackSlot.{_camel(item.name)}")
    return result


def _decode_track_slots(
    value: object,
    *,
    mode: Literal["circular", "linear"],
    path: str,
) -> tuple[str | CircularTrackSlot | LinearTrackSlot, ...] | None:
    if value is None:
        return None
    raw_slots = _array(value, path=path)
    expected_kind = "circularTrackSlot" if mode == "circular" else "linearTrackSlot"
    cls = CircularTrackSlot if mode == "circular" else LinearTrackSlot
    result: list[str | CircularTrackSlot | LinearTrackSlot] = []
    for index, raw in enumerate(raw_slots):
        if isinstance(raw, str):
            result.append(raw)
            continue
        slot_path = f"{path}[{index}]"
        slot = _object(raw, path=slot_path, required={"kind"}, exact=False)
        if slot["kind"] != expected_kind:
            raise CanonicalRequestDecodingError(
                f"Unsupported track slot kind at {slot_path}: {slot['kind']!r}."
            )
        field_map = {_camel(item.name): item.name for item in fields(cls)}
        _require_exact_fields(
            slot, path=slot_path, required={"kind", *field_map.keys()}
        )
        scalar_fields = (
            {"radius", "width", "spacing"}
            if mode == "circular"
            else {"height", "spacing"}
        )
        kwargs = {
            name: (
                _decode_scalar_if_needed(slot[key], path=f"{slot_path}.{key}")
                if name in scalar_fields
                else slot[key]
            )
            for key, name in field_map.items()
        }
        if str(kwargs.get("renderer", "")).strip().lower() == "annotations":
            params = dict(kwargs.get("params") or {})
            if isinstance(params.get("style_override"), MappingABC):
                params["style_override"] = _decode_annotation_style(
                    params["style_override"], path=f"{slot_path}.params.style_override"
                )
            kwargs["params"] = params
        result.append(cls(**kwargs))
    return tuple(result)


def _decode_scalar_if_needed(value: object, *, path: str) -> object:
    if not isinstance(value, MappingABC) or "unit" not in value:
        return value
    scalar = _object(value, path=path, required={"value", "unit"})
    result = ScalarSpec(value=scalar["value"], unit=scalar["unit"])
    _validate_dataclass_contract(result, path=path, error="decode")
    return result


def _encode_assembly_output(value: object) -> dict[str, Any]:
    if not isinstance(value, OutputOptions):
        raise CanonicalRequestEncodingError("diagramOptions.output must be OutputOptions.")
    return {
        "outputPrefix": value.output_prefix,
        "legend": value.legend,
        "plotTitlePosition": value.plot_title_position,
    }


def _decode_assembly_output(value: object) -> OutputOptions:
    payload = _object(
        value,
        path="renderRequest.diagramOptions.output",
        required={"outputPrefix", "legend", "plotTitlePosition"},
    )
    result = OutputOptions(
        output_prefix=payload["outputPrefix"],
        legend=payload["legend"],
        plot_title_position=payload["plotTitlePosition"],
    )
    _validate_dataclass_contract(result, path="diagramOptions.output", error="decode")
    return result


def _encode_resource_matrix(
    name: str,
    value: object,
    *,
    table: bool,
    resources: _ResourceBuilder,
) -> list[list[dict[str, str] | None]]:
    result: list[list[dict[str, str] | None]] = []
    for row_index, row in enumerate(_sequence(value, name=name), start=1):
        encoded_row: list[dict[str, str] | None] = []
        for column_index, item in enumerate(_sequence(row, name=name), start=1):
            if item is None:
                encoded_row.append(None)
                continue
            resource_name = f"{name}-{row_index}-{column_index}"
            encoded_row.append(
                _table_ref(resource_name, item, resources=resources)
                if table
                else _file_ref(resource_name, item, resources=resources)
            )
        result.append(encoded_row)
    return result


def _decode_resource_matrix(
    value: object,
    *,
    name: str,
    table: bool,
    resource_paths: Mapping[str, str | Path],
) -> tuple[tuple[DataFrame | str | None, ...], ...]:
    rows = _array(value, path=f"renderRequest.diagramOptions.{_camel(name)}")
    result: list[tuple[DataFrame | str | None, ...]] = []
    for row_index, raw_row in enumerate(rows, start=1):
        row = _array(
            raw_row,
            path=f"renderRequest.diagramOptions.{_camel(name)}[{row_index - 1}]",
        )
        decoded_row: list[DataFrame | str | None] = []
        for column_index, item in enumerate(row, start=1):
            if item is None:
                decoded_row.append(None)
                continue
            resource_name = f"{name}-{row_index}-{column_index}"
            if table:
                decoded_row.append(
                    _decode_table_ref(
                        item, name=resource_name, resource_paths=resource_paths
                    )
                )
            else:
                decoded_row.append(
                    str(
                        _decode_file_ref(
                            item, name=resource_name, resource_paths=resource_paths
                        )
                    )
                )
        result.append(tuple(decoded_row))
    return tuple(result)


def _encode_comparisons(
    options: DiagramOptions,
    *,
    mode: Literal["circular", "linear"],
    resources: _ResourceBuilder,
) -> list[dict[str, Any]]:
    if mode == "circular":
        return []
    result: list[dict[str, Any]] = []
    for index, comparison in enumerate(options.linear_comparisons or (), start=1):
        ref = _table_ref(
            f"comparison-explicit-{index}", comparison.matches, resources=resources
        )
        result.append(
            {
                "kind": "precomputedProteinComparison",
                "resourceId": ref["resourceId"],
                "encoding": "canonicalTsv",
                "queryRecordIndex": comparison.query_record_index,
                "subjectRecordIndex": comparison.subject_record_index,
            }
        )
    for index, path in enumerate(options.blast_files or (), start=1):
        resource_id = resources.add_path(
            f"comparison-nucleotide-{index}", kind="nucleotide-blast", value=path
        )
        result.append(
            {
                "kind": "nucleotideBlast",
                "resourceId": resource_id,
                "queryRecordIndex": index - 1,
                "subjectRecordIndex": index,
            }
        )
    for index, table in enumerate(options.protein_comparisons or (), start=1):
        ref = _table_ref(f"comparison-protein-{index}", table, resources=resources)
        result.append(
            {
                "kind": "precomputedProteinComparison",
                "resourceId": ref["resourceId"],
                "encoding": "canonicalTsv",
                "queryRecordIndex": index - 1,
                "subjectRecordIndex": index,
            }
        )
    if options.orthogroups is not None:
        resource_id = _typed_json_resource(
            "comparison-orthogroups",
            kind="orthogroup-result",
            value_kind="orthogroupResult",
            value=options.orthogroups,
            resources=resources,
        )
        result.append(
            {
                "kind": "orthogroupResult",
                "resourceId": resource_id,
                "encoding": "canonicalJson",
            }
        )
    if options.collinearity_blocks is not None:
        value_kind = (
            "result"
            if isinstance(options.collinearity_blocks, CollinearityResult)
            else "blocks"
        )
        resource_id = _typed_json_resource(
            "comparison-collinearity",
            kind="collinearity-result",
            value_kind=value_kind,
            value=options.collinearity_blocks,
            resources=resources,
        )
        result.append(
            {
                "kind": "collinearityResult",
                "resourceId": resource_id,
                "encoding": "canonicalJson",
                "valueKind": value_kind,
            }
        )
    if any(
        not _same_default(getattr(options, name), getattr(_DEFAULT_OPTIONS, name))
        for name in _PIPELINE_FIELDS
    ) or options.protein_comparison_pairs is not None:
        settings = {
            _camel(name): _encode_pipeline_value(name, getattr(options, name))
            for name in _PIPELINE_FIELDS
            if name != "protein_blastp_mode"
        }
        result.append(
            {
                "kind": "generatedProteinComparison",
                "mode": options.protein_blastp_mode,
                "pairs": [
                    {
                        "queryRecordIndex": int(pair[0]),
                        "subjectRecordIndex": int(pair[1]),
                    }
                    for pair in (options.protein_comparison_pairs or ())
                ],
                "settings": settings,
            }
        )
    return result


def _decode_comparisons(
    value: object,
    *,
    mode: Literal["circular", "linear"],
    schema: int,
    resource_paths: Mapping[str, str | Path],
) -> dict[str, Any]:
    comparisons = _array(value, path="renderRequest.comparisons")
    if mode == "circular" and comparisons:
        raise CanonicalRequestDecodingError(
            "A circular canonical request cannot contain linear comparisons."
        )
    blast_files: list[str] = []
    protein_tables: list[DataFrame] = []
    explicit_comparisons: list[LinearComparison] = []
    result: dict[str, Any] = {}
    singleton_kinds: set[str] = set()
    for index, raw in enumerate(comparisons):
        path = f"renderRequest.comparisons[{index}]"
        item = _object(raw, path=path, required={"kind"}, exact=False)
        kind = item["kind"]
        if kind == "nucleotideBlast":
            required = {"kind", "resourceId"}
            if schema >= 2:
                required |= {"queryRecordIndex", "subjectRecordIndex"}
            _require_exact_fields(item, path=path, required=required)
            resource_path = _resolve_resource(
                item["resourceId"],
                path=f"{path}.resourceId",
                resource_paths=resource_paths,
            )
            if schema >= 2:
                try:
                    table = read_csv(
                        resource_path,
                        sep="\t",
                        comment="#",
                        names=COMPARISON_COLUMNS,
                    )
                except Exception as exc:
                    raise CanonicalRequestDecodingError(
                        f"Could not decode BLAST resource for {path}."
                    ) from exc
                query_index = _non_negative_index(
                    item["queryRecordIndex"], f"{path}.queryRecordIndex"
                )
                subject_index = _non_negative_index(
                    item["subjectRecordIndex"], f"{path}.subjectRecordIndex"
                )
                if query_index == len(blast_files) and subject_index == query_index + 1:
                    blast_files.append(str(resource_path))
                else:
                    explicit_comparisons.append(
                        LinearComparison(query_index, subject_index, table)
                    )
            else:
                blast_files.append(str(resource_path))
        elif kind == "precomputedProteinComparison":
            required = {"kind", "resourceId", "encoding"}
            if schema >= 2:
                required |= {"queryRecordIndex", "subjectRecordIndex"}
            _require_exact_fields(item, path=path, required=required)
            if item["encoding"] != "canonicalTsv":
                raise CanonicalRequestDecodingError(
                    f"Unsupported protein comparison encoding at {path}."
                )
            table = _read_canonical_table(
                _resolve_resource(
                    item["resourceId"],
                    path=f"{path}.resourceId",
                    resource_paths=resource_paths,
                ),
                context=path,
            )
            if schema >= 2:
                query_index = _non_negative_index(
                    item["queryRecordIndex"], f"{path}.queryRecordIndex"
                )
                subject_index = _non_negative_index(
                    item["subjectRecordIndex"], f"{path}.subjectRecordIndex"
                )
                if query_index == len(protein_tables) and subject_index == query_index + 1:
                    protein_tables.append(table)
                else:
                    explicit_comparisons.append(
                        LinearComparison(query_index, subject_index, table)
                    )
            else:
                protein_tables.append(table)
        elif kind == "orthogroupResult":
            _check_singleton_kind(kind, singleton_kinds, path=path)
            _require_exact_fields(
                item, path=path, required={"kind", "resourceId", "encoding"}
            )
            if item["encoding"] != "canonicalJson":
                raise CanonicalRequestDecodingError(
                    f"Unsupported orthogroup encoding at {path}."
                )
            result["orthogroups"] = _read_typed_json_resource(
                item["resourceId"],
                value_kind="orthogroupResult",
                expected=OrthogroupResult,
                path=path,
                resource_paths=resource_paths,
            )
        elif kind == "collinearityResult":
            _check_singleton_kind(kind, singleton_kinds, path=path)
            _require_exact_fields(
                item,
                path=path,
                required={"kind", "resourceId", "encoding", "valueKind"},
            )
            if item["encoding"] != "canonicalJson" or item["valueKind"] not in {
                "result",
                "blocks",
            }:
                raise CanonicalRequestDecodingError(
                    f"Unsupported collinearity representation at {path}."
                )
            expected = CollinearityResult if item["valueKind"] == "result" else tuple[CollinearityBlock, ...]
            result["collinearity_blocks"] = _read_typed_json_resource(
                item["resourceId"],
                value_kind=item["valueKind"],
                expected=expected,
                path=path,
                resource_paths=resource_paths,
            )
        elif kind == "generatedProteinComparison":
            _check_singleton_kind(kind, singleton_kinds, path=path)
            result.update(_decode_pipeline(item, path=path, schema=schema))
        else:
            raise CanonicalRequestDecodingError(
                f"Unsupported comparison kind at {path}: {kind!r}."
            )
    if blast_files:
        result["blast_files"] = tuple(blast_files)
    if protein_tables:
        result["protein_comparisons"] = tuple(protein_tables)
    if explicit_comparisons:
        result["linear_comparisons"] = tuple(explicit_comparisons)
    return result


def _encode_pipeline_value(name: str, value: object) -> Any:
    if name == "collinearity_params":
        if value is None:
            return None
        if isinstance(value, CollinearityParameters):
            kind = "standard"
        elif isinstance(value, LosslessCollinearityParameters):
            kind = "lossless"
        else:
            raise CanonicalRequestEncodingError(
                "Unsupported collinearity parameter object."
            )
        value.validate()
        return {
            "kind": kind,
            "parameters": {
                _camel(item.name): _json_value(
                    getattr(value, item.name), path=f"collinearityParams.{item.name}"
                )
                for item in fields(value)
            },
        }
    return _json_value(value, path=f"comparisons.settings.{_camel(name)}")


def _decode_pipeline(
    item: Mapping[str, Any], *, path: str, schema: int
) -> dict[str, Any]:
    required = {"kind", "mode", "settings"}
    if schema >= 2:
        required.add("pairs")
    _require_exact_fields(item, path=path, required=required)
    settings = _object(item["settings"], path=f"{path}.settings")
    setting_fields = tuple(name for name in _PIPELINE_FIELDS if name != "protein_blastp_mode")
    field_map = {_camel(name): name for name in setting_fields}
    _require_exact_fields(settings, path=f"{path}.settings", required=set(field_map))
    mode = item["mode"]
    result = {"protein_blastp_mode": mode}
    if schema >= 2:
        pairs = _array(item["pairs"], path=f"{path}.pairs")
        decoded_pairs: list[tuple[int, int]] = []
        for index, raw_pair in enumerate(pairs):
            pair_path = f"{path}.pairs[{index}]"
            pair = _object(
                raw_pair,
                path=pair_path,
                required={"queryRecordIndex", "subjectRecordIndex"},
            )
            decoded_pairs.append(
                (
                    _non_negative_index(
                        pair["queryRecordIndex"], f"{pair_path}.queryRecordIndex"
                    ),
                    _non_negative_index(
                        pair["subjectRecordIndex"], f"{pair_path}.subjectRecordIndex"
                    ),
                )
            )
        # Early schema-2 writers also stored the derived row-adjacent search pairs
        # used by collinear mode here.  The public option is pairwise-only; current
        # collinear rendering derives those pairs from the saved layout instead.
        result["protein_comparison_pairs"] = (
            tuple(decoded_pairs) if decoded_pairs and mode == "pairwise" else None
        )
    for key, name in field_map.items():
        raw = settings[key]
        result[name] = (
            _decode_collinearity_params(raw, path=f"{path}.settings.{key}")
            if name == "collinearity_params"
            else raw
        )
    return result


def _decode_collinearity_params(value: object, *, path: str) -> object:
    if value is None:
        return None
    payload = _object(value, path=path, required={"kind", "parameters"})
    kind = payload["kind"]
    cls: type[CollinearityParameters] | type[LosslessCollinearityParameters]
    if kind == "standard":
        cls = CollinearityParameters
    elif kind == "lossless":
        cls = LosslessCollinearityParameters
    else:
        raise CanonicalRequestDecodingError(
            f"Unsupported collinearity parameter kind at {path}: {kind!r}."
        )
    parameters = _object(payload["parameters"], path=f"{path}.parameters")
    field_map = {_camel(item.name): item.name for item in fields(cls)}
    _require_exact_fields(
        parameters, path=f"{path}.parameters", required=set(field_map)
    )
    result = cls(**{name: parameters[key] for key, name in field_map.items()})
    _validate_dataclass_contract(result, path=path, error="decode")
    result.validate()
    return result


def _decode_output(
    value: object,
    *,
    output_directory: str | Path,
) -> RenderOutputRequest:
    payload = _object(
        value,
        path="renderRequest.output",
        required={"prefix", "formats", "overwrite", "interactiveMetadataPolicy"},
    )
    formats = payload["formats"]
    if not isinstance(formats, list):
        raise CanonicalRequestDecodingError(
            "renderRequest.output.formats must be an array."
        )
    return RenderOutputRequest(
        output_prefix=payload["prefix"],
        output_directory=output_directory,
        formats=tuple(formats),
        overwrite=payload["overwrite"],
        interactive_metadata_policy=payload["interactiveMetadataPolicy"],
    )


def _table_ref(
    name: str,
    value: object,
    *,
    resources: _ResourceBuilder,
) -> dict[str, str]:
    if not isinstance(value, DataFrame):
        raise CanonicalRequestEncodingError(f"{name} must be a pandas DataFrame.")
    if len(value.columns) == 0:
        raise CanonicalRequestEncodingError(
            f"DataFrame resource {name!r} must have at least one column."
        )
    try:
        text = value.to_csv(sep="\t", index=False, lineterminator="\n")
    except Exception as exc:
        raise CanonicalRequestEncodingError(
            f"Could not encode DataFrame resource {name!r} as canonical TSV."
        ) from exc
    resource_id = _resource_id(name)
    resources.add_bytes(
        resource_id,
        kind="canonical-tsv",
        name=f"{resource_id}.tsv",
        content=text.encode("utf-8"),
    )
    return {"resourceId": resource_id, "representation": "canonicalTsv"}


def _file_ref(
    name: str,
    value: object,
    *,
    resources: _ResourceBuilder,
) -> dict[str, str]:
    resource_id = _resource_id(name)
    resources.add_path(resource_id, kind=_resource_id(name), value=value)
    return {"resourceId": resource_id, "representation": "file"}


def _decode_table_ref(
    value: object,
    *,
    name: str,
    resource_paths: Mapping[str, str | Path],
) -> DataFrame:
    ref = _resource_ref(value, path=f"resource reference {name}", representation="canonicalTsv")
    return _read_canonical_table(
        _resolve_resource(
            ref["resourceId"],
            path=f"resource reference {name}.resourceId",
            resource_paths=resource_paths,
        ),
        context=name,
    )


def _decode_file_ref(
    value: object,
    *,
    name: str,
    resource_paths: Mapping[str, str | Path],
) -> Path:
    ref = _resource_ref(value, path=f"resource reference {name}", representation="file")
    return _resolve_resource(
        ref["resourceId"],
        path=f"resource reference {name}.resourceId",
        resource_paths=resource_paths,
    )


def _resource_ref(
    value: object,
    *,
    path: str,
    representation: str,
) -> Mapping[str, Any]:
    ref = _object(
        value, path=path, required={"resourceId", "representation"}
    )
    if ref["representation"] != representation:
        raise CanonicalRequestDecodingError(
            f"Unsupported resource representation at {path}: {ref['representation']!r}."
        )
    return ref


def _read_canonical_table(file_path: Path, *, context: str) -> DataFrame:
    try:
        return read_csv(file_path, sep="\t")
    except Exception as exc:
        raise CanonicalRequestDecodingError(
            f"Could not decode canonical TSV resource for {context}."
        ) from exc


def _typed_json_resource(
    resource_id: str,
    *,
    kind: str,
    value_kind: str,
    value: object,
    resources: _ResourceBuilder,
) -> str:
    body = {
        "schema": 1,
        "kind": value_kind,
        "value": _encode_typed_tree(value),
    }
    try:
        content = json.dumps(
            body,
            ensure_ascii=False,
            sort_keys=True,
            separators=(",", ":"),
            allow_nan=False,
        ).encode("utf-8")
    except (TypeError, ValueError) as exc:
        raise CanonicalRequestEncodingError(
            f"Could not encode typed resource {resource_id!r}."
        ) from exc
    resources.add_bytes(
        resource_id,
        kind=kind,
        name=f"{resource_id}.json",
        content=content,
    )
    return resource_id


def _read_typed_json_resource(
    resource_id: object,
    *,
    value_kind: str,
    expected: object,
    path: str,
    resource_paths: Mapping[str, str | Path],
) -> Any:
    resource_path = _resolve_resource(
        resource_id, path=f"{path}.resourceId", resource_paths=resource_paths
    )
    try:
        raw = json.loads(resource_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise CanonicalRequestDecodingError(
            f"Could not decode canonical JSON resource for {path}."
        ) from exc
    payload = _object(
        raw, path=f"{path} resource", required={"schema", "kind", "value"}
    )
    if payload["schema"] != 1 or payload["kind"] != value_kind:
        raise CanonicalRequestDecodingError(
            f"Canonical JSON resource metadata does not match {path}."
        )
    return _decode_typed_tree(payload["value"], expected, path=f"{path}.value")


def _encode_typed_tree(value: object) -> Any:
    if is_dataclass(value) and not isinstance(value, type):
        cls = type(value)
        if cls.__name__ not in _TYPED_TREE_CLASSES:
            raise CanonicalRequestEncodingError(
                f"Unsupported typed resource value: {cls.__name__}."
            )
        return {
            "type": cls.__name__,
            "fields": {
                _camel(item.name): _encode_typed_tree(getattr(value, item.name))
                for item in fields(value)
            },
        }
    if isinstance(value, MappingABC):
        if not all(isinstance(key, str) for key in value):
            raise CanonicalRequestEncodingError(
                "Canonical typed resource mappings require string keys."
            )
        return {key: _encode_typed_tree(item) for key, item in value.items()}
    if isinstance(value, SequenceABC) and not isinstance(value, (str, bytes, bytearray)):
        return [_encode_typed_tree(item) for item in value]
    return _json_value(value, path="typed resource")


def _decode_typed_tree(value: object, hint: object, *, path: str) -> Any:
    if hint is Any or hint is object:
        return value
    origin = get_origin(hint)
    args = get_args(hint)
    if origin in {types.UnionType, Union}:
        if value is None and type(None) in args:
            return None
        tagged_type = value.get("type") if isinstance(value, MappingABC) else None
        for branch in args:
            if branch is type(None):
                continue
            if isinstance(branch, type) and branch.__name__ == tagged_type:
                return _decode_typed_tree(value, branch, path=path)
            branch_origin = get_origin(branch)
            if branch_origin in {tuple, list, Sequence, SequenceABC} and isinstance(value, list):
                return _decode_typed_tree(value, branch, path=path)
            if _matches_type(value, branch):
                return _decode_typed_tree(value, branch, path=path)
        raise CanonicalRequestDecodingError(
            f"Typed resource value at {path} does not match its union contract."
        )
    if isinstance(hint, type) and is_dataclass(hint):
        tagged = _object(value, path=path, required={"type", "fields"})
        if tagged["type"] != hint.__name__ or hint.__name__ not in _TYPED_TREE_CLASSES:
            raise CanonicalRequestDecodingError(
                f"Unsupported typed resource object at {path}: {tagged['type']!r}."
            )
        raw_fields = _object(tagged["fields"], path=f"{path}.fields")
        hints = get_type_hints(hint)
        field_map = {_camel(item.name): item.name for item in fields(hint)}
        _require_exact_fields(raw_fields, path=f"{path}.fields", required=set(field_map))
        kwargs = {
            name: _decode_typed_tree(
                raw_fields[key], hints[name], path=f"{path}.fields.{key}"
            )
            for key, name in field_map.items()
        }
        return hint(**kwargs)
    if origin in {tuple, list, Sequence, SequenceABC}:
        raw = _array(value, path=path)
        item_hint = args[0] if args else Any
        decoded = [
            _decode_typed_tree(item, item_hint, path=f"{path}[{index}]")
            for index, item in enumerate(raw)
        ]
        return tuple(decoded) if origin is tuple else decoded
    if origin in {dict, Mapping, MappingABC}:
        raw = _object(value, path=path)
        key_hint, value_hint = args if len(args) == 2 else (str, Any)
        if key_hint is not str:
            raise CanonicalRequestDecodingError(
                f"Unsupported typed resource mapping key contract at {path}."
            )
        return {
            key: _decode_typed_tree(item, value_hint, path=f"{path}.{key}")
            for key, item in raw.items()
        }
    if origin is Literal:
        if value not in args:
            raise CanonicalRequestDecodingError(
                f"Typed resource literal at {path} has unsupported value {value!r}."
            )
        return value
    if not _matches_type(value, hint):
        raise CanonicalRequestDecodingError(
            f"Typed resource value at {path} has an invalid type."
        )
    return value


def _resolve_resource(
    resource_id: object,
    *,
    path: str,
    resource_paths: Mapping[str, str | Path],
) -> Path:
    if not isinstance(resource_id, str) or not _RESOURCE_ID_RE.fullmatch(resource_id):
        raise CanonicalRequestDecodingError(f"Invalid resource ID at {path}.")
    if resource_id not in resource_paths:
        raise CanonicalRequestDecodingError(
            f"Canonical resource {resource_id!r} referenced at {path} is missing."
        )
    raw_path = resource_paths[resource_id]
    if not isinstance(raw_path, (str, Path)) or not str(raw_path).strip():
        raise CanonicalRequestDecodingError(
            f"Materialized path for resource {resource_id!r} is invalid."
        )
    materialized = Path(str(raw_path))
    if not materialized.is_file():
        raise CanonicalRequestDecodingError(
            f"Materialized resource is not a file: {resource_id!r}."
        )
    return materialized


def _object(
    value: object,
    *,
    path: str,
    required: set[str] | frozenset[str] | None = None,
    exact: bool = True,
) -> Mapping[str, Any]:
    if not isinstance(value, MappingABC) or not all(
        isinstance(key, str) for key in value
    ):
        raise CanonicalRequestDecodingError(f"{path} must be an object.")
    if required is not None:
        missing = set(required) - set(value)
        if missing:
            raise CanonicalRequestDecodingError(
                f"Missing required field(s) at {path}: {', '.join(sorted(missing))}."
            )
        if exact:
            _require_exact_fields(value, path=path, required=set(required))
    return value


def _require_exact_fields(
    value: Mapping[str, Any],
    *,
    path: str,
    required: set[str],
) -> None:
    missing = required - set(value)
    if missing:
        raise CanonicalRequestDecodingError(
            f"Missing required field(s) at {path}: {', '.join(sorted(missing))}."
        )
    unknown = set(value) - required
    if unknown:
        raise CanonicalRequestDecodingError(
            f"Unknown field(s) at {path}: {', '.join(sorted(unknown))}."
        )


def _array(value: object, *, path: str) -> list[Any]:
    if not isinstance(value, list):
        raise CanonicalRequestDecodingError(f"{path} must be an array.")
    return value


def _sequence(value: object, *, name: str) -> Sequence[Any]:
    if not isinstance(value, SequenceABC) or isinstance(value, (str, bytes, bytearray)):
        raise CanonicalRequestEncodingError(f"{name} must be a sequence.")
    return value


def _json_value(value: object, *, path: str) -> Any:
    if value is None or isinstance(value, (str, bool)):
        return value
    if isinstance(value, int) and not isinstance(value, bool):
        return value
    if isinstance(value, float):
        if not math.isfinite(value):
            raise CanonicalRequestEncodingError(f"{path} must be finite.")
        return value
    if isinstance(value, MappingABC):
        if not all(isinstance(key, str) for key in value):
            raise CanonicalRequestEncodingError(f"{path} requires string mapping keys.")
        return {
            key: _json_value(item, path=f"{path}.{key}")
            for key, item in value.items()
        }
    if isinstance(value, SequenceABC) and not isinstance(value, (str, bytes, bytearray)):
        return [
            _json_value(item, path=f"{path}[{index}]")
            for index, item in enumerate(value)
        ]
    raise CanonicalRequestEncodingError(
        f"{path} contains unsupported value type {type(value).__name__}."
    )


def _validate_dataclass_contract(value: object, *, path: str, error: str) -> None:
    if not is_dataclass(value) or isinstance(value, type):
        _raise_contract_error(error, f"{path} must be a dataclass value.")
    hints = get_type_hints(type(value))
    for item in fields(value):
        raw = getattr(value, item.name)
        hint = hints.get(item.name, Any)
        if not _matches_type(raw, hint):
            _raise_contract_error(
                error,
                f"{path}.{_camel(item.name)} does not match its typed contract.",
            )
        if is_dataclass(raw) and not isinstance(raw, type):
            _validate_dataclass_contract(
                raw, path=f"{path}.{_camel(item.name)}", error=error
            )


def _matches_type(value: object, hint: object) -> bool:
    if hint is Any or hint is object:
        return True
    origin = get_origin(hint)
    args = get_args(hint)
    if origin in {types.UnionType, Union}:
        return any(_matches_type(value, branch) for branch in args)
    if origin is Literal:
        return value in args and any(type(value) is type(item) for item in args if item == value)
    if origin in {tuple, list, Sequence, SequenceABC}:
        if not isinstance(value, SequenceABC) or isinstance(value, (str, bytes, bytearray)):
            return False
        if not args:
            return True
        item_hint = args[0]
        return all(_matches_type(item, item_hint) for item in value)
    if origin in {dict, Mapping, MappingABC}:
        if not isinstance(value, MappingABC):
            return False
        if len(args) != 2:
            return True
        return all(
            _matches_type(key, args[0]) and _matches_type(item, args[1])
            for key, item in value.items()
        )
    if hint is float:
        return (
            isinstance(value, (int, float))
            and not isinstance(value, bool)
            and math.isfinite(float(value))
        )
    if hint is int:
        return isinstance(value, int) and not isinstance(value, bool)
    if hint is bool:
        return isinstance(value, bool)
    if hint is str:
        return isinstance(value, str)
    if hint is type(None):
        return value is None
    return isinstance(hint, type) and isinstance(value, hint)


def _raise_contract_error(kind: str, message: str) -> None:
    if kind == "encode":
        raise CanonicalRequestEncodingError(message)
    raise CanonicalRequestDecodingError(message)


def _same_default(value: object, default: object) -> bool:
    if value is default:
        return True
    if isinstance(value, DataFrame) or isinstance(default, DataFrame):
        return False
    try:
        return bool(value == default)
    except (TypeError, ValueError):
        return False


def _resource_id(value: str) -> str:
    normalized = re.sub(r"[^a-z0-9]+", "-", value.lower()).strip("-")
    if not normalized or not normalized[0].isalpha():
        normalized = f"resource-{normalized}".rstrip("-")
    return normalized


def _camel(value: str) -> str:
    head, *tail = value.split("_")
    return head + "".join(part[:1].upper() + part[1:] for part in tail)


def _integer(value: object, path: str) -> int:
    if not isinstance(value, int) or isinstance(value, bool):
        raise CanonicalRequestDecodingError(f"{path} must be an integer.")
    return value


def _optional_integer(value: object, path: str) -> int | None:
    return None if value is None else _integer(value, path)


def _non_negative_index(value: object, path: str) -> int:
    index = _integer(value, path)
    if index < 0:
        raise CanonicalRequestDecodingError(f"{path} must be non-negative.")
    return index


def _boolean(value: object, path: str) -> bool:
    if not isinstance(value, bool):
        raise CanonicalRequestDecodingError(f"{path} must be a boolean.")
    return value


def _optional_string(value: object, path: str) -> str | None:
    if value is not None and not isinstance(value, str):
        raise CanonicalRequestDecodingError(f"{path} must be a string or null.")
    return value


def _required_string(value: object, path: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise CanonicalRequestDecodingError(f"{path} must be a non-empty string.")
    return value.strip()


def _check_singleton_kind(kind: str, seen: set[str], *, path: str) -> None:
    if kind in seen:
        raise CanonicalRequestDecodingError(
            f"Duplicate singleton comparison kind at {path}: {kind}."
        )
    seen.add(kind)


__all__ = [
    "CANONICAL_REQUEST_SCHEMA",
    "UNKNOWN_FIELD_POLICY",
    "CanonicalRequestCodecError",
    "CanonicalRequestDecodingError",
    "CanonicalRequestEncodingError",
    "CanonicalRequestResource",
    "EncodedCanonicalRequest",
    "decode_canonical_request",
    "encode_canonical_request",
]
