"""CLI- and session-independent request value objects.

These models describe unresolved or materialized record inputs and render intent.
They intentionally do not parse session documents or expose CLI argument names.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from enum import Enum
from pathlib import Path, PureWindowsPath
from typing import Literal, NamedTuple, Sequence, TypeAlias, cast

from Bio.SeqRecord import SeqRecord

from gbdraw.core.record_metadata import WINDOW_OF_ANNOTATION
from gbdraw.exceptions import ValidationError
from gbdraw.features.source import SourceFeatureCatalog
from gbdraw.io.record_select import RecordSelector
from gbdraw.io.regions import RecordWindow, RegionSpec
from gbdraw.layout.similarity_alignment import (
    SimilarityAlignmentPlan,
    SimilarityAlignmentReference,
    SimilarityAlignmentReferenceError,
)
from gbdraw.render.export_limits import checked_png_dpi
from gbdraw.render.formats import parse_output_formats
from gbdraw.render.output_paths import is_windows_reserved_filename_component

from .options import (
    CircularDiagramOptions,
    CircularMultiRecordOptions,
    LinearDiagramOptions,
    LinearMultiRecordOptions,
    resolve_circular_diagram_options,
    resolve_linear_diagram_options,
)


@dataclass(frozen=True)
class _LegacySimilarityAlignment:
    """Reader-only alignment value carried only by supported old schemas."""

    target: str
    source_schema: int

    def __post_init__(self) -> None:
        if not isinstance(self.target, str) or not self.target.strip() or "\0" in self.target:
            raise ValidationError(
                "Legacy similarity alignment target must be non-empty text without NUL."
            )
        if (
            isinstance(self.source_schema, bool)
            or not isinstance(self.source_schema, int)
            or self.source_schema < 1
        ):
            raise ValidationError(
                "Legacy similarity alignment source schema must be a positive integer."
            )
        object.__setattr__(self, "target", self.target.strip())


def _with_legacy_similarity_alignment(
    request: "LinearDiagramRequest",
    legacy: _LegacySimilarityAlignment,
) -> "LinearDiagramRequest":
    """Attach old-schema state without exposing it through current construction."""

    if not isinstance(request, LinearDiagramRequest) or not isinstance(
        legacy,
        _LegacySimilarityAlignment,
    ):
        raise ValidationError("Legacy similarity alignment attachment is invalid.")
    object.__setattr__(request, "_legacy_similarity_alignment", legacy)
    return request


def _materialized_path(value: str | Path, *, field_name: str) -> Path:
    if not isinstance(value, (str, Path)):
        raise ValidationError(f"{field_name} must identify a materialized file.")
    raw = str(value).strip()
    if not raw or raw in {".", ".."}:
        raise ValidationError(f"{field_name} must identify a materialized file.")
    path = Path(raw)
    if not path.name:
        raise ValidationError(f"{field_name} must identify a materialized file.")
    return path


@dataclass(frozen=True)
class GenBankInputSource:
    """One materialized GenBank input file."""

    path: str | Path

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "path",
            _materialized_path(self.path, field_name="GenBank path"),
        )


@dataclass(frozen=True)
class GffFastaInputSource:
    """One materialized GFF3/FASTA input pair."""

    gff_path: str | Path
    fasta_path: str | Path

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "gff_path",
            _materialized_path(self.gff_path, field_name="GFF3 path"),
        )
        object.__setattr__(
            self,
            "fasta_path",
            _materialized_path(self.fasta_path, field_name="FASTA path"),
        )


@dataclass(frozen=True)
class FastaInputSource:
    """One materialized FASTA file given without GFF3: records with sequence and no features."""

    fasta_path: str | Path

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "fasta_path",
            _materialized_path(self.fasta_path, field_name="FASTA path"),
        )


@dataclass(frozen=True)
class InMemoryRecordSource:
    """One already-parsed sequence record.

    ``source_features``: the catalog of the record's original source, when the
    record is a planned (cropped or reversed) record of an earlier plan.
    """

    record: SeqRecord
    source_features: SourceFeatureCatalog | None = field(default=None, kw_only=True, repr=False)

    def __post_init__(self) -> None:
        if not isinstance(self.record, SeqRecord):
            raise ValidationError("In-memory record source must contain a SeqRecord.")
        if self.source_features is not None and not isinstance(self.source_features, SourceFeatureCatalog):
            raise ValidationError("source_features must be a SourceFeatureCatalog.")


RecordInputSource: TypeAlias = (
    GenBankInputSource | GffFastaInputSource | FastaInputSource | InMemoryRecordSource
)


class RecordCardinality(str, Enum):
    """How many displayed records one input source is allowed to resolve."""

    EXACTLY_ONE = "exactly_one"
    FIRST = "first"
    ALL = "all"


@dataclass(frozen=True)
class RecordPresentation:
    """Per-record presentation and placement values."""

    label: str | None = None
    subtitle: str | None = None
    reverse_complement: bool = False
    grid_row: int | None = None
    grid_column: int | None = None

    def __post_init__(self) -> None:
        for field_name in ("label", "subtitle"):
            value = getattr(self, field_name)
            if value is not None and not isinstance(value, str):
                raise ValidationError(f"Record {field_name} must be a string or None.")
            if isinstance(value, str):
                object.__setattr__(self, field_name, value.strip() or None)
        if not isinstance(self.reverse_complement, bool):
            raise ValidationError("reverse_complement must be a boolean.")
        for field_name in ("grid_row", "grid_column"):
            value = getattr(self, field_name)
            if value is not None and (
                not isinstance(value, int) or isinstance(value, bool) or value < 1
            ):
                raise ValidationError(f"{field_name} must be a positive integer or None.")
        if self.grid_column is not None and self.grid_row is None:
            raise ValidationError("grid_column requires grid_row.")


@dataclass(frozen=True)
class RecordDisplayOptions:
    """Per-record topology override and optional 1-based source display start."""

    is_circular: bool | None = None
    start_coordinate: int | None = None

    def __post_init__(self) -> None:
        if self.is_circular is not None and not isinstance(self.is_circular, bool):
            raise ValidationError("is_circular must be a boolean or None.")
        if self.start_coordinate is not None:
            if (
                not isinstance(self.start_coordinate, int)
                or isinstance(self.start_coordinate, bool)
                or self.start_coordinate < 1
            ):
                raise ValidationError("start_coordinate must be a positive integer or None.")
            if self.is_circular is False:
                raise ValidationError("An explicit display start requires a circular record.")


# The defaults of ``RecordInput``'s presentation and display: an input left at them gets a new value of its own.
_NEW_PRESENTATION = RecordPresentation()
_NEW_DISPLAY = RecordDisplayOptions()


class _NoRegion:
    """The default of ``RecordInput``'s ``region``: no region was passed, so ``windows`` stand."""


_NO_REGION = _NoRegion()


class _HeldWindows(tuple[RecordWindow, ...]):
    """The windows a ``RecordInput`` holds, which ``dataclasses.replace`` passes back with any region given to it.

    A region given with windows of this type replaces them (OV-2769); a region
    given with windows the caller built conflicts with them.
    """


@dataclass(frozen=True, init=False)
class RecordInput:
    """A record source plus selection, windows, and presentation metadata.

    ``windows`` are the spans of the record drawn, each its own unit of a
    Linear row (D-01); none draws the whole record. A record drawn as one unit
    is placed by its presentation's ``grid_column``; the windows of a record
    with several are placed by their own ``column``. ``region`` is not a
    field: a region given at construction, with its optional record selector,
    is read as one window (D-09); read back, ``region`` is that one window
    named by the record's selector, ``None`` without windows, and an error
    with several (D-20). Windows of one record may overlap (D-18).
    A ``region`` passed explicitly, ``None`` included, is the record's one span
    and replaces the windows the record held, so
    ``dataclasses.replace(record, region=None)`` draws the record whole and
    ``replace(record, region=other)`` crops it again, as before windows
    (OV-2769); without ``region``, ``replace`` copies ``windows``. A region
    given with windows of the caller's own is refused: give one or the other.
    """

    source: RecordInputSource
    selector: RecordSelector | None = None
    presentation: RecordPresentation = field(default_factory=RecordPresentation)
    record_key: str | None = None
    cardinality: RecordCardinality = RecordCardinality.EXACTLY_ONE
    display: RecordDisplayOptions = field(default_factory=RecordDisplayOptions)
    windows: Sequence[RecordWindow] = ()

    def __init__(
        self,
        source: RecordInputSource,
        selector: RecordSelector | None = None,
        region: RegionSpec | None | _NoRegion = _NO_REGION,
        presentation: RecordPresentation = _NEW_PRESENTATION,
        record_key: str | None = None,
        cardinality: RecordCardinality = RecordCardinality.EXACTLY_ONE,
        display: RecordDisplayOptions = _NEW_DISPLAY,
        windows: Sequence[RecordWindow] = (),
    ) -> None:
        # The positional order of the fields before windows (D-09): ``region`` stays third. Each input left at its
        # default gets its own presentation and display, as the fields' default factories give.
        presentation = RecordPresentation() if presentation is _NEW_PRESENTATION else presentation
        display = RecordDisplayOptions() if display is _NEW_DISPLAY else display
        # A region passed explicitly is the record's one span: it replaces the windows the record held, which
        # ``replace`` passes back (OV-2769); windows the caller built conflict with it.
        conflict = not isinstance(region, _NoRegion) and bool(windows) and not isinstance(windows, _HeldWindows)
        if isinstance(region, _NoRegion):
            region = None
        else:
            windows = ()
        for name, value in (
            ("source", source), ("selector", selector), ("presentation", presentation),
            ("record_key", record_key), ("cardinality", cardinality), ("display", display), ("windows", windows),
        ):
            object.__setattr__(self, name, value)
        if conflict:
            raise ValidationError(
                f"{self.refusal_name()}: give a region or windows, not both.", diagnostic={"code": "INPUT_INVALID"}
            )
        self._validate(region)

    @property
    def region(self) -> RegionSpec | None:
        """The one window as the region that gives it, named by the record's selector; ``None`` without windows."""

        if not self.windows:
            return None
        if len(self.windows) > 1:
            raise ValidationError(
                "This record has several windows; read RecordInput.windows instead of region.",
                diagnostic={"code": "INPUT_INVALID"},
            )
        span = self.windows[0].region()
        if self.selector is None:
            return span
        return replace(
            span,
            raw=f"{self.selector.raw}:{span.raw}",
            record_id=self.selector.record_id,
            record_index=self.selector.record_index,
        )

    def _validate(self, region: RegionSpec | None) -> None:
        if not isinstance(
            self.source,
            (GenBankInputSource, GffFastaInputSource, FastaInputSource, InMemoryRecordSource),
        ):
            raise ValidationError(
                "Record source must be GenBank, GFF3/FASTA, FASTA, or an in-memory SeqRecord."
            )
        if not isinstance(self.cardinality, RecordCardinality):
            raise ValidationError(
                "Record cardinality must be a RecordCardinality value."
            )
        if not isinstance(self.display, RecordDisplayOptions):
            raise ValidationError("Record display has an unsupported type.")
        if not isinstance(self.presentation, RecordPresentation):
            raise ValidationError("Record presentation has an unsupported type.")
        if self.record_key is not None:
            if not isinstance(self.record_key, str) or not self.record_key.strip():
                raise ValidationError("record_key must be a non-empty string or None.")
            object.__setattr__(self, "record_key", self.record_key.strip())
        if self.selector is not None:
            if not isinstance(self.selector, RecordSelector):
                raise ValidationError("Record selector has an unsupported type.")
            if (self.selector.record_id is None) == (self.selector.record_index is None):
                raise ValidationError("Record selector must identify one record ID or index.")
            if self.selector.record_index is not None and self.selector.record_index < 0:
                raise ValidationError("Record selector index must be non-negative.")
            if self.cardinality is not RecordCardinality.EXACTLY_ONE:
                raise ValidationError(
                    "A record selector identifies one record and therefore requires "
                    "RecordCardinality.EXACTLY_ONE."
                )
        if region is not None:
            if not isinstance(region, RegionSpec):
                raise ValidationError("Record region has an unsupported type.")
            if region.file_selector is not None:
                raise ValidationError(
                    "A materialized RecordInput region must not contain a file selector."
                )
            # A start after the end is a span across the origin (``:origin``); the crop needs a circular record.
            if region.start < 1 or region.end < 1:
                raise ValidationError("Record region coordinates are invalid.")
            if region.record_index is not None and region.record_index < 0:
                raise ValidationError("Record region selector index must be non-negative.")
            region_selects_record = (
                region.record_id is not None or region.record_index is not None
            )
            if self.selector is not None and region_selects_record:
                raise ValidationError(
                    "Specify the record selector either on RecordInput or in its region, not both."
                )
            if (
                region_selects_record
                and self.cardinality is not RecordCardinality.EXACTLY_ONE
            ):
                raise ValidationError(
                    "A selector-qualified region identifies one record and therefore "
                    "requires RecordCardinality.EXACTLY_ONE."
                )
            # The boundary adapter of D-09: a region is one window; a record
            # selector in it selects the record.
            if region_selects_record:
                object.__setattr__(self, "selector", region.record_selector())
            object.__setattr__(self, "windows", (RecordWindow.from_region(region),))
        self._validate_windows()

    @property
    def is_reversed(self) -> bool:
        """Whether the record is drawn reverse-complemented: by its presentation, or by its windows, all reversed."""

        return self.presentation.reverse_complement or (
            bool(self.windows) and all(window.reverse_complement for window in self.windows)
        )

    def oriented(self, reverse: bool) -> RecordInput:
        """This input drawn in the absolute direction ``reverse``: on its windows when it has them, else on the record."""

        if self.windows:
            return replace(
                self,
                windows=tuple(replace(window, reverse_complement=reverse) for window in self.windows),
                presentation=replace(self.presentation, reverse_complement=False),
            )
        return replace(self, presentation=replace(self.presentation, reverse_complement=reverse))

    def refusal_name(self) -> str:
        """How a refusal names this record (PD-OI-122): its selector, else its key, and its file."""

        source = self.source
        if isinstance(source, GenBankInputSource):
            file_name: str | None = Path(source.path).name
        elif isinstance(source, (GffFastaInputSource, FastaInputSource)):
            file_name = Path(source.fasta_path).name
        else:
            file_name = None
        name = (self.selector.raw if self.selector is not None else self.record_key) or (
            str(source.record.id) if isinstance(source, InMemoryRecordSource) else None
        )
        label = f"Record {name!r}" if name else "A record"
        return f"{label} ({file_name})" if file_name else label

    def _validate_windows(self) -> None:
        if isinstance(self.windows, (str, bytes)) or not isinstance(self.windows, Sequence):
            raise ValidationError("Record windows must be a sequence of RecordWindow values.", diagnostic={"code": "INPUT_INVALID"})
        windows = _HeldWindows(self.windows)
        if not all(isinstance(window, RecordWindow) for window in windows):
            raise ValidationError("Record windows must contain RecordWindow values.", diagnostic={"code": "INPUT_INVALID"})
        object.__setattr__(self, "windows", windows)
        if not windows:
            return
        if self.selector is None and self.cardinality is RecordCardinality.ALL:
            raise ValidationError(
                "A selectorless per-input region cannot be combined with "
                "RecordCardinality.ALL; use request record_options for "
                "collection-level regions.",
                diagnostic={"code": "INPUT_INVALID"},
            )
        # A record drawn as one unit is placed by its presentation; a record with several windows places each
        # window by its own column (D-08).
        if len(windows) == 1 and windows[0].column is not None:
            raise ValidationError("A record with one window is placed by its presentation's grid_column.", diagnostic={"code": "INPUT_INVALID"})
        if len(windows) > 1 and self.presentation.grid_column is not None:
            raise ValidationError("A record with several windows places each window by its own column.", diagnostic={"code": "INPUT_INVALID"})
        if self.presentation.reverse_complement and (
            len(windows) > 1 or windows[0].reverse_complement
        ):
            raise ValidationError(
                f"{self.refusal_name()}: reverse complement is set in both the region and record presentation."
                if len(windows) == 1
                else f"{self.refusal_name()} has several windows: reverse-complement each window, not the record.",
                diagnostic={"code": "INPUT_INVALID"},
            )


@dataclass(frozen=True)
class RecordCollectionOptions:
    """Transforms applied after all unresolved inputs have expanded."""

    regions: Sequence[RegionSpec] = ()
    labels: Sequence[str] = ()
    subtitles: Sequence[str] = ()

    def __post_init__(self) -> None:
        for field_name in ("regions", "labels", "subtitles"):
            value = getattr(self, field_name)
            if isinstance(value, (str, bytes)) or not isinstance(value, Sequence):
                raise ValidationError(
                    f"Record collection {field_name} must be a sequence."
                )
        regions = tuple(self.regions)
        if not all(isinstance(region, RegionSpec) for region in regions):
            raise ValidationError(
                "Record collection regions must contain RegionSpec values."
            )
        object.__setattr__(self, "regions", regions)
        for field_name in ("labels", "subtitles"):
            raw_values = tuple(getattr(self, field_name))
            if not all(isinstance(value, str) for value in raw_values):
                raise ValidationError(
                    f"Record collection {field_name} must contain strings."
                )
            values = tuple(value.strip() for value in raw_values)
            object.__setattr__(self, field_name, values)


def _normalized_formats(formats: Sequence[str] | str) -> tuple[str, ...]:
    if not isinstance(formats, (str, Sequence)):
        raise ValidationError("formats must be a string or sequence of strings.")
    return parse_output_formats(formats)


def _validate_interactive_metadata_policy(
    policy: object,
    formats: Sequence[str],
) -> None:
    if policy not in {"auto", "required", "omit"}:
        raise ValidationError(
            "interactive_metadata_policy must be 'auto', 'required', or 'omit'."
        )
    if policy == "required" and "interactive_svg" not in formats:
        raise ValidationError(
            "interactive metadata can be required only when interactive_svg is requested."
        )


def _contains_ascii_control(value: str) -> bool:
    return any(ord(character) < 32 or ord(character) == 127 for character in value)


# The longest name written from a prefix is the CLI Session sidecar
# "<prefix>.gbdraw-session.json", staged as ".<name>.<random>.tmp" (34 more
# bytes). 200 bytes keeps it within the 255-byte file name limit of common
# filesystems.
# A File card (RecordCardinality.ALL) that expands names its records "<card key>:<position>".
EXPANDED_RECORD_KEY_SEPARATOR = ":"


def expanded_record_key(card_key: str, position: int) -> str:
    """The key of the record at ``position`` (from 1) of the file a File card draws whole, when it draws several."""

    return f"{card_key}{EXPANDED_RECORD_KEY_SEPARATOR}{position}"


class InputPlace(NamedTuple):
    """Where one record input stands in a Linear grid before its records expand: its row, its column (a record drawn
    whole or in one window), the columns of its windows (a record drawn in several), and how many records it draws."""

    row: int | None
    column: int | None
    window_columns: tuple[int | None, ...]
    records: int

    @classmethod
    def at(cls, row: int | None, column: int | None, window_columns: Sequence[int | None], records: int) -> InputPlace:
        """The place of an input at ``row`` and ``column`` whose windows stand at ``window_columns`` (a record drawn
        in one window stands at the input's column) and which draws ``records`` records."""

        return cls(row, column, tuple(window_columns) if len(window_columns) > 1 else (), records)

    @classmethod
    def of(cls, record_input: RecordInput, records: int) -> InputPlace:
        """The place of ``record_input``, which draws ``records`` records."""

        presentation = record_input.presentation
        return cls.at(
            presentation.grid_row, presentation.grid_column, [window.column for window in record_input.windows], records
        )

    def record_column(self, index: int) -> int | None:
        """The column of the input's record ``index`` (from 0, in file order)."""

        return None if self.column is None else self.column + index

    def drawn_columns(self) -> tuple[int | None, ...]:
        """The column of each record or window the input draws, in drawing order: its windows' columns when it draws
        several, else one column per record (``record_column``)."""

        return self.window_columns or tuple(self.record_column(index) for index in range(self.records))


# A record with several windows draws each as a unit keyed "<record key>~w<window number>" (D-01); a record with one
# window draws one unit with the record's own key.
WINDOW_UNIT_KEY_SEPARATOR = "~w"


def window_unit_key(record_key: str, window_index: int, window_count: int) -> str:
    """The key of the drawn unit of one window of the record ``record_key`` (the one owner of that key)."""

    return record_key if window_count == 1 else f"{record_key}{WINDOW_UNIT_KEY_SEPARATOR}{window_index + 1}"


def drawn_record_identity(record_input: RecordInput) -> str | None:
    """The key of the record ``record_input`` draws: its own key, or, for a drawn unit that is one of several windows
    of a record (a projected request), that record's key (``WINDOW_OF_ANNOTATION``, D-01).

    Record translations and a Similarity alignment plan are keyed by it: the windows of a record move as one (D-50).
    """

    source = record_input.source
    if isinstance(source, InMemoryRecordSource):
        window_of = (source.record.annotations or {}).get(WINDOW_OF_ANNOTATION)
        if window_of is not None:
            return str(window_of)
    return record_input.record_key
_OUTPUT_PREFIX_MAX_BYTES = 200


@dataclass(frozen=True)
class RenderOutputRequest:
    """Output destination and format policy for a typed render request.

    ``embed_fonts``: the SVG and interactive SVG outputs carry a subset of each
    font they use that gbdraw does not bundle. Like ``overwrite``, it belongs
    to the caller and a Session does not store it.
    """

    output_prefix: str = "out"
    output_directory: str | Path | None = None
    formats: Sequence[str] | str = ("svg",)
    overwrite: bool = False
    interactive_metadata_policy: Literal["auto", "required", "omit"] = "auto"
    resolve_prefix_from_first_record: bool = False
    embed_fonts: bool = True
    # PNG pixels per inch (dpi / 96 per SVG pixel); None is the default, 300 (DEFAULT_PNG_DPI in
    # gbdraw/render/export_limits.py). Like ``overwrite``, it is the
    # caller's: a Session's request never stores it. A Session render sets it from ``--dpi``, else from the
    # PNG DPI a Web Session saved (``read_session_png_dpi`` in gbdraw/session.py).
    dpi: float | None = None

    def __post_init__(self) -> None:
        if not isinstance(self.output_prefix, str):
            raise ValidationError("output_prefix must be a string.")
        raw_prefix = str(self.output_prefix)
        prefix = raw_prefix.strip()
        if (
            not prefix
            or Path(prefix).name != prefix
            or PureWindowsPath(prefix).name != prefix
            or prefix in {".", ".."}
            or _contains_ascii_control(raw_prefix)
            or is_windows_reserved_filename_component(prefix)
        ):
            raise ValidationError(
                "output_prefix must be a valid portable filename prefix "
                "without directories.",
                diagnostic={
                    "code": "INPUT_INVALID",
                    "field": "output_prefix",
                    "reason": "FILENAME" if prefix else "REQUIRED",
                },
            )
        if len(prefix.encode("utf-8", "surrogatepass")) > _OUTPUT_PREFIX_MAX_BYTES:
            raise ValidationError(
                f"output_prefix must be at most {_OUTPUT_PREFIX_MAX_BYTES} bytes in UTF-8.",
                diagnostic={
                    "code": "INPUT_INVALID",
                    "field": "output_prefix",
                    "reason": "FILENAME_LENGTH",
                },
            )
        object.__setattr__(self, "output_prefix", prefix)

        if self.output_directory is not None:
            if not isinstance(self.output_directory, (str, Path)):
                raise ValidationError("output_directory must be a path or None.")
            raw_directory = str(self.output_directory)
            directory = raw_directory.strip()
            if not directory:
                raise ValidationError("output_directory must not be empty.")
            if _contains_ascii_control(raw_directory):
                raise ValidationError(
                    "output_directory must not contain control characters."
                )
            object.__setattr__(self, "output_directory", Path(directory))

        normalized_formats = _normalized_formats(self.formats)
        object.__setattr__(self, "formats", normalized_formats)

        if not isinstance(self.overwrite, bool):
            raise ValidationError("overwrite must be a boolean.")
        if not isinstance(self.embed_fonts, bool):
            raise ValidationError(
                "embed_fonts must be a boolean.",
                diagnostic={"code": "INPUT_INVALID", "reason": "BOOLEAN"},
            )
        if not isinstance(self.resolve_prefix_from_first_record, bool):
            raise ValidationError(
                "resolve_prefix_from_first_record must be a boolean."
            )
        checked_png_dpi(self.dpi)
        _validate_interactive_metadata_policy(
            self.interactive_metadata_policy,
            normalized_formats,
        )


@dataclass(frozen=True)
class CircularBatchOutputPolicy:
    """Resolve one collision-free output name per displayed Circular record."""

    output_prefix: str | Path | None = None
    formats: Sequence[str] | str = ("svg",)
    overwrite: bool = False
    interactive_metadata_policy: Literal["auto", "required", "omit"] = "auto"
    embed_fonts: bool = True
    dpi: float | None = None

    def __post_init__(self) -> None:
        checked_png_dpi(self.dpi)
        if self.output_prefix is not None:
            if not isinstance(self.output_prefix, (str, Path)):
                raise ValidationError(
                    "Circular batch output_prefix must be a path or None."
                )
            raw_prefix = str(self.output_prefix)
            prefix = raw_prefix.strip()
            if (
                not prefix
                or prefix in {".", ".."}
                or _contains_ascii_control(raw_prefix)
            ):
                raise ValidationError(
                    "Circular batch output_prefix must identify a filename prefix."
                )
            path = Path(prefix)
            if not path.name:
                raise ValidationError(
                    "Circular batch output_prefix must identify a filename prefix."
                )
            object.__setattr__(self, "output_prefix", path)
        normalized_formats = _normalized_formats(self.formats)
        object.__setattr__(self, "formats", normalized_formats)
        if not isinstance(self.overwrite, bool):
            raise ValidationError("overwrite must be a boolean.")
        if not isinstance(self.embed_fonts, bool):
            raise ValidationError(
                "embed_fonts must be a boolean.",
                diagnostic={"code": "INPUT_INVALID", "reason": "BOOLEAN"},
            )
        _validate_interactive_metadata_policy(
            self.interactive_metadata_policy,
            normalized_formats,
        )


def _request_records(records: Sequence[RecordInput]) -> tuple[RecordInput, ...]:
    try:
        normalized = tuple(records)
    except TypeError as exc:
        raise ValidationError("A diagram request requires a record input sequence.") from exc
    if not normalized:
        raise ValidationError("A diagram request requires at least one record input.")
    if not all(isinstance(record, RecordInput) for record in normalized):
        raise ValidationError("Diagram request records must be RecordInput values.")
    keys = [record.record_key for record in normalized if record.record_key is not None]
    if len(set(keys)) != len(keys):
        raise ValidationError("Diagram request record keys must be unique.")
    return normalized


def _circular_request_records(records: Sequence[RecordInput]) -> tuple[RecordInput, ...]:
    normalized = _request_records(records)
    windowed = next((record for record in normalized if len(record.windows) > 1), None)
    if windowed is not None:
        # Windows of one record share a Linear row (D-01); a Circular record is drawn whole or as one span.
        raise ValidationError(
            f"{windowed.refusal_name()} has several windows; a Circular diagram draws one window of a record at "
            "most, so several windows are Linear.",
            diagnostic={"code": "INPUT_INVALID"},
        )
    return normalized


def _validate_circular_placements(
    records: Sequence[RecordInput],
    *,
    layout: CircularMultiRecordOptions | None,
) -> None:
    if any(
        record.cardinality is RecordCardinality.ALL
        for record in records
    ) and any(
        record.presentation.grid_row is not None for record in records
    ):
        raise ValidationError(
            "RecordCardinality.ALL cannot use per-input grid placement; "
            "placement is resolved only after record expansion."
        )
    placements = [record.presentation for record in records]
    has_row = any(item.grid_row is not None for item in placements)
    has_column = any(item.grid_column is not None for item in placements)
    if not has_row:
        return
    if layout is None:
        raise ValidationError("Grid placement requires a circular multi-record request.")
    if any(item.grid_row is None for item in placements):
        raise ValidationError("If one circular record has grid_row, every record must have it.")
    if has_column and any(item.grid_column is None for item in placements):
        raise ValidationError("If one circular record has grid_column, every record must have it.")
    occupied = [
        (item.grid_row, item.grid_column)
        for item in placements
        if item.grid_column is not None
    ]
    if len(set(occupied)) != len(occupied):
        raise ValidationError("Circular record grid placements must be unique.")
    if layout.multi_record_positions:
        raise ValidationError(
            "Specify circular placement in RecordPresentation or layout, not both."
        )


def _validate_linear_placements(
    records: Sequence[RecordInput],
    *,
    layout: LinearMultiRecordOptions | None,
) -> None:
    placements = [record.presentation for record in records]
    # The cells a record takes before its records expand (``InputPlace``): one per window of a record with several
    # (D-08), else its own. A File card's column is the first of the consecutive columns its records take when they
    # expand (``numbered_places``, OV-3072).
    cells = [
        (place.row, column)
        for place in (InputPlace.of(record, 1) for record in records)
        for column in place.drawn_columns()
    ]
    has_row = any(item.grid_row is not None for item in placements)
    has_column = any(column is not None for _row, column in cells)
    if not has_row:
        return
    if layout is None:
        raise ValidationError("Grid placement requires a Linear layout.")
    if any(item.grid_row is None for item in placements):
        raise ValidationError("If one Linear record has grid_row, every record must have it.")
    if has_column and any(column is None for _row, column in cells):
        raise ValidationError("If one Linear record has grid_column, every record must have it.")
    occupied = [cell for cell in cells if cell[1] is not None]
    if len(set(occupied)) != len(occupied):
        raise ValidationError("Linear record grid placements must be unique.")
    if layout.multi_record_positions:
        raise ValidationError(
            "Specify Linear placement in RecordPresentation or layout, not both."
        )


def _resolve_circular_request_options(
    options: object,
    layout: object = None,
) -> CircularDiagramOptions:
    if not isinstance(options, CircularDiagramOptions):
        raise ValidationError(
            "Circular request options must be CircularDiagramOptions."
        )
    one_circle = isinstance(layout, CircularMultiRecordOptions) and layout.record_layout == "one_circle"
    return resolve_circular_diagram_options(options, one_circle=one_circle)


@dataclass(frozen=True)
class CircularDiagramRequest:
    """Possibly unresolved record inputs for one Circular diagram.

    ``grouping`` is normalized eagerly so the request states whether the one
    diagram is a regular single-record render or a grid.  A grid may contain
    one record; separate-diagram collections use :class:`CircularBatchRequest`.
    Request planning expands cardinality and returns an exact-one materialized
    projection suitable for canonical exact-one rendering.
    """

    records: Sequence[RecordInput]
    options: CircularDiagramOptions = field(default_factory=CircularDiagramOptions)
    layout: CircularMultiRecordOptions | None = None
    output: RenderOutputRequest = field(default_factory=RenderOutputRequest)
    grouping: Literal["single", "grid"] | None = None
    record_options: RecordCollectionOptions = field(
        default_factory=RecordCollectionOptions
    )

    def __post_init__(self) -> None:
        records = _circular_request_records(self.records)
        object.__setattr__(self, "records", records)
        object.__setattr__(
            self,
            "options",
            _resolve_circular_request_options(self.options, self.layout),
        )
        if self.layout is not None and not isinstance(
            self.layout, CircularMultiRecordOptions
        ):
            raise ValidationError("Circular request layout has an unsupported type.")
        grouping = self.grouping
        if grouping is None:
            grouping = (
                "grid"
                if self.layout is not None
                or len(records) > 1
                or any(
                    record.cardinality is RecordCardinality.ALL
                    for record in records
                )
                else "single"
            )
        if grouping not in {"single", "grid"}:
            raise ValidationError("Circular request grouping must be 'single' or 'grid'.")
        if grouping == "single" and len(records) != 1:
            raise ValidationError(
                "A single Circular request requires exactly one record."
            )
        if (
            grouping == "single"
            and records[0].cardinality is RecordCardinality.ALL
        ):
            raise ValidationError(
                "A single Circular request cannot use RecordCardinality.ALL."
            )
        if grouping == "single" and self.layout is not None:
            raise ValidationError(
                "A single Circular request cannot define a grid layout."
            )
        if grouping == "grid" and self.layout is None:
            object.__setattr__(self, "layout", CircularMultiRecordOptions())
        object.__setattr__(self, "grouping", grouping)
        if not isinstance(self.output, RenderOutputRequest):
            raise ValidationError("Circular request output has an unsupported type.")
        if not isinstance(self.record_options, RecordCollectionOptions):
            raise ValidationError(
                "Circular request record_options has an unsupported type."
            )
        _validate_circular_placements(records, layout=self.layout)


@dataclass(frozen=True)
class CircularBatchRequest:
    """Possibly unresolved inputs for separate one-record Circular diagrams."""

    records: Sequence[RecordInput]
    options: CircularDiagramOptions = field(default_factory=CircularDiagramOptions)
    outputs: Sequence[RenderOutputRequest] = ()
    output_policy: CircularBatchOutputPolicy | None = None
    record_options: RecordCollectionOptions = field(
        default_factory=RecordCollectionOptions
    )

    def __post_init__(self) -> None:
        records = _circular_request_records(self.records)
        object.__setattr__(self, "records", records)
        object.__setattr__(
            self,
            "options",
            _resolve_circular_request_options(self.options),
        )
        if isinstance(self.outputs, (str, bytes)) or not isinstance(
            self.outputs,
            Sequence,
        ):
            raise ValidationError(
                "Circular batch outputs must be a sequence of RenderOutputRequest values."
            )
        outputs = tuple(self.outputs)
        if bool(outputs) == bool(self.output_policy):
            raise ValidationError(
                "Circular batch requires either explicit outputs or one output_policy."
            )
        if outputs and len(outputs) != len(records):
            raise ValidationError(
                "Circular batch requires one resolved output request per "
                "single-record input."
            )
        if not all(isinstance(output, RenderOutputRequest) for output in outputs):
            raise ValidationError(
                "Circular batch outputs must be RenderOutputRequest values."
            )
        resolved_targets = {
            (
                Path(output.output_directory or "."),
                output.output_prefix,
            )
            for output in outputs
        }
        if len(resolved_targets) != len(outputs):
            raise ValidationError(
                "Circular batch output prefixes must be unique within each directory."
            )
        if outputs and any(
            record.cardinality is RecordCardinality.ALL for record in records
        ):
            raise ValidationError(
                "Circular batch inputs using RecordCardinality.ALL require an "
                "output_policy so outputs can be resolved after expansion."
            )
        if self.output_policy is not None and not isinstance(
            self.output_policy,
            CircularBatchOutputPolicy,
        ):
            raise ValidationError(
                "Circular batch output_policy has an unsupported type."
            )
        if not isinstance(self.record_options, RecordCollectionOptions):
            raise ValidationError(
                "Circular batch record_options has an unsupported type."
            )
        if any(record.presentation.grid_row is not None for record in records):
            raise ValidationError(
                "Circular batch records cannot define grid placement."
            )
        object.__setattr__(self, "outputs", outputs)

    @property
    def grouping(self) -> Literal["batch"]:
        return "batch"


@dataclass(frozen=True)
class LinearDiagramRequest:
    """Possibly unresolved record inputs and options for a Linear render."""

    records: Sequence[RecordInput]
    options: LinearDiagramOptions = field(default_factory=LinearDiagramOptions)
    layout: LinearMultiRecordOptions | None = None
    similarity_alignment: (
        SimilarityAlignmentPlan | SimilarityAlignmentReference | None
    ) = None
    output: RenderOutputRequest = field(default_factory=RenderOutputRequest)
    record_options: RecordCollectionOptions = field(
        default_factory=RecordCollectionOptions
    )
    _legacy_similarity_alignment: _LegacySimilarityAlignment | None = field(
        default=None,
        init=False,
        repr=False,
    )

    def __post_init__(self) -> None:
        records = _request_records(self.records)
        object.__setattr__(self, "records", records)
        if not isinstance(self.options, LinearDiagramOptions):
            raise ValidationError(
                "Linear request options must be LinearDiagramOptions."
            )
        object.__setattr__(
            self,
            "options",
            resolve_linear_diagram_options(self.options),
        )
        if self.layout is not None and not isinstance(self.layout, LinearMultiRecordOptions):
            raise ValidationError("Linear request layout has an unsupported type.")
        if self.similarity_alignment is not None and not isinstance(
            self.similarity_alignment,
            (SimilarityAlignmentPlan, SimilarityAlignmentReference),
        ):
            raise ValidationError(
                "Linear request similarity_alignment has an unsupported type."
            )
        if self.similarity_alignment is not None and self.options.anchor_feature is not None:
            # Both decide where records sit along their rows.
            raise ValidationError(
                "anchor_feature cannot be combined with a similarity alignment.",
                diagnostic={"code": "INPUT_INVALID"},
            )
        if (
            isinstance(self.similarity_alignment, SimilarityAlignmentReference)
            and (
                self.options.losat_search is None
                or self.options.losat_search.program != "losatp"
                or self.options.losat_search.losatp_mode != "similarity_groups"
            )
        ):
            raise SimilarityAlignmentReferenceError(
                self.similarity_alignment,
                "requires the Similarity groups analysis that resolves it "
                "(losat_search program 'losatp' with losatp_mode 'similarity_groups').",
            )
        if not isinstance(self.output, RenderOutputRequest):
            raise ValidationError("Linear request output has an unsupported type.")
        if not isinstance(self.record_options, RecordCollectionOptions):
            raise ValidationError(
                "Linear request record_options has an unsupported type."
            )
        if any(record.presentation.grid_row is not None for record in records) and self.layout is None:
            raise ValidationError(
                "Grid placement requires LinearMultiRecordOptions; without it, grid "
                "placement is supported only by circular multi-record requests."
            )
        _validate_linear_placements(records, layout=self.layout)
        translations = (
            self.layout.record_translations if self.layout is not None else ()
        )
        plan = (
            self.similarity_alignment
            if isinstance(self.similarity_alignment, SimilarityAlignmentPlan)
            else None
        )
        if translations or plan is not None:
            # None keys are rejected just below.
            # The windows of one record share its key (D-50).
            record_keys = cast(
                "tuple[str, ...]", tuple(dict.fromkeys(drawn_record_identity(record) for record in records))
            )
            if any(record_key is None for record_key in record_keys):
                raise ValidationError(
                    "Similarity alignment and record translations require stable record keys."
                )
            expected = set(record_keys)
            actual = {translation.record_key for translation in translations}
            if translations and actual != expected:
                missing = sorted(expected - actual)
                unknown = sorted(actual - expected)
                if any(
                    key.startswith(f"{card}{EXPANDED_RECORD_KEY_SEPARATOR}")
                    for card in unknown
                    for key in record_keys
                ):
                    raise ValidationError(
                        "anchor_feature cannot align a File card that holds several records; "
                        "split the card so each record has its own position.",
                        diagnostic={
                            "code": "INPUT_INVALID",
                            "field": "anchor_feature",
                            "reason": "ANCHOR_RECORD_CARDS",
                        },
                    )
                raise ValidationError(
                    "Linear record translation coverage differs from the displayed "
                    f"records (missing={missing}, unknown={unknown})."
                )
            if plan is not None:
                if not translations:
                    raise ValidationError(
                        "A similarity alignment plan requires one base translation per record."
                    )
                plan.validate_record_coverage(record_keys)


DiagramRequest: TypeAlias = (
    CircularDiagramRequest | CircularBatchRequest | LinearDiagramRequest
)


__all__ = [
    "CircularBatchRequest",
    "CircularBatchOutputPolicy",
    "CircularDiagramRequest",
    "DiagramRequest",
    "GenBankInputSource",
    "GffFastaInputSource",
    "FastaInputSource",
    "InMemoryRecordSource",
    "LinearDiagramRequest",
    "RecordCardinality",
    "RecordCollectionOptions",
    "RecordDisplayOptions",
    "RecordInput",
    "RecordInputSource",
    "RecordPresentation",
    "RenderOutputRequest",
]
