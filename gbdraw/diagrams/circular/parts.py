"""The parts one Circular record is assembled from, grouped by what they describe.

The API builds them for each record (``gbdraw.api.diagram``): the record and
its data frames, the requested track slots and region annotations, and how the
record's definition is presented. ``CircularPlannedTracks`` is what annotation
planning makes of the track inputs; the drawing reads only it, never the
requested slots.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Mapping, NamedTuple, Sequence

from Bio.SeqRecord import SeqRecord
from pandas import DataFrame

from ...analysis.conservation import DrawnConservationRing
from ...analysis.depth_tracks import DepthTrackData
from ...analysis.feature_density import FeatureDensityLayer
from ...annotations import AnnotationOptions, ResolvedAnnotationBundle, ResolvedAnnotationTrack
from ...layout.record_coordinates import RecordDisplayTransform
from ...render.groups.circular.definition import CircularDefinitionProfile
from ...tracks import CircularTrackSlot

if TYPE_CHECKING:
    from gbdraw.api.requests import RecordPresentation


@dataclass(frozen=True, slots=True)
class CircularRecordData:
    """The record and the data the API computed for it.

    ``gc_df`` is the GC frame of the request's dinucleotide;
    ``dinucleotide_dataframes`` holds one frame per dinucleotide a slot draws,
    read by its content and skew slots alike (``skew_df``).
    ``depth_track_count_value`` is the number of depth tracks the request
    declares (``None``: count ``depth_tracks``).
    """

    gb_record: SeqRecord
    gc_df: DataFrame
    depth_df: DataFrame | None = None
    depth_tracks: Sequence[DepthTrackData] | None = None
    depth_track_count_value: int | None = None
    # The rings the record draws, outermost first (``resolve_drawn_conservation_rings``).
    conservation_rings: Sequence[DrawnConservationRing] | None = None
    conservation_min_identity: float = 0.0
    dinucleotide_dataframes: dict[str, DataFrame] | None = None
    # The record's feature density layers by track id.
    feature_density: Mapping[str, FeatureDensityLayer] | None = None


@dataclass(frozen=True, slots=True)
class CircularTrackInputs:
    """The record's requested track slots and region annotations, before annotation planning.

    ``annotation_record_index`` is the record's index on its canvas: its
    record key, its annotation rows and its definition group read it.
    """

    circular_track_slots: list[CircularTrackSlot] | None = None
    # True when the reader gave no stack: ``circular_track_slots`` is then None or the list the
    # engine built (rings, depth tracks; a grid's records get the grid's), and region sets get
    # their slots on it as on the default list (OV-2820).
    implicit_track_stack: bool = False
    circular_track_axis_index: int | None = None
    annotations: AnnotationOptions | ResolvedAnnotationBundle | None = None
    annotation_record_index: int = 0


@dataclass(frozen=True, slots=True)
class CircularCircleUnit:
    """One record drawn on a circle it shares with other records.

    ``transform`` places it there (``RecordDisplayTransform.circle``); ``record_index`` is its index among the
    figure's records, which its record key, its feature inputs and its groups read.
    """

    record: CircularRecordData
    transform: RecordDisplayTransform
    record_index: int


@dataclass(frozen=True, slots=True)
class CircularOneCircle:
    """Several records drawn one after another on one circle (``gbdraw.layout.record_units.circle_units``).

    ``units`` are the drawn records in drawing order; ``hidden_count`` records shorter than ``min_length_bp``
    are not drawn, and the definition says how many there are. ``display_length`` is the circle's length in
    display bases, records and gaps together. ``record_band``: a band outside the tracks names the records.
    """

    units: tuple[CircularCircleUnit, ...]
    display_length: int
    hidden_count: int
    min_length_bp: int
    record_band: bool


class CircularPlannedTracks(NamedTuple):
    """What annotation planning made of one record's track inputs, and what the drawing reads of them.

    ``slots``, ``annotations`` and ``annotation_track_layouts`` are the
    planned ones; ``circular_track_axis_index`` and ``annotation_record_index``
    are the track inputs' (``CircularTrackInputs``).
    """

    slots: list[CircularTrackSlot] | None
    annotations: ResolvedAnnotationBundle
    annotation_track_layouts: dict[str, ResolvedAnnotationTrack]
    circular_track_axis_index: int | None
    annotation_record_index: int


@dataclass(frozen=True, slots=True)
class CircularPresentationInputs:
    """How the record is shown: its names, its definition, its centre and its display transform.

    ``species``, ``strain`` and ``plot_title`` are the names the definition
    draws (``None``: none given). ``layout_definition_profile`` is the
    definition the tracks are placed around when it differs from the drawn
    one; ``definition_record_count`` and ``definition_position`` keep the
    values a single-record figure draws unless a caller sets them.
    """

    species: str | None
    strain: str | None
    plot_title: str | None
    definition_position: str = "center"
    definition_profile: CircularDefinitionProfile = "full"
    layout_definition_profile: CircularDefinitionProfile | None = None
    definition_group_id: str | None = None
    definition_record_count: int = 1
    center_reserved_radius: float | None = None
    record_transform: RecordDisplayTransform | None = None
    record_presentation: RecordPresentation | None = None
