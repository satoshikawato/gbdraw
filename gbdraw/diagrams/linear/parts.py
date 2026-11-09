"""The parts a Linear diagram is assembled from, grouped by what they describe.

The API builds them once per figure (``gbdraw.api.diagram``): the records and
how they are shown, the requested tracks, the comparisons between them, and
the alignment that moves records sideways. The configurators travel as one
``DiagramConfigurators`` beside the canvas and comparison configurators, and
the title, size and padding as one ``CompositionInputs``.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Sequence

from Bio.SeqRecord import SeqRecord

from ...analysis.depth_tracks import DepthTrackSpec
from ...analysis.feature_density import FeatureDensityTrack
from ...analysis.protein_colinearity import OrthogroupGraphResult, OrthogroupResult
from ...annotations import AnnotationOptions, ResolvedAnnotationBundle
from ...layout.record_coordinates import RecordDisplayTransform
from ...layout.record_placement import FeatureAnchorPlacement
from ...layout.similarity_alignment import SimilarityAlignmentPlan
from ...legend.table import PairwiseIdentityLegend
from ...linear_comparison import LinearComparison, ProteinComparisonFrames
from ...tracks import LinearTrackSlot

if TYPE_CHECKING:
    from gbdraw.api.options import LinearMultiRecordOptions
    from gbdraw.api.requests import RecordPresentation


@dataclass(frozen=True, slots=True)
class LinearRecordInputs:
    """The records in drawing order, how each is shown, and where the rows put them.

    ``record_transforms`` and ``record_presentations`` hold one entry per
    record when given; ``linear_layout`` holds the row positions and the
    record gap (``None``: one record per row, the default gap).
    """

    records: list[SeqRecord]
    record_transforms: Sequence[RecordDisplayTransform] | None = None
    record_presentations: Sequence[RecordPresentation] | None = None
    linear_layout: LinearMultiRecordOptions | None = None


@dataclass(frozen=True, slots=True)
class LinearTrackInputs:
    """The requested track slots, depth tracks and region annotations, before annotation planning."""

    record_depth_tracks: list[list[DepthTrackSpec]] | None = None
    # The feature density tracks: a stack the reader did not write draws each; a written stack those it names.
    feature_density_tracks: tuple[FeatureDensityTrack, ...] = ()
    linear_track_slots: list[LinearTrackSlot] | None = None
    linear_track_axis_index: int | None = None
    annotations: AnnotationOptions | ResolvedAnnotationBundle | None = None


@dataclass(frozen=True, slots=True)
class LinearComparisonInputs:
    """The comparisons between records: files, explicit and protein comparisons, orthogroups.

    ``pairwise_identity_legend`` is the identity legend the protein
    comparisons call for, used when the drawn comparisons give none.
    """

    blast_files: list[str] | None
    linear_comparisons: list[LinearComparison] | None = None
    protein_comparisons: ProteinComparisonFrames | None = None
    orthogroups: OrthogroupResult | OrthogroupGraphResult | None = None
    pairwise_identity_legend: PairwiseIdentityLegend | None = None


@dataclass(frozen=True, slots=True)
class LinearAlignmentInputs:
    """What moves records sideways: a similarity alignment, its anchor centres, a feature anchor."""

    similarity_alignment: SimilarityAlignmentPlan | None = None
    alignment_anchor_centers: Sequence[float | None] | None = None
    feature_anchor: FeatureAnchorPlacement | None = None
