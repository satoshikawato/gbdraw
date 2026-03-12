from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Literal, Mapping, Sequence

from ...features.objects import FeatureObject
from ...features.tracks import arrange_feature_tracks
from ...tracks import CircularCustomRule, TrackSpec, iter_compiled_custom_rules


AnnularTrackKind = Literal["features", "gc_content", "gc_skew", "analysis", "custom"]


@dataclass(frozen=True)
class CircularFeatureTrackEntry:
    id: str
    kind: Literal["features", "custom"]
    show: bool
    feature_types: tuple[str, ...]
    caption: str | None
    strand_mode: Literal["all", "positive", "negative"]
    rules: tuple[CircularCustomRule, ...]
    track_spec: TrackSpec


@dataclass(frozen=True)
class CircularAnalysisTrackEntry:
    id: str
    kind: Literal["gc_content", "gc_skew", "analysis"]
    show: bool
    metric: Literal["content", "skew"]
    dinucleotide: str
    caption: str
    track_spec: TrackSpec


def _default_annular_track_specs(*, show_gc: bool, show_skew: bool) -> dict[str, TrackSpec]:
    return {
        "features": TrackSpec(id="features", kind="features", mode="circular", show=True),
        "gc_content": TrackSpec(id="gc_content", kind="gc_content", mode="circular", show=bool(show_gc)),
        "gc_skew": TrackSpec(id="gc_skew", kind="gc_skew", mode="circular", show=bool(show_skew)),
    }


def build_ordered_annular_track_specs(
    track_specs: Sequence[TrackSpec] | None,
    *,
    show_gc: bool,
    show_skew: bool,
) -> list[TrackSpec]:
    annular_specs = [
        track_spec
        for track_spec in (track_specs or [])
        if str(track_spec.kind) in {"features", "gc_content", "gc_skew", "analysis", "custom"}
    ]
    defaults = _default_annular_track_specs(show_gc=show_gc, show_skew=show_skew)
    explicit_kinds = {str(track_spec.kind) for track_spec in annular_specs}

    if not annular_specs:
        return [defaults["features"], defaults["gc_content"], defaults["gc_skew"]]

    ordered_specs: list[TrackSpec] = []
    if "features" not in explicit_kinds:
        ordered_specs.append(defaults["features"])
    ordered_specs.extend(annular_specs)
    if "gc_content" not in explicit_kinds:
        ordered_specs.append(defaults["gc_content"])
    if "gc_skew" not in explicit_kinds:
        ordered_specs.append(defaults["gc_skew"])
    return ordered_specs


def build_analysis_track_entries(
    annular_track_specs: Sequence[TrackSpec],
    *,
    default_dinucleotide: str,
) -> list[CircularAnalysisTrackEntry]:
    entries: list[CircularAnalysisTrackEntry] = []
    default_nt = str(default_dinucleotide or "GC").strip().upper() or "GC"

    for track_spec in annular_track_specs:
        kind = str(track_spec.kind)
        if kind not in {"gc_content", "gc_skew", "analysis"}:
            continue

        params = dict(track_spec.params or {})
        if kind == "gc_content":
            metric: Literal["content", "skew"] = "content"
            dinucleotide = default_nt
            caption = f"{dinucleotide} content"
        elif kind == "gc_skew":
            metric = "skew"
            dinucleotide = default_nt
            caption = f"{dinucleotide} skew"
        else:
            metric = str(params.get("metric") or "content").strip().lower()  # type: ignore[assignment]
            dinucleotide = str(params.get("dinucleotide") or default_nt).strip().upper() or default_nt
            caption = str(params.get("caption") or "").strip()

        entries.append(
            CircularAnalysisTrackEntry(
                id=str(track_spec.id),
                kind=kind,  # type: ignore[arg-type]
                show=bool(track_spec.show),
                metric=metric,
                dinucleotide=dinucleotide,
                caption=caption,
                track_spec=track_spec,
            )
        )

    return entries


def build_feature_track_entries(
    annular_track_specs: Sequence[TrackSpec],
    *,
    fallback_feature_types: Sequence[str],
) -> list[CircularFeatureTrackEntry]:
    entries: list[CircularFeatureTrackEntry] = []
    for track_spec in annular_track_specs:
        kind = str(track_spec.kind)
        if kind not in {"features", "custom"}:
            continue
        params = dict(track_spec.params or {})
        if kind == "features":
            feature_types = tuple(
                str(feature_type)
                for feature_type in list(params.get("feature_types") or fallback_feature_types)
                if str(feature_type or "").strip()
            )
            entries.append(
                CircularFeatureTrackEntry(
                    id=str(track_spec.id),
                    kind="features",
                    show=bool(track_spec.show),
                    feature_types=feature_types,
                    caption=None,
                    strand_mode="all",
                    rules=tuple(),
                    track_spec=track_spec,
                )
            )
            continue

        entries.append(
            CircularFeatureTrackEntry(
                id=str(track_spec.id),
                kind="custom",
                show=bool(track_spec.show),
                feature_types=tuple(str(feature_type) for feature_type in list(params.get("feature_types") or [])),
                caption=str(params.get("caption") or "").strip() or None,
                strand_mode=str(params.get("strand_mode") or "all"),  # type: ignore[arg-type]
                rules=iter_compiled_custom_rules(track_spec),
                track_spec=track_spec,
            )
        )
    return entries


def _matches_rule(feature_object: FeatureObject, rule: CircularCustomRule) -> bool:
    qualifiers = dict(getattr(feature_object, "qualifiers", {}) or {})
    for key, values in qualifiers.items():
        if str(key).lower() != str(rule.qualifier).lower():
            continue
        if isinstance(values, (list, tuple, set)):
            normalized_values = [str(value) for value in values if value is not None]
        else:
            normalized_values = [str(values)]
        return any(rule.regex.search(value) for value in normalized_values)
    return False


def _matches_custom_track(entry: CircularFeatureTrackEntry, feature_object: FeatureObject) -> bool:
    if feature_object.feature_type not in set(entry.feature_types):
        return False
    feature_strand = str(getattr(feature_object, "strand", "undefined") or "undefined")
    if entry.strand_mode == "positive" and feature_strand != "positive":
        return False
    if entry.strand_mode == "negative" and feature_strand != "negative":
        return False
    return all(_matches_rule(feature_object, rule) for rule in entry.rules)


def assign_features_to_feature_tracks(
    base_feature_dict: Mapping[str, FeatureObject],
    feature_track_entries: Sequence[CircularFeatureTrackEntry],
    *,
    separate_strands: bool,
    resolve_overlaps: bool,
    split_overlaps_by_strand: bool,
    genome_length: int,
) -> dict[str, dict[str, FeatureObject]]:
    rendered_track_dicts: dict[str, dict[str, FeatureObject]] = {
        entry.id: {} for entry in feature_track_entries if entry.show
    }
    built_in_features = next(
        (entry for entry in feature_track_entries if entry.kind == "features"),
        None,
    )
    ordered_custom_tracks = [
        entry for entry in feature_track_entries if entry.kind == "custom" and entry.show
    ]

    unclaimed_ids: list[str] = []
    for feature_id, feature_object in base_feature_dict.items():
        feature_object.circular_track_id = None
        feature_object.circular_track_center_factor = 1.0
        feature_object.feature_track_id = 0
        unclaimed_ids.append(str(feature_id))

    still_unclaimed: list[str] = []
    for feature_id in unclaimed_ids:
        feature_object = base_feature_dict[feature_id]
        claimed = False
        for entry in ordered_custom_tracks:
            if not _matches_custom_track(entry, feature_object):
                continue
            feature_object.circular_track_id = entry.id
            rendered_track_dicts.setdefault(entry.id, {})[feature_id] = feature_object
            claimed = True
            break
        if not claimed:
            still_unclaimed.append(feature_id)

    if built_in_features is not None and built_in_features.show:
        built_in_feature_types = set(built_in_features.feature_types)
        rendered_track_dicts.setdefault(built_in_features.id, {})
        for feature_id in still_unclaimed:
            feature_object = base_feature_dict[feature_id]
            if feature_object.feature_type not in built_in_feature_types:
                continue
            feature_object.circular_track_id = built_in_features.id
            rendered_track_dicts[built_in_features.id][feature_id] = feature_object

    for entry in feature_track_entries:
        if not entry.show:
            continue
        feature_subset = rendered_track_dicts.get(entry.id, {})
        if not feature_subset:
            continue
        arranged_subset = arrange_feature_tracks(
            feature_subset,
            separate_strands,
            resolve_overlaps,
            split_overlaps_by_strand=split_overlaps_by_strand,
            genome_length=genome_length,
        )
        rendered_track_dicts[entry.id] = arranged_subset
    return rendered_track_dicts


def summarize_feature_objects(
    feature_dict: Mapping[str, FeatureObject] | None,
) -> tuple[list[str], set[tuple[str, str]], set[str]]:
    if not feature_dict:
        return [], set(), set()

    feature_types: list[str] = []
    seen_feature_types: set[str] = set()
    used_color_rules: set[tuple[str, str]] = set()
    default_used_features: set[str] = set()

    for feature_object in feature_dict.values():
        feature_type = str(getattr(feature_object, "feature_type", "") or "")
        if feature_type and feature_type not in seen_feature_types:
            seen_feature_types.add(feature_type)
            feature_types.append(feature_type)
        color = str(getattr(feature_object, "color", "") or "")
        caption = getattr(feature_object, "color_caption", None)
        if caption:
            used_color_rules.add((str(caption), color))
        elif feature_type:
            default_used_features.add(feature_type)
    return feature_types, used_color_rules, default_used_features


__all__ = [
    "CircularAnalysisTrackEntry",
    "CircularFeatureTrackEntry",
    "assign_features_to_feature_tracks",
    "build_analysis_track_entries",
    "build_feature_track_entries",
    "build_ordered_annular_track_specs",
    "summarize_feature_objects",
]
