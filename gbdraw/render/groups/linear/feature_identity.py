"""Exact rendered identities shared by Linear features, labels, and matches."""

from collections import Counter
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field

from ....features.ids import compute_feature_object_hash, make_linear_rendered_feature_id
from ....layout.linear import VerticalBand
from ....svg.ids import instance_svg_id
from ...drawers.linear.features import FeatureDrawer


@dataclass(frozen=True)
class LinearFeatureDomIndex:
    """Source-bound DOM identities and final feature attachment geometry."""

    by_source_index: Mapping[tuple[int, int], str]
    by_view_id: Mapping[tuple[int, str], tuple[str, ...]]
    attachment_bands: Mapping[str, VerticalBand] = field(default_factory=dict)


def build_linear_feature_dom_index(
    feature_dicts: Sequence[Mapping[str, object]],
) -> LinearFeatureDomIndex:
    """Index the exact DOM IDs produced from the prepared feature layers."""

    record_count = len(feature_dicts)
    by_source_index: dict[tuple[int, int], str] = {}
    by_view_id: dict[tuple[int, str], tuple[str, ...]] = {}
    for record_index, feature_dict in enumerate(feature_dicts):
        features = list(feature_dict.values())
        stable_ids = [
            str(FeatureDrawer.get_feature_data_id(feature) or "")
            for feature in features
        ]
        stable_id_counts = Counter(stable_id for stable_id in stable_ids if stable_id)
        mutable_by_view_id: dict[str, list[str]] = {}
        for feature, stable_id in zip(features, stable_ids, strict=True):
            if not stable_id:
                continue
            rendered_id = make_linear_rendered_feature_id(
                record_index=record_index,
                stable_feature_id=stable_id,
                record_count=record_count,
            )
            if not rendered_id:
                continue
            source_index = getattr(feature, "source_feature_index", None)
            if stable_id_counts[stable_id] > 1 and source_index is not None:
                rendered_id = instance_svg_id(rendered_id, source_index)
            # Comparison rows name a drawn feature by its drawn-location hash
            # (`*_view_feature_svg_id`), which tells the row's frame.
            view_id = compute_feature_object_hash(feature)
            if view_id:
                mutable_by_view_id.setdefault(view_id, []).append(rendered_id)
            if source_index is None:
                continue
            key = (record_index, int(source_index))
            if key in by_source_index:
                raise ValueError(
                    "Prepared linear features contain a duplicate source feature index."
                )
            by_source_index[key] = rendered_id
        for view_id, rendered_ids in mutable_by_view_id.items():
            by_view_id[(record_index, view_id)] = tuple(rendered_ids)
    return LinearFeatureDomIndex(by_source_index, by_view_id)


__all__ = ["LinearFeatureDomIndex", "build_linear_feature_dom_index"]
