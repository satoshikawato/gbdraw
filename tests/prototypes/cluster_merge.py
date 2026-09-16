"""S07 test-only merge oracle, frozen from S07 baseline (S06 9e6c6e74).

Uses only unchanged endpoint ordering, compatibility and block construction.
Production must never import this exhaustive scan/copy implementation.
"""
from __future__ import annotations
from typing import Sequence
from gbdraw.analysis.collinearity import (
    CollinearityAnchor, CollinearityBlock, LosslessCollinearityParameters,
    _path_sorted_anchors, _lossless_anchors_are_compatible,
    _lossless_block_from_anchors, _final_block_sort_key,
)


def _lossless_conflicts_between_clusters(
    left: CollinearityBlock,
    right: CollinearityBlock,
    anchors: Sequence[CollinearityAnchor],
) -> int:
    left_path = _path_sorted_anchors(left.anchors, left.orientation)
    right_path = _path_sorted_anchors(right.anchors, right.orientation)
    if not left_path or not right_path:
        return 0
    left_end = left_path[-1]
    right_start = right_path[0]
    query_min = min(int(left_end.query_order), int(right_start.query_order))
    query_max = max(int(left_end.query_order), int(right_start.query_order))
    subject_min = min(int(left_end.subject_order), int(right_start.subject_order))
    subject_max = max(int(left_end.subject_order), int(right_start.subject_order))
    cluster_anchors = {*left.anchors, *right.anchors}
    return sum(
        1
        for anchor in anchors
        if anchor not in cluster_anchors
        and query_min < int(anchor.query_order) < query_max
        and subject_min < int(anchor.subject_order) < subject_max
    )


def _lossless_clusters_can_merge(
    left: CollinearityBlock,
    right: CollinearityBlock,
    *,
    anchors: Sequence[CollinearityAnchor],
    params: LosslessCollinearityParameters,
) -> bool:
    if left.kind != "cluster" or right.kind != "cluster":
        return False
    if left.orientation != right.orientation:
        return False
    left_path = _path_sorted_anchors(left.anchors, left.orientation)
    right_path = _path_sorted_anchors(right.anchors, right.orientation)
    if not left_path or not right_path:
        return False
    if not _lossless_anchors_are_compatible(
        left_path[-1],
        right_path[0],
        orientation=left.orientation,
        params=params,
    ):
        return False
    return _lossless_conflicts_between_clusters(left, right, anchors) <= int(
        params.max_conflicts
    )


def _merge_lossless_clusters(
    blocks: Sequence[CollinearityBlock],
    *,
    anchors: Sequence[CollinearityAnchor],
    params: LosslessCollinearityParameters,
) -> tuple[CollinearityBlock, ...]:
    merged: list[CollinearityBlock] = []
    singleton_blocks = [block for block in blocks if block.kind == "singleton"]
    cluster_blocks = [block for block in blocks if block.kind == "cluster"]
    for block in sorted(cluster_blocks, key=_final_block_sort_key):
        if not merged or not _lossless_clusters_can_merge(
            merged[-1],
            block,
            anchors=anchors,
            params=params,
        ):
            merged.append(block)
            continue
        previous = merged[-1]
        merged[-1] = _lossless_block_from_anchors(
            block_id=previous.block_id,
            pair=(previous.query_record_index, previous.subject_record_index),
            orientation=previous.orientation,
            anchors=(*previous.anchors, *block.anchors),
        )
    return tuple(sorted((*merged, *singleton_blocks), key=_final_block_sort_key))
