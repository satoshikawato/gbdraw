"""Sequence / genome analysis helpers (internal)."""

from .collinearity import (
    CollinearityAnchor,
    CollinearityAnchorMode,
    CollinearityBlock,
    CollinearityParameters,
    CollinearityResult,
    CollinearitySearchScope,
    LosslessCollinearityParameters,
    build_native_collinearity_blocks,
    build_orthogroup_collinearity_blocks,
    convert_collinearity_blocks_to_comparisons,
    iter_collinearity_search_pairs,
    normalize_collinearity_anchor_mode,
    normalize_collinearity_search_scope,
)
from .depth import depth_df, read_depth_tsv

__all__ = [
    "CollinearityAnchor",
    "CollinearityAnchorMode",
    "CollinearityBlock",
    "CollinearityParameters",
    "CollinearityResult",
    "CollinearitySearchScope",
    "LosslessCollinearityParameters",
    "build_native_collinearity_blocks",
    "build_orthogroup_collinearity_blocks",
    "convert_collinearity_blocks_to_comparisons",
    "iter_collinearity_search_pairs",
    "normalize_collinearity_anchor_mode",
    "normalize_collinearity_search_scope",
    "depth_df",
    "read_depth_tsv",
]


