# Overlap Resolution for Circular Mode

This document describes the `--resolve_overlaps` feature implemented for circular genome diagrams in gbdraw.

## Overview

The `--resolve_overlaps` option addresses the issue of overlapping features in circular genome diagrams, which is particularly useful for visualizing plasmid constructs where features often overlap.

Related: [GitHub Issue #68](https://github.com/satoshikawato/gbdraw/issues/68)

## Usage

```bash
gbdraw-circular --gbk plasmid.gb --resolve_overlaps --track_type middle
```

### Requirements

The `--resolve_overlaps` option is **only effective** when:
- `--separate_strands` is **not** used (i.e., `strandedness=False`)

### Displacement Direction by Track Type

| Track Type | Displacement Direction |
|------------|------------------------|
| `spreadout` | Outward (away from center) |
| `middle` | Outward (away from center) |
| `tuckin` | Inward (toward center) |

### Example Commands

```bash
# spreadout: overlapping features pushed outward
gbdraw-circular --gbk plasmid.gb --resolve_overlaps --track_type spreadout

# middle: overlapping features pushed outward
gbdraw-circular --gbk plasmid.gb --resolve_overlaps --track_type middle

# tuckin: overlapping features pushed inward (toward center)
gbdraw-circular --gbk plasmid.gb --resolve_overlaps --track_type tuckin

# NOT effective: separate_strands enabled (resolve_overlaps ignored)
gbdraw-circular --gbk genome.gb --resolve_overlaps --track_type middle --separate_strands
```

## Technical Implementation

### Origin-Spanning Feature Detection

The overlap detection algorithm correctly handles **origin-spanning features** (features that cross the 0/total_length boundary in circular genomes). This is critical for accurate plasmid visualization.

Example: A feature with coordinates `join(9500..10000,1..500)` in a 10,000 bp genome is correctly identified as spanning the origin.

### Track Assignment

Features are assigned to tracks based on:
1. **Occupied length** (longer features are prioritized)
2. **Start position**

The algorithm iteratively assigns features to the lowest available track number where no overlap occurs.

### Radial Positioning

When `resolve_overlaps` is enabled:
- Track 0: Default position (on the axis)
- Track 1, 2, ...: Pushed outward by `track_offset = track_id × cds_ratio × 1.2`

## Files Modified

| File | Changes |
|------|---------|
| `gbdraw/features/tracks.py` | Added origin-spanning detection to `get_feature_ends()`, `check_feature_overlap()`, and `arrange_feature_tracks()` |
| `gbdraw/features/factory.py` | Pass `genome_length` to `arrange_feature_tracks()` |
| `gbdraw/layout/circular.py` | Added `track_id` parameter to `calculate_feature_position_factors_circular()` |
| `gbdraw/svg/circular_features.py` | Added `track_id` to all path generation functions |
| `gbdraw/drawers/circular/features.py` | Pass `track_id` from `FeatureObject` to path generators |
| `gbdraw/labels/placement_circular.py` | Use `track_id` for label positioning |
| `gbdraw/drawers/circular/labels.py` | Use `track_id` for embedded label positioning |
| `gbdraw/circular.py` | Added `--resolve_overlaps` CLI option |

## API Changes

### `calculate_feature_position_factors_circular()`

```python
def calculate_feature_position_factors_circular(
    total_length: int,
    strand: str,
    track_ratio: float,
    cds_ratio: float,
    offset: float,
    track_type: str = "tuckin",
    strandedness: bool = True,
    track_id: int = 0,  # NEW: Track number for overlap resolution
) -> list[float]:
```

### `check_feature_overlap()`

```python
def check_feature_overlap(
    a: dict, 
    b: dict, 
    separate_strands: bool, 
    genome_length: Optional[int] = None  # NEW: Required for origin-spanning detection
) -> bool:
```

### `arrange_feature_tracks()`

```python
def arrange_feature_tracks(
    feature_dict: Dict[str, FeatureObject],
    separate_strands: bool,
    resolve_overlaps: bool,
    genome_length: Optional[int] = None,  # NEW: Required for origin-spanning detection
) -> Dict[str, FeatureObject]:
```

## Limitations

- Only works without `--separate_strands` (i.e., when strands are not separated)
- If both `--separate_strands` and `--resolve_overlaps` are specified, a warning is displayed and `--resolve_overlaps` is ignored
- Ticks and outer labels do not automatically adjust to the outermost feature track (may be addressed in future updates)

## Backward Compatibility

All changes are backward compatible. Existing code and command-line usage will work without modification. The `track_id` defaults to `0`, which produces the same output as before.

