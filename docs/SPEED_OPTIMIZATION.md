# Speed Optimization Recommendations

These recommendations focus on improving runtime while preserving identical outputs.

## Goals and constraints
- Output must not change (SVG/PNG/PDF/EPS/PS). No visual or ordering differences.
- Deterministic behavior and stable element ordering are required for reference tests.
- Avoid side effects that change floating point rounding or font metrics.

## Python (CLI + core rendering)
### High-impact, low-risk
- Share per-record feature dictionaries across legend, labels, and record groups.
- Reuse color rule matches when constructing feature objects and legends.
- Preprocess color tables and label filtering once per run and reuse.

### High-impact, medium-risk
- Compute GC/GC-skew DataFrame once per record and reuse for both tracks.
- Replace sliding-window GC/skew counting with an O(n) rolling counter (avoid slicing + count).
- Replace pandas row iteration with vectorized arrays where path generation iterates over rows.

### Medium impact
- Reduce repeated font metric calls by memoizing label text measurement per font/size.
- Avoid repeated `GbdrawConfig.from_dict()` by passing a pre-parsed config object.
- Optimize overlap resolution with interval indexing (careful to keep stable track assignment).

## Web app (Pyodide + SVG editing)
### High-impact, low-risk
- Preindex extracted features by `svg_id` and rules by qualifier/type to avoid O(n^2) scans.
- Batch DOM updates and serialize SVG once per edit cycle (debounce with requestAnimationFrame).
- Cache `getBBox()` results inside layout passes to avoid repeated layout thrashing.

### Medium impact
- Reuse extracted features returned by Python run instead of re-parsing GenBank files.
- Group DOM queries (single `querySelectorAll`) and reuse node lists across operations.
- Avoid repeated text measurement and transform parsing when legend entries are unchanged.

## Cross-cutting profiling and validation
### Profiling plan
- Python: `python -m cProfile -o prof.out ...` and review with `snakeviz` or `pstats`.
- Web: Chrome Performance panel; record generation + legend edits; watch layout thrash.

### Output stability checks
- Python: run reference SVG comparison tests; verify checksum equality.
- Web: export SVG before/after optimizations and diff with a normalized SVG tool.
- Verify label layout and legend ordering on representative genomes (small, medium, large).

## Risk notes
- Any change to text measurement or layout math can alter SVG output.
- DOM updates that change element order can affect serialization and reference diffs.
- Vectorization must preserve per-point ordering in generated paths.

## Suggested rollout
1) Apply low-risk reuse and caching in Python core.
2) Add profiling and compare against reference outputs.
3) Apply web UI batching and indexing.
4) Re-run output comparisons and manual spot checks.
