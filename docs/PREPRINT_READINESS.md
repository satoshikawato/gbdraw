# Preprint Readiness Checklist

This document captures the work needed to prepare a gbdraw preprint and the
accompanying codebase. It is intentionally scoped to documentation, stability,
and reproducibility tasks.

Baseline version:
- Display version: 0.9.0-beta
- Package version (PEP 440): 0.9.0b0

## 1) Stability and API polish (should-fix)

Goal: Improve library usage and testability by avoiding process exits in
non-CLI modules.

Tasks:
- Replace sys.exit() in library modules with exceptions.
  - Files: gbdraw/io/genome.py, gbdraw/io/colors.py, gbdraw/labels/filtering.py,
    gbdraw/circular.py, gbdraw/linear.py
  - Keep sys.exit() only in CLI entry points (gbdraw/cli.py).
  - Add clear exception types/messages that callers can handle.

- Decide on config.toml failure mode.
  - File: gbdraw/config/toml.py
  - Recommended: raise a clear exception (e.g., ConfigNotFoundError) when the
    packaged config is missing, and let the CLI catch and format the error.
  - Alternative: embed a built-in default config and log a warning.

Acceptance criteria:
- Library use does not terminate the process on input errors.
- CLI still exits with non-zero status on fatal errors.
- Missing config yields a clear, actionable error message.

## 2) Documentation consistency (Quickstart and Tutorials)

Goal: Ensure recent linear options are documented consistently.

Options to cover:
- --record_id
- --reverse_complement
- --region

Files to update:
- docs/QUICKSTART.md
- docs/TUTORIALS/1_Customizing_Plots.md
- docs/TUTORIALS/2_Comparative_Genomics.md
- docs/TUTORIALS/3_Advanced_Customization.md

Recommended content:
- A short "Linear advanced input selection" section with examples:
  - Selecting specific records by ID or index.
  - Reverse complement per input file.
  - Cropping regions with record selectors.

Acceptance criteria:
- Each tutorial that mentions linear mode includes at least one of the three
  options and links to the CLI reference.

## 3) Nice-to-have features (roadmap alignment)

Goal: Track future work in docs/ROADMAP.md.

Short-term:
- Tick label configuration (for example, per Mb).
- Feature-wise SVG grouping for easier post-editing.

Long-term:
- Label-overlap handling (especially for dense genomes).
- Custom tracks (read depth, motifs, etc.).
- PAF/MAF support for pairwise matches.

Acceptance criteria:
- docs/ROADMAP.md explicitly lists these items with time horizons.

## 4) Preprint-specific polish

Goal: Improve citation, reproducibility, and performance reporting.

Tasks:
- Add CITATION.cff at the repository root.
  - Include title, authors, version, and DOI (once available).
- Update docs/ABOUT.md to point to the preprint DOI (once available).
- Add a reproducibility section.
  - Include exact command lines, environment info, and reference outputs.
  - Note the exact gbdraw version used.
  - Consider a dedicated doc (for example docs/REPRODUCIBILITY.md).
- Add a small benchmarking section vs. related tools (OGDRAW/CGView/etc.).
  - Document performance on typical genome sizes.
  - Include hardware, runtime, and output size.

Acceptance criteria:
- CITATION.cff exists and is referenced from docs/ABOUT.md.
- Reproducibility and benchmarking are clearly documented with commands and
  versioned results.

## 5) Suggested sequencing

1. Implement stability/API changes.
2. Update docs and tutorials.
3. Refresh ROADMAP.md.
4. Add citation, reproducibility, and benchmarking materials.

