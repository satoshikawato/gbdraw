[Home](./DOCS.md) | [Installation](./INSTALL.md) | [Quickstart](./QUICKSTART.md) | [Tutorials](./TUTORIALS/TUTORIALS.md) | [Gallery](./GALLERY.md) | [FAQ](./FAQ.md) | [ABOUT](./ABOUT.md)

# Reproducibility Guide

This document captures the minimal information required to reproduce figures,
tests, and performance measurements for the gbdraw preprint. Fill in the TODOs
with your exact environment and commands.

## Version

- Display version: 0.9.0-beta
- Package version (PEP 440): 0.9.0b0

## Environment (fill in)

- OS: Ubuntu 24.04.3 LTS (WSL2, Linux 5.15.167.4-microsoft-standard-WSL2)
- CPU: Intel(R) Core(TM) i9-14900HX (16C/32T)
- RAM: 31 GiB
- Python: 3.13.3
- Environment manager: micromamba (conda-compatible), env: base

## Installation

```bash
# Example (conda/mamba)
mamba create -n gbdraw -c conda-forge -c bioconda gbdraw
conda activate gbdraw
gbdraw --version
```

```bash
# Example (source)
git clone https://github.com/satoshikawato/gbdraw.git
cd gbdraw
python -m pip install -e ".[dev]"
gbdraw --version
```

## Command lines used in the paper (fill in)

Provide the exact commands used to generate each figure. Example template:

```bash
# Figure 1: Circular genome diagram
gbdraw circular --gbk <INPUT_GB> -o <OUTPUT_PREFIX>

# Figure 2: Linear comparison
gbdraw linear --gbk <INPUT1_GB> <INPUT2_GB> -b <BLAST_OUT> -o <OUTPUT_PREFIX>
```

If you use advanced linear options, include the exact selectors:

```bash
gbdraw linear --gbk <INPUT_GB> --record_id <RECORD_ID> --reverse_complement true \
  --region <RECORD_ID:START-END> -o <OUTPUT_PREFIX>
```

## Web app reproducibility

- Wheel filename: gbdraw-0.9.0b0-py3-none-any.whl (matches gbdraw/web/js/config.js)
- Hosted URL (if used): https://gbdraw.app/
- Local run:
  ```bash
  gbdraw gui
  ```

## Tests

```bash
# Fast tests (skip slow)
pytest tests/ -v -m "not slow"

# Full test suite
pytest tests/ -v
```

Record the test results and commit hash:
- Commit: 4a18ea9ae14b44789b99db43123d244754c3ef9f
- Command: pytest tests/ -v -m "not slow"
- Environment: Python 3.13.3, pytest 9.0.2
- Outcome: 53 passed, 2 skipped, 2 deselected (1m 29s)

## Reference outputs

If diagram output changes, update `tests/reference_outputs/` and record:
- Which tests changed: None (no reference outputs updated)
- Rationale for change: N/A
- Regeneration command(s): N/A

## Performance benchmarks (fill in)

Provide at least one typical genome size and runtime. Example template:

| Task | Input size | Hardware | Runtime | Output size |
| --- | --- | --- | --- | --- |
| Circular diagram | 13.03 Mb (NC_010162.gb) | Intel(R) Core(TM) i9-14900HX (16C/32T), 31 GiB RAM | 5.20 s | 3.31 MiB (SVG) |
| Linear comparison | 0.306 Mb + 0.287 Mb (MjeNMV.gb, MelaMJNV.gb) | Intel(R) Core(TM) i9-14900HX (16C/32T), 31 GiB RAM | 4.68 s | 1.52 MiB (SVG) |

[Home](./DOCS.md) | [Installation](./INSTALL.md) | [Quickstart](./QUICKSTART.md) | [Tutorials](./TUTORIALS/TUTORIALS.md) | [Gallery](./GALLERY.md) | [FAQ](./FAQ.md) | [ABOUT](./ABOUT.md)
