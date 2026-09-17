# S09 integration follow-up — archived pandas dtypes

The maintainer requested push and merge into `dev` after the S09 handoff and
explanation of its retained compatibility paths. The reviewed content was
committed as `604364ad39e5e612fc2a987021a9b8368ba0eef6` and published in
[PR #537](https://github.com/satoshikawato/gbdraw/pull/537).
Its 375 source and three changed-test Git blobs match the completed S09 audit.
The earlier reports remain immutable records of the pre-publication state.

## CI finding and correction

The first completed [PR test run](https://github.com/satoshikawato/gbdraw/actions/runs/35172032566)
passed Gallery, recipes, Web contracts, Web smoke, lint and architecture checks.
Core reported 5 failures, 5,652 passes and 17 skips. The failures were the three
`test_s01_full_selector_paths_and_archived_oracle` cases and two
`test_gallery_science_matches_archived_s04_bytes` cases.

Those archived byte/digest oracles include DataFrame dtypes from pandas 2.
CI uses pandas 3.0.5, which infers `str` columns where the archived environment
inferred `object`. The two test functions now build and execute their archived
cases inside `pd.option_context("future.infer_string", False)`. The setting is
restored when each context exits. Exact bytes/digests, field order, path IDs,
numeric values, and dtype assertions remain unchanged. Other tests exercise
the dependency's default setting.

This correction changes two test files only, plus this explanatory document.
Production, dependency constraints, the benchmark canonicalizer, historical
oracles, reference outputs, and earlier evidence are unchanged. It introduces
no production owner/path, compatibility branch or Product outcome. The S09
production-source aggregate remains
`2786f88d13149749167f4f3fc5fa12c8d6bdba5e611ab06aa3c99c5d232daa6d`.

## Verification

Both complete affected modules passed all 34 tests under pandas 2.3.0 and
pandas 3.0.5 on Python 3.13.3 / NumPy 2.2.5. The pandas 3 check used a private
package directory; no shared environment was changed. Ruff and `git diff --check`
passed. The archived oracles were read, never regenerated.

The five previously failing cases also passed with CI's pandas 3.0.5,
NumPy 2.4.6 and Biopython 1.88, installed into private directories, on local
Python 3.13.3. The remote PR run separately covers Python 3.11.

```bash
python -m pytest tests/test_ortholog_path_contract.py tests/test_lossless_ortholog_paths.py -q --tb=short
PYTHONPATH=.venv/s09-pandas3:. python -m pytest tests/test_ortholog_path_contract.py tests/test_lossless_ortholog_paths.py -q --tb=short
```

The final PR checks bind this test-only correction to the actual merge head.
The scientific and browser evidence reuse limits in S09 remain in force.
