# S09 integration follow-up — exact oracles across Python and pandas versions

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

## Python sum boundary

The [second PR run](https://github.com/satoshikawato/gbdraw/actions/runs/35172948564)
passed the three synthetic cases but retained two Gallery digest failures:
5,655 core tests passed and 17 were skipped. Python 3.12 changed
[floating-point `sum()`](https://docs.python.org/3/library/functions.html#sum).
The archived reports were produced on Python 3.13, whereas the PR runs 3.11.
Using a sequential sum in a temporary diagnostic reproduced both remaining CI
digests exactly; no such patch is part of production or the regression tests.

The unchanged S04 commit `dc79915a7d30abe1e1cdea9b592b330bfef45696` was
exported to an isolated directory and run independently on Python 3.10.21,
3.11.16 and 3.13.3. The current code was run on the same interpreters and inputs.
Full canonical result bytes match between S04 and current on every interpreter.
The S04 Python 3.10 and 3.11 outputs also match one another exactly. Across
Python 3.11 and 3.13, only float values differ: 161 leaves for Collinear and
140 for Orthogroup. Structure, order, identifiers, other scalar values, and
DataFrame dtypes all match. [Comparison evidence](data/s09-merge-ci/comparison.json)
records every differing leaf; the tests use no numeric tolerance.

The Gallery test now selects the independently measured S04 Python 3.11
digests for Python before 3.12, retaining the original archived digests on
3.12 and later. This adds exact environment-specific expectations without
rewriting any historical report or accepting current output as its own oracle.
Production and all prior S09 source hashes remain unchanged.

After both corrections, all 34 tests in the two affected modules passed on
actual Python 3.10.21 (pandas 2.3.3 / NumPy 2.2.6) and Python 3.11.16
(pandas 3.0.5 / NumPy 2.4.6), both with Biopython 1.88. See the
[3.10 log](data/s09-merge-ci/python310-tests.log) and
[3.11 log](data/s09-merge-ci/python311-tests.log). The 207 historical Python
source files were verified against the S04 Git blobs before accepting its
interpreter-specific digests.

[The diagnostic script](data/s09-merge-ci/gallery_oracles.py) records dependencies,
source and fixture hashes and exact semantic digests without timing or searches.
For each interpreter, run it against an export of S04 and against this tree:

```bash
git archive dc79915a7d30abe1e1cdea9b592b330bfef45696 gbdraw | tar -x -C /tmp/gbdraw-s09-s04-baseline
python docs/internal/collinear_similarity_performance_plan_2026-09-15/results/data/s09-merge-ci/gallery_oracles.py --source-root /tmp/gbdraw-s09-s04-baseline --output /tmp/baseline.json
python docs/internal/collinear_similarity_performance_plan_2026-09-15/results/data/s09-merge-ci/gallery_oracles.py --source-root . --output /tmp/current.json
```

The final PR checks bind these test-only corrections to the actual merge head.
The scientific and browser evidence reuse limits in S09 remain in force.
