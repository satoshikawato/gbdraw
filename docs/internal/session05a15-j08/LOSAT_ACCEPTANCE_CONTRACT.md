# LOSAT acceptance schema contract

Newly saved Sessions must be checked against the current writer, while released
input fixtures retain their historical Session and request versions. Session,
request, raw protein cache and derived protein cache are independent namespaces;
a version change in one does not promote the others.

`tests/web/losat-session-schema-contract.test.mjs` enforces this boundary in normal
PR Web contract discovery. It compares the existing Python and Web writer
definitions, the Node LOSAT acceptance declarations, and the actual declarations
loaded from the Python fallback. Loading that fallback does not start a browser.
The test introduces no runtime schema owner or persisted representation.

When a writer changes, all four schema namespaces must agree before merge. Update only
the affected current-output expectation. Do not change the released fixture's
version, bytes or oracle, relax cache assertions, or accept multiple output
versions to make the contract pass. The existing full browser acceptance remains
responsible for migration, cancellation/error rollback, TSV bytes and geometry.
Both adapters use the same browser operations in
`tests/web/helpers/losat-cache-render-boundary.mjs` for cancellation, renderer
failure injection and render settlement. Their bodies are copied unchanged from
the accepted Node adapter; the superseded Python polling and old Worker payload
probe are removed.

The test is a top-level `tests/web/*.test.mjs` file selected by the existing PR
Web-contract job. PR smoke remains 10 tests. CI tiers, Node browser assertions,
timeouts and budgets are unchanged. The Python fallback is aligned with those
accepted assertions: migrate through Generate before Save, retain legacy
candidates, keep derived payloads in runtime memory, rebuild them on replay, and
preserve uploaded TSV bytes through the existing resource-backing reader.
Derived-reference resolution assertions inspect the live runtime cache, while
saved documents must contain no derived entries. This removes its obsolete legacy
pre-Generate Save expectation and the contradictory requirement to persist derived payloads.

## Incident and validation

PR #519 was merged at `4b7e4b62a2e994c1f07c679b5c6902c18bd11155`, tree
`1d211c001e7bb4d846854cb79875c0e5aa0395c5`. Its required LOSAT job failed at
`assertCurrentSessionBoundary`: expected Session 41, received 42. The first run
and both retries failed for the same reason. The two preceding dev LOSAT jobs
passed; this record does not classify their other CI failures as LOSAT failures.

The new contract first ran against unchanged dev code and acceptance adapters:
one comparison passed and two failed. It identified the Node Session expectation
41 versus writer 42 and the Python fallback request expectation 6 versus writer
7. The correction advances the Node output expectation to 42 and makes the
Python fallback import the existing canonical request constant.

The identical contract file, SHA-256
`d43892f2ed5df3be4ec1dc8cde1fe8a2c6f01b77c965954e5aa8f918ff8063ff`, then passed
all three comparisons. The unchanged Node browser acceptance program passed both
tests. Existing PR Web contract discovery passed 444 tests across 97 files.
Ruff and whitespace checks passed. No production file or released fixture changed.

Commands, from the repository root:

```sh
node --test tests/web/losat-session-schema-contract.test.mjs
PYTHONPATH=. python tests/run_losat_cache_browser_acceptance.py
PYTHONPATH=. python tests/run_losat_cache_browser_acceptance.py --python
```

The initial Python browser fallback stopped at its pre-Generate legacy Session
Save; unchanged dev reproduced the same failure. Its subsequent run exposed a
local installed-package mismatch (40 versus source 42). Final verification pins
Python imports to the candidate source with `PYTHONPATH`. These initial runs
are not passing fallback evidence. Final browser acceptance passes: Node 2 tests;
Python 22,567 assertions. The Python run includes strict derived-reference checks against the live cache and
uses the current resource-backing reader for the final uploaded TSV bytes.

Local logs: `/tmp/gbdraw-j08-cache-contract-evidence`. Exact-dev J08/J16 results
and GitHub gate snapshots: `/tmp/gbdraw-j08-exact-dev-evidence`. On the resulting
#519 dev, all 212 functional tests, Python 3.10–3.12 core
checks (3,756 passed and 17 skipped each), performance (13), Gallery (102),
recipes (188), and Gallery readiness passed. The LOSAT version mismatch was
the sole failing required job, so Dev staging / gate failed. The corrected
candidate is not exact-dev acceptance: that gate must pass after this correction
is reviewed, integrated and validated on its resulting dev commit. No machine
gate is human approval.
