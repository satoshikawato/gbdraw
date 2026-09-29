# S05 resumed runtime reproduction

Use the existing `fix/issue-597-input-session-20260926` checkout when its tree and ownership permit.
This run reused `/tmp/gbdraw-issue597-S05.ujOWyl`; `/tmp/issue597-S05-resume-evidence` owns all
fixtures, LOSAT cache, npm cache, measurement servers, logs and generated downloads. Shared checkouts
and other running sessions were not modified. Raw outputs stay outside git; the resume JSON manifests
record their hashes. Pinned S01/S04/S05 evidence and recipes remain unchanged.

## Sources and prerequisites

Start `7209d7e2ad70ac71eaaab5409e75f2f7aab126b3`; actual remote work branch had that SHA.
Actual dev was `d313b70b9f97c2c1d70f9ae885edbead80b62021`, whose object was missing locally.
Only dev was fetched, then integrated as `c11d4854a899349cd66863c859f6f871724ba3ab`.
The runtime measurement HEAD is that integration commit plus the separately recorded runtime file hashes.
Product receipts PD-OI-044/045 and the two exact PR #620 Worker permissions were independently verified.
A late remote check found dev `252986d096011fcf1a0f5564e940480d3b92844d` (PR #621,
only revision-23 Issue #619 decision text). The unchanged S05 decisions and runtime were verified.
That verified dev is merged separately after the implementation commit into a clean tree; measurement
source remains the integration HEAD and runtime hashes above. Final integration/remote proof is at handoff.

Run the unchanged pinned readiness recipe and base checker, preserving their outputs with new names:

```bash
python docs/internal/issue-597-input-session-implementation-20260926/evidence/S05_READINESS.py
node tools/check-web-change-budget.mjs --base origin/dev
```

Readiness is narrower than runtime admission. The protected guard dependency is documented in
[S05_GUARD_DELIVERY.md](S05_GUARD_DELIVERY.md); its patch is inert and is not applied here.

## Owned fixture and measurements

Prepare the ignored local browser wheel with `python tools/prepare_browser_wheel.py --no-build-isolation`.
Install browser tooling only in this checkout/environment when absent; this run installed
`@playwright/test@1.61.0` with `--no-save --package-lock=false --ignore-scripts` and its own npm cache.
The precise environment and fingerprints are in the resume manifests.

```bash
python tools/prepare_issue597_s01_fixture.py \
  --output /tmp/issue597-S05-resume-evidence/fixture --threads 8
```

This unchanged recipe selects the same twelve distinct chromosomes from six pinned GBFFs and
computes all 132 nonself plus twelve self protein comparisons. Positive large import uses gzip.
The expanded document exceeds the 200 MiB plain file cap; it is never a plain positive case.

After fixture generation and this session's other heavy checks finish, run these **sequentially**:

```bash
python docs/internal/issue-597-input-session-implementation-20260926/evidence/S05_MEASURE.py \
  --mode transport --repetitions 3 \
  --fixture /tmp/issue597-S05-resume-evidence/fixture/real-full-pairwise.gbdraw-session.json.gz \
  --output /tmp/issue597-S05-resume-evidence/transport
python docs/internal/issue-597-input-session-implementation-20260926/evidence/S05_MEASURE.py \
  --mode pipeline --repetitions 3 \
  --fixture /tmp/issue597-S05-resume-evidence/fixture/real-full-pairwise.gbdraw-session.json.gz \
  --output /tmp/issue597-S05-resume-evidence/pipeline
python docs/internal/issue-597-input-session-implementation-20260926/evidence/S05_MEASURE.py \
  --mode limits --repetitions 1 \
  --fixture /tmp/issue597-S05-resume-evidence/fixture/real-full-pairwise.gbdraw-session.json.gz \
  --output /tmp/issue597-S05-resume-evidence/limits
python docs/internal/issue-597-input-session-implementation-20260926/evidence/S05_MEASURE.py \
  --mode journey --repetitions 1 \
  --fixture /tmp/issue597-S05-resume-evidence/fixture/real-full-pairwise.gbdraw-session.json.gz \
  --output /tmp/issue597-S05-resume-evidence/journey
```

Each invocation creates its own loopback HTTP server, Chromium browser and CDP port. It reuses S01's
100 ms heartbeat, long-task observer, precise main heap, independent 100 ms process RSS sampling,
CDP used/backing heap sampling, overlapping codec stream spans, and stage definitions. Worker codec
observation and constructor/post stacks are test-served wrappers; their hashes are recorded. They
never log contents, choose a transport, pause production, or add an ACK/flush protocol. The probe
route smoke check must confirm `timings.codecProbe` before accepted measurements; the first smoke
failed that observation check and was corrected before measurement.

`completionMs` excludes the same 200 ms observer settlement interval used by S01 pipeline. Full graph
comparison and reply serialization occur after the responsiveness interval. Native structured-clone
wire/copy bytes and native receive-task duration are null: no direct observation exists. Zero explicit
transferred-buffer bytes reflects production's absent transfer list, not zero native copying.
CDP/RSS peaks are sampled lower bounds; missing worker isolates stay unknown, and sampling errors
remain in the raw evidence. Performance runs are sequential within this session; the shared host is
not certified idle while other sessions run.

For each saved pipeline document, use the unchanged semantic verifier. Add `--replay` once per dataset:

```bash
python tools/verify_issue597_s01_outputs.py \
  --source gbdraw/web/gallery/sessions/vibrio-harveyi-group-collinear.gbdraw-session.json.gz \
  --saved /tmp/issue597-S05-resume-evidence/pipeline/vibrio-0-saved.gbdraw-session.json.gz \
  --output /tmp/issue597-S05-resume-evidence/replay-vibrio --replay
python tools/verify_issue597_s01_outputs.py \
  --source /tmp/issue597-S05-resume-evidence/fixture/real-full-pairwise.gbdraw-session.json.gz \
  --saved /tmp/issue597-S05-resume-evidence/pipeline/real-0-saved.gbdraw-session.json.gz \
  --output /tmp/issue597-S05-resume-evidence/replay-real --replay
```

Do not use the historical verifier's `--fresh-load` to assess the new runtime: that observer asserts
no Worker of any kind. S05's fresh-load trace distinguishes import JS from diagram Python and retains
**Python count zero** as its acceptance; a nonzero Python count remains a failure. Vibrio's generated
journey is independently replayed. Real Generate fails on candidate and trusted dev with
`Sequence #1: Missing GenBank file.` Empty active file bindings are preserved. Its subsequently
saved retained preview is not evidence of successful Generate.
S01 feasibility, S05 lifecycle correctness, and full-pipeline performance are separate conclusions.

The validation manifest lists exact focused and required commands, exits and failures. No test timeout,
budget, mapped contract, reference output, or authority is changed to obtain a PASS.

Additional CLI SVG diagnosis retains three strict failures. Browser-saved previews omit root
`baseProfile`; CLI adds `full`. Only that root attribute differs. An additional comparison removing
only that attribute from disposable trees passes for all remaining SVG geometry/text/style/metadata.
The unchanged verifier, test contracts and references are not edited. Reproduce the diagnosis:

```python
import sys
import xml.etree.ElementTree as ET
from pathlib import Path
sys.path.insert(0, 'tools')
from measure_issue597_s01 import document
from verify_issue597_s01_outputs import compare_svg
from tests.utils.svg_compare import compare_svgs, parse_svg

# Also repeat for real-0 and generated Vibrio paths in S05-resume-semantic.json.
saved = Path('/tmp/issue597-S05-resume-evidence/pipeline/vibrio-0-saved.gbdraw-session.json.gz')
cli = Path('/tmp/issue597-S05-resume-evidence/replay-vibrio-0/cli-replay.svg')
expected = document(saved)['results'][0]['content']
print(compare_svg(expected, cli))  # Retain strict FAIL.
a, b = parse_svg(expected), parse_svg(cli.read_text())
differences = {k: (a.attrib.get(k), b.attrib.get(k))
               for k in a.attrib.keys() | b.attrib.keys()
               if a.attrib.get(k) != b.attrib.get(k)}
assert differences == {'baseProfile': (None, 'full')}
a.attrib.pop('baseProfile', None)
b.attrib.pop('baseProfile', None)
print(compare_svgs(ET.tostring(a, encoding='unicode'), ET.tostring(b, encoding='unicode')))
```
