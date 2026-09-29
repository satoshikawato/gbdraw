# S01 reproduction and metric definitions

Measurement source: `9584ed904c7d1ee93f5a07ea4f3a2a465b061724` plus the S01 test-only recipes/probes.
Run from a new independent checkout of the same work branch. Do not reuse a shared server, browser profile, output directory, LOSAT cache, or Python installation writes.
The recipes bind their own loopback server to port 0 and create/close their own Chromium contexts. Browser network probes abort and record external requests.

## Inputs and tools

Python 3.13.3; Node 26.8.2; Chromium 149.0.7827.55; Python Playwright 1.61.0;
Biopython 1.85; psutil 7.0.0; websocket-client 1.8.0. Node `@playwright/test` was absent.
No shared environment installation was performed. Source, input, authority and tool fingerprints are in
[S01-artifacts.json](./S01-artifacts.json), [S01-namespaces-and-source.json](./S01-namespaces-and-source.json),
[S01-real-fixture.json](./S01-real-fixture.json) and [S01-authority.json](./S01-authority.json).

The biological recipe selects the first two chromosomes from each of six pinned real GenBank assemblies,
checks unique record IDs **and sequence hashes**, and writes one combined GenBank file.
These are 12 real chromosomes from six assemblies, not 12 independent species/genomes or synthetic duplicates.
The record table places two chromosomes per row in a Linear diagram. No Circular multi-file collection or new bindings format is introduced.
LOSAT 0.1.0 is prepared using the repository release lock in the recipe's **subprocess-local** output cache.
Its executable SHA must equal the release lock. All directed pairs are actually searched; non-empty render comparisons alone are not proof of a full cache.
The recipe asserts 12 records, at least 25,000 biological catalog features, the exact 132 non-self directed cache pairs plus 12 self pairs (144 raw entries), writer 44/request 8, and existing 200/512 MiB limits. `--verify-only` can recheck saved outputs without repeating LOSAT searches.
Timestamp and output-directory metadata make regenerated gzip bytes non-identical; verify biological source/sequence fingerprints, counts, pair coverage and canonical semantics instead of claiming deterministic gzip bytes.

## Commands

Use an environment that already has the above dependencies, or install only into an independently owned environment.
Preparation/measurement commands write outputs to the chosen `/tmp` directory, never Gallery/reference assets.
The wheel is a generated, ignored local test asset; do not stage it or refresh cache bust for this measurement.

```bash
ISSUE597_S01_OUTPUT=$(mktemp -d /tmp/issue597-S01-reproduce.XXXXXX)
python tools/prepare_browser_wheel.py --no-build-isolation
python tools/inspect_issue597_s01_contracts.py > "$ISSUE597_S01_OUTPUT/authority.json"
python tools/inspect_issue597_s01_namespaces.py > "$ISSUE597_S01_OUTPUT/namespaces.json"
python tools/characterize_issue597_s01.py --output "$ISSUE597_S01_OUTPUT/discovery.json"
python tools/prepare_issue597_s01_fixture.py --output "$ISSUE597_S01_OUTPUT/fixture" --threads 8
python tools/measure_issue597_s01.py --mode pipeline --repetitions 3 \
  --fixture "$ISSUE597_S01_OUTPUT/fixture/real-full-pairwise.gbdraw-session.json.gz" \
  --output "$ISSUE597_S01_OUTPUT/pipeline"
python tools/measure_issue597_s01.py --mode transport --repetitions 3 \
  --fixture "$ISSUE597_S01_OUTPUT/fixture/real-full-pairwise.gbdraw-session.json.gz" \
  --output "$ISSUE597_S01_OUTPUT/transport"
python tools/verify_issue597_s01_outputs.py \
  --source "$ISSUE597_S01_OUTPUT/fixture/real-full-pairwise.gbdraw-session.json.gz" \
  --saved "$ISSUE597_S01_OUTPUT/pipeline/real-0-saved.gbdraw-session.json.gz" \
  --output "$ISSUE597_S01_OUTPUT/replay-real" --replay --fresh-load
python tools/verify_issue597_s01_outputs.py \
  --source gbdraw/web/gallery/sessions/vibrio-harveyi-group-collinear.gbdraw-session.json.gz \
  --saved "$ISSUE597_S01_OUTPUT/pipeline/vibrio-0-saved.gbdraw-session.json.gz" \
  --output "$ISSUE597_S01_OUTPUT/replay-vibrio" --replay --fresh-load
python docs/internal/issue-597-input-session-implementation-20260926/authority-candidates/verify_candidate.py privileged "$PWD" --inert
```

Run generation, pipeline and transport **sequentially** for performance evidence.
Each measurement mode always includes the existing Vibrio, plus the additional real fixture if supplied.
Repeat the semantic verifier for saved repetitions 1 and 2. The optional replay checks the existing Python reader/materializer/request boundary and CLI SVG conversion, separately from browser Save SVG equivalence.
A failed recipe asserts/exits nonzero. Measurement exceptions are preserved in the metrics JSON, not labelled PASS.

## Measurement meanings and limits

- Heartbeat: a main-thread 100 ms timer from operation start through terminal observation; p95 is sorted sample at `floor(n*0.95)` capped at `n-1`. Report count, p95 and max. A low p95 with a max above 500 ms still fails. Targets are p95 ≤250 ms **and** max ≤500 ms for every repetition and dataset.
- Long tasks: browser PerformanceObserver entries ≥50 ms in the observation interval. Main receiver callback duration excludes native structured-clone deserialization before the callback; heartbeat/long tasks include that stop.
- Pipeline wall includes automation/terminal-observer overhead. Lifecycle hooks delimit parse, preflight/sub-validation, canonical projection, restore gap, SVG admission, preview mount, save preparation/projection/compression. Preflight sub-stages are nested and must not be added to their enclosing stage.
- File read and decompression are streaming spans and byte/chunk counts, **overlapping** each other and decode. TextDecoder CPU is observed separately. Decode CPU is nested in the combined read/decode/decompress stage; these are not independent additive wall stages. Main TextEncoder encode/encodeInto CPU and output bytes report re-encoding work, including existing hashing/projection if invoked.
- Native application Load constructs no JS import Worker; any recorded diagram-generation Worker is an existing Python/helper baseline activity. Its File/stream/decoded-byte counters are observed; structured-clone reply bytes are not applicable. Native allocation/string-copy counts and internal engine copies are **unknown**, not zero.
- Transport-only: File/Blob is posted to one disposable JS Worker. The unchanged session-file codec enforces gzip magic, fatal UTF-8 and file/expanded limits; JSON.parse runs only in Worker. Replies never enter application validation/Result/History. Complete-document semantic hashes distinguish JSON types and normalise integral numeric spellings (CLI `1.0` / browser `1`), separately from byte SHA. Comparison and reserialization occur **after** the responsiveness interval.
- Whole reply: one structured-cloned plain graph; actual transfer-list bytes 0. Native wire size/copy bytes and maximum clone message size are unknown (`null`). Expanded JSON size is a logical reference, not measured wire size.
- Bounded comparison: conservative data clone budget 128 KiB, arrays ≤64 entries/message, large strings in transferred Uint16Array buffers ≤256 KiB with detach assertions and acknowledgement backpressure. Each JSON code unit, including lone surrogates, survives. Bounds/estimates count data or backing buffers and exclude protocol paths/metadata; they are not native wire-byte measurements. Main reconstructs a private plain graph, never full-document text/JSON.parse or Vue/DOM proxies.
- Transport observation adds an artificial 120 ms post-parse sampling pause and 120 ms final flush. These are in wall, excluded from Worker read/parse/transfer stage times. The two comparison modes are **not** a runtime size-based fallback. Implementation must choose only the result's selected mode, using existing known Session sections rather than introducing general RPC.
- Memory: precise performance.memory before/after and 100 ms timer samples; independent CDP Runtime.getHeapUsage on main and Worker isolates; independent 100 ms RSS sum of this Chromium process and descendants. Report each definition separately. CDP backingStorageSize observes external string/ArrayBuffer backing storage; field maxima occur at different times and cannot be summed as one simultaneous peak. RSS includes shared pages and is not PSS. Busy-main CDP timeouts are retained. Observed maxima are lower bounds; no unseen peak is certified as measured.
- Worker activity: the constructors array is cumulative within each fresh context; Save can include a Worker created by Load. Compare phase counts before claiming a new Save constructor. `ready.records` reflects Linear UI entries/files, while `ready.requestRecords` counts canonical records; one multi-record GenBank input can have one UI entry and 12 canonical records.
- Environment: i9-14900HX, 32 logical CPUs, 33,518,579,712 bytes RAM, WSL2 Linux 5.15.167.4/glibc 2.39. Only this session's own compute was kept sequential during final performance runs; other sessions were not stopped and the host was not certified idle.
- Existing Vibrio numerical budgets remain save wall <22,341 ms, observed heap delta <1,898,125,842 bytes, heartbeat max <1,000 ms or ≤2,875.8 ms. The original Node Playwright spec was not run without its dependency; numerical/semantic targeted evidence is distinct from its unexecuted UI assertions.

The real expanded document exceeds the 200 MiB plain-JSON file cap; the measured real positive input is gzip under the 512 MiB expanded cap. Real plain-JSON positive Load/performance was not measured or claimed. Auxiliary JSON/gzip transport and unchanged codec tests cover the plain branch. The optional real fresh-Load verifier currently exits 1 for the observed Python/helper Worker construction; its semantic/reader/CLI checks pass, while the Worker-free assertion remains FAIL.

Full-pipeline targets, transport feasibility, memory sampling completeness and final Product acceptance are separate statuses. A successful codec-only experiment does not resolve synchronous preflight, restore, SVG sanitation or DOM mount.
