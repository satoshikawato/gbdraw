# Documentation capture

This directory contains the executable browser journeys listed in `config.py`. `run_all.py` is the cross-surface evidence orchestrator: it runs those GUI journeys and directly calls the CLI and Python recipe APIs in `docs/recipes/`. It does not reimplement recipes or launch nested runner subprocesses.

For GUI scenarios, the runner starts the packaged web app on `127.0.0.1`,
uploads accession- and checksum-verified frozen test copies through visible
controls, writes every manifest-owned screenshot, and validates the downloaded
diagrams and evidence files. Those frozen copies are an offline regression
mechanism; public documentation directs readers to the authoritative sequence
databases.

Each scenario uses a fresh browser context and blocks every request that does
not target its temporary loopback server. The flows do not seed browser
storage, change biological record boundaries, or contact a remote service.
T-GUI-08 and H-GUI-08 upload all five Hepatoplasmataceae records and execute
their LOSATP Collinear searches from an empty cache. They select **Infer
orthogroups with self-comparisons** and clear **Max target seqs**, matching the
Gallery request and the CLI/Python recipes. Fresh Web Collinear settings cap raw
hits at `5` and leave inference off; without inference the runs execute 8 or 20
searches instead of 13 or 25 and draw 462 instead of 500 Collinear match
elements. H-GUI-14 reloads only the
current-format session downloaded earlier in the same raw-input journey, and
it does so in a fresh context. Multi-record Circular examples contain
unchanged complete natural records. Lambda comparisons use complete Lambda and
DE3 genomes. Protein Similarity-group examples use five complete linear BGC
records; the Collinear examples use five complete Hepatoplasmataceae genomes.
The download buttons are clicked through the UI; each file is checked for its
name, structure, biological endpoints, and static-SVG safety where applicable.

## Environment

Install the development dependencies and the pinned browser:

```bash
python -m pip install -e ".[dev]"
python -m playwright install chromium
```

The authoritative environment is:

- Python Playwright 1.61.0
- Node Playwright 1.61.1 (the paired JavaScript test environment)
- Chrome for Testing 149.0.7827.55, Playwright Chromium revision v1228
- 1440 × 900 viewport at device scale 1
- `en-US`, UTC, light color scheme, reduced motion, and vendored fonts

The runner rejects a different Python Playwright or Chromium version. The Node version is recorded here and in `config.py` so both browser-test surfaces remain explicit; this Python flow does not invoke Node.

## Regenerate and check

Regenerate the core GUI tier and every implemented CLI and Python recipe:

```bash
python docs/capture/run_all.py --tier core
```

`--scenario all` is the default. `--tier` filters only GUI captures; CLI and Python recipes are small, deterministic documentation contracts and all of them run for an `all` invocation. Use `--tier extended` or `--tier nightly` to add the corresponding GUI captures.

Regenerate one scenario from any surface:

```bash
python docs/capture/run_all.py --scenario T-GUI-01 --tier core
python docs/capture/run_all.py --scenario T-CLI-01
python docs/capture/run_all.py --scenario T-PY-01
```

Focused GUI commands keep their existing tier requirement:

```bash
python docs/capture/run_all.py --scenario T-GUI-02 --tier core
python docs/capture/run_all.py --scenario H-GUI-01 --tier core
python docs/capture/run_all.py --scenario T-GUI-03 --tier extended
python docs/capture/run_all.py --scenario H-GUI-02 --tier extended
python docs/capture/run_all.py --scenario H-GUI-03 --tier extended
python docs/capture/run_all.py --scenario H-GUI-04 --tier extended
python docs/capture/run_all.py --scenario H-GUI-05 --tier extended
python docs/capture/run_all.py --scenario H-GUI-06 --tier extended
python docs/capture/run_all.py --scenario T-GUI-04 --tier extended
python docs/capture/run_all.py --scenario H-GUI-07 --tier extended
python docs/capture/run_all.py --scenario H-GUI-08 --tier extended
python docs/capture/run_all.py --scenario H-GUI-09 --tier extended
python docs/capture/run_all.py --scenario H-GUI-10 --tier extended
python docs/capture/run_all.py --scenario H-GUI-11 --tier extended
python docs/capture/run_all.py --scenario H-GUI-12 --tier extended
python docs/capture/run_all.py --scenario H-GUI-13 --tier extended
python docs/capture/run_all.py --scenario H-GUI-14 --tier extended
python docs/capture/run_all.py --scenario H-GUI-15 --tier core
```

Regenerate every selected artifact without publishing changes, comparing GUI pixels and recipe outputs with the committed files:

```bash
python docs/capture/run_all.py --tier core --check
```

Check one scenario from any surface:

```bash
python docs/capture/run_all.py --scenario T-GUI-01 --tier core --check
python docs/capture/run_all.py --scenario T-CLI-01 --check
python docs/capture/run_all.py --scenario T-PY-01 --check
```

Focused GUI checks keep their existing tier requirement:

```bash
python docs/capture/run_all.py --scenario T-GUI-02 --tier core --check
python docs/capture/run_all.py --scenario H-GUI-01 --tier core --check
python docs/capture/run_all.py --scenario T-GUI-03 --tier extended --check
python docs/capture/run_all.py --scenario H-GUI-02 --tier extended --check
python docs/capture/run_all.py --scenario H-GUI-03 --tier extended --check
python docs/capture/run_all.py --scenario H-GUI-04 --tier extended --check
python docs/capture/run_all.py --scenario H-GUI-05 --tier extended --check
python docs/capture/run_all.py --scenario H-GUI-06 --tier extended --check
python docs/capture/run_all.py --scenario T-GUI-04 --tier extended --check
python docs/capture/run_all.py --scenario H-GUI-07 --tier extended --check
python docs/capture/run_all.py --scenario H-GUI-08 --tier extended --check
python docs/capture/run_all.py --scenario H-GUI-09 --tier extended --check
python docs/capture/run_all.py --scenario H-GUI-10 --tier extended --check
python docs/capture/run_all.py --scenario H-GUI-11 --tier extended --check
python docs/capture/run_all.py --scenario H-GUI-12 --tier extended --check
python docs/capture/run_all.py --scenario H-GUI-13 --tier extended --check
python docs/capture/run_all.py --scenario H-GUI-14 --tier extended --check
python docs/capture/run_all.py --scenario H-GUI-15 --tier core --check
```

GUI tiers are cumulative. The first GUI run can take a few minutes while the packaged Python diagram worker starts. Worker and generation waits are bounded at three minutes each. CLI and Python recipes may still be run through their standalone runners for surface-specific development; `run_all.py` is the authoritative whole-documentation entry point.

After regeneration, inspect every published Tutorial image directory at normal document width. Pixel comparison detects stale committed images. The semantic checks separately verify accessions, complete sequence lengths and topology, feature and comparison counts, track ordering and axes, popup metadata, exported evidence, legend placement, and static-SVG safety.

## Linear Result and draft checkpoints

The Web reference is the existing public owner for generation status, Auto,
Definition Lock, live editing, review directions and Reset. Its small worked
comparison is verified independently of Gallery sessions. The specialized
browser screenshot inventory in `tools/reproduce_examples_manifest.py` points
to its regeneration command; it does not add a Tutorial scenario or change
that registry's reviewed screenshot budgets:

```bash
python tools/prepare_browser_wheel.py
python docs/capture/verify_linear_live_edit.py --capture
```

The first command prepares the ignored browser wheel. The second is the
one-command reproduction of the reference image and its complete UI journey.
Omit `--capture` to verify without replacing the committed image; use `--output`
to choose an evidence directory (default `/tmp/gbdraw-linear-live-edit-evidence`).
It uses the existing pinned browser environment, fresh contexts, network-blocking
capture owner and random loopback server; no login, storage seed, or running
external service is needed. It reloads only the Session saved earlier in this
same raw-input journey.

`linear-live-edit-source-verification.json` binds both complete GenBank records
to the exact NCBI `gbwithparts` requests, UTC retrieval time, byte size and SHA-256.
DE3 uses its existing byte-identical tutorial mirror. Lambda uses the unchanged
original fetched for this capture, stored losslessly as gzip under `fixtures/`;
the receipt records both compressed and original bytes/hashes. The older tutorial
copy's taxonomy differs and is not used. Offline runs validate both pinned mirrors
and record identity before uploading the exact decompressed original bytes. The
versioned accessions do not promise a frozen annotation table: a new upstream revision must be reviewed and recorded
before replacing this evidence. Readers obtain originals through the links in
the reference page. The existing support TSV has the exact contents `CDS<TAB>gene`
plus a newline, published on that page.

The harness writes each checkpoint's current Result, mounted/exported SVG,
settings draft, committed request, feature identity and History observation.
It compares the mandatory sanitizer's output after fresh Load, retains both
source hashes in the saved resources, verifies pending Export without Generate,
and checks complete Generate → Undo → Redo restoration. Naming a new Session
adds the existing document-title History change; all other intent fields are
compared unchanged. Camera framing and capture-only search hiding change no
artifact or History. The public crop contains the current comparison,
metadata, gene labels, legend and ribbons; QA JSON/SVG/Session/log files remain
outside the repository.

Selector findings in this flow: the Auto status marker and generated SVG need
CSS locators because they have no separately named capture surface; the Feature
list search and popup label input have placeholders without associated labels.
Other operations use roles, labels, titles or existing test IDs. Imported camera
and search helpers retain their established canvas/search CSS selectors. These
are capture findings; the runtime's authority, selectors, hit checks and timeouts
are unchanged. Run this focused command in documentation CI to detect future
changes to the workflow.

## Alignment directions and Reset

`T-GUI-04` regenerates the finished five-BGC example from original inputs, then
uses the real review to reverse only its minority reference. It exercises the
combined and positions-only Reset controls, Undo between scopes, and restores
the Keep artifact before the SVG download. The Keep figure and group popup
are captured before the optional Reset checks, so a notification from those
checks does not obscure the finished comparison. Run:

```sh
python docs/capture/run_all.py --scenario T-GUI-04 --tier extended
```

`bgc-source-verification.json` records direct byte comparisons with all five
version-5 MIBiG GenBank downloads, including UTC retrieval timestamps, URLs,
sizes and source/mirror SHA-256. The BGC harness checks these pins before using
its offline copies. Public readers obtain inputs from MIBiG, as the Tutorial
specifies. No login, stored authentication or seed data is needed. `config.py`
owns the existing viewport/theme/locale; `web_server.py` selects a free loopback
port. The existing screenshot density remains unchanged.

## Preview layout, search rows, and current Result (Issue 599)

The [Web app reference](../REFERENCE/web-app.md#preview-search-and-editor) is
the public owner of this procedure. Its Session and export continuations live
in the existing [Session](../REFERENCE/session-and-request-compatibility.md#what-a-session-preserves)
and [export](../REFERENCE/output-formats-and-export.md#static-interactive-and-raster-output)
references. No public screenshot is replaced by this prose update.

Reproduce the real Gallery Circular, Linear, and two-output batch journey from
the retained checkout with the ignored, source-matched browser wheel:

```bash
python tools/prepare_browser_wheel.py
python tests/web/decoration-continuity.playwright.py --evidence-dir /tmp/issue599-s04-continuity
python tests/web/layout-affordance.playwright.py --evidence-dir /tmp/issue599-s04-affordance
```

These committed scripts load the existing Gallery sessions, use real pointer
drag on eligible items, Generate twice, change color/font/title/legend settings,
exercise Undo/Redo/Reset, save and load a Session, and download the current
Result. The second script also uses keyboard and touch for **Layout edit** and
checks Preview-only hints are absent from saved SVG, Interactive SVG, PNG, PDF,
and Session output. Both scripts write `observations.json`, real downloads, and
browser screenshots into their chosen evidence directories. The docked search
and toolbar matrix, drawer transitions, focus, zoom, short viewport, and soft
keyboard viewport are reproduced by:

```bash
GBDRAW_WEB_TEST_PORT=46211 PYTHONPATH="$PWD" npx playwright test --config=playwright.functional.config.js tests/web/preview-navigation.playwright.spec.js tests/web/right-drawer.playwright.spec.js --grep 'docked search and controls|background pan preserves|preview pan leaves|zoom controls do not overlap|mobile overlays preserve mode|compact Editor scroll recovery' --workers=1 --retries=0
```

These are loopback-only regression captures. The source Gallery sessions are
stable local fixtures for executable QA; public readers use their own input and
save the Session created from that input. Review generated figures at readable
scale before using them as publication examples.
