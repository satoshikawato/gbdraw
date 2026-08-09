# Documentation capture

This directory contains the implemented executable browser journeys listed in `config.py`. `run_all.py` is the cross-surface documentation orchestrator: it runs those GUI journeys and directly calls the canonical CLI and Python recipe APIs in `docs/recipes/`. It does not reimplement recipes or launch nested runner subprocesses.

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
their LOSATP Collinear searches from an empty cache. H-GUI-14 reloads only the
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

After regeneration, inspect every directory named by the manifest at normal document width. Pixel comparison detects stale committed images. The semantic checks separately verify accessions, complete sequence lengths and topology, feature and comparison counts, track ordering and axes, popup metadata, exported evidence, legend placement, and static-SVG safety.
