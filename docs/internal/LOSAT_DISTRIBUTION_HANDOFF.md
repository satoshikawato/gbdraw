# LOSAT native distribution handoff

Base: `f536dc758bbf917ecfa46aa5489b7b01f29f83f5` (`origin/dev`, 2026-09-17).
Work branch: `setup-losat-v010`. LOSAT v0.1.0 is published at candidate
`6bfb1b09b6cb9451fa771e687c82cbb860e8c779`; the public lock is generated from
verified fixed-URL downloads. Hosted public-install acceptance passed on all
four native targets and real Rosetta.

The public pin was generated with:

```sh
python tools/lock_losat_release.py --candidate-sha 6bfb1b09b6cb9451fa771e687c82cbb860e8c779
```

The tool consumes public `RC-HANDOFF.json`, `SHA256SUMS`, native metadata and
all four native archives at fixed v0.1.0 URLs. It writes the lightweight lock
only after their candidate/size/hash identities agree. It does not reproduce
LOSAT certification. Review the generated lock, then exercise setup, offline
cache reuse, version and representative BLASTP search on all four platforms.
Before publication, the empty lock intentionally made setup fail before downloading.

After reviewing and committing the public lock, the manual
`losat-distribution.yml` workflow builds/installs a wheel on the four supported
runners and invokes `tools/verify_losat_installation.py` outside the checkout.
An additional Apple Silicon job installs x64 Python and requires its own
`sysctl.proc_translated` result to prove Rosetta execution. Native jobs require
native execution. macOS reports preserve the installed binary's extended
attributes, quarantine presence, and `codesign --display` output without
altering them. This covers setup's actual download path; it does not claim
browser-downloaded archives bypass Gatekeeper or that Developer ID signing or
notarization was performed. Hosted Rosetta evidence confirms translated x86_64
Python via sysctl.proc_translated. All three macOS setup results had no extended
attributes/quarantine. arm64 has a linker-generated ad-hoc signature; x64 is
unsigned. Developer ID signing and notarization were not performed.
The process check follows [Apple's Rosetta detection procedure](https://developer.apple.com/documentation/apple-silicon/about-the-rosetta-translation-environment).
The eight execution-mode boundary tests pass locally. The driver retains each
search's exact argv, input FASTA files, return code and consumer text-mode
stdout/stderr, plus managed/explicit result tables. Text stdout is labelled as
such; it is not a replacement for LOSAT's raw process-byte certification.
The driver requires a fresh cache and exact pinned LOSAT fixture-source SHA,
downloads through the normal setup CLI, then forbids network calls while it
checks the release BLASTP smoke and Pairwise / Similarity groups / Collinear.
The initial work-branch push registered the workflow without changing main/dev.
Subsequent runs use manual dispatch with `platforms=all` or `platforms=macos`.
Run 35211378017 passed Linux and Windows. The macOS acceptance harness required
the native read-only `xattr` command; Rosetta requires explicit x86_64 Python
launches for package installation, the verifier and the setup subprocess.
Only the three failed macOS targets were rerun; installer code was unchanged.

## Public acceptance results

- [Linux x64 and Windows x64](https://github.com/satoshikawato/gbdraw/actions/runs/35211378017)
  passed at `f28097c8dd8dcf24cc3c18a1d01f2215d593389a`. This overall run failed
  because of the three macOS harness issues described above; it is not claimed
  as an all-green campaign.
- [macOS arm64, Intel and Rosetta](https://github.com/satoshikawato/gbdraw/actions/runs/35212810672)
  passed at `f1ddef6f4e3e2af857343dfe11706ac4fe098b29`. Between the two measured
  commits, only the acceptance workflow/verifier, its tests and this document
  changed. Installed production code and the public lock were identical.
- Each of five environments performed initial public setup, then reused cache
  with network calls forbidden. Each retained 19 actual search invocations and
  their input/consumer stdout/stderr plus six managed/explicit result tables.
  The BLASTP smoke is the canonical 475-row result; Pairwise and Similarity groups
  each returned six matching rows, Collinear one. The runtime source is LOSAT
  `6bfb1b09b6cb9451fa771e687c82cbb860e8c779` on every target.
- Published lock SHA-256:
  `72da7f7ac6371788e862093686a49a9abccf3142fc6f8f22fb78ce8e92cfde04`.
  Final wheel/sdist inspection confirms this lock is packaged and native
  executables/user caches are absent. Final lock read from a clean installed
  wheel passed; no Python package was uploaded.

The raw results, both full workflow logs, failure history and independent review
are retained in the LOSAT workspace under
`artifacts/v010-distribution-20260917/gbdraw-public-acceptance-complete.tar.gz`.
Subsequent documentation-only handoff edits do not relabel these measured commits.

Architecture: new installation/cache capability has one owner
(`gbdraw/losat_setup.py`) and one explicit network entry (`setup_losat`);
searches use only `managed_losat`. CLI dispatch and existing protein resolver
are adapters. No alternate search implementation or persisted-format migrator
is introduced (OE/PE/CB remain zero for this new capability). Existing bundled
source/PATH and explicit NCBI policies are retained; native package payloads are
removed from the one package-data owner. The release lock is generated data;
Command-line reference owns user setup/error/cache behavior (`keep` disposition).

Native argv uses current LOSAT CLI v2. The semantic search-cache keys retain
their existing spelling for saved-session/Web interoperability; no cache
migration or browser argument change is needed.

Local validation (Linux x64, Python 3.13.3): 261 focused tests passed, one
pre-existing test skipped; 256 CLI/collinearity integration tests passed.
Ruff passed for production code and the changed tests/tools. Wheel and sdist
build successfully and contain the release lock without native binaries/cache.
An installed wheel used a fixture-transport archive containing the real Rust
1.92 candidate binary (`0dd6c2b83f437a940d5fcbea5ca8a247934b48f3`): setup,
offline reuse and all three protein comparison modes passed; managed/explicit
outputs agree. The release BLASTP fixture produced 475 rows and SHA-256
`fd4b010800e32ce6c823cb38b42a10b7845f3342edae892acccc8f554f9edf34`.
This is local candidate evidence, not four-platform/public-URL acceptance.
The updated acceptance driver also passed a separate installed-wheel Linux run
with fixture transport: 19 real search invocations and six comparison tables
were retained, and all saved consumer-stdout hashes were rechecked. The durable
archive is in the LOSAT workspace at
`artifacts/v010-distribution-20260917/gbdraw-acceptance-driver-recorded.tar.gz`.

Commit title: Add explicit setup and verified offline caching for native LOSAT

Summary: Install pinned release archives atomically, preserve explicit runtime
selection, verify cache identity before use, and keep native binaries out of
platform-independent Python packages. Public release pinning and hosted public-install execution are complete.
Python package publication and merging this work branch are separate actions.
