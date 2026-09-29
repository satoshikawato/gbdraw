# Issue #597 S07 result — 2026-09-28

Status: **S07 public documentation, owner-generated screenshots, recipe fixes and
verification completed locally.** Several FAIL and unmeasured items remain and
are listed below; none is reported as PASS. Production code was not changed in
S07. S08 and BUG-01 were not started.

Raw evidence root (persistent, not `/tmp`):
`/mnt/c/Users/genom/GitHub/gbdraw/.worktrees/issue597-S05-recovery-evidence-20260928/S07/`
(below: `S07/`). Case-sensitive scratch, browser profiles and pytest base
directories used the native Linux directory `/home/kawato/gbdraw-issue597-s07-scratch/`.
No work, evidence or generated artifact was stored in `/tmp`.

## Intake, refs and preservation

- Checkout: `.worktrees/issue597-S05-recovered-20260928`, branch
  `fix/issue-597-input-session-20260926`, upstream the same-named remote.
- HEAD and same-named remote at intake: `ee6b6207845cc861ae428b4890a091509e1736b6`.
  Pending `MERGE_HEAD`: `007388567222638b707fbb16fe82dbeba61551c9` (PR #634 on dev).
- Classification at intake: **193 staged / 60 unstaged / 22 untracked**
  (`S07/intake/status-porcelain.txt`, `S07/intake/fingerprint-intake.json`).
  The staged and unstaged bytes are preserved in `S07/intake/staged.patch` and
  `S07/intake/unstaged.patch`; the untracked list is `S07/intake/untracked.txt`.
  All 22 untracked paths are Issue #597 S05/S06 evidence, results, session
  prompts, or three S05/S06 Node tests.
- Actual remote `dev` was `57cef3ba47f4b7790a9145f2ce6988422a00710e` (PR #640),
  not the previously recorded `ff5b58d`. It was fetched into
  `refs/remotes/origin/dev` only; no merge, reset, switch or clean was run.
  The advance is the Issue #619 boundary follow-up. Its runtime change to
  `services/config.js` restores the saved active mode when loading current Web
  Sessions; the S07 pages therefore do not describe which mode Load opens. The
  Product Contract (revision 26, PD-OI-044/045 unchanged) is byte-identical to
  that `dev`. Guard and authority files (`.github/workflows/test.yml`,
  `web-base-policy.yml`, `tools/web-*.json`, ratchet policies, checker) equal
  that `dev`.
- The shared checkout and Issue #619 work were not edited. S05/S06/S06
  follow-up changes were preserved without rewriting.

## Authority

PD-OI-044 (`diagram-generation.circular-transform-discoverability`) and
PD-OI-045 (`web.session-operation-consistency`), Contract revision 26, are the
active authority. S07 documents existing implemented behavior and corrects
public text that promised more than the implementation delivers. Product
Impact classification: `IMPLEMENT_EXISTING_AUTHORITY` (documentation corrected
to authority and measured behavior; no outcome chosen). Architecture: not
architecture-bearing for runtime; the documentation capture harness gained one
shared helper that replaces duplicated inline code (`open_ancestor_details`)
and one status assertion helper. OE/PE/CB unchanged.

## Public-page decision gate

Recorded before writing; `UNCONFIRMED` by a human because the session is
non-interactive.

| Reader question | Existing owner | Evidence | Disposition | Resulting owner |
| --- | --- | --- | --- | --- |
| How do I know which records a Circular upload contains, and what if inspection fails? | `docs/REFERENCE/web-app.md` (one sentence) | `S07/observe-discovery`, `S07/observe-followup` | merge | `web-app.md` § Circular source records and one-record settings |
| When can I crop, reverse or retitle one Circular record, and why is it unavailable? | none | `S07/observe-discovery`, `S07/observe-record-choice`, Node `circular-record-presentation` | merge | same section |
| What happens while a Session saves or loads; limits; failures | `web-app.md` § Main workflow (two sentences) | `S07/observe-sessions`, `S07/observe-errors-*`, `S07/observe-followup` | merge | `web-app.md` § Save and Load Sessions |
| Why does a loaded Circular Session say Records not inspected? | none | `S07/observe-sessions/observation-small.json` | merge | FAQ entry linking to `web-app.md` |
| Why are controls unavailable during Save/Load? | none | `S07/observe-sessions/observation-browse.json` | merge | FAQ entry linking to `web-app.md` |
| Rotation controls after Load | `web-app.md` § Rotate a record (stale button name) | `S07/observe-cli-load`, `S07/probe-rotation-button-*` | keep page, correct text | `web-app.md` |
| Does Load start the diagram engine? | `session-and-request-compatibility.md` (over-broad promise) | `S07/observe-cli-load`, `S07/observe-sessions/observation-cli.json` | keep page, correct text | links to `web-app.md` |
| First Circular Tutorial Step 1 | `first-circular-genome-diagram.md` | T-GUI-01 capture | keep, update text and images | same |
| Session Tutorial Step 6 | `create-and-resume-an-interactive-figure.md` | T-GUI-09 capture | keep, update text and images | same |

No public page was added, merged away or deleted.

## Public changes and their evidence

| Documented statement | Evidence (current source) |
| --- | --- |
| Circular upload is inspected at once; status `Inspecting source records…` then `1 source record(s) inspected`; rotation row before Generate; Python diagram Worker 0 | `S07/observe-discovery/observation.json` (1440 and 390 px) |
| GFF3 alone shows `Upload both GFF3 and FASTA files to inspect records.` | Node `record-display-discovery` "incomplete GFF pair" PASS |
| Failure offers `Retry source inspection`; Replace/Remove continue; Result unchanged | `S07/observe-followup/observation.json` (1440 and 390 px) |
| Record choices: one record → `Automatic (only record)`; several → `All records (separate diagrams)` gives 2 results; one listed record gives 1 result | `S07/observe-record-choice/observation.json` |
| Multi-Record Canvas on in fresh pages and after Reset Settings; `Record` disabled with the grid reason; `Show Multi-Record Canvas setting` moves focus to the checkbox | `S07/observe-discovery`, `S07/observe-record-choice`, `S07/observe-keyboard-390` |
| One-record section opens on upload / Record choice / canvas off; manual close survives unrelated edits; no automatic record or grouping change | `S07/observe-discovery`, `S07/observe-record-choice`, Node `circular-record-presentation` 5 PASS |
| Save title prompt, `<title>.gbdraw-session.json.gz`, repeat-download confirmation, Cancel saves nothing | `S07/observe-sessions/observation-small.json` |
| Compressed Session > 50 MiB asks to continue; Cancel saves nothing | `S07/observe-sessions/observation-cli.json` (`Compressed session size is 127.3 MB. Continue?`, 0 downloads) |
| Saving/Loading status; Generate, file inputs, mode, settings, Undo/Redo, Reset Settings, other Session button unavailable; scroll, pan/zoom, feature search available; double Save click → one download | `S07/observe-sessions/observation-{small,browse}.json`; Node `session-operation-consistency` PASS |
| Save/Load unavailable during Generate with `Generating diagram. Retry after generation finishes.` | `S07/observe-followup/observation.json` |
| 200 MiB plain / 512 MiB expanded rejection; Save needs gzip compression; Load needs Workers; failure shows **Operation error** and keeps Result and Undo history | `S07/observe-errors-current/observation.json`, `S07/observe-followup/observation.json` |
| Loaded Circular Session: `Records not inspected`, Inspect lists records without changing Result or Undo; Generate inspects first | `S07/observe-sessions/observation-small.json`, `S07/observe-followup/observation.json`, Node `record-display-discovery` "saved preview Generate inspects" PASS |
| Loading a saved preview does not start LOSATP; the diagram engine starts only to validate saved settings with no web control (command-line configuration) | Web Circular/Linear Sessions: 0 diagram Workers; CLI Circular Session and the reconstructed real Linear CLI Session: 1 Worker with only `validateConfigOverrides` (`S07/observe-cli-load`, `S07/observe-sessions/observation-cli.json`) |

The heartbeat max ≤500 ms miss, native structured-clone bytes (UNAVAILABLE)
and the lost original gzip remain internal limits and are not described as
resolved in any public page.

Changed public files: `docs/REFERENCE/web-app.md`,
`docs/REFERENCE/session-and-request-compatibility.md`, `docs/FAQ.md`,
`docs/TUTORIALS/GUI/first-circular-genome-diagram.md`,
`docs/TUTORIALS/GUI/create-and-resume-an-interactive-figure.md`, and the six
images each under `docs/images/t-gui-01/` and `docs/images/t-gui-09/`.
Gallery tutorial recipes: `gbdraw/web/gallery/tutorials/tobacco-chloroplast.json`,
`vibrio-harveyi-group-collinear.json` (capture metadata only). Internal:
`docs/internal/WEB_GALLERY_OPERATION_SCREENSHOT_REGISTER.md`.

## Generators, captures and regeneration evidence

- Owner generator: `python docs/capture/run_all.py --scenario T-GUI-01 --tier core`
  and `--scenario T-GUI-09 --tier extended`, with `TMPDIR` on native Linux.
  Logs: `S07/capture/T-GUI-0{1,9}-capture-final.log`; previous committed
  images: `S07/capture/before/`; published hashes: `S07/capture/published-image-sha256.txt`.
- Capture harness corrections, each required to regenerate these two scenarios:
  1. `Layout` is closed by default on both this branch and `dev` 57cef3b, so
     T-GUI-01 and the shared `apply_finished_human_settings` could not select
     **Track Preset** (`S07/probe-track-preset-*.json`). The shared
     `open_ancestor_details` helper now opens the containing disclosure; it
     replaces the identical inline loop in `gui_annotated_chloroplast.py`.
  2. Interactive SVG exports embed catalog schema 4 since `4d9fbfd2`; the GUI
     export validator still required schema 3. It now requires schema 4, like
     the CLI/Python recipes already did; the source-pin test was updated to the
     same value (`tests/test_gui_interactive_capture_contracts.py`).
  3. `expect_circular_source_status` makes the upload screenshots wait for
     `1 source record(s) inspected` and binds the Session Tutorial to
     `Records not inspected` after Load and the inspected status after Generate.
  4. The Generation status row reduces preview height; at 60% the bottom
     `tRNA-Asp` label sat behind the preview toolbar
     (`S07/capture/t-gui-01-04-bottom-compare.png`). Both scenarios now use
     50%; T-GUI-01 reuses `_reset_finished_preview_viewport` and
     `_frame_finished_preview_with_legend` instead of its fixed pan, so the
     whole record and the legend are verified inside the canvas.
- Visual review at readable scale: all 12 images; the finished figures keep
  every label, both GC tracks, the six-entry legend and the record metadata;
  no toolbar overlap. Input images show the new Source records status and
  one-record section.
- Reproduction: the final owner capture was followed by `--check` for both scenarios, which regenerated and matched all 12 committed images.
- Gallery: the two recipes clicked the removed Circular button or a Linear
  button that is no longer rendered after Load. Corrected recipes ran with
  `tools/capture_gallery_tutorial_screenshots.py --example … --operation …`
  and passed every declared assertion (`S07/gallery/`). The existing media
  bitmaps were kept (`Keep` in the register): the tobacco recapture truncates
  label cells and the Vibrio recapture differs only by an unrelated label.
  `examples.json`, Gallery Sessions and `artifact-manifest.json` reference the
  tutorial by path only and needed no regeneration.
- Not regenerated, with reason: T-GUI-05 (`--check` stale, all 5 images),
  T-GUI-06 and T-GUI-10 (flows fail on closed disclosures before capture),
  T-GUI-12 (`Custom Track Slots` button not found), H-GUI-16 (`Import TSV` not
  found). Their text did not change in S07; their failures reproduce the
  pre-existing closed-disclosure pattern and are an S08 entry. H-GUI images
  other than H-GUI-16 are not referenced by public pages.
- `examples/gbdraw_social_preview.png`, `tests/reference_outputs/`, `dist/`,
  `gbdraw.egg-info/`, the browser wheel and cache-bust were not changed. The
  existing ignored wheel (SHA-256 `4382ac20…`) was reused; Python sources are
  unchanged since S06 follow-up.

## Verification

| Command / check | Result |
| --- | --- |
| Python Playwright GUI observations (`S07/tools/observe_*.py`, probes) | Completed; findings below |
| Node `@playwright/test` 1.61.1 (resolved from the parent checkout) with native `TMPDIR`, `playwright.functional.config.js --workers=1`: `record-display-discovery` | **8 PASS / 1 FAIL** (F1) |
| Same, `circular-record-presentation`, `session-operation-consistency`, `session-save-lifecycle`, `session-loading-feedback`, `session-import-worker`, `settings-only-session`, `gallery-tutorial` | **45 PASS / 1 FAIL** (F4) |
| Docs/Gallery contract pytest (14 files, `-m "not slow"`) | first 207 PASS / 2 FAIL (alt-text pin, schema pin), after correction **209 PASS** |
| `pytest tests/test_output_comparison.py::TestOutputComparison` (read-only) | **16 PASS**, references unchanged |
| `ruff check gbdraw/` and `ruff check docs/capture` | PASS |
| `node tools/check-web-change-budget.mjs --base origin/dev` (57cef3b) | Gate **PASS** / Review **REQUIRED**, 0 cycles |
| `capture_gallery_tutorial_screenshots.py --example … --check` (tobacco, Vibrio) | PASS |
| `python docs/capture/run_all.py --scenario T-GUI-01 --tier core --check` and `--scenario T-GUI-09 --tier extended --check` after the final capture | **PASS** both: "all committed screenshots match a fresh capture" (`S07/capture/T-GUI-0{1,9}-check-final.log`); the independent `S07/tools/capture_to_dir.py` runs also reproduced the committed T-GUI-01 input images |

Node Playwright launched once browser profiles were on native Linux; the S06
`ILL_ILLOPN`/`SIGTRAP` launch failures were caused by profiles on `/mnt/c`
(`S07/observe-discovery-attempt1-ntfs-profile-ILL_ILLOPN.log`).

## Remaining FAIL, unmeasured items and reuse conditions

- **F1 (current branch, user-visible):** a Circular source inspection failure
  renders the whole normalized error object as raw JSON in **Source records**
  (`S07/observe-followup/discovery-error-{1440,390}.png`). Reproduce: load
  `HmmtDNA_basic_circular`, upload a file containing `invalid source`. Node
  `record-display-discovery` expects the filename there and fails. Not fixed
  (production out of S07 scope). Public text says only that the section
  "reports the failure".
- **F2 (also on `dev` 57cef3b):** Session size limits and missing browser APIs
  surface as generic `UNKNOWN` **Operation error** guidance
  (`S07/observe-errors-{current,dev-57cef3b}/`). Public text names the limits
  and the generic error without claiming a specific message.
- **F3 (also on `dev`):** after any operation error on a page without a
  Result, the Generation status shows `Invalid settings · Canonical resource
  record-1-genbank is missing.`
- **F4 (inherited test):** `session-save-lifecycle` "Save Session is
  single-flight…" passes its join/one-download assertions, then fails because a
  compression failure returns `{status: 'error'}` without the expected `error`
  object.
- **F5:** T-GUI-05/06/10/12 and H-GUI-16 captures are stale or cannot run (see
  above).
- **F6 (resolved for the published state, retained as history):** at 60%,
  T-GUI-09 `--check` twice reported raster anti-aliasing noise on preview
  images (e.g. 2,042 pixels, max channel delta 80, visually identical;
  `S07/capture/t-gui-09-06-diff.png`). After the 50% framing change the final
  capture and `--check` matched. Thresholds were not changed.
- **F7 (observation):** after editing a one-record field, programmatic
  `focus()` on **Save Session** loses focus; real Shift+Tab traversal reaches
  Save and Enter saves at 1440 and 390 px (`S07/probe-shift-tab-save`).
- **F8:** a `git merge-tree` simulation of the complete working tree against
  `dev` 57cef3b reports content conflicts in `gbdraw/web/index.html`,
  `app/run-analysis.js`, `services/config.js`, `services/error-normalization.js`,
  `tests/web/feature-selection.test.mjs` and
  `tests/web/right-drawer.playwright.spec.js` (`S07/merge-sim/`). Resolving
  them requires a dev integration merge, which this session did not perform.
- Inherited and unchanged: heartbeat max ≤500 ms FAIL (user-allowed), native
  structured-clone bytes UNAVAILABLE, original gzip lost, S05 broad pytest
  interruption, full Web Node and full pytest suites not rerun in S07.
- Reuse: S07 GUI observations hold only for the recorded runtime fingerprint
  (`S07/final/fingerprint-final.json`); any runtime or fixture change requires
  rerunning `S07/tools/observe_*.py`.

## Delivery

- `3f74d507ccd5af737987a71837a18323c9405cb2` concludes the pending merge of
  `0073885` with the inherited S05, S06 and S06 follow-up work.
- `6c8baf3e8bd50cbab883cc2b96665770e80ac173` contains the S07 changes and this
  result. Both were pushed non-force to the same-named branch; local and remote
  matched. Post-commit policy: Gate PASS / Review REQUIRED
  (`S07/gates/web-policy-postcommit.log`).
- PR [#641](https://github.com/satoshikawato/gbdraw/pull/641) targets `dev`
  and reports `CONFLICTING`/`DIRTY` against `57cef3b`, so it was not merged.
  Handoff: `S07/final/HANDOFF.md`.
- Next-session prompt:
  [S08_RESUME_AFTER_S07_20260928.md](../sessions/S08_RESUME_AFTER_S07_20260928.md).

## S08 entry

1. In a clean checkout of the pushed branch, integrate `dev` (57cef3b or
   later) and resolve F8 in the six files, preserving S04/S05 Session
   exclusivity and dev's Issue #599/#601/#619 changes; rerun the Node specs
   listed above and the architecture checker.
2. Resolve F1 (render the error summary, not the object) and F4, or record
   explicit decisions; rerun `record-display-discovery` and
   `session-save-lifecycle`.
3. Repair the closed-disclosure docs capture flows (F5) with
   `open_ancestor_details`, regenerate T-GUI-05/06/10/12 and H-GUI-16 with
   `docs/capture/run_all.py`, and run their `--check`.
4. Then complete the S08 acceptance matrix (D-01–D-04, S-01–S-05, A-01, W-01).
