# Similarity Group alignment with record-owned orientation — master plan

Status: planning committed; implementation not started.
Implementation branch: `fix/similarity-alignment-orientation-owner-20260926`
Branch base at creation: `origin/dev@22cbcca96f2ef5e20bb45fe397ba4232fb158574`, fetched 2026-09-26.
Product authority: receipts approved on 2026-09-26 and serialized as commit `61eec6c5` on
`authority/similarity-alignment-orientation-owner-20260926` (contract revision 16); merge into
`origin/dev` is pending. See [00_DECISION_PACK.md](00_DECISION_PACK.md) and Session 00.

## 1. Background and problem

gbdraw draws Linear multi-record genome diagrams. Protein comparisons between
records produce **Similarity Groups** (orthogroup-like sets of homologous
features). In the Web app, a user can pick one exact feature of a group as the
**reference** and choose **Align…**. The Python resolver
(`gbdraw/layout/similarity_alignment.py`) chooses one anchor feature in each
other record, or marks that record as skipped. The renderer then shifts each
target record horizontally so that its anchor center lines up with the
reference anchor center. A resolved plan is applied immediately. An ambiguous
plan, or the explicit **Review alignment options…** action, opens the
**Select alignment anchors** review. In the review the user can Select or Skip
a candidate per target and tick a per-target **Match reference direction**
checkbox.

Match reference direction never produces a visible effect in the Web app:

- If the target anchor lies on the strand opposite to the reference anchor,
  Apply always fails with `Alignment generation failed. Review the draft and
  retry Apply.`
- If the strands already agree, the checkbox changes nothing.

The underlying error is `ValueError: Comparison source feature index conflicts
with its view feature ID.`, raised in
`gbdraw/render/groups/linear/pairwise_match.py`. The Web hides it because the
alignment controller replaces the error that `runAnalysis` publishes with a
generic message. The template also appends a second copy of the retry sentence.

### Root cause

Record orientation has two owners:

1. the record's **Reverse complement** setting (`linearSeqs[].region_reverse`,
   projected to `presentation.reverseComplement` or `region.reverseComplement`
   in the render request); and
2. the alignment plan (`records[].orientationPolicy` and
   `records[].effectiveReverseComplement`), which Python applies late in
   `materialize_similarity_alignment_display()` (`gbdraw/api/record_planning.py`).

Before sending the render request, the browser projects protein-comparison rows
into displayed coordinates. The projection reads owner 1 through
`buildRegionSpec()` in `gbdraw/web/js/app/run-analysis.js`. Each row carries a
`view_feature_svg_id` derived from displayed coordinates. Apply
(`applyPlan()` in `gbdraw/web/js/app/similarity-alignment.js`) keeps owner 1 at
the pre-alignment orientation and puts the reversal only in owner 2. Python
reverses the record after the rows were projected, so the renderer finds rows
for the old orientation and raises. A manual Reverse works because the whole
pipeline reads owner 1. Existing tests missed the defect:

- the Python rendering tests exercise `MATCH_REFERENCE` without browser-projected
  comparison rows;
- the CLI never produces `match_reference`;
- the only browser journey disables comparisons and expects `0 reversed`.

There is a second partial path: `runAnalysisInternal()` patches
`canonicalStateOverride.linearRecordOrientations` into the serialized request
only after the comparison projection has already read the live state.

### Reproduction on the branch base

The first command below fails with the error above. The second command sets a
plain record Reverse without re-projecting and fails the same way. This proves
that the saved comparison resources are tied to the orientation at projection
time.

~~~bash
python - <<'EOF'
import copy, json, tempfile
from pathlib import Path
from gbdraw.web_support.request_render import render_embedded_canonical_web_request
session = json.loads(Path(
    "gbdraw/web/gallery/sessions/BGC0000708-BGC0000713.gbdraw-session.json"
).read_text())
for variant in ("plan-reversal", "record-reverse-without-reprojection"):
    request = copy.deepcopy(session["renderRequest"])
    if variant == "plan-reversal":
        for decision in request["layout"]["similarityAlignment"]["records"]:
            if decision["recordKey"] == "record-2":
                decision.update(rationale="user_selected",
                                orientationPolicy="match_reference",
                                effectiveReverseComplement=True)
    else:
        request["layout"]["similarityAlignment"] = None
        request["records"][1]["presentation"]["reverseComplement"] = True
    try:
        render_embedded_canonical_web_request(
            request, resources=session["resources"],
            workspace=Path(tempfile.mkdtemp()) / "render")
        print(variant, "OK")
    except Exception as error:
        print(variant, type(error).__name__, error)
EOF
~~~

In the Web app on the branch base, open the Gallery example
`BGC0000708-BGC0000713`. In the first record, click `CAG38712.1` (livA) and
choose **Review alignment options…**. Enable Match reference direction for
*Streptomyces fradiae* ATCC 10745, whose anchor neoA is on the opposite strand,
and click **Apply**.

## 2. Branches and gates

- Sessions 01–04 use exactly `fix/similarity-alignment-orientation-owner-20260926`.
  - It was created with no upstream from the base above, and the planning
    commit was pushed to `origin/fix/similarity-alignment-orientation-owner-20260926`.
  - Before every commit or push, verify the current branch, its upstream, and
    the push target. Push only to that same-named remote branch.
  - Never rebase, replace, or recreate this branch. Bring in newer `origin/dev`
    by merge.
- Session 00 is the only exception. Product authority is serialized on a
  separate authority-only branch, `authority/similarity-alignment-orientation-owner-20260926`,
  created from the then-latest `origin/dev`. That branch contains no runtime
  code, tests, or plan files.
- Sessions 02–04 change dependent behavior. They start only after the approved
  authority is merged into `origin/dev` and `origin/dev` is merged into the fix
  branch. An unmerged authority candidate or this plan is not runtime authority.
- Session 01 corrects error reporting under existing authority: `PD-OI-031`
  revision 3 already requires actionable errors. Its classification is
  `IMPLEMENT_EXISTING_AUTHORITY`, so it may run while authority is pending.
- Pushing this branch does not authorize a pull request, merge, deploy, or tag.
  Each of those needs its own authorization.

## 3. Selected behavior (acceptance semantics)

The Product Decision Owner selected this outcome on 2026-09-26. The exact
receipts are proposed in the Decision Pack.

1. **Orientation belongs to records.** Each record's Reverse complement
   setting is the only source of orientation. An alignment plan stores the
   reference, per-record anchor or Skip, and rationale only. Features, labels,
   annotations, comparison ribbons, anchor centers, and the rev indication all
   follow the record setting.
2. **Default Align is unchanged.**
   - A Python-resolved plan is applied immediately.
   - Every orientation is preserved.
   - The summary reports `0 reversed`.
3. **The review has one Match reference direction checkbox.**
   - It sits directly below the exact-reference card and is unchecked each time
     a review opens. Editing it is local and starts no Worker job.
   - It is disabled, with the reason `All selected anchors already face the
     reference direction.`, when no selected target has an opposite strand.
   - A status line under the checkbox names what Apply will do:
     - checked: `Apply reverses N record(s): <record labels>.`
     - unchecked: `Record directions stay unchanged.`
     - any selected target with an unknown strand:
       `Unchanged because a strand is unknown: <record labels>.`
   - Each target row replaces `Effective orientation: …` with one line:
     `Direction: same as reference`, `Direction: opposite to reference`
     (plus ` — reversed on Apply` when the checkbox is on), or
     `Direction: unknown strand — unchanged`. The line always reflects the
     currently selected candidate.
   - Per-target orientation checkboxes are removed. Skip rows keep their
     current position and orientation.
4. **Apply** runs the existing batch Python validation. When the checkbox is
   on, Apply flips the Reverse setting of every aligned target whose selected
   anchor is `opposite` to the reference. It then generates once with the
   validated plan and those orientations, and commits plan, orientations, and
   Result as one History entry. The summary's `N reversed` counts the records
   whose orientation this Apply changed.
5. **Records that never move or reverse:**
   - The reference record never moves or reverses.
   - Skipped, missing, and unusable targets keep position and orientation.
   - Targets with an unknown strand keep their orientation.
   - Every target keeps its vertical position.
   - Alignment is exact and idempotent. Opening a review again after a
     successful match shows every matched target as `same as reference`.
6. **A manual Reverse after alignment keeps the active plan.** No
   `Alignment cleared` notice appears. The next generation aligns the same
   anchors in the new orientation. Other clearing triggers are unchanged:
   manual record movement, source replacement, crop change, and selector
   change.
7. **Reset Align** clears the plan and restores the record positions from
   immediately before that Align. It does not change any record orientation,
   including orientation set by Match reference direction. Notice:
   `Alignment reset: record positions restored; record directions unchanged.`
8. **Undo** restores the complete prior artifact, including orientation and
   any prior plan. Failed, canceled, stale, or superseded work commits nothing
   and changes no record orientation.
9. **Failures show the underlying error.** The review and the global error
   banner both show the underlying failure summary once. The draft, including
   the checkbox state, is retained for retry.
10. **Sessions and exports.**
    - Save and load round-trip the plan anchors and record orientations.
    - Regeneration reproduces the same figure.
    - Exports contain no transient review UI.
11. **The CLI is unchanged.** `--align_orthogroup_feature` preserves
    orientation; no new CLI option is added.

## 4. Architecture

### 4.1 Ownership before and after

| Responsibility | Before | After |
| --- | --- | --- |
| Record orientation value | `region_reverse` and plan `effectiveReverseComplement` | `region_reverse` only (projected once to the request) |
| Applying orientation to the drawing | `resolve_record_inputs()` and late reversal in `materialize_similarity_alignment_display()` | `resolve_record_inputs()` only |
| Displayed-coordinate projection of comparison rows | Browser, from `region_reverse` | Unchanged |
| Strand fact per candidate | `_orientation_result()` combining policy and effect | `strand_relation` (`same`, `opposite`, `unknown`) computed once in Python |
| Rule "reverse only when known opposite" | `_orientation_result()` | One pure function in `similarity-alignment.js` applied to Python facts |
| Anchor-center translation | `_final_record_translations()` | Unchanged; runs after orientation, so it works for any orientation |
| Orientation during one generation | Live state for projection; late patch for serialization | One run-local orientation input read by every stage |
| Plan clearing on edits | `clearForMutation()` and `setManualOrientation()` | `clearForMutation()` only |
| Error shown after a failed generation | `runAnalysis` publishes it; the controller overwrites it | `runAnalysis` publishes and returns it; the controller displays the same error |

### 4.2 Canonical workflow

~~~text
feature popup or Similarity Groups drawer (exact reference)
  -> one controller: gbdraw/web/js/app/similarity-alignment.js
  -> Worker helper -> resolve_similarity_alignment(): anchors, rationale, strand facts
  -> resolved + normal Align: applyPlan(plan, no orientation change)
     ambiguous or explicit review: local draft
        (Select/Skip per target + matchReferenceDirection)
        -> Apply: batch Python validation -> validated plan
        -> matchedOrientations(response, request, enabled)
  -> runAnalysis(canonicalStateOverride {
        similarityAlignmentPlan, linearRecordTranslations,
        linearRecordOrientations })
     orientation override applied once at run entry; LOSAT display projection,
     canonical serialization, and request construction read the same value
  -> typed request -> resolve_record_inputs() applies orientation
  -> project_similarity_alignment_centers() -> _final_record_translations()
  -> SVG
  -> success: plan + orientations + Result commit as one History entry
     failure/cancel/stale: nothing changes; draft retained
~~~

### 4.3 Contracts after the change

Python (`gbdraw/layout/similarity_alignment.py`):

- `AlignmentRecordChoice(record_key, kind, anchor)`
- `AlignmentRecordDecision(record_key, status, rationale, anchor)`
- `AlignmentStrandRelation`: `SAME`, `OPPOSITE`, `UNKNOWN`. An unknown strand
  on either side gives `UNKNOWN`.
- `AlignmentReviewCandidate(candidate, usable, direct_evidence, strand_relation)`
- `SimilarityAlignmentCandidate` without `effective_reverse_complement`
- Removed: `AlignmentOrientationPolicy`, `AlignmentOrientationEffect`,
  `_orientation_result()`, and the `gbdraw.api` export of
  `AlignmentOrientationPolicy`. Add no replacement public export unless a
  documented public example needs one.

Python render (`gbdraw/api/record_planning.py`):

- `project_similarity_alignment_centers(collection, plan) -> tuple[float | None, ...]`
  replaces `materialize_similarity_alignment_display()`. It validates record
  coverage and projects each non-skipped anchor center. It never modifies
  records. Callers use the resolved collection unchanged.

Worker JSON and the persisted plan:

- Choice: `{recordKey, kind, anchor}`.
- Candidate: the current fields, with `orientation` replaced by
  `strandRelation`.
- Decision and plan record: `{recordKey, status, rationale, anchor}`
  (response decisions keep `kind`, `reviewReason`, and `candidates`).
- The plan and helper `schema` integer stays `2`; see section 5.

Web state:

- The draft gains `matchReferenceDirection: boolean`.
- The controller action `setMatchReferenceDirection(enabled)` replaces
  `setOrientation(recordKey, policy)`.
- `setManualOrientation()` and the app-setup helpers
  `setLinearRecordOrientation()` and `linearRecordOrientationValue()` are
  removed. The record Reverse complement checkbox binds
  `v-model="seq.region_reverse"`, as it does on `origin/main`.
- `orientationsFromRequest(request, plan)`, `requestWithOrientations()`, and
  the orientation part of `baseline()` and `installBaseState()` are removed.
  Record drag still materializes translations only.

### 4.4 Design principles

- **Single responsibility:** each part owns one thing.
  - The plan owns anchors and positions; records own orientation.
  - Python owns eligibility, strand facts, and validation.
  - The controller owns the workflow; `runAnalysis` owns generation.
- **Open/closed:** any later displayed-coordinate projection follows record
  orientation automatically, with no alignment-specific code.
- **Interface segregation:** the controller's public API loses the orientation
  operations it no longer owns.
- **KISS:** one checkbox, one orientation path, and no synchronization rule
  between two stores.
- **DRY:**
  - Match reference direction and manual Reverse share the same orientation
    path and History transaction.
  - One error normalizer (`services/error-normalization.js`) serves the banner
    and the review.
- **YAGNI:** none of the following is added:
  - per-record orientation policy;
  - an orientation baseline for Reset;
  - a CLI option or feature flag;
  - a compatibility reader or new schema version;
  - a Python re-projection of comparison rows.
- **Net deletion:** `AGENTS.md` treats every added line as debt. Production code
  should shrink; justify any net growth in the handoff.

### 4.5 Prohibited approaches

- Re-projecting comparison rows in Python after a plan-driven reversal: this
  adds a second projection path.
- Weakening or removing the ID consistency check in `pairwise_match.py`.
- Adding any second orientation store, for example "reversed by alignment"
  flags or a Reset orientation baseline.
- Changing anchor resolution, recommendations, CLI behavior, or the Worker
  lifecycle.
- Adding an alignment-specific history.
- Editing the historical plan directories
  `docs/internal/issue-586-similarity-alignment-ux/`,
  `docs/internal/issue-586-alignment-followup-20260925/`, or
  `docs/internal/ISSUE_561_SIMILARITY_ALIGNMENT_MASTER_PLAN_2026-09-22.md`.

## 5. Persisted-format compatibility

The typed alignment plan (`layout.similarityAlignment`) and
`gbdraw/layout/similarity_alignment.py` are absent from `origin/main` and tag
`0.13.0`; only the legacy `alignOrthogroupFeature` string exists there. Under
the persisted-format rules in `CLAUDE.md`:

- Change the plan field set in place and keep `schema: 2`. Do not advertise a
  branch-only intermediate version or add a reader for the removed fields. A
  Session with the removed fields fails with an explicit validation error.
- Regenerate the branch-owned Gallery artifacts from current code with the
  unfiltered owner command `python tools/refresh_gallery_sessions.py`, then run
  `python tools/gallery_artifact_manifest.py` to verify.
- The legacy `alignOrthogroupFeature` reader (`PD-OI-030`) stays. It stops
  emitting orientation fields, and the legacy record orientation stays in the
  record.
- Update every current document that describes the removed fields.

## 6. Change inventory

| Area | Files | Session |
| --- | --- | --- |
| Product authority | `docs/internal/OPTION_INTEGRITY_PRODUCT_CONTRACT.md` (authority-only branch) | 00 |
| Error reporting | `gbdraw/web/js/app/run-analysis.js` (error outcomes), `gbdraw/web/js/app/similarity-alignment.js` (`applyPlan`, `publishError`), `gbdraw/web/index.html` (duplicated retry sentence) | 01 |
| Python domain and render | `gbdraw/layout/similarity_alignment.py`, `gbdraw/api/record_planning.py`, `gbdraw/api/request_render.py`, `gbdraw/api/session_compat.py`, `gbdraw/api/__init__.py`, `gbdraw/session_request_codec.py`, `gbdraw/web_support/similarity_alignment.py` | 02 |
| Web contract readers | `similarity-alignment.js` (validators, `reviewRows`, `editRow`, `planChoices`), `services/session-request.js`, `services/legacy-similarity-alignment.js`, `index.html` (per-target checkbox removed) | 02 |
| Orientation workflow and UI | `similarity-alignment.js` (draft option, `matchedOrientations`, `applyPlan`, `resetAlignment`, summary, inspector), `run-analysis.js` (single run-entry orientation input), `app-setup.js`, `index.html` | 03 |
| Tests | `tests/test_similarity_alignment*.py`, `tests/test_session_request_codec.py`, `tests/test_session_compat.py`, `tests/test_api_session.py`, `tests/test_documentation_reference_contracts.py`, `tests/web/similarity-alignment-actions.test.mjs`, `tests/web/session-request.test.mjs`, `tests/web/run-analysis-simple-path.test.mjs`, `tests/web/gallery-session-publication.test.mjs`, `tests/web/similarity-alignment-ui.playwright.spec.js` | 01–03 |
| Documentation | `docs/REFERENCE/python-api.md` (executable example S07-PY-01), `docs/REFERENCE/typed-requests.md`, `docs/REFERENCE/session-and-request-compatibility.md`, `docs/SESSION_COMPATIBILITY.md` in 02; `docs/REFERENCE/web-app.md`, `docs/RELEASE_NOTES_0.14.0.md`, Gallery tutorials if affected in 04 | 02, 04 |
| Generated Gallery artifacts | `gbdraw/web/gallery/` session, SVGs, thumbnails, `examples.json`, manifest | 02 (refresh), 04 (final refresh if needed) |

Update a document in the same session that makes it inaccurate. Session 04
performs the final consistency sweep.

## 7. Sessions

| Session | Instruction file | Scope | Start condition |
| --- | --- | --- | --- |
| 00 | [SESSION_00_AUTHORITY.md](SESSION_00_AUTHORITY.md) | Verify and publish the approved authority-only serialization of `PD-OI-027/028/029/031/034`; approval and serialization are done | Now |
| 01 | [SESSION_01_ERROR_REPORTING.md](SESSION_01_ERROR_REPORTING.md) | Underlying error shown once; in-browser confirmation of the root cause | Now (independent of authority) |
| 02 | [SESSION_02_PLAN_CONTRACT.md](SESSION_02_PLAN_CONTRACT.md) | Orientation removed from plan and render; strand facts; contract readers; Gallery refresh; API docs | Authority merged into `origin/dev` and then into the fix branch; Session 01 committed |
| 03 | [SESSION_03_ORIENTATION_WORKFLOW.md](SESSION_03_ORIENTATION_WORKFLOW.md) | One Match option; orientation through record state; manual Reverse keeps plan; Reset | Session 02 committed |
| 04 | [SESSION_04_DOCS_GALLERY_ACCEPTANCE.md](SESSION_04_DOCS_GALLERY_ACCEPTANCE.md) | Public docs, Gallery, integrated acceptance, full gates, record | Sessions 01–03 committed |

Every instruction file is self-contained. At the end of Sessions 00–03, print
the complete next-session file in a copyable code block for the user to paste
into a new agent session. Session 00 hands off to Session 01. Session 01 names
the authority gate before handing off to Session 02.

Between Sessions 02 and 03, the branch intentionally has no Match reference
direction control; the branch is not merged in that state.

## 8. Verification

Focused evidence:

- Python: `tests/test_similarity_alignment.py`,
  `tests/test_similarity_alignment_rendering.py`,
  `tests/test_similarity_alignment_web_adapter.py`,
  `tests/test_session_request_codec.py`, `tests/test_session_compat.py`,
  `tests/test_api_session.py`, `tests/test_documentation_reference_contracts.py`,
  and `tests/test_gallery_session_semantics.py`.
- Web unit: `node --test tests/web/similarity-alignment-actions.test.mjs
  tests/web/session-request.test.mjs tests/web/run-analysis-simple-path.test.mjs
  tests/web/gallery-session-publication.test.mjs`.
- Browser: `tests/web/similarity-alignment-ui.playwright.spec.js` and
  `tests/web/linear-multi-record.playwright.spec.js`.
  - Check both Playwright installations as described in `CLAUDE.md`. If Node's
    `@playwright/test` is missing, run equivalent targeted checks with Python
    Playwright.
  - Rerun a Chromium sandbox failure with the documented escalation.
- Required regression journey, with comparisons enabled, on the Gallery
  example `BGC0000708-BGC0000713`:
  1. Reference `CAG38712.1`, then **Review alignment options…**, then enable
     Match reference direction.
  2. Apply succeeds. The *S. fradiae* ATCC 10745 record is reversed, its
     comparison ribbons are drawn, and anchor centers match the reference
     within 0.5 px.
  3. Reopening the review shows `same as reference` for that record.
  4. Undo restores orientation and positions. Reset Align restores positions
     and keeps the reversal.
  5. A manual Reverse keeps the plan and re-aligns on Generate.
  6. Save and load round-trip the result.

Final gates, run after focused failures are fixed:

- `ruff check gbdraw/`
- `node tests/web/architecture-contracts.test.mjs`
- `node tools/check-web-change-budget.mjs --base origin/dev`
- `python tools/update_cli_reference_help.py --check`
- `python -m pytest tests/ -v -m "not slow"`
- `python -m pytest tests/test_output_comparison.py::TestOutputComparison -v`
- `python -m build`
- `git diff --check`

Prepare the gitignored browser wheel with `python tools/prepare_browser_wheel.py`
when wheel-dependent checks need it. Allow at least 30 minutes for the full
pytest run and monitor it incrementally. Never update
`tests/reference_outputs/` to silence a comparison. The CLI orientation
behavior is unchanged, so no reference SVG change is expected.

## 9. Exit criteria and records

The work is complete only when all of the following hold:

- every behavior in section 3 is demonstrated;
- the regression journey passes with real rendering;
- the authority is on the runtime base;
- every required gate passes;
- production, tests, documentation, and generated artifacts have been reviewed
  as separate diffs.

Architecture evidence is the concise non-increasing form of
[the architecture ratchet](../ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md):

- **Owner excess (OE) decreases:** record orientation goes from two owners to
  one.
- **Path excess (PE) decreases:** applying orientation goes from two paths to
  one, and so does orientation during a run.
- **Compatibility burden (CB) is unchanged:** no reader is added or removed.

If an exception condition applies, provide the full owner, path, and
compatibility sets instead. Record exact commits, commands, results, measured
geometry, visual observations, authority SHAs, and remaining risk in section 10.

The planning commit contains only this directory.

## 10. Acceptance record

Session 04 appends the integrated acceptance record here.
