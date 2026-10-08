---
name: love-me-love-my-docs
description: >-
  Generate or renovate reproducible user documentation for Web, CLI, Python,
  or mobile workflows. Capture real GUI screenshots with executable scripts,
  run documented commands from clean directories, execute literal code
  examples, and bind every visible result to regeneration evidence. Use for
  user manuals, user guides, Tutorials, procedural docs, onboarding docs,
  documentation with screenshots, or requests that mention
  love-me-love-my-docs or /love-me-love-my-docs.
---

# Love Me, Love My Docs

Manuals rot when screenshots are pasted, commands are never rerun, or code
examples are treated as decoration. This skill makes documentation evidence a
build artifact: committed automation performs the real workflow on the surface
being taught, produces the documented result, and regenerates the manual.

## The Prime Directive (family rule)

> **Every visible result is reproducible and every documented action was
> performed on the surface being taught.** GUI images come from committed
> capture automation, CLI commands run in clean directories, and literal
> Python examples execute unchanged. A step the harness cannot perform is a
> finding, not a gap to hide with prose.

## Progress checklist

Copy this into your response and check items off:

```
Docs Progress:
- [ ] Step 1: Frame — audience, language(s), surfaces, output format
- [ ] Step 2: Flow census — user journeys mined from routes/screens; Public-page Decision Gate passed
- [ ] Step 3: Demo data — stable public inputs or safe seeded data; zero private data
- [ ] Step 4: Smoke proof — one workflow runs end-to-end and its result renders in the docs
- [ ] Step 5: Execution harness — GUI capture and/or clean CLI/Python execution per documented flow
- [ ] Step 6: Evidence run — screenshots, commands, code, downloads, and results generated completely
- [ ] Step 7: Manual written — page decisions applied with inputs, actions, expected results, and useful media
- [ ] Step 8: Verify + report — no broken images, rot pre-mortem passed, regeneration command documented
```

## Step 1 — Frame

- **Audience:** end users / admins / both — separate manuals if both
  (mixed audiences make unusable docs).
- **Language(s):** write in the product's language first; if multiple,
  captures may need per-locale runs (the harness parameterizes locale).
- **Surfaces:** select the real interfaces the reader will use. Use Playwright
  for Web, Maestro for explicitly requested mobile work, a subprocess runner
  for CLI, and literal-snippet execution for Python. Do not use evidence from
  one surface to claim another surface works. See
  [references/capture-web.md](references/capture-web.md),
  [references/capture-mobile.md](references/capture-mobile.md), and
  [references/execute-cli-python.md](references/execute-cli-python.md).
- **Output:** MkDocs Material site (recommended — beautiful by default,
  searchable), plain Markdown in `docs/manual/`, or PDF. One choice.

## Step 2 — Flow census

Mine user journeys and evidence scenarios from evidence, not memory: routes,
navigation menus, screen registries, CLI entry points, public APIs, and
existing documentation manifests (`file:line` each). Record the reader goal,
surfaces exercised, screens or commands touched, and any internal/debug flows
that should be skipped. Treat scenarios as verification units, not as a public
page inventory.

**Public-page Decision Gate** — before writing a capture script or public
page, present a compact table with one row per reader question or affected
page: reader question, route/navigation evidence, existing public owner,
evidence scenarios and surfaces, `keep` / `merge` / `delete` / `new`
disposition, and resulting owner. Default to `merge`; use `new` only when no
existing owner can answer the distinct reader question clearly. One or many
evidence scenarios may support one public page, and a scenario may support no
public page. Surface-specific evidence proves that surface works; it does not
require a surface-specific page unless the reader journey is materially
different. Use
[references/manual-structure.md](references/manual-structure.md) for the page
ownership rules. Ask for confirmation **once**. No harness code may exist
before this gate passes. If the user cannot respond (headless/CI run), proceed
and mark the decision table `UNCONFIRMED` in the final report.

## Step 3 — Demo data hygiene

Screenshots outlive databases. Before any capture:

- A dedicated demo account (name like "Somchai Demo", not a real person)
  and seeded content that looks real but is fictional.
- **Zero real user data, emails, tokens, or keys in any frame** — treat
  every screenshot as public forever.
- Consistent state: the same seed produces the same screens, so re-runs
  diff cleanly. Seed script lives next to the capture script.
- For sequence-based tools, make the reader obtain original sequence records
  from an authoritative public database by accession. A public procedural page
  must not use a repository-bundled sequence file or a prebuilt project/session
  as its reader-facing input. Give the database page or API, format choice,
  exact save name, and identity checks. A page may reload only a session the
  reader created earlier in that same page
  from the original inputs.
- Do not assume that an accession's nucleotide version also freezes its
  feature table. When rendered features or protein searches depend on an
  annotation revision, link an official database revision-history snapshot or
  another authoritative annotation release that readers can download, and
  record that exact request in the evidence.
- A frozen repository copy may support deterministic offline automation only
  after a mirror-verification record captures the authoritative database URL
  or exact API request, versioned accession, requested format, UTC retrieval
  date, source byte size and SHA-256, mirror byte size and SHA-256, and the
  direct comparison result. Byte-preserving copies require equal hashes. If a
  documented deterministic normalization is necessary, retain the original
  source hash and executable transformation, then compare the mirror hash with
  a freshly derived output. Accession, length, topology, or parser checks are
  useful additional guards but do not replace the byte comparison.
- Mark a mirror without that evidence `legacy-unverified`. Do not use it to
  regenerate public screenshots or other visible results until it is fetched
  again from the authoritative source and compared or rebuilt. If the exact
  source cannot be recovered, report the blocked evidence path or adopt a new
  authoritative input and regenerate the documentation; never grandfather the
  mirror. Once verified, routine offline runs need only check the pinned mirror
  hash and identity metadata and must still keep the mirror out of the reader's
  acquisition path.
- Distinguish files readers download, create, generate, and compare. Show the
  complete contents or reproducible derivation for support tables and other
  non-sequence inputs.
- **Production gate (mechanical):** seeding writes rows and the harness
  clicks real buttons ("Publish", "Delete") — pointed at production, it
  mutates production. That is a ONE-WAY action. Before seeding one row or
  scripting one click, check the target base URL host: `localhost`,
  `127.0.0.1`, or a `.local`/`.test` domain passes automatically; **any
  other host requires the user to confirm it by name, once**. No
  confirmation available (headless run) = do not seed, do not capture —
  fall back to the degraded mode in Step 5.

## Step 4 — Smoke proof (one workflow end-to-end)

Before writing every flow, prove the thinnest slice works. This is where the
operational unknowns live, and finding them on flow 1 of 1 is cheaper than
finding them after all pages are drafted:

1. Boot the app (document the exact command) and verify the base URL
   responds.
2. Log the demo account in once; save the storage state to
   `docs/capture/auth.json`.
3. Verify one seeded entity is visible on a real page.
4. Capture ONE screenshot through the harness skeleton and render ONE public
   page that references it, in the chosen output format.

For a CLI or Python-only manual, replace items 1–4 with one clean-directory
run of the exact documented command or literal program, validate the named
output, and render that generated result in one public page.

**Exit criterion (mechanical):** one image file exists on disk and one
rendered public page displays it. Until then, no second capture script
may be written. Any numbered item that fails is a named finding (app
won't boot, auth broken, seed invisible) — report it; do not script
around it.

## Step 5 — Build the capture harness

- **Web:** Playwright for Python — patterns in
  [references/capture-web.md](references/capture-web.md): stored auth
  state, fixed viewport/theme, wait-for-stable strategies, element
  highlighting before the shot, per-locale parameterization. **Selector
  policy:** prefer `get_by_role`, then `get_by_label` for a form control with
  an associated accessible label, then a stable `get_by_test_id`. Label
  locators are compliant and do not imply a missing test ID. Every bare-text
  or CSS selector that survives is counted and listed in the final report as
  an app finding (missing accessible semantics or a stable test ID). No running
  app, or Playwright not installable = degraded mode (mirror of the
  mobile one): still commit the harness and seed scripts (they encode the
  steps), emit a manual capture checklist with the exact filenames to
  shoot, mark every unfilled image slot with a visible TODO placeholder,
  and say so honestly in the report.
- **Mobile (optional):** Maestro YAML flows with `takeScreenshot` —
  patterns, alternatives (Fastlane snapshot/screengrab, raw simctl/adb),
  and the install-security rule (auditable channels only — **never
  `curl | bash`**) in
  [references/capture-mobile.md](references/capture-mobile.md). No
  simulator/emulator available = degraded mode: generate the flow files +
  a manual capture checklist, and say so honestly.
- **CLI and Python:** use the execution and artifact rules in
  [references/execute-cli-python.md](references/execute-cli-python.md). Run
  commands from a new temporary directory with only declared inputs. Execute
  the literal published Python block rather than a parallel test-only example.
  Prefer generated diagrams and parsed output checks over terminal screenshots.
- Organize scripts and flows by independently runnable evidence scenario, not
  by public page. Several scenarios may support one page, and a scenario may
  support no page. Separate GUI, CLI, and Python harnesses do not by themselves
  require separate public pages. Name screenshots
  `<scenario-id>/<step-number>-<slug>.png` — the filename IS the step order.
- Commit the harness to `docs/capture/`. It is product code now.

## Step 6 — Capture run

Run every selected harness. Triage failures rather than skipping them: a GUI
failure may expose a missing accessible selector; a CLI/Python failure may
expose an undeclared input, stale option, wrong working-directory assumption,
or output mismatch. Re-run until the evidence set is complete and consistent.

## Step 7 — Write the manual

Structure and style per
[references/manual-structure.md](references/manual-structure.md). Apply the
decision table: extend the resulting owner, merge duplicated material, delete
pages only after preserving unique content or evidence, and create a page only
for a distinct unanswered reader question. Do not clone prose for each
surface merely because its execution evidence is separate. Procedural pages
identify the starting state, every input and filename, the performed actions,
the expected result, verification, and troubleshooting. Write in user language
("select **Publish**"), never developer language ("trigger the POST endpoint").

## Step 8 — Verify and report

- Every image referenced exists; every GUI flow, CLI recipe, and Python example
  runs green end-to-end; the regeneration command is documented in the
  manual's own README (`python docs/capture/run_all.py`, a recipe runner, or
  `maestro test flows/`).
- **Rot pre-mortem** — assume the manual rotted three months from now;
  the known causes are checked mechanically, not pondered:
  - **Selectors:** grep committed scripts for bare text selectors
    (`:has-text`, `text=`) and CSS selectors (`#`, `.`, `[`); count must
    be zero or each one listed as an app finding with `file:line`.
  - **Auth:** the `auth.json` re-mint procedure is documented in the
    manual's README (stored auth state expires).
  - **Seed drift:** the seed script is committed next to the harness and
    `run_all` invokes it (or its README documents the seed command as
    step one).
  - **Environment:** base URL, viewport, theme, and locale are pinned
    constants in one place — never repeated per script.
- Report: pages kept, merged, deleted, or added; screenshots generated; flows
  that failed automation (app findings); and the one-command regeneration
  story.
- Offer the standing suggestion: wire the capture run into CI so UI
  changes that break the manual fail loudly instead of rotting silently.

## When things go wrong

| Situation | Response |
|-----------|----------|
| App won't boot or base URL unreachable | Named finding in Step 4 — report it; do not script around it. Fall back to degraded mode: commit harness skeleton + manual capture checklist with exact filenames, mark image slots with TODO placeholders. |
| Playwright not installable or no running app available | Degraded mode (Step 5): commit harness + seed scripts (they encode steps), emit manual capture checklist, mark unfilled image slots with visible TODO, report honestly. |
| Capture script fails on a specific flow mid-run | Triage (Step 6): missing test-id/accessibility-label in app (finding — report it) or flow changed since census (update census). Never skip; re-run until complete. |
| CLI or Python example works only from the repository root | Declare and copy every input into a new temporary directory, execute the literal published block there, and fix the prose or recipe. Do not preserve a hidden checkout-relative dependency. |
| Generated figure or file differs from the page | Treat the executable recipe as the authority, inspect the semantic difference, then update the prose/artifact together. Do not hand-edit a generated result. |
| Stored auth state (`auth.json`) expires during regeneration | Re-mint procedure documented in manual's README (Step 8 rot pre-mortem). Log demo account in manually once, save storage state, document the exact command. |
| Screenshots inconsistent (different viewport/theme/locale) | Environment constants not pinned (Step 8 rot pre-mortem). Pin BASE, VIEWPORT, locale, color_scheme in one place; never repeat per script. |
| Manual regeneration command undocumented or fails | Step 8 exit criterion failed. Document full chain in manual's README: boot app, seed, re-mint auth if expired, run capture. Test end-to-end before reporting complete. |
