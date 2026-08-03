---
name: love-me-love-my-docs
description: >-
  Auto-generates a beautiful user manual with REAL screenshots captured by
  executable scripts: censuses the app's user flows from routes and
  navigation, builds a Playwright (Python) capture harness for web (the
  primary focus; mobile via Maestro is an optional path), seeds safe demo
  data so no real user data appears in images, writes the manual chapter-per-flow with a
  numbered step per screenshot, and renders it beautifully (MkDocs Material
  site, Markdown, or PDF). The capture scripts are committed, so the manual
  regenerates when the UI changes instead of rotting. Use when the user
  asks to generate a user manual, user guide, documentation with
  screenshots, onboarding docs, or mentions love-me-love-my-docs or
  /love-me-love-my-docs.
license: MIT
argument-hint: "[project-root] [platform: web/mobile] [output: mkdocs/markdown/pdf] [language]"
---

# Love Me, Love My Docs

Manuals rot for one reason: screenshots are pasted, not generated. The day
the UI changes, every image lies. This skill makes screenshots a build
artifact — a committed capture script walks the real app and shoots every
step — so regenerating the whole manual is one command, forever.

## The Prime Directive (family rule)

> **Every screenshot is reproducible and every documented step was
> actually performed.** Images come from the capture script (committed to
> the repo), never from hand-cropped one-offs. Steps come from flows the
> script successfully walked — a step the harness couldn't execute is a
> finding about the app, not a thing to paper over with prose.

## Progress checklist

Copy this into your response and check items off:

```
Docs Progress:
- [ ] Step 1: Frame — audience, language(s), platform, output format
- [ ] Step 2: Flow census — user journeys mined from routes/screens; Chapter Plan Gate passed
- [ ] Step 3: Demo data — target passed the production gate; safe seeded account; zero real user data
- [ ] Step 4: Smoke capture — ONE screenshot end-to-end (boot → auth → seed visible → shot → rendered page)
- [ ] Step 5: Capture harness — Playwright (web) / Maestro (mobile) script per flow, stable selectors only
- [ ] Step 6: Capture run — screenshots generated, consistent and complete
- [ ] Step 7: Manual written — chapter per flow, step per screenshot, rendered beautifully
- [ ] Step 8: Verify + report — no broken images, rot pre-mortem passed, regeneration command documented
```

## Step 1 — Frame

- **Audience:** end users / admins / both — separate manuals if both
  (mixed audiences make unusable docs).
- **Language(s):** write in the product's language first; if multiple,
  captures may need per-locale runs (the harness parameterizes locale).
- **Platform:** **web is the primary focus** → Playwright. Mobile is an
  optional path (Maestro — see
  [references/capture-mobile.md](references/capture-mobile.md)) — take it
  only when the user explicitly asks for a mobile manual.
- **Output:** MkDocs Material site (recommended — beautiful by default,
  searchable), plain Markdown in `docs/manual/`, or PDF. One choice.

## Step 2 — Flow census

Mine the user journeys from evidence, not memory: routes, navigation
menus, screen registries (`file:line` each). Classify:

- **Core flows** — what 80% of users do (login, create X, publish, search)
  → each becomes a chapter.
- **Secondary flows** — settings, profile, exports → grouped chapters.
- **Admin flows** — separate manual or section, per Step 1.
- **Skip** — internal/debug routes, with a note.

**Chapter Plan Gate** — before writing a single capture script, present a
compact brief in chat: one numbered line per chapter (flow, route/nav
evidence `file:line`, screens touched, screenshot count estimate) plus the
skip list. 10–20 lines total; ask for confirmation **once**. Changing the
chapter list here costs one message; changing it after the harness exists
costs a rewrite. No harness code may exist before this gate passes. If the
user cannot respond (headless/CI run), proceed and mark the chapter plan
`UNCONFIRMED` in the final report.

## Step 3 — Demo data hygiene

Screenshots outlive databases. Before any capture:

- A dedicated demo account (name like "Somchai Demo", not a real person)
  and seeded content that looks real but is fictional.
- **Zero real user data, emails, tokens, or keys in any frame** — treat
  every screenshot as public forever.
- Consistent state: the same seed produces the same screens, so re-runs
  diff cleanly. Seed script lives next to the capture script.
- **Production gate (mechanical):** seeding writes rows and the harness
  clicks real buttons ("Publish", "Delete") — pointed at production, it
  mutates production. That is a ONE-WAY action. Before seeding one row or
  scripting one click, check the target base URL host: `localhost`,
  `127.0.0.1`, or a `.local`/`.test` domain passes automatically; **any
  other host requires the user to confirm it by name, once**. No
  confirmation available (headless run) = do not seed, do not capture —
  fall back to the degraded mode in Step 5.

## Step 4 — Smoke capture (one screenshot end-to-end)

Before writing per-flow scripts, prove the thinnest slice works. This is
where every operational unknown lives — dev-server boot, auth, base URL,
seed visibility — and finding them on flow 1 of 1 is cheap; finding them
on flow 12 of 12 is a rewrite:

1. Boot the app (document the exact command) and verify the base URL
   responds.
2. Log the demo account in once; save the storage state to
   `docs/capture/auth.json`.
3. Verify one seeded entity is visible on a real page.
4. Capture ONE screenshot through the harness skeleton and render ONE
   chapter page that references it, in the chosen output format.

**Exit criterion (mechanical):** one image file exists on disk and one
rendered chapter page displays it. Until then, no second capture script
may be written. Any numbered item that fails is a named finding (app
won't boot, auth broken, seed invisible) — report it; do not script
around it.

## Step 5 — Build the capture harness

- **Web:** Playwright for Python — patterns in
  [references/capture-web.md](references/capture-web.md): stored auth
  state, fixed viewport/theme, wait-for-stable strategies, element
  highlighting before the shot, per-locale parameterization. **Selector
  policy:** committed scripts use `get_by_test_id` / `get_by_role` —
  every bare text or CSS selector that survives is counted and listed in
  the final report as an app finding (missing `data-testid`s). No running
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
- One script/flow per chapter; screenshots named
  `<chapter>/<step-number>-<slug>.png` — the filename IS the step order.
- Commit the harness to `docs/capture/`. It is product code now.

## Step 6 — Capture run

Run the harness. Every failure is triaged, not skipped: a step that can't
be automated is either a missing test-id/accessibility-label in the app (a
finding — report it) or a flow that changed since the census (update the
census). Re-run until the set is complete and consistent (same viewport,
same theme, same locale per set).

## Step 7 — Write the manual

Structure and style per
[references/manual-structure.md](references/manual-structure.md): chapter
per flow — goal, prerequisites, numbered steps (one screenshot each, with
a caption saying what to notice), expected result, troubleshooting.
Written in user language ("click **Publish**"), never developer language
("trigger the POST endpoint"). Render the chosen format; MkDocs Material
config included when that's the choice.

## Step 8 — Verify and report

- Every image referenced exists; every capture script runs green
  end-to-end; the regeneration command is documented in the manual's own
  README (`python docs/capture/run_all.py` or `maestro test flows/`).
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
- Report: chapters written, screenshots generated, flows that failed
  automation (app findings), and the one-command regeneration story.
- Offer the standing suggestion: wire the capture run into CI so UI
  changes that break the manual fail loudly instead of rotting silently.

## When things go wrong

| Situation | Response |
|-----------|----------|
| App won't boot or base URL unreachable | Named finding in Step 4 — report it; do not script around it. Fall back to degraded mode: commit harness skeleton + manual capture checklist with exact filenames, mark image slots with TODO placeholders. |
| Playwright not installable or no running app available | Degraded mode (Step 5): commit harness + seed scripts (they encode steps), emit manual capture checklist, mark unfilled image slots with visible TODO, report honestly. |
| Capture script fails on a specific flow mid-run | Triage (Step 6): missing test-id/accessibility-label in app (finding — report it) or flow changed since census (update census). Never skip; re-run until complete. |
| Stored auth state (`auth.json`) expires during regeneration | Re-mint procedure documented in manual's README (Step 8 rot pre-mortem). Log demo account in manually once, save storage state, document the exact command. |
| Screenshots inconsistent (different viewport/theme/locale) | Environment constants not pinned (Step 8 rot pre-mortem). Pin BASE, VIEWPORT, locale, color_scheme in one place; never repeat per script. |
| Manual regeneration command undocumented or fails | Step 8 exit criterion failed. Document full chain in manual's README: boot app, seed, re-mint auth if expired, run capture. Test end-to-end before reporting complete. |
