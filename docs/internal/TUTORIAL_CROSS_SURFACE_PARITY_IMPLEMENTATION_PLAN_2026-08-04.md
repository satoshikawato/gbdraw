# Tutorial cross-surface parity implementation plan

Date: 2026-08-04  
Status: complete  
Scope: the nine projects listed as migration backlog in `docs/scenarios/manifest.json`

## Objective

Complete the GUI, CLI, and Python API variants for every existing Tutorial
project. Each variant must build the project's declared biological figure from
the same public fixtures and preserve record order, comparisons, feature
selection, labels, colors, track geometry, and finished SVG semantics.

## Starting census

The approved manifest has 16 Tutorial chapters. Two projects already have all
three variants. The remaining work is:

- three GUI chapters;
- five CLI chapters;
- nine Python chapters;
- one correction to the existing first-linear CLI recipe so its layout matches
  the GUI figure.

The completed census is 33 Tutorial chapters and 84 chapters overall.

## Scenario allocation

### GUI

- `T-GUI-10`: mitochondrial feature presentation;
- `T-GUI-11`: Majanivirus table-driven comparison;
- `T-GUI-12`: quantitative Hepatoplasma map.

### CLI

- `T-CLI-07`: Lambda-DE3 LOSATN comparison;
- `T-CLI-08`: aminoglycoside BGC Similarity groups;
- `T-CLI-09`: mitochondrial comparison rings;
- `T-CLI-10`: Hepatoplasmataceae Collinear map;
- `T-CLI-11`: interactive session handoff.

### Python API

- `T-PY-03`: first linear genome;
- `T-PY-04`: Lambda-DE3 LOSATN comparison;
- `T-PY-05`: aminoglycoside BGC Similarity groups;
- `T-PY-06`: mitochondrial comparison rings;
- `T-PY-07`: Hepatoplasmataceae Collinear map;
- `T-PY-08`: interactive session handoff;
- `T-PY-09`: mitochondrial feature presentation;
- `T-PY-10`: Majanivirus table-driven comparison;
- `T-PY-11`: quantitative Hepatoplasma map.

## Implementation phases

1. Freeze the current documentation contracts with focused tests and inspect
   each source Tutorial's executable settings.
2. Add the five CLI and nine Python recipes, keeping the literal code in the
   public pages as the executable authority.
3. Add the three GUI visible-interface captures and their owned screenshots.
4. Extend the manifest, fixture ownership, surface indexes, root project table,
   and recipe dispatchers.
5. Run every new recipe in a clean directory, run browser captures from a fresh
   context, and compare each project's durable SVG semantics.
6. Mark a project verified only after all three variants pass their executable
   checks. Keep any unproved project as migration rather than inferring parity.

## Required evidence

- scenario and Tutorial contract tests pass;
- each CLI recipe exits successfully and produces only declared files;
- each Python code block executes unchanged and produces only declared files;
- each GUI capture uses the visible controls, owns its screenshots, and reports
  a clean offline browser run;
- standard SVG safety and biological-identity checks pass;
- published non-browser SVGs retain the declared record order, comparison
  endpoints, track order, labels, and geometry;
- public prose passes the repository AI-writing-pattern audit.

## Definition of done

- no `tutorial_projects.*.variants` value is null;
- `legacy_migration_projects` is empty and all 11 projects are `verified`;
- the root Tutorial table contains links for all 33 variants;
- GUI, CLI, and Python indexes list their manifest chapters once and in order;
- all focused verification is green and the final diff contains no unrelated
  production or generated-package changes.

## Completion evidence

- `docs/scenarios/manifest.json` contains 33 Tutorial chapters and 84 chapters
  overall. All 11 Tutorial projects are `verified`, every project has GUI, CLI,
  and Python variants, and `legacy_migration_projects` is empty.
- The five new CLI recipes (`T-CLI-07` through `T-CLI-11`) and nine new Python
  recipes (`T-PY-03` through `T-PY-11`) passed `--check` from clean temporary
  directories. The corrected `T-CLI-02`, `T-CLI-03`, and paired `T-PY-09`
  outputs also matched their published SVGs.
- Fresh Playwright runs matched both committed screenshots for `T-GUI-10`,
  `T-GUI-11`, and `T-GUI-12`. The comparison permits at most 100 pixels of
  one-channel-unit Chromium antialiasing noise; larger or stronger changes
  fail. Their downloaded SVG validators reported 37, 706, and 576 feature
  elements, respectively; `T-GUI-11` also retained its exact 80-plus-2
  comparison-link set.
- The focused documentation, fixture, and browser-capture suite passed all 117
  tests. Ruff passed for the changed Python capture, recipe, and contract files.
- The public Tutorial prose completed the `avoid-ai-writing` audit with no
  prose flags after excluding Markdown table syntax, code comments, and CLI
  option prefixes. Audit status: `DONE`.
