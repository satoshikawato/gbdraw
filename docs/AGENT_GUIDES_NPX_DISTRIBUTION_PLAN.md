# gbdraw Agent Guides Distribution Plan

This document defines a staged plan for publishing gbdraw AI-agent guidance.
The first milestone is to lock the repository-local prerequisites that make the
guidance safe to generate and test. The second milestone is a reliable,
versioned `agent-guides/` tree in this repository. The `npx` installer is a
later milestone and must remain a thin distribution wrapper around the canonical
guidance.

## Current Constraints

- The repository currently has a `.codex` file, not a `.codex/` directory. At
  the time this plan was revised, that file was an empty regular file. A
  default Codex install target of `.codex/skills/gbdraw/` will fail in this
  repository until that file is removed, migrated, or intentionally preserved
  with a different local smoke-test target.
- `docs/CLI_Reference.md` claims to mirror current command help, but it has
  already drifted from the installed CLI version. Agent guidance must not add
  another hand-maintained copy of the same CLI surface without validation.
- gbdraw's CLI is broad enough that a small set of examples is not sufficient:
  guidance must distinguish circular, linear, GUI, Python API, GFF3+FASTA,
  BLAST outfmt 6/7, circular conservation rings, depth tracks, label filtering,
  custom track slots, and LOSATP protein comparison workflows.
- `print --agent <name>` cannot rely on local relative references unless it
  prints a bundle. A printed Codex adapter must be useful as a standalone file.
- OpenClaw install conventions are not yet specified well enough to test. The
  v1 package must not advertise an OpenClaw adapter until its target path,
  filename, conflict handling, and smoke tests are defined.

## Goals

- Keep canonical agent guidance in the gbdraw repository and version it with
  CLI/API changes.
- Resolve the repository's `.codex` file policy before running local Codex
  install smoke tests.
- Make generated CLI help the source of truth before adding agent guidance that
  contains command options or defaults.
- Make the Codex `SKILL.md` useful even when the reference files are missing or
  when the skill is printed to stdout.
- Provide task-oriented references that help agents choose the right gbdraw
  mode and options instead of copying stale examples.
- Validate generated or copied CLI help so docs, references, and package
  templates do not silently drift.
- Keep v1 scoped to Codex, Claude, and generic guidance until OpenClaw behavior
  is fully specified.
- Add an `npx` installer only after the canonical guidance is reviewable and
  useful without npm.

## Non-Goals

- Do not distribute the Python `gbdraw` package through npm.
- Do not run `gbdraw`, BLAST, LOSATP, `pip`, `conda`, or other external tools
  from the installer.
- Do not silently edit existing `CLAUDE.md`, `AGENTS.md`, `.codex`, or other
  agent configuration files.
- Do not require a JavaScript build step for the gbdraw web UI.
- Do not ship an OpenClaw adapter in v1.
- Do not make the npm package or its `doctor` command execute Python, `gbdraw`,
  BLAST, LOSATP, `pip`, `conda`, or shell validation commands.

## Implementation Strategy

### Phase 0: Lock repository prerequisites

Complete these items before writing `SKILL.md` or npm templates:

1. Decide and document the repository-local `.codex` policy.
   - Preferred if the file is unused: remove the empty `.codex` file and allow
     `.codex/skills/gbdraw/` in this repository.
   - If the file is intentional: keep it, but local Codex smoke tests must use a
     different output directory and must still verify that `.codex` as a file is
     a hard conflict for normal user installs.
2. Add a CLI reference helper and check mode:
   - `python tools/update_cli_reference.py`
   - `python tools/update_cli_reference.py --check`
3. Refresh `docs/CLI_Reference.md` from current command help.
4. Ensure the same generated help can populate or check
   `agent-guides/references/cli.md` once `agent-guides/` exists.
5. Do not add command examples with option defaults to agent guidance until the
   CLI help check exists.

### Phase 1: Canonical repository guidance

Create and review a canonical `agent-guides/` tree first:

```text
agent-guides/
├── core/
│   └── gbdraw-core.md
├── references/
│   ├── api.md
│   ├── cli.md
│   ├── tables.md
│   ├── troubleshooting.md
│   └── workflows.md
├── codex/
│   └── gbdraw/
│       └── SKILL.md
├── claude/
│   └── CLAUDE-gbdraw.md
└── generic/
    └── gbdraw-agent-prompt.md
```

Phase 1 should not add `packages/agent-guides-npm/`. The first review target is
the accuracy and usefulness of the guidance itself.

OpenClaw is deferred from v1. Add it only after a separate specification defines
the adapter filename, install target, merge/conflict policy, `print` behavior,
and smoke tests.

### Phase 2: npx distribution wrapper

After Phase 1 is stable, add the npm package:

```text
packages/
└── agent-guides-npm/
    ├── package.json
    ├── bin/
    │   └── gbdraw-agent-guides.mjs
    ├── scripts/
    │   └── sync-templates.mjs
    ├── templates/
    │   ├── core/
    │   ├── references/
    │   ├── codex/
    │   ├── claude/
    │   └── generic/
    └── test/
        └── smoke.mjs
```

The npm package must be generated from or checked against `agent-guides/`.
Manual template copying is acceptable only if a sync check fails the release
when `templates/` differs from the canonical source.

## `.codex` Conflict Policy

Before enabling `init --agent codex`, decide how this repository should handle
the existing `.codex` file:

- If the file is unused, remove it and allow `.codex/skills/gbdraw/`.
- If the file is intentional, keep it and choose a different default target for
  this repository's local smoke tests.
- In all cases, the installer must treat an existing `.codex` file in a user
  repository as a hard conflict. It must not delete, replace, or rename that
  file, even with `--force`.

Local repository tests must cover both cases:

- successful Codex install into a temporary directory where `.codex` is absent
  or already a directory
- refused Codex install into a temporary directory where `.codex` is a file

The `doctor` command should report this case explicitly:

```text
Cannot install Codex skill: .codex exists and is a file.
Expected a directory for .codex/skills/gbdraw/.
```

## Codex SKILL.md Requirements

`agent-guides/codex/gbdraw/SKILL.md` should be concise, but not merely a router.
It must contain enough information to run safely without reference files:

- Mode selection:
  - single circular genome diagram: `gbdraw circular`
  - multiple linear genome comparison: `gbdraw linear`
  - browser workflow: `gbdraw gui`
  - automation/pipeline workflow: `gbdraw.api`
- Input selection:
  - GenBank/DDBJ flatfiles use `--gbk`
  - GFF3 requires paired FASTA: `--gff ... --fasta ...`
  - do not mix `--gbk` with `--gff --fasta`
- Comparison selection:
  - linear nucleotide comparison uses existing BLAST outfmt 6/7 via `-b` or
    `--blast`
  - circular conservation rings use existing BLAST outfmt 6/7 via
    `--conservation_blast`
  - LOSATP/protein blastp modes may execute external tools and require explicit
    user confirmation before running
- Output selection:
  - default to SVG
  - PNG/PDF/EPS/PS require CairoSVG and its platform dependencies
- Safety checks before execution:
  - confirm input files and output prefix
  - check whether the output files already exist
  - confirm external tools before BLAST/LOSATP/protein comparison workflows
  - prefer existing comparison tables over generating new ones automatically
- Verification:
  - run the smallest relevant command first when possible
  - inspect produced SVG or export files
  - use targeted pytest checks when changing code or documented behavior

Reference files can add depth, examples, table formats, and troubleshooting, but
they must not be required for these baseline decisions.

Keep the boundary explicit:

- `SKILL.md` should be a compact decision and safety document, not an option
  catalog.
- It may include a few skeletal command patterns with placeholders, but it must
  not copy full help output or hand-maintain defaults.
- Full command help belongs in generated `references/cli.md`.
- Workflow examples, table schemas, and troubleshooting belong in the
  corresponding reference files.

## Reference Design

`references/workflows.md` should be a task decision tree, not only a command
gallery. It should cover at least:

- single circular plot from GenBank/DDBJ
- single circular plot from GFF3+FASTA
- multi-record circular canvas
- circular conservation rings from precomputed BLAST
- linear plot without comparison data
- linear BLAST comparison from precomputed outfmt 6/7
- linear LOSATP/protein blastp comparison
- region cropping and record selection
- depth tracks and multiple logical depth tracks
- non-SVG export with CairoSVG
- GUI workflow and local server behavior
- Python API workflow

`references/cli.md` should be generated or verified from command help, not
hand-maintained as a second source of truth. Narrative examples belong in
`references/workflows.md`.

`references/tables.md` should define TSV schemas for color tables, default
colors, label whitelist/blacklist, qualifier priority, label overrides, feature
visibility, and any comparison/depth inputs that agents are expected to prepare
or validate.

`references/troubleshooting.md` should include common validation failures:
missing paired FASTA for GFF3, BLAST outfmt mismatch, CairoSVG missing for
non-SVG export, output overwrite risk, unsupported record selectors, and
external LOSATP availability.

## CLI and Documentation Synchronization

Add a documentation helper in Phase 0, before writing agent guidance that
contains command options or defaults:

```bash
python tools/update_cli_reference.py
python tools/update_cli_reference.py --check
```

The helper should capture:

- `python -m gbdraw.cli --help`
- `python -m gbdraw.cli circular --help`
- `python -m gbdraw.cli linear --help`
- `python -m gbdraw.cli --version`

The check mode should fail when generated help differs from
`docs/CLI_Reference.md` or `agent-guides/references/cli.md`. This turns CLI help
into the source that agent guidance and published docs can safely mirror.

This helper is a repository development and CI tool. It may execute
`python -m gbdraw.cli` in the development environment. It must not be called by
the npm package, by `npx @gbdraw/agent-guides ...`, or by the npm `doctor`
command.

## npx User Experience

Phase 2 command shape:

```bash
npx @gbdraw/agent-guides@latest init --agent codex
npx @gbdraw/agent-guides@latest init --agent claude
npx @gbdraw/agent-guides@latest init --agent generic
npx @gbdraw/agent-guides@latest init --agent all
npx @gbdraw/agent-guides@latest print --agent codex
npx @gbdraw/agent-guides@latest print --agent codex --bundle
npx @gbdraw/agent-guides@latest doctor
```

Expected behavior:

- `init --agent codex` writes a Codex skill under `.codex/skills/gbdraw/` only
  when `.codex` is absent or already a directory.
- `init --agent claude` writes `CLAUDE-gbdraw.md` and prints merge guidance when
  `CLAUDE.md` already exists.
- `init --agent generic` writes `gbdraw-agent-prompt.md`.
- `init --agent all` writes the supported v1 adapters and shared references,
  subject to the same conflict checks. For v1, supported adapters are Codex,
  Claude, and generic.
- `print --agent <name>` prints a standalone adapter.
- `print --agent <name> --bundle` prints the adapter plus references in a clear
  multi-file format.
- `doctor` validates template consistency and reports existing local agent files.
  It must not run Python, `gbdraw`, BLAST, LOSATP, package managers, or shell
  validators.

The installer should require explicit flags for normal overwrites:

```bash
npx @gbdraw/agent-guides@latest init --agent codex --force
```

`--force` must not override hard conflicts such as `.codex` being a file.

## npx Package Requirements

Use only Node.js built-ins for v1:

- `node:fs/promises`
- `node:path`
- `node:url`
- `node:process`
- `node:crypto` for template checksums if needed

Command parser requirements:

- Accept `init`, `print`, `doctor`, `--help`, and `--version`.
- Accept `--agent codex|claude|generic|all`.
- Accept `--out <dir>` with default `process.cwd()`.
- Accept `--force` for ordinary overwrite conflicts.
- Exit nonzero on unknown commands, unknown agents, missing templates, template
  drift, hard conflicts, and refused overwrites.
- Print every file written and every conflict detected.

Minimal package metadata:

```json
{
  "name": "@gbdraw/agent-guides",
  "version": "0.1.0",
  "description": "Agent instructions for using gbdraw with Codex, Claude, and other AI coding agents.",
  "license": "MIT",
  "type": "module",
  "bin": {
    "gbdraw-agent-guides": "./bin/gbdraw-agent-guides.mjs"
  },
  "files": [
    "bin",
    "templates",
    "scripts",
    "test"
  ],
  "engines": {
    "node": ">=18"
  },
  "keywords": [
    "gbdraw",
    "bioinformatics",
    "genome-diagram",
    "ai-agent",
    "codex",
    "claude"
  ]
}
```

## Safety Model

Separate installer safety from agent-operation safety.

Installer safety:

- Never execute `gbdraw`, BLAST, LOSATP, `pip`, `conda`, or shell commands.
- Never download files after npm has installed the package.
- Never overwrite files unless `--force` is present.
- Never treat hard conflicts as forceable.
- Print a summary of every write and conflict.
- Keep `doctor` limited to template checks, expected file/path checks, and local
  conflict reporting.

Agent-operation safety:

- Before proposing a command, identify input files, output prefix, and expected
  output files.
- Check whether outputs already exist before running examples.
- Do not run BLAST, LOSATP, protein comparison modes, package installs, or GUI
  browser launch steps without explicit user confirmation.
- Prefer SVG for first-pass validation and mention CairoSVG for non-SVG export.
- Keep code changes separate from generated diagram output unless the user asks
  for both.

## Testing

Phase 0 checks:

- The `.codex` policy decision is recorded before Codex install smoke tests are
  enabled.
- `python tools/update_cli_reference.py --check` passes after
  `docs/CLI_Reference.md` is refreshed.

Phase 1 checks:

- Markdown link/path checks for files under `agent-guides/`.
- Generated `agent-guides/references/cli.md` is checked by the Phase 0 CLI
  helper once `agent-guides/` exists.
- Review `SKILL.md` as a standalone printed file.

Phase 2 checks:

- Temporary-directory smoke tests for `init`, `print`, `print --bundle`,
  overwrite refusal, `--force`, `.codex` file conflict, and `doctor`.
- Template sync check between `agent-guides/` and
  `packages/agent-guides-npm/templates/`.
- `npm pack --dry-run` check confirming only expected files are included.
- A negative check that npm `doctor` does not invoke Python, `gbdraw`, BLAST,
  LOSATP, package managers, or shell commands.

## Milestones

1. Resolve or document the repository's existing `.codex` file policy.
2. Add CLI help generation/checking and refresh `docs/CLI_Reference.md`.
3. Add `agent-guides/` with self-contained Codex skill and task-oriented
   references.
4. Add README documentation linking to canonical `agent-guides/`.
5. Review guidance quality using real gbdraw circular, linear, GUI, and API
   tasks.
6. Add the dependency-free npm wrapper for Codex, Claude, and generic adapters,
   plus template sync tests.
7. Publish `@gbdraw/agent-guides@0.1.0` only after Phase 2 smoke tests pass.
8. Define OpenClaw adapter behavior in a separate future plan before adding it
   to the package.

## README Update

After Phase 1, add repository-local documentation:

````markdown
## Use gbdraw with AI agents

gbdraw provides versioned agent instructions under
[`agent-guides/`](./agent-guides/). Start with the adapter for your agent, such
as the Codex skill at `agent-guides/codex/gbdraw/SKILL.md`.

The guidance files do not install or run gbdraw. They describe safe command,
GUI, and API workflows for agents helping with gbdraw tasks.
````

After Phase 2, add the optional npx install command:

````markdown
You can also copy the guidance into another repository with:

```bash
npx @gbdraw/agent-guides@latest init --agent codex
npx @gbdraw/agent-guides@latest print --agent claude
```

The installer only writes instruction files. It does not install or run gbdraw.
````
