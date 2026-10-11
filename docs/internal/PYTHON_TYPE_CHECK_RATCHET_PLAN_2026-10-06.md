# Python type-check ratchet plan

Status: proposed, 2026-10-06. The Owner approved the approach on 2026-10-06:
a pinned type checker runs in CI over the `gbdraw` package as a ratchet, no
file gains a new type error, and each file's error count is recorded in a
baseline that only goes down.
Baseline: `dev` at `c9c3b2e5`. Every measurement below was taken at
`01694c3a`; no file under `gbdraw/` changed between the two.
Scope: `gbdraw/**/*.py` except `gbdraw/web/`, `pyproject.toml`, the `core`
and `core-pr` jobs of `.github/workflows/test.yml`,
`tests/test_type_check_ratchet.py`, root `CLAUDE.md`, `.gitignore`.

## 1. Objective

The Python sources carry type annotations almost everywhere: frozen
dataclasses, `NamedTuple`, `Literal` option domains, and typed request and
session contracts. Nothing checks them. CI runs Ruff only, and the repository
has no mypy or pyright configuration. This plan adds a checker as a ratchet,
in the same spirit as the Web rule R14 "Typed boundaries":

1. A checker at a pinned version runs in CI.
2. No file gains a new type error.
3. Each file's error count is recorded in a baseline that only goes down.
4. The contract modules (slices S1a and S1b) have no type debt: no errors and
   no `# type: ignore`.

Bringing the whole package to zero is not a completion condition of this
plan. Section 10 gives the outlook.

## 2. Evidence

### 2.1 Tooling and CI

| Question | Finding | Evidence |
| --- | --- | --- |
| Is a checker configured? | No. No `[tool.mypy]`, `mypy.ini`, or `pyrightconfig.json`. The `Lint` job installs `ruff==0.15.12` and runs `ruff check gbdraw/`. | `pyproject.toml`; `.github/workflows/test.yml:1019-1023` |
| Existing suppressions | 412 `# type: ignore` comments in `gbdraw/`. 331 of them are `# type: ignore[reportMissingImports]` on third-party imports, a pyright rule name, added with the refactoring of 2026-01-09 (#100). The CLAUDE.md convention "`type: ignore` comments for BioPython (missing type stubs)" describes them. | `grep -rn "type: *ignore" gbdraw` |
| Which jobs run unmarked pytest tests? | Only `Core PR` (pull requests, Python 3.11) and `Core` (dev pushes, 3.10/3.11/3.12). Both install `-e ".[export]" pytest pytest-timeout pytest-xdist Pillow` and select `not slow and not (recipe or gallery or browser)`. Every other pytest job selects a marker (`recipe`, `gallery`, `browser`, `slow`, acceptance surfaces). | `.github/workflows/test.yml:173-180,207-214` |
| Can the `Lint` job host a pytest guard? | Not without installing the package: `tests/conftest.py` imports `gbdraw.session_io`, which imports pandas and BioPython. | `tests/conftest.py:16` |
| Does `Core PR` run for every change that can alter a checked file? | Yes. A `.py` path under `gbdraw/` outside `gbdraw/web/` classifies as `python-core`, `renderer`, `session-persistence`, `losat-integration`, `packaging` (`gbdraw/_build_support.py`), or `full` (the default for an unmatched path), and all of these route `core-pr`, as do `tests-only` and `ci-only`. The classes without `core-pr` (`documentation`, `metadata`, `policy-documentation`, `web-runtime`, `gallery`) cannot change a checked file. A dev push of any non-light class runs every dev job. | `tools/ci-impact-policy.mjs:42-56,158-210` |
| Is a workflow change governed? | Yes. `.github/workflows/test.yml` is CI authority, so changing it is a GOVERNANCE pull request. The Web Gate's co-change block concerns production runtime paths, which it defines as `gbdraw/web/index.html`, `gbdraw/web/js/`, and `gbdraw/web/vendor/`. Python sources are not among them. | `docs/internal/WEB_CHANGE_POLICY.md:256-270`; `tools/check-web-change-budget.mjs:333-335,1830-1836` |
| Are annotations ever runtime data? | Yes, in three ways. (1) `get_type_hints` reads dataclass field annotations: `gbdraw/config/modify.py` for `GbdrawConfig` and its nested dataclasses (config overrides), and `_validate_dataclass_contract` in `gbdraw/session_request_codec.py` for request options, layouts, outputs, track slots, and other decoded or encoded dataclasses. A field annotation therefore changes behavior, and a name used in it must exist at run time. (2) `tests/test_public_contract.py` hashes the signatures, annotation strings included, of every name in `gbdraw.__all__`. (3) Only 143 of the 226 checked files have `from __future__ import annotations`. In the other 83 (51 of them annotated, 36 with type debt), annotations are evaluated when the module is imported. | `gbdraw/config/modify.py:76-95`; `gbdraw/session_request_codec.py:4579`; `tests/test_public_contract.py` |
| Are there platform branches? | `gbdraw/losat_setup.py:155-171` uses `msvcrt.locking` under an `os.name` check, which a checker cannot narrow. `gbdraw/comparisons/losat_runtime.py:216-221` uses `sys.platform`, which it can. | source |

### 2.2 Checkers compared

Both checkers ran over `gbdraw/` without `gbdraw/web/`, targeting Python
3.10. "CI-like environment" is Python 3.11 with `pip install -e ".[dev]"` on
2026-10-06: pandas 3.0.6, numpy 2.4.6, biopython 1.88. The local base
environment has pandas 2.3.0, numpy 2.2.5, and biopython 1.85.

| Checker and setting | Errors | Files with errors | Wall time, memory |
| --- | ---: | ---: | --- |
| pyright 1.1.414 basic, base environment | 1,087 | 112 | 14 s, 1.5 GB |
| pyright 1.1.414 basic, CI-like environment | 1,014 | 110 | 14 s, 1.5 GB |
| pyright 1.1.414 basic, CI-like + `pandas-stubs` 3.0.5.260914 | 922 | 106 | 14 s |
| pyright 1.1.414 basic, no site-packages (empty venv) | 739 | 78 | 13 s, 1.0 GB |
| mypy 2.4.0 defaults, CI-like environment | 738 | 103 | 4.6 s cold, 0.1 s warm, 360 MB |
| mypy 2.4.0 defaults, CI-like + `pandas-stubs` | 774 | — | — |
| mypy 2.4.0, `--no-site-packages` | 648 | 94 | 3.1 s cold, 250 MB |
| mypy 2.4.0, `--no-site-packages --check-untyped-defs` (the configuration in 4.3) | 670 | 95 | 3-7 s cold, 240 MB |

Findings:

- **Counts that read installed packages depend on their versions.** The same
  pyright run gives 1,087 or 1,014 errors depending on which pandas, numpy,
  and BioPython are installed. `pyproject.toml` pins none of them, and CI
  installs the newest. A per-file baseline built on such counts could fail an
  unrelated pull request on the day a library releases.
- **Without site-packages the count is reproducible.** The mypy configuration
  in 4.3 produced byte-identical output on Python 3.10, 3.11 (with every
  dependency installed), 3.12, and 3.13 (with none installed).
- **What third-party typing would add.** BioPython 1.85 and later ship
  `py.typed`. With site-packages, its types add 94 mypy errors, mostly `len()`
  of `SeqRecord.seq` (`Seq | MutableSeq | None`) and uses of `record.id`
  (`str | None`); 4 errors exist only without them. numpy, which also ships
  `py.typed`, adds none in mypy. pandas ships no `py.typed`, so it is `Any`
  either way. `pandas-stubs` changes the counts in opposite directions in the
  two checkers: mypy +36, pyright -92.
- **Counting granularity.** pyright reports a `**dict[str, object]` argument
  once per parameter it fills. In `gbdraw/interface.py`, pyright reports 122
  errors and mypy 13 for the same call sites.
- **Integration.** mypy is a Python package with compiled wheels for 3.10 and
  later, so `pip` installs it. pyright needs Node. Its PyPI wrapper downloads
  the npm package the first time it runs, and the npm package would move
  Python checking into the jobs that run `npm ci`.

### 2.3 Error composition (configuration in 4.3)

670 errors by mypy code: `arg-type` 269, `union-attr` 100, `assignment` 91,
`attr-defined` 86, `no-redef` 30, `call-overload` 28, `var-annotated` 24,
`index` 12, `return-value` 10, others 20.

Grouped from the messages:

| Cause | Errors |
| --- | ---: |
| A value typed `object`, usually from a `dict[str, object]` bag passed with `**` or read back from it | 149 |
| `None` not narrowed | 105 |
| A non-`None` union not narrowed | 57 |
| A variable annotated twice in two branches (`no-redef`) | 30 |
| A variable that needs an annotation (`var-annotated`) | 24 |
| Other mismatches (`Sequence` where `list` is declared, `str` where a `Literal` is declared, wrong attribute names, and so on) | 305 |

Of the 412 `# type: ignore` comments, 352 suppress nothing in this
configuration (`--warn-unused-ignores`): the 331 import comments and 21
others. The other 60 suppress a live error.

Cases that the slices must examine, not yet confirmed as defects:

- `gbdraw/interface.py` passes a title position typed
  `Literal['none', 'center', 'top', 'bottom'] | None` to
  `CircularOutputOptions`, which declares `Literal['none', 'top', 'bottom']`.
- `gbdraw/interface.py:1142-1197` reads `.drawing`, `.annotation_warnings`,
  and `.feature_identity_notices` from
  `PreparedDiagramRequest | PreparedCircularBatchRequest`; the batch type has
  none of them.
- `gbdraw/diagrams/linear/assemble.py:2539` sets
  `canvas_config.height_below_final_record` on a `LinearCanvasConfigurator`,
  which does not declare it, and `gbdraw/diagrams/linear/builders.py:446-448`
  reads it back: one layer writes an undeclared attribute into another
  layer's object.
- `gbdraw/api/request_render.py` passes `**dict[str, object]` bags to
  `read_feature_placement_table`, `read_feature_override_table`, and the
  diagram builders (about 40 errors in one file).

### 2.4 Type debt and slices

Type debt is defined in D3: a file's mypy errors plus its `# type: ignore`
comments. At `c9c3b2e5` the debt is 1,082 in 138 of 226 files. T1 removes the
352 dead comments, which leaves 730 in 107 files.

| PR | Slice | Modules | Lines | Errors | `# type: ignore` (dead) | Debt at T0 (files) | Debt after T1 (files) |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| S1a | Persisted contracts | 22 | 15,954 | 111 | 29 (19) | 140 (12) | 121 (10) |
| S1b | Public API | 11 | 12,535 | 122 | 67 (67) | 189 (11) | 122 (6) |
| S2 | Entry points | 15 | 7,864 | 59 | 9 (6) | 68 (7) | 62 (7) |
| S3 | Assembly | 39 | 17,498 | 106 | 109 (105) | 215 (27) | 110 (19) |
| S4 | Rendering | 61 | 13,873 | 79 | 100 (100) | 179 (35) | 79 (24) |
| S5 | Features and I/O | 50 | 18,897 | 86 | 35 (22) | 121 (26) | 99 (23) |
| S6 | Analysis and comparisons | 16 | 13,066 | 67 | 60 (32) | 127 (12) | 95 (10) |
| S7 | Web support | 12 | 5,713 | 40 | 3 (1) | 43 (8) | 42 (8) |
| | Total | 226 | 105,400 | 670 | 412 (352) | 1,082 (138) | 730 (107) |

Module lists:

- S1a: `session_request_codec.py`, `session_io.py`, `session.py`,
  `linear_comparison.py`, `api/session_compat.py`, `config/` (10),
  `mode_profiles.py`, `exceptions.py`, `tracks/` (5). Largest after T1:
  `session_request_codec.py` 51, `session_io.py` 28,
  `linear_comparison.py` 17, `api/session_compat.py` 16.
- S1b: `api/` except `session_compat.py`. Largest: `api/request_render.py`
  86, `api/diagram.py` 18, `api/record_planning.py` 9, `api/options.py` 6.
- S2: `interface.py`, `cli.py`, `cli_utils/`, `circular.py`, `linear.py`,
  `__init__.py`, `version.py`, `crop_genbank.py`, `losat_setup.py`,
  `_build_support.py`, `_web_assets.py`, `data/__init__.py`. Largest:
  `interface.py` 16, `linear.py` 14, `circular.py` 14, `losat_setup.py` 9.
- S3: `diagrams/`, `canvas/`, `configurators/`, `layout/`. Largest:
  `diagrams/circular/radial_layout.py` 25, `diagrams/circular/assemble.py`
  14, `configurators/legend.py` 14, `diagrams/linear/assemble.py` 13.
- S4: `render/`, `svg/`, `definition_line_styles.py`. Largest:
  `render/groups/circular/legend.py` 20, `render/interactive_svg.py` 15,
  `render/interactive_context.py` 7.
- S5: `features/`, `io/`, `labels/`, `legend/`, `core/`, `annotations/`.
  Largest: `features/placement.py` 18, `io/cli_tables.py` 14,
  `features/overrides.py` 10.
- S6: `analysis/`, `comparisons/`. Largest:
  `analysis/protein_colinearity.py` 52, `analysis/ortholog_paths.py` 14.
- S7: `web_support/`. Largest: `web_support/feature_catalog.py` 15,
  `web_support/error_adapter.py` 8.

Python already declares the persisted key sets that the Web typedefs mirror:
`_TOP_LEVEL_FIELDS_V5` (`gbdraw/session_request_codec.py:197`) and
`CURRENT_SESSION_TOP_LEVEL_FIELDS` (`gbdraw/session_io.py:62`). The parity
tests in `tests/test_session_request_codec.py` compare them with the Web
typedefs.

## 3. Decisions

All are Owner-delegated: each recommended option is adopted under the
standing instruction of 2026-09-29. Section 9 lists the ones the Owner may
want to override before T0.

- **D1 mypy 2.4.0, exact pin.** It is reproducible without site-packages
  (2.2), installs with `pip`, runs in 3-7 s and about 250 MB, and counts
  close to one error per site. pyright was the alternative. Its strengths are
  that it matches Pylance in VS Code and infers through unannotated code. Its
  counts depend on installed packages unless it is pointed at an empty
  virtual environment. It needs Node, and it takes 13-14 s and 1-1.5 GB. The
  pin lives in a new `typecheck` extra (`typecheck = ["mypy==2.4.0"]`), and
  `dev` includes it as `gbdraw[typecheck]`. A checker upgrade is its own pull
  request, because a new version can add diagnostics; one that raises an
  entry needs the Owner's approval.
- **D2 Configuration.** `[tool.mypy]` sets `python_version = "3.10"` (the
  oldest supported version), `platform = "linux"`, `files = ["gbdraw"]`,
  `exclude = ["^gbdraw/web/"]`, `no_site_packages = true`, and
  `check_untyped_defs = true`. One override sets `ignore_missing_imports` for
  the seven third-party top-level modules imported by `gbdraw/`: `BCBio`,
  `Bio`, `fontTools`, `numpy`, `pandas`, `svgwrite`, and `tomli`. Listing
  them, instead of setting `ignore_missing_imports` globally, keeps a
  mistyped first-party import an error. All other options keep mypy's
  defaults; `strict_optional` stays on. `check_untyped_defs` costs 22 errors
  and makes sure that an unannotated new function is still checked.
  `platform` makes the result independent of the machine; the `os.name`
  branch in `losat_setup.py` costs 4 errors.
- **D3 Baseline shape: a count map of type debt, compared exactly.** A file's
  type debt is its mypy errors plus its `# type: ignore` comments.
  `TYPE_DEBT_BASELINE` in `tests/test_type_check_ratchet.py` maps each file
  with debt to that number. The guard fails when a file's debt is above its
  entry, a file that is not listed has debt (a new file included), a file
  with an entry has less debt than recorded ("lower its entry to N" or
  "remove its entry"), or an entry names a file that no longer exists.
  - Counting the comments closes the obvious bypass: a new error hidden by a
    new `# type: ignore` still raises the file's debt.
  - Exact comparison keeps the baseline equal to the code. Every fix lowers
    the entry in the same pull request, so a later change cannot spend the
    slack. The comparison can be exact only because the count is
    reproducible (D2).
  - The map is a literal in the guard, like `UNCHECKED_MODULES` in R14, so
    every change to it is a reviewed line in a diff. No tool or CI step
    rewrites it.
  - A count map starts protecting every file on the first day. The Web's
    "unchecked module" list would leave 138 files unprotected until their
    slice lands, while their errors are spread across the package.
- **D4 The guard runs in `Core PR` and `Core`.** It is an ordinary pytest
  test without a marker. Those two jobs install `.[export,typecheck]`
  instead of `.[export]`. This is the only workflow change, and it lands in
  T0. `Lint` was the alternative, but it cannot import the package (2.1).
  Running the guard through the existing `.[dev]` jobs would need a `recipe`
  or `browser` marker, which would mislabel it. The guard costs about 6 s per
  job and runs three times on a dev push, once per Python version.
- **D5 Third-party packages are `Any`.** `no_site_packages` hides installed
  packages, including `py.typed` packages and stub packages. T1 removes the
  352 comments that suppress nothing, including all 331 import comments; the
  60 live ones remain as debt and go with their slices. `pandas-stubs` is not
  added. It would have no effect under `no_site_packages`, and it changes the
  counts in opposite directions in the two checkers (2.2). Checking calls into
  numpy, pandas, and BioPython against pinned stubs is a candidate for the
  next period (Section 10).
- **D6 Scope.** `gbdraw/**/*.py` except `gbdraw/web/` (only
  `gbdraw/web/__init__.py` there; R14 covers the Web JavaScript). `tests/`
  and `tools/` are out of scope. Tests build ill-typed inputs on purpose to
  exercise validation, and use fixtures and monkeypatching that would need
  many casts. Tools are scripts outside the wheel, several of them tied to
  one issue (`characterize_issue597_s01.py`, `measure_issue597_s01.py`).
  Both import `gbdraw`, so editors still see the typed signatures. `tests/`
  is a candidate for the next period.
- **D7 Annotation-only changes.** A slice pull request changes typing only,
  apart from the fixes its body names. Allowed:
  - parameter, return, and variable annotations;
  - `typing.cast`;
  - imports from `typing`, `typing_extensions`, and `collections.abc`, and
    `from __future__ import annotations` (added to a file that lacks it
    before its annotations change);
  - `if TYPE_CHECKING:` blocks, in files that have the `__future__` import;
  - `Protocol` and `TypedDict` classes, `TypeVar`, and `@overload` stubs;
  - `: TypeAlias` added to an existing alias assignment;
  - removing a `# type: ignore` comment.

  These count as code changes and are named in the body: every class-body
  annotation (dataclass and `NamedTuple` fields, which `get_type_hints` reads
  in `gbdraw/config/modify.py` and `gbdraw/session_request_codec.py`), and a
  new type alias, which is module-level code. A name imported only under
  `TYPE_CHECKING` must not appear in a class-body annotation, because
  `get_type_hints` would raise `NameError`. The annotations of names in
  `gbdraw.__all__` are part of the public contract snapshot
  (`tests/test_public_contract.py`); changing one is a public API change
  (D10). Appendix A has the check; a slice pull request records its output.
- **D8 Narrowing `None` is a runtime change.** If the `None` can arrive from
  a public input (the Python API, the CLI, a Session, or a table file), it is
  a defect `PY-xx`: a failing test comes first, then a fix that raises the
  appropriate `GbdrawError` subclass. If the `None` cannot occur by
  construction, the first choice is to remove `None` from the declared type,
  which is annotation-only. Otherwise `assert x is not None` states the
  invariant (the package already uses 26 such asserts), and the body names
  it.
- **D9 A missing or different mypy fails the guard** with "run
  `pip install -e '.[typecheck]'`". The guard does not skip; a guard that
  skips does not protect anything.
- **D10 Design problems found while typing** (Owner instruction,
  2026-10-06). Examples: an argument nobody reads, one layer reading another
  layer's internals, the same value derived in two places, callers that
  disagree about whether `None` is returned, or an untyped `dict` used as a
  record.
  - A behavior-preserving fix of a few lines goes into the slice, named in
    the body as `F-xx`.
  - A larger one is its own STANDARD pull request: moving a layer, replacing
    a `dict` with a dataclass or `TypedDict`, or changing an API shape. If it
    is architecture-bearing, it follows
    `docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md`. A change to the
    shape of a public API (what `gbdraw` or `gbdraw.api` exports) is recorded
    in `CHANGELOG.md`; if it affects users, the recommendation goes to the
    Owner first.
  - A behavior change is a defect `PY-xx` with a failing test first. If SVG
    output changes, `tests/reference_outputs/` is regenerated only after the
    difference is reviewed (root `CLAUDE.md`, "Updating Reference Outputs").
- **D11 Rule changes.** The guard compares `[tool.mypy]` with
  `EXPECTED_MYPY_CONFIG` and rejects the suppressions that hide more than one
  counted line (4.2, assertion 3). Relaxing the
  configuration or raising an entry expands authority. It goes in its own
  pull request, which waits for the Owner's approval and is not
  auto-merged. Lowering or removing entries is a contraction and may ship
  with code changes.

## 4. Rule and guard (landed by T0)

### 4.1 Root `CLAUDE.md`

In "Quick Commands":

```bash
# Type-check gbdraw/ (pinned mypy from the typecheck extra) and compare with the baseline
pytest tests/test_type_check_ratchet.py
mypy
```

In "Coding Conventions", the Patterns line "`type: ignore` comments for
BioPython (missing type stubs)" is removed, and a new subsection follows
Patterns:

```markdown
### Type checking
- mypy, pinned in the `typecheck` extra and configured in `[tool.mypy]`,
  checks `gbdraw/` except `gbdraw/web/`. Third-party packages are `Any` to it,
  so their imports need no `# type: ignore`.
- A file's type debt is its mypy errors plus its `# type: ignore` comments.
  `TYPE_DEBT_BASELINE` in `tests/test_type_check_ratchet.py` only goes down:
  lower or remove a file's entry when its debt falls. A file without an
  entry, including a new file, has no debt. Fix the type, use
  `typing.cast`, or fix the code instead of adding `# type: ignore`. Raising
  an entry or relaxing `[tool.mypy]` needs the Owner's approval.
- Class-body annotations are runtime data. `get_type_hints` reads dataclass
  fields in `gbdraw/config/modify.py` (config overrides) and
  `gbdraw/session_request_codec.py` (canonical request validation), and
  `tests/test_public_contract.py` hashes the signatures of `gbdraw.__all__`.
  Treat a change to either as a code change, and use only names importable at
  run time in them, not names imported under `TYPE_CHECKING`.
```

### 4.2 Guard (`tests/test_type_check_ratchet.py`)

1. **Pinned version.** `[project.optional-dependencies] typecheck` holds
   exactly one `mypy==X` entry, and the installed mypy is version X.
   Otherwise: "run `pip install -e '.[typecheck]'`".
2. **Configuration.** `[tool.mypy]` equals `EXPECTED_MYPY_CONFIG`.
3. **No suppression outside the count.** Under `gbdraw/` (outside
   `gbdraw/web/`), no comment starts with `# mypy:` (per-file
   configuration), no `# type: ignore` stands on a line of its own (at the top
   of a module it silences the whole module), and no `no_type_check` appears
   (it silences a whole function).
4. **Debt.** Assertion 1 first, then one run of `sys.executable -m mypy --config-file pyproject.toml
   --cache-dir <os.devnull> --output json` from the repository root. Exit
   code 0, or 1 with JSON lines, is accepted; anything else fails with mypy's
   output. Errors outside the
   checked file set fail. `# type: ignore` comments are counted with
   `tokenize`, so a string that contains the text does not count. Then D3
   applies, and the keys of `TYPE_DEBT_BASELINE` must be sorted. A failure
   names the file, the debt split into errors and comments, the new error
   lines, and the change to make.

A cold run without a cache keeps the result independent of earlier runs; it
costs 3-7 s.

### 4.3 `pyproject.toml`

```toml
[project.optional-dependencies]
export = ["cairosvg"]
typecheck = ["mypy==2.4.0"]
dev = [
    "build",
    "cairosvg",
    "gbdraw[typecheck]",
    # ... unchanged
]

[tool.mypy]
python_version = "3.10"
platform = "linux"
files = ["gbdraw"]
exclude = ["^gbdraw/web/"]
no_site_packages = true
check_untyped_defs = true

[[tool.mypy.overrides]]
module = ["BCBio.*", "Bio.*", "fontTools.*", "numpy.*", "pandas.*", "svgwrite.*", "tomli"]
ignore_missing_imports = true
```

`pip install --dry-run --ignore-installed` confirms that `.[dev]` and
`.[export,typecheck]` both resolve `mypy 2.4.0`.

### 4.4 Workflow

In `.github/workflows/test.yml`, the "Install core test dependencies" step of
`core` and `core-pr` becomes:

```yaml
run: python -m pip install -e ".[export,typecheck]" pytest pytest-timeout pytest-xdist Pillow
```

## 5. PR sequence

Common rules for every slice pull request (S1a-S7):

- Change class STANDARD. "This is not architecture-bearing", unless a named
  F-xx fix moves an owner or a path, in which case that fix is its own pull
  request (D10).
- Fix the slice's errors with annotations and casts (D7), name every other
  change (D8, D10), and lower or remove the slice's entries in
  `TYPE_DEBT_BASELINE`. The target is that no module of the slice keeps an
  entry. An error that needs a larger design change keeps its entry, and the
  body names the F-xx pull request that will remove it. For S1a and S1b the
  target is required (completion condition 4).
- Verification:
  ```bash
  pip install -e ".[dev]"   # once per environment; brings mypy 2.4.0
  python -m pytest tests/test_type_check_ratchet.py
  mypy   # full list, for working through a file
  python <scratch>/annotation_only.py origin/dev HEAD   # Appendix A; expected: no runtime change, or exactly the named fixes
  python -m pytest tests/ -m "not slow and not (recipe or gallery or browser)" -n 4   # under the machine lock
  python -m pytest tests/test_output_comparison.py::TestOutputComparison   # when a named fix can reach rendering
  ruff check gbdraw/
  ```
- After merge, the next slice rebases onto `dev`. Two slices edit different
  lines of the same literal, so an adjacent-line conflict is resolved by
  keeping both sides' lower numbers and rerunning the guard.

### P0: this plan (STANDARD, documentation)

`docs/internal/PYTHON_TYPE_CHECK_RATCHET_PLAN_2026-10-06.md`.

### T0: checker, guard, and baseline (GOVERNANCE)

Files:

- `pyproject.toml`: the `typecheck` extra, `gbdraw[typecheck]` in `dev`,
  `[tool.mypy]` (4.3).
- `.github/workflows/test.yml`: the install line of `core` and `core-pr`
  (4.4).
- `tests/test_type_check_ratchet.py` (new): the guard (4.2), with
  `TYPE_DEBT_BASELINE` holding the debt of the base: 1,082 in 138 files at
  `c9c3b2e5`.
- Root `CLAUDE.md`: 4.1.
- `.gitignore`: `.mypy_cache/`.

No file under `gbdraw/` changes. The baseline records the untouched base and
relaxes nothing; from then on each entry can only go down. The R14 T0 (#844)
registered its starting list the same way.

Verification:

```bash
pip install -e ".[dev]"
python -m pytest tests/test_type_check_ratchet.py        # 4 passed
python -m pytest tests/ -m "not slow and not (recipe or gallery or browser)" -n 4
node tools/ci-impact.mjs classify --base origin/dev --head HEAD   # ci-only + packaging + tests-only
# Negative checks in a scratch commit, then drop it:
#   a new module with a type error -> "is above its baseline 0", with the error line
#   `# type: ignore[assignment]` added to gbdraw/exceptions.py -> "is above its baseline 0"
#   check_untyped_defs = false -> assertion 2 fails; assertion 4 lists files now below their entries
#   `# mypy: ignore-errors`, a top-of-module `# type: ignore`, and `@typing.no_type_check` -> assertion 3 names each
#   mypy uninstalled -> "run `pip install -e '.[typecheck]'`"
```

Proposed title: `Type-check gbdraw with a pinned mypy and a shrink-only debt baseline`

### T1: remove `# type: ignore` comments that suppress nothing (STANDARD)

Remove the 352 comments that `mypy --warn-unused-ignores` reports, including
all 331 `# type: ignore[reportMissingImports]` comments on third-party
imports. Comments are not in the AST, so Appendix A reports no change. The
debt goes from 1,082 in 138 files to 730 in 107. `warn_unused_ignores` is not
turned on: every comment already counts as debt.

### S1a-S7: slices (STANDARD)

| PR | Slice (2.4) | Debt after T1 | Notes |
| --- | --- | ---: | --- |
| S1a | Persisted contracts | 121 | `_TOP_LEVEL_FIELDS_V5` and `CURRENT_SESSION_TOP_LEVEL_FIELDS` keep their names and members unless the same pull request updates `tests/test_session_request_codec.py` and the Web typedefs. `CURRENT_SESSION_TOP_LEVEL_FIELDS` still contains `files`, which Session 40 and later reject and the writer no longer writes; removing it is F-xx in its own pull request. `config/models/` follows the runtime-annotation exception of D7. |
| S1b | Public API | 122 | `api/request_render.py` carries most of the `**dict[str, object]` bags; replacing them with typed parameters or a `TypedDict` is F-xx. Public signatures exported by `gbdraw.api` keep their shape unless the Owner agrees (D10). |
| S2 | Entry points | 62 | Examine the `plot_title_position` Literal mismatch and the batch-request union in `interface.py` (2.3). `losat_setup.py`: replace the `os.name` check with `sys.platform == "win32"`, a named fix that keeps Windows behavior. |
| S3 | Assembly | 110 | `height_below_final_record` (2.3) is F-xx. `no-redef` errors are fixed by annotating once before the branches. |
| S4 | Rendering | 79 | Any named fix that can reach SVG output runs `TestOutputComparison`. |
| S5 | Features and I/O | 99 | |
| S6 | Analysis and comparisons | 95 | `analysis/protein_colinearity.py` has 52. |
| S7 | Web support | 42 | These modules run in Pyodide for the Web app; a named fix there also runs the Web contracts locally. |

Order: S1a, S1b, then S2-S7 in order. A slice whose files no in-flight pull
request touches may go earlier.

### Final report to the Owner

The debt trajectory (T0, T1, after each slice), the PY-xx defects and F-xx
design findings with their pull requests, the Owner-delegated choices, and
the recommendation for the next period.

## 6. Interaction with in-flight work

| In-flight | Overlap | What waits |
| --- | --- | --- |
| Web typed boundaries, phase 2 (`docs/internal/WEB_TYPED_BOUNDARIES_PHASE2_PLAN_2026-10-06.md`) | The parity tests on `_TOP_LEVEL_FIELDS_V5` and `CURRENT_SESSION_TOP_LEVEL_FIELDS` | S1a does not change those sets (6.1). |
| Web owner-coupling Phase E (#857, #859) | none: Web JavaScript and Web authority only | — |
| Open Python pull requests | none on 2026-10-06 | — |
| Any Python pull request that merges after T0 | The guard applies: new files have no debt, and changed files cannot exceed their entries. A pull request opened before T0 is updated onto `dev` before it merges, so that the guard runs on it. | — |

### 6.1 Shared contract sets

S1a may annotate `_TOP_LEVEL_FIELDS_V5` and `CURRENT_SESSION_TOP_LEVEL_FIELDS`
(for example as `frozenset[str]`) but does not rename them or change their
members. A change to either is a separate pull request that also updates
`tests/test_session_request_codec.py` and the Web typedef that mirrors it.

## 7. Acceptance

- T0: the guard runs in `Core PR` and `Core` in under 15 s, the negative
  checks fail with the documented messages, and the baseline records 1,082
  in 138 files at the T0 base.
- T1: 730 in 107 files; Appendix A reports no runtime change.
- After S1a and S1b: no S1 module has an entry (completion condition 4).
- Every slice pull request passes Appendix A apart from the fixes its body
  names, and the fast suite passes. `tests/reference_outputs/` is unchanged
  unless a named PY-xx says otherwise.
- The final report gives the trajectory and the next-period recommendation.

## 8. Rollback

- A slice pull request is fixed forward. Reverting one would raise entries,
  which is an expansion that needs the Owner (D11).
- To retire the mechanism, revert T0 in a GOVERNANCE pull request: the guard,
  the configuration, the extra, the workflow line, and the CLAUDE.md text.
  The annotations stay, and they do no harm.
- If mypy 2.4.0 misbehaves, pin another version in one GOVERNANCE pull
  request with re-measured entries. Raised entries need the Owner's approval.

## 9. Owner decisions

None blocks. The Owner may want to override these delegated choices before T0:

1. **D1 checker.** mypy 2.4.0, which is reproducible without site-packages,
   installs with `pip`, and takes 3-7 s. The alternative is pyright 1.1.414,
   which matches Pylance and infers more, but needs Node, takes 13-14 s, and
   gives counts that depend on installed packages unless it is pointed at an
   empty virtual environment.
2. **D3 debt counts `# type: ignore`, compared exactly.** The alternative,
   errors only with "at most", is simpler, but it lets a suppression or an
   unrecorded fix pass silently.
3. **D4 the workflow change.** `core` and `core-pr` install the `typecheck`
   extra. The alternative is the `Lint` job, which would then have to
   install the package and pytest.
4. **D5 third-party packages as `Any`.** This gives up the BioPython types
   (94 errors that the checker would otherwise report) in exchange for counts
   that do not change when a library releases.

## 10. Non-goals and the next period

- Not in this plan: bringing the whole package to zero, stricter flags,
  typed third-party packages, `tests/`, `tools/`, and pyright or Pylance
  settings.
- The composition in 2.3 suggests three kinds of work after T1:
  annotation-only fixes (the 30 `no-redef` and 24 `var-annotated` errors and
  many of the 305 mismatches, such as `Sequence` against `list`), named small
  fixes (the 105 `None` and 57 other union narrowings), and F-xx design pull
  requests (most of the 149 errors from `object`-typed bags, concentrated in
  `api/request_render.py`, `interface.py`, and `render/interactive_svg.py`).
  The final report gives the actual split. If S1a-S7 meet their targets, the
  entries left are the ones tied to open F-xx pull requests.
- Next-period candidates, measured at `01694c3a` on top of 4.3:

  | Additional flag | Errors |
  | --- | ---: |
  | `disallow_incomplete_defs` | +126 |
  | `disallow_untyped_defs` | +179 |
  | `warn_return_any` | +45 |
  | `strict` | +963 |

  Unlike `tsc`, mypy accepts per-module options, so a stricter flag can
  start with the S1 modules through `[[tool.mypy.overrides]]` and spread
  slice by slice. A tightening raises entries only by the errors the new flag
  adds, so it is an authority change for the Owner, as in D11.
- Typed third-party packages: a second program with exact pins of numpy,
  `pandas-stubs`, and biopython, whose debt is counted separately.
- `tests/`, after the package reaches zero.

## Appendix A: annotation-only check

Save as `annotation_only.py` outside the repository, then run it in the
checkout: `python annotation_only.py origin/dev HEAD`. It compares the
normalized AST of every changed `gbdraw/**/*.py` file (outside `gbdraw/web/`)
and prints a unified diff of the normalized code of each file whose runtime
code changed. It exits 0 when there is none.

```python
"""Report runtime code changes in gbdraw/ between two revisions, ignoring typing.

Usage: python annotation_only.py [BASE] [HEAD]   (defaults: origin/dev HEAD; run in the checkout)

Both sides of every changed ``gbdraw/**/*.py`` file (except ``gbdraw/web/``) are
parsed and normalized, then compared with ``ast.dump``. Normalization removes
what has no runtime effect under ``from __future__ import annotations``:

- parameter and return annotations;
- annotations of assignments; a bare ``x: T`` inside a function body;
- ``cast(T, x)`` and ``typing.cast(T, x)`` (replaced by ``x``);
- ``if TYPE_CHECKING:`` blocks (the ``else`` branch is kept);
- imports from ``typing``, ``typing_extensions``, ``collections.abc``;
- classes deriving from ``Protocol`` or ``TypedDict``, ``@overload`` stubs,
  and ``T = TypeVar(...)``. ``X: TypeAlias = v`` counts as ``X = v``.

Normalization applies only when the HEAD side of the file has
``from __future__ import annotations``; without it, annotations are evaluated at
import time, so any typing change in that file is reported. Add the import first.

Class-body annotations always stay: they are runtime data. Dataclass fields
are read with ``get_type_hints`` by ``gbdraw/config/modify.py`` (config
overrides) and ``gbdraw/session_request_codec.py`` (canonical request
validation), and they decide ``ClassVar``/``InitVar``/``KW_ONLY`` and
``NamedTuple`` fields. Comments, including ``# type: ignore``, are not in the AST.
Exit 0 when every changed file is annotation-only; otherwise print a unified
diff of the normalized code per file and exit 1.
"""

from __future__ import annotations

import ast
import difflib
import subprocess
import sys

TYPING_MODULES = {"typing", "typing_extensions", "collections.abc", "__future__"}


def _name(node: ast.AST) -> str:
    if isinstance(node, ast.Name):
        return node.id
    if isinstance(node, ast.Attribute):
        return node.attr
    if isinstance(node, ast.Subscript):
        return _name(node.value)
    if isinstance(node, ast.Call):
        return _name(node.func)
    return ""


class Normalize(ast.NodeTransformer):
    def __init__(self, path: str) -> None:
        self.scopes: list[str] = ["module"]

    def _body(self, statements: list[ast.stmt]) -> list[ast.stmt]:
        out: list[ast.stmt] = []
        for statement in statements:
            result = self.visit(statement)
            if result is None:
                continue
            out.extend(result if isinstance(result, list) else [result])
        return out or [ast.Pass()]

    def _scoped(self, node, kind):
        self.scopes.append(kind)
        node.body = self._body(node.body)
        for field in ("orelse", "finalbody"):
            if getattr(node, field, None):
                setattr(node, field, self._body(getattr(node, field)))
        if hasattr(node, "handlers"):
            node.handlers = [self.visit(h) for h in node.handlers]
        self.scopes.pop()
        return node

    def visit_FunctionDef(self, node):
        if any(_name(d) == "overload" for d in node.decorator_list):
            return None
        node.decorator_list = [self.visit(d) for d in node.decorator_list]
        node.returns = None
        node.type_comment = None
        for arg in [*node.args.posonlyargs, *node.args.args, *node.args.kwonlyargs,
                    node.args.vararg, node.args.kwarg]:
            if arg is not None:
                arg.annotation = None
                arg.type_comment = None
        node.args.defaults = [self.visit(d) for d in node.args.defaults]
        node.args.kw_defaults = [d if d is None else self.visit(d) for d in node.args.kw_defaults]
        return self._scoped(node, "function")

    visit_AsyncFunctionDef = visit_FunctionDef

    def visit_ClassDef(self, node):
        if any(_name(b) in {"Protocol", "TypedDict"} for b in node.bases):
            return None
        node.bases = [self.visit(b) for b in node.bases]
        node.keywords = [self.visit(k) for k in node.keywords]
        node.decorator_list = [self.visit(d) for d in node.decorator_list]
        return self._scoped(node, "class")

    def visit_If(self, node):
        if _name(node.test) == "TYPE_CHECKING":
            return self._body(node.orelse) if node.orelse else None
        node.test = self.visit(node.test)
        node.body = self._body(node.body)
        node.orelse = self._body(node.orelse) if node.orelse else []
        return node

    def visit_For(self, node):
        node.target = self.visit(node.target)
        node.iter = self.visit(node.iter)
        return self._scoped(node, self.scopes[-1])

    visit_AsyncFor = visit_For

    def visit_While(self, node):
        node.test = self.visit(node.test)
        return self._scoped(node, self.scopes[-1])

    def visit_With(self, node):
        node.items = [self.visit(i) for i in node.items]
        return self._scoped(node, self.scopes[-1])

    visit_AsyncWith = visit_With

    def visit_Try(self, node):
        return self._scoped(node, self.scopes[-1])

    def visit_ExceptHandler(self, node):
        if node.type is not None:
            node.type = self.visit(node.type)
        node.body = self._body(node.body)
        return node

    def visit_ImportFrom(self, node):
        return None if node.module in TYPING_MODULES else node

    def visit_Import(self, node):
        node.names = [a for a in node.names if a.name not in TYPING_MODULES]
        return node if node.names else None

    def visit_AnnAssign(self, node):
        if self.scopes[-1] == "class":
            node.value = None if node.value is None else self.visit(node.value)
            return node
        if node.value is None:
            return None
        return ast.Assign(targets=[node.target], value=self.visit(node.value), lineno=0)

    def visit_Assign(self, node):
        if isinstance(node.value, ast.Call) and _name(node.value.func) in {"TypeVar", "ParamSpec", "TypeVarTuple"}:
            return None
        return self.generic_visit(node)

    def visit_Call(self, node):
        if _name(node.func) == "cast" and len(node.args) == 2:
            return self.visit(node.args[1])
        return self.generic_visit(node)


def _has_future_annotations(source: str | None) -> bool:
    return source is not None and any(
        isinstance(node, ast.ImportFrom) and node.module == "__future__"
        and any(alias.name == "annotations" for alias in node.names)
        for node in ast.parse(source).body
    )


def _normalized(path: str, source: str | None, strip: bool) -> str:
    if source is None:
        return ""
    tree = ast.parse(source)
    if strip:
        tree = Normalize(path).visit(tree)
        ast.fix_missing_locations(tree)
    return ast.unparse(tree)


def _show(revision: str, path: str) -> str | None:
    result = subprocess.run(["git", "show", f"{revision}:{path}"], capture_output=True, text=True)
    return result.stdout if result.returncode == 0 else None


def main(base: str = "origin/dev", head: str = "HEAD") -> int:
    paths = subprocess.run(
        ["git", "diff", "--name-only", base, head, "--", "gbdraw"],
        capture_output=True, text=True, check=True,
    ).stdout.split()
    paths = [p for p in paths if p.endswith(".py") and not p.startswith("gbdraw/web/")]
    changed = 0
    for path in paths:
        head_source = _show(head, path)
        strip = _has_future_annotations(head_source)
        before = _normalized(path, _show(base, path), strip)
        after = _normalized(path, head_source, strip)
        if before != after:
            changed += 1
            sys.stdout.writelines(difflib.unified_diff(
                before.splitlines(True), after.splitlines(True),
                f"{base}:{path}", f"{head}:{path}", n=2))
            print()
    print(f"annotation-only: {len(paths) - changed} of {len(paths)} changed files; runtime code changed: {changed}")
    return 1 if changed else 0


if __name__ == "__main__":
    sys.exit(main(*sys.argv[1:]))
```
