---
name: maintain-python-api
description: Maintain and extend gbdraw's public Python API while preserving contract snapshots, CLI behavior, rendering outputs, and error semantics. Use when modifying gbdraw/api, public exports or signatures, DiagramOptions or other public dataclasses, interactive SVG builders, table readers, render/export APIs, high-level diagram builders, Python API documentation, or related contract and regression tests.
---

# Maintain the Python API

Keep `gbdraw.api` thin, truthful, and independent of CLI parsing. Reuse the real
implementation owner, make failures observable to library callers, and prove that
public contracts and default diagrams remain stable.

## Inspect the contract first

1. Read `CLAUDE.md`, `docs/PYTHON_API_IMPROVEMENT_PLAN.md`,
   `docs/PYTHON_API.md`, and any related ADR.
2. Inspect `git status` before editing. Preserve unrelated and pre-existing work.
3. Inventory the affected exports, signatures, dataclass fields/defaults, callers,
   documentation examples, and contract snapshots.
4. Trace each public option through the adapter and assembler to the function that
   actually consumes it. Do not assume successful validation means successful
   forwarding.

## Preserve dependency direction

- Put shared behavior in a CLI-independent owner and make both CLI and API call it.
- Keep `gbdraw.api` to re-exports and small validation or translation adapters.
- Do not import CLI parsers or private CLI helpers from the public API.
- Reuse existing table validators, renderers, and config models instead of copying
  their logic.
- Add a public option class only when it reduces the combined caller and
  implementation complexity. Keep low-level `assemble_*` functions as compatibility
  and advanced-use entry points.

## Apply repository-specific patterns

### Options and typed config

- Normalize an option once, then test that the normalized value reaches its real
  downstream owner. Include aliases and the unchanged default in forwarding tests.
- Treat `LabelsFilteringConfig.raw` as payload, not ordinary nested dataclass state.
  A naive `dataclasses.asdict()` adds another `raw` layer and can hide label
  DataFrames from downstream code. Preserve `filtering.as_dict()` losslessly when
  applying typed config overrides.
- Accept either an in-memory DataFrame or a file path for a logical table, never
  both. Raise `ValidationError` for ambiguous input.
- Compare DataFrame contents and downstream behavior, not object identity.

### Render and export

- Keep the library API strict even if the CLI intentionally warns and skips an
  unavailable converter.
- Return only paths that exist when the function returns. Wrap conversion failures
  in a typed `GbdrawError` subclass and preserve the original exception as cause.
- Test stale outputs and overwrite behavior so an old file cannot make a failed
  conversion look successful.
- Keep browser/Pyodide conversion distinct from successful local file creation.

### Interactive SVG

- Own context construction in a CLI-independent module.
- Re-export the context type, enrichment function, and builder through `gbdraw.api`.
- Add context as an optional keyword-only argument to byte/file rendering APIs.
- Verify context-free interactive output still works and static SVG bytes do not
  change.
- Compare CLI and API metadata schema for feature and match popups.

### Public surface and sessions

- Update `gbdraw.api.__all__` and `tests/fixtures/public_contract.json` together.
  Build the new contract with the test helper, then review and apply only intended
  changes; never accept a broad snapshot rewrite blindly.
- Re-export public table models and readers from `gbdraw.api.io` while sharing the
  existing CLI-compatible validation owner.
- Do not expose session replay as CLI argument strings. Require CLI-independent
  typed request models and pure conversion first; otherwise record the non-public
  decision in an ADR.

## Keep documentation executable

- Add capability-oriented recipes to `docs/PYTHON_API.md` using public imports.
- Make optional-tool examples deterministic with fixtures or precomputed data.
- Remember that `tests/test_api_library_usage.py` executes every Python code block
  sequentially in one namespace. Make each block valid in that execution model.
- Explain intentional low-level or non-public boundaries instead of implying that
  an unstable workflow is supported.

## Verify without contaminating evidence

1. Run focused regression and contract tests for the changed boundary.
2. Run the targeted command in `docs/PYTHON_API_IMPROVEMENT_PLAN.md`.
3. Run `python -m pytest tests/ -v -m "not slow"`.
4. Lint touched Python files first, then run the repository lint gate. Distinguish
   new failures from a recorded pre-existing baseline.
5. Check the working tree and public contract diff after every broad validation run.

Treat reference generation as a write operation. `TestGenerateReferences` may
rewrite tracked SVG fixtures before comparison, which can invalidate the regression
evidence. Record the initial status of `tests/reference_outputs/`, prefer
comparison-only tests, and run generation only when intentionally refreshing an
approved geometry change. Never discard a reference file that was dirty before the
test run, and never mechanically accept reference changes for a non-geometry fix.

## Finish the change

- Confirm every new public symbol imports from `gbdraw.api`.
- Confirm defaults retain previous behavior and non-default values reach the owner.
- Confirm error types and returned artifacts match actual outcomes.
- Confirm CLI behavior remains intentional and uses the shared owner.
- Confirm docs, ADRs, contract snapshots, and tests describe the implemented state.
- Report targeted/full test counts, lint scope, reference-output status, and any
  pre-existing failure separately.
