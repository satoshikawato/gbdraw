# Web Pyodide Worker Cancellation Roadmap

## Goal

Make web generation cancellation real.

The UI must remain responsive while gbdraw is assembling and rendering a diagram,
and a cancelled run must stop promptly instead of waiting for a synchronous
main-thread Pyodide call to return.

## Current State

The current progress and cancel MVP already handles much of the LOSAT phase:

- `generationProgress` exists in `gbdraw/web/js/state.js`.
- The processing overlay shows phases, LOSAT pair progress, and a Cancel button.
- `run-analysis.js` creates an `AbortController` for manual generation runs.
- `services/losat.js` accepts an `AbortSignal` and terminates the LOSAT Worker
  pool on cancellation.
- Completed LOSAT pair results are cached immediately from `onCompleted`.

The remaining weak points are:

- `run_gbdraw_wrapper(...)` still runs synchronously on the main thread through
  Pyodide. While this call is active, the Vue UI cannot reliably process a
  Cancel click.
- The sequential LOSAT fallback path can still run a pair on the main thread.
  If that path is used, an active pair cannot be interrupted promptly.
- Pyodide is used in several places, not only for final drawing. Existing direct
  calls include palette loading, record listing, FASTA extraction, feature
  extraction, legend entry generation, and definition rerendering.

## Target Architecture

Keep the UI on the browser main thread. Move Python and gbdraw execution into a
dedicated Pyodide Worker.

```text
main thread
  Vue UI
  generation overlay
  Cancel button
  SVG preview DOM
  LOSAT cache state
  LOSAT Worker pool
      |
      | postMessage
      v
pyodide-gbdraw-worker
  Pyodide runtime
  gbdraw wheel
  Python helper functions
  run_gbdraw_wrapper()
  feature extraction helpers
```

The LOSAT Worker pool should stay separate from the Pyodide drawing Worker for
the first implementation. The main thread can continue to build LOSAT jobs,
cache completed TSV output, and pass final virtual BLAST file contents to the
Pyodide Worker when drawing starts.

Avoid nested Workers in the first pass. They add browser compatibility and asset
resolution complexity without being necessary for responsive cancellation.

## Cancellation Model

Use two cancellation levels.

### Soft Cancel

When available, use Pyodide interrupts:

- Pyodide runs inside a Worker.
- The main thread creates a `SharedArrayBuffer`.
- The Worker calls `pyodide.setInterruptBuffer(interruptBuffer)`.
- On Cancel, the main thread writes `2` into the buffer, which is SIGINT.
- Python should receive `KeyboardInterrupt`.

This requires a cross-origin isolated page. In practice that means the hosted app
and `gbdraw gui` server need compatible HTTP headers such as COOP and COEP.

Reference:
https://pyodide.org/en/latest/usage/keyboard-interrupts.html

### Hard Cancel

If soft cancel is unavailable or does not stop the run within a short grace
period, terminate the Pyodide Worker.

Hard cancel is the reliable fallback. It discards the current Pyodide runtime,
so the app must recreate and reinitialize the Worker before the next Python
operation.

Suggested behavior:

- Cancel clicked during drawing.
- If `SharedArrayBuffer` is available, request SIGINT.
- Wait 2 to 5 seconds.
- If the run has not completed or aborted, terminate the Worker.
- Show a runtime restart message while the Worker is recreated.

## Worker Ownership

### Main Thread Owns

- Vue state and DOM.
- Uploaded `File` objects and conversion to transferable buffers.
- Generation progress state.
- Cancel button state.
- LOSAT cache and LOSAT cache info.
- LOSAT Worker pool.
- SVG preview assignment and sanitization.
- Export actions.

### Pyodide Worker Owns

- Loading Pyodide.
- Installing local Python wheels and gbdraw.
- Running `PYTHON_HELPERS`.
- Writing run-specific files into Pyodide FS.
- Calling Python helper functions.
- Calling `run_gbdraw_wrapper()`.
- Returning structured results and structured errors.

The Worker must not reach into DOM state. Every action should be expressed as a
message with explicit input and output.

## Message Protocol

Use request IDs for every command. Ignore stale responses on the main thread.

Suggested command types:

- `init`
- `getPalettes`
- `listGenbankRecords`
- `extractFirstFasta`
- `extractFeatures`
- `generateLegendEntrySvg`
- `regenerateDefinitionSvgs`
- `runGbdraw`
- `interrupt`
- `dispose`

Suggested response event types:

- `ready`
- `status`
- `progress`
- `result`
- `error`
- `cancelled`

Errors should be structured:

```javascript
{
  type: "PythonError" | "JsError" | "AbortError",
  summary: "...",
  details: [
    { label: "STDERR", text: "..." },
    { label: "STDOUT", text: "..." },
    { label: "Traceback", text: "..." }
  ]
}
```

## Stepwise Roadmap

### Phase 0: Record the Current Contract

Deliverables:

- Update the existing generation progress checklist so it reflects what is
  already implemented.
- Add a short manual test matrix for the current behavior.
- Confirm the current app still runs after the existing progress/cancel MVP.

Acceptance criteria:

- The project has one clear source of truth for current progress/cancel status.
- Known limitations are documented as limitations, not TODO ambiguity.

### Phase 1: Add a Python Runner Interface

Introduce a main-thread service interface that hides direct Pyodide access.

Initial implementation can still call the current main-thread Pyodide runtime.
The point of this phase is to make call sites depend on a runner API instead of
raw `getPyodide()`, `runPython()`, `globals.get(...)`, and `FS.writeFile(...)`.

Candidate API:

```javascript
pythonRunner.init()
pythonRunner.getPalettes()
pythonRunner.listGenbankRecords({ file })
pythonRunner.extractFirstFasta({ file, selector, region, reverse })
pythonRunner.extractFeatures({ files, mode, selectors, regions, reverseFlags })
pythonRunner.generateLegendEntrySvg(payload)
pythonRunner.regenerateDefinitionSvgs(payload)
pythonRunner.runGbdraw(payload, { signal })
pythonRunner.cancelRun(requestId)
pythonRunner.dispose()
```

Call sites to migrate behind the interface:

- `gbdraw/web/js/app/pyodide.js`
- `gbdraw/web/js/app/run-analysis.js`
- `gbdraw/web/js/app/results.js`
- `gbdraw/web/js/app/legend/entry-actions.js`

Acceptance criteria:

- Direct Pyodide access is localized to one adapter module.
- No user-visible behavior changes are required in this phase.
- Existing offline asset checks still pass.

### Phase 2: Remove Uncancellable LOSAT Fallback from the UI Path

Keep the sequential runner for non-UI fallback or debugging, but do not use it
for browser generation unless the user explicitly accepts the non-cancellable
behavior.

Implementation options:

- Add `allowSequentialFallback: false` when `run-analysis.js` calls
  `runLosatPairsParallel(...)`.
- If Worker construction or Worker execution setup fails in the UI path, surface
  a structured error instead of running LOSAT on the main thread.
- Keep abort checks before and after loading LOSAT assets.

Acceptance criteria:

- A browser LOSAT run never starts a long active pair on the main thread by
  default.
- Cancelling a LOSAT run either terminates active Workers or stops before any
  pair begins.
- Completed pair cache behavior remains unchanged.

### Phase 3: Add the Pyodide Worker Shell

Add a dedicated Worker file, for example:

- `gbdraw/web/js/workers/pyodide-gbdraw-worker.js`
- `gbdraw/web/js/services/pyodide-runner.js`

The Worker should support `init` first:

- Load the local Pyodide runtime.
- Load micropip.
- Install local dependency wheels.
- Install the local gbdraw wheel.
- Run `PYTHON_HELPERS`.
- Send status events that can replace the current `loadingStatus` updates.

The main-thread runner service should create the Worker, send `init`, track
readiness, and recreate the Worker after hard cancel.

Acceptance criteria:

- Pyodide can initialize inside the Worker.
- Startup status still appears in the UI.
- The main thread no longer needs to own the Pyodide object for initialization.
- `tools/verify_gui_offline.py check-assets` includes the new Worker file.

### Phase 4: Move Final Drawing into the Worker

Move the final `run_gbdraw_wrapper(...)` execution behind
`pythonRunner.runGbdraw(...)`.

Payload should contain:

- mode
- CLI args
- files as transferable buffers or text
- virtual BLAST files from completed LOSAT results
- any generated TSV contents currently written into Pyodide FS

The Worker should:

- Reset or overwrite run-specific FS paths.
- Write incoming files to Pyodide FS.
- Call `run_gbdraw_wrapper(...)`.
- Return the JSON result string parsed as structured JS data, or return a
  structured error.

Main thread behavior:

- Set phase to `Rendering diagram`.
- Keep the overlay interactive.
- On Cancel, hard terminate the Worker if soft cancel is not available yet.
- Do not assign returned results if the request ID is stale or cancelled.

Acceptance criteria:

- Clicking Cancel during drawing is processed immediately by the UI.
- A cancelled drawing run does not assign SVG results.
- A hard-cancelled Worker is recreated before the next Python operation.
- Normal drawing output is unchanged.

### Phase 5: Add Soft Cancel with Pyodide Interrupts

Add optional interrupt support.

Main thread:

- Detect `globalThis.crossOriginIsolated`.
- Detect `typeof SharedArrayBuffer === "function"`.
- Create an interrupt buffer when supported.
- Send it to the Pyodide Worker during initialization.
- On Cancel, write SIGINT value `2` into the buffer.

Worker:

- Call `pyodide.setInterruptBuffer(interruptBuffer)`.
- Map Python `KeyboardInterrupt` to a structured `AbortError`.
- Clear or replace the interrupt buffer before each run to avoid stale signals.

Server and hosting:

- Add or verify COOP and COEP headers for `gbdraw gui`.
- Add or verify equivalent headers in deployment config.
- Confirm CDN and local assets are compatible with cross-origin isolation.

Acceptance criteria:

- In cross-origin isolated mode, many Python drawing cancellations finish via
  `KeyboardInterrupt` without killing the Worker.
- If soft cancel does not complete within the grace period, hard cancel still
  terminates the Worker.
- In non-isolated mode, cancellation still works through hard cancel.

### Phase 6: Move Remaining Python Helpers Behind the Worker

After final drawing works in the Worker, migrate the other Python helper calls.

Targets:

- Palette loading.
- GenBank record listing.
- FASTA extraction fallback.
- Feature extraction for the color editor.
- Dynamic legend entry SVG generation.
- Linear definition SVG rerendering.

Acceptance criteria:

- Main-thread modules no longer call Pyodide directly.
- The only Python access path is the runner service.
- Worker restart after hard cancel does not break subsequent helper operations.

### Phase 7: Add Regression and Browser Verification

Automated checks:

- Offline asset check includes new Worker assets.
- Smoke test confirms Pyodide Worker initialization.
- Smoke test generates a normal circular diagram.
- Smoke test generates a normal linear diagram with existing BLAST input.
- LOSAT cancellation test with three or more sequences:
  - cancel after at least one pair completes
  - completed TSV remains cached
  - unfinished TSV is not cached
  - rerun reuses completed cache hits
- Drawing cancellation test:
  - start a heavy render
  - click Cancel during `Rendering diagram`
  - assert overlay exits
  - assert no stale SVG result is assigned
  - assert the next run can initialize or reuse a healthy Worker

Manual checks:

- Same tests in a browser with `SharedArrayBuffer` available.
- Same tests in a browser or server mode without cross-origin isolation.
- Confirm export actions still work after a hard-cancelled run.

## Implementation Notes

- Prefer transferring `ArrayBuffer` objects to the Worker instead of sending
  `File` objects and re-reading them there.
- Keep request IDs on both sides. Cancellation and stale response handling should
  be explicit.
- Do not make fine-grained Python drawing progress part of the first Worker
  migration. Coarse phase progress is enough until Python-side hooks are added.
- Do not let Worker restart erase main-thread LOSAT cache state.
- Treat Worker termination as a normal cancellation path, not an exceptional app
  crash.
- Keep the direct main-thread Pyodide adapter during migration until every call
  site has moved behind the runner interface.

## Open Decisions

- Whether `gbdraw gui` should always serve COOP/COEP headers, or only enable
  them when local assets are guaranteed to be compatible.
- The hard-cancel grace period after a soft cancel request.
- Whether Worker initialization should happen eagerly on page load or lazily on
  first Python operation after the UI is mounted.
- Whether the sequential LOSAT fallback should be kept behind a developer flag
  only, or exposed to users with a clear warning.

