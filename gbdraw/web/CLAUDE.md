# Web application maintenance

Read the repository-level `CLAUDE.md` first. This file describes the intended
Web architecture and its ownership boundaries. It deliberately avoids line
numbers and implementation-size snapshots because those become stale quickly.

## Non-negotiable properties

- The Web UI is a single-page app with no JavaScript build step.
- `index.html` owns markup, local styles, and templates. ES modules live under
  `gbdraw/web/js/`.
- All runtime assets are served from the same origin. Do not add a CDN or other
  runtime network dependency.
- Genome inputs and rendered results stay in the browser. Generation must not
  upload user data.
- Generated SVG is sanitized before insertion into the preview.
- The browser wheel is generated and gitignored. Prepare it when packaging or
  wheel-dependent tests need it; never edit it by hand.

## Runtime data flow

The canonical generation path is:

```text
index.html
  -> app.js and app/app-setup.js
  -> reactive state
  -> services/session-request.js builds one canonical schema-5 render request
  -> app/run-analysis.js validates and orchestrates optional LOSAT work
  -> services/diagram-generation.js dispatches the request and resources
  -> workers/diagram-generation-worker.js owns Pyodide rendering
  -> typed request decoding and render_request()
  -> sanitized preview, interactive editing, export, and session save
```

Do not introduce a second argv-shaped generation contract. Fresh generation,
saved sessions, and replay must converge on the typed render-request boundary.
JavaScript owns browser state and resource bytes; Python owns request validation,
planning, loading, diagram assembly, and rendering.

The application may have two Pyodide runtimes:

- The module worker owns diagram rendering and generation-time feature
  extraction. Heavy rendering must not return to the main thread.
- The main-thread runtime may serve small UI helpers that have not yet moved,
  such as palette data or isolated compatibility helpers. It is not an
  alternative render path.

LOSAT workers are separate from the diagram-generation worker. Threaded LOSAT
is available only when the page is cross-origin isolated; keep a single-thread
fallback.

## Module ownership

| Area | Owner |
|---|---|
| Vue mount and exports | `js/app.js` |
| Composition and dependency wiring | `js/app/app-setup.js` |
| Reactive state and computed values | `js/state.js` |
| Generate-button orchestration | `js/app/run-analysis.js` |
| Canonical request/session projection | `js/services/session-request.js` |
| Save/load coordination | `js/services/config.js` |
| Render-worker client and lifecycle | `js/services/diagram-generation.js` |
| Pyodide typed rendering | `js/workers/diagram-generation-worker.js` |
| LOSAT dispatch and workers | `js/services/losat.js`, `js/workers/losat-*` |
| SVG/PNG/PDF downloads | `js/services/export.js` |
| Preview SVG sanitization | `js/state.js` |
| Feature-editor entry point | `js/app/feature-editor.js` |
| Feature-editor helpers | `js/app/feature-editor/` |
| Legend entry point | `js/app/legend.js` |
| Legend helpers | `js/app/legend/` |
| Legend/diagram positioning | `js/app/legend-layout.js`, `js/app/legend-layout/` |
| Gallery data, media, and tutorials | `gallery/` |

Keep top-level `create*` entry points in `js/app/*.js`. Put larger,
single-purpose helpers in the matching subfolder instead of growing another
general utility module.

## Request and session boundary

`services/session-request.js` is the projection boundary between reactive UI
state and the persisted/rendered model. It owns:

- canonical request schema and resource descriptors;
- Circular `single`, `grid`, and `batch` grouping;
- Linear record/group topology;
- output-prefix and per-batch output projection;
- explicit track-slot projection;
- the current `ui.layoutPreferences` representation.

`services/config.js` coordinates user-facing save and load actions. It must not
grow a parallel model of render fields. Compatibility migrations are reader
concerns: normalize supported old data once, then use the current model. Do not
write retired fields into new sessions.

Explicit track slots are authoritative when enabled. Legacy flat controls are
compatibility inputs, not a second source of truth. Preserve empty positions in
per-record depth and comparison inputs when those positions carry alignment
meaning.

See `docs/SESSION_COMPATIBILITY.md` for accepted versions and migration limits.

## Modes and output topology

- Circular `single` renders one record.
- Circular `grid` renders several records in one figure; a one-record grid is
  valid.
- Circular `batch` renders one output per record and has one resolved output
  target per item.
- Linear renders one ordered multi-record layout. Record groups and explicit
  comparison pairs affect topology and must survive session round trips.

The preview and download layers consume render results; they must not infer a
different grouping from the number of uploaded files.

## Assets, CSP, and privacy

Vue, Pyodide, styling, fonts, icons, DOMPurify, and export libraries are vendored
under `gbdraw/web/vendor/` or otherwise packaged locally. Keep the Content
Security Policy in `index.html` aligned with actual same-origin needs. Adding an
external host to the CSP is not a substitute for vendoring a dependency.

Do not log genome sequence, full uploaded file contents, or generated
comparison rows. Error messages should identify the failed stage and safe
resource label without exposing private data.

Generated SVG must pass through the shared sanitization profile before preview
or interactive editing. Event-handler attributes, scripts, foreign content,
and unsafe URL schemes remain forbidden.

## Changing a setting

Before adding or changing a control, trace its complete lifecycle:

1. Declare its reactive state and default in the appropriate state/setup owner.
2. Bind the visible control and its help text.
3. Project it once into the canonical request.
4. Decode and validate it in the typed Python request layer if it is new.
5. Preserve it through session save/load when it is user-owned state.
6. Add focused state/request/session tests and a browser assertion when visual
   behavior changes.
7. Update Gallery tutorial captures when the control is part of a documented
   workflow.

Do not add a diagram argv builder. Project each control once at the canonical
request boundary instead of duplicating it in configuration and session modules.

## Gallery ownership

`gbdraw/web/gallery/index.html` is generated from
`gbdraw/web/gallery/examples.json`. Tutorial instructions live under
`gbdraw/web/gallery/tutorials/`; screenshots and thumbnails live under
`gbdraw/web/gallery/media/` and `thumbnails/`.

Use `tools/build_web_gallery.py` for generated gallery markup. Use
`tools/capture_gallery_tutorial_screenshots.py` for declarative tutorial
captures. Data-dependent captures must load the example's own session, declare
the expected app state, and prove that the controls or data identity named by
the instruction are visible in the final crop.

Read `docs/skills/web-gallery-screenshot-maintenance/SKILL.md` before editing
Gallery tutorials or screenshots.

## Local build and verification

```bash
gbdraw gui
python tools/prepare_browser_wheel.py
python tools/prepare_browser_wheel.py --refresh-cache-bust
python -m build
pytest tests/ -v -m "not slow"
```

Refresh the cache-bust token only when preparing a deployable bundle.

When browser verification matters, check both Playwright installations:

```bash
command -v playwright && playwright --version
python -c "from playwright.sync_api import sync_playwright; print('python playwright ok')"
node -e "console.log(require.resolve('@playwright/test'))"
```

The JavaScript specs under `tests/web/` require Node's `@playwright/test`. If it
is unavailable, use Python Playwright for focused browser checks. In an agent
sandbox, Chromium may fail with `sandbox_host_linux.cc ... Operation not
permitted`; rerun the same check with the required sandbox escalation.

For offline packaging, test that the wheel version matches `pyproject.toml`,
all runtime URLs are local, CSP allows the required worker/runtime behavior, and
the app reaches a ready state without external network access.

## Debugging principles

- Treat request validation failures as boundary errors, not reasons to bypass
  the typed request path.
- Keep worker errors structured enough for the main thread to distinguish
  initialization, resource transfer, validation, LOSAT, render, and export
  failures.
- Revoke object URLs and terminate superseded workers so repeated generation
  does not leak browser memory.
- Preserve cancellation and stale-result guards when changing async code.
- Verify both Circular and Linear modes after shared state or request changes.

Related project documentation:

- `CLAUDE.md`
- `docs/TYPED_API.md`
- `docs/SESSION_COMPATIBILITY.md`
- `docs/SVG_SEMANTIC_HOOKS.md`
