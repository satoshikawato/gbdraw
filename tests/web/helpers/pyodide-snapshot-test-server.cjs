// The Pyodide snapshot that a test server may serve to the diagram Worker
// (CD-17, CD-18, CI cost 2026-10-16). The Worker restores the snapshot when
// one is served and takes its cold boot (Pyodide, micropip, the gbdraw wheel)
// on any restore error. The test server and the Playwright config import the
// path from here.

/** The URL prefix under which a test server serves a Pyodide snapshot. */
const PYODIDE_SNAPSHOT_TEST_PATH = '/__pyodide-snapshot__/';

/**
 * Makes the page's diagram Worker take its cold boot: every request under
 * {@link PYODIDE_SNAPSHOT_TEST_PATH} fails, so no snapshot can be restored.
 * A case that asserts the cold boot (the wheel request, its failure) calls
 * this before the app boots the Worker; a Worker that already booted keeps
 * its runtime.
 * @param {import('@playwright/test').Page} page
 * @returns {Promise<void>}
 */
const forceColdPyodideBoot = (page) => page.route(
  (url) => url.pathname.startsWith(PYODIDE_SNAPSHOT_TEST_PATH),
  (route) => route.abort()
);

module.exports = { PYODIDE_SNAPSHOT_TEST_PATH, forceColdPyodideBoot };
