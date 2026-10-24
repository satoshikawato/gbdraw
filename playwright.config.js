// @ts-check
const { resolve } = require('node:path');
const { defineConfig, devices } = require('@playwright/test');
const { PYODIDE_SNAPSHOT_TEST_PATH } = require('./tests/web/helpers/pyodide-snapshot-test-server.cjs');
const serverPort = Number(process.env.GBDRAW_WEB_TEST_PORT || 4173);
const serverUrl = `http://127.0.0.1:${serverPort}`;
// A directory that tools/build_pyodide_snapshot.mjs wrote: the server then also
// serves it, and the diagram Worker restores it instead of installing its
// packages (docs/internal/SELECTIVE_CI.md, "Pyodide snapshot").
const pyodideSnapshotDir = process.env.GBDRAW_TEST_PYODIDE_SNAPSHOT_DIR;

// The variant configurations in tests/web/playwright/ spread this object, and
// Playwright resolves relative paths against the loaded config file, so these
// paths are absolute.
module.exports = defineConfig({
  testDir: resolve(__dirname, 'tests', 'web'),
  outputDir: resolve(__dirname, 'test-results'),
  testMatch: /.*\.playwright\.spec\.js/,
  fullyParallel: true,
  forbidOnly: Boolean(process.env.CI),
  retries: process.env.CI ? 2 : 0,
  workers: process.env.CI ? 1 : undefined,
  reporter: process.env.CI ? 'github' : 'list',
  // In CI Playwright reads the pull request's whole diff for the HTML report,
  // which no reporter here uses; a diff over Node's string limit (a Gallery
  // regeneration) crashed the runner before any test ran (OV-1742).
  captureGitInfo: { diff: false },
  use: {
    baseURL: serverUrl,
    trace: 'on-first-retry'
  },
  webServer: {
    command: pyodideSnapshotDir
      ? `python3 tools/serve_pyodide_snapshot.py --port ${serverPort} --snapshot-dir ${JSON.stringify(resolve(pyodideSnapshotDir))} --snapshot-path ${PYODIDE_SNAPSHOT_TEST_PATH}`
      : `python3 -m http.server ${serverPort} --bind 127.0.0.1`,
    cwd: __dirname,
    url: serverUrl,
    reuseExistingServer: !process.env.CI,
    stderr: 'ignore'
  },
  projects: [
    {
      name: 'chromium',
      use: { ...devices['Desktop Chrome'] }
    }
  ]
});
