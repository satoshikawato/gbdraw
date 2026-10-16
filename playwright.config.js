// @ts-check
const { resolve } = require('node:path');
const { defineConfig, devices } = require('@playwright/test');
const serverPort = Number(process.env.GBDRAW_WEB_TEST_PORT || 4173);
const serverUrl = `http://127.0.0.1:${serverPort}`;

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
  use: {
    baseURL: serverUrl,
    trace: 'on-first-retry'
  },
  webServer: {
    // PROTOTYPE (proto/pyodide-snapshot, CI measurement only): with
    // GBDRAW_TEST_PYODIDE_SNAPSHOT_DIR set, the server also serves that
    // snapshot and config.js points the Worker at it.
    command: process.env.GBDRAW_TEST_PYODIDE_SNAPSHOT_DIR
      ? `python3 tools/proto-pyodide-snapshot/serve.py . "${process.env.GBDRAW_TEST_PYODIDE_SNAPSHOT_DIR}" ${serverPort}`
      : `python3 -m http.server ${serverPort} --bind 127.0.0.1`,
    env: process.env.GBDRAW_TEST_PYODIDE_SNAPSHOT_DIR ? { PROTO_SNAPSHOT_URL: '/__snapshot__/' } : {},
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
