// @ts-check
// Playwright config for the manual audit sweeps in tools/audit/.
// These sweeps are not CI checks. Run them with:
//   GBDRAW_WEB_TEST_PORT=<port> npx playwright test -c tools/audit/playwright.audit.config.cjs <spec>
const { resolve } = require('node:path');
const { defineConfig, devices } = require('@playwright/test');
const { AUDIT_OUT } = require('./helpers/audit-common.cjs');

const repositoryRoot = resolve(__dirname, '..', '..');
const serverPort = Number(process.env.GBDRAW_WEB_TEST_PORT || 4173);
const serverUrl = `http://127.0.0.1:${serverPort}`;

module.exports = defineConfig({
  testDir: __dirname,
  testMatch: /.*\.audit\.spec\.js/,
  fullyParallel: false,
  workers: 1,
  retries: 0,
  reporter: 'line',
  outputDir: resolve(AUDIT_OUT, 'playwright'),
  use: {
    baseURL: serverUrl,
    acceptDownloads: true
  },
  webServer: {
    command: `python3 -m http.server ${serverPort} --bind 127.0.0.1`,
    cwd: repositoryRoot,
    url: serverUrl,
    reuseExistingServer: true,
    stderr: 'ignore'
  },
  projects: [
    {
      name: 'chromium',
      use: { ...devices['Desktop Chrome'] }
    }
  ]
});
