// S00 baseline runner. Serve one snapshot per origin first, for example
//   (cd <snapshot> && python3 -m http.server 4302 --bind 127.0.0.1)
// then run from the worktree root:
//   S00_TARGET=dev S00_BASE_URL=http://127.0.0.1:4302 S00_OUT_DIR=<dir> \
//   S00_VNIG_SESSION=<session> npx playwright test -c docs/internal/gui-remediation-20260928/results/s00/playwright.s00.config.js
'use strict';
const { defineConfig, devices } = require('@playwright/test');

module.exports = defineConfig({
  testDir: __dirname,
  testMatch: /s00-baseline\.spec\.js/,
  fullyParallel: false,
  workers: 1,
  retries: 0,
  reporter: 'list',
  use: {
    ...devices['Desktop Chrome'],
    viewport: { width: 1440, height: 900 },
    baseURL: process.env.S00_BASE_URL,
    trace: 'off',
    // S00_DIAG=1 (diagnosis only, never a budget pass): precise heap sizes per frame.
    launchOptions: process.env.S00_DIAG === '1' ? { args: ['--enable-precise-memory-info'] } : {}
  }
});
