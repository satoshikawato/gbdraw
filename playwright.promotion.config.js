// @ts-check
const { defineConfig } = require('@playwright/test');
const baseConfig = require('./playwright.config.js');

// Promotion-only checks (docs/internal/WEB_PERIODIC_AUDIT.md). No PR or dev
// workflow runs this configuration, and the shared configurations match
// `.playwright.spec.js`, which `.promotion.spec.js` is not.
module.exports = defineConfig({
  ...baseConfig,
  testMatch: /.*\.promotion\.spec\.js/,
  testIgnore: [],
  fullyParallel: false,
  retries: 0,
  workers: 1,
  reporter: 'list',
  timeout: 3_600_000,
  use: {
    ...baseConfig.use,
    trace: 'off'
  }
});
