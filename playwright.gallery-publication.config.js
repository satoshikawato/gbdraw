// @ts-check
const { defineConfig } = require('@playwright/test');
const baseConfig = require('./playwright.config.js');

module.exports = defineConfig({
  ...baseConfig,
  testDir: './tests/web/contracts',
  testMatch: 'gallery-publication-parity.serial.spec.js',
  fullyParallel: false,
  retries: 0,
  workers: 1,
  timeout: 600_000,
  use: {
    ...baseConfig.use,
    trace: 'retain-on-failure'
  }
});
