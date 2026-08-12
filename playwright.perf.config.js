// @ts-check
const { defineConfig } = require('@playwright/test');
const baseConfig = require('./playwright.config.js');

module.exports = defineConfig({
  ...baseConfig,
  testMatch: /.*performance\.playwright\.spec\.js/,
  testIgnore: [],
  fullyParallel: false,
  retries: 0,
  workers: 1,
  use: {
    ...baseConfig.use,
    trace: 'retain-on-failure'
  }
});
