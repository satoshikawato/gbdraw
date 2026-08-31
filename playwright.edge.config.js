// @ts-check
const { defineConfig, devices } = require('@playwright/test');
const performanceConfig = require('./playwright.perf.config.js');

module.exports = defineConfig({
  ...performanceConfig,
  projects: [
    {
      name: 'msedge',
      use: {
        ...devices['Desktop Chrome'],
        ...performanceConfig.use,
        channel: 'msedge'
      }
    }
  ]
});
