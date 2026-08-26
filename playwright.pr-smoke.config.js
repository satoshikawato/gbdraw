// @ts-check
const { defineConfig } = require('@playwright/test');
const baseConfig = require('./playwright.config.js');

module.exports = defineConfig({
  ...baseConfig,
  testIgnore: [
    /.*performance\.playwright\.spec\.js/,
    /losat-cache-migration\.playwright\.spec\.js/
  ],
  grep: /@pr-smoke/,
  retries: 0,
  workers: 1,
  reporter: 'list',
  use: {
    ...baseConfig.use,
    trace: 'retain-on-failure'
  },
  projects: baseConfig.projects.filter(({ name }) => name === 'chromium')
});
