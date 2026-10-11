// @ts-check
const { defineConfig } = require('@playwright/test');
const baseConfig = require('../../../playwright.config.js');

module.exports = defineConfig({
  ...baseConfig,
  // Run two browser workers per functional shard.
  workers: process.env.CI ? 2 : baseConfig.workers,
  reporter: process.env.CI
    ? [['list'], ['github']]
    : baseConfig.reporter,
  testIgnore: [
    /.*performance\.playwright\.spec\.js/,
    /losat-cache-migration\.playwright\.spec\.js/
  ]
});
