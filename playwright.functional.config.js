// @ts-check
const { defineConfig } = require('@playwright/test');
const baseConfig = require('./playwright.config.js');

module.exports = defineConfig({
  ...baseConfig,
  testIgnore: [
    /.*performance\.playwright\.spec\.js/,
    /losat-cache-migration\.playwright\.spec\.js/
  ]
});
