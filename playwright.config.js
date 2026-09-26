// @ts-check
const { defineConfig, devices } = require('@playwright/test');
const serverPort = Number(process.env.GBDRAW_WEB_TEST_PORT || 4173);
const serverUrl = `http://127.0.0.1:${serverPort}`;

module.exports = defineConfig({
  testDir: './tests/web',
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
    command: `python3 -m http.server ${serverPort} --bind 127.0.0.1`,
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
