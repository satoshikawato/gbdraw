// The input summary (gbdraw/api/input_summary.py): a record with no feature to draw still draws, and the status
// area shows the warning the command line prints; Run Info lists the summary lines. Python owns every word.
const { readFileSync } = require('node:fs');
const { join } = require('node:path');
const { test, expect } = require('@playwright/test');
const { generateAndWaitForResult, openApp } = require('./helpers/app-lifecycle.cjs');

const ZERO_FEATURES = () => readFileSync(join(process.cwd(), 'tests/fixtures/input_variants/zero_features_record.gb'));

test('a record with no feature to draw: the status area says so and Run Info lists the summary', async ({ page }) => {
  test.setTimeout(300000);
  await openApp(page);
  await page.getByLabel('GenBank/DDBJ File', { exact: true })
    .setInputFiles({ name: 'zero_features_record.gb', mimeType: 'text/plain', buffer: ZERO_FEATURES() });
  await page.waitForFunction(() => window.__GBDRAW_APP__.circularRecordList.length === 1);
  await generateAndWaitForResult(page);

  const notice = page.getByTestId('input-summary-notice');
  await expect(notice).toBeVisible();
  await expect(notice.locator('p').first()).toHaveText('Part of the input was not drawn.');
  await expect(notice.locator('[data-input-warning="no_drawable_features"]')).toHaveText(
    'Record EMPTY001.1 (zero_features_record.gb) has no feature to draw: it has no features.'
  );

  await page.evaluate(() => { window.__GBDRAW_APP__.resultPanelTab = 'run-info'; });
  await expect(page.locator('[data-run-info-input-summary] li')).toHaveText([
    'Record EMPTY001.1 (zero_features_record.gb), record 1 of 1, 1,000 bp: drew 0 features.'
  ]);
});
