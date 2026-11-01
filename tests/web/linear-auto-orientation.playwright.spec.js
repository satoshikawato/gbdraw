// "Orient records from comparisons" (Linear) asks Python, through the diagram
// Worker, what `--auto_orient` decides from the comparisons of the last
// Generate, and writes the decided Reverse complement checkboxes as one History
// step; the next Generate draws them. tests/test_web_auto_orientation.py shows
// that the helper and the CLI planner decide with the same owner.
const { test, expect } = require('@playwright/test');
const { readFileSync } = require('node:fs');
const path = require('node:path');
const { evaluateWithRetainedPromise, generateAndWaitForResult } = require('./helpers/app-lifecycle.cjs');
const { openFresh, settle } = require('./helpers/audit-browser.cjs');

const root = path.resolve(process.env.GBDRAW_REPO || process.cwd());
const examples = path.join(root, 'examples');

// The tblastx table with MelaMJNV's coordinates mirrored: its hits against MjeNMV are mostly reverse.
const mirroredTable = () => {
  const subject = readFileSync(path.join(examples, 'MelaMJNV.gb'), 'utf8');
  const length = Number(/^LOCUS\s+\S+\s+(\d+)\s+bp/m.exec(subject)[1]);
  const rows = readFileSync(path.join(examples, 'MjeNMV.MelaMJNV.tblastx.out'), 'utf8').split('\n')
    .filter((line) => line && !line.startsWith('#'))
    .map((line) => {
      const fields = line.split('\t');
      fields[8] = String(length + 1 - Number(fields[8]));
      fields[9] = String(length + 1 - Number(fields[9]));
      return fields.join('\t');
    });
  return { subject, table: `${rows.join('\n')}\n` };
};
// The table as tblastx wrote it: MelaMJNV's hits are mostly forward in its source frame.
const forwardTable = () => readFileSync(path.join(examples, 'MjeNMV.MelaMJNV.tblastx.out'), 'utf8');
const useTable = (page, blast, name) => page.evaluate(async ({ blast, name }) => {
  const app = window.__GBDRAW_APP__;
  app.linearComparisonPlan.edges.splice(0, app.linearComparisonPlan.edges.length, {
    id: 'auto-orient-edge', queryUid: app.linearSeqs[0].uid, subjectUid: app.linearSeqs[1].uid,
    included: true, fileActive: true, losatFilenameActive: false, source: 'upload',
    file: new File([blast], name, { type: 'text/plain' }), losatFilename: ''
  });
  await window.Vue.nextTick();
}, { blast, name });

const reverseFlags = (page) => page.evaluate(() => window.__GBDRAW_APP__.linearSeqs.map((seq) => seq.region_reverse));
// The orientation the last Generate drew: the records of its committed request.
const drawnReverseFlags = (page) => page.evaluate(async () => {
  const { getCommittedCanonicalSession } = await import('./js/services/config.js');
  return getCommittedCanonicalSession().renderRequest.records
    .map((record) => Boolean(record.presentation?.reverseComplement));
});

test('Orient records from comparisons checks the reversed record as one Undo step, and Generate draws it', async ({ page }) => {
  test.setTimeout(420_000);
  await openFresh(page);
  const { subject, table } = mirroredTable();
  await page.evaluate(async ({ first, second }) => {
    const app = window.__GBDRAW_APP__;
    app.selectDrawing('linear');
    await window.Vue.nextTick();
    if (app.linearSeqs.length < 2) app.addLinearSeq();
    app.setLinearSeqPrimaryFile(0, 'gb', new File([first], 'MjeNMV.gb', { type: 'text/plain' }));
    app.setLinearSeqPrimaryFile(1, 'gb', new File([second], 'MelaMJNV.gb', { type: 'text/plain' }));
    app.linearComparisonPlan.mode = 'adjacent';
    app.linearComparisonPlan.defaultSource = 'upload';
    app.form.show_labels_linear = 'none';
    await window.Vue.nextTick();
  }, { first: readFileSync(path.join(examples, 'MjeNMV.gb'), 'utf8'), second: subject });
  await useTable(page, table, 'MjeNMV.MelaMJNV.mirrored.tsv');
  await settle(page);

  const button = page.getByTestId('linear-auto-orientation-button');
  await button.evaluate((element) => {
    for (let details = element.closest('details'); details; details = details.parentElement?.closest('details')) {
      details.open = true;
    }
  });
  // Before a Generate there are no comparisons to read.
  await expect(button).toBeDisabled();
  await expect(button).toHaveAttribute('aria-describedby', 'help-linear-auto-orientation');

  await generateAndWaitForResult(page);
  expect(await drawnReverseFlags(page)).toEqual([false, false]);
  await expect(button).toBeEnabled();
  await button.click();
  const status = page.getByTestId('linear-auto-orientation-status');
  await expect(status).toContainText('Generate to draw the new orientation');
  expect(await reverseFlags(page)).toEqual([false, true]);
  await expect(page.getByLabel('Reverse complement for sequence 2', { exact: true })).toBeChecked();
  expect(await page.evaluate(() => window.__GBDRAW_HISTORY__.undoLabel())).toBe('Orient records from comparisons');

  await evaluateWithRetainedPromise(page, () => window.__GBDRAW_HISTORY__.undo());
  expect(await reverseFlags(page)).toEqual([false, false]);
  await evaluateWithRetainedPromise(page, () => window.__GBDRAW_HISTORY__.redo());
  expect(await reverseFlags(page)).toEqual([false, true]);

  // A second press before Generate reads the same comparisons against changed orientation.
  await button.click();
  await expect(status).toContainText('changed since the last Generate');
  expect(await reverseFlags(page)).toEqual([false, true]);

  await generateAndWaitForResult(page);
  expect(await drawnReverseFlags(page)).toEqual([false, true]);
  // Drawn as reversed, the comparisons give every record its current orientation.
  await button.click();
  await expect(status).toContainText('give every record its current orientation');
  expect(await reverseFlags(page)).toEqual([false, true]);

  // With hits mostly forward in its source frame, a press unchecks the record it reversed.
  await useTable(page, forwardTable(), 'MjeNMV.MelaMJNV.tblastx.out');
  await generateAndWaitForResult(page);
  await button.click();
  await expect(status).toContainText('Set Reverse complement from the comparisons');
  expect(await reverseFlags(page)).toEqual([false, false]);
});
