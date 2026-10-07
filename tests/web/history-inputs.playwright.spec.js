const { test, expect } = require('@playwright/test');
const { readFileSync, writeFileSync } = require('node:fs');
const { resolve } = require('node:path');
const { gunzipSync } = require('node:zlib');
const { openApp, reveal } = require('./helpers/app-lifecycle.cjs');

const seed = resolve('gbdraw/web/gallery/sessions/HmmtDNA_basic_circular.gbdraw-session.json');
const conservationSession = resolve('tests/fixtures/sessions/synthetic_conservation.gbdraw-session.json.gz');
const sessionInput = 'input[type="file"][accept^=".json,"]';

const preparePage = async (page) => {
  await page.addInitScript(() => {
    window.__HISTORY_INPUT_EVENTS__ = [];
    window.__GBDRAW_TEST_HOOKS__ = {
      onSessionLifecycleEvent: (event) => window.__HISTORY_INPUT_EVENTS__.push(event)
    };
  });
  page.on('dialog', (dialog) => dialog.accept());
  await openApp(page);
};

const loadSession = async (page, path) => {
  await page.locator(sessionInput).setInputFiles(path);
  await page.waitForFunction(() => (
    window.__HISTORY_INPUT_EVENTS__.some((event) => event.name === 'interactiveReady')
    && !window.__GBDRAW_APP__.sessionImportPending
  ));
};

const generate = async (page) => {
  const before = await page.evaluate(() => window.__HISTORY_INPUT_EVENTS__
    .filter((event) => event.name === 'generate.processing-cleared').length);
  await page.getByRole('button', { name: 'Generate Diagram', exact: true }).click();
  await page.waitForFunction((count) => window.__HISTORY_INPUT_EVENTS__
    .filter((event) => event.name === 'generate.processing-cleared').length > count, before,
  { timeout: 180_000 });
  expect(await page.evaluate(() => window.__GBDRAW_APP__.errorLog?.summary || '')).toBe('');
  return page.evaluate(() => window.__GBDRAW_APP__.results[0].content);
};

const expectHistory = async (page, undo, redo, label, options = {}) => {
  await expect.poll(() => page.evaluate(() => {
    const history = window.__GBDRAW_HISTORY__;
    return [history.getUndoCount(), history.getRedoCount(), history.undoLabel()];
  }), options).toEqual([undo, redo, label]);
};

for (const inputMethod of ['keyboard', 'pointer']) {
  test(`Label Mode ${inputMethod} edit has one Undo step and survives generation and fresh Load`, async ({
    page, context, browser, baseURL
  }, testInfo) => {
    test.setTimeout(300_000);
    const externalRequests = [];
    const allowLocal = (route) => {
      if (new URL(route.request().url()).origin === new URL(baseURL).origin) {
        return route.continue();
      }
      externalRequests.push(route.request().url());
      return route.abort();
    };
    await context.route('**/*', allowLocal);
    await preparePage(page);
    await loadSession(page, seed);
    const originalSvg = await generate(page);
    const baseline = await page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount());
    await expectHistory(page, baseline, 0, 'Generate diagram');

    const summary = page.locator('summary[aria-label="Labels"]');
    await summary.press('Enter');
    const select = page.locator('#circular-label-mode');
    await expect(select).toHaveValue('out');
    if (inputMethod === 'keyboard') {
      await summary.focus();
      // The Labels and Label Mode help tips are tab stops before the select.
      await page.keyboard.press('Tab');
      await expect(summary.getByRole('button', { name: 'Help', exact: true })).toBeFocused();
      await page.keyboard.press('Tab');
      await expect(page.locator('button[aria-describedby="help-label-label-mode"]')).toBeFocused();
      await page.keyboard.press('Tab');
      await expect(select).toBeFocused();
      await page.keyboard.press('ArrowDown');
    } else {
      // Open the native popup through its pointer admission path.
      await select.click();
      await page.keyboard.press('ArrowDown');
      await page.keyboard.press('Enter');
    }
    await page.keyboard.press('Tab');
    await expect(select).toHaveValue('both');
    await testInfo.attach('after-select', {
      body: JSON.stringify(await page.evaluate(() => ({
        labels: window.__GBDRAW_APP__.form.labels_mode,
        undo: window.__GBDRAW_HISTORY__.getUndoCount(),
        undoLabel: window.__GBDRAW_HISTORY__.undoLabel()
      }))),
      contentType: 'application/json'
    });
    await expectHistory(page, baseline + 1, 0, 'Change setting');
    expect(await page.evaluate(() => window.__GBDRAW_APP__.results[0].content)).toBe(originalSvg);

    await page.getByRole('button', { name: 'Undo', exact: true }).click();
    await expect(select).toHaveValue('out');
    await expectHistory(page, baseline, 1, 'Generate diagram');
    await page.getByRole('button', { name: 'Redo', exact: true }).click();
    await expect(select).toHaveValue('both');
    await expectHistory(page, baseline + 1, 0, 'Change setting');
    const bothSvg = await generate(page);
    expect(bothSvg).not.toBe(originalSvg);
    await expectHistory(page, baseline + 2, 0, 'Generate diagram');

    // Undo the artifact replacement, then the form edit; neither may consume both.
    await page.getByRole('button', { name: 'Undo', exact: true }).click();
    await expectHistory(page, baseline + 1, 1, 'Change setting');
    await expect(select).toHaveValue('both');
    await page.getByRole('button', { name: 'Undo', exact: true }).click();
    await expect(select).toHaveValue('out');
    await expectHistory(page, baseline, 2, 'Generate diagram');

    const downloadPromise = page.waitForEvent('download');
    await page.getByRole('button', { name: 'Save Session', exact: true }).click();
    const download = await downloadPromise;
    const savedPath = testInfo.outputPath('restored.gbdraw-session.json.gz');
    await download.saveAs(savedPath);
    const saved = JSON.parse(gunzipSync(readFileSync(savedPath)));
    expect(saved.config.form.labels_mode).toBe('out');
    expect(await generate(page)).toBe(originalSvg);

    const freshContext = await browser.newContext({ baseURL });
    try {
      await freshContext.route('**/*', allowLocal);
      const freshPage = await freshContext.newPage();
      await preparePage(freshPage);
      await loadSession(freshPage, savedPath);
      expect(await freshPage.evaluate(() => window.__GBDRAW_APP__.form.labels_mode)).toBe('out');
      expect(await generate(freshPage)).toBe(originalSvg);
    } finally {
      await freshContext.close();
    }
    expect(externalRequests).toEqual([]);
  });
}

test('Linear File removal choices are atomic, undoable, and preserve one slot', async ({ page }) => {
  test.setTimeout(120_000);
  await page.setViewportSize({ width: 390, height: 844 });
  await preparePage(page);
  await page.getByRole('button', { name: 'Linear', exact: true }).click();
  await page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    const first = new File(['LOCUS A\n'], 'first.gbk', { type: 'text/plain', lastModified: 1 });
    const last = new File(['LOCUS B\n'], 'last.gbk', { type: 'text/plain', lastModified: 2 });
    app.linearSeqs[0].gb = first;
    app.linearSeqs[0].region_record_id = 'A1';
    app.linearSeqs[0].definition = 'First one';
    app.linearSeqs[0].depth = [new File(['1\t10\n'], 'first.tsv')];
    app.addLinearSeq();
    app.linearSeqs[1].gb = first;
    app.linearSeqs[1].region_record_id = 'A2';
    app.linearSeqs[1].definition = 'First two';
    app.addLinearSeq();
    app.linearSeqs[2].gb = last;
    app.linearSeqs[2].region_record_id = 'B1';
  });

  const sources = page.locator('[data-linear-source-card]');
  const dialog = page.getByRole('dialog', { name: 'Clear or delete File?' });
  await expect(sources).toHaveCount(2);
  const baseline = await page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount());
  const firstRemove = sources.first().getByRole('button', { name: /Remove$/ });

  await firstRemove.click();
  expect(await page.evaluate(() => ({
    dialog: { ...window.__GBDRAW_APP__.linearSourceRemovalDialog },
    target: window.__GBDRAW_APP__.linearSourceRemovalTarget?.uid || null
  }))).toEqual({
    dialog: { open: true, sourceUid: expect.any(String), origin: 'card' },
    target: expect.any(String)
  });
  await expect(dialog).toBeVisible();
  await expect(dialog).toContainText('first.gbk');
  await expect(dialog).toContainText('2 records');
  const clearFile = dialog.getByRole('button', { name: 'Clear file only', exact: true });
  const cancel = dialog.getByRole('button', { name: 'Cancel', exact: true });
  await expect(clearFile).toBeFocused();
  await page.keyboard.press('Shift+Tab');
  await expect(cancel).toBeFocused();
  await page.keyboard.press('Tab');
  await expect(clearFile).toBeFocused();
  await cancel.click();
  await expect(dialog).toHaveCount(0);
  await expect(firstRemove).toBeFocused();
  expect(await page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount())).toBe(baseline);

  await firstRemove.focus();
  await firstRemove.press('Enter');
  await expect(dialog).toBeVisible();
  await page.keyboard.press('Escape');
  await expect(dialog).toHaveCount(0);
  expect(await page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount())).toBe(baseline);

  await firstRemove.click();
  await dialog.getByRole('button', { name: 'Clear file only', exact: true }).click();
  await expect(sources).toHaveCount(2);
  await expect(sources.first().getByRole('button', { name: 'Choose GenBank / DDBJ File', exact: true })).toBeFocused();
  expect(await page.evaluate(() => window.__GBDRAW_APP__.linearSeqs.map((sequence) => ({
    uid: sequence.uid,
    file: sequence.gb?.name || null,
    selector: sequence.region_record_id,
    definition: sequence.definition,
    depth: (sequence.depth || []).map((file) => file?.name || null)
  })))).toEqual([
    expect.objectContaining({ file: null, selector: '', definition: '', depth: [null] }),
    expect.objectContaining({ file: 'last.gbk', selector: 'B1' })
  ]);
  await expectHistory(page, baseline + 1, 0, 'Clear Linear File');
  await page.getByRole('button', { name: 'Undo', exact: true }).click();
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.linearSeqs.map(
    (sequence) => sequence.region_record_id
  ))).toEqual(['A1', 'A2', 'B1']);
  await expect.poll(() => page.evaluate(() => [
    window.__GBDRAW_HISTORY__.getUndoCount(),
    window.__GBDRAW_HISTORY__.getRedoCount()
  ])).toEqual([baseline, 1]);
  await page.getByRole('button', { name: 'Redo', exact: true }).click();
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.linearSeqs.map(
    (sequence) => sequence.region_record_id
  ))).toEqual(['', 'B1']);
  await expectHistory(page, baseline + 1, 0, 'Clear Linear File');
  await page.getByRole('button', { name: 'Undo', exact: true }).click();

  const globalRemove = page.getByRole('button', { name: 'Remove last sequence', exact: true });
  await globalRemove.click();
  const globalDialog = page.getByRole('dialog', { name: 'Remove last File?' });
  await expect(globalDialog).toBeVisible();
  await expect(globalDialog).toContainText('last.gbk');
  await globalDialog.getByRole('button', { name: 'Cancel', exact: true }).click();
  expect(await page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount())).toBe(baseline);
  await globalRemove.click();
  await globalDialog.getByRole('button', { name: 'Delete card', exact: true }).click();
  await expect(sources).toHaveCount(1);
  await expect(sources.first().getByRole('button', { name: 'Choose GenBank / DDBJ File', exact: true })).toBeFocused();
  await expectHistory(page, baseline + 1, 0, 'Delete Linear File');
  await page.getByRole('button', { name: 'Undo', exact: true }).click();
  await expect(sources).toHaveCount(2);

  await page.getByRole('group', { name: 'Linear File actions' })
    .getByRole('button', { name: 'Add sequence', exact: true }).click();
  await expect(sources).toHaveCount(3);
  const beforeBlankRemoval = await page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount());
  await globalRemove.click();
  await expect(globalDialog).toHaveCount(0);
  await expect(sources).toHaveCount(2);
  await expectHistory(page, beforeBlankRemoval + 1, 0, 'Delete Linear File');
  expect(await page.locator('[data-linear-record-list]').evaluate(
    (element) => element.scrollWidth <= element.clientWidth + 1
  )).toBe(true);
});

// R11 input matrix (SE-02, SE-03, SE-04): control type x input means x focus
// state. Every changed control adds exactly one step, and Undo is LIFO.
const controlState = (page) => page.evaluate(() => {
  const app = window.__GBDRAW_APP__;
  return {
    rich: app.richFeaturePopup,
    input: app.cInputType,
    prefix: app.form.prefix,
    mode: app.mode,
    labels: app.form.labels_mode
  };
});
const historyCounts = (page) => page.evaluate(() => [
  window.__GBDRAW_HISTORY__.getUndoCount(),
  window.__GBDRAW_HISTORY__.getRedoCount()
]);

test('each discrete control change records one step for every input means and focus state', async ({ page }) => {
  test.setTimeout(300_000);
  const errors = [];
  page.on('pageerror', (error) => errors.push(error.message));
  await preparePage(page);
  await loadSession(page, seed);
  const [baseline] = await historyCounts(page);
  const states = [await controlState(page)];
  // Each change is one expected step, applied in order.
  const expectSteps = async (name, ...changes) => {
    changes.forEach((changed) => states.push({ ...states.at(-1), ...changed }));
    await expect.poll(() => controlState(page), { message: name }).toEqual(states.at(-1));
    await expect.poll(() => historyCounts(page), { message: name })
      .toEqual([baseline + states.length - 1, 0]);
  };

  const rich = page.locator('label.option-label', { hasText: 'Rich Feature Popup' }).first();
  await reveal(rich);
  const richBox = rich.locator('input[type="checkbox"]');
  const genbankLabel = page.locator('label', { hasText: 'GenBank' }).first();
  const prefix = page.locator('#output-prefix');
  const initial = states[0];

  // No focused text field: label text, Space, radio Arrow, and radio label text.
  await rich.locator('span').first().click({ position: { x: 4, y: 4 } });
  await expectSteps('checkbox label text', { rich: !initial.rich });
  await richBox.focus();
  await page.keyboard.press('Space');
  await expectSteps('checkbox Space', { rich: initial.rich });
  await genbankLabel.locator('input[type="radio"]').focus();
  await page.keyboard.press('ArrowRight');
  await expectSteps('radio Arrow', { input: 'gff' });
  await genbankLabel.click({ position: { x: 40, y: 8 } });
  await expectSteps('radio label text', { input: 'gb' });

  // A focused text field, then a typed one, before a checkbox or a mode button.
  await prefix.click();
  await richBox.click();
  await expectSteps('checkbox after a focused text field', { rich: !initial.rich });
  await prefix.click();
  await prefix.press('End');
  await page.keyboard.type('X');
  await richBox.click();
  await expectSteps('checkbox after typed text', { prefix: `${initial.prefix}X` }, { rich: initial.rich });
  await prefix.click();
  await prefix.press('End');
  await page.keyboard.type('Y');
  await page.getByRole('button', { name: 'Linear', exact: true }).click();
  await expectSteps('mode button after typed text', { prefix: `${initial.prefix}XY` }, { mode: 'linear' });

  // Undo returns through every recorded state in reverse order.
  const undo = page.getByRole('button', { name: 'Undo', exact: true });
  for (let index = states.length - 2; index >= 0; index -= 1) {
    await undo.click();
    await expect.poll(() => controlState(page), { message: `Undo to state ${index}` }).toEqual(states[index]);
    await expect.poll(() => historyCounts(page))
      .toEqual([baseline + index, states.length - 1 - index]);
  }
  expect(errors).toEqual([]);
});

test('History shortcuts work on a focused select, and a text field keeps native Undo', async ({ page }) => {
  test.setTimeout(300_000);
  await preparePage(page);
  await loadSession(page, seed);
  const [baseline] = await historyCounts(page);
  await page.locator('summary[aria-label="Labels"]').click();
  const select = page.locator('#circular-label-mode');
  await expect(select).toHaveValue('out');
  await select.focus();
  await page.keyboard.press('ArrowDown');
  await expect(select).toHaveValue('both');
  await expect.poll(() => historyCounts(page)).toEqual([baseline + 1, 0]);

  await page.keyboard.press('Control+z');
  await expect(select).toHaveValue('out');
  await expect.poll(() => historyCounts(page)).toEqual([baseline, 1]);
  await page.keyboard.press('Control+Shift+z');
  await expect(select).toHaveValue('both');
  await expect.poll(() => historyCounts(page)).toEqual([baseline + 1, 0]);
  await page.keyboard.press('Control+z');
  await expect(select).toHaveValue('out');
  await page.keyboard.press('Control+y');
  await expect(select).toHaveValue('both');
  await page.keyboard.press('Control+z');
  await expect(select).toHaveValue('out');
  await expect(select).toBeFocused();
  // Leaving the select does not record the Undo as a new edit.
  await page.keyboard.press('Tab');
  await expect.poll(() => historyCounts(page)).toEqual([baseline, 1]);

  const prefix = page.locator('#output-prefix');
  const before = await prefix.inputValue();
  await prefix.click();
  await prefix.press('End');
  await page.keyboard.type('abc');
  await expect(prefix).toHaveValue(`${before}abc`);
  await page.keyboard.press('Control+z');
  await expect(prefix).not.toHaveValue(`${before}abc`);
  await expect(prefix).toBeFocused();
  expect(await historyCounts(page)).toEqual([baseline, 1]);
});

// B21 (R11): a file chosen in a hidden input that a button or a label opens
// records one step, whether the picker, the input, or a script sets the file.
// An Add Seq step commits after the ring reader answers (B22); the first read
// starts the diagram worker.
const ringReadTimeout = { timeout: 120_000 };
const ringRows = (page) => page.evaluate(() => ({
  files: (window.__GBDRAW_APP__.files.c_conservation_fastas || []).map((file) => file?.name || null),
  rows: window.__GBDRAW_APP__.circularConservation.series.map((entry) => entry.label)
}));
const fastaFile = (testInfo, name) => {
  const path = testInfo.outputPath(name);
  writeFileSync(path, `>${name}\nACGTACGTACGTACGTACGTACGTACGTACGT\n`);
  return path;
};

test('each Circular LOSAT ring file added through Add Seq is one Undo step', async ({ page }, testInfo) => {
  test.setTimeout(180_000);
  await preparePage(page);
  await loadSession(page, conservationSession);
  const [baseline] = await historyCounts(page);
  const initial = await ringRows(page);
  const addSeq = page.locator('button', { hasText: 'Add Seq' });
  await reveal(addSeq);
  const ringInput = page.locator('input[type="file"][accept^=".fa,.fas,.fasta,.fna,.ffn,.gb"]');
  const states = [initial];
  const expectAdded = async (name) => {
    states.push({ files: [...states.at(-1).files, name], rows: [...states.at(-1).rows, name.replace(/\.fasta$/, '')] });
    await expect.poll(() => ringRows(page), { message: name }).toEqual(states.at(-1));
    await expectHistory(page, baseline + states.length - 1, 0, 'Change uploaded file', ringReadTimeout);
  };

  const [chooser] = await Promise.all([page.waitForEvent('filechooser'), addSeq.click()]);
  await chooser.setFiles(fastaFile(testInfo, 'ring-e.fasta'));
  await expectAdded('ring-e.fasta');
  await ringInput.setInputFiles(fastaFile(testInfo, 'ring-f.fasta'));
  await expectAdded('ring-f.fasta');
  await ringInput.setInputFiles({
    name: 'ring-g.fasta', mimeType: 'text/plain', buffer: Buffer.from('>ring-g\nACGTACGTACGT\n')
  });
  await expectAdded('ring-g.fasta');

  const undo = page.getByRole('button', { name: 'Undo', exact: true });
  for (let index = states.length - 2; index >= 0; index -= 1) {
    await undo.click();
    await expect.poll(() => ringRows(page), { message: `Undo to state ${index}` }).toEqual(states[index]);
    await expect.poll(() => historyCounts(page)).toEqual([baseline + index, states.length - 1 - index]);
  }
  await page.getByRole('button', { name: 'Redo', exact: true }).click();
  await expect.poll(() => ringRows(page)).toEqual(states[1]);
  await expectHistory(page, baseline + 1, 2, 'Change uploaded file');
});

test('an uploaded BLAST row comparison sequence chosen through its label is one Undo step', async ({ page }, testInfo) => {
  test.setTimeout(120_000);
  await preparePage(page);
  const upload = page.locator('input[type="radio"][value="upload"]');
  await reveal(upload);
  await upload.check();
  const blast = testInfo.outputPath('ring-hits.tsv');
  writeFileSync(blast, 'q1\ts1\t99.0\t100\t1\t0\t1\t100\t1\t100\t1e-50\t200\n');
  await page.locator('input[type="file"][aria-label="BLAST outfmt 6/7 files"]').setInputFiles(blast);
  const sequenceSources = () => page.evaluate(() => (
    window.__GBDRAW_APP__.files.c_conservation_sequence_sources || []
  ).map((file) => file?.name || null));
  const companion = page.locator('label', { hasText: 'Comparison sequence (optional)' })
    .locator('input[type="file"]');
  await expect(companion).toHaveCount(1);
  await expect.poll(sequenceSources).toEqual([]);
  const [baseline] = await historyCounts(page);

  await companion.setInputFiles(fastaFile(testInfo, 'ring-subject.fasta'));
  await expect.poll(sequenceSources).toEqual(['ring-subject.fasta']);
  await expectHistory(page, baseline + 1, 0, 'Change uploaded file');
  await page.getByRole('button', { name: 'Undo', exact: true }).click();
  await expect.poll(sequenceSources).toEqual([]);
  await expect.poll(() => historyCounts(page)).toEqual([baseline, 1]);
  await page.getByRole('button', { name: 'Redo', exact: true }).click();
  await expect.poll(sequenceSources).toEqual(['ring-subject.fasta']);
  await expectHistory(page, baseline + 1, 0, 'Change uploaded file');
});

// B22 (R11, D12): the Python ring reader names a GenBank or DDBJ ring row after
// the file is chosen. That label belongs to the row's Add Seq step, so Redo
// restores it and the add is still one step.
const ringFlatFile = ({ definition }) => {
  const sequence = 'ACGTTGCAAC'.repeat(12);
  const origin = [];
  for (let start = 0; start < sequence.length; start += 60) {
    const chunk = sequence.slice(start, start + 60).match(/.{1,10}/g).join(' ');
    origin.push(`${String(start + 1).padStart(9)} ${chunk}`);
  }
  return [
    `LOCUS       RING1            ${String(sequence.length).padStart(11)} bp    DNA     linear   UNK 01-JAN-1980`,
    `DEFINITION  ${definition}`,
    'ACCESSION   RING1',
    'VERSION     RING1',
    'KEYWORDS    .',
    'SOURCE      Synthetic ring organism',
    '  ORGANISM  Synthetic ring organism',
    '            Unclassified.',
    'FEATURES             Location/Qualifiers',
    'ORIGIN',
    ...origin,
    '//',
    ''
  ].join('\n');
};

test('a GenBank or DDBJ ring added through Add Seq keeps its record label through Undo and Redo', async ({ page }, testInfo) => {
  test.setTimeout(180_000);
  await preparePage(page);
  await loadSession(page, conservationSession);
  const [baseline] = await historyCounts(page);
  const initial = await ringRows(page);
  await reveal(page.locator('button', { hasText: 'Add Seq' }));
  const ringInput = page.locator('input[type="file"][accept^=".fa,.fas,.fasta,.fna,.ffn,.gb"]');
  const added = (state, name, label) => ({ files: [...state.files, name], rows: [...state.rows, label] });
  const expectState = async (state, undo, redo, message) => {
    await expect.poll(() => ringRows(page), { message, ...ringReadTimeout }).toEqual(state);
    await expect.poll(() => historyCounts(page), { message, ...ringReadTimeout }).toEqual([baseline + undo, redo]);
  };

  // GenBank through the input: the DEFINITION label is part of the one step.
  const genbankPath = testInfo.outputPath('ring-genbank.gbk');
  writeFileSync(genbankPath, ringFlatFile({ definition: 'synthetic ring genbank.' }));
  const genbank = added(initial, 'ring-genbank.gbk', 'synthetic ring genbank');
  await ringInput.setInputFiles(genbankPath);
  await expectState(genbank, 1, 0, 'GenBank added');
  await expectHistory(page, baseline + 1, 0, 'Change uploaded file');
  await page.getByRole('button', { name: 'Undo', exact: true }).click();
  await expectState(initial, 0, 1, 'GenBank add undone');
  await page.getByRole('button', { name: 'Redo', exact: true }).click();
  await expectState(genbank, 1, 0, 'GenBank add redone with its DEFINITION');

  // DDBJ set by a script (buffer): an empty DEFINITION falls back to the organism.
  const ddbj = added(genbank, 'ring-ddbj.ddbj', 'Synthetic ring organism');
  await ringInput.setInputFiles({
    name: 'ring-ddbj.ddbj', mimeType: 'text/plain', buffer: Buffer.from(ringFlatFile({ definition: '.' }))
  });
  await expectState(ddbj, 2, 0, 'DDBJ added');
  await expectHistory(page, baseline + 2, 0, 'Change uploaded file');

  // A label typed before the reader answers wins, in the same step.
  const typed = added(ddbj, 'ring-typed.gbk', 'Typed ring');
  await page.evaluate(async (text) => {
    const { state } = await import('/gbdraw/web/js/state.js');
    const file = new File([text], 'ring-typed.gbk', { type: 'text/plain', lastModified: 0 });
    window.__GBDRAW_APP__.addCircularConservationComparisonFile({ target: { files: [file], value: '' } });
    state.activeDrawing().circularConservation.series[state.activeDrawing().circularConservation.series.length - 1].label = 'Typed ring';
  }, ringFlatFile({ definition: 'synthetic ring typed.' }));
  await expectState(typed, 3, 0, 'typed ring added');
  await expectHistory(page, baseline + 3, 0, 'Change uploaded file');

  const states = [initial, genbank, ddbj, typed];
  for (let index = states.length - 2; index >= 0; index -= 1) {
    await page.getByRole('button', { name: 'Undo', exact: true }).click();
    await expectState(states[index], index, states.length - 1 - index, `Undo to state ${index}`);
  }
  for (let index = 1; index < states.length; index += 1) {
    await page.getByRole('button', { name: 'Redo', exact: true }).click();
    await expectState(states[index], index, states.length - 1 - index, `Redo to state ${index}`);
  }
});

const replayState = (page) => page.evaluate(() => {
  const app = window.__GBDRAW_APP__;
  const names = (files) => (files || []).map((file) => file?.name || null);
  return {
    source: app.files.c_conservation_blasts_source,
    blasts: names(app.files.c_conservation_blasts),
    fastas: names(app.files.c_conservation_fastas),
    sequenceSources: names(app.files.c_conservation_sequence_sources),
    rows: app.circularConservation.series.map((entry) => [entry.label, entry.fileName]),
    slots: (app.adv.circular_track_slots || [])
      .filter((slot) => slot?.renderer === 'sequence_conservation')
      .map((slot) => [slot.id, slot.params?.series_key, slot.params?.label])
  };
});

// B23: the custom stack's ring slots name the LOSAT-cache replay's BLAST rows,
// so Undo restores those rows with the ring rows and the next Generate matches.
test('Undo of a ring change in a LOSAT-cache replay restores its rows and Generate result', async ({ page }, testInfo) => {
  test.setTimeout(300_000);
  await preparePage(page);
  await loadSession(page, conservationSession);
  const customStack = page.locator('label', { hasText: 'Use custom stack' }).locator('input[type="checkbox"]');
  await reveal(customStack);
  await customStack.check();
  const baseline = await generate(page);
  const loaded = await replayState(page);
  expect(loaded.source).toBe('losat-cache');
  expect(loaded.blasts).toHaveLength(3);
  expect(loaded.slots).toHaveLength(3);
  const [steps] = await historyCounts(page);
  const undo = page.getByRole('button', { name: 'Undo', exact: true });
  const removeSeries = page.locator('button[title="Remove series"]');
  await reveal(removeSeries.first());
  for (const index of [0, 2]) {
    await removeSeries.nth(index).click();
    await expect.poll(async () => (await replayState(page)).rows.length).toBe(2);
    await expectHistory(page, steps + 1, 0, 'Change setting');
    await undo.click();
    await expect.poll(() => replayState(page), { message: `Undo Remove series ${index}` }).toEqual(loaded);
  }
  await page.locator('input[type="file"][accept^=".fa,.fas,.fasta,.fna,.ffn,.gb"]')
    .setInputFiles(fastaFile(testInfo, 'ring-e.fasta'));
  await expectHistory(page, steps + 1, 0, 'Change uploaded file', ringReadTimeout);
  await undo.click();
  await expect.poll(() => replayState(page), { message: 'Undo Add Seq' }).toEqual(loaded);
  expect(await generate(page)).toBe(baseline);
});
