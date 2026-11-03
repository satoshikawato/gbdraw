// OV-1949: a per-feature edit names its feature by the record key of the view
// it was made in (Circular canvas `record-<k>`, one Circular record drawn alone
// `circular-<recordId>-#<k>`, a Linear card drawing one record `<uid>`, and
// `<uid>:<k>` for record k of the card's file while the card draws another
// record). Each Generate moves the rows of every record whose key changes, so
// a "This feature only" fill, a label edit, Feature visibility Off and a stroke
// made in one view keep applying in the others, live and after Generate, and
// after Undo/Redo of that Generate. Duplicates stay separate (PD-OI-069 rev 2),
// a Session written by 0.14.2 keeps its edits, and Run Info's Source recipe
// names the CLI key.
const { test, expect } = require('@playwright/test');
const { readFileSync } = require('node:fs');
const { spawnSync } = require('node:child_process');
const { join, resolve } = require('node:path');
const { openApp, openLinearFileCards } = require('./helpers/app-lifecycle.cjs');
const { generate } = require('./helpers/mode-transition.cjs');
const { showDrawing } = require('./helpers/drawing-bar.cjs');

const repoRoot = resolve(process.env.GBDRAW_REPO || process.cwd());
const twoRecords = ['BGC0000708', 'BGC0000709']
  .map((name) => readFileSync(join(repoRoot, `tests/test_inputs/${name}.gbk`), 'utf8')).join('\n');
const session0142 = join(repoRoot, 'tests/fixtures/sessions/record-identity-circular-single.v48.gbdraw-session.json.gz');
// Records 1 and 2 of the batch fixture; record 2 gets a second CDS at the
// coordinates of its first, so the two have one hash (IDs `<hash>~<n>`).
const batchRecords = readFileSync(join(repoRoot, 'tests/fixtures/web_batch_two_records.gb'), 'utf8')
  .split(/^\/\/\s*$/m).map((chunk) => chunk.trim()).filter((chunk) => chunk.startsWith('LOCUS')).map((chunk) => `${chunk}\n//\n`);
const identicalTwins = [batchRecords[0], batchRecords[1].replace(
  '     CDS             301..600\n',
  '     CDS             301..600\n                     /locus_tag="TWIN_A"\n     CDS             301..600\n'
)].join('');

const FILL = '#c026d3';
const STROKE = '#00aa00';
const LABEL = 'RI_LABEL_ONE';
const LABEL_LIVE = 'RI_LABEL_TWO';
// Record 2 (BGC0000709) of the file, by view.
const CANVAS_KEY = 'record-2';
const SINGLE_KEY = 'circular-BGC0000709-#2';

// No edit, rule preparation or label rerender is running: a label rerender
// remounts the Result, and the remount resets the label text a popup shows.
const settled = (page) => page.waitForFunction(() => {
  const app = window.__GBDRAW_APP__;
  return !window.__GBDRAW_HISTORY__.mutationPending() && !app.ruleMatchingPending && !app.labelReflowProcessing && !app.processing;
}, null, { timeout: 120_000 });

const historySettled = (page) => expect.poll(() => page.evaluate(() => {
  const history = window.__GBDRAW_HISTORY__;
  return !history.restoring.value && !history.traversalPending() && !history.mutationPending();
}), { timeout: 120_000 }).toBe(true);

const startCircular = async (page, content = twoRecords, form = {}) => {
  page.on('dialog', (dialog) => dialog.accept());
  await openApp(page);
  await showDrawing(page, { kind: 'circular' });
  await page.evaluate(({ content, form }) => {
    const app = window.__GBDRAW_APP__;
    app.setCircularInputType('gb');
    // The File transition runs the record discovery the record list waits for (R10).
    app.setCircularSourceFile('c_gb', new File([content], 'two-records.gbk', { type: 'text/plain', lastModified: 1 }));
    Object.assign(app.form, { multi_record_canvas: true, suppress_gc: true, suppress_skew: true, labels_mode: 'out' }, form);
    app.sessionTitle = 'record-identity-moves';
  }, { content, form });
  await page.waitForFunction(() => window.__GBDRAW_APP__?.circularRecordList?.length === 2, null, { timeout: 180_000 });
};

// The two-record file chosen in `files` Linear cards: each choice expands into
// one card per record (record 1 keeps the card, record 2 gets a card of its
// own). Returns the card uids, in card order.
const startLinear = async (page, files = 1) => {
  page.on('dialog', (dialog) => dialog.accept());
  await openApp(page);
  await showDrawing(page, { kind: 'linear' });
  await page.evaluate(async () => {
    const app = window.__GBDRAW_APP__;
    app.lInputType = 'gb';
    await app.setLinearComparisonGlobalAction('none');
    Object.assign(app.form, { show_gc: false, show_skew: false, show_depth: false, show_labels_linear: 'all' });
    app.sessionTitle = 'record-identity-moves-linear';
  });
  for (let file = 0; file < files; file += 1) {
    await page.evaluate(({ content, file: number }) => {
      const app = window.__GBDRAW_APP__;
      if (number > 0) app.addLinearSeq();
      app.setLinearSeqPrimaryFile(app.linearSeqs.length - 1, 'gb',
        new File([content], `two-records-${number + 1}.gbk`, { type: 'text/plain', lastModified: 1 }));
    }, { content: twoRecords, file });
    await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.linearSeqs.map((seq) => seq.region_record_id)),
      { timeout: 180_000 }).toEqual(Array(file + 1).fill(['BGC0000708', 'BGC0000709']).flat());
  }
  return page.evaluate(() => window.__GBDRAW_APP__.linearSeqs.map((seq) => seq.uid));
};

// The card's Record selector, found by the card's uid, in its opened File
// card, record list and Record options.
const selectLinearRecord = async (page, uid, recordId) => {
  await openLinearFileCards(page);
  await page.evaluate((key) => {
    document.querySelectorAll(`[data-linear-source-records], [data-linear-record-options="${CSS.escape(key)}"]`)
      .forEach((details) => { details.open = true; });
  }, uid);
  await page.locator(`select[data-record-selector-uid="${uid}"]`).selectOption(recordId);
  await expect.poll(() => page.evaluate((key) => window.__GBDRAW_APP__.linearSeqs
    .find((seq) => seq.uid === key).region_record_id, uid)).toBe(recordId);
};

// Four drawn CDS of the record keyed `recordKey`, one per edit, so hiding one
// leaves the others drawn.
const pickTargets = (page, recordKey) => page.evaluate(async (key) => {
  const app = window.__GBDRAW_APP__;
  const svg = app.svgContainer.querySelector('svg');
  const { getFeatureElements } = await import('/gbdraw/web/js/app/feature-editor/svg-actions.js');
  const ids = app.extractedFeatures
    .filter((feature) => feature.recordKey === key && feature.type === 'CDS'
      && getFeatureElements(svg, String(feature.svg_id)).length > 0)
    .map((feature) => feature.biologicalFeatureId);
  if (ids.length < 4) {
    throw new Error(`record ${key} has ${ids.length} drawn CDS; records: ${[...new Set(app.extractedFeatures.map((f) => f.recordKey))]}`);
  }
  return { fill: ids[0], label: ids[1], hide: ids[2], stroke: ids[3] };
}, recordKey);

const openPopup = async (page, recordKey, featureId) => {
  await page.evaluate(({ recordKey: key, featureId: id }) => {
    const app = window.__GBDRAW_APP__;
    const feature = app.extractedFeatures.find((candidate) => candidate.recordKey === key && candidate.biologicalFeatureId === id);
    if (!feature) throw new Error(`no feature ${id} in record ${key}`);
    app.openFeatureEditorFromList(feature, { clientX: 220, clientY: 220 });
  }, { recordKey, featureId });
  const popup = page.getByRole('dialog', { name: /Feature details:/ });
  await expect(popup).toBeVisible();
  await settled(page);
  return popup;
};

const closePopup = async (page) => {
  await page.getByRole('button', { name: 'Close feature popup', exact: true }).click();
  await expect(page.getByRole('dialog', { name: /Feature details:/ })).toHaveCount(0);
};

// The narrow scope: "This feature only" (fill, stroke) or "This feature"
// (Feature visibility).
const chooseScope = async (page, title, button = { name: 'This feature only', exact: true }) => {
  const scope = page.getByRole('dialog', { name: title });
  await expect(scope).toBeVisible();
  await settled(page);
  await scope.getByRole('button', button).click();
  await expect(scope).toHaveCount(0, { timeout: 120_000 });
  await settled(page);
};

const editFill = async (page, recordKey, featureId) => {
  const popup = await openPopup(page, recordKey, featureId);
  await popup.getByLabel('Feature fill color').first().fill(FILL);
  await chooseScope(page, 'Color Change Scope');
  await closePopup(page);
};

const editLabel = async (page, recordKey, featureId, text) => {
  const popup = await openPopup(page, recordKey, featureId);
  const input = popup.getByLabel('Label text', { exact: true });
  await input.fill(text);
  await popup.getByLabel('Label visibility', { exact: true }).selectOption('on');
  await expect(input).toHaveValue(text);
  await popup.getByRole('button', { name: 'Apply Label', exact: true }).click();
  await settled(page);
  await expect.poll(() => page.evaluate(({ key, id }) => window.__GBDRAW_APP__.featureOverrides[JSON.stringify([key, id])]?.labelText ?? null,
    { key: recordKey, id: featureId }), 'the label edit is written to the feature row').toBe(text);
  await closePopup(page);
};

const editHide = async (page, recordKey, featureId) => {
  const popup = await openPopup(page, recordKey, featureId);
  await popup.getByLabel('Feature visibility', { exact: true }).selectOption('off');
  await chooseScope(page, 'Feature Visibility Scope', { name: /^This feature One rule/ });
  await closePopup(page);
};

const editStroke = async (page, recordKey, featureId) => {
  const popup = await openPopup(page, recordKey, featureId);
  const mode = popup.getByLabel('Feature stroke color mode', { exact: true });
  if (await mode.inputValue() !== 'color') {
    await mode.selectOption('color');
    await settled(page);
  }
  await popup.getByLabel('Feature stroke color', { exact: true }).fill(STROKE);
  await chooseScope(page, 'Stroke Change Scope');
  await closePopup(page);
};

const editFour = async (page, recordKey) => {
  const targets = await pickTargets(page, recordKey);
  await editFill(page, recordKey, targets.fill);
  await editLabel(page, recordKey, targets.label, LABEL);
  await editHide(page, recordKey, targets.hide);
  await editStroke(page, recordKey, targets.stroke);
  return targets;
};

// The displayed Result, read for the four target features. With `recordKey`,
// each is found by the SVG ID its feature of that record has in this Result
// (Linear: `<feature id>_record_<n>`); without it, by its biological feature
// ID (Circular).
const shown = (page, targets, recordKey = null) => page.evaluate(async ({ ids: { fill, label, hide, stroke }, key }) => {
  const app = window.__GBDRAW_APP__;
  const svg = app.svgContainer.querySelector('svg');
  const { getFeatureElements, getFeatureFillElements } = await import('/gbdraw/web/js/app/feature-editor/svg-actions.js');
  const drawn = (element) => element.getAttribute('display') !== 'none' && element.closest('[display="none"]') === null;
  const svgId = (id) => (key === null ? id : String(app.extractedFeatures
    .find((feature) => feature.recordKey === key && feature.biologicalFeatureId === id)?.svg_id ?? ''));
  const elements = (id) => (svgId(id) ? getFeatureElements(svg, svgId(id)) : []);
  return {
    records: [...new Set([...svg.querySelectorAll('[data-gbdraw-record-id]')]
      .map((element) => element.getAttribute('data-gbdraw-record-id')))].sort(),
    fill: [...new Set((svgId(fill) ? getFeatureFillElements(svg, svgId(fill)) : []).map((element) => element.getAttribute('fill')))],
    label: (svgId(label) ? [...svg.querySelectorAll(`text[data-label-feature-id="${CSS.escape(svgId(label))}"]`)] : [])
      .filter(drawn).map((element) => element.textContent.trim()),
    hidden: elements(hide).filter(drawn).length === 0,
    stroke: [...new Set(elements(stroke).map((element) => element.getAttribute('stroke')))]
  };
}, { ids: targets, key: recordKey });

const allShown = (records, labelText = LABEL) => ({
  records, fill: [FILL], label: [labelText], hidden: true, stroke: [STROKE]
});

// The record keys of the drawing's draft rows naming the target features.
const drafts = (page, targets) => page.evaluate(({ fill, label, hide, stroke }) => {
  const app = window.__GBDRAW_APP__;
  const overrides = Object.values(app.featureOverrides);
  const keyOf = (id, field) => overrides.filter((row) => row.biologicalFeatureId === id && row[field] !== null)
    .map((row) => row.recordKey);
  return {
    fill: app.manualSpecificRules.filter((rule) => rule.qual === 'hash' && String(rule.val).includes(fill))
      .map((rule) => String(rule.val).replace(/^\^|\$$/g, '').split('@').slice(1).join('@')),
    label: keyOf(label, 'labelText'),
    hide: keyOf(hide, 'featureVisibility'),
    stroke: Object.keys(app.featureStrokeOverrides).filter((key) => key.endsWith(`\0${stroke}`))
      .map((key) => key.slice(0, key.indexOf('\0')))
  };
}, targets);

const keyedUnder = (key) => ({ fill: [key], label: [key], hide: [key], stroke: [key] });

// The fill of the feature `featureId` of the record keyed `recordKey`, read
// through the SVG ID that feature has in the displayed Result (a duplicate
// record's features have IDs of their own).
const fillOf = (page, recordKey, featureId) => page.evaluate(async ({ key, id }) => {
  const app = window.__GBDRAW_APP__;
  const svg = app.svgContainer.querySelector('svg');
  const { getFeatureFillElements } = await import('/gbdraw/web/js/app/feature-editor/svg-actions.js');
  const feature = app.extractedFeatures.find((candidate) => candidate.recordKey === key && candidate.biologicalFeatureId === id);
  if (!feature) return `no feature ${id} in record ${key}`;
  return [...new Set(getFeatureFillElements(svg, String(feature.svg_id)).map((element) => element.getAttribute('fill')))];
}, { key: recordKey, id: featureId });

// Canvas ON, the four edits on record 2, then canvas OFF with record 2 drawn
// alone and Generate.
const editOnCanvasThenDrawRecord2Alone = async (page) => {
  await startCircular(page);
  await generate(page);
  const targets = await editFour(page, CANVAS_KEY);
  expect(await shown(page, targets)).toEqual(allShown(['BGC0000708', 'BGC0000709']));
  expect(await drafts(page, targets)).toEqual(keyedUnder(CANVAS_KEY));
  await page.evaluate(() => { window.__GBDRAW_APP__.form.multi_record_canvas = false; });
  await page.getByLabel('Circular record', { exact: true }).selectOption('BGC0000709');
  await generate(page);
  return targets;
};

const undoOrRedo = async (page, action) => {
  await page.getByRole('button', { name: action, exact: true }).click();
  await historySettled(page);
};

const loadSession = async (page, path) => {
  page.on('dialog', (dialog) => dialog.accept());
  await openApp(page);
  await page.locator('input[accept^=".json,"]').first().setInputFiles(path);
  await expect.poll(() => page.evaluate(() => !window.__GBDRAW_APP__.sessionImportPending
    && window.__GBDRAW_APP__.extractedFeatures.length > 0), { timeout: 180_000 }).toBe(true);
};

// The text files of the reproducibility bundle that Run info downloads.
const bundleTexts = async (page, testInfo) => {
  await page.getByRole('button', { name: /Run info/i }).click();
  const pending = page.waitForEvent('download');
  await page.getByRole('button', { name: /Download reproducibility files/ }).click();
  const path = testInfo.outputPath('reproducibility.zip');
  await (await pending).saveAs(path);
  const listed = spawnSync('python3', ['-c', [
    'import json, sys, zipfile',
    'z = zipfile.ZipFile(sys.argv[1])',
    'print(json.dumps({n: z.read(n).decode("utf-8", "replace") for n in z.namelist() if not n.endswith((".gb", ".gbk"))}))'
  ].join('\n'), path], { encoding: 'utf8' });
  expect(listed.status, listed.stderr).toBe(0);
  return JSON.parse(listed.stdout);
};

const forcedLabelRerenders = (page) => page.evaluate(async () => (await import('./js/state.js')).state.labelReflowForceRequestSeq.value);

test.describe('record identity moves (OV-1949)', () => {
  test.describe.configure({ timeout: 420_000 });

  test('Circular: canvas edits apply to record 2 drawn alone, live, and back on the canvas', async ({ page }) => {
    const targets = await editOnCanvasThenDrawRecord2Alone(page);
    expect(await shown(page, targets), 'record 2 drawn alone shows the canvas edits').toEqual(allShown(['BGC0000709']));
    expect(await drafts(page, targets)).toEqual(keyedUnder(SINGLE_KEY));

    // The Result Python drew already shows the moved Label visibility, so
    // opening a popup asks for no label rerender (which would remount the
    // Result and reset label text being typed in the popup).
    const rerenders = await forcedLabelRerenders(page);
    await openPopup(page, SINGLE_KEY, targets.label);
    await closePopup(page);
    expect(await forcedLabelRerenders(page), 'no label rerender after a Generate that moved the rows').toBe(rerenders);

    // Live: the label row moved, so a popup edit updates the Result at once.
    await editLabel(page, SINGLE_KEY, targets.label, LABEL_LIVE);
    await expect.poll(async () => (await shown(page, targets)).label, 'the live label edit shows').toEqual([LABEL_LIVE]);

    await page.evaluate(() => { window.__GBDRAW_APP__.form.multi_record_canvas = true; });
    await generate(page);
    expect(await shown(page, targets), 'back on the canvas').toEqual(allShown(['BGC0000708', 'BGC0000709'], LABEL_LIVE));
    expect(await drafts(page, targets)).toEqual(keyedUnder(CANVAS_KEY));
  });

  test('Circular: Undo of the view-switch Generate restores the canvas Result and rows; Redo restores record 2 alone', async ({ page }) => {
    const targets = await editOnCanvasThenDrawRecord2Alone(page);
    await undoOrRedo(page, 'Undo');
    expect(await shown(page, targets), 'Undo shows the canvas Result').toEqual(allShown(['BGC0000708', 'BGC0000709']));
    expect(await drafts(page, targets), 'Undo puts the rows back under the canvas key').toEqual(keyedUnder(CANVAS_KEY));
    await undoOrRedo(page, 'Redo');
    expect(await shown(page, targets), 'Redo shows record 2 alone with the edits').toEqual(allShown(['BGC0000709']));
    expect(await drafts(page, targets)).toEqual(keyedUnder(SINGLE_KEY));
  });

  test('Linear: edits on a card drawing record 2 wait while it draws record 1 and come back with record 2', async ({ page }) => {
    const [, uid] = await startLinear(page);
    await generate(page);
    const targets = await editFour(page, uid);
    expect(await shown(page, targets, uid)).toEqual(allShown(['BGC0000708', 'BGC0000709']));
    expect(await drafts(page, targets)).toEqual(keyedUnder(uid));

    await selectLinearRecord(page, uid, 'BGC0000708');
    await generate(page);
    const record1 = await page.evaluate(({ fill, stroke, label }) => {
      const content = window.__GBDRAW_APP__.currentResults()[0].content.toLowerCase();
      return { fill: content.includes(`fill="${fill}"`), stroke: content.includes(`stroke="${stroke}"`), label: content.includes(label.toLowerCase()) };
    }, { fill: FILL, stroke: STROKE, label: LABEL });
    expect(record1, 'record 1 shows none of the record 2 edits').toEqual({ fill: false, stroke: false, label: false });
    expect((await shown(page, targets)).records).toEqual(['BGC0000708']);
    expect(await drafts(page, targets), 'the rows wait under record 2 of the card\'s file').toEqual(keyedUnder(`${uid}:2`));

    await selectLinearRecord(page, uid, 'BGC0000709');
    await generate(page);
    expect(await shown(page, targets, uid), 'record 2 again').toEqual(allShown(['BGC0000708', 'BGC0000709']));
    expect(await drafts(page, targets)).toEqual(keyedUnder(uid));

    // Live: the moved label row takes a popup edit at once.
    await editLabel(page, uid, targets.label, LABEL_LIVE);
    await expect.poll(async () => (await shown(page, targets, uid)).label, 'the live label edit shows').toEqual([LABEL_LIVE]);
  });

  // GT B11 (PD-OI-104 rev 5, review 2 M4): a live Depth hide replaces the
  // committed Session's resources; the record identities the Generate read stay
  // with them, so the next Generate moves the rows by those identities and never
  // reads the committed records again against the current sources (the fallback).
  test('GT B11: after a live Depth hide the view-switch Generate moves the rows by the identities of the last Generate', async ({ page }) => {
    await page.addInitScript(() => {
      const fallbacks = window.__RI_FALLBACKS__ = [];
      window.__GBDRAW_TEST_HOOKS__ = {
        onStructuralMetric: (metric) => { if (metric.name === 'recordIdentityFallbackCount') fallbacks.push(metric); }
      };
    });
    await startCircular(page);
    const depth = ['BGC0000708', 'BGC0000709']
      .flatMap((id) => [0, 1, 2, 3].map((index) => `${id}\t${index * 9000 + 1}\t${10 + index}`)).join('\n');
    await page.evaluate((text) => {
      const app = window.__GBDRAW_APP__;
      app.setCircularDepthFile(0, new File([text], 'depth.tsv', { type: 'text/tab-separated-values' }));
      app.setTrackShown('show_depth', true);
    }, depth);
    await generate(page);
    const targets = await pickTargets(page, CANVAS_KEY);
    await editFill(page, CANVAS_KEY, targets.fill);
    const depthGroups = () => page.evaluate(async () => {
      const drawnGroups = await import('/gbdraw/web/js/services/drawn-groups.js');
      const svg = window.__GBDRAW_APP__.svgContainer.querySelector('svg');
      return drawnGroups.findTrackGroups(svg, drawnGroups.DEPTH_GROUPS).filter((group) => group.getAttribute('display') !== 'none').length;
    });
    expect(await depthGroups(), 'the Generate draws Depth').toBeGreaterThan(0);
    // A History step, as the Show Depth checkbox makes one: the live apply runs from it.
    await page.evaluate(() => window.__GBDRAW_HISTORY__.runUndoable('Show Depth',
      () => window.__GBDRAW_APP__.setTrackShown('show_depth', false)));
    await settled(page);
    await expect.poll(depthGroups, { message: 'the live hide removes Depth', timeout: 120_000 }).toBe(0);
    await page.evaluate(() => { window.__GBDRAW_APP__.form.multi_record_canvas = false; });
    await page.getByLabel('Circular record', { exact: true }).selectOption('BGC0000709');
    await generate(page);
    expect(await page.evaluate(() => window.__RI_FALLBACKS__.length), 'no Generate read the record identities again').toBe(0);
    expect(await fillOf(page, SINGLE_KEY, targets.fill), 'record 2 drawn alone keeps the canvas fill').toEqual([FILL]);
    expect((await drafts(page, targets)).fill).toEqual([SINGLE_KEY]);
  });

  test('identical features: a "This feature only" fill stays on its own copy in every view', async ({ page }) => {
    // No labels: the twin CDS share one place, so their labels would overlap.
    await startCircular(page, identicalTwins, { labels_mode: 'none' });
    await generate(page);
    const twins = await page.evaluate(() => {
      const features = window.__GBDRAW_APP__.extractedFeatures
        .filter((feature) => feature.recordKey === 'record-2' && /~\d+$/.test(feature.biologicalFeatureId));
      return features.map((feature) => feature.biologicalFeatureId).sort();
    });
    expect(twins).toHaveLength(2);
    expect(new Set(twins.map((id) => id.split('~')[0])).size).toBe(1);
    await editFill(page, 'record-2', twins[1]);
    const fills = async (recordKey) => ({ edited: await fillOf(page, recordKey, twins[1]), other: (await fillOf(page, recordKey, twins[0])).includes(FILL) });
    expect(await fills('record-2')).toEqual({ edited: [FILL], other: false });

    await page.evaluate(() => { window.__GBDRAW_APP__.form.multi_record_canvas = false; });
    await page.getByLabel('Circular record', { exact: true }).selectOption('TESTB');
    await generate(page);
    expect(await fills('circular-TESTB-#2'), 'record 2 drawn alone').toEqual({ edited: [FILL], other: false });
    expect((await drafts(page, { fill: twins[1] })).fill).toEqual(['circular-TESTB-#2']);

    await page.evaluate(() => { window.__GBDRAW_APP__.form.multi_record_canvas = true; });
    await generate(page);
    expect(await fills('record-2'), 'back on the canvas').toEqual({ edited: [FILL], other: false });
  });

  test('Linear duplicates: the same file in two cards keeps a "This feature only" fill on its own card', async ({ page }) => {
    // Cards: file 1 record 1, file 1 record 2 (`first`), file 2 record 1, file 2 record 2 (`twin`).
    const [, first, , twin] = await startLinear(page, 2);
    await generate(page);
    const { fill } = await pickTargets(page, first);
    await editFill(page, first, fill);
    await generate(page);
    expect({ first: await fillOf(page, first, fill), twin: (await fillOf(page, twin, fill)).includes(FILL) },
      'the edited card is filled, its twin is not').toEqual({ first: [FILL], twin: false });
    expect((await drafts(page, { fill })).fill).toEqual([first]);

    await selectLinearRecord(page, first, 'BGC0000708');
    await generate(page);
    expect((await fillOf(page, twin, fill)).includes(FILL), 'the twin stays unfilled').toBe(false);
    expect((await drafts(page, { fill })).fill).toEqual([`${first}:2`]);

    await selectLinearRecord(page, first, 'BGC0000709');
    await generate(page);
    expect({ first: await fillOf(page, first, fill), twin: (await fillOf(page, twin, fill)).includes(FILL) },
      'back on record 2, only the edited card is filled').toEqual({ first: [FILL], twin: false });
  });

  test('a 0.14.2 Session of record 2 drawn alone keeps its edits on the canvas', async ({ page }) => {
    const provenance = JSON.parse(readFileSync(session0142.replace('.v48.gbdraw-session.json.gz', '.provenance.json'), 'utf8'));
    const { edits } = provenance;
    const targets = { fill: edits.fill.featureId, label: edits.label.featureId, hide: edits.hide.featureId, stroke: edits.stroke.featureId };
    await loadSession(page, session0142);
    expect(await drafts(page, targets)).toEqual(keyedUnder(provenance.recordKey));
    await page.evaluate(() => { window.__GBDRAW_APP__.form.multi_record_canvas = true; });
    await generate(page);
    expect(await shown(page, targets), 'the canvas draws the 0.14.2 edits on record 2').toEqual({
      records: ['BGC0000708', 'BGC0000709'], fill: [edits.fill.color], label: [edits.label.text], hidden: true, stroke: [edits.stroke.color]
    });
    expect(await drafts(page, targets)).toEqual(keyedUnder(CANVAS_KEY));
  });

  test('Run info after the view switch: the Source recipe names the fill by the CLI record key', async ({ page }, testInfo) => {
    const targets = await editOnCanvasThenDrawRecord2Alone(page);
    const recipe = await page.evaluate(() => {
      const source = window.__GBDRAW_APP__.lastRunInfo?.sourceRecipe || {};
      return { available: source.available, reason: source.unavailableReason || '', command: source.command || '' };
    });
    expect(recipe.available, recipe.reason).toBe(true);
    expect(recipe.command).toContain(' -t ');
    const table = Object.values(await bundleTexts(page, testInfo)).find((text) => /\thash\t/.test(text));
    expect(table, 'the bundle holds the specific color table').toBeTruthy();
    const row = table.split('\n').find((line) => line.includes(targets.fill));
    expect(row, table).toMatch(new RegExp(`^CDS\thash\t\\^?${targets.fill}@record-1\\$?\t${FILL}\t`, 'i'));
  });
});
