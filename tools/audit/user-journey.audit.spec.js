// User journeys for the dev to main promotion (docs/internal/WEB_PERIODIC_AUDIT.md,
// Procedure step 4). Each journey is one independent test that walks a user flow
// on realistic inputs (Gallery Sessions, Session fixtures written by main and by
// the 0.13.0 release, and the tutorial GenBank file) and checks, after every
// step, the common oracles of tools/audit/helpers/journey-evidence.cjs: no page
// error, console error, or unhandled rejection, and no busy indicator left.
//
// Run all journeys, or one with --grep "J5 ":
//   GBDRAW_WEB_TEST_PORT=<port> npx playwright test -c tools/audit/playwright.audit.config.cjs user-journey
// Evidence: $GBDRAW_AUDIT_OUT/journeys/<journey>/ (a screenshot per step and
// steps.json) and journeys/index.html (the contact sheet, failures first).
// The run also prints the capabilities changed since origin/main that no
// journey covers (COVERAGE below); the promotion PR lists them as uncovered.
const { execFileSync } = require('node:child_process');
const { readFileSync, rmSync, statSync, writeFileSync } = require('node:fs');
const { join } = require('node:path');
const { pathToFileURL } = require('node:url');
const { test, expect } = require('@playwright/test');
const A = require('./helpers/audit-common.cjs');
const { createJourney, writeContactSheet } = require('./helpers/journey-evidence.cjs');

const repo = (path) => join(A.REPO_ROOT, path);
const {
  CURRENT_SESSION_VERSION, assertOperationHealth, diffUserOwnedState, evaluateWithRetainedPromise,
  openApp, reveal, snapshotUserOwnedState
} = A.helpers();
const { expectLiveEqualsGenerate, semanticSnapshot, settleLive } = require(repo('tests/web/helpers/live-generate-parity.cjs'));
const { colorLegendRow, generate, history, renameRow } = require(repo('tests/web/helpers/live-generate-parity-steps.cjs'));
const { compareSvgFiles } = require(repo('tests/web/helpers/svg-semantic-compare.cjs'));

// Journey -> the capabilities (tools/ci-impact-policy.mjs classifyChanges) it
// exercises. A changed capability that no journey lists is printed as uncovered.
const COVERAGE = Object.freeze({
  J1: ['web-runtime', 'python-core', 'renderer'],
  J2: ['web-runtime', 'python-core', 'renderer'],
  J3: ['session-persistence', 'gallery', 'web-runtime'],
  J4: ['web-runtime'],
  J5: ['web-runtime'],
  // Comparison owners through the Session's cached LOSATP results; LOSAT
  // execution itself is the release tier's `losat-cache-browser-acceptance`.
  J6: ['losat-integration', 'gallery', 'session-persistence', 'web-runtime'],
  J7: ['web-runtime', 'session-persistence']
});
// Capabilities with no user flow to walk.
const NOT_USER_FACING = Object.freeze(['none', 'metadata', 'documentation', 'policy-documentation', 'tests-only', 'ci-only']);

const HMMT_GBK = repo('tests/test_inputs/HmmtDNA.gbk');
const GALLERY = (name) => repo(`gbdraw/web/gallery/sessions/${name}`);
const FIXTURE = (name) => repo(`tests/fixtures/sessions/${name}`);
const PHONE = { width: 390, height: 844 };

const changedCapabilities = async () => {
  let diff;
  try {
    diff = execFileSync('git', ['diff', '--name-status', '--no-renames', 'origin/main...HEAD'], {
      cwd: A.REPO_ROOT, encoding: 'utf8', maxBuffer: 64 * 1024 * 1024
    });
  } catch (error) {
    return { error: `git diff origin/main...HEAD failed; fetch origin main first (${String(error?.message || error).split('\n')[0]})` };
  }
  const changes = diff.split('\n').filter(Boolean).map((line) => {
    const [status, ...paths] = line.split('\t');
    return { status, paths };
  });
  if (!changes.length) return { changed: [], paths: {} };
  const { classifyChanges } = await import(pathToFileURL(repo('tools/ci-impact-policy.mjs')).href);
  const classified = classifyChanges(changes);
  const paths = {};
  for (const { path, impact } of classified.paths) (paths[impact] ||= []).push(path);
  return { changed: [...classified.capabilities], paths };
};

test.beforeAll(async () => {
  const { changed, paths, error } = await changedCapabilities();
  const covered = new Set(Object.values(COVERAGE).flat());
  const uncovered = error ? [] : changed.filter((capability) => !covered.has(capability) && !NOT_USER_FACING.includes(capability));
  const report = error
    ? { error }
    : {
      changed,
      covered: changed.filter((capability) => covered.has(capability)),
      uncovered,
      uncoveredPaths: Object.fromEntries(uncovered.map((capability) => [capability, paths[capability] || []])),
      notUserFacing: changed.filter((capability) => NOT_USER_FACING.includes(capability)),
      coverage: COVERAGE
    };
  A.writeEvidence(A.outDir('journeys'), 'coverage.json', report);
  // eslint-disable-next-line no-console
  console.log(error
    ? `user journeys: changed capabilities unknown: ${error}`
    : `user journeys: changed capabilities since origin/main: ${changed.join(', ') || 'none'}\n`
      + `  uncovered (hand look with a time limit, or an Owner waiver): ${uncovered.join(', ') || 'none'}\n`
      + uncovered.map((capability) => `    ${capability}: ${report.uncoveredPaths[capability].join(', ')}\n`).join('')
      + `  not user-facing: ${report.notUserFacing.join(', ') || 'none'}`);
});

test.afterAll(() => {
  writeContactSheet(A.outDir('journeys'));
});

// One journey: a fresh evidence folder, the step recorder, and the pages it opened.
const journey = (id, title, minutes, body) => test(`${id} ${title}`, async ({ page, browser, baseURL }) => {
  test.setTimeout(minutes * 60_000);
  const dir = A.outDir('journeys', id);
  rmSync(dir, { recursive: true, force: true });
  A.outDir('journeys', id);
  const steps = createJourney({ id, title, dir });
  const contexts = [];
  const freshPage = async (label, viewport = null) => {
    const context = await browser.newContext({ baseURL, acceptDownloads: true, ...(viewport ? { viewport } : {}) });
    contexts.push(context);
    return steps.watch(await context.newPage(), label);
  };
  try {
    await body({ page, steps, dir, freshPage });
    steps.finish();
  } finally {
    steps.close();
    for (const context of contexts) await context.close().catch(() => {});
  }
});

const displayedSvg = (page) => page.evaluate(() => {
  const app = window.__GBDRAW_APP__;
  return String(app.results?.[app.selectedResultIndex]?.content || '');
});

// The loaded preview, then Generate from the loaded draft: the two Results must
// be the same drawing (the Gallery publication parity comparison).
const generateEqualsLoadedPreview = async (page, dir, name) => {
  const loaded = await displayedSvg(page);
  expect(loaded, 'the Session shows a saved preview').toContain('<svg');
  await generate(page);
  const generated = await displayedSvg(page);
  const loadedPath = join(dir, `${name}-loaded.svg`);
  const generatedPath = join(dir, `${name}-generated.svg`);
  writeFileSync(loadedPath, loaded);
  writeFileSync(generatedPath, generated);
  const comparison = compareSvgFiles(loadedPath, generatedPath, { repoRoot: A.REPO_ROOT });
  expect(comparison.status, `Generate differs from the loaded preview:\n${comparison.report.slice(0, 4000)}`).toBe(0);
  return { loadedChars: loaded.length, generatedChars: generated.length };
};

const loadSession = async (page, file) => {
  const loaded = await A.loadSession(page, file);
  expect(loaded.error, JSON.stringify(loaded.error)).toBeNull();
  await settleLive(page);
  return { ms: loaded.ms, dialogs: loaded.dialogs };
};

const uploadCircular = async (page, path, { viaChooser = false } = {}) => {
  if (viaChooser) {
    const chooser = page.waitForEvent('filechooser');
    await page.getByRole('button', { name: 'Choose GenBank/DDBJ File', exact: true }).click();
    await (await chooser).setFiles(path);
  } else {
    await page.getByLabel('GenBank/DDBJ File', { exact: true }).setInputFiles(path);
  }
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.circularRecordList.length), { timeout: 60_000 })
    .toBeGreaterThan(0);
  await settleLive(page);
};

// The control lies inside the viewport width once scrolled to, and takes a click there.
const expectReachable = async (page, locator, name) => {
  await locator.scrollIntoViewIfNeeded();
  const box = await locator.boundingBox();
  const width = page.viewportSize().width;
  expect(box, `${name} is laid out`).not.toBeNull();
  expect(box.x >= -1 && box.x + box.width <= width + 1, `${name} lies inside the ${width} px viewport: ${JSON.stringify(box)}`)
    .toBe(true);
  await locator.click({ trial: true });
};

const expectNoHorizontalOverflow = async (page, label) => {
  const overflow = await A.measureOverflow(page);
  expect(overflow.overflowX, `${label}: the page scrolls horizontally: ${JSON.stringify(overflow.offenders)}`).toBe(false);
  return { scrollWidth: overflow.scrollWidth, clientWidth: overflow.clientWidth };
};

const parsesAsSvg = (page, source) => page.evaluate((text) => {
  const doc = new DOMParser().parseFromString(text, 'image/svg+xml');
  return !doc.querySelector('parsererror') && doc.documentElement.localName === 'svg';
}, source);

const exportFile = async (page, button, path) => {
  const control = page.getByRole('button', { name: button, exact: true });
  const pending = page.waitForEvent('download', { timeout: 120_000 });
  await control.click();
  const download = await pending;
  await download.saveAs(path);
  const bytes = statSync(path).size;
  expect(bytes, `${button} export is not empty`).toBeGreaterThan(0);
  return { file: download.suggestedFilename(), bytes };
};

// Starts an edit in the page and answers each choice dialog it opens by clicking
// the button `choose(title, texts)` names (default: the first that is not Cancel).
const runEdit = async (page, start, argument, choose = (title, texts) => texts.findIndex((text) => text !== 'Cancel')) => {
  await page.evaluate(`(() => {
    const edit = window.__JOURNEY_EDIT__ = { settled: false };
    edit.promise = Promise.resolve((${start.toString()})(${JSON.stringify(argument)}));
    edit.promise.then(() => { edit.settled = true; }, (error) => { edit.error = String(error?.message || error); edit.settled = true; });
  })()`);
  const answered = [];
  const openDialog = page.locator('[role="dialog"][aria-modal="true"]:visible');
  for (let round = 0; round < 6; round += 1) {
    await page.waitForFunction(() => window.__JOURNEY_EDIT__.settled
      || [...document.querySelectorAll('[role="dialog"][aria-modal="true"]')].some((dialog) => dialog.getClientRects().length > 0),
    null, { timeout: 120_000 });
    if (!await openDialog.count()) break;
    const dialog = openDialog.first();
    const title = (await dialog.locator('h2, h3').first().innerText()).trim();
    const texts = (await dialog.getByRole('button').allInnerTexts()).map((text) => text.replace(/\s+/g, ' ').trim());
    const index = choose(title, texts);
    if (index < 0 || index >= texts.length) throw new Error(`no choice for "${title}": ${JSON.stringify(texts)}`);
    answered.push(`${title}: ${texts[index]}`);
    await dialog.getByRole('button').nth(index).click();
    // The dialog stays open, marked busy, while the choice applies.
    await page.waitForFunction((shown) => ![...document.querySelectorAll('[role="dialog"][aria-modal="true"]')]
      .some((element) => element.getClientRects().length > 0 && element.querySelector('h2, h3')?.textContent.trim() === shown),
    title, { timeout: 120_000 });
  }
  const error = await evaluateWithRetainedPromise(page, async () => {
    await window.__JOURNEY_EDIT__.promise.catch(() => {});
    return window.__JOURNEY_EDIT__.error || null;
  });
  if (error) throw new Error(`the edit failed: ${error}`);
  await settleLive(page);
  return answered;
};

// Opens the feature popup of the `index`-th feature of `type` from the Features
// list, as its Edit button does.
const openPopup = (page, type, index) => evaluateWithRetainedPromise(page, async (target) => {
  const app = window.__GBDRAW_APP__;
  const feature = app.filteredFeatures.filter((item) => item.type === target.type)[target.index];
  if (!feature) throw new Error(`no ${target.type} feature #${target.index}`);
  await app.openFeatureEditorFromList(feature, null);
  await window.Vue.nextTick();
  if (!app.clickedFeature) throw new Error(`no popup for ${target.type} #${target.index}`);
  return [feature.type, feature.start, feature.end, feature.locus_tag || feature.product || ''].join('|');
}, { type, index });

const closePopup = async (page) => {
  const close = page.getByRole('button', { name: 'Close feature popup', exact: true });
  if (await close.isVisible()) await close.click();
  await settleLive(page);
};

const legendCaptions = (page) => page.evaluate(() => window.__GBDRAW_APP__.legendEntries.map((entry) => String(entry.caption || '')));

const historyDepth = (page) => page.evaluate(() => ({
  undo: window.__GBDRAW_HISTORY__.getUndoCount(),
  redo: window.__GBDRAW_HISTORY__.getRedoCount()
}));

const readDrawings = (page, paths) => page.evaluate(async (fieldPaths) => {
  const { state } = await import('/gbdraw/web/js/state.js');
  const unwrap = (value) => (value && typeof value === 'object' && value.__v_isRef ? value.value : value);
  const at = (drawing, path) => path.split('.').reduce((value, key) => unwrap(value)?.[key], drawing);
  const pick = (drawing) => Object.fromEntries(fieldPaths.map((path) => [path, JSON.parse(JSON.stringify(at(drawing, path) ?? null))]));
  return { mode: state.mode.value, circular: pick(state.drawings.circular), linear: pick(state.drawings.linear) };
}, paths);

const showMode = async (page, mode) => {
  await page.getByRole('button', { name: mode === 'linear' ? 'Linear' : 'Circular', exact: true }).click();
  await page.waitForFunction((value) => window.__GBDRAW_APP__?.mode === value, mode);
  await settleLive(page);
};

const pairwiseMatches = (svg) => (svg.match(/data-gbdraw-pairwise-match-id=/g) || []).length;
// The data-identity of each drawn comparison match.
const matchIdentities = (svg) => (svg.match(/<[^>]*data-gbdraw-pairwise-match-id=[^>]*>/g) || [])
  .map((tag) => Number(tag.match(/data-identity="([^"]*)"/)?.[1]))
  .filter(Number.isFinite);

// J1 and J7: upload, Generate with the button, and every export. At phone width
// each step also checks the controls it uses are reachable and nothing overflows.
const inputGenerateExport = async (page, steps, dir, { phone = false } = {}) => {
  const check = async (label, controls = []) => {
    if (!phone) return {};
    for (const [name, locator] of controls) await expectReachable(page, locator, name);
    return expectNoHorizontalOverflow(page, label);
  };
  await steps.step('open the app', async () => {
    await openApp(page);
    return check('open');
  }, { limitMs: 240_000 });
  await steps.step('upload HmmtDNA.gbk (Circular)', async () => {
    await check('before upload', [['Choose GenBank/DDBJ File', page.getByRole('button', { name: 'Choose GenBank/DDBJ File', exact: true })]]);
    await uploadCircular(page, HMMT_GBK, { viaChooser: true });
    return check('after upload');
  });
  await steps.step('Generate with the Generate button', async () => {
    const button = page.getByRole('button', { name: 'Generate Diagram', exact: true });
    await check('before Generate', [['Generate Diagram', button]]);
    await button.click();
    await page.waitForFunction(() => {
      const app = window.__GBDRAW_APP__;
      return !app.processing && String(app.results?.[0]?.content || '').includes('<svg');
    }, null, { timeout: 300_000 });
    await settleLive(page);
    await assertOperationHealth(page, { operation: 'Generate (button)' });
    expect(await page.evaluate(() => window.__GBDRAW_APP__.errorLog), 'Generate left an alert').toBeFalsy();
    return check('after Generate');
  }, { limitMs: 360_000 });
  const exported = {};
  for (const [button, name, verify] of [
    ['SVG', 'static.svg', async (path) => expect(await parsesAsSvg(page, readFileSync(path, 'utf8')), 'the SVG parses').toBe(true)],
    ['Interactive SVG', 'interactive.svg', async (path) => expect(await parsesAsSvg(page, readFileSync(path, 'utf8')), 'the SVG parses').toBe(true)],
    ['PNG', 'diagram.png', (path) => expect(readFileSync(path).subarray(0, 8).toString('hex'), 'PNG signature').toBe('89504e470d0a1a0a')],
    ['PDF', 'diagram.pdf', (path) => expect(readFileSync(path).subarray(0, 5).toString('latin1'), 'PDF header').toBe('%PDF-')]
  ]) {
    await steps.step(`export ${button}`, async () => {
      await check(`before ${button} export`, [[button, page.getByRole('button', { name: button, exact: true })]]);
      exported[button] = await exportFile(page, button, join(dir, name));
      await verify(join(dir, name));
      return exported[button];
    }, { limitMs: 180_000 });
  }
  return exported;
};

// The exported interactive SVG opened alone, as a reader opens the file: no page
// error, and its feature popup script runs on a click.
const openInteractiveAlone = async (steps, dir, freshPage, page) => {
  await steps.step('open the exported interactive SVG alone', async () => {
    const viewer = await freshPage('interactive SVG');
    steps.show(viewer);
    await viewer.goto(pathToFileURL(join(dir, 'interactive.svg')).href);
    const features = viewer.locator('[data-gbdraw-interactive-feature="true"]');
    await expect.poll(() => features.count(), { timeout: 30_000 }).toBeGreaterThan(0);
    const count = await features.count();
    await features.first().dispatchEvent('click');
    await viewer.waitForTimeout(500);
    return { interactiveFeatures: count };
  });
  steps.show(page);
};

journey('J1', 'Input -> Generate -> export', 15, async ({ page, steps, dir, freshPage }) => {
  await steps.watch(page);
  await inputGenerateExport(page, steps, dir);
  await openInteractiveAlone(steps, dir, freshPage, page);
});

journey('J2', 'Mode switch with per-mode settings', 15, async ({ page, steps }) => {
  await steps.watch(page);
  const fields = ['adv.label_font_size', 'adv.axis_stroke_width'];
  let initial;
  let circular;
  let linear;
  const shown = async () => ({
    drawings: await readDrawings(page, fields),
    result: await semanticSnapshot(page),
    legend: await legendCaptions(page)
  });
  await steps.step('open the app and upload HmmtDNA.gbk (Circular)', async () => {
    await openApp(page);
    await uploadCircular(page, HMMT_GBK);
    initial = await readDrawings(page, fields);
    return initial;
  }, { limitMs: 240_000 });
  await steps.step('Circular: label font size 13, Generate', async () => {
    await page.evaluate(() => { window.__GBDRAW_APP__.adv.label_font_size = 13; });
    return generate(page);
  }, { limitMs: 300_000 });
  await steps.step('Circular: rename the Legend row tRNA to "transfer RNA"', async () => {
    await renameRow(page, 'tRNA', 'transfer RNA');
    await settleLive(page);
    circular = await shown();
    expect(circular.legend).toContain('transfer RNA');
    return circular.drawings;
  });
  await steps.step('switch to Linear: its own file, axis stroke width 2, Generate', async () => {
    await showMode(page, 'linear');
    await page.evaluate(async (text) => {
      window.__GBDRAW_APP__.setLinearSeqPrimaryFile(0, 'gb', new File([text], 'HmmtDNA.gbk', { type: 'text/plain', lastModified: 1000 }));
      await window.Vue.nextTick();
    }, readFileSync(HMMT_GBK, 'utf8'));
    await settleLive(page);
    await page.evaluate(() => { window.__GBDRAW_APP__.adv.axis_stroke_width = 2; });
    await generate(page);
    linear = await shown();
    expect(linear.drawings.linear['adv.axis_stroke_width']).toBe(2);
    expect(linear.drawings.linear['adv.label_font_size'], 'the Circular font size stays in Circular')
      .toBe(initial.linear['adv.label_font_size']);
    expect(linear.drawings.circular['adv.axis_stroke_width'], 'the Linear stroke width stays in Linear')
      .toBe(initial.circular['adv.axis_stroke_width']);
    return linear.drawings;
  }, { limitMs: 300_000 });
  await steps.step('back to Circular: settings, Result, and Legend edit unchanged', async () => {
    await showMode(page, 'circular');
    const back = await shown();
    expect(back.drawings).toEqual({ ...circular.drawings, linear: linear.drawings.linear });
    expect(back.result, 'the Circular Result').toEqual(circular.result);
    expect(back.legend, 'the Circular Legend rows').toEqual(circular.legend);
  });
  await steps.step('back to Linear: settings and Result unchanged', async () => {
    await showMode(page, 'linear');
    const back = await shown();
    expect(back.drawings).toEqual({ ...linear.drawings, mode: 'linear' });
    expect(back.result, 'the Linear Result').toEqual(linear.result);
    expect(back.legend, 'the Linear Legend rows').toEqual(linear.legend);
  });
});

const J3_SESSIONS = [
  { key: 'v30-bgc', name: '0.13.0 Session 30 (BGC0000708-BGC0000713, Linear)', file: FIXTURE('BGC0000708-BGC0000713.v30.gbdraw-session.json.gz') },
  { key: 'v44-two-mode', name: 'main-written Session 44 (two-mode project)', file: FIXTURE('two-mode-project.v44.gbdraw-session.json.gz') },
  { key: 'gallery-hmmt', name: 'Gallery HmmtDNA_basic_circular (Circular)', file: GALLERY('HmmtDNA_basic_circular.gbdraw-session.json') },
  { key: 'gallery-bgc', name: 'Gallery BGC0000708-BGC0000713 (Linear, LOSATP comparisons)', file: GALLERY('BGC0000708-BGC0000713.gbdraw-session.json') }
];

journey('J3', 'Session round trip across versions', 30, async ({ steps, dir, freshPage }) => {
  // One failed Session does not stop the others; the journey fails at the end.
  for (const session of J3_SESSIONS) {
    const first = await freshPage(`${session.key} first`);
    steps.show(first);
    const savedPath = join(dir, `${session.key}.saved.gbdraw-session.json.gz`);
    const step = (name, body, limitMs = 600_000) => steps.step(`${session.key}: ${name}`, body, { limitMs, soft: true });
    const passed = () => steps.record.steps[steps.record.steps.length - 1].result === 'passed';
    let savedState = null;
    await step(`Load ${session.name}`, async () => {
      await openApp(first);
      return loadSession(first, session.file);
    });
    if (passed()) {
      // A mismatch here (the OV-278 class) does not stop the round trip below.
      await step('Generate equals the loaded preview', () => generateEqualsLoadedPreview(first, dir, `${session.key}-first`));
      await step('Save (writes the current Session version)', async () => {
        const state = await snapshotUserOwnedState(first);
        const result = await A.saveSession(first, savedPath);
        expect(result.doc.version, 'saved Session version').toBe(CURRENT_SESSION_VERSION);
        savedState = state;
        return { file: result.suggested, version: result.doc.version };
      });
    }
    if (savedState) {
      const second = await freshPage(`${session.key} reopened`);
      steps.show(second);
      await step('Load the saved file in a fresh context: user state equal', async () => {
        await openApp(second);
        const loaded = await loadSession(second, savedPath);
        const after = await snapshotUserOwnedState(second);
        const changes = diffUserOwnedState(savedState.state, after.state);
        expect(changes, JSON.stringify(changes).slice(0, 4000)).toEqual([]);
        return loaded;
      });
      await step('Generate equals the loaded preview again', () => generateEqualsLoadedPreview(second, dir, `${session.key}-reopened`));
      await second.context().close();
    }
    await first.context().close();
  }
});

journey('J4', 'Legend editing and History', 15, async ({ page, steps }) => {
  await steps.watch(page);
  const shown = async () => ({ result: await semanticSnapshot(page), legend: await legendCaptions(page) });
  let start;
  let floor;
  await steps.step('Load the Gallery Session HmmtDNA_basic_circular and Generate', async () => {
    await openApp(page);
    await loadSession(page, GALLERY('HmmtDNA_basic_circular.gbdraw-session.json'));
    await generate(page);
    start = await shown();
    floor = (await historyDepth(page)).undo;
    return { legend: start.legend, undoDepth: floor };
  }, { limitMs: 600_000 });
  const rowIndex = async (caption) => {
    const index = (await legendCaptions(page)).indexOf(caption);
    expect(index, `Legend row "${caption}"`).toBeGreaterThanOrEqual(0);
    return index;
  };
  const edits = [
    ['rename the row tRNA to "transfer RNA"', () => renameRow(page, 'tRNA', 'transfer RNA')],
    ['color the row rRNA #7b2cbf', () => colorLegendRow(page, 'rRNA', '#7b2cbf')],
    ['stroke the row CDS #e63946, width 2', async () => {
      const index = await rowIndex('CDS');
      await page.evaluate((row) => window.__GBDRAW_APP__.setLegendEntryStrokeColorValue(row, '#e63946'), index);
      await settleLive(page);
      await page.evaluate((row) => window.__GBDRAW_APP__.updateLegendEntryStrokeWidth(row, 2), index);
    }],
    ['delete the row GC content', async () => {
      const index = await rowIndex('GC content');
      await evaluateWithRetainedPromise(page, (row) => window.__GBDRAW_APP__.deleteLegendEntry(row), index);
    }],
    ['sort the rows Z to A', () => page.evaluate(() => window.__GBDRAW_APP__.sortLegendEntries('desc'))]
  ];
  for (const [name, edit] of edits) {
    await steps.step(name, async () => {
      await edit();
      await settleLive(page);
      return { legend: await legendCaptions(page) };
    });
  }
  let end;
  await steps.step('Undo every step: the drawing equals the start', async () => {
    end = await shown();
    let undone = 0;
    while ((await historyDepth(page)).undo > floor) {
      await history(page, 'undo');
      undone += 1;
    }
    const back = await shown();
    expect(back.legend, 'Legend rows after Undo').toEqual(start.legend);
    expect(back.result, 'drawing after Undo').toEqual(start.result);
    return { undone };
  });
  await steps.step('Redo every step: the drawing equals the end', async () => {
    let redone = 0;
    while ((await historyDepth(page)).redo > 0) {
      await history(page, 'redo');
      redone += 1;
    }
    const forward = await shown();
    expect(forward.legend, 'Legend rows after Redo').toEqual(end.legend);
    expect(forward.result, 'drawing after Redo').toEqual(end.result);
    return { redone };
  });
  await steps.step('live = Generate', () => expectLiveEqualsGenerate(page, { label: 'Legend edits after Undo and Redo' })
    .then(({ tolerated }) => ({ tolerated: tolerated.length })), { limitMs: 300_000 });
});

journey('J5', 'Feature editing and Reset', 20, async ({ page, steps }) => {
  await steps.watch(page);
  await steps.step('open the app, upload HmmtDNA.gbk (Circular), Generate', async () => {
    await openApp(page);
    await uploadCircular(page, HMMT_GBK);
    return generate(page);
  }, { limitMs: 360_000 });
  const COLORS = ['#2a9d8f', '#e63946', '#7b2cbf', '#f4a261', '#264653'];
  // The scope choice a Color Change Scope button names; "Apply to all "X" (n)"
  // is the rule or the Legend row the dialog offers.
  const scopeKind = (text) => [
    [/^This feature only/, 'single'], [/^Apply to all label /, 'displayLabel'],
    [/^Apply to all source label /, 'annotationLabel'], [/^Use existing /, 'useExisting'], [/^Apply to all "/, 'group']
  ].find(([pattern]) => pattern.test(text))?.[1] || text;
  // Each fill takes a scope choice no earlier fill took, until the dialog offers
  // no untaken choice. The choices a later dialog may stop offering go first.
  const PREFERENCE = ['useExisting', 'displayLabel', 'annotationLabel', 'single', 'group'];
  const rank = (text) => (PREFERENCE.includes(scopeKind(text)) ? PREFERENCE.indexOf(scopeKind(text)) : PREFERENCE.length);
  const taken = new Set();
  let untaken = true;
  for (let fill = 0; untaken && fill < COLORS.length; fill += 1) {
    await steps.step(`popup fill on CDS #${fill + 1} with an untaken scope choice; live = Generate`, async () => {
      const feature = await openPopup(page, 'CDS', fill);
      let offered = [];
      const answered = await runEdit(page, (color) => window.__GBDRAW_APP__.updateClickedFeatureColor(color), COLORS[fill],
        (title, texts) => {
          if (!/Change Scope/.test(title)) return texts.findIndex((text) => text !== 'Cancel');
          offered = texts.filter((text) => text !== 'Cancel');
          const pick = [...offered].sort((left, right) => rank(left) - rank(right))
            .find((text) => !taken.has(scopeKind(text))) || offered[0];
          taken.add(scopeKind(pick));
          untaken = offered.some((text) => !taken.has(scopeKind(text)));
          return texts.indexOf(pick);
        });
      if (!offered.length) untaken = false;
      await closePopup(page);
      const { tolerated } = await expectLiveEqualsGenerate(page, { label: `fill, ${answered.join('; ') || 'no scope dialog'}` });
      return { feature, offered, answered, tolerated: tolerated.length };
    }, { limitMs: 300_000 });
  }
  await steps.step('popup: visibility Off on tRNA #1; live = Generate', async () => {
    const feature = await openPopup(page, 'tRNA', 0);
    const answered = await runEdit(page, async () => {
      const app = window.__GBDRAW_APP__;
      app.clickedFeature.featureVisibility = 'off';
      return app.updateClickedFeatureVisibility('off');
    }, null);
    await closePopup(page);
    const { tolerated } = await expectLiveEqualsGenerate(page, { label: `visibility Off, ${answered.join('; ') || 'no dialog'}` });
    return { feature, answered, tolerated: tolerated.length };
  }, { limitMs: 300_000 });
  await steps.step('popup: label text on CDS #1; live = Generate', async () => {
    const feature = await openPopup(page, 'CDS', 0);
    const answered = await runEdit(page, async (text) => {
      const app = window.__GBDRAW_APP__;
      app.clickedFeature.labelText = text;
      return app.updateClickedFeatureLabelText();
    }, 'NADH dehydrogenase 1, edited');
    await closePopup(page);
    const { tolerated } = await expectLiveEqualsGenerate(page, { label: `label text, ${answered.join('; ') || 'no dialog'}` });
    return { feature, answered, tolerated: tolerated.length };
  }, { limitMs: 300_000 });
  await steps.step('Reset Settings (confirm)', async () => {
    await page.getByRole('button', { name: 'Reset Settings', exact: true }).click();
    await settleLive(page);
    return { dialogs: steps.record.dialogs.slice(-1) };
  });
  // PD-OI-066: right after Reset Settings, the shown Result and the Legend list
  // equal what the next Generate draws from the reset draft (OV-287, OV-289, OV-290).
  await steps.step('after Reset Settings: shown Result = Generate', async () => {
    const legendBefore = await legendCaptions(page);
    const { tolerated } = await expectLiveEqualsGenerate(page, { label: 'Reset Settings' });
    expect(await legendCaptions(page), 'the Legend list after Reset equals the list after Generate').toEqual(legendBefore);
    return { tolerated: tolerated.length };
  }, { limitMs: 300_000 });
});

journey('J6', 'Comparisons', 20, async ({ page, steps, dir }) => {
  await steps.watch(page);
  let before;
  await steps.step('Load the Gallery Session BGC0000708-BGC0000713 (Linear, LOSATP pairs)', async () => {
    await openApp(page);
    return loadSession(page, GALLERY('BGC0000708-BGC0000713.gbdraw-session.json'));
  }, { limitMs: 600_000 });
  await steps.step('Generate equals the loaded preview', async () => {
    const outcome = await generateEqualsLoadedPreview(page, dir, 'bgc');
    before = pairwiseMatches(await displayedSvg(page));
    expect(before, 'drawn comparison matches').toBeGreaterThan(0);
    return { ...outcome, matches: before };
  }, { limitMs: 600_000 });
  // Minimum identity is a Generate-time setting. No live-generate-parity spec
  // pins a comparison setting, and docs/REFERENCE/web-app.md ("The operation
  // labels below state when each kind of edit reaches the Result") lists no
  // comparison filter under Live edit: a setting outside that row stays in the
  // settings draft and "the Result keeps its applied settings" until the next
  // successful Generate. So the edit must leave the shown Result unchanged, and
  // Generate must apply it: a threshold above the lowest drawn identity removes
  // at least one match.
  let edited;
  await steps.step('raise Minimum identity with its control: the shown Result keeps its applied settings', async () => {
    const svg = await displayedSvg(page);
    const identities = matchIdentities(svg).sort((left, right) => left - right);
    expect(identities.length, 'drawn matches with an identity').toBe(before);
    const median = identities[Math.floor(identities.length / 2)];
    const target = Math.ceil(median) > identities[0] ? Math.ceil(median) : Math.floor(identities[0]) + 1;
    const current = Number(await page.evaluate(() => window.__GBDRAW_APP__.adv.identity));
    const shownBefore = await semanticSnapshot(page);
    const input = await reveal(page.getByLabel('Linear comparison minimum identity', { exact: true }));
    await expect(input, 'the Minimum identity control is reachable').toBeVisible();
    await input.fill(String(target));
    await input.press('Tab');
    await settleLive(page);
    expect(Number(await page.evaluate(() => window.__GBDRAW_APP__.adv.identity))).toBe(target);
    expect(pairwiseMatches(await displayedSvg(page)), 'drawn matches before Generate').toBe(before);
    expect(await semanticSnapshot(page), 'the shown Result before Generate').toEqual(shownBefore);
    const removable = identities.filter((value) => value < target).length;
    edited = { target, removable };
    return { identity: [current, target], lowestDrawn: identities[0], removable };
  }, { limitMs: 120_000 });
  await steps.step('Generate applies Minimum identity', async () => {
    await generate(page);
    const after = pairwiseMatches(await displayedSvg(page));
    expect(after, `matches at minimum identity ${edited.target} (was ${before}; ${edited.removable} drawn below it)`).toBeLessThan(before);
    return { matches: [before, after] };
  }, { limitMs: 600_000 });
  await steps.step('shown Result = Generate', () => expectLiveEqualsGenerate(page, { label: 'after the Minimum identity Generate' })
    .then(({ tolerated }) => ({ tolerated: tolerated.length })), { limitMs: 600_000 });
});

journey('J7', 'Phone width', 20, async ({ page, steps, dir }) => {
  await page.setViewportSize(PHONE);
  await steps.watch(page);
  await inputGenerateExport(page, steps, dir, { phone: true });
  await steps.step('popup fill on CDS #1 at phone width', async () => {
    await openPopup(page, 'CDS', 0);
    const popup = page.locator('.feature-popup');
    await expectReachable(page, popup.getByLabel('Feature fill color', { exact: true }), 'Feature fill color');
    await expectReachable(page, page.getByRole('button', { name: 'Close feature popup', exact: true }), 'Close feature popup');
    const answered = await runEdit(page, (color) => window.__GBDRAW_APP__.updateClickedFeatureColor(color), '#7b2cbf',
      (title, texts) => {
        const index = texts.findIndex((text) => text === 'This feature only');
        return index >= 0 ? index : texts.findIndex((text) => text !== 'Cancel');
      });
    await closePopup(page);
    return { answered, ...(await expectNoHorizontalOverflow(page, 'after the popup edit')) };
  }, { limitMs: 300_000 });
  await steps.step('Save Session at phone width', async () => {
    await expectReachable(page, page.getByRole('banner').getByRole('button', { name: 'Save Session', exact: true }), 'Save Session');
    const saved = await A.saveSession(page, join(dir, 'phone.gbdraw-session.json.gz'));
    expect(saved.doc.version).toBe(CURRENT_SESSION_VERSION);
    return { file: saved.suggested, ...(await expectNoHorizontalOverflow(page, 'after Save')) };
  }, { limitMs: 300_000 });
});
