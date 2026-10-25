const { expect, test } = require('@playwright/test');
const fs = require('node:fs/promises');
const { gunzipSync } = require('node:zlib');
const { answerSessionTitle, assertOperationHealth, openApp, readErrorSignature } = require('./app-lifecycle.cjs');

const seeds = {
  circular: 'gbdraw/web/gallery/sessions/HmmtDNA_basic_circular.gbdraw-session.json',
  linear: 'gbdraw/web/gallery/sessions/lambda_basic_linear.gbdraw-session.json'
};

const load = async (browser, file = seeds.circular, viewport = { width: 1600, height: 1000 }) => {
  const context = await browser.newContext({ viewport });
  const page = await context.newPage();
  page.externalRequests = [];
  await context.route('**/*', route => {
    if (new URL(route.request().url()).hostname === '127.0.0.1') return route.continue();
    page.externalRequests.push(route.request().url());
    return route.abort();
  });
  await page.addInitScript(() => {
    window.__MODE_EVENTS__ = [];
    window.__GBDRAW_TEST_HOOKS__ = {
      onSessionLifecycleEvent: event => window.__MODE_EVENTS__.push(event)
    };
  });
  page.on('dialog', dialog => dialog.dismiss());
  await openApp(page);
  await page.locator('input[accept^=".json,"]').setInputFiles(file);
  await expect.poll(() => page.evaluate(() => !window.__GBDRAW_APP__.sessionImportPending
    && window.__GBDRAW_APP__.extractedFeatures.length > 0), { timeout: 180000 }).toBe(true);
  return page;
};

// G-G(3): options name the audit ID of a known defect that trips a health check.
const generate = async (page, health = {}) => {
  const key = await page.evaluate(async () => (await import('./js/state.js')).state.resultGenerationKey.value);
  await page.getByRole('button', { name: 'Generate Diagram', exact: true }).click();
  await expect.poll(() => page.evaluate(async () => {
    const { state } = await import('./js/state.js');
    return { key: state.resultGenerationKey.value, processing: state.processing.value, error: state.errorLog.value };
  }), { timeout: 180000 }).toEqual({ key: key + 1, processing: false, error: null });
  await assertOperationHealth(page, { operation: 'Generate', ...health });
};

const switchMode = async (page, mode) => {
  await page.getByRole('button', { name: mode === 'circular' ? 'Circular' : 'Linear', exact: true }).click();
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.mode)).toBe(mode);
  await expect.poll(() => page.evaluate(async () =>
    (await import('./js/state.js')).state.circularRecordDiscovery.status)).not.toBe('loading');
  // Settle the Vue render and its asynchronous preview binder, without a time delay.
  await page.evaluate(async () => { await window.Vue.nextTick(); });
};

const popup = async (page, index = 0) => {
  if (!await page.locator('.right-drawer').isVisible()) await page.locator('.drawer-toggle').click();
  await page.locator('.right-drawer').getByRole('button', { name: 'Edit', exact: true }).nth(index).click();
  await expect(page.locator('.feature-popup')).toBeVisible();
  return page.evaluate(() => {
    const feature = window.__GBDRAW_APP__.clickedFeature;
    return { id: feature.id, featureId: feature.svg_id, sourceText: feature.labelSourceText };
  });
};

const closeEditor = async page => {
  if (await page.locator('.feature-popup').isVisible()) {
    await page.getByRole('button', { name: 'Close feature popup', exact: true }).click();
  }
  if (await page.locator('.right-drawer').isVisible()) await page.locator('.drawer-toggle').click();
};

// `title`: the Session has none, so Save asks for one; the dialog is answered
// with it ('' saves an untitled file). Leave it out for a titled Session.
const download = async (page, button, path, { title, ...health } = {}) => {
  const control = button === 'Save Session' ? page.getByRole('banner') : page;
  const errorSignatureBefore = await readErrorSignature(page);
  const pending = page.waitForEvent('download');
  await control.getByRole('button', { name: button, exact: true }).click();
  if (title !== undefined) await answerSessionTitle(page, title);
  // A large file asks first (`Download large ...?`), at any time before the
  // download starts; the download action answers it.
  const large = page.locator('[data-choice-prompt] [data-choice-prompt-action="download"]');
  const asked = large.waitFor({ state: 'visible', timeout: 0 }).then(() => 'asked', () => null);
  const first = await Promise.race([pending, asked]);
  if (first === 'asked') await large.click();
  const download = first === 'asked' ? await pending : first;
  await download.saveAs(path);
  if (button === 'Save Session') {
    await assertOperationHealth(page, { operation: 'Save Session', errorSignatureBefore, ...health });
  }
  return fs.readFile(path);
};

// Legend rows the editor added (`featureIds: []`, captions Python's inventory
// lacks), as a Session saved while Add legend item existed holds them; R15-2
// retired that action, so a Session is their only source. The page's Session
// is saved; `rows` ([caption, color]) are written into the shown drawing's
// Legend entries and, as that writer drew them, into its Results (a copy of
// the last row of each feature Legend group, owned by the editor, a unit lower
// per row so the editor lists them last, in order); the Session is loaded into
// the page and generated, which lays the rows out as Python does.
const loadEditorLegendRows = async (page, rows) => {
  const path = test.info().outputPath(`editor-legend-rows-${rows.map(([caption]) => caption).join('-')
    .replace(/[^A-Za-z0-9-]+/g, '_').slice(0, 60)}-${Date.now()}.gbdraw-session.json`);
  // Saved under a title of its own, so a later Save of the page's title is no
  // second download of one file name (which asks first); loaded with the title.
  const { mode, title } = await page.evaluate(async () => {
    const { state } = await import('./js/state.js');
    const shown = state.sessionTitle.value || 'editor-legend-rows';
    state.sessionTitle.value = `${shown} editor rows`;
    return { mode: state.mode.value, title: shown };
  });
  const bytes = await download(page, 'Save Session', path);
  const saved = JSON.parse((bytes[0] === 0x1f ? gunzipSync(bytes) : bytes).toString('utf8'));
  saved.title = title;
  expect(saved.renderRequest.mode, 'the saved Results are the shown drawing\'s').toBe(mode);
  saved.modes[mode].editorState.legend.entries.push(...rows.map(([caption, color]) => (
    { caption, originalCaption: caption, color, featureIds: [] })));
  saved.results = await page.evaluate(async ({ results, added }) => {
    const { getAllFeatureLegendGroups, getLegendEntrySwatch } = await import('./js/services/legend-svg.js');
    return results.map((result) => {
      const svg = new DOMParser().parseFromString(result.content, 'image/svg+xml').documentElement;
      getAllFeatureLegendGroups(svg).forEach((group) => {
        const last = [...group.querySelectorAll('g[data-legend-key]')].at(-1);
        added.forEach(([caption, color], index) => {
          const row = last.cloneNode(true);
          [row, ...row.querySelectorAll('*')].forEach((element) => [...element.attributes]
            .filter(({ name }) => name === 'display' || name.startsWith('data-gbdraw-base-'))
            .forEach(({ name }) => element.removeAttribute(name)));
          row.setAttribute('data-legend-key', caption);
          row.setAttribute('data-legend-owner', 'direct-editor');
          row.setAttribute('transform', `translate(0,${index + 1}) ${row.getAttribute('transform') || ''}`.trim());
          getLegendEntrySwatch(row).setAttribute('fill', color);
          row.querySelector('text').textContent = caption;
          last.parentNode.appendChild(row);
        });
      });
      return { ...result, content: new XMLSerializer().serializeToString(svg) };
    });
  }, { results: saved.results, added: rows });
  await fs.writeFile(path, JSON.stringify(saved));
  await page.locator('input[accept^=".json,"]').setInputFiles(path);
  await expect.poll(() => page.evaluate(captions => {
    const app = window.__GBDRAW_APP__;
    return !app.sessionImportPending && captions.every(caption => app.legendEntries.some(entry => entry.caption === caption));
  }, rows.map(([caption]) => caption)), { timeout: 180000 }).toBe(true);
  await generate(page);
};

const snapshot = page => page.evaluate(async () => {
  const { state: s } = await import('./js/state.js');
  const ingestion = await import('./js/services/svg-result-ingestion.js');
  const { getCommittedCanonicalRenderRequest } = await import('./js/services/config.js');
  // A mode without its own Result shows the empty Preview (E1).
  const root = s.svgContainer.value?.querySelector('svg') || null;
  const result = window.__GBDRAW_APP__.currentResults()[s.selectedResultIndex.value] || null;
  return {
    mode: s.mode.value, generation: s.resultGenerationKey.value,
    // Per-feature edits by identity key, one map per edited field.
    ...Object.fromEntries([['labels', 'labelText'], ['labelSources', 'labelSourceText'], ['visibility', 'labelVisibility'],
      ['featureVisibility', 'featureVisibility']].map(([name, field]) => [name, Object.fromEntries(
      Object.entries(s.activeDrawing().featureOverrides).filter(([, row]) => row[field] !== null).map(([key, row]) => [key, row[field]]))])),
    bulkLabels: { ...s.activeDrawing().labelTextBulkOverrides }, visibilityRules: [...s.activeDrawing().featureVisibilityManualRules],
    colors: { ...s.activeDrawing().featureColorOverrides }, rules: s.activeDrawing().manualSpecificRules.map(rule => ({ ...rule, fromFile: Boolean(rule.fromFile) })),
    featureCount: s.extractedFeatures.value.length,
    resultIdentity: result ? ingestion.getCommittedSvgResultRuntimeIdentity(result) : null,
    markedMounted: result ? ingestion.isCommittedSvgResultMounted(result) : false,
    result: result?.content ?? null, payload: s.svgContent.value, mounted: root?.outerHTML ?? null,
    sameRoot: Boolean(root) && root === window.__MODE_EDITED_ROOT__,
    mountEvents: window.__MODE_EVENTS__.filter(e => e.name === 'preview.mount-observed'),
    request: getCommittedCanonicalRenderRequest()
  };
});

const semantics = (page, content) => page.evaluate(content => {
  const root = new DOMParser().parseFromString(content, 'image/svg+xml').documentElement;
  return [...root.querySelectorAll('[data-gbdraw-feature-id]')].map(element => ({
    id: element.getAttribute('data-gbdraw-feature-id'),
    part: element.getAttribute('data-gbdraw-feature-part'),
    d: element.getAttribute('d'), fill: element.getAttribute('fill'),
    stroke: element.getAttribute('stroke'), display: element.getAttribute('display')
  }));
}, content);

module.exports = { seeds, load, generate, switchMode, popup, closeEditor, download, loadEditorLegendRows, snapshot, semantics };
