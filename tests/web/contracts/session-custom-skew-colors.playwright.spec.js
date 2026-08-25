const { test, expect } = require('@playwright/test');
const { readFileSync } = require('node:fs');
const { join, resolve } = require('node:path');
const { gunzipSync } = require('node:zlib');
const { openApp } = require('../helpers/app-lifecycle.cjs');

test.describe.configure({ retries: 0 });

const repoRoot = resolve(process.env.GBDRAW_REPO || process.cwd());
const sourceSessionPath = join(
  repoRoot,
  'gbdraw',
  'web',
  'gallery',
  'sessions',
  'HmmtDNA_ATskew.gbdraw-session.json'
);
const sourceSession = JSON.parse(readFileSync(sourceSessionPath, 'utf8'));
const AT_SLOT_ID = 'a_skew_2';
const AT_POSITIVE = '#deaf6e';
const AT_NEGATIVE = '#7294e3';

const expectedAtParams = {
  nt: 'AT',
  positive_color: AT_POSITIVE,
  negative_color: AT_NEGATIVE,
  legend_label: 'AT skew'
};

const installRequestObserver = async (page) => {
  await page.addInitScript(() => {
    window.__GBDRAW_CUSTOM_SKEW_REQUESTS__ = [];
    const NativeWorker = window.Worker;
    window.Worker = new Proxy(NativeWorker, {
      construct(target, args) {
        const worker = Reflect.construct(target, args, target);
        if (!String(args[0] || '').includes('diagram-generation-worker.js')) return worker;
        const nativePostMessage = worker.postMessage.bind(worker);
        worker.postMessage = (message, transfer) => {
          if (message?.type === 'run' && message?.payload?.request) {
            window.__GBDRAW_CUSTOM_SKEW_REQUESTS__.push(
              JSON.parse(JSON.stringify(message.payload.request))
            );
          }
          if (transfer === undefined) return nativePostMessage(message);
          return nativePostMessage(message, transfer);
        };
        return worker;
      }
    });
  });
};

const collectErrors = (page) => {
  const errors = [];
  page.on('pageerror', (error) => errors.push(`page: ${error.message}`));
  page.on('console', (message) => {
    if (message.type() === 'error') errors.push(`console: ${message.text()}`);
  });
  return errors;
};

const loadSessionThroughUi = async (page, sessionPath) => {
  const chooserPromise = page.waitForEvent('filechooser');
  await page.getByRole('button', { name: 'Load Session' }).click();
  const chooser = await chooserPromise;
  const dialogPromise = page.waitForEvent('dialog', { timeout: 180_000 });
  await chooser.setFiles(sessionPath);
  const dialog = await dialogPromise;
  expect(dialog.message()).toBe('Session loaded successfully!');
  await dialog.accept();
  await expect.poll(() => page.evaluate(() => ({
    mode: window.__GBDRAW_APP__?.mode,
    processing: window.__GBDRAW_APP__?.processing,
    slotIds: (window.__GBDRAW_APP__?.adv?.circular_track_slots || [])
      .map((slot) => slot.id)
  })), { timeout: 180_000 }).toEqual({
    mode: 'circular',
    processing: false,
    slotIds: ['features', 'gc_content', 'gc_skew', AT_SLOT_ID, 'ticks']
  });
};

const openCustomTrackSlots = async (page) => {
  const panel = page.locator('#circular-custom-track-slots-panel');
  if (!await panel.isVisible()) {
    await page.getByRole('button', { name: 'Custom Track Slots' }).click();
  }
  await expect(panel).toBeVisible();
};

const expectExplicitAtControls = async (page) => {
  await openCustomTrackSlots(page);
  const group = page.getByRole('group', { name: `Circular track slot ${AT_SLOT_ID}` });
  await expect(group).toBeVisible();
  await expect(group.getByLabel('Track dinucleotide')).toHaveValue('AT');
  await expect(group.getByLabel('Track legend label')).toHaveValue('AT skew');
  await expect(group.locator('input[type="color"][title="Positive skew color"]'))
    .toHaveValue(AT_POSITIVE);
  await expect(group.locator('input[type="color"][title="Negative skew color"]'))
    .toHaveValue(AT_NEGATIVE);
  await expect(group.getByText('Custom +', { exact: true })).toBeVisible();
  await expect(group.getByText('Custom -', { exact: true })).toBeVisible();
  return group;
};

const activeAtParams = (page) => page.evaluate((slotId) => {
  const slot = window.__GBDRAW_APP__.adv.circular_track_slots
    .find((candidate) => candidate.id === slotId);
  return JSON.parse(JSON.stringify(slot?.params || null));
}, AT_SLOT_ID);

const requestEvidence = (request) => {
  const slots = request?.diagramOptions?.tracks?.circularTrackSlots || [];
  const at = slots.find((slot) => slot.id === AT_SLOT_ID);
  const gc = slots.find((slot) => slot.id === 'gc_skew');
  return {
    schema: request?.schema,
    at: at && {
      id: at.id,
      renderer: at.renderer,
      params: at.params
    },
    gc: gc && {
      id: gc.id,
      renderer: gc.renderer,
      params: gc.params
    }
  };
};

const expectRenderedAtColors = (evidence) => {
  expect(evidence.result).toEqual(evidence.mounted);
  expect(evidence.result.atFills).toContain(AT_POSITIVE);
  expect(evidence.result.atFills).toContain(AT_NEGATIVE);
  expect(evidence.result.gcFills).toContain(evidence.paletteSkew.positive);
  expect(evidence.result.gcFills).toContain(evidence.paletteSkew.negative);
  expect(evidence.result.atFills).not.toEqual(evidence.result.gcFills);
  expect(evidence.result.atLegend).toEqual({
    positive: AT_POSITIVE,
    negative: AT_NEGATIVE
  });
};

const generateThroughUi = async (page) => {
  const requestCount = await page.evaluate(() => (
    window.__GBDRAW_CUSTOM_SKEW_REQUESTS__.length
  ));
  await page.getByRole('button', { name: 'Generate Diagram' }).click();
  await expect.poll(() => page.evaluate((previousCount) => ({
    processing: window.__GBDRAW_APP__.processing,
    requestCount: window.__GBDRAW_CUSTOM_SKEW_REQUESTS__.length,
    error: window.__GBDRAW_APP__.errorLog
  }), requestCount), { timeout: 300_000 }).toEqual({
    processing: false,
    requestCount: requestCount + 1,
    error: null
  });
  return page.evaluate(() => ({
    request: JSON.parse(JSON.stringify(window.__GBDRAW_CUSTOM_SKEW_REQUESTS__.at(-1))),
    svg: (() => {
      const normalizeColor = (value) => String(value || '').trim().toLowerCase();
      const inspectSvg = (svg) => {
        if (!svg) throw new Error('Expected an SVG Result.');
        const trackFills = (slotId) => {
          const group = svg.querySelector(`[data-gbdraw-slot-id="${CSS.escape(slotId)}"]`);
          if (!group) throw new Error(`Missing rendered track slot ${slotId}.`);
          return [...new Set([group, ...group.querySelectorAll('[fill]')]
            .map((element) => normalizeColor(element.getAttribute('fill')))
            .filter((value) => value && value !== 'none'))].sort();
        };
        const legendFill = (label) => {
          const text = [...svg.querySelectorAll('text')]
            .find((element) => String(element.textContent || '').trim() === label);
          if (!text) throw new Error(`Missing legend label ${label}.`);
          let group = text.parentElement;
          while (group && group !== svg) {
            const fill = [...group.querySelectorAll('[fill]')]
              .map((element) => normalizeColor(element.getAttribute('fill')))
              .find((value) => value && value !== 'none');
            if (fill) return fill;
            group = group.parentElement;
          }
          throw new Error(`Missing legend fill for ${label}.`);
        };
        return {
          atFills: trackFills('a_skew_2'),
          gcFills: trackFills('gc_skew'),
          atLegend: {
            positive: legendFill('AT skew (+)'),
            negative: legendFill('AT skew (-)')
          }
        };
      };
      const app = window.__GBDRAW_APP__;
      const selected = app.results?.[app.selectedResultIndex];
      const resultSvg = new DOMParser().parseFromString(
        String(selected?.content || ''),
        'image/svg+xml'
      ).documentElement;
      const mountedSvg = app.svgContainer?.querySelector?.('svg') || null;
      return {
        result: inspectSvg(resultSvg),
        mounted: inspectSvg(mountedSvg),
        paletteSkew: {
          positive: normalizeColor(app.currentColors?.skew_high),
          negative: normalizeColor(app.currentColors?.skew_low)
        }
      };
    })()
  }));
};

const changePaletteThroughUi = async (page) => {
  const colorsPanel = page.locator('details').filter({
    has: page.locator('summary[aria-label="Colors"]')
  });
  if (!await colorsPanel.getAttribute('open')) {
    await colorsPanel.locator('summary[aria-label="Colors"]').click();
  }
  const palette = page.getByLabel('Palette', { exact: true });
  await expect(palette).toBeVisible();
  const current = await palette.inputValue();
  const options = await palette.locator('option').evaluateAll((elements) => (
    elements.map((element) => element.value)
  ));
  const next = options.find((value) => value && value !== current);
  expect(next).toBeTruthy();
  await palette.selectOption(next);
  await expect(palette).toHaveValue(next);
  return next;
};

const saveSessionThroughUi = async (page) => {
  const downloadPromise = page.waitForEvent('download', { timeout: 180_000 });
  page.once('dialog', (dialog) => dialog.accept('custom-skew-round-trip'));
  await page.getByRole('button', { name: 'Save Session' }).click();
  const download = await downloadPromise;
  const path = await download.path();
  expect(path).toBeTruthy();
  const saved = JSON.parse(gunzipSync(readFileSync(path)).toString('utf8'));
  return { path, saved };
};

test('explicit AT-skew colors survive schema-5 Load, Generate, Save, fresh Load, and continued Generate', async ({
  browser,
  page
}) => {
  test.setTimeout(600_000);
  expect(sourceSession.renderRequest.schema).toBe(5);
  expect(requestEvidence(sourceSession.renderRequest).at.params).toEqual(expectedAtParams);

  const initialErrors = collectErrors(page);
  await installRequestObserver(page);
  await openApp(page);
  await loadSessionThroughUi(page, sourceSessionPath);
  await expectExplicitAtControls(page);
  expect(await activeAtParams(page)).toEqual(expectedAtParams);

  const firstGenerate = await generateThroughUi(page);
  expect(requestEvidence(firstGenerate.request)).toEqual({
    schema: 6,
    at: {
      id: AT_SLOT_ID,
      renderer: 'dinucleotide_skew',
      params: expectedAtParams
    },
    gc: {
      id: 'gc_skew',
      renderer: 'dinucleotide_skew',
      params: {}
    }
  });
  expectRenderedAtColors(firstGenerate.svg);

  await changePaletteThroughUi(page);
  expect(await activeAtParams(page)).toEqual(expectedAtParams);
  const paletteGenerate = await generateThroughUi(page);
  expect(requestEvidence(paletteGenerate.request).at.params).toEqual(expectedAtParams);
  expectRenderedAtColors(paletteGenerate.svg);

  const { path: savedPath, saved } = await saveSessionThroughUi(page);
  expect(saved.renderRequest.schema).toBe(6);
  expect(requestEvidence(saved.renderRequest).at.params).toEqual(expectedAtParams);
  expect(
    saved.config.adv.circular_track_slots.find((slot) => slot.id === AT_SLOT_ID).params
  ).toEqual(expectedAtParams);

  const freshContext = await browser.newContext();
  const freshPage = await freshContext.newPage();
  const freshErrors = collectErrors(freshPage);
  await installRequestObserver(freshPage);
  await openApp(freshPage);
  await loadSessionThroughUi(freshPage, savedPath);
  const freshAtGroup = await expectExplicitAtControls(freshPage);
  expect(await activeAtParams(freshPage)).toEqual(expectedAtParams);

  const freshGenerate = await generateThroughUi(freshPage);
  expect(requestEvidence(freshGenerate.request).at.params).toEqual(expectedAtParams);
  expectRenderedAtColors(freshGenerate.svg);

  const continuedGenerate = await generateThroughUi(freshPage);
  expect(requestEvidence(continuedGenerate.request).at.params).toEqual(expectedAtParams);
  expectRenderedAtColors(continuedGenerate.svg);

  await freshAtGroup.getByTitle('Inherit positive skew color').click();
  await expect(freshAtGroup.getByText('Inherited +', { exact: true })).toBeVisible();
  const inheritedPositive = await freshAtGroup
    .locator('input[type="color"][title="Positive skew color"]')
    .inputValue();
  expect(inheritedPositive).not.toBe(AT_POSITIVE);
  const clearedParams = await activeAtParams(freshPage);
  expect(clearedParams).toEqual({
    nt: 'AT',
    negative_color: AT_NEGATIVE,
    legend_label: 'AT skew'
  });
  const inheritedGenerate = await generateThroughUi(freshPage);
  expect(requestEvidence(inheritedGenerate.request).at.params).toEqual(clearedParams);
  expect(inheritedGenerate.svg.result.atFills).toContain(inheritedPositive);
  expect(inheritedGenerate.svg.result.atLegend.positive).toBe(inheritedPositive);

  expect(initialErrors).toEqual([]);
  expect(freshErrors).toEqual([]);
  await freshContext.close();
});
