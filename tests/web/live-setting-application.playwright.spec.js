// G1 live settings (R1 (c), PD-OI-102, PD-OI-103; PD-OI-066 scenario revision
// 2): a Legend or Plot Title change redraws the shown Result through the
// automatic rerender, which renders the committed Session plus the changed
// live values and commits them only when it succeeds and the committed
// Session has not changed since it started. Controls are selected by CSS.
const { test, expect } = require('@playwright/test');
const { readFileSync } = require('node:fs');
const { evaluateWithRetainedPromise, generateAndWaitForResult, settleLive } = require('./helpers/app-lifecycle.cjs');
const { load } = require('./helpers/mode-transition.cjs');
const { expectLiveEqualsGenerate, semanticSnapshot } = require('./helpers/live-generate-parity.cjs');

const TOBACCO = 'gbdraw/web/gallery/sessions/tobacco-chloroplast.gbdraw-session.json';
const HEPA = 'gbdraw/web/gallery/sessions/hepatoplasmataceae_orthogroup.gbdraw-session.json.gz';
const LAMBDA = 'gbdraw/web/gallery/sessions/lambda_basic_linear.gbdraw-session.json';
const LEGEND_FONT = 'input[aria-label="Legend font size"]';
const PLOT_TITLE = 'input[aria-label="Plot Title"]';
const LEGEND_FONT_OVERRIDE = 'objects.legends.font_size.short';
const DEFINITION_FONT_OVERRIDE = 'objects.definition.circular.font_size';
const LINEAR_DEFINITION_FONT_OVERRIDE = 'objects.definition.linear.font_size.short';
// An apply takes seconds; a base without live rows fails each case at this bound.
const APPLY_TIMEOUT = 120_000;

const generateNow = async (page) => {
  await generateAndWaitForResult(page);
  await settleLive(page);
};

// What the shown Result draws, the committed Session, and the draft.
const facts = (page) => page.evaluate(async () => {
  const app = window.__GBDRAW_APP__;
  const { state } = await import('./js/state.js');
  const { getCommittedCanonicalRenderRequest } = await import('./js/services/config.js');
  const request = getCommittedCanonicalRenderRequest();
  const svg = app.svgContainer?.querySelector('svg');
  const sizes = svg
    ? [...new Set([...svg.querySelectorAll('g[data-legend-key] text')]
        .map((text) => String(parseFloat(getComputedStyle(text).fontSize))))]
    : [];
  return {
    mode: app.mode,
    generationKey: state.resultGenerationKey.value,
    legendFontSizes: sizes,
    overrides: request?.diagramOptions?.configOverrides ?? null,
    plotTitle: request?.diagramOptions?.plotTitle ?? null,
    diagramOptions: JSON.stringify(request?.diagramOptions ?? null),
    content: app.results?.[app.selectedResultIndex]?.content ?? null,
    resultCount: Array.isArray(app.results) ? app.results.length : 0,
    draftLegendFont: app.adv.legend_font_size ?? null,
    reflowError: app.labelReflowLastError ? JSON.stringify(app.labelReflowLastError) : null,
    reflowing: Boolean(app.labelReflowProcessing)
  };
});

// Opens the collapsed cards that hold a control, then commits a value as the
// user does: type it and leave the field (its `change` ends the step).
const commitValue = async (page, selector, value) => {
  const input = page.locator(selector);
  for (const details of await input.locator('xpath=ancestor::details').all()) {
    if (await details.getAttribute('open') === null) await details.locator(':scope > summary').press('Enter');
  }
  await input.fill(String(value));
  await input.press('Tab');
};

const nextLegendFont = (sizes) => Math.round(Number(sizes[0] || 16)) + 5;

const appliedLegendFont = async (page, size) => {
  await expect.poll(async () => (await facts(page)).legendFontSizes, { timeout: APPLY_TIMEOUT }).toEqual([String(size)]);
  await settleLive(page);
};

// Holds the first diagram reply after it is installed until `releaseGate`,
// and counts every reply (the Worker's renders).
const installGate = (page) => page.evaluate(() => {
  const gate = { calls: 0, entered: false };
  gate.released = new Promise((resolve) => { gate.release = resolve; });
  window.__G1_GATE__ = gate;
  window.__GBDRAW_TEST_HOOKS__.beforeDiagramGenerationResponse = async () => {
    gate.calls += 1;
    if (gate.calls !== 1) return;
    gate.entered = true;
    await gate.released;
  };
});
const gateEntered = (page) => expect.poll(
  () => page.evaluate(() => window.__G1_GATE__.entered), { timeout: APPLY_TIMEOUT }
).toBe(true);
const releaseGate = (page) => page.evaluate(() => window.__G1_GATE__.release());
const gateCalls = (page) => page.evaluate(() => window.__G1_GATE__.calls);
const removeGate = (page) => page.evaluate(() => { delete window.__GBDRAW_TEST_HOOKS__.beforeDiagramGenerationResponse; });

for (const [name, file] of [['tobacco', TOBACCO], ['hepa', HEPA]]) {
  test(`${name}: a Legend font size change redraws the Result without Generate and equals a Generate from the same draft`, async ({ browser }) => {
    test.setTimeout(900_000);
    const page = await load(browser, file);
    try {
      await generateNow(page);
      const before = await facts(page);
      const size = nextLegendFont(before.legendFontSizes);
      await commitValue(page, LEGEND_FONT, size);
      await appliedLegendFont(page, size);
      const after = await facts(page);
      expect(after.generationKey, 'no Generate ran').toBe(before.generationKey);
      expect(after.reflowError).toBeNull();
      expect(after.overrides[LEGEND_FONT_OVERRIDE], 'the committed Session names the applied value').toBe(size);
      await expectLiveEqualsGenerate(page, { label: `${name} Legend font size` });
    } finally {
      await page.context().close();
    }
  });
}

test('lambda: a Plot Title change redraws the Linear Result and equals a Generate from the same draft', async ({ browser }) => {
  test.setTimeout(600_000);
  const page = await load(browser, LAMBDA);
  try {
    await generateNow(page);
    const before = await facts(page);
    await commitValue(page, PLOT_TITLE, 'G1 live title');
    await expect.poll(async () => (await facts(page)).content?.includes('G1 live title'), { timeout: APPLY_TIMEOUT }).toBe(true);
    await settleLive(page);
    const after = await facts(page);
    expect(after.generationKey).toBe(before.generationKey);
    expect(after.plotTitle).toBe('G1 live title');
    await expectLiveEqualsGenerate(page, { label: 'Linear Plot Title' });
  } finally {
    await page.context().close();
  }
});

// OV-1071: a select changed without a pointer or focus first (assistive or
// scripted input) records its History step, so its live value applies.
test('tobacco: a Legend position chosen without a pointer is one History step and redraws the Circular Result', async ({ browser }) => {
  test.setTimeout(600_000);
  const page = await load(browser, TOBACCO);
  try {
    await generateNow(page);
    const side = () => page.evaluate(() => JSON.parse(window.__GBDRAW_APP__.svgContainer.querySelector('svg')
      .getAttribute('data-gbdraw-composition')).legendSide);
    const steps = await page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount());
    const position = page.locator('#legend-position');
    for (const details of await position.locator('xpath=ancestor::details').all()) {
      if (await details.getAttribute('open') === null) await details.locator(':scope > summary').press('Enter');
    }
    const before = await side();
    const current = await position.inputValue();
    const next = (await position.locator('option').evaluateAll((options) => options.map((option) => option.value)))
      .find((value) => value !== current && value !== 'none');
    await position.selectOption(next);
    await expect.poll(async () => JSON.parse((await facts(page)).diagramOptions).output.legend, { timeout: APPLY_TIMEOUT })
      .toBe(next);
    await settleLive(page);
    expect(await side()).not.toBe(before);
    expect(await page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount())).toBe(steps + 1);
    await expectLiveEqualsGenerate(page, { label: 'Circular legend position' });
  } finally {
    await page.context().close();
  }
});

test('a live change is one History step; Undo and Redo return the Result after their apply', async ({ browser }) => {
  test.setTimeout(900_000);
  const page = await load(browser, TOBACCO);
  try {
    await generateNow(page);
    const original = await facts(page);
    const originalSnapshot = await semanticSnapshot(page);
    const steps = await page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount());
    const size = nextLegendFont(original.legendFontSizes);
    await commitValue(page, LEGEND_FONT, size);
    await appliedLegendFont(page, size);
    expect(await page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount())).toBe(steps + 1);
    const liveSnapshot = await semanticSnapshot(page);

    await evaluateWithRetainedPromise(page, () => window.__GBDRAW_HISTORY__.undo());
    await expect.poll(async () => (await facts(page)).legendFontSizes, { timeout: APPLY_TIMEOUT })
      .toEqual(original.legendFontSizes);
    await settleLive(page);
    const undone = await facts(page);
    expect(undone.diagramOptions).toBe(original.diagramOptions);
    expect(await semanticSnapshot(page)).toEqual(originalSnapshot);

    await evaluateWithRetainedPromise(page, () => window.__GBDRAW_HISTORY__.redo());
    await appliedLegendFont(page, size);
    expect((await facts(page)).overrides[LEGEND_FONT_OVERRIDE]).toBe(size);
    expect(await semanticSnapshot(page)).toEqual(liveSnapshot);
    expect(await page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount())).toBe(steps + 1);
  } finally {
    await page.context().close();
  }
});

test('a failed apply keeps the Result and the committed Session, says Live edit failed, and keeps the draft value', async ({ browser }) => {
  test.setTimeout(600_000);
  const page = await load(browser, TOBACCO);
  try {
    await generateNow(page);
    const before = await facts(page);
    await page.evaluate(() => {
      window.__GBDRAW_TEST_HOOKS__.beforeDiagramGenerationResponse = () => { throw new Error('G1 forced failure'); };
    });
    const size = nextLegendFont(before.legendFontSizes);
    await commitValue(page, LEGEND_FONT, size);
    await expect.poll(async () => (await facts(page)).reflowError, { timeout: APPLY_TIMEOUT }).not.toBeNull();
    await settleLive(page);
    await expect(page.locator('[data-live-application-feedback]')).toContainText('Live edit failed');
    const after = await facts(page);
    expect(after.content).toBe(before.content);
    expect(after.diagramOptions).toBe(before.diagramOptions);
    expect(after.draftLegendFont).toBe(size);
    await removeGate(page);
  } finally {
    await page.context().close();
  }
});

test('with an Applies on Generate value pending, a live change shows only the live value', async ({ browser }) => {
  test.setTimeout(900_000);
  const page = await load(browser, TOBACCO);
  try {
    await generateNow(page);
    const before = await facts(page);
    const definitionFont = before.overrides[DEFINITION_FONT_OVERRIDE];
    const draftDefinitionFont = await page.evaluate(() => {
      const value = window.__GBDRAW_APP__.adv.def_font_size;
      window.__GBDRAW_APP__.adv.def_font_size = 31;
      return value ?? null;
    });
    const size = nextLegendFont(before.legendFontSizes);
    await commitValue(page, LEGEND_FONT, size);
    await appliedLegendFont(page, size);
    const after = await facts(page);
    expect(after.overrides[DEFINITION_FONT_OVERRIDE], 'the pending value is not committed').toEqual(definitionFont);
    expect(after.overrides[LEGEND_FONT_OVERRIDE]).toBe(size);
    // The Result equals a Generate of the committed Session plus the live value.
    await page.evaluate((value) => { window.__GBDRAW_APP__.adv.def_font_size = value; }, draftDefinitionFont);
    await expectLiveEqualsGenerate(page, { label: 'a live change with a pending value' });
  } finally {
    await page.context().close();
  }
});

// 6.3: an Undo of a Generate step during a rerender makes it stale; it reruns
// from the restored committed Session.
test('an Undo of a Generate step during a rerender is not admitted, and the rerender runs from the restored Session', async ({ browser }) => {
  test.setTimeout(900_000);
  const page = await load(browser, TOBACCO);
  try {
    await generateNow(page);
    const first = await facts(page);
    const firstSnapshot = await semanticSnapshot(page);
    await page.evaluate(() => { window.__GBDRAW_APP__.adv.def_font_size = 31; });
    await generateNow(page);
    expect((await facts(page)).diagramOptions).not.toBe(first.diagramOptions);
    await installGate(page);
    await page.evaluate(async () => {
      const { state } = await import('./js/state.js');
      state.labelReflowForceRequestSeq.value += 1;
    });
    await gateEntered(page);
    await evaluateWithRetainedPromise(page, () => window.__GBDRAW_HISTORY__.undo());
    await releaseGate(page);
    await settleLive(page);
    await removeGate(page);
    const after = await facts(page);
    expect(after.reflowError).toBeNull();
    expect(after.diagramOptions).toBe(first.diagramOptions);
    expect(await semanticSnapshot(page)).toEqual(firstSnapshot);
  } finally {
    await page.context().close();
  }
});

test('a Generate during an apply wins: the apply is not admitted and no empty rerender follows', async ({ browser }) => {
  test.setTimeout(900_000);
  const page = await load(browser, TOBACCO);
  try {
    await generateNow(page);
    const before = await facts(page);
    await installGate(page);
    const size = nextLegendFont(before.legendFontSizes);
    await commitValue(page, LEGEND_FONT, size);
    await gateEntered(page);
    const generated = evaluateWithRetainedPromise(page, async () => (await window.__GBDRAW_APP__.runAnalysis())?.status);
    await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.processing)).toBe(true);
    await releaseGate(page);
    expect(await generated).toBe('ok');
    await settleLive(page);
    const after = await facts(page);
    expect(after.generationKey).toBe(before.generationKey + 1);
    expect(after.legendFontSizes).toEqual([String(size)]);
    expect(after.overrides[LEGEND_FONT_OVERRIDE]).toBe(size);
    expect(await gateCalls(page), 'the held apply and the Generate; no rerender after it').toBe(2);
    await removeGate(page);
  } finally {
    await page.context().close();
  }
});

test('Save waits for an apply, and Save then Load gives the same Result and request', async ({ browser }, testInfo) => {
  test.setTimeout(900_000);
  const page = await load(browser, TOBACCO);
  try {
    await generateNow(page);
    const before = await facts(page);
    await installGate(page);
    const size = nextLegendFont(before.legendFontSizes);
    await commitValue(page, LEGEND_FONT, size);
    await gateEntered(page);
    const save = page.locator('button[title^="Save the current Result"]');
    await expect(save).toBeDisabled();
    expect(await page.evaluate(async () => (await import('./js/state.js')).sessionOperationAvailability('save')?.reason))
      .toBe('Updating diagram. Retry after the update finishes.');
    await releaseGate(page);
    await appliedLegendFont(page, size);
    await removeGate(page);
    const applied = await facts(page);
    const file = testInfo.outputPath('g1-live.gbdraw-session.json');
    const pending = page.waitForEvent('download');
    await save.click();
    const large = page.locator('[data-choice-prompt] [data-choice-prompt-action="download"]');
    const first = await Promise.race([pending, large.waitFor({ state: 'visible' }).then(() => 'asked', () => null)]);
    if (first === 'asked') await large.click();
    await (first === 'asked' ? await pending : first).saveAs(file);
    expect(readFileSync(file).length).toBeGreaterThan(0);
    const loaded = await load(browser, file);
    try {
      const reloaded = await facts(loaded);
      expect(reloaded.content).toBe(applied.content);
      expect(reloaded.diagramOptions).toBe(applied.diagramOptions);
    } finally {
      await loaded.context().close();
    }
  } finally {
    await page.context().close();
  }
});

// Linear: in Circular, Reset turns Multi-Record Canvas on (its default), and a
// pending Multi-Record Canvas change holds every automatic rerender
// (`shouldDeferCircularPreviewUpdates`, OV-1070).
test('Reset Settings applies the live rows and leaves the rest pending', async ({ browser }) => {
  test.setTimeout(900_000);
  const page = await load(browser, LAMBDA);
  try {
    await generateNow(page);
    const size = nextLegendFont((await facts(page)).legendFontSizes);
    await commitValue(page, LEGEND_FONT, size);
    await appliedLegendFont(page, size);
    await page.evaluate(() => { window.__GBDRAW_APP__.adv.def_font_size = 31; });
    const before = await facts(page);
    await page.locator('button[aria-label="Reset Settings"]').click();
    await page.locator('[data-choice-prompt] [data-choice-prompt-action="reset"]').click();
    await expect.poll(async () => (await facts(page)).legendFontSizes, { timeout: APPLY_TIMEOUT }).not.toEqual([String(size)]);
    await settleLive(page);
    const after = await facts(page);
    expect(after.generationKey).toBe(before.generationKey);
    expect(after.draftLegendFont).toBeNull();
    expect(after.overrides[LEGEND_FONT_OVERRIDE], 'the default leaves the override out').toBeUndefined();
    expect(after.overrides[LINEAR_DEFINITION_FONT_OVERRIDE], 'an Applies on Generate value stays pending')
      .toEqual(before.overrides[LINEAR_DEFINITION_FONT_OVERRIDE]);
  } finally {
    await page.context().close();
  }
});

test('a live change in a mode without a Result applies nothing until that mode has one', async ({ browser }) => {
  test.setTimeout(900_000);
  const page = await load(browser, TOBACCO);
  try {
    await generateNow(page);
    const circular = await facts(page);
    await page.locator('button[aria-label="Linear"][aria-pressed]').click();
    await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.mode)).toBe('linear');
    await settleLive(page);
    const linear = await facts(page);
    expect(linear.resultCount).toBe(0);
    await commitValue(page, LEGEND_FONT, 23);
    await settleLive(page);
    const after = await facts(page);
    expect(after.reflowing).toBe(false);
    expect(after.resultCount).toBe(0);
    expect(after.generationKey).toBe(linear.generationKey);
    await page.locator('button[aria-label="Circular"][aria-pressed]').click();
    await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.mode)).toBe('circular');
    await settleLive(page);
    const back = await facts(page);
    expect(back.content).toBe(circular.content);
    expect(back.diagramOptions).toBe(circular.diagramOptions);
    // The shown mode's Result applies its own change.
    const size = nextLegendFont(back.legendFontSizes);
    await commitValue(page, LEGEND_FONT, size);
    await appliedLegendFont(page, size);
  } finally {
    await page.context().close();
  }
});
