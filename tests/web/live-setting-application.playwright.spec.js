// G1 live settings (R1 (c), PD-OI-102, PD-OI-103; PD-OI-066 scenario revision
// 2): a Legend or Plot Title change redraws the shown Result through the
// automatic rerender, which renders the committed Session plus the changed
// live values and commits them only when it succeeds and the committed
// Session has not changed since it started. Controls are selected by CSS.
const { test, expect } = require('@playwright/test');
const { readFileSync, writeFileSync } = require('node:fs');
const { evaluateWithRetainedPromise, generateAndWaitForResult, reveal, settleLive } = require('./helpers/app-lifecycle.cjs');
const { load } = require('./helpers/mode-transition.cjs');
const { expectLiveEqualsGenerate, semanticSnapshot } = require('./helpers/live-generate-parity.cjs');
const { open, addColorRule, addVisibilityRule } = require('./helpers/live-generate-parity-steps.cjs');

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

// Holds the diagram replies `held` (1-based, counted from its install) until
// `releaseGate`, and counts every reply (the Worker's renders).
const installGate = (page, held = [1]) => page.evaluate((held) => {
  const gate = { calls: 0, entered: [], releases: {} };
  const released = Object.fromEntries(held.map((call) => [call, new Promise((resolve) => { gate.releases[call] = resolve; })]));
  window.__G1_GATE__ = gate;
  window.__GBDRAW_TEST_HOOKS__.beforeDiagramGenerationResponse = async () => {
    gate.calls += 1;
    const call = gate.calls;
    if (!released[call]) return;
    gate.entered.push(call);
    await released[call];
  };
}, held);
const gateEntered = (page, call = 1) => expect.poll(
  () => page.evaluate((call) => window.__G1_GATE__.entered.includes(call), call), { timeout: APPLY_TIMEOUT }
).toBe(true);
const releaseGate = (page, call = 1) => page.evaluate((call) => window.__G1_GATE__.releases[call](), call);
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

// Review F1: a Session whose request has no `configOverrides`; a live value creates them, as a Generate from the same draft
// sends them. The Gallery Session holds its command line's overrides, so the
// case removes them.
test('tobacco: a loaded Session applies a Legend font size with no Generate', async ({ browser }, testInfo) => {
  test.setTimeout(600_000);
  const session = JSON.parse(readFileSync(TOBACCO, 'utf8'));
  delete session.renderRequest.diagramOptions.configOverrides;
  const file = testInfo.outputPath('no-overrides.gbdraw-session.json');
  writeFileSync(file, JSON.stringify(session));
  const page = await load(browser, file);
  try {
    await settleLive(page);
    const before = await facts(page);
    expect(before.overrides, 'the Gallery request has no config overrides').toBeNull();
    const size = nextLegendFont(before.legendFontSizes);
    await commitValue(page, LEGEND_FONT, size);
    await appliedLegendFont(page, size);
    const after = await facts(page);
    expect(after.reflowError).toBeNull();
    expect(after.generationKey).toBe(before.generationKey);
    expect(after.overrides[LEGEND_FONT_OVERRIDE]).toBe(size);
  } finally {
    await page.context().close();
  }
});

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

// OV-1116: Enter commits a value and the field keeps focus. Each value committed
// with Enter in the still-focused field is one History step and one apply, and
// leaving the field with the committed value adds no step.
test('tobacco: each value committed with Enter in a focused field is one step and redraws the Result', async ({ browser }) => {
  test.setTimeout(900_000);
  const page = await load(browser, TOBACCO);
  try {
    await generateNow(page);
    const undoCount = () => page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount());
    // The Gallery Session draws no plot title; draw it so its edits show.
    const position = page.locator('select[aria-label="Plot Title Position"]');
    for (const details of await position.locator('xpath=ancestor::details').all()) {
      if (await details.getAttribute('open') === null) await details.locator(':scope > summary').press('Enter');
    }
    await position.selectOption('top');
    await expect.poll(async () => JSON.parse((await facts(page)).diagramOptions).output?.plotTitlePosition ?? null, {
      timeout: APPLY_TIMEOUT
    }).toBe('top');
    await settleLive(page);
    for (const [selector, values] of [
      [PLOT_TITLE, ['Enter A', 'Enter B', 'Enter C']],
      ['input[aria-label="Plot Title Font Size"]', ['40', '24', '36']],
      [LEGEND_FONT, ['20', '28', '16']]
    ]) {
      const input = page.locator(selector);
      for (const details of await input.locator('xpath=ancestor::details').all()) {
        if (await details.getAttribute('open') === null) await details.locator(':scope > summary').press('Enter');
      }
      await input.focus();
      for (const value of values) {
        const before = await facts(page);
        const steps = await undoCount();
        await input.fill(value);
        await input.press('Enter');
        await expect.poll(async () => (await facts(page)).diagramOptions, {
          timeout: APPLY_TIMEOUT, message: `${selector} = ${value} applies`
        }).not.toBe(before.diagramOptions);
        await settleLive(page);
        expect((await facts(page)).content, `${selector} = ${value} redraws`).not.toBe(before.content);
        expect(await undoCount(), `${selector} = ${value} is one step`).toBe(steps + 1);
      }
      const steps = await undoCount();
      await input.press('Tab');
      await settleLive(page);
      expect(await undoCount(), `leaving ${selector} adds no step`).toBe(steps);
    }
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

// Review F7 (d): a failed value waits for its next change or Generate; a
// later rerender does not apply it (R1 (c)).
test('a failed apply does not run again', async ({ browser }) => {
  test.setTimeout(600_000);
  const page = await load(browser, TOBACCO);
  try {
    await generateNow(page);
    const before = await facts(page);
    await page.evaluate(() => {
      window.__G1_FAILURES__ = 0;
      window.__GBDRAW_TEST_HOOKS__.beforeDiagramGenerationResponse = () => {
        window.__G1_FAILURES__ += 1;
        throw new Error('G1 forced failure');
      };
    });
    const size = nextLegendFont(before.legendFontSizes);
    await commitValue(page, LEGEND_FONT, size);
    await expect.poll(async () => (await facts(page)).reflowError, { timeout: APPLY_TIMEOUT }).not.toBeNull();
    await settleLive(page);
    expect(await page.evaluate(() => window.__G1_FAILURES__)).toBe(1);
    await removeGate(page);
    // An editor rerender after it draws the committed Session without the failed value.
    await page.evaluate(async () => {
      const { state } = await import('./js/state.js');
      state.labelReflowForceRequestSeq.value += 1;
    });
    await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.labelReflowProcessing)).toBe(true);
    await settleLive(page);
    const after = await facts(page);
    expect(after.reflowError).toBeNull();
    expect(after.legendFontSizes).toEqual(before.legendFontSizes);
    expect(after.diagramOptions).toBe(before.diagramOptions);
    expect(after.draftLegendFont).toBe(size);
  } finally {
    await page.context().close();
  }
});

// Review F6: the label binder of a live apply reports to it after the apply
// committed its value (the committed revision it checks is the new one).
test('a forced label the live apply does not draw says Live edit failed', async ({ browser }) => {
  test.setTimeout(600_000);
  const page = await load(browser, TOBACCO);
  try {
    await generateNow(page);
    // Label visibility On for a labelled feature: the binder checks its label.
    const forced = await evaluateWithRetainedPromise(page, async () => {
      const app = window.__GBDRAW_APP__;
      const label = app.svgContainer.querySelector('svg text[data-label-feature-id]');
      const featureId = label.getAttribute('data-label-feature-id');
      await app.openFeatureEditorFromList(app.extractedFeatures.find((item) => item.svg_id === featureId), null);
      await window.Vue.nextTick();
      app.clickedFeature.labelVisibility = 'on';
      await app.updateClickedFeatureLabelText();
      app.clickedFeature = null;
      return label.textContent;
    });
    await settleLive(page);
    expect((await facts(page)).reflowError).toBeNull();
    // The apply's Result is mounted without that label. Each Result gets a
    // new preview surface, whose markup v-preview-markup sets through
    // `replaceChildren`; the binder reads it in the same flush.
    await page.evaluate((text) => {
      const replaceChildren = HTMLElement.prototype.replaceChildren;
      window.__G1_RESTORE_MOUNT__ = () => { HTMLElement.prototype.replaceChildren = replaceChildren; };
      HTMLElement.prototype.replaceChildren = function replaceWithoutLabel(...nodes) {
        if (this.classList.contains('gbdraw-preview-surface')) {
          for (const node of nodes) {
            node.querySelectorAll?.('text').forEach((label) => { if (label.textContent === text) label.remove(); });
          }
        }
        return replaceChildren.apply(this, nodes);
      };
    }, forced);
    const size = nextLegendFont((await facts(page)).legendFontSizes);
    await commitValue(page, LEGEND_FONT, size);
    await expect.poll(async () => (await facts(page)).reflowError, { timeout: APPLY_TIMEOUT }).not.toBeNull();
    await page.evaluate(() => window.__G1_RESTORE_MOUNT__());
    await settleLive(page);
    await expect(page.locator('[data-live-application-feedback]')).toContainText('Live edit failed');
    expect((await facts(page)).overrides[LEGEND_FONT_OVERRIDE], 'the apply itself committed').toBe(size);
  } finally {
    await page.context().close();
  }
});

// Review F3: a checkpoint step (the Bakta preset changes the Legend rows)
// applies the Legend sizes it sets in its own rerender; Undo and Redo give
// the Result of the action.
test('the Bakta preset applies its Legend sizes; Undo and Redo give the same Results', async ({ browser }) => {
  test.setTimeout(600_000);
  const page = await load(browser, TOBACCO);
  try {
    await generateNow(page);
    const original = await facts(page);
    const originalSnapshot = await semanticSnapshot(page);
    expect(original.legendFontSizes).not.toEqual(['12']);
    const preset = page.locator('select[aria-label="Preset scheme"]');
    for (const details of await preset.locator('xpath=ancestor::details').all()) {
      if (await details.getAttribute('open') === null) await details.locator(':scope > summary').press('Enter');
    }
    await preset.selectOption('bakta');
    await evaluateWithRetainedPromise(page, () => window.__GBDRAW_APP__.applySpecificRulePreset());
    await appliedLegendFont(page, 12);
    const applied = await facts(page);
    expect(applied.generationKey).toBe(original.generationKey);
    expect(applied.reflowError).toBeNull();
    expect(applied.overrides[LEGEND_FONT_OVERRIDE]).toBe(12);
    const appliedSnapshot = await semanticSnapshot(page);

    await evaluateWithRetainedPromise(page, () => window.__GBDRAW_HISTORY__.undo());
    await expect.poll(async () => (await facts(page)).legendFontSizes, { timeout: APPLY_TIMEOUT })
      .toEqual(original.legendFontSizes);
    await settleLive(page);
    expect(await semanticSnapshot(page)).toEqual(originalSnapshot);

    await evaluateWithRetainedPromise(page, () => window.__GBDRAW_HISTORY__.redo());
    await appliedLegendFont(page, 12);
    expect(await semanticSnapshot(page)).toEqual(appliedSnapshot);
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

// 6.3, review F7 (a): a value committed while a Generate runs applies after
// it (the replay). An Undo of that Generate step while the apply is held
// makes it stale by the committed revision alone: nothing of it is committed
// into the restored Session, and it runs again from that Session.
test('an Undo of a Generate step during an apply commits nothing from it; it runs again from the restored Session', async ({ browser }) => {
  test.setTimeout(900_000);
  const page = await load(browser, TOBACCO);
  try {
    await generateNow(page);
    const first = await facts(page);
    const firstSnapshot = await semanticSnapshot(page);
    await page.evaluate(() => { window.__GBDRAW_APP__.adv.def_font_size = 31; });
    await installGate(page, [1, 2]);
    const generated = evaluateWithRetainedPromise(page, async () => (await window.__GBDRAW_APP__.runAnalysis())?.status);
    await gateEntered(page, 1);
    const size = nextLegendFont(first.legendFontSizes);
    await commitValue(page, LEGEND_FONT, size);
    await releaseGate(page, 1);
    expect(await generated).toBe('ok');
    // The replay applies the value the Generate did not draw; hold it.
    await gateEntered(page, 2);
    expect(await page.evaluate(() => window.__GBDRAW_HISTORY__.undoLabel())).toMatch(/Generate/i);
    await evaluateWithRetainedPromise(page, () => window.__GBDRAW_HISTORY__.undo());
    const restored = await facts(page);
    expect(restored.diagramOptions).toBe(first.diagramOptions);
    await releaseGate(page, 2);
    await settleLive(page);
    await removeGate(page);
    const after = await facts(page);
    expect(after.reflowError).toBeNull();
    const drawnFont = after.draftLegendFont === null ? undefined : Number(after.draftLegendFont);
    // The Result shows the restored Session plus the restored draft's live value.
    expect(after.overrides?.[LEGEND_FONT_OVERRIDE]).toBe(drawnFont);
    expect(after.overrides?.[DEFINITION_FONT_OVERRIDE]).toEqual(first.overrides?.[DEFINITION_FONT_OVERRIDE]);
    if (drawnFont === undefined) {
      expect(after.diagramOptions).toBe(first.diagramOptions);
      expect(await semanticSnapshot(page)).toEqual(firstSnapshot);
    } else {
      expect(await gateCalls(page), 'the stale apply and its rerun').toBe(3);
      expect(after.legendFontSizes).toEqual([String(drawnFont)]);
    }
  } finally {
    await page.context().close();
  }
});

// Review F7 (b): values committed while a Generate runs. One changed and
// changed back is equal to the Generate's: the replay renders nothing. One
// the Generate did not draw is applied by the replay.
test('a value committed during a Generate: the replay skips an equal value and applies a new one', async ({ browser }) => {
  test.setTimeout(900_000);
  const page = await load(browser, TOBACCO);
  try {
    await generateNow(page);
    const before = await facts(page);
    const size = nextLegendFont(before.legendFontSizes);
    const draft = await page.locator(LEGEND_FONT).inputValue().catch(() => '');
    await installGate(page);
    let generated = evaluateWithRetainedPromise(page, async () => (await window.__GBDRAW_APP__.runAnalysis())?.status);
    await gateEntered(page);
    await commitValue(page, LEGEND_FONT, size);
    await commitValue(page, LEGEND_FONT, draft);
    await releaseGate(page);
    expect(await generated).toBe('ok');
    await settleLive(page);
    expect(await gateCalls(page), 'the Generate alone').toBe(1);
    expect((await facts(page)).legendFontSizes).toEqual(before.legendFontSizes);

    await installGate(page);
    generated = evaluateWithRetainedPromise(page, async () => (await window.__GBDRAW_APP__.runAnalysis())?.status);
    await gateEntered(page);
    await commitValue(page, LEGEND_FONT, size);
    await releaseGate(page);
    expect(await generated).toBe('ok');
    await appliedLegendFont(page, size);
    expect(await gateCalls(page), 'the Generate and the replay').toBe(2);
    const after = await facts(page);
    expect(after.reflowError).toBeNull();
    expect(after.overrides[LEGEND_FONT_OVERRIDE]).toBe(size);
    await removeGate(page);
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

// OV-1070 (review F2): in Circular, Reset turns Multi-Record Canvas on (its
// default) while the Result is single-canvas; that pending change no longer
// holds the automatic rerender.
for (const [name, file, definitionOverride] of [
  ['lambda', LAMBDA, LINEAR_DEFINITION_FONT_OVERRIDE], ['tobacco', TOBACCO, DEFINITION_FONT_OVERRIDE]
]) {
test(`${name}: Reset Settings applies the live rows and leaves the rest pending`, async ({ browser }) => {
  test.setTimeout(900_000);
  const page = await load(browser, file);
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
    expect(after.overrides[definitionOverride], 'an Applies on Generate value stays pending')
      .toEqual(before.overrides[definitionOverride]);
  } finally {
    await page.context().close();
  }
});
}

// OV-1074: Reset Settings turns Multi-Record Canvas on and fills the
// multi-record Legend and Plot Title positions. A single-canvas Result keeps
// its own canvas's positions. Reset holds the rerender over its checkpoint
// (OV-1076), so the removed color rule and the live values draw in one render.
test('Reset Settings on a batch Circular Result keeps its positions and draws once', async ({ browser }) => {
  test.setTimeout(600_000);
  const page = await (await browser.newContext({ viewport: { width: 1600, height: 1000 } })).newPage();
  try {
    await open(page, { mode: 'circular', results: 'batch', reflow: 'off' });
    await addColorRule(page, { feat: 'CDS', qual: 'locus_tag', val: '_0001$', color: '#c83366', cap: 'CDS' });
    await generateNow(page);
    const before = JSON.parse((await facts(page)).diagramOptions).output;
    expect(before).toEqual({ legend: 'left', plotTitlePosition: 'none' });
    await installGate(page, []);
    await page.locator('button[aria-label="Reset Settings"]').click();
    await page.locator('[data-choice-prompt] [data-choice-prompt-action="reset"]').click();
    await settleLive(page);
    const after = await facts(page);
    expect(await page.evaluate(() => window.__GBDRAW_APP__.form.multi_record_canvas), 'the draft canvas is multi-record').toBe(true);
    expect(JSON.parse(after.diagramOptions).output, 'the committed single-canvas positions').toEqual(before);
    expect(await gateCalls(page), 'one rerender for the removed rule').toBe(1);
    expect(after.reflowError).toBeNull();
  } finally {
    await removeGate(page);
    await page.context().close();
  }
});

// OV-1075: a live request for a value the running editor rerender does not
// apply, equal to the committed one, does not make that rerender stale: it
// commits, and the next iteration finds nothing to write.
test('a live request during an editor rerender for a value equal to the committed one draws once', async ({ browser }) => {
  test.setTimeout(600_000);
  const page = await load(browser, TOBACCO);
  try {
    await generateNow(page);
    const before = await facts(page);
    const size = nextLegendFont(before.legendFontSizes);
    const draft = await page.locator(LEGEND_FONT).inputValue().catch(() => '');
    await installGate(page);
    await page.evaluate(async () => {
      const { state } = await import('./js/state.js');
      state.labelReflowForceRequestSeq.value += 1;
    });
    await gateEntered(page);
    await commitValue(page, LEGEND_FONT, size);
    await commitValue(page, LEGEND_FONT, draft);
    await releaseGate(page);
    await settleLive(page);
    const after = await facts(page);
    expect(await gateCalls(page), 'the editor rerender alone').toBe(1);
    expect(after.generationKey).toBe(before.generationKey);
    expect(after.legendFontSizes).toEqual(before.legendFontSizes);
    expect(after.diagramOptions).toBe(before.diagramOptions);
    expect(after.reflowError).toBeNull();
  } finally {
    await removeGate(page);
    await page.context().close();
  }
});

// OV-1076: one Reset Settings asks for the editor rerender twice: its label
// projection for a feature the Result does not draw, then the rule owner for
// the removed color rule, with the same editor state. One render draws both.
for (const [mode, results, rule, hidden] of [
  ['circular', 'batch', '_0001$', '_0004$'], ['linear', 'single', 'FL1$', 'FL2$']
]) {
test(`${mode}: Reset Settings that shows a hidden feature and removes a color rule draws once`, async ({ browser }) => {
  test.setTimeout(600_000);
  const page = await (await browser.newContext({ viewport: { width: 1600, height: 1000 } })).newPage();
  try {
    await open(page, { mode, results, reflow: 'off' });
    await addColorRule(page, { feat: 'CDS', qual: 'locus_tag', val: rule, color: '#c83366', cap: 'CDS' });
    await addVisibilityRule(page, { recordId: '*', featureType: 'CDS', qualifier: 'locus_tag', value: hidden, action: 'off' });
    await generateNow(page);
    await installGate(page, []);
    await page.locator('button[aria-label="Reset Settings"]').click();
    await page.locator('[data-choice-prompt] [data-choice-prompt-action="reset"]').click();
    await settleLive(page);
    expect(await gateCalls(page), 'one rerender for the shown feature and the removed rule').toBe(1);
    expect((await facts(page)).reflowError).toBeNull();
  } finally {
    await removeGate(page);
    await page.context().close();
  }
});
}

// Review F2: the layout rows apply live only to a multi-record canvas Result.
// OV-1077: during that Multi-Record Canvas change a Legend or Plot Title
// position applies on Generate; their sizes and the title text stay live.
test('with a single-canvas Result, the Multi-Record layout rows say Applies on Generate', async ({ browser }) => {
  test.setTimeout(600_000);
  const page = await load(browser, TOBACCO);
  try {
    await generateNow(page);
    const positionNote = /the position applies on Generate after the Multi-Record Canvas change/;
    for (const id of ['legend', 'plot-title']) await expect(page.locator(`[data-setting="${id}"]`).getByText(positionNote)).toHaveCount(0);
    await page.evaluate(() => { window.__GBDRAW_APP__.form.multi_record_canvas = true; });
    const note = page.locator('[data-setting="multi-record-layout"] > span').first();
    await expect(note).toHaveText(/Applies on Generate/);
    await expect(note).toHaveAttribute('data-applies-on-generate', 'multi-record-layout');
    for (const id of ['legend', 'plot-title']) await expect(page.locator(`[data-setting="${id}"]`).getByText(positionNote)).toHaveCount(1);
  } finally {
    await page.context().close();
  }
});

// Review F4 (design 6.7): Load applies nothing. A value an apply left pending
// (its render was canceled after the reply) does not reach a loaded Session
// whose draft holds another value of that setting.
test('a Load after a canceled apply applies none of the loaded draft', async ({ browser }, testInfo) => {
  test.setTimeout(900_000);
  const page = await load(browser, TOBACCO);
  try {
    const gallery = await facts(page);
    await generateNow(page);
    const generated = await facts(page);
    const size = nextLegendFont(generated.legendFontSizes);
    await page.evaluate(() => {
      window.__GBDRAW_TEST_HOOKS__.beforeDiagramGenerationResponse = async () => {
        delete window.__GBDRAW_TEST_HOOKS__.beforeDiagramGenerationResponse;
        (await import('./js/services/diagram-generation.js')).cancelDiagramGeneration();
      };
    });
    await commitValue(page, LEGEND_FONT, size);
    await expect.poll(() => page.evaluate(() => window.__GBDRAW_TEST_HOOKS__.beforeDiagramGenerationResponse === undefined),
      { timeout: APPLY_TIMEOUT }).toBe(true);
    await settleLive(page);
    const canceled = await facts(page);
    expect(canceled.reflowError).toBeNull();
    expect(canceled.diagramOptions, 'the canceled apply committed nothing').toBe(generated.diagramOptions);
    // The Gallery Session with a draft Legend font size its Result does not
    // draw (and the config overrides object a Generate writes).
    const session = JSON.parse(readFileSync(TOBACCO, 'utf8'));
    session.modes.circular.config.adv.legend_font_size = size + 7;
    session.renderRequest.diagramOptions.configOverrides = {};
    const file = testInfo.outputPath('pending-legend.gbdraw-session.json');
    writeFileSync(file, JSON.stringify(session));
    await page.locator('input[accept^=".json,"]').setInputFiles(file);
    await expect.poll(() => page.evaluate(() => !window.__GBDRAW_APP__.sessionImportPending
      && window.__GBDRAW_APP__.extractedFeatures.length > 0), { timeout: 180_000 }).toBe(true);
    await settleLive(page);
    const loaded = await facts(page);
    expect(Number(loaded.draftLegendFont)).toBe(size + 7);
    expect(loaded.legendFontSizes).toEqual(gallery.legendFontSizes);
    // An editor rerender draws the loaded committed Session.
    await page.evaluate(async () => {
      const { state } = await import('./js/state.js');
      state.labelReflowForceRequestSeq.value += 1;
    });
    await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.labelReflowProcessing)).toBe(true);
    await settleLive(page);
    const after = await facts(page);
    expect(after.reflowError).toBeNull();
    expect(after.diagramOptions).toBe(loaded.diagramOptions);
    expect(after.overrides[LEGEND_FONT_OVERRIDE]).toBeUndefined();
    expect(after.legendFontSizes).toEqual(gallery.legendFontSizes);
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

// Opens the collapsed cards that hold a control, then clicks it.
const clickControl = async (page, selector) => {
  const control = page.locator(selector);
  for (const details of await control.locator('xpath=ancestor::details').all()) {
    if (await details.getAttribute('open') === null) await details.locator(':scope > summary').press('Enter');
  }
  await control.click();
};

// PD-OI-103 G1: the Linear Lock Definition Column redraws the shown Result.
test('lambda: Lock Definition Column redraws the Linear Result without Generate and equals a Generate from the same draft', async ({ browser }) => {
  test.setTimeout(600_000);
  const page = await load(browser, LAMBDA);
  const lockKey = 'canvas.linear.keep_definition_left_aligned';
  try {
    await generateNow(page);
    const before = await facts(page);
    const steps = await page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount());
    await clickControl(page, 'input[type="checkbox"][aria-describedby="linear-definition-lock-help"]');
    await expect.poll(async () => (await facts(page)).overrides?.[lockKey], { timeout: APPLY_TIMEOUT })
      .toBe(!before.overrides[lockKey]);
    await settleLive(page);
    const after = await facts(page);
    expect(after.generationKey, 'no Generate ran').toBe(before.generationKey);
    expect(after.reflowError).toBeNull();
    expect(after.content).not.toBe(before.content);
    expect(await page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount())).toBe(steps + 1);
    await expect(page.locator('[data-setting="definition-column"] > span.apply-note')).toHaveText(/Live edit/);
    await expectLiveEqualsGenerate(page, { label: 'Linear Lock Definition Column' });
  } finally {
    await page.context().close();
  }
});

// OV-1078 (D-80): Reset Settings turns Multi-Record Canvas on, so the reset
// single-canvas positions wait; turning it off again draws the shown canvas,
// and they apply without Generate.
test('tobacco: Reset Settings, then Multi-Record Canvas off, applies the reset positions without Generate', async ({ browser }) => {
  test.setTimeout(600_000);
  const page = await load(browser, TOBACCO);
  const output = async () => JSON.parse((await facts(page)).diagramOptions).output;
  try {
    await generateNow(page);
    const before = await facts(page);
    const shown = await output();
    await page.locator('button[aria-label="Reset Settings"]').click();
    await page.locator('[data-choice-prompt] [data-choice-prompt-action="reset"]').click();
    await settleLive(page);
    expect(await page.evaluate(() => window.__GBDRAW_APP__.form.multi_record_canvas)).toBe(true);
    expect(await output(), 'the positions wait while the draft draws the multi-record canvas').toEqual(shown);
    const reset = await facts(page);
    await clickControl(page, 'input[type="checkbox"][data-circular-canvas-setting]');
    const draft = await page.evaluate(() => ({
      legend: window.__GBDRAW_APP__.form.legend, plotTitlePosition: window.__GBDRAW_APP__.adv.plot_title_position
    }));
    expect(draft, 'Reset changed the single-canvas positions').not.toEqual(shown);
    await expect.poll(output, { timeout: APPLY_TIMEOUT }).toEqual(draft);
    await settleLive(page);
    const after = await facts(page);
    expect(after.generationKey, 'no Generate ran').toBe(before.generationKey);
    expect(after.reflowError).toBeNull();
    expect(after.content, 'the Result draws them').not.toBe(reset.content);
  } finally {
    await page.context().close();
  }
});

// S4-Q5 (design 6.6): an editor edit made while a live setting apply runs
// waits for that apply, then applies; the Result then equals a Generate from
// the same draft.
const COLOR_RULE_PATTERN = 'input[aria-label="Color rule 1 pattern"]';
const VISIBILITY_RULE = (field) => `[aria-label="Visibility rule 1 ${field}"]`;

// OV-1133: the pattern edit waited for the apply and was then dropped as
// stale, with its text left in the field as "Not applied".
test('tobacco: a color rule pattern committed during a Legend font size apply is drawn and equals a Generate (OV-1133)', async ({ browser }) => {
  test.setTimeout(900_000);
  const page = await load(browser, TOBACCO);
  try {
    await generateNow(page);
    const pattern = page.locator(COLOR_RULE_PATTERN);
    await reveal(pattern);
    const before = await facts(page);
    const value = `${await pattern.inputValue()}|psbA`;
    const size = nextLegendFont(before.legendFontSizes);
    await installGate(page);
    await commitValue(page, LEGEND_FONT, size);
    await gateEntered(page);
    await pattern.fill(value);
    await pattern.press('Tab');
    await releaseGate(page);
    await appliedLegendFont(page, size);
    const rule = await page.evaluate(() => {
      const app = window.__GBDRAW_APP__;
      const row = app.manualSpecificRules[0];
      return { val: row.val, draft: app.specificRulePatternDraft(row) };
    });
    expect(rule, 'the edit is committed and no draft stays').toEqual({ val: value, draft: null });
    await expect(pattern).toHaveValue(value);
    expect((await facts(page)).reflowError).toBeNull();
    await expectLiveEqualsGenerate(page, { label: 'color rule pattern during a Legend font size apply' });
  } finally {
    await removeGate(page);
    await page.context().close();
  }
});

// OV-1134: the edit shows psbA, which the apply's Result does not draw (the
// apply started with the old rule). The display of that Result superseded
// the edit's projection, and the rerender that must draw psbA never ran.
test('tobacco: a visibility rule edit during a Legend font size apply equals a Generate (OV-1134)', async ({ browser }) => {
  test.setTimeout(900_000);
  const page = await load(browser, TOBACCO);
  try {
    // The rule controls are enabled once a Result is shown.
    await generateNow(page);
    const add = page.locator('button[aria-label="Add feature visibility rule"]');
    await reveal(add);
    await add.click();
    await reveal(page.locator(VISIBILITY_RULE('value regex')));
    await page.locator(VISIBILITY_RULE('feature type')).selectOption('CDS');
    await commitValue(page, VISIBILITY_RULE('qualifier'), 'gene');
    await commitValue(page, VISIBILITY_RULE('value regex'), '^psbA$');
    await page.locator(VISIBILITY_RULE('action')).selectOption('off');
    await settleLive(page);
    const size = nextLegendFont((await facts(page)).legendFontSizes);
    await commitValue(page, LEGEND_FONT, size);
    // The edit as the field commits it, at once once the apply runs.
    const editedDuringApply = await evaluateWithRetainedPromise(page, async (selector) => {
      const app = window.__GBDRAW_APP__;
      const deadline = Date.now() + 60_000;
      while (!app.labelReflowProcessing && Date.now() < deadline) await new Promise((resolve) => requestAnimationFrame(resolve));
      const input = document.querySelector(selector);
      input.value = '^psbB$';
      input.dispatchEvent(new Event('change', { bubbles: true }));
      return app.labelReflowProcessing;
    }, VISIBILITY_RULE('value regex'));
    expect(editedDuringApply, 'the edit is made while the apply runs').toBe(true);
    await appliedLegendFont(page, size);
    expect(await page.evaluate(() => window.__GBDRAW_APP__.featureVisibilityManualRules[0].value)).toBe('^psbB$');
    expect((await facts(page)).reflowError).toBeNull();
    await expectLiveEqualsGenerate(page, { label: 'visibility rule during a Legend font size apply' });
  } finally {
    await page.context().close();
  }
});
