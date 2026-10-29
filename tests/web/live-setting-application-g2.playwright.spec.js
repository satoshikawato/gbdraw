// G2 live settings (PD-OI-103 scenario revision 2, PD-OI-059 scenario revision
// 3, PD-OI-066 scenario revision 3): each G2 unit redraws the shown Result
// through the automatic rerender, and the live Result equals a Generate from
// the same draft (R3). Values that still wait show "Applies on Generate" on
// their own item and leave the Result as it is until Generate. Controls are
// selected by CSS; each case is one -g selectable unit and mode.
const { test, expect } = require('@playwright/test');
const { readFileSync, writeFileSync } = require('node:fs');
const { generateAndWaitForResult, reveal, settleLive } = require('./helpers/app-lifecycle.cjs');
const { load } = require('./helpers/mode-transition.cjs');
const { expectLiveEqualsGenerate } = require('./helpers/live-generate-parity.cjs');
const { shownSavedDrawing } = require('./helpers/saved-session.cjs');

const SESSIONS = 'gbdraw/web/gallery/sessions';
const TOBACCO = `${SESSIONS}/tobacco-chloroplast.gbdraw-session.json`;
const LAMBDA = `${SESSIONS}/lambda_basic_linear.gbdraw-session.json`;
const BGC = `${SESSIONS}/BGC0000708-BGC0000713.gbdraw-session.json`;
const VIBRIO = `${SESSIONS}/vibrio-harveyi-group-collinear.gbdraw-session.json.gz`;
// An apply takes seconds; a base without the G2 rows fails each live case at this bound.
const APPLY_TIMEOUT = 120_000;
// How long a waiting value is given to (wrongly) redraw the Result.
const WAIT_WINDOW = 5_000;

const generateNow = async (page) => {
  await generateAndWaitForResult(page);
  await settleLive(page);
};

// The committed request and the shown Result.
const facts = (page) => page.evaluate(async () => {
  const app = window.__GBDRAW_APP__;
  const { state } = await import('./js/state.js');
  const { getCommittedCanonicalRenderRequest } = await import('./js/services/config.js');
  return {
    generationKey: state.resultGenerationKey.value,
    request: JSON.parse(JSON.stringify(getCommittedCanonicalRenderRequest() ?? null)),
    content: app.currentResults()[app.selectedResultIndex]?.content ?? null,
    reflowError: app.labelReflowLastError ? JSON.stringify(app.labelReflowLastError) : null
  };
});

// Commits a value as the user does: type it and leave the field.
const commitValue = async (locator, value) => {
  await reveal(locator);
  await locator.scrollIntoViewIfNeeded();
  await locator.fill(String(value));
  await locator.press('Tab');
};

// A color picked without a pointer (a keyboard-opened picker): `fill` sends
// input and change, which open and close one History step (OV-2243).
const commitColor = async (locator, value) => {
  await reveal(locator);
  await locator.scrollIntoViewIfNeeded();
  await locator.fill(value);
};

const numberAfter = (value, fallback, step) => (String(value).trim() === '' ? fallback : Number(value) + step);

// Changes one control, waits until the committed Session holds the value
// (`applied`), then: no Generate ran, the Result was redrawn, and a Generate
// from the same draft draws the same Result, semantically and byte for byte.
const expectLiveUnit = async (page, label, change, applied) => {
  await settleLive(page);
  const before = await facts(page);
  await change();
  await expect.poll(async () => applied((await facts(page)).request), { timeout: APPLY_TIMEOUT }).toBe(true);
  await settleLive(page);
  const live = await facts(page);
  expect(live.generationKey, 'no Generate ran').toBe(before.generationKey);
  expect(live.reflowError).toBeNull();
  expect(live.content, 'the live apply redrew the Result').not.toBe(before.content);
  await expectLiveEqualsGenerate(page, { label });
  expect((await facts(page)).content, 'the Generate draws the live Result byte for byte').toBe(live.content);
};

// The Result and the committed request stay as they are for a waiting value.
const expectWaits = async (page, before) => {
  await page.waitForTimeout(WAIT_WINDOW);
  await settleLive(page);
  const after = await facts(page);
  expect(after.content, 'the Result waits for Generate').toBe(before.content);
  expect(after.request, 'the committed Session waits for Generate').toEqual(before.request);
};

const noteOf = (scope) => scope.locator('.apply-note').first();

// Turns "Use custom stack" on in a mode; true when it was off.
const enableCustomStack = async (page, mode) => {
  const toggle = page.locator(`button[aria-controls="${mode}-custom-track-slots-panel"]`);
  await reveal(toggle);
  await toggle.scrollIntoViewIfNeeded();
  if (await toggle.getAttribute('aria-expanded') !== 'true') await toggle.click();
  const box = toggle.locator('xpath=ancestor::div[.//input[@type="checkbox"]/parent::label][1]')
    .locator('label input[type="checkbox"]').first();
  if (await box.isChecked()) return false;
  await box.click();
  await settleLive(page);
  return true;
};

// The Linear custom-stack rows: id, renderer, and Height, by the row's slot id.
const linearSlotRows = (page) => page.evaluate(() => [...document.querySelectorAll('input[aria-label^="Linear track slot "][aria-label$=" Height"]')]
  .map((input) => {
    const id = input.getAttribute('aria-label').replace(/^Linear track slot /, '').replace(/ Height$/, '');
    const renderer = document.querySelector(`select[aria-label="Linear track renderer ${CSS.escape(id)}"]`);
    return { id, renderer: renderer?.value ?? null, height: input.value };
  }));
const slotHeight = (page, id) => page.locator(`input[aria-label="Linear track slot ${id} Height"]`);
const slotScope = (page, id) => page.locator('[data-setting="track-geometry"]').filter({ has: slotHeight(page, id) });

const run = (name, file, body, timeout = 600_000) => test(name, async ({ browser }, testInfo) => {
  test.setTimeout(timeout);
  const page = await load(browser, typeof file === 'function' ? file(testInfo) : file);
  try {
    await settleLive(page);
    await body(page, testInfo);
  } finally {
    await page.context().close();
  }
});

// Titles definition rows: Default font size in each mode.
run('definition-text circular: a Default font size change is live and equals a Generate', TOBACCO, async (page) => {
  const input = page.locator('input[aria-label="Default font size"]');
  await reveal(input);
  const size = numberAfter(await input.inputValue(), 20, 4);
  await expectLiveUnit(page, 'Circular definition font size', () => commitValue(input, size),
    (request) => request.diagramOptions.configOverrides?.['objects.definition.circular.font_size'] === size);
});
run('definition-text linear: a Default font size change is live and equals a Generate', LAMBDA, async (page) => {
  const input = page.locator('input[aria-label="Default font size"]');
  await reveal(input);
  const size = numberAfter(await input.inputValue(), 20, 4);
  await expectLiveUnit(page, 'Linear definition font size', () => commitValue(input, size),
    (request) => request.diagramOptions.configOverrides?.['objects.definition.linear.font_size.short'] === size);
});

// A Linear row's Size input stores its text; the live apply sends the number a
// Generate sends (OV-2244: it sent the text and the apply failed validation).
run('definition-text linear: a Name / Species font size change is live and equals a Generate', VIBRIO, async (page) => {
  const input = page.locator('input[aria-label="Name / Species font size"]');
  await reveal(input);
  const size = numberAfter(await input.inputValue(), 20, 4);
  await expectLiveUnit(page, 'Linear Name / Species font size', () => commitValue(input, size),
    (request) => request.diagramOptions.configOverrides?.['objects.definition.linear.line_styles.name.font_size'] === size);
}, 900_000);

run('species circular: a Species change is live and equals a Generate', TOBACCO, async (page) => {
  const input = page.locator('#circular-species');
  await reveal(input);
  const species = `${await input.inputValue()} G2`;
  await expectLiveUnit(page, 'Circular species', () => commitValue(input, species),
    (request) => request.diagramOptions.species === species);
});

// Lambda's Session has GC off; the Linear case turns it on and generates first.
for (const [mode, file] of [['circular', TOBACCO], ['linear', LAMBDA]]) {
  run(`dinucleotide-window ${mode}: a Window change is live and equals a Generate`, file, async (page) => {
    if (mode === 'linear') {
      const gc = page.locator('input[aria-describedby="help-linear-show-gc-content"]');
      await reveal(gc);
      if (!(await gc.isChecked())) await gc.click();
      await generateNow(page);
    }
    const input = page.locator('input[aria-label="Window"]');
    await reveal(input);
    const current = await input.inputValue();
    const windowSize = current === '' ? 3000 : Number(current) * 2;
    await expectLiveUnit(page, `${mode} window`, () => commitValue(input, windowSize),
      (request) => request.diagramOptions.window === windowSize);
  });
}

// Region Annotations: a committed item's Fill. Tobacco's plastome_tint items
// draw their fill; lambda has no set, so one is added and generated first.
run('annotation-style circular: an annotation item Fill change is live and equals a Generate', TOBACCO, async (page) => {
  const fill = page.locator('input[aria-label="Annotation fill color"]').nth(4);
  await expectLiveUnit(page, 'Circular annotation fill', () => commitColor(fill, '#ff0000'),
    (request) => request.diagramOptions.annotations.sets[1].annotations[0].style?.fill === '#ff0000');
});
run('annotation-style linear: an annotation item Fill change is live and equals a Generate', LAMBDA, async (page) => {
  const section = page.locator('[data-setting="region-annotations"]');
  await reveal(section.locator('button', { hasText: 'Add set' }).first());
  await section.locator('button', { hasText: 'Add set' }).first().click();
  await section.locator('button', { hasText: 'Coordinates' }).first().click();
  await commitValue(page.locator('input[aria-label="Annotation end (inclusive)"]').first(), 20000);
  const target = page.locator('select[aria-label="Annotation target record"]').first();
  if (await target.count() > 0) {
    const value = await target.evaluate((el) => [...el.options].find((o) => o.value !== '' && !o.disabled)?.value ?? null);
    if (value !== null) await target.selectOption(value);
  }
  await generateNow(page);
  const fill = page.locator('input[aria-label="Annotation fill color"]').first();
  await expectLiveUnit(page, 'Linear annotation fill', () => commitColor(fill, '#ff0000'),
    (request) => request.diagramOptions.annotations?.sets?.[0]?.annotations?.[0]?.style?.fill === '#ff0000');
}, 900_000);

run('comparison-height linear: a Match height change is live and equals a Generate', BGC, async (page) => {
  const input = page.locator('input[aria-label="Comparison match height"]');
  await reveal(input);
  const height = numberAfter(await input.inputValue(), 80, 20);
  await expectLiveUnit(page, 'Linear comparison match height', () => commitValue(input, height),
    (request) => request.diagramOptions.configOverrides?.['canvas.linear.comparison_height'] === height);
});

// Custom track stacks: Circular Center Reserved Radius; Linear Height of a
// committed row that is neither features nor Depth.
run('track-geometry circular: a Center Reserved Radius change is live and equals a Generate', TOBACCO, async (page) => {
  if (await enableCustomStack(page, 'circular')) await generateNow(page);
  const input = page.locator('input[aria-label="Center Reserved Radius"]');
  const current = parseFloat(await input.inputValue());
  const radius = Number.isFinite(current) && current > 0 ? Math.round(current * 2) : 60;
  await expectLiveUnit(page, 'Circular center reserved radius', () => commitValue(input, radius),
    (request) => request.diagramOptions.tracks.centerReservedRadius === radius);
});
run('track-geometry linear: a custom-stack row Height change is live and equals a Generate', LAMBDA, async (page) => {
  if (await enableCustomStack(page, 'linear')) await generateNow(page);
  const row = (await linearSlotRows(page)).find((slot) => slot.renderer !== 'features' && slot.renderer !== 'depth');
  expect(row, 'the stack has a row that is neither features nor Depth').toBeTruthy();
  await expect(noteOf(slotScope(page, row.id))).toHaveText(/Live edit/);
  const before = JSON.stringify((await facts(page)).request.diagramOptions.tracks.linearTrackSlots);
  const current = parseFloat(row.height);
  const height = Number.isFinite(current) && current > 0 ? Math.round(current * 2) : 60;
  await expectLiveUnit(page, 'Linear slot height', () => commitValue(slotHeight(page, row.id), height),
    (request) => JSON.stringify(request.diagramOptions.tracks.linearTrackSlots) !== before);
}, 900_000);

// Record names: the Circular record label (a single-record Result), and a Linear
// record's name edited right after Session Load (no Generate first).
run('record-names circular: a record label change is live and equals a Generate', TOBACCO, async (page) => {
  // The label is editable once the source is inspected and its one record is
  // selected; the selection applies on Generate.
  const input = page.locator('input[aria-label="Circular record label"]');
  const inspect = page.locator('[data-circular-inspect]');
  await reveal(inspect);
  await inspect.click();
  const record = page.locator('select[aria-label="Circular record"]');
  await expect(record.locator('option[value="NC_001879.2"]')).toHaveCount(1, { timeout: APPLY_TIMEOUT });
  await reveal(record);
  await record.selectOption('NC_001879.2');
  await expect(input).toBeEnabled();
  await generateNow(page);
  await expectLiveUnit(page, 'Circular record label', () => commitValue(input, 'G2 record label'),
    (request) => JSON.stringify(request.records).includes('G2 record label'));
});
run('record-names linear: a record name edited after Session Load is live and equals a Generate', BGC, async (page) => {
  const card = page.locator('[data-linear-source-card]').first();
  await expect(noteOf(card.locator('[data-setting="record-names"]'))).toHaveText(/Live edit/);
  const input = page.locator('input[aria-label="Definition for sequence 1"]');
  await expectLiveUnit(page, 'Linear record name', () => commitValue(input, 'G2 record name'),
    (request) => JSON.stringify(request.records).includes('G2 record name'));
});

// Values that wait: an annotation item with a Legend label (O-8), a features
// row's Height, and a record name in a card whose File was replaced.
const tobaccoWithLegendItem = (testInfo) => {
  const session = JSON.parse(readFileSync(TOBACCO, 'utf8'));
  const irb = (sets) => sets.find((set) => set.id === 'plastome_tint').annotations.find((item) => item.id === 'irb_tint');
  const drawing = shownSavedDrawing(session);
  irb(drawing.renderRequest.diagramOptions.annotations.sets).legendLabel = 'Inverted repeat';
  irb(drawing.config.annotationSets).legendLabel = 'Inverted repeat';
  const file = testInfo.outputPath('legend-item.gbdraw-session.json');
  writeFileSync(file, JSON.stringify(session));
  return file;
};
run('annotation-style waits: an item with a Legend label says Applies on Generate and keeps the Result', tobaccoWithLegendItem, async (page) => {
  const rows = page.locator('[data-annotation-row]');
  const scopeOf = (index) => rows.nth(index).locator('[data-setting="annotation-style"]');
  await reveal(rows.nth(5));
  // plastome_tint: lsc_tint (row 5) has no Legend label, irb_tint (row 6) has one.
  await expect(noteOf(scopeOf(4).first())).toHaveText(/Live edit/);
  await expect(noteOf(scopeOf(5).first())).toHaveText(/Applies on Generate/);
  const before = await facts(page);
  await commitColor(page.locator('input[aria-label="Annotation fill color"]').nth(5), '#ff0000');
  await expectWaits(page, before);
  await generateNow(page);
  expect((await facts(page)).request.diagramOptions.annotations.sets[1].annotations[1].style.fill).toBe('#ff0000');
});
run('track-geometry waits: a features row Height says Applies on Generate and keeps the Result', LAMBDA, async (page) => {
  if (await enableCustomStack(page, 'linear')) await generateNow(page);
  const slots = await linearSlotRows(page);
  const features = slots.find((slot) => slot.renderer === 'features');
  const other = slots.find((slot) => slot.renderer !== 'features' && slot.renderer !== 'depth');
  expect(features, 'the stack has a features row').toBeTruthy();
  await expect(noteOf(slotScope(page, features.id))).toHaveText(/Applies on Generate/);
  if (other) await expect(noteOf(slotScope(page, other.id))).toHaveText(/Live edit/);
  const before = await facts(page);
  const current = parseFloat(features.height);
  await commitValue(slotHeight(page, features.id), Number.isFinite(current) && current > 0 ? Math.round(current * 2) : 60);
  await expectWaits(page, before);
}, 900_000);
run('record-names waits: a record name in a card whose File was replaced keeps the Result', BGC, async (page, testInfo) => {
  const session = JSON.parse(readFileSync(BGC, 'utf8'));
  const resource = session.resources['record-1-genbank'];
  const file = testInfo.outputPath(resource.name);
  writeFileSync(file, Buffer.from(resource.data, 'base64'));
  const card = page.locator('[data-linear-source-card]').first();
  await reveal(card.locator('input[type="file"]').first());
  await card.locator('input[type="file"]').first().setInputFiles(file);
  await settleLive(page);
  await expect(noteOf(card.locator('[data-setting="record-names"]'))).toHaveText(/Applies on Generate/);
  const before = await facts(page);
  await commitValue(page.locator('input[aria-label="Definition for sequence 1"]'), 'G2 replaced name');
  await expectWaits(page, before);
});
