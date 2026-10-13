// Layout edit drags (PD-OI-053: only a real drag is an operation; R1, R11).
// A press that moves nothing leaves the stored Result, History, and the
// mounted SVG as they were (OV-608). A drag that moves is one step, leaves no
// inline opacity behind, and Undo restores the stored Result of before it.
// The pointer shows `grab` on a target, `grabbing` while it is pressed, and
// `grab` after, also on a Legend row and a record definition (OV-603, OV-604);
// a record definition that no press drags shows no `grab` (OV-625).
const { test, expect } = require('@playwright/test');
const { BATCH_FIXTURE, HMMT_SESSION, loadSessionFile, openFresh, openWithGenBank, settle } = require('./helpers/audit-browser.cjs');
const { generateAndWaitForResult } = require('./helpers/app-lifecycle.cjs');

test.describe.configure({ retries: 0 });

const LINEAR_SESSION = 'gbdraw/web/gallery/sessions/BGC0000708-BGC0000713.gbdraw-session.json';
const DRAG = { x: 24, y: 12 };

// Pans the first drawn element of `selector` to the canvas centre and returns
// a point where the pointer hits it.
const pointOn = async (page, selector) => {
  const centre = () => page.evaluate((query) => {
    const element = [...document.querySelectorAll(`.gbdraw-preview-surface ${query}`)]
      .find((candidate) => candidate.getBoundingClientRect().width > 0);
    if (!element) return null;
    const box = element.getBoundingClientRect();
    // The middle first, then outwards: another element may cover a part.
    for (let step = 0; step <= 40; step += 1) {
      const x = box.x + box.width * (0.5 + (step % 2 ? 1 : -1) * Math.ceil(step / 2) * 0.01);
      const y = box.y + box.height / 2;
      if (element.contains(document.elementFromPoint(x, y))) return { x, y, hit: true };
    }
    return { x: box.x + box.width / 2, y: box.y + box.height / 2, hit: false };
  }, selector);
  const found = await centre();
  expect(found, selector).not.toBeNull();
  await page.evaluate(({ x, y }) => {
    const app = window.__GBDRAW_APP__;
    const view = app.canvasContainerRef.getBoundingClientRect();
    app.canvasPan.x += view.x + view.width / 2 - x;
    app.canvasPan.y += view.y + view.height / 2 - y;
  }, found);
  // The pan is a CSS transition: wait until the element rests under a point.
  let point = null;
  await expect.poll(async () => {
    const previous = point;
    point = await centre();
    return Boolean(previous && point.hit && previous.x === point.x && previous.y === point.y);
  }, { message: `${selector} rests under the pointer` }).toBe(true);
  return point;
};

// Keeps the stored Results and the mounted SVG in the page (they are large)
// and returns the History depth.
const remember = (page) => page.evaluate(() => {
  const app = window.__GBDRAW_APP__;
  window.__LAYOUT_PRESS_BEFORE__ = {
    results: app.results.map((result) => result.content),
    mounted: app.svgContainer.querySelector('svg').outerHTML
  };
  return window.__GBDRAW_HISTORY__.getUndoCount();
});

const observe = (page) => page.evaluate(() => {
  const app = window.__GBDRAW_APP__;
  const before = window.__LAYOUT_PRESS_BEFORE__;
  const firstDifference = (left, right) => {
    if (left === right) return null;
    let index = 0;
    while (index < left.length && left[index] === right[index]) index += 1;
    return { index, before: left.slice(Math.max(0, index - 80), index + 80), after: right.slice(Math.max(0, index - 80), index + 80) };
  };
  const results = app.results.map((result) => result.content);
  const svg = app.svgContainer.querySelector('svg');
  return {
    undo: window.__GBDRAW_HISTORY__.getUndoCount(),
    storedDifference: firstDifference(before.results.join('\n'), results.join('\n')),
    // An emptied `style=""` (a feature the pointer passed over) is no inline style.
    mountedDifference: firstDifference(before.mounted.replaceAll(' style=""', ''), svg.outerHTML.replaceAll(' style=""', '')),
    storedStyles: results.join('\n').match(/\sstyle="[^"]*"/g) || [],
    mountedOpacityStyles: [...svg.querySelectorAll('[data-gbdraw-composition-role]')]
      .filter((element) => element.style.opacity !== '')
      .map((element) => `${element.id}: ${element.getAttribute('style')}`)
  };
});

const layoutEditOn = async (page) => {
  const toggle = page.getByRole('button', { name: 'Toggle layout edit mode' });
  await toggle.click();
  await expect(toggle).toHaveAttribute('aria-pressed', 'true');
  await settle(page);
};

// The cursor the pointer shows: that of the element under it.
const cursorAt = (page, point) => page.evaluate(({ x, y }) => (
  getComputedStyle(document.elementFromPoint(x, y)).cursor
), point);

// A target shows `grab`; a press on it starts a Layout edit drag, which shows
// `grabbing` (PD-OI-053; a Legend row too, OV-603).
const pressDown = async (page, point) => {
  expect(await cursorAt(page, point), 'a Layout edit target shows grab').toBe('grab');
  await page.mouse.down();
  expect(await cursorAt(page, point), 'the press starts a Layout edit drag').toBe('grabbing');
};

const history = (page, method) => page.evaluate(async (name) => {
  await window.__GBDRAW_APP__[name]();
}, method);

const expectPressAndDrag = async (page, targets) => {
  // With Layout edit off a Legend row is the button of its popup (OV-603).
  const row = await pointOn(page, targets.legend);
  await page.mouse.move(row.x, row.y);
  expect(await cursorAt(page, row)).toBe('pointer');
  await layoutEditOn(page);

  for (const [kind, selector] of Object.entries(targets)) {
    await test.step(`${kind}: a press that moves nothing changes nothing`, async () => {
      const point = await pointOn(page, selector);
      await page.mouse.move(point.x, point.y);
      const undo = await remember(page);
      await pressDown(page, point);
      await page.mouse.up();
      await settle(page);
      expect(await observe(page)).toMatchObject({ undo, storedDifference: null, mountedDifference: null });

      // The same after the pointer left the press point and came back to it.
      await pressDown(page, point);
      await page.mouse.move(point.x + DRAG.x, point.y + DRAG.y, { steps: 4 });
      await page.mouse.move(point.x, point.y, { steps: 4 });
      await page.mouse.up();
      await settle(page);
      expect(await observe(page)).toMatchObject({ undo, storedDifference: null, mountedDifference: null });
    });
  }

  const drag = async (selector) => {
    const point = await pointOn(page, selector);
    await page.mouse.move(point.x, point.y);
    const undo = await remember(page);
    await pressDown(page, point);
    await page.mouse.move(point.x + DRAG.x, point.y + DRAG.y, { steps: 6 });
    await page.mouse.up();
    await settle(page);
    expect(await cursorAt(page, { x: point.x + DRAG.x, y: point.y + DRAG.y }), 'after the drag').toBe('grab');
    return undo;
  };

  // The first edit of a Result replaces Python's serialization with the
  // browser's, so the byte comparisons start from an edited Result.
  await drag(Object.values(targets)[0]);
  await history(page, 'undoHistory');
  await settle(page);

  for (const [kind, selector] of Object.entries(targets)) {
    await test.step(`${kind}: a drag is one step without a leftover style, and Undo restores the Result`, async () => {
      const undo = await drag(selector);
      const dragged = await observe(page);
      expect(dragged.undo).toBe(undo + 1);
      expect(dragged.storedDifference).not.toBeNull();
      expect(dragged.storedStyles).toEqual([]);
      expect(dragged.mountedOpacityStyles).toEqual([]);

      await history(page, 'undoHistory');
      await settle(page);
      expect(await observe(page)).toMatchObject({ undo, storedDifference: null, storedStyles: [] });
      await history(page, 'redoHistory');
      await settle(page);
      expect((await observe(page)).undo).toBe(undo + 1);
      await history(page, 'undoHistory');
      await settle(page);
    });
  }
};

test('Circular: an unmoved press commits nothing and a drag is one clean step', async ({ page }) => {
  test.setTimeout(600_000);
  await openFresh(page);
  await loadSessionFile(page, HMMT_SESSION);
  await page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    app.form.plot_title = 'Layout edit press';
    app.adv.plot_title_position = 'top';
  });
  await generateAndWaitForResult(page);
  await settle(page);

  await expectPressAndDrag(page, {
    diagram: 'g[data-gbdraw-slot-renderer="ticks"] text',
    legend: '#legend text',
    'plot title': '#plot_title text'
  });
});

// A record definition inside a record group is a handle of that record (OV-604).
test('Circular Multi-Record Canvas: a record definition shows the drag of its record', async ({ page }) => {
  test.setTimeout(600_000);
  await openWithGenBank(page, BATCH_FIXTURE, () => {
    const app = window.__GBDRAW_APP__;
    app.form.multi_record_canvas = true;
    app.adv.circular_grouping_intent = 'grid';
  });
  await generateAndWaitForResult(page);
  await settle(page);
  await layoutEditOn(page);

  const point = await pointOn(page, '[data-gbdraw-role="record-definition"] text');
  await page.mouse.move(point.x, point.y);
  await pressDown(page, point);
  await page.mouse.up();
  expect(await cursorAt(page, point), 'after the press').toBe('grab');
});

test('Linear: an unmoved press commits nothing and a drag is one clean step', async ({ page }) => {
  test.setTimeout(600_000);
  await openFresh(page);
  await loadSessionFile(page, LINEAR_SESSION);

  await expectPressAndDrag(page, {
    // The axis of a record: with several records a drag moves one record.
    record: 'g[data-record-key] > line:not([data-label-feature-id])',
    legend: '#legend text',
    'plot title': '#plot_title text',
    'Linear scale': '#length_bar text'
  });

  // A record definition is no drag target here: it shows no grab, and a
  // press on it moves nothing and opens nothing (OV-625).
  const definition = await pointOn(page, '[data-gbdraw-role="record-definition"] text');
  await page.mouse.move(definition.x, definition.y);
  expect(await cursorAt(page, definition)).toBe('default');
  const undo = await remember(page);
  await page.mouse.down();
  await page.mouse.move(definition.x + DRAG.x, definition.y + DRAG.y, { steps: 6 });
  await page.mouse.up();
  await settle(page);
  expect(await observe(page)).toMatchObject({ undo, storedDifference: null, mountedDifference: null });
  await expect(page.locator('[data-preview-item-popup]')).toHaveCount(0);
});
