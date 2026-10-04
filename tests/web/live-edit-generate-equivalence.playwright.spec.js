// G-A (Web GUI audit 2026-09-30; R1, R2, PD-OI-059/060): a live editor edit
// draws what Generate draws from the same draft, and Save, a fresh Load and
// Generate keep it (IN-01, GE-02, PV-10 class). The fixtures carry a Circular
// grid, and a Linear crop and reverse complement. The comparison covers what an
// edit changes (see editProjection), not label or legend geometry. Edits that
// apply on Generate (title, definition, global stroke) are covered by
// history-generated-authority.
const { test, expect } = require('@playwright/test');
const { readFileSync } = require('node:fs');
const { generateAndWaitForResult } = require('./helpers/app-lifecycle.cjs');
const { BATCH_FIXTURE, editFeature, openFresh, openWithGenBank, settle } = require('./helpers/audit-browser.cjs');
const { download, load } = require('./helpers/mode-transition.cjs');

test.describe.configure({ retries: 0 });

const [TESTA, TESTB] = readFileSync(BATCH_FIXTURE, 'utf8').split(/^\/\/\s*$/m)
  .map((chunk) => chunk.trim()).filter((chunk) => chunk.startsWith('LOCUS')).map((chunk) => `${chunk}\n//\n`);

const CIRCULAR = 'Circular grid';
const LINEAR = 'Linear crop and reverse complement';
const MODES = [
  {
    name: CIRCULAR,
    open: (page) => openWithGenBank(page, BATCH_FIXTURE, () => {
      Object.assign(window.__GBDRAW_APP__.form, { labels_mode: 'out', multi_record_canvas: true });
    })
  },
  {
    name: LINEAR,
    open: async (page) => {
      await openFresh(page);
      await page.getByRole('button', { name: 'Linear', exact: true }).click();
      await page.waitForFunction(() => window.__GBDRAW_APP__?.mode === 'linear');
      await page.evaluate(async (records) => {
        const app = window.__GBDRAW_APP__;
        while (app.linearSeqs.length < records.length) app.addLinearSeq();
        records.forEach((text, index) => app.setLinearSeqPrimaryFile(index, 'gb', new File(
          [text], `record-${index}.gb`, { type: 'text/plain', lastModified: 1000 + index }
        )));
        await window.Vue.nextTick();
      }, [TESTA, TESTB]);
      await settle(page);
      await page.evaluate(() => {
        const app = window.__GBDRAW_APP__;
        app.setLinearRecordCrop(app.linearSeqs[0], 'region_start', 201);
        app.setLinearRecordCrop(app.linearSeqs[0], 'region_end', 3800);
        app.linearSeqs[1].region_reverse = true;
        app.form.show_labels_linear = 'all';
      });
      await settle(page);
    }
  }
];

const legendCaption = (page, caption) => page.evaluate((wanted) => {
  const index = window.__GBDRAW_APP__.legendEntries.findIndex((entry) => entry.caption === wanted);
  if (index < 0) throw new Error(`no legend entry ${wanted}`);
  return index;
}, caption);

// [name, apply, { legend }]
const EDITS = [
  // Generate turns a This feature only color into a specific rule and renames the
  // default caption (docs/REFERENCE/palettes-feature-rules-labels-shapes-and-tracks.md).
  ['feature fill, this feature only', (page) => editFeature(page, 'TESTA_0001', { fill: '#c83366' }),
    { legend: false }],
  ['label text', (page) => editFeature(page, 'TESTB_0002', { labelText: 'LIVE_EDIT_LABEL' })],
  ['label hidden', (page) => editFeature(page, 'TESTA_0005', { labelVisibility: 'off' })],
  ['feature hidden', (page) => editFeature(page, 'TESTB_0006', { visibility: 'off' })],
  ['legend entry color', (page) => page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    const index = app.legendEntries.findIndex((entry) => entry.caption);
    if (index < 0) throw new Error('no legend entry');
    return app.updateLegendEntryColor(index, '#2a9d8f');
  })],
  ['legend rename', async (page) => {
    const index = await legendCaption(page, 'tRNA');
    await page.evaluate((entry) => window.__GBDRAW_APP__.renameLegendEntry(entry, 'transfer RNA'), index);
  }],
  ['legend order', async (page) => {
    if (!await page.evaluate(() => window.__GBDRAW_APP__.showRightDrawer)) await page.locator('.drawer-toggle').click();
    await page.locator('.right-drawer').getByRole('button', { name: 'Legend' }).click();
    await page.locator('.right-drawer').getByTitle('Sort Z-A', { exact: true }).click();
  }],
  // N-06 (PD-OI-042): a rule captioned like the generated `other proteins` row
  // draws its own `other proteins [#hex]` row; that row then edits the rule.
  ['rule caption naming a generated row', (page) => page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    Object.assign(app.newSpecRule, {
      feat: 'CDS', qual: 'locus_tag', val: '^TESTA_0006$', color: '#e63946', cap: 'other proteins'
    });
    return app.addSpecificRule();
  })],
  ['suffixed legend row color', async (page) => {
    const index = await legendCaption(page, 'other proteins [#e63946]');
    await page.evaluate((entry) => window.__GBDRAW_APP__.updateLegendEntryColor(entry, '#7b2cbf'), index);
  }]
];

// What an edit changes: paint and visibility of each drawn feature part, the
// visible label texts, the legend caption -> swatch color map, and the legend
// reading order. Label and legend geometry is left out (positions only order
// the captions): the live layout measures rendered text, Python uses its own
// font metrics (the PV-10 measurement, PD-OI-084).
const editProjection = (page) => page.evaluate(() => {
  const app = window.__GBDRAW_APP__;
  // The displayed SVG; Save and export serialize it.
  const root = app.svgContainer.querySelector('svg');
  const hidden = (element) => {
    for (let node = element; node && node.nodeType === 1; node = node.parentNode) {
      const style = node.getAttribute('style') || '';
      if (node.getAttribute('display') === 'none' || /display\s*:\s*none/.test(style)
        || node.getAttribute('visibility') === 'hidden') return true;
    }
    return false;
  };
  const paint = (element, name) => {
    const inline = (element.getAttribute('style') || '').match(new RegExp(`(?:^|;)\\s*${name}\\s*:\\s*([^;]+)`));
    return (inline ? inline[1] : element.getAttribute(name) || '').trim().toLowerCase();
  };
  const features = [...root.querySelectorAll('[data-gbdraw-feature-id]')]
    .filter((element) => !hidden(element) && element.localName === 'path')
    .map((element) => [element.getAttribute('data-gbdraw-feature-id'), element.getAttribute('data-gbdraw-feature-part') || '',
      paint(element, 'fill'), paint(element, 'stroke')].join(' '))
    .sort();
  const entries = [...root.querySelectorAll('g[data-legend-key]')].filter((entry) => !hidden(entry)).map((entry) => {
    const text = entry.querySelector('text');
    const swatch = [...entry.querySelectorAll('path,rect')].find((element) => !['', 'none'].includes(paint(element, 'fill')));
    const box = text?.getBoundingClientRect();
    return { caption: text?.textContent || entry.getAttribute('data-legend-key'), fill: swatch ? paint(swatch, 'fill') : '',
      x: box?.x || 0, y: box?.y || 0 };
  });
  const legend = entries.map(({ caption, fill }) => `${caption} = ${fill}`).sort();
  // Reading order of the captions (rows top to bottom, then left to right).
  const legendOrder = [...entries].sort((left, right) => (Math.abs(left.y - right.y) < 2 ? left.x - right.x : left.y - right.y))
    .map(({ caption }) => caption);
  const labels = [...root.querySelectorAll('text')]
    .filter((element) => !hidden(element) && !element.closest('g[data-legend-key]'))
    .map((element) => [element.closest('[data-label-feature-id]')?.getAttribute('data-label-feature-id')
      || element.closest('[data-gbdraw-feature-id]')?.getAttribute('data-gbdraw-feature-id'), element.textContent.trim()])
    .filter(([feature, text]) => feature && text).map((pair) => pair.join(' ')).sort();
  return { features, labels, legend, legendOrder: [legendOrder.join(' | ')] };
});
const changes = (before, after) => Object.fromEntries(Object.keys(before).map((key) => [key, {
  removed: before[key].filter((item) => !after[key].includes(item)).slice(0, 6),
  added: after[key].filter((item) => !before[key].includes(item)).slice(0, 6)
}]).filter(([, { removed, added }]) => removed.length || added.length));

for (const mode of MODES) {
  test(`${mode.name}: each live edit draws what Generate draws, and Save and Load keep it`, async ({ page, browser }, testInfo) => {
    test.setTimeout(900_000);
    await mode.open(page);
    await generateAndWaitForResult(page);
    await settle(page);
    let previous = await editProjection(page);
    for (const [name, apply, { legend = true } = {}] of EDITS) {
      await test.step(name, async () => {
        await apply(page);
        await settle(page);
        const live = await editProjection(page);
        const overrides = () => page.evaluate(() => JSON.stringify({
          colors: window.__GBDRAW_APP__.featureColorOverrides,
          featureOverrides: window.__GBDRAW_APP__.featureOverrides
        }));
        expect.soft(changes(previous, live), `${mode.name}: ${name} changed the live Result; overrides ${await overrides()}`)
          .not.toEqual({});
        await generateAndWaitForResult(page);
        await settle(page);
        const generated = await editProjection(page);
        const differ = changes(live, generated);
        if (!legend) {
          delete differ.legend;
          delete differ.legendOrder;
        }
        expect.soft(differ, `${mode.name}: ${name}: live Result and Generate differ`).toEqual({});
        previous = generated;
      });
    }
    const saved = testInfo.outputPath('live-edits.gbdraw-session.json');
    await download(page, 'Save Session', saved);
    const reloaded = await load(browser, saved);
    await settle(reloaded);
    await generateAndWaitForResult(reloaded);
    await settle(reloaded);
    expect(changes(previous, await editProjection(reloaded)), `${mode.name}: Save, Load and Generate`).toEqual({});
    await reloaded.context().close();
  });
}
