// PD-OI-066 (R1, R3, R10): changes of track data and mode switches leave the
// live Result equal to the next Generate: Undo and Redo of a Depth source
// removal (OV-66), Legend styles that follow the captions of track data
// (OV-65), a decoration moved on one mode's Result (OV-104), and the track
// toggles of the mode not shown. The matrix of edit kinds is in
// live-generate-parity.playwright.spec.js.
const { test, expect } = require('@playwright/test');
const { readFileSync } = require('node:fs');
const { reveal } = require('./helpers/app-lifecycle.cjs');
const { openWithGenBank } = require('./helpers/audit-browser.cjs');
const { expectLiveEqualsGenerate, semanticSnapshot, settleLive } = require('./helpers/live-generate-parity.cjs');
const {
  SINGLE_FIXTURE, open, generate, history, DEPTH_TSV, colorLegendRow, openCanvas, renameRow
} = require('./helpers/live-generate-parity-steps.cjs');

test.describe.configure({ retries: 0 });

// OV-66: Undo of a Depth source removal brings back the Depth track and its
// tick text, and Redo hides them again; the live Result equals what Generate
// draws. A Depth source History step restores files, which suppresses the
// track-visibility watcher, so the step projects the visibility itself.
const DEPTH_CASES = {
  'circular, clearing the Depth file': {
    mode: 'circular',
    add: (app, file) => app.setCircularDepthFile(0, file),
    remove: (app) => app.setCircularDepthFile(0, null)
  },
  'circular, removing the Depth track': {
    mode: 'circular',
    add: (app, file) => app.setCircularDepthFile(0, file),
    remove: (app) => app.removeCircularDepthTrack(0)
  },
  'linear, clearing the Depth file': {
    mode: 'linear',
    add: (app, file) => app.setLinearDepthFile(app.linearSeqs[0], 0, file),
    remove: (app) => app.setLinearDepthFile(app.linearSeqs[0], 0, null)
  }
};
for (const [name, { mode, add, remove }] of Object.entries(DEPTH_CASES)) {
  test(`Undo and Redo of a Depth source removal match Generate (${name})`, async ({ page }) => {
    test.setTimeout(180_000);
    await open(page, { mode, results: 'single', reflow: 'off' });
    const inStep = async (label, change, ...args) => {
      await page.evaluate(async ({ stepLabel, source, values }) => {
        const run = new Function('app', 'text', `return (${source})(app, text && new File([text], 'depth.tsv', { type: 'text/tab-separated-values' }));`);
        await window.__GBDRAW_HISTORY__.runUndoable(stepLabel, () => run(window.__GBDRAW_APP__, values[0]));
      }, { stepLabel: label, source: change.toString(), values: args });
      await settleLive(page);
    };
    await inStep('Change uploaded file', add, DEPTH_TSV);
    await generate(page);
    await inStep('Remove Depth', remove, null);
    await history(page, 'undo');
    await expectLiveEqualsGenerate(page, { label: `${name}: Undo` });
    await history(page, 'redo');
    await expectLiveEqualsGenerate(page, { label: `${name}: Redo` });
  });
}

// OV-65: a Legend color on a row named only by a track's data (an annotation set,
// a depth file) follows the caption: a data change retires the styles of the
// captions the data no longer names, in the History step of the change, so
// Undo brings back the data and the style (RETIRING_CASES, further below).
// A region annotation with a legend label draws a Legend row from its set; the
// slot of the set is added through the track slot control.
const addAnnotationRow = async (page, label) => {
  await page.evaluate((legendLabel) => {
    const app = window.__GBDRAW_APP__;
    const set = app.addAnnotationSet('regions');
    const annotation = app.addCoordinateAnnotation(set, { start: 100, end: 400 });
    annotation.legendLabel = legendLabel;
  }, label);
  await settleLive(page);
  await generate(page);
};

// OV-65: Legend styles follow the captions that track data names. Each data
// change below runs in one History step, as its control's does.
const inHistoryStep = (page, label, body, arg) => page.evaluate(
  async ({ stepLabel, source, value }) => {
    const change = new Function('app', 'value', `return (${source})(app, value);`);
    await window.__GBDRAW_HISTORY__.runUndoable(stepLabel, () => change(window.__GBDRAW_APP__, value));
  },
  { stepLabel: label, source: body.toString(), value: arg }
);

const legendStyleOf = (page, caption) => page.evaluate(async (target) => {
  const { state } = await import('/gbdraw/web/js/state.js');
  return {
    color: state.activeDrawing().legendColorOverrides[target] ?? null,
    stroke: state.activeDrawing().legendStrokeOverrides[target] ?? null
  };
}, caption);

const undo = async (page) => {
  await page.evaluate(() => window.__GBDRAW_APP__.undoHistory());
  await settleLive(page);
};

// Types a new text into a label field of an opened section, and leaves the
// field, as the reader does; the History step commits when the field loses focus.
const typeIntoLabel = async (page, section, label, text) => {
  await page.evaluate((name) => {
    document.querySelectorAll('details > summary').forEach((summary) => {
      if (summary.textContent.trim().startsWith(name) || summary.getAttribute('aria-label') === name) {
        summary.parentElement.open = true;
      }
    });
  }, section);
  const field = page.getByLabel(label, { exact: true });
  await field.fill(text);
  await field.blur();
  await settleLive(page);
};

const addDepthFile = async (page, name, mode = 'circular') => {
  await inHistoryStep(page, 'Change uploaded file', (app, { text, fileName, linear }) => {
    const file = new File([text], fileName, { type: 'text/tab-separated-values' });
    if (linear) app.setLinearDepthFile(app.linearSeqs[0], 0, file);
    else app.setCircularDepthFile(0, file);
  }, { text: DEPTH_TSV, fileName: name, linear: mode === 'linear' });
  await settleLive(page);
};
const depthCaption = (page) => page.evaluate(() => (
  window.__GBDRAW_APP__.legendEntries.find((entry) => /depth/i.test(entry.caption))?.caption ?? null
));
// The setup of a Depth case: a file `depth.tsv` drawn once, with a Legend color on its row.
const colorDepthRow = (mode = 'circular') => async (page) => {
  await addDepthFile(page, 'depth.tsv', mode);
  await generate(page);
  await colorLegendRow(page, await depthCaption(page));
};
// OV-87: the setup of a renamed case: a row drawn once, renamed in the Legend
// and colored under its new name; `generateAfter` draws the renamed row once.
const renameAndColorRow = (draw, from, to, { generateAfter = false } = {}) => async (page) => {
  await draw(page);
  await renameRow(page, typeof from === 'function' ? await from(page) : from, to);
  await settleLive(page);
  await colorLegendRow(page, to);
  if (generateAfter) await generate(page);
};
const drawDepthFile = (mode = 'circular') => async (page) => {
  await addDepthFile(page, 'depth.tsv', mode);
  await generate(page);
};
const legendNames = (page) => page.evaluate(() => (
  window.__GBDRAW_APP__.legendEntries.map((entry) => `${entry.originalCaption}=>${entry.caption}`)
));
const depthFileName = (page) => page.evaluate(() => {
  const app = window.__GBDRAW_APP__;
  return (app.mode === 'linear' ? app.linearSeqs[0]?.depth?.[0] : app.files.c_depth?.[0]?.[0])?.name ?? null;
});

const ANNOTATION_TSV = (legendLabel) => [
  'set_id\tid\tmark\tstart\tend\tlegend_label',
  `regions\tregion_1\thighlight\t100\t400\t${legendLabel}`
].join('\n');

const drawnLegendCaptions = async (page) => (await semanticSnapshot(page)).legend.map(({ caption }) => caption);
const NO_STYLE = { color: null, stroke: null };
const COLORED = '#7b2cbf';

// Each case runs a data change that removes a caption's rows as one History
// step. The style is retired in that step; Undo brings back the data and the
// style, and the live Result equals Generate. The same change made again on
// the restored data retires the style again, and Generate succeeds and draws
// no row of the data. (One case per row: the open and setup dominate a case's
// time, and the second change starts from the data the setup made.)
const RETIRING_CASES = [
  {
    name: 'removing an annotation set',
    caption: async () => 'Region X',
    setup: async (page) => {
      await addAnnotationRow(page, 'Region X');
      await colorLegendRow(page, 'Region X');
    },
    change: (page) => inHistoryStep(page, 'Delete set', (app) => app.removeAnnotationSet(app.annotationSets[0])),
    restored: (page) => page.evaluate(() => window.__GBDRAW_APP__.annotationSets.length === 1)
  },
  {
    name: 'removing the depth file',
    caption: async (page) => depthCaption(page),
    setup: colorDepthRow(),
    change: (page) => inHistoryStep(page, 'Change uploaded file', (app) => app.setCircularDepthFile(0, null)),
    restored: async (page) => await depthFileName(page) === 'depth.tsv'
  },
  // A Depth track removal and a file of another name reach the same rule.
  {
    name: 'removing the Depth track',
    caption: async (page) => depthCaption(page),
    setup: colorDepthRow(),
    change: (page) => inHistoryStep(page, 'Remove Depth', (app) => app.removeCircularDepthTrack(0)),
    restored: async (page) => await depthFileName(page) === 'depth.tsv'
  },
  {
    name: 'replacing the depth file with a file of another name',
    caption: async (page) => depthCaption(page),
    setup: colorDepthRow(),
    change: (page) => addDepthFile(page, 'coverage.tsv'),
    restored: async (page) => await depthFileName(page) === 'depth.tsv',
    drawn: ['coverage']
  },
  {
    name: 'removing the Depth track in Linear',
    mode: 'linear',
    caption: async (page) => depthCaption(page),
    setup: colorDepthRow('linear'),
    change: (page) => inHistoryStep(page, 'Remove Depth', (app) => app.removeLinearDepthTrack(0)),
    restored: async (page) => await depthFileName(page) === 'depth.tsv'
  },
  {
    name: 'replacing the depth file with a file of another name in Linear',
    mode: 'linear',
    caption: async (page) => depthCaption(page),
    setup: colorDepthRow('linear'),
    change: (page) => addDepthFile(page, 'coverage.tsv', 'linear'),
    restored: async (page) => await depthFileName(page) === 'depth.tsv',
    drawn: ['coverage']
  },
  {
    name: 'replacing the annotation data with another legendLabel',
    caption: async () => 'Region X',
    setup: async (page) => {
      await addAnnotationRow(page, 'Region X');
      await colorLegendRow(page, 'Region X');
    },
    change: (page) => inHistoryStep(page, 'Import annotations', async (app, tsv) => {
      await app.importAnnotationTableFile({ target: { files: [new File([tsv], 'annotations.tsv')], value: '' } });
    }, ANNOTATION_TSV('Region Y')),
    restored: (page) => page.evaluate(() => window.__GBDRAW_APP__.annotationSets[0]?.annotations[0]?.legendLabel === 'Region X'),
    drawn: ['Region Y']
  },
  // OV-67: a label edit is a data change of the caption. The edits go through
  // the label fields a reader types into.
  {
    name: 'editing the legend label of an annotation set',
    caption: async () => 'Region X',
    setup: async (page) => {
      await page.evaluate(() => {
        const app = window.__GBDRAW_APP__;
        const set = app.addAnnotationSet('regions');
        app.addCoordinateAnnotation(set, { start: 100, end: 400 });
        set.legendLabel = 'Region X';
      });
      await settleLive(page);
      await generate(page);
      await colorLegendRow(page, 'Region X');
    },
    change: (page) => typeIntoLabel(page, 'Region Annotations', 'Set legend label', 'Region Z'),
    restored: (page) => page.evaluate(() => window.__GBDRAW_APP__.annotationSets[0]?.legendLabel === 'Region X'),
    drawn: ['Region Z']
  },
  {
    name: 'editing the legend title of a Depth series',
    caption: async (page) => depthCaption(page),
    setup: colorDepthRow(),
    change: (page) => typeIntoLabel(page, 'Depth TSV tracks', 'Depth legend title', 'Coverage'),
    restored: (page) => page.evaluate(() => window.__GBDRAW_APP__.adv.depth_tracks[0]?.label === 'depth'),
    drawn: ['Coverage']
  },
  // OV-87: names follow the caption as styles do. The data change also retires
  // the Legend rename of the row and the styles stored under the new name.
  {
    name: 'removing the depth file of a row renamed in the Legend',
    caption: async () => 'Coverage',
    renamedFrom: 'depth',
    setup: renameAndColorRow(drawDepthFile(), depthCaption, 'Coverage'),
    change: (page) => inHistoryStep(page, 'Change uploaded file', (app) => app.setCircularDepthFile(0, null)),
    restored: async (page) => await depthFileName(page) === 'depth.tsv'
  },
  {
    name: 'removing the Depth track of a row renamed in the Legend',
    caption: async () => 'Coverage',
    renamedFrom: 'depth',
    setup: renameAndColorRow(drawDepthFile(), depthCaption, 'Coverage'),
    change: (page) => inHistoryStep(page, 'Remove Depth', (app) => app.removeCircularDepthTrack(0)),
    restored: async (page) => await depthFileName(page) === 'depth.tsv'
  },
  {
    name: 'replacing the depth file of a row renamed in the Legend with a file of another name',
    caption: async () => 'Coverage',
    renamedFrom: 'depth',
    setup: renameAndColorRow(drawDepthFile(), depthCaption, 'Coverage'),
    change: (page) => addDepthFile(page, 'cov2.tsv'),
    restored: async (page) => await depthFileName(page) === 'depth.tsv',
    drawn: ['cov2']
  },
  {
    name: 'removing the Depth track of a row renamed in the Legend in Linear',
    mode: 'linear',
    caption: async () => 'Coverage',
    renamedFrom: 'depth',
    setup: renameAndColorRow(drawDepthFile('linear'), depthCaption, 'Coverage'),
    change: (page) => inHistoryStep(page, 'Remove Depth', (app) => app.removeLinearDepthTrack(0)),
    restored: async (page) => await depthFileName(page) === 'depth.tsv'
  },
  {
    name: 'removing the depth file of a renamed row drawn once in Linear',
    mode: 'linear',
    caption: async () => 'Coverage',
    renamedFrom: 'depth',
    setup: renameAndColorRow(drawDepthFile('linear'), depthCaption, 'Coverage', { generateAfter: true }),
    change: (page) => inHistoryStep(page, 'Change uploaded file', (app) => app.setLinearDepthFile(app.linearSeqs[0], 0, null)),
    restored: async (page) => await depthFileName(page) === 'depth.tsv'
  },
  {
    name: 'removing an annotation set whose row is renamed in the Legend',
    caption: async () => 'Region Y',
    renamedFrom: 'Region X',
    setup: renameAndColorRow((page) => addAnnotationRow(page, 'Region X'), 'Region X', 'Region Y'),
    change: (page) => inHistoryStep(page, 'Delete set', (app) => app.removeAnnotationSet(app.annotationSets[0])),
    restored: (page) => page.evaluate(() => window.__GBDRAW_APP__.annotationSets.length === 1)
  }
];

test.describe('OV-65 Legend styles follow the captions of track data', () => {
  test.beforeEach(() => { test.setTimeout(180_000); });

  for (const { name, mode = 'circular', caption, renamedFrom = null, setup, change, restored, drawn = [] } of RETIRING_CASES) {
    test(`${name} retires the styles of its rows, Undo restores them as Generate draws, and Generate succeeds`, async ({ page }) => {
      await openCanvas(page, mode, null);
      await setup(page);
      const row = await caption(page);
      expect(row, 'row caption').toBeTruthy();
      const changeRetires = async () => {
        await change(page);
        await settleLive(page);
        expect(await legendStyleOf(page, row)).toEqual(NO_STYLE);
        if (renamedFrom) {
          expect(await legendStyleOf(page, renamedFrom)).toEqual(NO_STYLE);
          expect(await legendNames(page), 'the rename is retired').not.toContain(`${renamedFrom}=>${row}`);
        }
      };
      await changeRetires();
      await undo(page);
      expect(await restored(page), 'data restored').toBe(true);
      expect((await legendStyleOf(page, row)).color).toBe(COLORED);
      if (renamedFrom) expect(await legendNames(page), 'the rename is restored').toContain(`${renamedFrom}=>${row}`);
      await expectLiveEqualsGenerate(page, { label: `${name}, Undo` });

      await changeRetires();
      await generate(page);
      const captions = await drawnLegendCaptions(page);
      expect(captions).not.toContain(row);
      if (renamedFrom) expect(captions).not.toContain(renamedFrom);
      for (const kept of drawn) expect(captions).toContain(kept);
    });
  }

  test('replacing the depth file with the label unchanged keeps the row style', async ({ page }) => {
    await openCanvas(page, 'circular', null);
    await addDepthFile(page, 'depth.tsv');
    await generate(page);
    const caption = await depthCaption(page);
    await colorLegendRow(page, caption);
    await addDepthFile(page, 'depth.tsv');
    expect((await legendStyleOf(page, caption)).color).toBe(COLORED);
    await expectLiveEqualsGenerate(page, { label: 'depth file replaced, same label' });
  });

  test('replacing the depth file with the label unchanged keeps the Legend rename and its style (OV-87)', async ({ page }) => {
    await openCanvas(page, 'circular', null);
    await renameAndColorRow(drawDepthFile(), depthCaption, 'Coverage')(page);
    await addDepthFile(page, 'depth.tsv');
    expect(await legendNames(page)).toContain('depth=>Coverage');
    expect((await legendStyleOf(page, 'Coverage')).color).toBe(COLORED);
    await expectLiveEqualsGenerate(page, { label: 'renamed depth row, file replaced with the same label' });
    expect(await drawnLegendCaptions(page)).toContain('Coverage');
  });

  test('replacing the annotation data with the legendLabel unchanged keeps the row style', async ({ page }) => {
    await openCanvas(page, 'circular', null);
    await addAnnotationRow(page, 'Region X');
    await colorLegendRow(page, 'Region X');
    await inHistoryStep(page, 'Import annotations', async (app, tsv) => {
      await app.importAnnotationTableFile({ target: { files: [new File([tsv], 'annotations.tsv')], value: '' } });
    }, ANNOTATION_TSV('Region X'));
    await settleLive(page);
    expect((await legendStyleOf(page, 'Region X')).color).toBe(COLORED);
    await expectLiveEqualsGenerate(page, { label: 'annotation data replaced, same legendLabel' });
  });

  // OV-68: the generic Track legend label of a slot renames a GC row. Python
  // lists the default GC caption among the rows it can produce (OV-63), so the
  // Legend color on the old caption stays stored, as for a switched-off track,
  // and Generate succeeds. The case guards that.
  test('renaming the GC content row through its Track legend label keeps Generate working', async ({ page }) => {
    await openCanvas(page, 'circular', null);
    await colorLegendRow(page, 'GC content');
    const panel = page.locator('button[aria-controls="circular-custom-track-slots-panel"]');
    await reveal(panel);
    if (await panel.getAttribute('aria-expanded') !== 'true') await panel.click({ timeout: 15_000 });
    await page.getByText('Use custom stack', { exact: true }).locator('input').check({ timeout: 15_000 });
    await settleLive(page);
    const field = page.getByRole('group', { name: 'Circular track slot gc_content', exact: true })
      .getByRole('textbox', { name: 'Track legend label' });
    await field.fill('Custom GC', { timeout: 15_000 });
    await field.blur();
    await settleLive(page);
    await generate(page);
    const captions = await drawnLegendCaptions(page);
    expect(captions).toContain('Custom GC');
    expect(captions).not.toContain('GC content');
  });
});

// OV-104 (PD-OI-052): a Legend or plot title moved on one mode's Result. Each
// mode keeps its own Result (E1), so the other mode's Generate reads only its
// own mode's Results: it succeeds and carries no move, and switching back shows
// the original Result, the same object, with its move.
test.describe('OV-104 a decoration moved on one mode\'s Result', () => {
  const MODE_NAMES = { circular: 'Circular', linear: 'Linear' };
  const showMode = async (page, target) => {
    await page.getByRole('button', { name: MODE_NAMES[target], exact: true }).click();
    await page.waitForFunction((value) => window.__GBDRAW_APP__?.mode === value, target);
    await settleLive(page);
  };
  const moveDecoration = async (page, role) => {
    await page.evaluate(async (targetRole) => {
      const app = window.__GBDRAW_APP__;
      app.layoutRepositionMode = true;
      await window.Vue.nextTick();
      const svg = app.svgContainer.querySelector('svg');
      const target = svg.querySelector(`[data-gbdraw-composition-role="${targetRole}"]`);
      if (!target) throw new Error(`No ${targetRole} to move`);
      const bounds = target.getBoundingClientRect();
      const x = bounds.left + Math.min(bounds.width / 2, 8);
      const y = bounds.top + Math.min(bounds.height / 2, 8);
      const mouse = (type, dx, dy, buttons) => new MouseEvent(type, {
        bubbles: true, cancelable: true, clientX: x + dx, clientY: y + dy, buttons, view: window
      });
      const frame = () => new Promise((resolve) => requestAnimationFrame(resolve));
      target.dispatchEvent(mouse('mousedown', 0, 0, 1));
      await frame();
      const moveTarget = targetRole === 'legend' ? svg : document;
      moveTarget.dispatchEvent(mouse('mousemove', 12, -9, 1));
      await frame();
      moveTarget.dispatchEvent(mouse('mouseup', 12, -9, 0));
    }, role);
    await expect.poll(() => page.evaluate(() => window.__GBDRAW_HISTORY__?.undoLabel?.() || ''))
      .toBe(role === 'title' ? 'Move plot title' : 'Move legend');
    await page.evaluate(() => { window.__GBDRAW_APP__.layoutRepositionMode = false; });
    await settleLive(page);
  };
  const displayedMove = (page, role) => page.evaluate(async (targetRole) => {
    const { compositionUserDeltas } = await import('/gbdraw/web/js/app/legend-layout/composition-actions.js');
    const { getCommittedSvgResultRuntimeIdentity } = await import('/gbdraw/web/js/services/svg-result-ingestion.js');
    const { state } = await import('/gbdraw/web/js/state.js');
    const svg = window.__GBDRAW_APP__.svgContainer.querySelector('svg');
    return {
      mode: state.generatedMode.value,
      identity: getCommittedSvgResultRuntimeIdentity(state.results.value[state.selectedResultIndex.value]),
      delta: compositionUserDeltas(svg)[targetRole]
    };
  }, role);
  const moved = (delta) => Array.isArray(delta) && delta.some((value) => Math.abs(value) > 0.1);
  const expectSameMove = (actual, expected) => {
    expect(actual).toHaveLength(2);
    actual.forEach((value, index) => expect(value).toBeCloseTo(expected[index], 5));
  };

  for (const { from, role } of [
    { from: 'circular', role: 'legend' },
    { from: 'linear', role: 'legend' },
    { from: 'circular', role: 'title' }
  ]) {
    const other = from === 'linear' ? 'circular' : 'linear';
    test(`a ${role} moved on the ${MODE_NAMES[from]} Result stays with it through a ${MODE_NAMES[other]} Generate`, async ({ page }) => {
      test.setTimeout(180_000);
      await open(page, { mode: 'circular', results: 'single', reflow: 'off' });
      await showMode(page, 'linear');
      await page.evaluate(async (text) => {
        window.__GBDRAW_APP__.setLinearSeqPrimaryFile(0, 'gb', new File([text], 'forced_label_underlay.gb', { type: 'text/plain', lastModified: 1000 }));
        await window.Vue.nextTick();
      }, readFileSync(SINGLE_FIXTURE, 'utf8'));
      await settleLive(page);
      await showMode(page, from);
      if (role === 'title') {
        await page.evaluate(() => {
          const app = window.__GBDRAW_APP__;
          app.form.plot_title = 'OV-104 title';
          app.adv.plot_title_position = 'top';
        });
      }
      await generate(page);
      await moveDecoration(page, role);
      const before = await displayedMove(page, role);
      expect(before.mode).toBe(from);
      expect(moved(before.delta), `the ${role} moved`).toBe(true);

      await showMode(page, other);
      await generate(page);
      const otherResult = await displayedMove(page, role);
      expect(otherResult.mode).toBe(other);
      expect(moved(otherResult.delta), `the ${MODE_NAMES[other]} Result carries no move`).toBe(false);

      await showMode(page, from);
      const returned = await displayedMove(page, role);
      expect(returned.mode).toBe(from);
      expect(returned.identity, 'the original Result is shown again').toBe(before.identity);
      expectSameMove(returned.delta, before.delta);

      await generate(page);
      const regenerated = await displayedMove(page, role);
      expect(regenerated.mode).toBe(from);
      expectSameMove(regenerated.delta, before.delta);
    });
  }
});

// OV-83: each mode keeps its own Result (E1). A track toggle of the current mode
// applies at its next Generate when the mode has no Result yet, and the other
// mode's Result waits in its slot unchanged: switching back shows the same
// Result, committed SVG and drawn tracks included.
const TRACK_TOGGLES = {
  circular: {
    'GC content': { form: 'suppress_gc', toggled: true, group: 'content' },
    'GC skew': { form: 'suppress_skew', toggled: true, group: 'skew' }
  },
  linear: {
    'GC content': { form: 'show_gc', toggled: false, group: 'content' },
    'GC skew': { form: 'show_skew', toggled: false, group: 'skew' }
  }
};

// The track group `group` is drawn in the committed Result: present and not
// hidden by display none.
const trackDrawn = (page, group) => page.evaluate((base) => {
  const app = window.__GBDRAW_APP__;
  const root = new DOMParser().parseFromString(app.results[app.selectedResultIndex].content, 'image/svg+xml').documentElement;
  // Circular names the groups gc_content_N and skew_N, Linear gc_content and gc_skew.
  const groups = [...root.querySelectorAll('g[id]')].filter((element) => (
    new RegExp(`^(gc_)?${base}(_\\d+)?$`).test(element.id)
  ));
  return groups.length > 0 && groups.every((element) => element.getAttribute('display') !== 'none');
}, group);

const displayedResult = (page) => page.evaluate(async () => {
  const { getCommittedSvgResultRuntimeIdentity } = await import('/gbdraw/web/js/services/svg-result-ingestion.js');
  const app = window.__GBDRAW_APP__;
  const result = app.results[app.selectedResultIndex] || null;
  return {
    mode: app.mode,
    count: app.results.length,
    identity: result ? getCommittedSvgResultRuntimeIdentity(result) : null,
    content: result?.content ?? null
  };
});

const switchTo = async (page, mode) => {
  await page.getByRole('button', { name: mode === 'circular' ? 'Circular' : 'Linear', exact: true }).click();
  await page.waitForFunction((wanted) => window.__GBDRAW_APP__?.mode === wanted, mode);
  await settleLive(page);
};

for (const [shown, current] of [['circular', 'linear'], ['linear', 'circular']]) {
  for (const [track, toggle] of Object.entries(TRACK_TOGGLES[current])) {
    test(`a ${current} ${track} toggle leaves the ${shown} Result in its slot and applies at the next ${current} Generate`, async ({ page }) => {
      test.setTimeout(240_000);
      await openWithGenBank(page, SINGLE_FIXTURE, () => { window.__GBDRAW_APP__.form.labels_mode = 'out'; });
      await switchTo(page, 'linear');
      await page.evaluate(async (text) => {
        const app = window.__GBDRAW_APP__;
        app.setLinearSeqPrimaryFile(0, 'gb', new File([text], 'forced_label_underlay.gb', { type: 'text/plain', lastModified: 1000 }));
        app.form.show_gc = true;
        app.form.show_skew = true;
        await window.Vue.nextTick();
      }, readFileSync(SINGLE_FIXTURE, 'utf8'));
      await settleLive(page);
      await switchTo(page, shown);
      await generate(page);
      expect(await trackDrawn(page, toggle.group), `${shown} Result draws ${track}`).toBe(true);
      const before = await displayedResult(page);
      await switchTo(page, current);
      expect(await displayedResult(page), `${current} has no Result yet`).toMatchObject({ mode: current, count: 0 });
      await page.evaluate(async ({ field, value }) => {
        window.__GBDRAW_APP__.form[field] = value;
        await window.Vue.nextTick();
      }, { field: toggle.form, value: toggle.toggled });
      await settleLive(page);
      await switchTo(page, shown);
      const after = await displayedResult(page);
      expect(after.identity, `the ${shown} Result after switching back`).toBe(before.identity);
      expect(after.content, `the ${shown} Result's committed SVG`).toBe(before.content);
      expect(await trackDrawn(page, toggle.group)).toBe(true);
      await switchTo(page, current);
      await generate(page);
      const other = Object.values(TRACK_TOGGLES[current]).find((item) => item !== toggle);
      expect(await trackDrawn(page, toggle.group), `the next ${current} Generate hides ${track}`).toBe(false);
      expect(await trackDrawn(page, other.group), `the next ${current} Generate keeps the other track`).toBe(true);
    });
  }
}
