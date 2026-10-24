// PD-OI-066 (R1, R3, R10): changes of track data and mode switches leave the
// live Result equal to the next Generate: Undo and Redo of a Depth source
// removal (OV-66), Legend styles that follow the captions of track data
// (OV-65), a decoration moved on one mode's Result (OV-104), and the track
// toggles of the mode not shown. The matrix of edit kinds is in
// live-generate-parity.playwright.spec.js.
const { expect } = require('@playwright/test');
const { readFileSync } = require('node:fs');
const { writeFile } = require('node:fs/promises');
const { gunzipSync } = require('node:zlib');
const { download } = require('./helpers/mode-transition.cjs');
const { evaluateWithRetainedPromise } = require('./helpers/app-lifecycle.cjs');
const { reveal, warmAppTest: test } = require('./helpers/app-lifecycle.cjs');
const { openWithGenBank } = require('./helpers/audit-browser.cjs');
const { expectLiveEqualsGenerate, semanticSnapshot, settleLive, showResult } = require('./helpers/live-generate-parity.cjs');
const { showDrawing } = require('./helpers/drawing-bar.cjs');
const {
  SINGLE_FIXTURE, open, generate, history, DEPTH_TSV, colorLegendRow, openCanvas, renameRow
} = require('./helpers/live-generate-parity-steps.cjs');

test.describe.configure({ retries: 0 });

// OV-66: Undo of a Depth source removal brings back the Depth track and its
// tick text, and Redo hides them again; the live Result equals what Generate
// draws. A Depth source History step restores files, which suppresses the
// track-visibility watcher, so the step projects the visibility itself.
// OV-1198 (PD-OI-104): the Legend row of the hidden track leaves the Result and
// the Legend list at once, and Undo, Redo, and adding the file again follow.
const BATCH_DEPTH_TSV = ['TESTA', 'TESTB']
  .flatMap((id) => [0, 1, 2, 3].map((index) => `${id}\t${index * 700 + 1}\t${10 + index}`)).join('\n');
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
  'circular, switching Show Depth off': {
    mode: 'circular',
    add: (app, file) => app.setCircularDepthFile(0, file),
    // The Show Depth checkbox's change handler.
    remove: (app) => app.setTrackShown('show_depth', false)
  },
  'circular batch, Result 2 shown, clearing the Depth file': {
    mode: 'circular',
    results: 'batch',
    shown: 1,
    tsv: BATCH_DEPTH_TSV,
    add: (app, file) => app.setCircularDepthFile(0, file),
    remove: (app) => app.setCircularDepthFile(0, null)
  },
  'linear, clearing the Depth file': {
    mode: 'linear',
    add: (app, file) => app.setLinearDepthFile(app.linearSeqs[0], 0, file),
    remove: (app) => app.setLinearDepthFile(app.linearSeqs[0], 0, null)
  }
};
for (const [name, { mode, results = 'single', shown = null, tsv = DEPTH_TSV, add, remove }] of Object.entries(DEPTH_CASES)) {
  test(`Undo and Redo of a Depth source removal match Generate (${name})`, async ({ page }) => {
    test.setTimeout(240_000);
    await open(page, { mode, results, reflow: 'off' });
    const inStep = async (label, change, ...args) => {
      await page.evaluate(async ({ stepLabel, source, values }) => {
        const run = new Function('app', 'text', `return (${source})(app, text && new File([text], 'depth.tsv', { type: 'text/tab-separated-values' }));`);
        await window.__GBDRAW_HISTORY__.runUndoable(stepLabel, () => run(window.__GBDRAW_APP__, values[0]));
      }, { stepLabel: label, source: change.toString(), values: args });
      await settleLive(page);
    };
    const listed = () => page.evaluate(() => window.__GBDRAW_APP__.legendEntries.map((entry) => entry.caption));
    const expectDepthRow = async (step, drawn) => {
      const live = await semanticSnapshot(page);
      expect(live.legend.some((row) => row.caption === 'depth'), `${name}: ${step}: Result row`).toBe(drawn);
      expect((await listed()).includes('depth'), `${name}: ${step}: Legend list row`).toBe(drawn);
    };
    await inStep('Change uploaded file', add, tsv);
    await generate(page);
    if (shown !== null) await showResult(page, shown);
    await expectDepthRow('generated', true);
    await inStep('Remove Depth', remove, null);
    await expectDepthRow('removed', false);
    await history(page, 'undo');
    await expectDepthRow('Undo', true);
    await expectLiveEqualsGenerate(page, { label: `${name}: Undo` });
    // A comparison's Generate drops the Redo branch, so Redo and the re-add
    // follow a new removal on the Result that draws the track.
    await inStep('Remove Depth', remove, null);
    await history(page, 'undo');
    await history(page, 'redo');
    await expectDepthRow('Redo', false);
    await inStep('Change uploaded file', add, tsv);
    await expectDepthRow('added again', true);
    await expectLiveEqualsGenerate(page, { label: `${name}: added again` });
    await inStep('Remove Depth', remove, null);
    await expectDepthRow('removed after Generate', false);
    await expectLiveEqualsGenerate(page, { label: `${name}: Remove` });
  });
}

// OV-1198, OV-1215 (PD-OI-104): a drawn track hidden or shown again through
// its checkbox (one History step): its Legend rows leave and return on the
// Result and in the Legend list at once, and Undo, Redo and Generate agree.
// The names the Legend list shows (its name inputs), in list order.
const shownListNames = (page, pattern) => page.evaluate((source) => [...document.querySelectorAll('input[aria-label^="Legend entry "][aria-label$=" name"]')]
  .map((input) => input.value).filter((value) => new RegExp(source, 'i').test(value)), pattern.source);
const clickTrackCheckbox = (page, name) => page.evaluate((label) => {
  const input = [...document.querySelectorAll('input[type=checkbox]')].find((each) => (
    each.getAttribute('aria-label') === label || each.closest('label')?.textContent.trim() === label));
  if (!input || input.disabled) throw new Error(`no enabled checkbox "${label}"`);
  input.click();
}, name);
const twoDepthSeries = (app, text) => {
  app.setCircularDepthFile(0, new File([text], 'depth.tsv', { type: 'text/tab-separated-values' }));
  app.addCircularDepthTrack();
  app.setCircularDepthFile(1, new File([text], 'depth.tsv', { type: 'text/tab-separated-values' }));
};
const oneDepthSeries = (app, text) => {
  app.setCircularDepthFile(0, new File([text], 'depth.tsv', { type: 'text/tab-separated-values' }));
};
// The drawn groups of a track kind (drawn-groups.js) the displayed Result shows.
const shownTrackGroups = (page, groups) => page.evaluate(async (name) => {
  const drawnGroups = await import('./js/services/drawn-groups.js');
  const svg = window.__GBDRAW_APP__.svgContainer?.querySelector('svg');
  return drawnGroups.findTrackGroups(svg, drawnGroups[name]).filter((group) => group.getAttribute('display') !== 'none').length;
}, groups);
const HIDDEN_TRACK_CASES = {
  'circular Hide GC Content': { mode: 'circular', box: 'Hide GC Content', rows: /^GC content/i, groups: 'DINUCLEOTIDE_CONTENT_GROUPS' },
  'circular Hide GC Skew': { mode: 'circular', box: 'Hide GC Skew', rows: /^GC skew/i, count: 2, groups: 'DINUCLEOTIDE_SKEW_GROUPS' },
  'linear Show GC Content': {
    mode: 'linear', box: 'Show GC Content', rows: /^GC content/i, form: { show_gc: true }, groups: 'DINUCLEOTIDE_CONTENT_GROUPS'
  },
  'linear Show GC Skew': {
    mode: 'linear', box: 'Show GC Skew', rows: /^GC skew/i, count: 2, form: { show_skew: true }, groups: 'DINUCLEOTIDE_SKEW_GROUPS'
  },
  'circular Show Depth, two series labelled depth': {
    mode: 'circular', box: 'Show Depth', rows: /^depth( \([12]\))?$/, count: 2, data: twoDepthSeries, groups: 'DEPTH_GROUPS'
  },
  // R2-M2: a batch Result selected after the hide shows neither the track nor its row.
  'circular batch Show Depth, Result 2 selected after the hide': {
    mode: 'circular', results: 'batch', box: 'Show Depth', rows: /^depth$/, data: oneDepthSeries, tsv: BATCH_DEPTH_TSV, groups: 'DEPTH_GROUPS'
  }
};
for (const [name, { mode, results = 'single', box, rows, count = 1, form = {}, data = null, tsv = DEPTH_TSV, groups }]
  of Object.entries(HIDDEN_TRACK_CASES)) {
  test(`A hidden drawn track's Legend rows leave and return with it (${name})`, async ({ page }) => {
    test.setTimeout(240_000);
    await open(page, { mode, results, reflow: 'off' });
    await page.evaluate(({ values, source, text }) => {
      const app = window.__GBDRAW_APP__;
      Object.assign(app.form, values);
      if (source) new Function('app', 'text', `return (${source})(app, text);`)(app, text);
    }, { values: form, source: data && data.toString(), text: tsv });
    await settleLive(page);
    await generate(page);
    const expectRows = async (step, drawn) => {
      const live = (await semanticSnapshot(page)).legend.filter((row) => rows.test(row.caption)).map((row) => row.caption);
      expect(live.length, `${name}: ${step}: Result rows`).toBe(drawn ? count : 0);
      // OV-1217 (PD-OI-095): the list names the rows as the Result draws them.
      expect(await shownListNames(page, rows), `${name}: ${step}: Legend list names`).toEqual(live);
      expect((await shownTrackGroups(page, groups)) > 0, `${name}: ${step}: track shown`).toBe(drawn);
    };
    const toggle = async () => { await clickTrackCheckbox(page, box); await settleLive(page); };
    await expectRows('generated', true);
    await toggle();
    await expectRows('hidden', false);
    if (results === 'batch') {
      await showResult(page, 1);
      await expectRows('hidden, Result 2 selected', false);
    }
    await history(page, 'undo');
    await expectRows('Undo', true);
    await expectLiveEqualsGenerate(page, { label: `${name}: Undo` });
    await toggle();
    await history(page, 'undo');
    await history(page, 'redo');
    await expectRows('Redo', false);
    await toggle();
    await expectRows('shown again', true);
    await expectLiveEqualsGenerate(page, { label: `${name}: shown again` });
    // A returned row takes Python's place, in the list as on the Result.
    const order = await page.evaluate(() => window.__GBDRAW_APP__.legendEntries.map((entry) => entry.originalCaption));
    expect(order, `${name}: shown again: Legend list order`).toEqual((await semanticSnapshot(page)).legend.map((row) => row.key));
    await toggle();
    await expectRows('hidden after Generate', false);
    await expectLiveEqualsGenerate(page, { label: `${name}: hidden` });
  });
}

// REVIEW H1 (OV-1198): a renamed Depth row hidden and shown again returns once,
// under its name, and leaves no waiting copy; renamed back, it keeps that name.
test('A renamed Depth row hidden and shown again returns once under its name', async ({ page }) => {
  test.setTimeout(240_000);
  await open(page, { mode: 'circular', results: 'single', reflow: 'off' });
  await drawDepthFile()(page);
  const rowsNow = () => page.evaluate(async () => {
    const { state } = await import('/gbdraw/web/js/state.js');
    const drawing = state.activeDrawing();
    const unref = (value) => (value && 'value' in value ? value.value : value) || [];
    return {
      listed: unref(drawing.legendEntries).map((entry) => entry.caption).filter((caption) => /depth|coverage/i.test(caption)),
      dormant: unref(drawing.dormantLegendEntries).length
    };
  });
  const hideAndShow = async () => {
    await clickTrackCheckbox(page, 'Show Depth');
    await settleLive(page);
    await clickTrackCheckbox(page, 'Show Depth');
    await settleLive(page);
  };
  await renameRow(page, 'depth', 'Coverage');
  await settleLive(page);
  await hideAndShow();
  expect(await rowsNow(), 'renamed, hidden, shown').toEqual({ listed: ['Coverage'], dormant: 0 });
  await renameRow(page, 'Coverage', 'depth');
  await settleLive(page);
  await hideAndShow();
  expect(await rowsNow(), 'renamed back, hidden, shown').toEqual({ listed: ['depth'], dormant: 0 });
  await expectLiveEqualsGenerate(page, { label: 'renamed back' });
});

// Brief 7 (OV-1198, OV-1216): the Depth row across hide/show, Session Save and Load, Reset, a
// second series, and a label equal to a feature type.
const depthRows = (page) => page.evaluate(() => {
  const app = window.__GBDRAW_APP__;
  return app.legendEntries.map((entry) => entry.caption).filter((caption) => /depth|coverage/i.test(caption));
});
const resultDepthRows = async (page) => (await semanticSnapshot(page)).legend
  .map((row) => row.caption).filter((caption) => /depth|coverage/i.test(caption));
const toggleDepth = async (page) => { await clickTrackCheckbox(page, 'Show Depth'); await settleLive(page); };
const saveSessionJson = async (page, name) => {
  const path = test.info().outputPath(`${name}-${Date.now()}.gbdraw-session.json`);
  await page.evaluate(async (title) => {
    const { state } = await import('./js/state.js');
    state.sessionTitle.value = title;
  }, name);
  const bytes = await download(page, 'Save Session', path);
  return JSON.parse((bytes[0] === 0x1f ? gunzipSync(bytes) : bytes).toString('utf8'));
};
// Counts of a caption among the saved Legend entries and every object key of the Session.
const sessionKeyCount = (saved, key) => {
  let count = 0;
  const walk = (value) => {
    if (Array.isArray(value)) value.forEach(walk);
    else if (value && typeof value === 'object') {
      for (const [name, each] of Object.entries(value)) {
        if (name === key) count += 1;
        walk(each);
      }
    }
  };
  walk(saved.modes);
  return count;
};

test('A Depth row renamed to Coverage, hidden, shown, saves its key once', async ({ page }) => {
  test.setTimeout(240_000);
  await open(page, { mode: 'circular', results: 'single', reflow: 'off' });
  await drawDepthFile()(page);
  await renameRow(page, 'depth', 'Coverage');
  await settleLive(page);
  await toggleDepth(page);
  await toggleDepth(page);
  expect(await depthRows(page), 'listed').toEqual(['Coverage']);
  const saved = await saveSessionJson(page, 'depth-key-once');
  const entries = saved.modes.circular.editorState.legend.entries
    .filter((entry) => /depth|coverage/i.test(`${entry.originalCaption} ${entry.caption}`));
  expect(entries.map((entry) => `${entry.originalCaption}=>${entry.caption}`), 'saved Legend entries').toEqual(['depth=>Coverage']);
  const text = JSON.stringify(saved.modes.circular.editorState);
  expect((text.match(/"depth"/g) || []).length, 'the key "depth" in the saved editor state').toBeLessThanOrEqual(1);
  expect(sessionKeyCount(saved, 'Coverage'), 'a "Coverage" object key').toBeLessThanOrEqual(1);
});

test('Two Depth series labelled depth draw depth (1) and depth (2), and Show Depth off hides both', async ({ page }) => {
  test.setTimeout(240_000);
  await open(page, { mode: 'circular', results: 'single', reflow: 'off' });
  await page.evaluate(({ text }) => {
    const app = window.__GBDRAW_APP__;
    const file = () => new File([text], 'depth.tsv', { type: 'text/tab-separated-values' });
    app.setCircularDepthFile(0, file());
    app.addCircularDepthTrack();
    app.setCircularDepthFile(1, file());
  }, { text: DEPTH_TSV });
  await settleLive(page);
  const request = await page.evaluate(async () => {
    const { getCommittedCanonicalRenderRequest } = await import('./js/services/config.js');
    const committed = getCommittedCanonicalRenderRequest();
    return JSON.stringify({
      depthTracks: committed.depthTracks ?? committed.tracks?.depth ?? null,
      slots: (committed.trackSlots ?? committed.circularTrackSlots ?? []).map((slot) => ({ r: slot.renderer, e: slot.enabled }))
    });
  });
  console.log(`OV-1217 committed request: ${request}`);
  await generate(page);
  // Python numbers equal captions `depth (1)`, `depth (2)`; the row set must equal Generate's.
  expect(await resultDepthRows(page), 'Generate rows').toEqual(['depth (1)', 'depth (2)']);
  // OV-1217 (PD-OI-095): the Legend list shows the names the Result draws; a row's key stays `depth`.
  expect(await shownListNames(page, /depth/i), 'Legend list names after Generate').toEqual(['depth (1)', 'depth (2)']);
  expect(await depthRows(page), 'Legend list keys after Generate').toEqual(['depth', 'depth (2)']);
  await toggleDepth(page);
  expect(await resultDepthRows(page), 'hidden: Result rows').toEqual([]);
  expect(await depthRows(page), 'hidden: listed rows').toEqual([]);
  await expectLiveEqualsGenerate(page, { label: 'two Depth series hidden' });
});

// R2-H1 (OV-1198): while the Linear custom stack is on, its enabled slots draw the tracks; the Show Depth
// and Show GC Content checkboxes it disables keep a stale false. Generate draws the track and its row,
// and the Result and the Legend list keep both (PD-OI-066).
const CUSTOM_STACK_CASES = {
  'an enabled Depth slot, Show Depth off': { depth: true, flag: 'show_depth', renderer: 'depth', rows: /^depth$/, groups: 'DEPTH_GROUPS' },
  'a GC content slot, Show GC Content off': {
    depth: false, flag: 'show_gc', renderer: 'dinucleotide_content', rows: /^GC content$/, groups: 'DINUCLEOTIDE_CONTENT_GROUPS'
  }
};
for (const [name, { depth, flag, renderer, rows, groups }] of Object.entries(CUSTOM_STACK_CASES)) {
  test(`A Linear custom stack draws its track and row from Generate (${name})`, async ({ page }) => {
    test.setTimeout(240_000);
    await open(page, { mode: 'linear', results: 'single', reflow: 'off' });
    if (depth) await drawDepthFile('linear')(page);
    // The checkbox's change handler, before the stack disables the checkbox.
    await page.evaluate((key) => window.__GBDRAW_APP__.setTrackShown(key, false), flag);
    await settleLive(page);
    const draft = await page.evaluate(({ key, kind }) => {
      const app = window.__GBDRAW_APP__;
      app.setLinearTrackSlotsEnabled(true);
      if (!app.adv.linear_track_slots.some((slot) => slot.renderer === kind && slot.enabled !== false)) app.addLinearTrackSlot(kind);
      return { flag: app.form[key], slots: app.adv.linear_track_slots.filter((slot) => slot.renderer === kind && slot.enabled !== false).length };
    }, { key: flag, kind: renderer });
    expect(draft, `${name}: the draft`).toEqual({ flag: false, slots: 1 });
    await settleLive(page);
    expect(await page.evaluate((key) => window.__GBDRAW_APP__.form[key], flag), `${name}: ${flag} before Generate`).toBe(false);
    await generate(page);
    const live = (await semanticSnapshot(page)).legend.filter((row) => rows.test(row.caption)).map((row) => row.caption);
    expect(live, `${name}: Result rows`).toHaveLength(1);
    expect(await shownListNames(page, rows), `${name}: Legend list names`).toEqual(live);
    expect(await shownTrackGroups(page, groups), `${name}: track shown`).toBeGreaterThan(0);
  });
}

// R3-H-A (OV-1198): a slot edit that stops a track being drawn hides the track and its Legend rows in
// one History step, as the track checkboxes do; Undo returns both (PD-OI-104).
const SLOT_EDIT_CASES = {
  'circular, Depth slot switched off': {
    mode: 'circular', rows: /^depth$/, groups: 'DEPTH_GROUPS',
    setup: (app, text) => {
      app.setCircularDepthFile(0, new File([text], 'depth.tsv', { type: 'text/tab-separated-values' }));
      app.setTrackShown('show_depth', true);
      app.setCircularTrackSlotsEnabled(true);
      return app.adv.circular_track_slots.filter((slot) => slot.renderer === 'depth' && slot.enabled !== false).length;
    },
    edit: (app) => app.setCircularTrackSlotEnabled(app.adv.circular_track_slots.find((slot) => slot.renderer === 'depth'), false)
  },
  'linear, custom stack switched off with Show GC Content off': {
    mode: 'linear', rows: /^GC content$/, groups: 'DINUCLEOTIDE_CONTENT_GROUPS',
    setup: (app) => {
      app.setTrackShown('show_gc', false);
      app.setLinearTrackSlotsEnabled(true);
      if (!app.adv.linear_track_slots.some((slot) => slot.renderer === 'dinucleotide_content' && slot.enabled !== false)) {
        app.addLinearTrackSlot('dinucleotide_content');
      }
      return app.adv.linear_track_slots.filter((slot) => slot.renderer === 'dinucleotide_content' && slot.enabled !== false).length;
    },
    edit: (app) => app.setLinearTrackSlotsEnabled(false)
  }
};
for (const [name, { mode, rows, groups, setup, edit }] of Object.entries(SLOT_EDIT_CASES)) {
  test(`A slot edit that stops drawing a track hides it and its Legend rows at once (${name})`, async ({ page }) => {
    test.setTimeout(240_000);
    await open(page, { mode, results: 'single', reflow: 'off' });
    const slots = await page.evaluate(({ source, text }) => new Function('app', 'text', `return (${source})(app, text);`)(
      window.__GBDRAW_APP__, text), { source: setup.toString(), text: DEPTH_TSV });
    expect(slots, `${name}: enabled slots`).toBe(1);
    await settleLive(page);
    await generate(page);
    const expectRows = async (step, drawn) => {
      const live = (await semanticSnapshot(page)).legend.filter((row) => rows.test(row.caption)).map((row) => row.caption);
      expect(live.length, `${name}: ${step}: Result rows`).toBe(drawn ? 1 : 0);
      expect(await shownListNames(page, rows), `${name}: ${step}: Legend list names`).toEqual(live);
      expect((await shownTrackGroups(page, groups)) > 0, `${name}: ${step}: track shown`).toBe(drawn);
    };
    await expectRows('generated', true);
    await inHistoryStep(page, 'Edit track slot', edit);
    await expectRows('slot edit', false);
    // One Undo returns the track and its row; the parity Generate comes last, as it is a History step of its own.
    await history(page, 'undo');
    await expectRows('Undo', true);
    await history(page, 'redo');
    await expectRows('Redo', false);
    await expectLiveEqualsGenerate(page, { label: `${name}: slot edit` });
  });
}

// R3-M-A (OV-1198): a Circular custom stack without an enabled Depth slot draws no Depth whatever the
// Show Depth checkbox it disables says, so a Depth series without a file does not stop Generate (the build used its own copy of the rule).
test('A Circular custom stack without a Depth slot generates with a Depth series that has no file', async ({ page }) => {
  test.setTimeout(240_000);
  await open(page, { mode: 'circular', results: 'single', reflow: 'off' });
  const draft = await page.evaluate((text) => {
    const app = window.__GBDRAW_APP__;
    // Series 1 has no file; series 2 has one.
    app.addCircularDepthTrack();
    app.setCircularDepthFile(1, new File([text], 'depth.tsv', { type: 'text/tab-separated-values' }));
    app.setTrackShown('show_depth', true);
    app.setCircularTrackSlotsEnabled(true);
    app.adv.circular_track_slots.filter((slot) => slot.renderer === 'depth').forEach((slot) => { slot.enabled = false; });
    return { showDepth: app.form.show_depth, depthSlots: app.adv.circular_track_slots.filter((slot) => slot.renderer === 'depth' && slot.enabled !== false).length };
  }, DEPTH_TSV);
  expect(draft, 'the draft').toEqual({ showDepth: true, depthSlots: 0 });
  await settleLive(page);
  await generate(page);
  expect(await page.evaluate(() => window.__GBDRAW_APP__.errorLog ?? null), 'Generate error').toBeFalsy();
  expect(await shownTrackGroups(page, 'DEPTH_GROUPS'), 'shown Depth groups').toBe(0);
});

// Reset Settings keeps the Depth file and turns Show Depth off, so the next Generate draws no Depth
// row. The list must equal that Generate, and the Result must equal the list at once, also after the
// automatic rerender the Reset requests, which draws the committed request again (OV-1198 H2).
const resetWithDepthKept = async (page) => {
  await open(page, { mode: 'circular', results: 'single', reflow: 'off' });
  await drawDepthFile()(page);
  await renameRow(page, 'depth', 'Coverage');
  await settleLive(page);
  await evaluateWithRetainedPromise(page, async () => { await window.__GBDRAW_APP__.resetSettings(); });
  await settleLive(page);
};
test('Reset Settings with the Depth files kept lists what the next Generate draws', async ({ page }) => {
  test.setTimeout(240_000);
  await resetWithDepthKept(page);
  const listed = await depthRows(page);
  await generate(page);
  expect(await depthRows(page), 'listed after the next Generate').toEqual(listed);
  expect(await resultDepthRows(page), 'Result after the next Generate').toEqual(listed);
});
test('Reset Settings with the Depth files kept shows on the Result what the list holds', async ({ page }) => {
  test.setTimeout(240_000);
  await resetWithDepthKept(page);
  expect(await resultDepthRows(page), 'Result after Reset').toEqual(await depthRows(page));
  expect(await depthRows(page), 'listed after Reset').toEqual([]);
  expect(await shownTrackGroups(page, 'DEPTH_GROUPS'), 'shown Depth groups').toBe(0);
});
// R3-H-B (OV-1198): the rerender the Reset requests hides the track the draft hides also when the
// track claims no Legend row: a Depth series labelled CDS in the CDS fill shares the CDS row (D-15-65).
test('Reset Settings with the Depth files kept hides a Depth track that claims no Legend row', async ({ page }) => {
  test.setTimeout(240_000);
  await open(page, { mode: 'circular', results: 'single', reflow: 'off' });
  const [cds] = (await semanticSnapshot(page)).legend.filter((row) => row.caption === 'CDS');
  await page.evaluate(({ text, fill }) => {
    const app = window.__GBDRAW_APP__;
    app.setCircularDepthFile(0, new File([text], 'depth.tsv', { type: 'text/tab-separated-values' }));
    Object.assign(app.adv.depth_tracks[0], { label: 'CDS', color: fill });
    app.adv.block_stroke_width = 0;
  }, { text: DEPTH_TSV, fill: cds.fill });
  await settleLive(page);
  await generate(page);
  expect(await shownTrackGroups(page, 'DEPTH_GROUPS'), 'shown Depth groups after Generate').toBeGreaterThan(0);
  // The renamed row makes the Reset ask for the rerender, as in the cases above.
  await renameRow(page, 'CDS', 'Coding');
  await settleLive(page);
  await evaluateWithRetainedPromise(page, async () => { await window.__GBDRAW_APP__.resetSettings(); });
  await settleLive(page);
  expect(await page.evaluate(() => window.__GBDRAW_APP__.form.show_depth), 'Show Depth after Reset').toBe(false);
  expect(await shownTrackGroups(page, 'DEPTH_GROUPS'), 'shown Depth groups after Reset').toBe(0);
});

test('A Session saved with Depth hidden loads without a Depth row in the list', async ({ page }) => {
  test.setTimeout(240_000);
  await open(page, { mode: 'circular', results: 'single', reflow: 'off' });
  await drawDepthFile()(page);
  await toggleDepth(page);
  expect(await depthRows(page), 'hidden before Save').toEqual([]);
  const saved = await saveSessionJson(page, 'depth-hidden');
  const path = test.info().outputPath(`depth-hidden-reload-${Date.now()}.gbdraw-session.json`);
  await writeFile(path, JSON.stringify(saved));
  await page.locator('input[accept^=".json,"]').setInputFiles(path);
  await expect.poll(() => page.evaluate(() => !window.__GBDRAW_APP__.sessionImportPending), { timeout: 180000 }).toBe(true);
  await settleLive(page);
  expect(await depthRows(page), 'listed after Load').toEqual([]);
  expect(await resultDepthRows(page), 'Result after Load').toEqual([]);
  await generate(page);
  expect(await depthRows(page), 'listed after Generate').toEqual([]);
});

// Residual risk 2: a Depth label equal to a feature-type palette key. Generate's colour and the live
// colour are recorded; the case pins that they agree.
test('A Depth label equal to a feature type draws the same colour live and from Generate', async ({ page }) => {
  test.setTimeout(240_000);
  await open(page, { mode: 'circular', results: 'single', reflow: 'off' });
  await drawDepthFile()(page);
  await renameRow(page, 'depth', 'tRNA');
  await settleLive(page);
  const swatch = async () => (await semanticSnapshot(page)).legend.filter((row) => row.caption === 'tRNA');
  const live = await swatch();
  await generate(page);
  const generated = await swatch();
  console.log(`CASE6 live=${JSON.stringify(live)} generated=${JSON.stringify(generated)}`);
  expect(live, 'live rows named tRNA').toEqual(generated);
});

// D-15-65 (OV-1198 review 2 M3): a Depth series labelled like a feature row, in
// that row's fill, shares the feature row and draws none of its own, so Show
// Depth off keeps the row on the Result and in the list, as Generate draws it.
test('A Depth series labelled CDS in the CDS fill shares the CDS row, and Show Depth off keeps it', async ({ page }) => {
  test.setTimeout(240_000);
  await open(page, { mode: 'circular', results: 'single', reflow: 'off' });
  const cdsRows = async () => (await semanticSnapshot(page)).legend.filter((row) => /^CDS( \(\d+\))?$/.test(row.caption));
  const [cds] = await cdsRows();
  await page.evaluate(({ text, fill }) => {
    const app = window.__GBDRAW_APP__;
    app.setCircularDepthFile(0, new File([text], 'depth.tsv', { type: 'text/tab-separated-values' }));
    Object.assign(app.adv.depth_tracks[0], { label: 'CDS', color: fill });
    // Python shares a row of the same name and swatch, stroke included: the feature blocks drop theirs.
    app.adv.block_stroke_width = 0;
  }, { text: DEPTH_TSV, fill: cds.fill });
  await settleLive(page);
  await generate(page);
  const expectCdsRow = async (step, depthShown) => {
    expect((await cdsRows()).map((row) => [row.caption, row.fill]), `${step}: Result rows`).toEqual([['CDS', cds.fill]]);
    expect(await shownListNames(page, /^CDS( \(\d+\))?$/), `${step}: Legend list names`).toEqual(['CDS']);
    expect((await shownTrackGroups(page, 'DEPTH_GROUPS')) > 0, `${step}: Depth shown`).toBe(depthShown);
  };
  await expectCdsRow('generated', true);
  await clickTrackCheckbox(page, 'Show Depth');
  await settleLive(page);
  await expectCdsRow('Show Depth off', false);
  await expectLiveEqualsGenerate(page, { label: 'Show Depth off' });
});

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
// style, and the live Result equals Generate. Generate after the change
// succeeds and draws no row of the data. A row without a rename runs both in
// one case: the change is made again after the Undo half, on data Generate
// drew. A row renamed in the Legend (OV-87) runs Generate right after the
// first change, on the setup's live edits with no Generate in between, in a
// case of its own; `generateAfter` draws the renamed row once first.
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
    const openWithRow = async (page) => {
      await openCanvas(page, mode, null);
      await setup(page);
      const row = await caption(page);
      expect(row, 'row caption').toBeTruthy();
      return row;
    };
    const changeRetires = async (page, row) => {
      await change(page);
      await settleLive(page);
      expect(await legendStyleOf(page, row)).toEqual(NO_STYLE);
      if (renamedFrom) {
        expect(await legendStyleOf(page, renamedFrom)).toEqual(NO_STYLE);
        expect(await legendNames(page), 'the rename is retired').not.toContain(`${renamedFrom}=>${row}`);
      }
    };
    const undoRestores = async (page, row) => {
      await undo(page);
      expect(await restored(page), 'data restored').toBe(true);
      // A renamed row's style is stored under its key (PR 3, D-25).
      expect((await legendStyleOf(page, renamedFrom || row)).color).toBe(COLORED);
      if (renamedFrom) expect(await legendNames(page), 'the rename is restored').toContain(`${renamedFrom}=>${row}`);
      await expectLiveEqualsGenerate(page, { label: `${name}, Undo` });
    };
    const generateDrawsNoRow = async (page, row) => {
      await generate(page);
      const captions = await drawnLegendCaptions(page);
      expect(captions).not.toContain(row);
      if (renamedFrom) expect(captions).not.toContain(renamedFrom);
      for (const kept of drawn) expect(captions).toContain(kept);
    };
    if (renamedFrom) {
      test(`${name} retires the styles of its rows, and Generate succeeds`, async ({ page }) => {
        const row = await openWithRow(page);
        await changeRetires(page, row);
        await generateDrawsNoRow(page, row);
      });
      test(`${name}: Undo restores the data and the style, and live equals Generate`, async ({ page }) => {
        const row = await openWithRow(page);
        await changeRetires(page, row);
        await undoRestores(page, row);
      });
      continue;
    }
    test(`${name} retires the styles of its rows, Undo restores them as Generate draws, and Generate succeeds`, async ({ page }) => {
      const row = await openWithRow(page);
      await changeRetires(page, row);
      await undoRestores(page, row);
      await changeRetires(page, row);
      await generateDrawsNoRow(page, row);
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
    expect((await legendStyleOf(page, 'depth')).color, 'the style stays under the row key').toBe(COLORED);
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
    const field = await reveal(page.getByRole('group', { name: 'Circular track slot gc_content', exact: true })
      .getByRole('textbox', { name: 'Track legend label', includeHidden: true }));
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
    await showDrawing(page, { kind: target });
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
  await showDrawing(page, { kind: mode });
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
