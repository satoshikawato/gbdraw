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
const { download, load } = require('./helpers/mode-transition.cjs');
const { evaluateWithRetainedPromise } = require('./helpers/app-lifecycle.cjs');
const { reveal, warmAppTest: test } = require('./helpers/app-lifecycle.cjs');
const { openWithGenBank } = require('./helpers/audit-browser.cjs');
const { expectLiveEqualsGenerate, semanticSnapshot, settleLive, showResult } = require('./helpers/live-generate-parity.cjs');
const {
  SINGLE_FIXTURE, open, generate, history, DEPTH_TSV, colorLegendRow, openCanvas, renameRow
} = require('./helpers/live-generate-parity-steps.cjs');

test.describe.configure({ retries: 0 });

// OV-66: Undo of a Depth source removal brings back the Depth track and its
// tick text, and Redo hides them again; the live Result equals what Generate
// draws. A Depth source History step restores files, which suppresses the
// track-visibility watcher, so the step projects the visibility itself.
// OV-1198 (PD-OI-104): the Legend row of the hidden track leaves the Result and
// the Legend list at once, and Undo and Redo follow at once. Adding the file
// again after a removal shows on Generate (PD-OI-104 rev 3).
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
    // PD-OI-104 rev 3: a Depth file added after a removal shows on Generate.
    await generate(page);
    if (shown !== null) await showResult(page, shown);
    await expectDepthRow('added again, after Generate', true);
    await expectLiveEqualsGenerate(page, { label: `${name}: added again, after Generate` });
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
    // The hide is a live redraw from the kept track parts (PD-OI-104 rev 5): wait for it.
    await settleLive(page);
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
      identity: getCommittedSvgResultRuntimeIdentity(window.__GBDRAW_APP__.currentResults()[state.selectedResultIndex.value]),
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
  const root = new DOMParser().parseFromString(app.currentResults()[app.selectedResultIndex].content, 'image/svg+xml').documentElement;
  // Circular names the groups gc_content_N and skew_N, Linear gc_content and gc_skew.
  const groups = [...root.querySelectorAll('g[id]')].filter((element) => (
    new RegExp(`^(gc_)?${base}(_\\d+)?$`).test(element.id)
  ));
  return groups.length > 0 && groups.every((element) => element.getAttribute('display') !== 'none');
}, group);

const displayedResult = (page) => page.evaluate(async () => {
  const { getCommittedSvgResultRuntimeIdentity } = await import('/gbdraw/web/js/services/svg-result-ingestion.js');
  const app = window.__GBDRAW_APP__;
  const result = app.currentResults()[app.selectedResultIndex] || null;
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

// PD-OI-104 rev 5 (the GT set): Show GC Content, Show GC Skew, Show Depth, Hide GC Content,
// Hide GC Skew, a GC, skew or Depth row's switch, and Use custom stack apply live, also for a
// track the Result does not draw yet, from the parts the last Generate built (D-199). Each case
// ends with the live Result equal to the next Generate (R3).
const depthFileOf = (text) => new File([text], 'depth.tsv', { type: 'text/tab-separated-values' });
const withTwoDepthSeries = (mode) => async (page) => page.evaluate(({ text, linear }) => {
  const app = window.__GBDRAW_APP__;
  const file = () => new File([text], 'depth.tsv', { type: 'text/tab-separated-values' });
  if (linear) {
    app.setLinearDepthFile(app.linearSeqs[0], 0, file());
    app.addLinearDepthTrack();
    app.setLinearDepthFile(app.linearSeqs[0], 1, file());
    app.setTrackShown('show_gc', true);
    app.setLinearTrackSlotsEnabled(true);
  } else {
    app.setCircularDepthFile(0, file());
    app.addCircularDepthTrack();
    app.setCircularDepthFile(1, file());
    app.setCircularTrackSlotsEnabled(true);
  }
}, { text: DEPTH_TSV, linear: mode === 'linear' });
const slotSwitch = (page, mode, renderer, nth, enabled) => inHistoryStep(page, 'Edit track slot', (app, value) => {
  const slots = value.linear ? app.adv.linear_track_slots : app.adv.circular_track_slots;
  const row = slots.filter((slot) => slot.renderer === value.renderer)[value.nth];
  if (!row) throw new Error(`no ${value.renderer} row ${value.nth}`);
  if (value.linear) app.setLinearTrackSlotEnabled(row, value.enabled);
  else app.setCircularTrackSlotEnabled(row, value.enabled);
}, { linear: mode === 'linear', renderer, nth, enabled });
// In a Circular custom stack, turning a Hide control on asks first; the test chooses to hide.
const hideInStack = async (page, box, choice) => {
  await clickTrackCheckbox(page, box);
  await page.getByRole('button', { name: choice, exact: true }).click();
  await settleLive(page);
};
// The text of the track toggles' scope (the Layout card may be closed).
const trackScopeText = (page) => page.evaluate(() => [...document.querySelectorAll('[data-setting="track-visibility"]')]
  .map((scope) => scope.textContent).join(' | '));

for (const mode of ['circular', 'linear']) {
  test(`GT B2: one of two Depth rows switched off and on applies live (${mode} stack)`, async ({ page }) => {
    test.setTimeout(300_000);
    await open(page, { mode, results: 'single', reflow: 'off' });
    await withTwoDepthSeries(mode)(page);
    await settleLive(page);
    await generate(page);
    const both = await shownTrackGroups(page, 'DEPTH_GROUPS');
    await slotSwitch(page, mode, 'depth', 1, false);
    await settleLive(page);
    expect(await shownTrackGroups(page, 'DEPTH_GROUPS'), `${mode}: one Depth row off`).toBeLessThan(both);
    await expectLiveEqualsGenerate(page, { label: `${mode}: second Depth row off` });
    // The row was off at that Generate: it turns on from the kept parts.
    await slotSwitch(page, mode, 'depth', 1, true);
    await settleLive(page);
    expect(await shownTrackGroups(page, 'DEPTH_GROUPS'), `${mode}: the row returns live`).toBe(both);
    await expectLiveEqualsGenerate(page, { label: `${mode}: second Depth row on` });
  });

  test(`GT B3/B4: a GC content row's switch hides and shows live (${mode} stack)`, async ({ page }) => {
    test.setTimeout(300_000);
    await open(page, { mode, results: 'single', reflow: 'off' });
    await page.evaluate((linear) => {
      const app = window.__GBDRAW_APP__;
      if (linear) {
        app.setTrackShown('show_gc', true);
        app.setLinearTrackSlotsEnabled(true);
        if (!app.adv.linear_track_slots.some((slot) => slot.renderer === 'dinucleotide_content')) app.addLinearTrackSlot('dinucleotide_content');
      } else app.setCircularTrackSlotsEnabled(true);
    }, mode === 'linear');
    await settleLive(page);
    await generate(page);
    await slotSwitch(page, mode, 'dinucleotide_content', 0, false);
    await settleLive(page);
    expect(await shownTrackGroups(page, 'DINUCLEOTIDE_CONTENT_GROUPS'), `${mode}: GC row off`).toBe(0);
    await expectLiveEqualsGenerate(page, { label: `${mode}: GC row off` });
    await slotSwitch(page, mode, 'dinucleotide_content', 0, true);
    await settleLive(page);
    expect(await shownTrackGroups(page, 'DINUCLEOTIDE_CONTENT_GROUPS'), `${mode}: GC row on`).toBeGreaterThan(0);
    await expectLiveEqualsGenerate(page, { label: `${mode}: GC row on` });
  });

  test(`GT B6/B7: Show Depth on after a Generate with it off shows live, through History, a mode switch, Save and Load (${mode})`, async ({ page }) => {
    test.setTimeout(400_000);
    await open(page, { mode, results: 'single', reflow: 'off' });
    await addDepthFile(page, 'depth.tsv', mode);
    await page.evaluate(() => window.__GBDRAW_APP__.setTrackShown('show_depth', false));
    await settleLive(page);
    await generate(page);
    expect(await shownTrackGroups(page, 'DEPTH_GROUPS'), `${mode}: generated without Depth`).toBe(0);
    await toggleDepth(page);
    expect(await shownTrackGroups(page, 'DEPTH_GROUPS'), `${mode}: shown live`).toBeGreaterThan(0);
    expect(await depthRows(page), `${mode}: Legend list`).toEqual(await resultDepthRows(page));
    await history(page, 'undo');
    expect(await shownTrackGroups(page, 'DEPTH_GROUPS'), `${mode}: Undo hides`).toBe(0);
    await history(page, 'redo');
    expect(await shownTrackGroups(page, 'DEPTH_GROUPS'), `${mode}: Redo shows`).toBeGreaterThan(0);
    const other = mode === 'circular' ? 'linear' : 'circular';
    await switchTo(page, other);
    await switchTo(page, mode);
    expect(await shownTrackGroups(page, 'DEPTH_GROUPS'), `${mode}: after a mode switch and back`).toBeGreaterThan(0);
    const saved = await saveSessionJson(page, `gt-depth-shown-${mode}`);
    const path = test.info().outputPath(`gt-depth-shown-${mode}-${Date.now()}.gbdraw-session.json`);
    await writeFile(path, JSON.stringify(saved));
    await page.locator('input[accept^=".json,"]').setInputFiles(path);
    await expect.poll(() => page.evaluate(() => !window.__GBDRAW_APP__.sessionImportPending), { timeout: 180000 }).toBe(true);
    await settleLive(page);
    expect(await resultDepthRows(page), `${mode}: Result after Load`).toEqual(await depthRows(page));
    expect(await shownTrackGroups(page, 'DEPTH_GROUPS'), `${mode}: Depth after Load`).toBeGreaterThan(0);
    await expectLiveEqualsGenerate(page, { label: `${mode}: Depth shown live, saved and loaded` });
  });
}

// A Session a Save wrote names its Depth files in `webFiles.bindings`: a live
// hide of the only Depth track clears them with its resources, and the show
// restores them. Otherwise the rerender fails its binding check, **Live edit
// failed**, and the track stays in the Result until the next Generate.
const LOADED_DEPTH_HIDES = {
  'Show Depth': { stack: false, turn: (page) => toggleDepth(page) },
  "a Depth row's switch": {
    stack: true,
    turn: async (page, mode, shown) => { await slotSwitch(page, mode, 'depth', 0, shown); await settleLive(page); }
  }
};
for (const mode of ['circular', 'linear']) {
  for (const [name, { stack, turn }] of Object.entries(LOADED_DEPTH_HIDES)) {
    test(`${name} after a Load of a saved Session hides and shows the only Depth track live, with Undo and Redo (${mode})`, async ({ page, browser }, testInfo) => {
      test.setTimeout(400_000);
      await open(page, { mode, results: 'single', reflow: 'off' });
      await addDepthFile(page, 'depth.tsv', mode);
      if (stack) {
        await page.evaluate((linear) => {
          const app = window.__GBDRAW_APP__;
          if (linear) app.setLinearTrackSlotsEnabled(true);
          else app.setCircularTrackSlotsEnabled(true);
        }, mode === 'linear');
      }
      await settleLive(page);
      await generate(page);
      const file = testInfo.outputPath('depth.gbdraw-session.json');
      const titled = await page.evaluate(async () => Boolean((await import('./js/state.js')).state.sessionTitle?.value));
      await download(page, 'Save Session', file, titled ? {} : { title: 'depth' });
      const loaded = await load(browser, file);
      const expectApplied = async (count, step) => {
        expect(await loaded.evaluate(async () => (await import('./js/state.js')).state.labelReflowLastError.value?.summary ?? null), `${mode}: ${step}: no Live edit failed`).toBeNull();
        const groups = await shownTrackGroups(loaded, 'DEPTH_GROUPS');
        if (count === 0) expect(groups, `${mode}: ${step}: the track is gone at once`).toBe(0);
        else expect(groups, `${mode}: ${step}: the track is drawn`).toBe(count);
      };
      try {
        const drawn = await shownTrackGroups(loaded, 'DEPTH_GROUPS');
        expect(drawn).toBeGreaterThan(0);
        await turn(loaded, mode, false);
        await expectApplied(0, 'off');
        await turn(loaded, mode, true);
        await expectApplied(drawn, 'on');
        await turn(loaded, mode, false);
        await expectApplied(0, 'off again');
        await history(loaded, 'undo');
        await expectApplied(drawn, 'Undo');
        await history(loaded, 'redo');
        await expectApplied(0, 'Redo');
        await expectLiveEqualsGenerate(loaded, { label: `${mode}: ${name} off after a Load` });
      } finally { await loaded.context().close(); }
    });
  }
}

test('GT B5: Use custom stack turned on after a Generate with it off applies live, with a spacer row (Linear)', async ({ page }) => {
  test.setTimeout(300_000);
  await open(page, { mode: 'linear', results: 'single', reflow: 'off' });
  await addDepthFile(page, 'depth.tsv', 'linear');
  await page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    app.setTrackShown('show_gc', true);
    app.setLinearTrackSlotsEnabled(true);
    app.addLinearTrackSlot('spacer');
    app.setLinearTrackSlotsEnabled(false);
  });
  await settleLive(page);
  await generate(page);
  const before = await page.evaluate(() => window.__GBDRAW_APP__.currentResults()[0].content);
  await inHistoryStep(page, 'Use custom stack', (app) => app.setLinearTrackSlotsEnabled(true));
  await settleLive(page);
  expect(await page.evaluate(() => window.__GBDRAW_APP__.currentResults()[0].content), 'the stack redraws live').not.toBe(before);
  await expectLiveEqualsGenerate(page, { label: 'Use custom stack on after a Generate with it off' });
});

test('GT B8: Undo of a Generate restores the parts it kept, so Show Depth shows live', async ({ page }) => {
  test.setTimeout(300_000);
  await open(page, { mode: 'circular', results: 'single', reflow: 'off' });
  await addDepthFile(page, 'depth.tsv');
  await page.evaluate(() => window.__GBDRAW_APP__.setTrackShown('show_depth', false));
  await settleLive(page);
  await generate(page);
  await inHistoryStep(page, 'Legend font size', (app) => { app.adv.legend_font_size = 17; });
  await settleLive(page);
  await generate(page);
  await history(page, 'undo');
  await toggleDepth(page);
  expect(await shownTrackGroups(page, 'DEPTH_GROUPS'), 'Depth shown live after the Undo of a Generate').toBeGreaterThan(0);
  await expectLiveEqualsGenerate(page, { label: 'Show Depth after the Undo of a Generate' });
});

test('GT B9 (PQ-4 A): in a Circular stack Hide GC Content off after a live hide waits for Generate, with the note', async ({ page }) => {
  test.setTimeout(300_000);
  await open(page, { mode: 'circular', results: 'single', reflow: 'off' });
  await page.evaluate(() => window.__GBDRAW_APP__.setCircularTrackSlotsEnabled(true));
  await settleLive(page);
  await generate(page);
  await hideInStack(page, 'Hide GC Content', 'Hide GC content');
  expect(await shownTrackGroups(page, 'DINUCLEOTIDE_CONTENT_GROUPS'), 'hidden live').toBe(0);
  await clickTrackCheckbox(page, 'Hide GC Content');
  await settleLive(page);
  expect(await shownTrackGroups(page, 'DINUCLEOTIDE_CONTENT_GROUPS'), 'Hide GC Content off waits').toBe(0);
  await expect.poll(() => trackScopeText(page), 'the note names the wait').toContain('Showing a track applies on Generate');
  await generate(page);
  expect(await shownTrackGroups(page, 'DINUCLEOTIDE_CONTENT_GROUPS'), 'Generate draws it').toBeGreaterThan(0);
});

test('GT B1: Hide GC Skew keeps an AT skew row in the Result and the Legend list (HmmtDNA_ATskew)', async ({ browser }) => {
  test.setTimeout(300_000);
  const page = await load(browser, 'gbdraw/web/gallery/sessions/HmmtDNA_ATskew.gbdraw-session.json');
  try {
    await settleLive(page);
    const atRows = async () => (await semanticSnapshot(page)).legend.filter((row) => /^AT skew/i.test(row.caption)).map((row) => row.caption);
    const drawnAt = await atRows();
    expect(drawnAt.length, 'the Session draws an AT skew row').toBeGreaterThan(0);
    await hideInStack(page, 'Hide GC Skew', 'Hide GC skew');
    expect((await semanticSnapshot(page)).legend.filter((row) => /^GC skew/i.test(row.caption)), 'GC skew rows hidden').toEqual([]);
    expect(await atRows(), 'AT skew rows stay on the Result').toEqual(drawnAt);
    expect(await shownListNames(page, /^AT skew/i), 'AT skew rows stay in the Legend list').toEqual(drawnAt);
    await expectLiveEqualsGenerate(page, { label: 'Hide GC Skew with an AT skew row' });
  } finally {
    await page.context().close();
  }
});

// PD-OI-104 scenario revision 6 (D-217, D-219): Hide GC Content in the tobacco
// Session, whose Track row 2 (plastome_tint) extends to the GC content track,
// asks; the choice that turns that row off with the hide applies live as one
// History step, equal to Generate; its Undo turns both back on in the draft.
test('GT B12 (rev 6): Hide GC Content with a row that extends to it turns that row off too (tobacco)', async ({ browser }) => {
  test.setTimeout(300_000);
  const page = await load(browser, 'gbdraw/web/gallery/sessions/tobacco-chloroplast.gbdraw-session.json');
  try {
    await settleLive(page);
    const rows = () => page.evaluate(() => Object.fromEntries(window.__GBDRAW_APP__.adv.circular_track_slots
      .filter((slot) => ['plastome_tint', 'gc_content'].includes(slot.id)).map((slot) => [slot.id, slot.enabled !== false])));
    expect(await rows()).toEqual({ plastome_tint: true, gc_content: true });
    expect(await shownTrackGroups(page, 'DINUCLEOTIDE_CONTENT_GROUPS'), 'the Session draws GC content').toBeGreaterThan(0);
    await clickTrackCheckbox(page, 'Hide GC Content');
    const dialog = page.getByRole('dialog');
    await expect(dialog).toContainText('Track row 2 (plastome_tint)');
    await dialog.getByRole('button', { name: 'Hide GC content and turn off Track row 2', exact: true }).click();
    await settleLive(page);
    expect(await rows(), 'both rows are off').toEqual({ plastome_tint: false, gc_content: false });
    expect(await shownTrackGroups(page, 'DINUCLEOTIDE_CONTENT_GROUPS'), 'hidden live').toBe(0);
    expect(await page.evaluate(() => window.__GBDRAW_APP__.labelReflowLastError ?? null), 'no refusal').toBeNull();
    await history(page, 'undo');
    expect(await rows(), 'Undo turns both on again').toEqual({ plastome_tint: true, gc_content: true });
    expect(await page.evaluate(() => window.__GBDRAW_APP__.form.suppress_gc), 'Undo clears Hide GC Content').toBe(false);
    await history(page, 'redo');
    expect(await rows(), 'Redo turns both off').toEqual({ plastome_tint: false, gc_content: false });
    await expectLiveEqualsGenerate(page, { label: 'Hide GC Content with Track row 2 turned off' });
  } finally {
    await page.context().close();
  }
});

test('GT B10 (OV-2237): on a Result loaded without a request a hide waits for Generate and the note says so', async ({ browser }) => {
  test.setTimeout(300_000);
  const page = await load(browser, 'tests/fixtures/sessions/BGC0000708-BGC0000713.v29.gbdraw-session.json.gz');
  try {
    await settleLive(page);
    const before = await page.evaluate(() => window.__GBDRAW_APP__.currentResults()[window.__GBDRAW_APP__.selectedResultIndex].content);
    const listed = await page.evaluate(() => window.__GBDRAW_APP__.legendEntries.map((entry) => entry.caption));
    await expect.poll(() => page.evaluate(() => [...document.querySelectorAll('[data-applies-on-generate="track-visibility"]')]
      .map((note) => note.textContent.trim())), 'the toggles say Applies on Generate').toContain('Applies on Generate');
    const box = await page.evaluate(() => (window.__GBDRAW_APP__.mode === 'linear' ? 'Show GC Content' : 'Hide GC Content'));
    await clickTrackCheckbox(page, box);
    await settleLive(page);
    expect(await page.evaluate(() => window.__GBDRAW_APP__.currentResults()[window.__GBDRAW_APP__.selectedResultIndex].content),
      'the Result waits for Generate').toBe(before);
    expect(await page.evaluate(() => window.__GBDRAW_APP__.legendEntries.map((entry) => entry.caption)), 'the Legend list follows the Result').toEqual(listed);
  } finally {
    await page.context().close();
  }
});

// Undo of a row switch's choice: the row switched off and the Overlay row that
// reaches it leave the Result and the committed request together, and Undo and
// Redo return both to what they were. A Save after the Undo holds the row.
const UNDO_SESSIONS = {
  circular: { session: 'gbdraw/web/gallery/sessions/tobacco-chloroplast.gbdraw-session.json', record: 'NC_001879.2', set: 'plastome_tint', tint: 'plastome_tint' },
  linear: { session: 'gbdraw/web/gallery/sessions/lambda_basic_linear.gbdraw-session.json', record: 'NC_001416.1', set: 'marks', tint: 'tint' }
};
// A stack with one Depth row, `mark` anchored on it, and a tint row that
// extends to the GC content row, generated.
const openReachingRows = async (browser, mode) => {
  const { session, record, set } = UNDO_SESSIONS[mode];
  const page = await load(browser, session);
  if (mode === 'linear') {
    await page.locator('[data-setting="region-annotations"] input[type=file]').first().setInputFiles({
      name: 'regions.tsv', mimeType: 'text/tab-separated-values',
      buffer: Buffer.from(`set_id\tid\tmark\trecord\tstart\tend\nmarks\tm1\thighlight\t${record}\t100\t900\n`)
    });
    await expect.poll(() => page.evaluate(async () => (await import('./js/state.js')).state.activeDrawing().annotationSets.length)).toBe(1);
    await page.evaluate(() => window.__GBDRAW_APP__.setLinearTrackSlotsEnabled(true));
  }
  await evaluateWithRetainedPromise(page, async ({ stackMode, tsv, setId }) => {
    const { state } = await import('./js/state.js');
    const app = window.__GBDRAW_APP__;
    const file = new File([tsv], 's1.tsv', { type: 'text/tab-separated-values' });
    await window.__GBDRAW_HISTORY__.runUndoable('Setup', () => {
      if (stackMode === 'circular') app.setCircularDepthFile(0, file);
      else app.setLinearDepthFile(app.linearSeqs[0], 0, file);
    });
    await window.__GBDRAW_HISTORY__.runUndoable('Setup', () => {
      const rows = state.activeDrawing().adv[`${stackMode}_track_slots`];
      const row = (id, params) => ({ id, renderer: 'annotations', enabled: true, side: 'overlay', z: params.to_slot ? -1 : 1,
        ...(stackMode === 'linear' ? { height: '', spacing: '' } : {}),
        params: { set_id: setId, marks: ['highlight'], ...params } });
      rows.push(row('mark', { anchor_slot: 'depth', layer: 'foreground' }));
      if (stackMode === 'linear') rows.push(row('tint', { anchor_slot: 'features', to_slot: 'gc_content', layer: 'underlay' }));
    });
    await app.runAnalysis();
  }, { stackMode: mode, tsv: [0, 1, 2, 3].map((index) => `${record}\t${index * 700 + 1}\t${10 + index}`).join('\n'), setId: set });
  await settleLive(page);
  return page;
};
const savedUndoSession = async (page, file) => {
  const titled = await page.evaluate(async () => Boolean((await import('./js/state.js')).state.sessionTitle?.value));
  await download(page, 'Save Session', file, titled ? {} : { title: 'undo' });
  return JSON.parse(gunzipSync(readFileSync(file)).toString('utf8'));
};
const committedSlots = (page, mode) => page.evaluate(async (stackMode) => {
  const { getCommittedCanonicalSession } = await import('./js/services/config.js');
  const tracks = getCommittedCanonicalSession().renderRequest.diagramOptions.tracks;
  return JSON.parse(JSON.stringify(tracks[`${stackMode}TrackSlots`]));
}, mode);
const overlayElements = (page, id) => page.evaluate((slotId) => window.__GBDRAW_APP__.svgContainer.querySelector('svg')
  .querySelectorAll(`[data-gbdraw-slot-id="${slotId}"] > *`).length, id);
const slotIds = (slots) => slots.map((slot) => (typeof slot === 'string' ? slot.split(':')[0] : slot.id));
const chooseToHide = async (page) => {
  await page.locator('[data-choice-prompt] [data-choice-prompt-action="hide"]').click();
  await settleLive(page);
};
for (const mode of ['circular', 'linear']) {
  for (const [target, reacher] of [['depth', 'mark'], ['gc_content', UNDO_SESSIONS[mode].tint]]) {
    for (const after of ['Generate', 'Load']) {
      test(`Undo and Redo of a row switch's choice return the Overlay row to the Result and the committed request (${mode}, ${target} row, after ${after})`, async ({ browser }, testInfo) => {
        test.setTimeout(420_000);
        let page = await openReachingRows(browser, mode);
        try {
          if (after === 'Load') {
            const file = testInfo.outputPath('reached.gbdraw-session.json');
            await savedUndoSession(page, file);
            await page.context().close();
            page = await load(browser, file);
          }
          const drawnBefore = await overlayElements(page, reacher);
          const slotsBefore = await committedSlots(page, mode);
          expect(drawnBefore, 'the Result draws the Overlay row').toBeGreaterThan(0);
          expect(slotIds(slotsBefore)).toContain(reacher);
          expect((await semanticSnapshot(page)).overlays).toContain(reacher);

          await page.evaluate(async ({ stackMode, id }) => {
            const { state } = await import('./js/state.js');
            const app = window.__GBDRAW_APP__;
            const row = state.activeDrawing().adv[`${stackMode}_track_slots`].find((each) => each.id === id);
            void app[stackMode === 'circular' ? 'setCircularTrackSlotEnabled' : 'setLinearTrackSlotEnabled'](row, false);
          }, { stackMode: mode, id: target });
          await chooseToHide(page);
          expect(await overlayElements(page, reacher), 'the choice takes it out').toBe(0);
          const slotsHidden = await committedSlots(page, mode);
          expect(slotIds(slotsHidden)).not.toContain(reacher);

          await history(page, 'undo');
          expect(await overlayElements(page, reacher), 'Undo draws it again').toBe(drawnBefore);
          expect(await committedSlots(page, mode), 'Undo gives the committed request back').toEqual(slotsBefore);
          const saved = await savedUndoSession(page, testInfo.outputPath('after-undo.gbdraw-session.json'));
          expect(slotIds(saved.renderRequest.diagramOptions.tracks[`${mode}TrackSlots`]), 'a Save after the Undo holds the row').toContain(reacher);

          await history(page, 'redo');
          expect(await overlayElements(page, reacher), 'Redo takes it out again').toBe(0);
          expect(await committedSlots(page, mode)).toEqual(slotsHidden);

          await history(page, 'undo');
          await expectLiveEqualsGenerate(page, { label: `${mode} ${target} after ${after}: Undo of the choice` });
        } finally { await page.context().close(); }
      });
    }
  }
}

// Undo of a live Hide GC Content or Hide GC Skew applies at once, with the
// rows that left with it, and no note stays; Redo hides again at once. Turning
// the hide off by hand still waits for Generate, with the note.
const noteOf = (page) => trackScopeText(page);
const committedTracks = (page) => page.evaluate(async () => {
  const { getCommittedCanonicalSession } = await import('./js/services/config.js');
  const { renderRequest } = getCommittedCanonicalSession();
  return JSON.parse(JSON.stringify({ tracks: renderRequest.diagramOptions.tracks, overrides: renderRequest.diagramOptions.configOverrides }));
});
const hideStep = async (page, box, ask) => {
  await clickTrackCheckbox(page, box);
  if (ask) await chooseToHide(page);
  else await settleLive(page);
};
const undoOnce = async (page) => { await history(page, 'undo'); await settleLive(page); };
const redoOnce = async (page) => { await history(page, 'redo'); await settleLive(page); };
const WAIT_NOTE = 'Showing a track applies on Generate';

test('Undo of Hide GC Content with the dialog\'s choice applies at once, with the row that extends to it; Save and Load after it hold them (tobacco, after Load)', async ({ browser }, testInfo) => {
  test.setTimeout(420_000);
  let page = await load(browser, UNDO_SESSIONS.circular.session);
  try {
    const before = await committedTracks(page);
    const tint = await overlayElements(page, 'plastome_tint');
    const groups = await shownTrackGroups(page, 'DINUCLEOTIDE_CONTENT_GROUPS');
    expect(tint).toBeGreaterThan(0);
    await hideStep(page, 'Hide GC Content', true);
    expect(await overlayElements(page, 'plastome_tint')).toBe(0);
    expect(await shownTrackGroups(page, 'DINUCLEOTIDE_CONTENT_GROUPS')).toBe(0);

    await undoOnce(page);
    expect(await overlayElements(page, 'plastome_tint'), 'Undo returns the row').toBe(tint);
    expect(await shownTrackGroups(page, 'DINUCLEOTIDE_CONTENT_GROUPS'), 'Undo returns the track').toBe(groups);
    expect(await committedTracks(page), 'the committed request is what it was').toEqual(before);
    expect(await noteOf(page), 'no note remains').not.toContain(WAIT_NOTE);

    await redoOnce(page);
    expect(await shownTrackGroups(page, 'DINUCLEOTIDE_CONTENT_GROUPS'), 'Redo hides again at once').toBe(0);
    expect(await overlayElements(page, 'plastome_tint')).toBe(0);
    await undoOnce(page);
    const saved = await savedUndoSession(page, testInfo.outputPath('after-undo.gbdraw-session.json'));
    expect(slotIds(saved.renderRequest.diagramOptions.tracks.circularTrackSlots)).toEqual(expect.arrayContaining(['gc_content', 'plastome_tint']));
    await expectLiveEqualsGenerate(page, { label: 'Undo of Hide GC Content with its choice' });
    await page.context().close();
    page = await load(browser, testInfo.outputPath('after-undo.gbdraw-session.json'));
    expect(await overlayElements(page, 'plastome_tint'), 'the loaded Result draws the row').toBe(tint);
    await expectLiveEqualsGenerate(page, { label: 'Load after the Undo' });
  } finally { await page.context().close(); }
});

for (const [box, groupsName, kind] of [['Hide GC Content', 'DINUCLEOTIDE_CONTENT_GROUPS', 'gc_content'], ['Hide GC Skew', 'DINUCLEOTIDE_SKEW_GROUPS', 'gc_skew']]) {
  test(`Undo of ${box} in a Circular stack applies at once after a Generate; a pending value stays out; off by hand still waits`, async ({ page }) => {
    test.setTimeout(420_000);
    await open(page, { mode: 'circular', results: 'single', reflow: 'off' });
    await page.evaluate(() => window.__GBDRAW_APP__.setCircularTrackSlotsEnabled(true));
    await generate(page);
    const before = await committedTracks(page);
    const groups = await shownTrackGroups(page, groupsName);
    expect(groups).toBeGreaterThan(0);
    // A pending row value stays in the draft through the Undo and out of the Result.
    await page.evaluate(async () => {
      const { state } = await import('./js/state.js');
      const rows = state.activeDrawing().adv.circular_track_slots;
      await window.__GBDRAW_HISTORY__.runUndoable('Pending', () => { rows.find((row) => row.id === 'gc_skew').width = '8'; });
    });
    await hideStep(page, box, true);
    expect(await page.evaluate(() => window.__GBDRAW_APP__.labelReflowLastError ?? null), 'no refusal').toBeNull();
    expect(await shownTrackGroups(page, groupsName)).toBe(0);
    await undoOnce(page);
    expect(await shownTrackGroups(page, groupsName), 'Undo returns the track').toBe(groups);
    expect(await committedTracks(page)).toEqual(before);
    expect(await noteOf(page), 'no note remains').not.toContain(WAIT_NOTE);
    expect(await page.evaluate(() => window.__GBDRAW_APP__.adv.circular_track_slots.find((row) => row.id === 'gc_skew').width), 'the value stays pending').toBe('8');
    await redoOnce(page);
    expect(await shownTrackGroups(page, groupsName), 'Redo hides again').toBe(0);
    await clickTrackCheckbox(page, box);
    await settleLive(page);
    expect(await shownTrackGroups(page, groupsName), 'off by hand waits').toBe(0);
    expect(slotIds((await committedTracks(page)).tracks.circularTrackSlots)).not.toContain(kind);
    await expect.poll(() => noteOf(page), 'the note names the wait').toContain('Showing a track applies on Generate');
    await generate(page);
    expect(await shownTrackGroups(page, groupsName)).toBe(groups);
  });

  test(`Undo of ${box} with Circular track widths applies at once (stack off, after a Generate)`, async ({ page }) => {
    test.setTimeout(420_000);
    await open(page, { mode: 'circular', results: 'single', reflow: 'off' });
    await page.evaluate(() => { window.__GBDRAW_APP__.adv.gc_content_width_circular = 40; });
    await generate(page);
    const before = await committedTracks(page);
    const groups = await shownTrackGroups(page, groupsName);
    expect(groups).toBeGreaterThan(0);
    await hideStep(page, box, false);
    expect(await shownTrackGroups(page, groupsName)).toBe(0);
    await undoOnce(page);
    expect(await shownTrackGroups(page, groupsName), 'Undo returns the track').toBe(groups);
    expect(await committedTracks(page)).toEqual(before);
    expect(await noteOf(page), 'no note remains').not.toContain(WAIT_NOTE);
    await expectLiveEqualsGenerate(page, { label: `Undo of ${box} with a track width` });
  });
}

// A later change of the hide ends what its Undo named, also inside that Undo's
// rerender. The rerender's reply is held, so each step lands inside it.
const CONTENT_GROUPS = 'DINUCLEOTIDE_CONTENT_GROUPS';
const SKEW_GROUPS = 'DINUCLEOTIDE_SKEW_GROUPS';
const generateWithTrackWidth = async (page) => {
  await open(page, { mode: 'circular', results: 'single', reflow: 'off' });
  await page.evaluate(() => { window.__GBDRAW_APP__.adv.gc_content_width_circular = 40; });
  await generate(page);
};
const holdNextRerender = (page) => page.evaluate(() => {
  const hold = { entered: false };
  const released = new Promise((resolve) => { hold.release = resolve; });
  window.__UNDO_HOLD__ = hold;
  window.__GBDRAW_TEST_HOOKS__.beforeDiagramGenerationResponse = async () => {
    if (hold.entered) return;
    hold.entered = true;
    await released;
  };
});
const insideHeldRerender = (page) => expect.poll(() => page.evaluate(() => window.__UNDO_HOLD__.entered), { timeout: 120_000 }).toBe(true);
const releaseRerender = async (page) => {
  await page.evaluate(() => {
    window.__UNDO_HOLD__.release();
    delete window.__GBDRAW_TEST_HOOKS__.beforeDiagramGenerationResponse;
  });
  await settleLive(page);
};
const step = (page, name) => page.evaluate((each) => window.__GBDRAW_HISTORY__[each](), name);
const undoWithHeldRerender = async (page) => {
  await holdNextRerender(page);
  await step(page, 'undo');
  await insideHeldRerender(page);
};
const handOffStates = [
  ['Redo, Redo', async (page) => { await step(page, 'redo'); await step(page, 'redo'); }],
  ['two clicks on the hide', async (page) => { await clickTrackCheckbox(page, 'Hide GC Content'); await clickTrackCheckbox(page, 'Hide GC Content'); }]
];
for (const [name, handOff] of handOffStates) {
  test(`Off by hand after a live hide keeps its wait when ${name} land inside the rerender of an Undo`, async ({ page }) => {
    test.setTimeout(420_000);
    await generateWithTrackWidth(page);
    await hideStep(page, 'Hide GC Content', false);
    await hideStep(page, 'Hide GC Content', false);
    await undoOnce(page);
    await undoWithHeldRerender(page);
    await handOff(page);
    await releaseRerender(page);
    expect(await page.evaluate(() => window.__GBDRAW_APP__.form.suppress_gc), 'the hide is off in the draft').toBe(false);
    expect(await shownTrackGroups(page, CONTENT_GROUPS), 'it waits for Generate').toBe(0);
    expect(await noteOf(page), 'the note names the wait').toContain(WAIT_NOTE);
  });
}

test('Hide, off by hand, Undo, Redo keeps the wait', async ({ page }) => {
  test.setTimeout(420_000);
  await generateWithTrackWidth(page);
  await hideStep(page, 'Hide GC Content', false);
  await hideStep(page, 'Hide GC Content', false);
  await undoOnce(page);
  await redoOnce(page);
  expect(await shownTrackGroups(page, CONTENT_GROUPS), 'it waits for Generate').toBe(0);
  expect(await noteOf(page), 'the note names the wait').toContain(WAIT_NOTE);
});

test('Undo, Undo inside the rerender returns both tracks, with no note while it runs; Undo, Redo, Undo inside it returns the track', async ({ page }) => {
  test.setTimeout(420_000);
  await generateWithTrackWidth(page);
  const before = await committedTracks(page);
  const groups = [await shownTrackGroups(page, CONTENT_GROUPS), await shownTrackGroups(page, SKEW_GROUPS)];
  await hideStep(page, 'Hide GC Content', false);
  await hideStep(page, 'Hide GC Skew', false);
  await undoWithHeldRerender(page);
  expect(await noteOf(page), 'the Undo is applying: nothing waits for Generate').not.toContain(WAIT_NOTE);
  await step(page, 'undo');
  expect(await noteOf(page), 'nor after the second Undo').not.toContain(WAIT_NOTE);
  await releaseRerender(page);
  expect([await shownTrackGroups(page, CONTENT_GROUPS), await shownTrackGroups(page, SKEW_GROUPS)], 'both tracks return').toEqual(groups);
  expect(await committedTracks(page)).toEqual(before);
  expect(await noteOf(page), 'no note remains').not.toContain(WAIT_NOTE);

  await hideStep(page, 'Hide GC Content', false);
  await undoWithHeldRerender(page);
  await step(page, 'redo');
  await step(page, 'undo');
  await releaseRerender(page);
  expect(await shownTrackGroups(page, CONTENT_GROUPS), 'the last Undo returns the track').toBe(groups[0]);
  expect(await committedTracks(page)).toEqual(before);
});

test('An Undo whose rerender fails keeps the Result and shows the wait; the next off by hand waits too', async ({ page }) => {
  test.setTimeout(420_000);
  await generateWithTrackWidth(page);
  await hideStep(page, 'Hide GC Content', false);
  await page.evaluate(() => { window.__GBDRAW_TEST_HOOKS__.beforeDiagramGenerationResponse = () => { throw new Error('forced rerender failure'); }; });
  await undoOnce(page);
  await page.evaluate(() => { delete window.__GBDRAW_TEST_HOOKS__.beforeDiagramGenerationResponse; });
  expect(await page.evaluate(() => window.__GBDRAW_APP__.labelReflowLastError), 'the failure is reported').not.toBeNull();
  expect(await shownTrackGroups(page, CONTENT_GROUPS), 'the Result stays').toBe(0);
  await expect.poll(() => noteOf(page), 'the note is back').toContain(WAIT_NOTE);
  await hideStep(page, 'Hide GC Content', false);
  await hideStep(page, 'Hide GC Content', false);
  expect(await shownTrackGroups(page, CONTENT_GROUPS), 'off by hand waits').toBe(0);
  await expect.poll(() => noteOf(page)).toContain(WAIT_NOTE);
});

const resetSettings = async (page) => {
  await evaluateWithRetainedPromise(page, async () => { await window.__GBDRAW_APP__.resetSettings(); });
};
const hideIsOn = (page) => page.evaluate(() => window.__GBDRAW_APP__.form.suppress_gc);

test('Redo of Reset Settings inside the rerender of an Undo keeps the wait: a restore that turns the hide on drops its name', async ({ page }) => {
  test.setTimeout(420_000);
  await generateWithTrackWidth(page);
  await hideStep(page, 'Hide GC Content', false);
  await resetSettings(page);
  await settleLive(page);
  await hideStep(page, 'Hide GC Content', false);
  await undoWithHeldRerender(page);
  await step(page, 'undo');
  expect(await hideIsOn(page), 'the Undo of Reset Settings turned the hide on').toBe(true);
  await step(page, 'redo');
  expect(await noteOf(page), 'the Redo of Reset Settings hands the hide off').toContain(WAIT_NOTE);
  await releaseRerender(page);
  expect(await hideIsOn(page)).toBe(false);
  expect(await shownTrackGroups(page, CONTENT_GROUPS), 'it waits for Generate').toBe(0);
  expect(await noteOf(page), 'the note names the wait').toContain(WAIT_NOTE);
});

test('An Undo that leaves its hide on names nothing', async ({ page }) => {
  test.setTimeout(420_000);
  await generateWithTrackWidth(page);
  await hideStep(page, 'Hide GC Content', false);
  await hideStep(page, 'Hide GC Content', false);
  await holdNextRerender(page);
  await clickTrackCheckbox(page, 'Hide GC Skew');
  await insideHeldRerender(page);
  await step(page, 'undo');
  await step(page, 'undo');
  expect(await hideIsOn(page), 'the second Undo leaves the hide on').toBe(true);
  expect(await page.evaluate(async () => {
    const { state } = await import('./js/state.js');
    return [...state.undoneTrackHides.value.get('circular')];
  })).toEqual(['skew']);
  await releaseRerender(page);
});

test('Undo, Undo through Reset Settings inside the rerender returns the track: a restore that leaves the hide off keeps the name', async ({ page }) => {
  test.setTimeout(420_000);
  await generateWithTrackWidth(page);
  const shown = await shownTrackGroups(page, CONTENT_GROUPS);
  await resetSettings(page);
  await settleLive(page);
  await hideStep(page, 'Hide GC Content', false);
  await undoWithHeldRerender(page);
  await step(page, 'undo');
  expect(await hideIsOn(page)).toBe(false);
  await releaseRerender(page);
  expect(await shownTrackGroups(page, CONTENT_GROUPS), 'the track returns').toBe(shown);
  expect(await noteOf(page), 'no note remains').not.toContain(WAIT_NOTE);
});

// A Generate waits for the running rerender; its Cancel ends that rerender too.
const generateAndCancel = async (page) => {
  await page.evaluate(() => { void window.__GBDRAW_APP__.runAnalysis(); });
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.processing)).toBe(true);
  await evaluateWithRetainedPromise(page, async () => { await window.__GBDRAW_APP__.cancelGeneration(); });
  await releaseRerender(page);
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.processing)).toBe(false);
  await settleLive(page);
};
const resultContent = (page) => page.evaluate(() => window.__GBDRAW_APP__.results[0].content);

test('Cancel of a Generate inside the rerender of an Undo leaves the Result and shows the wait; the next off by hand waits too', async ({ page }) => {
  test.setTimeout(420_000);
  await generateWithTrackWidth(page);
  await hideStep(page, 'Hide GC Content', false);
  const tracks = await committedTracks(page);
  const result = await resultContent(page);
  expect(result).toContain('<svg');
  await undoWithHeldRerender(page);
  await generateAndCancel(page);
  expect(await hideIsOn(page)).toBe(false);
  expect(await shownTrackGroups(page, CONTENT_GROUPS), 'the track stays out').toBe(0);
  await expect.poll(() => noteOf(page), 'the note is back').toContain(WAIT_NOTE);
  expect(await committedTracks(page), 'the committed tracks are unchanged').toEqual(tracks);
  expect((await resultContent(page)) === result, 'the Result is unchanged').toBe(true);
  await hideStep(page, 'Hide GC Content', false);
  await hideStep(page, 'Hide GC Content', false);
  expect(await shownTrackGroups(page, CONTENT_GROUPS), 'off by hand waits').toBe(0);
  await expect.poll(() => noteOf(page)).toContain(WAIT_NOTE);
});

test('Cancel of a Generate inside the rerender of an Undo keeps the names when another Undo queued an iteration', async ({ page }) => {
  test.setTimeout(420_000);
  await generateWithTrackWidth(page);
  const groups = [await shownTrackGroups(page, CONTENT_GROUPS), await shownTrackGroups(page, SKEW_GROUPS)];
  await hideStep(page, 'Hide GC Content', false);
  await hideStep(page, 'Hide GC Skew', false);
  await undoWithHeldRerender(page);
  await step(page, 'undo');
  await generateAndCancel(page);
  expect([await shownTrackGroups(page, CONTENT_GROUPS), await shownTrackGroups(page, SKEW_GROUPS)], 'both tracks return').toEqual(groups);
  expect(await noteOf(page), 'no note remains').not.toContain(WAIT_NOTE);
});
