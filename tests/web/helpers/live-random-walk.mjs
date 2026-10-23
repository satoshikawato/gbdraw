// A seeded random walk over the live edit kinds (PD-OI-066, gbdraw/web/CLAUDE.md
// R3): after every step, the displayed Result must equal the Result the next
// Generate draws from the same draft (expectLiveEqualsGenerate, minus
// tests/web/contracts/live-generate-parity-allowed.json), and the step must
// leave the slices of the drawing it did not edit as they were (M6 isolation:
// the sha256 of each slice's `capture`, as tests/web/per-mode-drawings.test.mjs
// hashes them in node). Reset Settings and Session Load change every drawing
// by design (PD-OI-090), so the isolation check skips those steps.
//
// Two specs run it: the PR spec tests/web/live-generate-random-walk.playwright.spec.js
// (fixed seed, 10 steps per fixture) and the promotion spec
// tests/web/live-generate-random-walk.promotion.spec.js (20 steps, seed
// 20261006 by default). Each fixture opens both drawings with a Result, so a
// drawing switch arrives at a Result. The kinds are stratified: `planWalk`
// gives every kind at least one step over the three fixtures, and draws the
// popup state (a feature popup opened and closed before the step, which binds
// the labels) and Auto Reflow for every step. A failure names the seed, the
// fixture, the step index, the steps so far, and the difference; the same seed
// and budget replay the same walk.
import { expect, test } from '@playwright/test';
import { readFileSync } from 'node:fs';
import { evaluateWithRetainedPromise } from './app-lifecycle.cjs';
import { expectLiveEqualsGenerate, settleLive, showResult } from './live-generate-parity.cjs';
import { download } from './mode-transition.cjs';
import { SINGLE_FIXTURE, generate, open } from './live-generate-parity-steps.cjs';

// The step kinds. Each run draws every kind at least once.
export const WALK_KINDS = Object.freeze([
  'fill', 'stroke', 'featureVisibility', 'labelText', 'labelVisibility',
  'visibilityRuleAdd', 'visibilityRuleEdit', 'visibilityRuleRemove', 'colorRuleAdd', 'colorRuleEdit', 'colorRuleRemove',
  'reflow', 'legendRename', 'legendColor', 'legendSort', 'legendSortDefault', 'legendDelete', 'legendRestore', 'legendAdd',
  'trackData', 'palette', 'undo', 'undoRedo', 'resultSwitch', 'reset', 'load', 'drawingSwitch'
]);

// The kinds that change every drawing by design (PD-OI-090 Reset, Save and Load).
const ALL_DRAWINGS_KINDS = Object.freeze(['reset', 'load']);

// The step kind that walks the path of each row of the R3 work allowlist
// (tests/web/helpers/live-work-allowlist.mjs), or null for a row whose path
// the walk does not draw: a Save that races a Result display, and Embedded
// Only label rendering with its Feature Edits TSV.
export const ALLOWLIST_WALK_KINDS = Object.freeze({
  'feature stroke (popup)': 'stroke',
  'Legend row stroke': 'stroke',
  'Legend row color (no rule)': 'legendColor',
  'History step (Undo of a Legend row color)': 'undo',
  'palette change': 'palette',
  'Apply to all on a palette row': 'fill',
  'palette switch, Keep my colors': 'palette',
  'Result display after a palette change and a stroke': 'resultSwitch',
  'Feature visibility rule add': 'visibilityRuleAdd',
  'color rule commit': 'colorRuleAdd',
  'Result display, nothing changed': 'resultSwitch',
  'Result display declined by Save, shown when Save ends': null,
  'Save while the displayed Result lacks paint': null,
  'Feature Edits TSV load (Embedded Only)': null,
  'Feature Edits TSV load, label Off (Embedded Only)': null,
  'popup label text (Embedded Only)': null,
  'Legend rename, row without features': 'legendRename',
  'color rule commit that splits a row': 'colorRuleAdd',
  'Legend rename, rule row one Result draws (whole row)': 'legendRename',
  'History step (Undo of the rule row rename)': 'undo',
  'History step (Redo of the rule row rename)': 'undoRedo',
  'Legend row color, rule row': 'legendColor',
  'Rules-panel color of a captioned rule': 'colorRuleEdit',
  'Legend rename, rule row both Results draw': 'legendRename',
  'Legend rename, palette row with features': 'legendRename',
  'Legend delete': 'legendDelete',
  'History step (Undo of a Legend delete)': 'undo',
  'Legend Restore': 'legendRestore',
  'Legend sort': 'legendSort',
  'mode switch back to a Result': 'drawingSwitch',
  'Session Load': 'load',
  'Legend rename onto a deleted row (Restore and merge)': 'legendRename',
  'Reset Settings': 'reset'
});

// The confirmed mismatches of a step kind (multi-drawing-project FINDINGS.md),
// expected until their fix: a difference line each one names.
/** @type {Readonly<Record<string, { id: string, matches: (line: string) => boolean }>>} */
const KNOWN_MISMATCHES = Object.freeze({});

// The fixtures: the shown drawing with its Result, and the other drawing's
// Result (Linear always draws the one-record file).
export const FIXTURES = Object.freeze([
  { name: 'Circular, one Result', mode: 'circular', results: 'single' },
  { name: 'Linear, one Result', mode: 'linear', results: 'single' },
  { name: 'Circular, two-Result batch', mode: 'circular', results: 'batch' }
]);
// The records of each file, for a Depth file over them.
const FILE_RECORDS = { single: ['FORCEDLBL'], batch: ['TESTA', 'TESTB'] };

// mulberry32: a 32-bit seeded generator; the same seed gives the same walk.
const createRandom = (seed) => {
  let state = seed >>> 0;
  const next = () => {
    state = (state + 0x6D2B79F5) >>> 0;
    let value = state;
    value = Math.imul(value ^ (value >>> 15), value | 1);
    value ^= value + Math.imul(value ^ (value >>> 7), value | 61);
    return ((value ^ (value >>> 14)) >>> 0) / 4294967296;
  };
  const int = (count) => Math.floor(next() * count);
  const pick = (items) => items[int(items.length)];
  const weighted = (entries) => {
    const total = entries.reduce((sum, [, weight]) => sum + weight, 0);
    let point = next() * total;
    for (const [value, weight] of entries) {
      point -= weight;
      if (point < 0) return value;
    }
    return entries[entries.length - 1][0];
  };
  const shuffle = (items) => {
    const copy = [...items];
    for (let index = copy.length - 1; index > 0; index -= 1) {
      const other = int(index + 1);
      [copy[index], copy[other]] = [copy[other], copy[index]];
    }
    return copy;
  };
  return { next, int, pick, weighted, shuffle };
};

// The steps after the kinds every run draws: the edits, weighted as the
// promotion walk drew them. A Result switch, a drawing switch, Reset, and
// Load are drawn once.
const FILLER_KINDS = [
  ['fill', 5], ['stroke', 2], ['featureVisibility', 2], ['labelText', 3], ['labelVisibility', 3],
  ['visibilityRuleAdd', 2], ['visibilityRuleEdit', 2], ['visibilityRuleRemove', 2],
  ['colorRuleAdd', 2], ['colorRuleEdit', 2], ['colorRuleRemove', 1], ['reflow', 1],
  ['legendRename', 2], ['legendColor', 2], ['legendSort', 2], ['legendSortDefault', 1],
  ['legendDelete', 1], ['legendRestore', 1], ['legendAdd', 1], ['trackData', 1], ['palette', 1],
  ['undo', 2], ['undoRedo', 2]
];
// Only the batch has two Results to switch between. Reset clears its batch
// settings and its drawing switch would leave it, so neither runs there.
const allowedIn = (kind, fixture) => {
  if (kind === 'resultSwitch') return fixture.results === 'batch';
  if (kind === 'reset' || kind === 'drawingSwitch') return fixture.results !== 'batch';
  return true;
};

/**
 * The kind schedule of a walk: per fixture, `steps` entries of
 * `{ kind, popup, reflow }`. Every kind of WALK_KINDS has at least one step.
 * @param {number} seed
 * @param {number} steps
 */
export const planWalk = (seed, steps) => {
  if (steps * FIXTURES.length < WALK_KINDS.length) {
    throw new Error(`a walk of ${steps} steps per fixture cannot draw all ${WALK_KINDS.length} kinds`);
  }
  const random = createRandom(seed ^ 0x5EED5EED);
  const slots = FIXTURES.map(() => /** @type {string[]} */ ([]));
  const fixturesOf = (kind) => FIXTURES.flatMap((fixture, index) => (allowedIn(kind, fixture) ? [index] : []));
  // The kinds with the fewest fixtures are placed first; each goes to an
  // allowed fixture with the fewest kinds.
  const kinds = random.shuffle(WALK_KINDS).sort((left, right) => fixturesOf(left).length - fixturesOf(right).length);
  for (const kind of kinds) {
    const open = fixturesOf(kind).filter((index) => slots[index].length < steps);
    if (!open.length) throw new Error(`no fixture has a step left for ${kind}`);
    const fewest = Math.min(...open.map((index) => slots[index].length));
    slots[random.pick(open.filter((index) => slots[index].length === fewest))].push(kind);
  }
  return slots.map((kindsOfFixture) => {
    while (kindsOfFixture.length < steps) kindsOfFixture.push(random.weighted(FILLER_KINDS));
    return random.shuffle(kindsOfFixture).map((kind) => ({ kind, popup: random.next() < 0.5, reflow: random.next() < 0.5 }));
  });
};

const COLORS = ['#c83366', '#123456', '#2a9d8f', '#e63946', '#7b2cbf', '#f4a261', '#264653'];
const LABEL_TEXTS = [
  'alpha protein, edited live',
  'x',
  'walk label',
  'a label text far too long to stay inside its own arc, drawn outside with leaders'
];
const CAPTIONS = ['Walk caption', 'transfer RNA', 'Renamed group', 'CDS'];
const ADDED_CAPTIONS = ['Walk group', 'Added row'];
const escapeRegex = (text) => text.replace(/[.*+?^${}()|[\]\\]/g, '\\$&');
const otherMode = (mode) => (mode === 'circular' ? 'linear' : 'circular');

// One drawing switch (the mode buttons call `selectDrawing`), settled.
const switchDrawing = async (page, drawingId) => {
  await settleLive(page);
  const outcome = await page.evaluate((id) => window.__GBDRAW_APP__.selectDrawing(id), drawingId);
  expect(outcome, `the switch to ${drawingId}`).toEqual({ status: 'ok' });
  await expect.poll(() => page.evaluate(async () => (await import('/gbdraw/web/js/state.js')).state.circularRecordDiscovery.status))
    .not.toBe('loading');
  await page.evaluate(async () => { await window.Vue.nextTick(); });
  await settleLive(page);
};

// The fixture: its drawing opened as the parity specs open it (labels drawn
// outside), then the other drawing with a Result, then its drawing shown.
const openFixture = async (page, fixture, reflow) => {
  await open(page, { mode: fixture.mode, results: fixture.results, reflow: reflow ? 'on' : 'off' });
  await switchDrawing(page, otherMode(fixture.mode));
  if (fixture.mode === 'circular') {
    await page.evaluate(async (text) => {
      const app = window.__GBDRAW_APP__;
      app.setLinearSeqPrimaryFile(0, 'gb', new File([text], 'forced_label_underlay.gb', { type: 'text/plain', lastModified: 1000 }));
      await window.Vue.nextTick();
    }, readFileSync(SINGLE_FIXTURE, 'utf8'));
    await settleLive(page);
    await page.evaluate(() => { window.__GBDRAW_APP__.form.show_labels_linear = 'all'; });
  } else {
    await page.getByLabel('GenBank/DDBJ File', { exact: true }).setInputFiles(SINGLE_FIXTURE);
    await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.circularRecordList.length), { timeout: 60_000 })
      .toBeGreaterThan(0);
    await page.evaluate(() => { window.__GBDRAW_APP__.form.labels_mode = 'out'; });
    await settleLive(page);
  }
  await generate(page);
  const otherFloor = await page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount());
  await switchDrawing(page, fixture.mode);
  return { [otherMode(fixture.mode)]: otherFloor };
};

// The file each drawing draws.
const fileOf = (fixture, drawingId) => (drawingId === 'circular' && fixture.mode === 'circular' ? fixture.results : 'single');

// The sha256 of each slice capture of each drawing, and the shown drawing.
const sliceHashes = (page) => page.evaluate(async () => {
  const { state } = await import('/gbdraw/web/js/state.js');
  const { DRAWING_SLICES } = await import('/tests/web/helpers/drawing-state.mjs');
  const hex = async (text) => [...new Uint8Array(await crypto.subtle.digest('SHA-256', new TextEncoder().encode(text)))]
    .map((byte) => byte.toString(16).padStart(2, '0')).join('');
  const drawings = {};
  for (const drawing of state.drawingList()) {
    const slices = {};
    for (const [name, slice] of Object.entries(DRAWING_SLICES)) slices[name] = await hex(JSON.stringify(slice.capture(drawing.slices[name])));
    drawings[drawing.id] = slices;
  }
  return { active: state.activeDrawing().id, drawings };
});
// The slices of `drawingId` whose captures differ between two hash sets.
const changedSlices = (before, after, drawingId) => Object.keys(before.drawings[drawingId])
  .filter((name) => before.drawings[drawingId][name] !== after.drawings[drawingId][name]);

// What the next step may use: the features of the displayed Result (sorted by
// a stable key), the editor lists, and the History and Result state. Undo is
// available above the drawing's floor, its History depth after its first
// Generate (or after a Load): Undo of that Generate leaves no Result to compare.
const readState = (page) => page.evaluate(async () => {
  const { state } = await import('/gbdraw/web/js/state.js');
  const app = window.__GBDRAW_APP__;
  const history = window.__GBDRAW_HISTORY__;
  const featureKey = (feature) => [feature.type, feature.start, feature.end, feature.locus_tag || '', feature.product || ''].join('|');
  const features = app.filteredFeatures.map((feature) => ({
    key: featureKey(feature), type: String(feature.type || ''), tag: String(feature.locus_tag || ''), product: String(feature.product || '')
  })).sort((left, right) => (left.key < right.key ? -1 : left.key > right.key ? 1 : 0));
  return {
    drawing: state.activeDrawing().id,
    features,
    resultCount: app.results.length,
    resultIndex: app.selectedResultIndex,
    reflow: Boolean(app.autoLabelReflowEnabled),
    visibilityRules: app.featureVisibilityManualRules.length,
    colorRules: app.manualSpecificRules.length,
    legend: app.legendEntries.map((entry) => String(entry.caption || '')),
    deletedLegend: app.deletedLegendEntries.length,
    palettes: app.paletteNames.filter((name) => name !== app.selectedPalette),
    depth: app.mode === 'linear' ? Boolean(app.linearSeqs[0]?.depth?.[0]) : Boolean(app.files.c_depth?.[0]?.[0]),
    undoDepth: history.getUndoCount(),
    canUndo: Boolean(history.canUndo())
  };
});

// The arguments of a step of `kind`, drawn from the generator and the state:
// { kind, ...arguments }. `picks` answers the dialogs the edit opens, in the
// order they open. A kind whose subject is missing first runs the step that
// makes one (`prepare`), in the same step.
const nextStep = (random, state, kind, { fixture, canUndo }) => {
  const features = state.features;
  const tagged = features.filter((feature) => feature.tag);
  const feature = () => random.pick(features).key;
  const rulePattern = () => {
    const target = random.pick(tagged);
    return { target, value: random.next() < 0.5 ? `^${escapeRegex(target.tag)}$` : `${escapeRegex(target.tag.slice(-3))}$` };
  };
  const index = (count) => random.int(count);
  const picks = [random.int(1000), random.int(1000), random.int(1000)];
  const visibilityRuleAdd = () => {
    const { target, value } = rulePattern();
    return {
      kind: 'visibilityRuleAdd',
      fields: { recordId: '*', featureType: target.type, qualifier: 'locus_tag', value, action: random.pick(['off', 'show']) }
    };
  };
  const colorRuleAdd = () => {
    const { target, value } = rulePattern();
    return {
      kind: 'colorRuleAdd',
      fields: { feat: target.type, qual: 'locus_tag', val: value, color: random.pick(COLORS), cap: random.pick(CAPTIONS) }
    };
  };
  const fill = () => ({ kind: 'fill', feature: feature(), color: random.pick(COLORS), picks });
  switch (kind) {
    case 'fill': return fill();
    case 'stroke': return random.next() < 0.5
      ? { kind, target: 'popup', feature: feature(), color: random.pick(COLORS), picks }
      : { kind, target: 'legend', index: index(state.legend.length), color: random.pick(COLORS), picks };
    case 'featureVisibility': return { kind, feature: feature(), value: random.pick(['on', 'off']), picks };
    case 'labelText': return { kind, feature: feature(), text: random.pick(LABEL_TEXTS), picks };
    case 'labelVisibility': return { kind, feature: feature(), value: random.pick(['on', 'off']), picks };
    case 'visibilityRuleAdd': return visibilityRuleAdd();
    case 'visibilityRuleEdit':
      return state.visibilityRules
        ? { kind, index: index(state.visibilityRules), field: 'action', value: random.pick(['off', 'show']) }
        : { kind, prepare: visibilityRuleAdd(), index: 0, field: 'action', value: random.pick(['off', 'show']) };
    case 'visibilityRuleRemove':
      return state.visibilityRules ? { kind, index: index(state.visibilityRules) } : { kind, prepare: visibilityRuleAdd(), index: 0 };
    case 'colorRuleAdd': return colorRuleAdd();
    case 'colorRuleEdit': {
      const field = random.pick(['color', 'cap']);
      const value = field === 'color' ? random.pick(COLORS) : random.pick(CAPTIONS);
      return state.colorRules
        ? { kind, index: index(state.colorRules), field, value, picks }
        : { kind, prepare: colorRuleAdd(), index: 0, field, value, picks };
    }
    case 'colorRuleRemove':
      return state.colorRules ? { kind, index: index(state.colorRules) } : { kind, prepare: colorRuleAdd(), index: 0 };
    case 'reflow': return { kind, enabled: !state.reflow };
    case 'legendRename': return { kind, index: index(state.legend.length), caption: random.pick(CAPTIONS), picks };
    case 'legendColor': return { kind, index: index(state.legend.length), color: random.pick(COLORS) };
    case 'legendSort': return { kind, direction: random.pick(['asc', 'desc']) };
    case 'legendSortDefault': return { kind };
    case 'legendDelete': return { kind, index: index(state.legend.length), picks };
    case 'legendRestore':
      return state.deletedLegend ? { kind, picks } : { kind, prepare: { kind: 'legendDelete', index: index(state.legend.length), picks }, picks };
    case 'legendAdd': return { kind, feature: feature(), caption: random.pick(ADDED_CAPTIONS), picks };
    case 'trackData': {
      // A track the Result lacks applies on Generate (R1), so a drawing
      // without a Depth file first gets one drawn; the step removes it.
      if (state.depth) return { kind };
      const records = FILE_RECORDS[fileOf(fixture, state.drawing)];
      const tsv = records.flatMap((record) => Array.from({ length: 4 }, (_, row) => `${record}\t${row * 700 + 1}\t${10 + row + random.int(20)}`))
        .join('\n');
      return { kind, prepare: { kind: 'depthFile', tsv } };
    }
    case 'palette': return { kind, palette: random.pick(state.palettes), picks };
    case 'undo': return canUndo ? { kind } : { kind, prepare: fill() };
    case 'undoRedo': return canUndo ? { kind } : { kind, prepare: fill() };
    case 'resultSwitch': {
      const others = [...Array(state.resultCount).keys()].filter((value) => value !== state.resultIndex);
      if (!others.length) throw new Error(`a Result switch needs two Results; the drawing shows ${state.resultCount}`);
      return { kind, index: random.pick(others) };
    }
    case 'drawingSwitch': return { kind, drawing: otherMode(state.drawing) };
    case 'reset':
    case 'load':
      return { kind };
    default: throw new Error(`unknown step kind ${kind}`);
  }
};

// One step through the actions the editor controls call. The dialogs the edit
// opens are answered from `picks`.
const applyStep = (page, step) => evaluateWithRetainedPromise(page, async (step) => {
  const app = window.__GBDRAW_APP__;
  const history = window.__GBDRAW_HISTORY__;
  const notes = [];
  window.__GBDRAW_WALK_NOTES__ = notes;
  const sleep = (ms) => new Promise((resolve) => setTimeout(resolve, ms));
  const featureKey = (feature) => [feature.type, feature.start, feature.end, feature.locus_tag || '', feature.product || ''].join('|');

  const run = async (current) => {
    let picked = 0;
    const choose = (choices) => choices[(current.picks?.[picked++ % current.picks.length] ?? 0) % choices.length];

    // Answers the dialogs the edit opens until the edit has finished and none is open.
    const answerDialogs = async (applied) => {
      // An answer can wait for a later dialog (a hidden label text answer opens
      // Label visibility On), so answers are started, not awaited, in the loop.
      const answers = [];
      let outstanding = 0;
      const answer = (result) => {
        outstanding += 1;
        const settled = Promise.resolve(result).finally(() => { outstanding -= 1; });
        answers.push(settled);
        settled.catch(() => {});
      };
      let finished = false;
      applied.then(() => { finished = true; }, () => { finished = true; });
      const deadline = Date.now() + 150_000;
      // A dialog is answered once, until it has closed.
      let answering = '';
      while (Date.now() < deadline) {
        const style = app.featureStyleScopeDialog;
        const visibility = app.featureVisibilityScopeDialog;
        const shown = ['featureStyleScopeDialog', 'featureVisibilityScopeDialog', 'labelOnDialog', 'labelTextScopeDialog',
          'hiddenLabelTextDialog', 'legendRenameDialog', 'paletteColorsDialog'].find((name) => app[name].show) || '';
        if (!shown) answering = '';
        if (shown && shown === answering) {
          await sleep(20);
          continue;
        }
        answering = shown;
        if (style.show) {
          const offered = ['single'];
          if (style.matchingRule) offered.push('rule');
          if (style.siblingCount > 0) offered.push('caption');
          // OV-648: the dialog offers a label scope that selects other features.
          if (style.displayLabelOffered) offered.push('displayLabel');
          if (style.annotationLabelOffered) offered.push('annotationLabel');
          if (style.kind === 'fill' && style.existingCaptionColor) offered.push('useExisting');
          const choice = choose(offered);
          notes.push(`style scope ${choice} of [${offered}]`);
          answer(app.handleFeatureStyleScopeChoice(choice));
        } else if (visibility.show) {
          const offered = (visibility.scopes || []).map((scope) => scope.id);
          const choice = choose(offered);
          notes.push(`visibility scope ${choice} of [${offered}]`);
          answer(app.handleFeatureVisibilityScopeChoice(choice));
        } else if (app.labelOnDialog.show) {
          const reason = app.labelOnDialog.reason;
          const choice = reason === 'hidden' ? choose(['show', 'keep']) : choose(['keep', 'cancel']);
          notes.push(`label On ${reason}: ${choice}`);
          answer(app.handleLabelOnChoice(choice));
        } else if (app.labelTextScopeDialog.show) {
          const choice = choose(['single', 'all']);
          notes.push(`label text scope ${choice}`);
          answer(app.handleLabelTextScopeChoice(choice));
        } else if (app.hiddenLabelTextDialog.show) {
          const choice = choose(['show', 'text_only', 'cancel']);
          notes.push(`hidden label text ${choice}`);
          answer(app.handleHiddenLabelTextChoice(choice));
        } else if (app.legendRenameDialog.show) {
          const choice = choose(app.legendRenameDialog.mode === 'target' ? ['merge', 'keep-both'] : ['single', 'group']);
          notes.push(`legend rename ${app.legendRenameDialog.mode} ${choice}`);
          answer(app.handleLegendRenameChoice(choice));
        } else if (app.paletteColorsDialog.show) {
          const choice = choose(app.paletteColorsDialog.kind === 'reset' ? ['reset', 'cancel'] : ['keep', 'palette']);
          notes.push(`palette colors ${choice}`);
          answer(app.handlePaletteColorsChoice(choice));
        } else if (finished && outstanding === 0) {
          break;
        } else {
          // Once the edit has settled, its promise wins the race at once: the
          // loop would then never yield a task, and an answer that waits for one
          // (a Worker reply) would never settle while the timers pile up until
          // the renderer runs out of memory (OV-330).
          await (finished ? sleep(20) : Promise.race([applied.catch(() => {}), sleep(20)]));
        }
      }
      await Promise.all([applied, ...answers]);
    };

    const withPopup = async (edit) => {
      const feature = app.filteredFeatures.find((item) => featureKey(item) === current.feature);
      if (!feature) throw new Error(`no feature ${current.feature}`);
      await app.openFeatureEditorFromList(feature, null);
      await window.Vue.nextTick();
      if (!app.clickedFeature) throw new Error(`no popup for ${current.feature}`);
      try {
        await edit();
      } finally {
        app.clickedFeature = null;
      }
    };

    if (current.prepare) {
      notes.push(`prepare ${current.prepare.kind}`);
      await run(current.prepare);
    }
    switch (current.kind) {
      case 'fill':
        await withPopup(() => answerDialogs(Promise.resolve(app.updateClickedFeatureColor(current.color))));
        break;
      case 'stroke':
        if (current.target === 'popup') {
          await withPopup(() => answerDialogs(Promise.resolve(app.setClickedFeatureStrokeColorValue(current.color))));
        } else {
          await answerDialogs(Promise.resolve(app.setLegendEntryStrokeColorValue(current.index, current.color)));
        }
        break;
      case 'featureVisibility':
        await withPopup(() => {
          app.clickedFeature.featureVisibility = current.value;
          return answerDialogs(Promise.resolve(app.updateClickedFeatureVisibility(current.value)));
        });
        break;
      case 'labelText':
        await withPopup(() => {
          app.clickedFeature.labelText = current.text;
          return answerDialogs(Promise.resolve(app.updateClickedFeatureLabelText()));
        });
        break;
      case 'labelVisibility':
        await withPopup(() => {
          app.clickedFeature.labelVisibility = current.value;
          return answerDialogs(Promise.resolve(app.updateClickedFeatureLabelText()));
        });
        break;
      case 'legendAdd':
        // A Legend name of its own on one feature adds its row.
        await withPopup(() => {
          app.clickedFeature.legendName = current.caption;
          return answerDialogs(Promise.resolve(app.handleLegendNameCommit()));
        });
        break;
      case 'visibilityRuleAdd': {
        await app.addFeatureVisibilityRule();
        const index = app.featureVisibilityManualRules.length - 1;
        for (const [field, value] of Object.entries(current.fields)) await app.setFeatureVisibilityRuleField(index, field, value);
        break;
      }
      case 'visibilityRuleEdit':
        await app.setFeatureVisibilityRuleField(current.index, current.field, current.value);
        break;
      case 'visibilityRuleRemove':
        await app.removeFeatureVisibilityRule(current.index);
        break;
      case 'colorRuleAdd':
        Object.assign(app.newSpecRule, current.fields);
        await app.addSpecificRule();
        break;
      case 'colorRuleEdit':
        await answerDialogs(Promise.resolve(app.setSpecificRuleField(current.index, current.field, current.value)));
        break;
      case 'colorRuleRemove':
        await app.removeSpecificRule(current.index);
        break;
      case 'legendRename':
        await answerDialogs(Promise.resolve(app.renameLegendEntry(current.index, current.caption)));
        break;
      case 'legendColor':
        await app.updateLegendEntryColor(current.index, current.color);
        break;
      case 'legendSort':
        await app.sortLegendEntries(current.direction);
        break;
      case 'legendSortDefault':
        await app.sortLegendEntriesByDefault();
        break;
      case 'legendDelete':
        await answerDialogs(Promise.resolve(app.deleteLegendEntry(current.index)));
        break;
      case 'legendRestore':
        await answerDialogs(Promise.resolve(app.restoreAllDeletedLegendEntries()));
        break;
      case 'depthFile':
      case 'trackData': {
        // The Depth file of the shown drawing's first track, set (`depthFile`)
        // or removed (`trackData`) in a History step, as the track-data parity
        // cases change it.
        const file = current.kind === 'depthFile' ? new File([current.tsv], 'walk-depth.tsv', { type: 'text/tab-separated-values' }) : null;
        await history.runUndoable(file ? 'Change uploaded file' : 'Remove Depth', () => (app.mode === 'linear'
          ? app.setLinearDepthFile(app.linearSeqs[0], 0, file)
          : app.setCircularDepthFile(0, file)));
        break;
      }
      case 'palette':
        app.setPaletteInstantPreview(true);
        await answerDialogs(Promise.resolve(app.selectPalette(current.palette)));
        break;
      case 'reset':
        await app.resetSettings();
        break;
      case 'undo':
        await history.undo();
        break;
      case 'undoRedo':
        await history.undo();
        await history.redo();
        break;
      case 'reflow':
        app.autoLabelReflowEnabled = current.enabled;
        break;
      default:
        throw new Error(`unknown step kind ${current.kind}`);
    }
  };
  await run(step);
  return notes;
}, step);

// A step that never settles (for example, a dialog no answer covers) fails with
// the dialogs that are open, instead of waiting for the test deadline.
const STEP_TIMEOUT_MS = 180_000;
const applyStepWithin = async (page, step) => {
  let timer;
  const expired = new Promise((_, reject) => {
    timer = setTimeout(async () => {
      const open = await page.evaluate(() => Object.entries(window.__GBDRAW_APP__)
        .filter(([name, value]) => /Dialog$/.test(name) && value?.show).map(([name]) => name)).catch((error) => [`(page unreadable: ${error?.message || error})`]);
      const answered = await page.evaluate(() => window.__GBDRAW_WALK_NOTES__).catch(() => []);
      reject(new Error(`the step did not settle within ${STEP_TIMEOUT_MS / 1000} s; open dialogs: [${open}]; answered: [${answered}]`));
    }, STEP_TIMEOUT_MS);
  });
  try {
    return await Promise.race([applyStep(page, step), expired]);
  } finally {
    clearTimeout(timer);
  }
};

// Save Session, then Load that file into the page.
const saveAndLoad = async (page, path) => {
  await download(page, 'Save Session', path, { title: '' });
  const shown = () => page.evaluate(async () => {
    const { getCommittedSvgResultRuntimeIdentity } = await import('/gbdraw/web/js/services/svg-result-ingestion.js');
    const app = window.__GBDRAW_APP__;
    return app.sessionImportPending ? null : getCommittedSvgResultRuntimeIdentity(app.results[app.selectedResultIndex]);
  });
  const before = await shown();
  await page.locator('input[accept^=".json,"]').setInputFiles(path);
  await expect.poll(async () => { const now = await shown(); return now !== null && now !== before; }, { timeout: 180_000 }).toBe(true);
  await settleLive(page);
};

// The draft settings the fixtures set that apply on Generate (R1). Reset
// Settings returns them to their defaults in the draft only, so the Result
// keeps drawing them until Generate; the reset step writes them back in both
// drawings, so its comparison covers what Reset Settings shows live.
const GENERATE_SETTINGS = Object.freeze({
  circular: [['form', 'labels_mode'], ['form', 'multi_record_canvas'], ['adv', 'circular_grouping_intent']],
  linear: [['form', 'show_labels_linear']]
});
const resetKeepingGenerateSettings = async (page) => {
  const kept = await page.evaluate(async (settings) => {
    const { state } = await import('/gbdraw/web/js/state.js');
    return Object.fromEntries(Object.entries(settings).map(([drawingId, paths]) => [
      drawingId, paths.map(([group, key]) => [group, key, window.Vue.unref(state.drawingOf(drawingId)[group])[key]])
    ]));
  }, GENERATE_SETTINGS);
  const notes = await applyStepWithin(page, { kind: 'reset' });
  await settleLive(page);
  const shown = await page.evaluate(() => window.__GBDRAW_APP__.mode);
  const writeBack = (drawingId) => page.evaluate((values) => {
    const app = window.__GBDRAW_APP__;
    for (const [group, key, value] of values) app[group][key] = value;
  }, kept[drawingId]);
  await writeBack(shown);
  await switchDrawing(page, otherMode(shown));
  await writeBack(otherMode(shown));
  await switchDrawing(page, shown);
  return [...notes, 'the settings that apply on Generate written back'];
};

// One step: the kinds that need the test runner (a Result or drawing switch,
// Save and Load, Reset Settings, a Generate inside the step) run here, the
// others in the page.
const runStep = async (page, step, sessionPath) => {
  switch (step.kind) {
    case 'resultSwitch': await showResult(page, step.index); return [];
    case 'drawingSwitch': await switchDrawing(page, step.drawing); return [];
    case 'load': await saveAndLoad(page, sessionPath); return [];
    case 'reset': return resetKeepingGenerateSettings(page);
    case 'trackData':
      if (!step.prepare) return applyStepWithin(page, step);
      await applyStepWithin(page, step.prepare);
      await generate(page);
      return ['prepare depthFile, generated', ...await applyStepWithin(page, { kind: step.kind })];
    default: return applyStepWithin(page, step);
  }
};

// The popup state of a step: a feature popup opened and closed before it, as
// a reader looking at the feature does, which binds the labels.
const bindLabels = (page, key) => evaluateWithRetainedPromise(page, async (target) => {
  const app = window.__GBDRAW_APP__;
  const featureKey = (feature) => [feature.type, feature.start, feature.end, feature.locus_tag || '', feature.product || ''].join('|');
  await app.openFeatureEditorFromList(app.filteredFeatures.find((item) => featureKey(item) === target), null);
  await window.Vue.nextTick();
  app.clickedFeature = null;
}, key).then(() => settleLive(page));

const describeStep = (index, { kind, popup, reflow }, step, notes = []) => {
  const { picks, ...rest } = step || { kind };
  return `#${index} popup ${popup ? 'opened' : 'not'}, Auto Reflow ${reflow ? 'on' : 'off'}: ${JSON.stringify(rest)}`
    + `${notes.length ? ` (${notes.join('; ')})` : ''}`;
};

/**
 * The walk tests, one per fixture, for the spec to declare (`test(title, body)`).
 * @param {{ seed: number, steps: number }} options
 */
export const randomWalkTests = ({ seed, steps }) => {
  const schedule = planWalk(seed, steps);
  return FIXTURES.map((fixture, fixtureIndex) => ({
    title: `random walk (seed ${seed}): ${fixture.name}`,
    body: async ({ page }) => {
      test.setTimeout(300_000 + steps * 240_000);
      // One generator per fixture, so a fixture replays alone with the same seed.
      const random = createRandom(seed + Math.imul(fixtureIndex + 1, 0x9E3779B1));
      const sequence = [];
      const mismatches = [];
      const isolation = [];
      const report = (index, difference) => (
        `seed ${seed}, fixture "${fixture.name}", step ${index}\n`
        + `steps so far:\n${sequence.map((line) => `  ${line}`).join('\n')}\n${difference}`
      );
      // The drawings a step must leave as they were: the other drawing, or both for a switch.
      const checkIsolation = (index, phase, before, after, drawings) => {
        for (const drawingId of drawings) {
          const changed = changedSlices(before, after, drawingId);
          if (changed.length) isolation.push(report(index, `${phase} changed the ${drawingId} drawing's slices: ${changed.join(', ')}`));
        }
      };

      const reflow = random.next() < 0.5;
      sequence.push(`#init ${fixture.mode} shown, Auto Reflow ${reflow ? 'on' : 'off'}`);
      // eslint-disable-next-line no-console
      console.log(`seed ${seed}, fixture "${fixture.name}": kinds ${schedule[fixtureIndex].map(({ kind }) => kind).join(', ')}`);
      /** @type {Record<string, number>} */
      const floors = await openFixture(page, fixture, reflow);
      expect(await page.evaluate(() => window.__GBDRAW_APP__.results.length)).toBe(fixture.results === 'batch' ? 2 : 1);

      for (const [offset, planned] of schedule[fixtureIndex].entries()) {
        const index = offset + 1;
        let step = null;
        let notes = [];
        let before;
        let afterStep;
        try {
          await page.evaluate((enabled) => { window.__GBDRAW_APP__.autoLabelReflowEnabled = enabled; }, planned.reflow);
          await settleLive(page);
          let state = await readState(page);
          floors[state.drawing] ??= state.undoDepth;
          if (planned.popup) {
            await bindLabels(page, random.pick(state.features).key);
            state = await readState(page);
          }
          step = nextStep(random, state, planned.kind, { fixture, canUndo: state.canUndo && state.undoDepth > floors[state.drawing] });
          before = await sliceHashes(page);
          notes = await runStep(page, step, test.info().outputPath(`walk-step-${index}.gbdraw-session.json`));
          await settleLive(page);
          afterStep = await sliceHashes(page);
        } catch (error) {
          sequence.push(describeStep(index, planned, step));
          throw new Error(`${report(index, `the step failed: ${error?.message || error}`)}`
            + `${mismatches.length ? `\nmismatches before it:\n${mismatches.join('\n')}` : ''}`
            + `${isolation.length ? `\nisolation failures before it:\n${isolation.join('\n')}` : ''}`);
        }
        sequence.push(describeStep(index, planned, step, notes));
        // eslint-disable-next-line no-console
        console.log(`  ${sequence[sequence.length - 1]}`);
        if (!ALL_DRAWINGS_KINDS.includes(step.kind)) {
          checkIsolation(index, 'the step', before, afterStep, Object.keys(before.drawings)
            .filter((drawingId) => step.kind === 'drawingSwitch' || drawingId !== before.active));
        }
        if (step.kind === 'load') for (const drawingId of Object.keys(floors)) delete floors[drawingId];
        const known = KNOWN_MISMATCHES[step.kind];
        let differences = [];
        try {
          await expectLiveEqualsGenerate(page, { label: `step ${index} (${step.kind})` });
        } catch (error) {
          // expectLiveEqualsGenerate fails with the unexpected differences as `actual`.
          if (!Array.isArray(error?.matcherResult?.actual)) throw error;
          differences = error.matcherResult.actual;
          // Generate has drawn the draft; the walk goes on from that state.
        }
        if (known) {
          // A known mismatch is expected, as test.fail expects it: the step
          // fails when the mismatch no longer reproduces.
          if (!differences.some(known.matches)) mismatches.push(report(index, `${known.id} no longer reproduces; remove its mark`));
          else test.info().annotations.push({ type: 'known mismatch', description: `${known.id}: step ${index} (${step.kind})` });
          differences = differences.filter((line) => !known.matches(line));
        }
        if (differences.length) {
          mismatches.push(report(index, `differences:\n${differences.map((line) => `  ${line}`).join('\n')}`));
        }
        // Generate draws the shown drawing only.
        const afterGenerate = await sliceHashes(page);
        checkIsolation(index, 'Generate', afterStep, afterGenerate, Object.keys(afterStep.drawings)
          .filter((drawingId) => drawingId !== afterStep.active));
      }
      const failures = [
        ...(mismatches.length ? [`${mismatches.length} of ${steps} steps differ from Generate:\n\n${mismatches.join('\n\n')}`] : []),
        ...(isolation.length ? [`${isolation.length} changes of a drawing a step did not edit:\n\n${isolation.join('\n\n')}`] : [])
      ];
      expect(failures.length, `${failures.join('\n\n')}\n`).toBe(0);
    }
  }));
};
