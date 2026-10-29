// A seeded random walk over the live edit kinds (PD-OI-066, gbdraw/web/CLAUDE.md
// R3): after every step, the displayed Result must equal the Result the next
// Generate draws from the same draft (expectLiveEqualsGenerate, minus
// tests/web/contracts/live-generate-parity-allowed.json), and the step must
// leave the slices of the drawing it did not edit as they were (M6 isolation:
// the sha256 of each slice's `capture`, as tests/web/per-mode-drawings.test.mjs
// hashes them in node). Reset Settings and Session Load change every drawing
// by design (PD-OI-090), so the isolation check skips those steps.
//
// The kind schedule is the plan of tests/web/helpers/live-random-walk-plan.mjs
// (its targets: every ordered pair of kinds, and every kind in each popup and
// Auto Reflow state; at promotion also adjacent pairs and every fixture). Each
// walk of the plan is one test that opens its fixture (both drawings with a
// Result, so a drawing switch arrives at a Result) and draws the step
// arguments from its own generator, so a walk replays alone with the same
// seed. Two specs run it: the PR spec
// tests/web/live-generate-random-walk.playwright.spec.js (fixed seed) and the
// promotion spec tests/web/live-generate-random-walk.promotion.spec.js (the
// run ID as the seed).
//
// Every walk attaches the steps it ran as data (`walk-steps.json`). When a
// step differs from Generate or changes a drawing it did not edit, the walk
// replays the recorded steps from a fresh fixture to find a short sequence
// that still fails (tests/web/helpers/live-random-walk-shrink.mjs; budget
// GBDRAW_RANDOM_WALK_SHRINK_REPLAYS), reports it as the minimal reproduction,
// and attaches it as `walk-shrunk.json`. GBDRAW_RANDOM_WALK_REPLAY=<file>
// makes the promotion spec replay such a file.
import { expect, test } from '@playwright/test';
import { readFileSync, writeFileSync } from 'node:fs';
import { evaluateWithRetainedPromise } from './app-lifecycle.cjs';
import { showDrawing } from './drawing-bar.cjs';
import { expectLiveEqualsGenerate, settleLive, showResult } from './live-generate-parity.cjs';
import { download } from './mode-transition.cjs';
import { SINGLE_FIXTURE, generate, open } from './live-generate-parity-steps.cjs';
import { BARRIER_KINDS, FIXTURES, createRandom, planWalk } from './live-random-walk-plan.mjs';
import { shrinkWalk } from './live-random-walk-shrink.mjs';

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

// The records of each file, for a Depth file over them.
const FILE_RECORDS = { single: ['FORCEDLBL'], batch: ['TESTA', 'TESTB'] };

const COLORS = ['#c83366', '#123456', '#2a9d8f', '#e63946', '#7b2cbf', '#f4a261', '#264653'];
const LABEL_TEXTS = [
  'alpha protein, edited live',
  'x',
  'walk label',
  'a label text far too long to stay inside its own arc, drawn outside with leaders'
];
// The label text short enough to be drawn inside the smallest feature. Two
// features with the same span are drawn on top of each other, so this text on
// both draws two labels on one place, and Generate rightly refuses the draft
// (PD-OI-097, LP-D52: "Labels do not fit", fix "separate overlapping
// features"; OV-2195). The walk gives a feature with a same-span twin the
// outside text 'walk label' instead; the random draws stay as they were.
const INSIDE_LABEL_TEXT = 'x';
const CAPTIONS = ['Walk caption', 'transfer RNA', 'Renamed group', 'CDS'];
const ADDED_CAPTIONS = ['Walk group', 'Added row'];
const escapeRegex = (text) => text.replace(/[.*+?^${}()|[\]\\]/g, '\\$&');
const otherMode = (mode) => (mode === 'circular' ? 'linear' : 'circular');

// One drawing switch as a user makes it (a tab, or the Linear offer when the project holds no Linear drawing), settled.
const switchDrawing = async (page, kind) => {
  await settleLive(page);
  await showDrawing(page, { kind });
  await expect.poll(() => page.evaluate(async () => (await import('/gbdraw/web/js/state.js')).state.drawingOf('circular').circularRecordDiscovery.status))
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
    key: featureKey(feature),
    span: `${feature.start}|${feature.end}`,
    type: String(feature.type || ''),
    tag: String(feature.locus_tag || ''),
    product: String(feature.product || '')
  })).sort((left, right) => (left.key < right.key ? -1 : left.key > right.key ? 1 : 0));
  return {
    drawing: state.activeDrawing().id,
    mode: state.activeDrawing().mode,
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
    case 'labelText': {
      const target = random.pick(features);
      const text = random.pick(LABEL_TEXTS);
      const stacked = features.some((other) => other.key !== target.key && other.span === target.span);
      return { kind, feature: target.key, text: text === INSIDE_LABEL_TEXT && stacked ? 'walk label' : text, picks };
    }
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
    case 'drawingSwitch': return { kind, mode: otherMode(state.mode) };
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
    case 'drawingSwitch': await switchDrawing(page, step.mode); return [];
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

/**
 * A step as a walk ran it: the data a replay applies again.
 * @typedef {object} RecordedStep
 * @property {number} index the step's position in its run (1-based)
 * @property {string} kind
 * @property {boolean} popup a feature popup was opened and closed before the step
 * @property {boolean} reflow Auto Reflow before the step
 * @property {string} drawing the drawing shown before the step
 * @property {string | null} popupFeature the key of the feature whose popup was opened
 * @property {Record<string, unknown> & { kind: string }} step the concrete step, with its `picks` and `prepare`
 * @property {string[]} notes what the run answered and did besides the step
 */
/**
 * The steps of one run, as `walk-steps.json` and `walk-shrunk.json` hold them.
 * @typedef {{ seed: number, title: string, fixture: string, reflow: boolean, steps: RecordedStep[], shrink?: object }} WalkRecord
 */
/** @typedef {{ phase: string, drawing: string, slices: string[] }} IsolationChange */
/** @typedef {{ index: number, kind: string, differences: string[], isolation: IsolationChange[] }} StepResult */
/**
 * Where the steps of a run come from: the walk draws them, a replay reads
 * them. Each returns null when the step cannot apply in `state`.
 * @typedef {{
 *   popupFeature: (offset: number, state: Awaited<ReturnType<typeof readState>>) => string | null,
 *   step: (offset: number, state: Awaited<ReturnType<typeof readState>>, context: { fixture: object, canUndo: boolean }) => { step: RecordedStep['step'] | null, reason: string }
 * }} StepSource
 */

// The walk draws the popup feature and the step arguments from its generator.
/** @returns {StepSource} */
const drawnSource = (random, planned) => ({
  popupFeature: (offset, state) => random.pick(state.features).key,
  step: (offset, state, context) => ({ step: nextStep(random, state, planned[offset].kind, context), reason: '' })
});

// Why a recorded step cannot apply in `state` (its subject is gone), or ''.
const whyNotApplicable = (step, state, { canUndo }) => {
  if (step.prepare) {
    const reason = whyNotApplicable(step.prepare, state, { canUndo });
    if (reason) return `its prepare step: ${reason}`;
  }
  const feature = () => (state.features.some(({ key }) => key === step.feature) ? '' : `no feature ${step.feature}`);
  const legendRow = () => (step.index < state.legend.length ? '' : `no Legend row ${step.index}`);
  switch (step.kind) {
    case 'fill': case 'featureVisibility': case 'labelText': case 'labelVisibility': case 'legendAdd': return feature();
    case 'stroke': return step.target === 'popup' ? feature() : legendRow();
    case 'legendRename': case 'legendColor': case 'legendDelete': return legendRow();
    case 'visibilityRuleEdit': case 'visibilityRuleRemove':
      return step.prepare || step.index < state.visibilityRules ? '' : `no visibility rule ${step.index}`;
    case 'colorRuleEdit': case 'colorRuleRemove': return step.prepare || step.index < state.colorRules ? '' : `no color rule ${step.index}`;
    case 'legendRestore': return step.prepare || state.deletedLegend ? '' : 'no deleted Legend row';
    case 'trackData': return step.prepare || state.depth ? '' : 'no Depth file';
    case 'palette': return state.palettes.includes(step.palette) ? '' : `palette ${step.palette} is selected`;
    case 'undo': case 'undoRedo': return step.prepare || canUndo ? '' : 'Undo is at its floor';
    case 'resultSwitch':
      return step.index < state.resultCount && step.index !== state.resultIndex ? '' : `Result ${step.index} is shown or absent`;
    case 'drawingSwitch': return step.mode !== state.mode ? '' : `the ${step.mode} drawing is shown`;
    default: return '';
  }
};

// A replay applies the recorded steps, drawing nothing new.
/** @param {RecordedStep[]} recorded @returns {StepSource} */
const recordedSource = (recorded) => ({
  popupFeature: (offset, state) => (state.features.some(({ key }) => key === recorded[offset].popupFeature) ? recorded[offset].popupFeature : null),
  step: (offset, state, context) => {
    const { drawing, step } = recorded[offset];
    const reason = drawing === state.drawing ? whyNotApplicable(step, state, context) : `the ${state.drawing} drawing is shown, not ${drawing}`;
    return reason ? { step: null, reason } : { step, reason: '' };
  }
});

// The slices of the drawings a phase must leave as they were that changed.
/** @returns {IsolationChange[]} */
const isolationChanges = (phase, before, after, drawings) => drawings
  .map((drawing) => ({ phase, drawing, slices: changedSlices(before, after, drawing) }))
  .filter(({ slices }) => slices.length);

/**
 * Opens the fixture and runs the planned steps from `source`, with the oracle
 * (one Generate) after each step. Returns the recorded steps, each step's
 * differences and isolation changes, and the step error or the step that could
 * not apply, if any; it stops there.
 * @param {import('@playwright/test').Page} page
 * @param {{ fixture: typeof FIXTURES[number], reflow: boolean, planned: readonly { kind: string, popup: boolean, reflow: boolean }[],
 *   source: StepSource, sessionPath: (index: number) => string, log?: (record: RecordedStep) => void }} options
 * @returns {Promise<{ records: RecordedStep[], results: StepResult[], error?: { index: number, message: string },
 *   inapplicable?: { index: number, reason: string } }>}
 */
const runWalk = async (page, { fixture, reflow, planned, source, sessionPath, log = () => {} }) => {
  /** @type {RecordedStep[]} */
  const records = [];
  /** @type {StepResult[]} */
  const results = [];
  /** @type {Record<string, number>} */
  const floors = await openFixture(page, fixture, reflow);
  expect(await page.evaluate(() => window.__GBDRAW_APP__.results.length)).toBe(fixture.results === 'batch' ? 2 : 1);
  for (const [offset, { kind, popup, reflow: stepReflow }] of planned.entries()) {
    const index = offset + 1;
    /** @type {RecordedStep} */
    const record = { index, kind, popup, reflow: stepReflow, drawing: '', popupFeature: null, step: { kind }, notes: [] };
    let before;
    let afterStep;
    try {
      await page.evaluate((enabled) => { window.__GBDRAW_APP__.autoLabelReflowEnabled = enabled; }, stepReflow);
      await settleLive(page);
      let state = await readState(page);
      floors[state.drawing] ??= state.undoDepth;
      record.drawing = state.drawing;
      if (popup) {
        record.popupFeature = source.popupFeature(offset, state);
        if (record.popupFeature === null) return { records, results, inapplicable: { index, reason: 'the popup feature is gone' } };
        await bindLabels(page, record.popupFeature);
        state = await readState(page);
      }
      const { step, reason } = source.step(offset, state, { fixture, canUndo: state.canUndo && state.undoDepth > floors[state.drawing] });
      if (!step) return { records, results, inapplicable: { index, reason } };
      record.step = step;
      before = await sliceHashes(page);
      record.notes = await runStep(page, step, sessionPath(index));
      await settleLive(page);
      afterStep = await sliceHashes(page);
    } catch (error) {
      records.push(record);
      return { records, results, error: { index, message: error?.message || String(error) } };
    }
    records.push(record);
    log(record);
    const isolation = BARRIER_KINDS.includes(record.step.kind) ? [] : isolationChanges('the step', before, afterStep,
      Object.keys(before.drawings).filter((drawingId) => record.step.kind === 'drawingSwitch' || drawingId !== before.active));
    if (record.step.kind === 'load') for (const drawingId of Object.keys(floors)) delete floors[drawingId];
    let differences = [];
    try {
      await expectLiveEqualsGenerate(page, { label: `step ${index} (${record.step.kind})` });
    } catch (error) {
      // expectLiveEqualsGenerate fails with the unexpected differences as `actual`.
      if (!Array.isArray(error?.matcherResult?.actual)) {
        return { records, results, error: { index, message: `the comparison failed: ${error?.message || error}` } };
      }
      differences = error.matcherResult.actual;
      // Generate has drawn the draft; the walk goes on from that state.
    }
    // Generate draws the shown drawing only.
    const afterGenerate = await sliceHashes(page);
    isolation.push(...isolationChanges('Generate', afterStep, afterGenerate,
      Object.keys(afterStep.drawings).filter((drawingId) => drawingId !== afterStep.active)));
    results.push({ index, kind: record.step.kind, differences, isolation });
  }
  return { records, results };
};

/** @param {RecordedStep} record */
const describeRecord = ({ index, popup, reflow, step, notes }) => {
  const { picks, ...rest } = step;
  return `#${index} popup ${popup ? 'opened' : 'not'}, Auto Reflow ${reflow ? 'on' : 'off'}: ${JSON.stringify(rest)}`
    + `${notes.length ? ` (${notes.join('; ')})` : ''}`;
};

// A step result without the known mismatches of its kind; a known mismatch
// that no longer reproduces is a difference of its own.
/** @param {StepResult} result */
const unknownDifferences = (result) => {
  const known = KNOWN_MISMATCHES[result.kind];
  if (!known) return result.differences;
  if (!result.differences.some(known.matches)) return [`${known.id} no longer reproduces; remove its mark`];
  return result.differences.filter((line) => !known.matches(line));
};
// The lines that identify a failure: a candidate reproduces it when it shows one of them.
/** @param {StepResult} result */
const failureLines = (result) => [
  ...unknownDifferences(result).map((line) => `difference: ${line}`),
  ...result.isolation.flatMap(({ phase, drawing, slices }) => slices.map((slice) => `isolation: ${phase} changed ${drawing}.${slice}`))
];

const attachJson = async (name, value) => {
  const path = test.info().outputPath(name);
  writeFileSync(path, `${JSON.stringify(value, null, 2)}\n`);
  await test.info().attach(name, { path, contentType: 'application/json' });
  return path;
};

/** The replay budget of the shrink: GBDRAW_RANDOM_WALK_SHRINK_REPLAYS, or `fallback`. */
export const shrinkReplaysFromEnv = (fallback) => {
  const value = Number.parseInt(process.env.GBDRAW_RANDOM_WALK_SHRINK_REPLAYS ?? `${fallback}`, 10);
  if (!Number.isInteger(value) || value < 0) throw new Error('GBDRAW_RANDOM_WALK_SHRINK_REPLAYS must be a non-negative integer');
  return value;
};

// The test time a run of `steps` steps may take, and the shrink replays after it.
const STEP_ALLOWANCE_MS = 240_000;
const OPEN_ALLOWANCE_MS = 300_000;
const runAllowance = (steps, shrinkReplays) => (OPEN_ALLOWANCE_MS + steps * STEP_ALLOWANCE_MS) * (1 + shrinkReplays);

/**
 * Fails the test with the run's step error, or with its differences and
 * isolation changes after shrinking the first failing step to a minimal
 * reproduction (`walk-shrunk.json`).
 * @param {import('@playwright/test').Page} page
 * @param {WalkRecord} walk
 * @param {Awaited<ReturnType<typeof runWalk>>} run
 * @param {{ header: string, shrinkReplays: number }} options
 */
const reportRun = async (page, walk, run, { header, shrinkReplays }) => {
  const sequence = [`#init ${FIXTURES.find(({ name }) => name === walk.fixture).mode} shown, Auto Reflow ${walk.reflow ? 'on' : 'off'}`,
    ...run.records.map(describeRecord)];
  const report = (index, difference) => `${header}, step ${index}\nsteps so far:\n`
    + `${sequence.slice(0, index + 1).map((line) => `  ${line}`).join('\n')}\n${difference}`;
  const failing = run.results.map((result) => ({ result, lines: failureLines(result) })).filter(({ lines }) => lines.length);
  const mismatches = failing.filter(({ lines }) => lines.some((line) => line.startsWith('difference: ')))
    .map(({ result }) => report(result.index, `differences:\n${unknownDifferences(result).map((line) => `  ${line}`).join('\n')}`));
  const isolation = failing.flatMap(({ result }) => result.isolation.map(({ phase, drawing, slices }) => (
    report(result.index, `${phase} changed the ${drawing} drawing's slices: ${slices.join(', ')}`))));
  if (run.error || run.inapplicable) {
    const { index } = run.error || run.inapplicable;
    const cause = run.error ? `the step failed: ${run.error.message}` : `the recorded step cannot apply: ${run.inapplicable.reason}`;
    throw new Error(`${report(index, cause)}`
      + `${mismatches.length ? `\nmismatches before it:\n${mismatches.join('\n')}` : ''}`
      + `${isolation.length ? `\nisolation failures before it:\n${isolation.join('\n')}` : ''}`);
  }
  for (const result of run.results) {
    const known = KNOWN_MISMATCHES[result.kind];
    if (known && result.differences.some(known.matches)) {
      test.info().annotations.push({ type: 'known mismatch', description: `${known.id}: step ${result.index} (${result.kind})` });
    }
  }
  if (!failing.length) return;

  // The first failing step, cut down to the steps it needs.
  const { result: first, lines } = failing[0];
  const original = new Set(lines);
  const fixture = FIXTURES.find(({ name }) => name === walk.fixture);
  let replay = 0;
  const reproduces = async (indices) => {
    replay += 1;
    const steps = indices.map((index) => walk.steps[index - 1]);
    try {
      const rerun = await runWalk(page, {
        fixture, reflow: walk.reflow, planned: steps, source: recordedSource(steps),
        sessionPath: (index) => test.info().outputPath(`shrink-${replay}-step-${index}.gbdraw-session.json`)
      });
      const last = rerun.results[steps.length - 1];
      const reproduced = Boolean(!rerun.error && !rerun.inapplicable && last && failureLines(last).some((line) => original.has(line)));
      // eslint-disable-next-line no-console
      console.log(`  shrink replay ${replay} [${indices}]: ${reproduced ? 'reproduces'
        : rerun.inapplicable ? `step ${rerun.inapplicable.index} cannot apply (${rerun.inapplicable.reason})`
          : rerun.error ? `step ${rerun.error.index} failed` : 'does not reproduce'}`);
      return reproduced;
    } catch (error) {
      // eslint-disable-next-line no-console
      console.log(`  shrink replay ${replay} [${indices}] failed to run: ${error?.message || error}`);
      return false;
    }
  };
  let shrinkReport = '';
  if (shrinkReplays > 0) {
    const shrunk = await shrinkWalk({
      prefix: walk.steps.slice(0, first.index - 1).map(({ index }) => index), failing: first.index, budget: shrinkReplays, reproduces
    });
    const steps = shrunk.indices.map((index) => walk.steps[index - 1]);
    const path = await attachJson('walk-shrunk.json', {
      ...walk, steps, shrink: { failingStep: first.index, replays: shrunk.replays, budget: shrinkReplays, finished: shrunk.finished, confirmed: shrunk.confirmed }
    });
    shrinkReport = `minimal reproduction of step ${first.index}: ${steps.length} of ${first.index} steps, ${shrunk.replays} replays, `
      + `${shrunk.finished ? 'search finished' : `search stopped at the budget of ${shrinkReplays} replays`}, `
      + `${shrunk.confirmed ? 'reproduced in a replay' : 'no shorter sequence reproduced in a replay'}:\n`
      + `${steps.map((step) => `  ${JSON.stringify(step)}`).join('\n')}\n`
      + `replay: GBDRAW_RANDOM_WALK_REPLAY=${path} npx playwright test --config=tests/web/playwright/promotion.config.js --workers=1\n\n`;
  }
  const failures = [
    ...(mismatches.length ? [`${mismatches.length} of ${run.results.length} steps differ from Generate:\n\n${mismatches.join('\n\n')}`] : []),
    ...(isolation.length ? [`${isolation.length} changes of a drawing a step did not edit:\n\n${isolation.join('\n\n')}`] : [])
  ];
  expect(failures.length, `${shrinkReport}${failures.join('\n\n')}\n`).toBe(0);
};

/**
 * The walk tests of a plan, one per walk, for the spec to declare
 * (`test(title, body)`).
 * @param {{ seed: number, level: import('./live-random-walk-plan.mjs').PlanLevel, walkLength?: number,
 *   required?: readonly import('./live-random-walk-plan.mjs').RequiredStep[], shrinkReplays: number }} options
 */
export const randomWalkTests = ({ seed, level, walkLength, required = [], shrinkReplays }) => {
  const { walks } = planWalk(seed, { level, walkLength, required });
  return walks.map(({ fixtureIndex, steps }, walkIndex) => {
    const fixture = FIXTURES[fixtureIndex];
    const title = `random walk (seed ${seed}): ${fixture.name}, walk ${walkIndex + 1} of ${walks.length}`;
    return {
      title,
      body: async ({ page }) => {
        test.setTimeout(runAllowance(steps.length, shrinkReplays));
        // One generator per walk, so a walk replays alone with the same seed.
        const random = createRandom(seed + Math.imul(walkIndex + 1, 0x9E3779B1));
        const reflow = random.next() < 0.5;
        // eslint-disable-next-line no-console
        console.log(`seed ${seed}, ${fixture.name}, walk ${walkIndex + 1}: kinds ${steps.map(({ kind }) => kind).join(', ')}`);
        const run = await runWalk(page, {
          fixture, reflow, planned: steps, source: drawnSource(random, steps),
          sessionPath: (index) => test.info().outputPath(`walk-step-${index}.gbdraw-session.json`),
          // eslint-disable-next-line no-console
          log: (record) => console.log(`  ${describeRecord(record)}`)
        });
        const walk = { seed, title, fixture: fixture.name, reflow, steps: run.records };
        await attachJson('walk-steps.json', walk);
        await reportRun(page, walk, run, { header: `seed ${seed}, fixture "${fixture.name}", walk ${walkIndex + 1}`, shrinkReplays });
      }
    };
  });
};

/** @param {string} path @returns {WalkRecord} */
const readWalkRecord = (path) => {
  const walk = JSON.parse(readFileSync(path, 'utf8'));
  if (!FIXTURES.some(({ name }) => name === walk.fixture) || !Array.isArray(walk.steps)) {
    throw new Error(`${path} is not a walk-steps.json or walk-shrunk.json file`);
  }
  return walk;
};
// Replays a recorded run: the same fixture, the same concrete steps, the oracle after every step.
const replayRun = (page, walk, log = () => {}) => runWalk(page, {
  fixture: FIXTURES.find(({ name }) => name === walk.fixture), reflow: walk.reflow, planned: walk.steps, source: recordedSource(walk.steps),
  sessionPath: (index) => test.info().outputPath(`replay-step-${index}.gbdraw-session.json`), log
});

/**
 * The test that replays a `walk-steps.json` or `walk-shrunk.json` file
 * (GBDRAW_RANDOM_WALK_REPLAY in the promotion spec).
 * @param {string} path
 * @param {{ shrinkReplays: number }} options
 */
export const replayWalkTest = (path, { shrinkReplays }) => {
  const walk = readWalkRecord(path);
  return {
    title: `random walk replay of ${path}: ${walk.fixture}, ${walk.steps.length} steps`,
    body: async ({ page }) => {
      test.setTimeout(runAllowance(walk.steps.length, shrinkReplays));
      // eslint-disable-next-line no-console
      const run = await replayRun(page, walk, (record) => console.log(`  ${describeRecord(record)}`));
      const replayed = { ...walk, steps: run.records };
      delete replayed.shrink;
      await attachJson('walk-steps.json', replayed);
      await reportRun(page, replayed, run, { header: `replay of ${path}`, shrinkReplays });
    }
  };
};

// The steps of the record-and-replay test: a popup edit with its dialogs, a
// rule, and Undo, on the Circular fixture.
const RECORD_REPLAY_STEPS = Object.freeze([
  { kind: 'fill', popup: true, reflow: false },
  { kind: 'colorRuleAdd', popup: false, reflow: true },
  { kind: 'undo', popup: false, reflow: false }
]);

/**
 * The test that keeps the replay path working on every pull request: it runs a
 * short walk, attaches its steps, replays them from the attached file, and
 * expects the same steps and no failure in either run.
 * @param {{ seed: number }} options
 */
export const recordReplayTest = ({ seed }) => ({
  title: `random walk (seed ${seed}): a recorded walk replays from its walk-steps.json`,
  body: async ({ page }) => {
    test.setTimeout(2 * runAllowance(RECORD_REPLAY_STEPS.length, 0));
    const fixture = FIXTURES[0];
    const random = createRandom(seed);
    const reflow = random.next() < 0.5;
    const run = await runWalk(page, {
      fixture, reflow, planned: RECORD_REPLAY_STEPS, source: drawnSource(random, RECORD_REPLAY_STEPS),
      sessionPath: (index) => test.info().outputPath(`walk-step-${index}.gbdraw-session.json`)
    });
    const walk = { seed, title: test.info().title, fixture: fixture.name, reflow, steps: run.records };
    const path = await attachJson('walk-steps.json', walk);
    await reportRun(page, walk, run, { header: `seed ${seed}, recorded walk`, shrinkReplays: 0 });
    const recorded = readWalkRecord(path);
    const replay = await replayRun(page, recorded);
    await reportRun(page, { ...recorded, steps: replay.records }, replay, { header: `replay of ${path}`, shrinkReplays: 0 });
    expect(replay.records.map(describeRecord)).toEqual(run.records.map(describeRecord));
  }
});
