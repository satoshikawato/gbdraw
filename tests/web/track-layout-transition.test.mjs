// R-3, Q3 (Owner, 2026-10-04), R10, R11: every edit of a draft feature-slot
// input runs through the feature placement owner's transition, which asks
// before it leaves a lane Feature placement undrawable. The guard lists every
// path that writes an input of draftPlacementTargets (form.track_type,
// form.linear_track_layout, form.separate_strands, adv.<mode>_track_slots*,
// and a stack row's renderer, enabled, side, or lane) and fails on a writer
// outside the transition that is not a listed restore or reconcile.
import assert from 'node:assert/strict';
import { readFileSync, readdirSync, statSync } from 'node:fs';
import { join, relative } from 'node:path';
import test from 'node:test';
import { createDefaultAdv, createDefaultForm } from '../../gbdraw/web/js/services/session-active-config-contract.js';
import {
  createFeaturePlacementActions, draftPlacementTargets, restoreTrackLayout, saveTrackLayout
} from '../../gbdraw/web/js/app/feature-editor/placement-actions.js';
import { createCircularTrackSlotEditor } from '../../gbdraw/web/js/app/circular-track-slots.js';
import { createLinearTrackSlotEditor } from '../../gbdraw/web/js/app/linear-track-slots.js';
import { trackLayoutWithVisibility } from '../../gbdraw/web/js/app/track-slot-edits.js';
import { withDrawings, withModeDrawings } from './helpers/drawing-state.mjs';
import { SOURCE_LOCALE, setUiLocale } from '../../gbdraw/web/js/utils/i18n.js';

const WEB = 'gbdraw/web';
const MODES = ['circular', 'linear'];
// A drawing's lane row; its mode is the drawing's (PR-1).
const lane = (name, side) => ({ recordKey: `${name}-record`, biologicalFeatureId: `${name}-cds`,
  placement: { kind: 'lane', side, level: 1 } });
const draftState = (mode = 'circular') => {
  const adv = createDefaultAdv(mode);
  // The simple controls differ from the default stacks, so Reset changes them.
  return { mode: { value: mode }, form: { ...createDefaultForm(), track_type: 'middle', linear_track_layout: 'above',
    separate_strands: false },
    adv: { ...adv, circular_track_slots_enabled: true, linear_track_slots_enabled: true }, featurePlacementOverrides: {},
    files: { c_depth: [], c_conservation_blasts: [], c_conservation_fastas: [] },
    circularConservation: { enabled: false, source: 'upload', labels: '', series: [] },
    annotationSets: [], linearSeqs: [], circularRecordList: { value: [] }, selectedResultIndex: { value: 0 },
    trackSlotResolvedGeometry: { value: null }, currentColors: { value: {} }, paletteDefinitions: { value: {} },
    selectedPalette: { value: '' } };
};
// The resolved inputs of draftPlacementTargets, without geometry or labels.
const project = ({ form, adv }) => JSON.stringify([form.track_type, form.linear_track_layout, form.separate_strands,
  ...MODES.map((mode) => [adv[`${mode}_track_slots_enabled`], adv[`${mode}_track_slots_axis_index`],
    adv[`${mode}_track_slots`].map((slot) => [slot.renderer, slot.enabled, slot.side,
      slot.params?.lane_direction ?? slot.params?.lanes ?? null])])]);
const everything = ({ form, adv }) => JSON.stringify({ form, adv });

test('a stack edit that drops a lane asks, Cancel keeps the draft, and Reset is one step (Q3)', async () => {
  const state = draftState('circular');
  // The Linear drawing holds its own lane; a Circular edit never reaches it (R2).
  const linear = draftState('linear');
  linear.featurePlacementOverrides.l = lane('linear', 'above');
  const linearBefore = everything(linear);
  const store = withModeDrawings(state, { circular: state, linear });
  const steps = [];
  const placement = createFeaturePlacementActions({ state: store, getCommittedRequest: () => null, isCurrentFeature: () => true,
    runUndoable: async (label, fn) => { steps.push(label); fn(); } });
  const editor = createCircularTrackSlotEditor({ state: store, changeTrackLayout: placement.changeTrackLayout });
  editor.normalizeCircularTrackSlots(state);
  Object.assign(state.featurePlacementOverrides, { c: lane('circular', 'outward'),
    m: { recordKey: 'circular-record', biologicalFeatureId: 'main-cds', placement: { kind: 'main' } } });
  const features = () => state.adv.circular_track_slots.findIndex((slot) => slot.renderer === 'features');
  const before = everything(state);
  const click = () => {
    const event = new Event('click');
    Object.defineProperty(event, 'currentTarget', { value: { tagName: 'BUTTON', textContent: '', title: 'Move outside Axis' } });
    return event;
  };
  assert.equal(editor.moveCircularTrackSlotOutside(features(), click()), false);
  assert.deepEqual({ ...placement.layoutChange },
    { open: true, count: 1, setting: 'Move outside Axis', value: '' });
  assert.equal(everything(state), before);
  await placement.resolveLayoutChange('cancel');
  assert.equal(placement.layoutChange.open, false);
  assert.equal(everything(state), before);
  assert.deepEqual(steps, []);

  editor.moveCircularTrackSlotOutside(features(), click());
  await placement.resolveLayoutChange('reset');
  assert.equal(state.adv.circular_track_slots[features()].params.lane_direction, 'outside');
  // Only the lost row goes; the Main row stays.
  assert.deepEqual(Object.keys(state.featurePlacementOverrides), ['m']);
  assert.deepEqual(steps, ['Change setting and reset Feature placements']);

  // Separate Strands of the Circular drawing is no input of a Circular lane,
  // and the Linear drawing has its own: the edit applies now and the control's
  // adapter records it.
  const strands = { type: 'checkbox', checked: true, getAttribute: () => 'Separate Strands' };
  assert.equal(placement.changeLayoutSetting({ target: strands }, 'separate_strands'), undefined);
  assert.equal(state.form.separate_strands, true);
  assert.equal(placement.layoutChange.open, false);
  assert.equal(steps.length, 1);
  assert.equal(everything(linear), linearBefore);
});

// OV-1418: the dialog names a checkbox's new value in the shown language.
test('a checkbox edit that drops a lane names its new value in the shown language', async () => {
  const linear = draftState('linear');
  linear.featurePlacementOverrides.l = lane('linear', 'above');
  const store = withModeDrawings(linear, { circular: draftState('circular'), linear });
  const placement = createFeaturePlacementActions({ state: store, getCommittedRequest: () => null, isCurrentFeature: () => true,
    runUndoable: async (label, fn) => { fn(); } });
  try {
    await setUiLocale('ja', async () => ({ default: { On: 'オン', Off: 'オフ' } }));
    const strands = { type: 'checkbox', checked: true, getAttribute: () => 'Separate Strands' };
    placement.changeLayoutSetting({ target: strands }, 'separate_strands');
    assert.deepEqual({ ...placement.layoutChange }, { open: true, count: 1, setting: 'Separate Strands', value: 'オン' });
    await placement.resolveLayoutChange('cancel');
  } finally {
    await setUiLocale(SOURCE_LOCALE);
  }
});

// Each editor export runs with every call below, on a stack whose feature row
// is on the Axis and on one whose feature row is off it.
const slotOf = (state, mode, renderer = 'features') => state.adv[`${mode}_track_slots`]
  .find((slot) => slot.renderer === renderer);
const CALLS = (state, mode) => {
  const slots = state.adv[`${mode}_track_slots`];
  const feature = slots.indexOf(slotOf(state, mode));
  return [[], [true], [false], ['spacer'], ['tuckin'],
    ...slots.flatMap((_, index) => [[index], [index, index - 1], [index, index + 1]]),
    ...[slots[feature], slots.at(-1)].flatMap((slot) => [[slot], [slot, false],
      [slot, 'inside'], [slot, 'outside'], [slot, 'overlay'], [slot, 'above'], [slot, 'below'], [slot, 'spacer'],
      [slot, 'width', '5'], [slot, 'positive_color', '#000000']])];
};
// Exports that may write the inputs outside the transition because they keep
// every mode's lanes: R10 reconciles (managed Depth and Comparison rows, the
// suppress controls) and in-place normalization.
const RECONCILE = new Set([
  'normalizeCircularTrackSlots', 'syncCircularConservationSlots', 'changeCircularDepthSources',
  'setCircularGcSuppressed', 'setCircularSkewSuppressed',
  'normalizeLinearTrackSlots', 'changeLinearDepthSources'
]);
const EDITORS = { circular: createCircularTrackSlotEditor, linear: createLinearTrackSlotEditor };
// The exports an editor routes through its `changeTrackLayout` port: each asks
// the port first, so a port that never applies the edit names them.
const transitionalExports = (mode) => {
  let asked = false;
  const editor = EDITORS[mode]({ state: draftState(mode), changeTrackLayout: () => { asked = true; } });
  return Object.keys(editor).filter((name) => {
    if (typeof editor[name] !== 'function') return false;
    asked = false;
    try { editor[name](); } catch { /* an unfit call */ }
    return asked;
  });
};
const STACK = {
  circular: { normalize: 'normalizeCircularTrackSlots', offAxis: 'moveCircularTrackSlotOutside' },
  linear: { normalize: 'normalizeLinearTrackSlots', offAxis: 'moveLinearTrackSlotAbove' }
};

for (const mode of MODES) {
  test(`${mode} stack editor writes the feature-slot inputs only through the transition (R10)`, () => {
    const fixture = (offAxis) => {
      // The fixture is the drawing the editor reads, so it names its mode.
      const state = Object.assign(draftState(mode), { mode });
      state.adv.linear_track_slots.push(
        { id: 'gc_content', renderer: 'dinucleotide_content', enabled: true, side: 'below', params: { nt: 'GC' } },
        { id: 'gap', renderer: 'spacer', enabled: true, side: 'below', height: '12px', params: {} });
      const windows = [];
      // The port applies every edit and records the inputs around it.
      const editor = EDITORS[mode]({ state: withDrawings(state), changeTrackLayout: (apply) => {
        const entry = project(state);
        try { return apply(); } finally { windows.push([entry, project(state)]); }
      } });
      editor[STACK[mode].normalize](state);
      if (offAxis) editor[STACK[mode].offAxis](state.adv[`${mode}_track_slots`].indexOf(slotOf(state, mode)));
      // The simple controls now differ from the stack, so Reset changes it.
      Object.assign(state.form, { track_type: 'tuckin', linear_track_layout: 'below' });
      windows.length = 0;
      return { state, editor, windows };
    };
    const transitional = transitionalExports(mode);
    assert.ok(transitional.length >= 10, 'the stack edits run through changeTrackLayout');
    const changedBy = new Set();
    for (const offAxis of [false, true]) {
      const { editor: probe } = fixture(offAxis);
      for (const name of Object.keys(probe).filter((key) => typeof probe[key] === 'function')) {
        for (let call = 0; call < CALLS(fixture(offAxis).state, mode).length; call += 1) {
          const { state, editor, windows } = fixture(offAxis);
          const args = CALLS(state, mode)[call];
          const lanes = laneSidesOf(state);
          const before = project(state);
          const whole = everything(state);
          const saved = saveTrackLayout(state);
          const label = `${mode} ${name} call ${call}${offAxis ? ' (feature row off the Axis)' : ''}`;
          try { editor[name](...args); } catch { /* an unfit call */ }
          if (transitional.includes(name)) {
            // Cancel restores every field the edit wrote (restoreTrackLayout).
            if (everything(state) !== whole) changedBy.add(name);
            restoreTrackLayout(state, saved);
            assert.equal(everything(state), whole, `${label}: Cancel leaves a field the edit wrote`);
          } else if (RECONCILE.has(name)) {
            assert.deepEqual(laneSidesOf(state), lanes, `${label} changes the lanes`);
          } else {
            const chain = [before, ...windows.flat(), project(state)];
            for (let index = 0; index < chain.length; index += 2) {
              assert.equal(chain[index + 1], chain[index], `${label} writes a feature-slot input outside the transition`);
            }
          }
        }
      }
    }
    // Every transition edit was exercised by at least one call.
    assert.deepEqual(transitional.filter((name) => !changedBy.has(name)), []);
  });
}

const laneSidesOf = (state) => MODES.map((mode) => draftPlacementTargets({ mode, form: state.form, adv: state.adv }));

// Every function outside the two stack editors that writes a feature-slot
// input, and why it may: the transition owner, a restore that installs a state
// as is (R11: Undo/Redo, Session load, Import, Reset Settings, which also
// clears placements), a reconcile that keeps the feature slot (R10: managed
// Depth rows, Generate's in-place canonical stack), or a payload migration or
// field list that writes no input. A new writer fails here.
const WRITERS = {
  'app/app-setup.js': { removeCircularDepthTrack: 'reconcile', removeLinearDepthTrack: 'reconcile' },
  'app/feature-editor/placement-actions.js': { changeLayoutSetting: 'transition', restoreTrackLayout: 'transition' },
  'app/run-analysis.js': { runAnalysisInternal: 'reconcile' },
  'services/config.js': { applyConfigData: 'restore',
    reconcileDepthTrackStateAfterSessionFiles: 'restore', restoreStoredNonCanonicalConfig: 'payload' },
  'services/reset.js': { resetSettings: 'restore' }
};
const FIELD = '(?:track_type|linear_track_layout|separate_strands|(?:circular|linear)_track_slots(?:_enabled|_axis_index)?)';
const WRITES = [
  new RegExp(`\\b(?:form|adv)\\.${FIELD}\\s*(?:=(?![=>])|\\+\\+|--)`, 'g'),
  /\badv\.(?:circular|linear)_track_slots(?:\[[^\]]*\])?\.(?:splice|push|pop|shift|unshift|sort|reverse|fill|copyWithin)\(/g,
  /\badv\.(?:circular|linear)_track_slots\[[^\]]+\]\s*=(?!=)/g,
  /\b(?:form|adv)\[[^\]]+\]\s*=(?!=)/g,
  /(?:replaceReactiveObject|replacePlainObject|Object\.assign)\(\s*(?:(?:state|drawing)\.)?(?:form|adv)\b/g,
  // The merge of stored settings into a drawing (`form` and `adv` included).
  /\binstallStoredSettings\(\s*(?:state|drawing)\b/g,
  // The settings slice (services/drawing/settings.js): Reset Settings, and an
  // install that gives `form` or `adv`, replace the whole member.
  /\bSETTINGS_SLICE\.(?:reset\(|install\([^;]*?\b(?:form|adv)\s*:)/g
];
const DECLARATION = /^\s*(?:export\s+)?(?:const|let|function|async\s+function)\s+([A-Za-z_$][\w$]*)/;
// The nearest declaration that encloses a line, by indentation.
const enclosingFunction = (lines, lineIndex) => {
  const own = lines[lineIndex].match(DECLARATION);
  if (own && /=>|function/.test(lines[lineIndex])) return own[1];
  let limit = lines[lineIndex].search(/\S/);
  for (let index = lineIndex - 1; index >= 0; index -= 1) {
    const indent = lines[index].search(/\S/);
    if (indent < 0 || indent >= limit) continue;
    // The end of a multi-line parameter list belongs to the declaration above it.
    const opener = /^\s*[})].*=>\s*\{?\s*$/.test(lines[index])
      ? lines.slice(0, index).findLast((line) => line.search(/\S/) === indent) : lines[index];
    const match = opener?.match(DECLARATION);
    if (match) return match[1];
    limit = indent;
  }
  return '(module)';
};
const sources = (directory) => readdirSync(directory).flatMap((name) => {
  const path = join(directory, name);
  return statSync(path).isDirectory() ? sources(path) : path.endsWith('.js') ? [path] : [];
});

test('every writer of a feature-slot input outside the stack editors is listed (R10, R11)', () => {
  const found = {};
  for (const path of sources(`${WEB}/js`)) {
    const file = relative(`${WEB}/js`, path);
    if (/^app\/(?:circular|linear)-track-slots\.js$/.test(file)) continue;
    const source = readFileSync(path, 'utf8');
    const lines = source.split('\n');
    for (const pattern of WRITES) {
      for (const match of source.matchAll(pattern)) {
        const name = enclosingFunction(lines, source.slice(0, match.index).split('\n').length - 1);
        (found[file] ||= {})[name] = WRITERS[file]?.[name] || 'unlisted: route it through changeTrackLayout';
      }
    }
  }
  assert.deepEqual(found, WRITERS);
});

// A handler is an `@event` attribute, which passes its DOM event as `$event`,
// or the `run` of an ActionMenu item in `:items`, which passes the item's click
// as its first parameter (the track slot Row actions menus): an arrow function
// (with or without parentheses or `async`) or a method. A `run` with no
// parameter, or a bare function reference, has no event to pass (event '').
const ITEM_RUN = /run(?:\s*:\s*(?:(?:async\s*)?(?:\(([^)]*)\)|(\w+))\s*=>\s*([^}]*)|([\w.$]+))|\s*\(([^)]*)\)\s*\{([^}]*)\})/g;
const templateHandlers = (html) => [
  ...[...html.matchAll(/@[\w.:-]+="([^"]*)"/g)].map((match) => ({ handler: match[1], index: match.index, event: '$event' })),
  ...[...html.matchAll(/:items="([^"]*)"/g)].flatMap((match) => [...match[1].matchAll(ITEM_RUN)]
    .map(([, params, param, body, reference, methodParams, methodBody]) => {
      if (reference !== undefined) return { handler: `${reference}()`, index: match.index, event: '' };
      const first = (list) => (list.split(',')[0] || '').trim();
      return methodBody === undefined
        ? { handler: body, index: match.index, event: param ?? first(params) }
        : { handler: methodBody, index: match.index, event: first(methodParams) };
    }))
];

test('the template scan reads every ActionMenu run, with or without a parameter', () => {
  const html = `<action-menu :items="[{ label: 'A', run: (event) => removeCircularTrackSlot(entry.index, event) },
    { label: 'B', run: () => removeCircularTrackSlot(entry.index) }, { label: 'C', run: removeCircularTrackSlot },
    { label: 'D', run: event => removeCircularTrackSlot(entry.index, event) },
    { label: 'E', run: async () => removeCircularTrackSlot(entry.index) },
    { label: 'F', run(event) { removeCircularTrackSlot(entry.index, event) } }]"></action-menu>`;
  assert.deepEqual(templateHandlers(html).map(({ handler, event }) => [handler.trim(), event]), [
    ['removeCircularTrackSlot(entry.index, event)', 'event'],
    ['removeCircularTrackSlot(entry.index)', ''],
    ['removeCircularTrackSlot()', ''],
    ['removeCircularTrackSlot(entry.index, event)', 'event'],
    ['removeCircularTrackSlot(entry.index)', ''],
    ['removeCircularTrackSlot(entry.index, event)', 'event']
  ]);
});

test('templates reach the feature-slot inputs only through the transition (R10, R11)', () => {
  const html = readFileSync(`${WEB}/index.html`, 'utf8');
  const inputs = '(?:form\\.(?:track_type|linear_track_layout|separate_strands)|adv\\.(?:circular|linear)_track_slots\\w*'
    + '|entry\\.slot\\.(?:enabled|renderer|side|params\\.(?:lane_direction|lanes)))';
  for (const [, model] of html.matchAll(/v-model(?:\.\w+)*="([^"]*)"/g)) {
    assert.doesNotMatch(model.trim(), new RegExp(`^${inputs}$`), `v-model="${model}" writes a feature-slot input`);
  }
  const names = new Set(['featurePlacementActions.changeLayoutSetting', ...MODES.flatMap(transitionalExports)]);
  // The whole tag of an attribute, read with quotes so `=>` does not end it.
  const tagAt = (index) => {
    let start = html.lastIndexOf('<', index);
    let quote = false;
    for (let cursor = start; cursor < html.length; cursor += 1) {
      if (html[cursor] === '"') quote = !quote;
      else if (html[cursor] === '>' && !quote) return html.slice(start, cursor + 1);
    }
    return html.slice(start);
  };
  const calls = [];
  for (const { handler, index, event } of templateHandlers(html)) {
    assert.doesNotMatch(handler, new RegExp(`(?:^|[^.\\w])${inputs}\\s*=(?!=)`),
      `${handler} assigns a feature-slot input`);
    for (const name of names) {
      const escaped = name.replace(/[.*+?^${}()|[\]\\]/g, '\\$&');
      for (const [call, args] of handler.matchAll(new RegExp(`(?:^|[^.\\w])${escaped}\\(([^)]*)\\)`, 'g'))) {
        calls.push(name);
        // The control names itself for Cancel and focus; its History adapter
        // records an edit that loses no lane (R11).
        assert.ok(event && args.split(',').map((arg) => arg.trim()).includes(event), `${call.trim()} passes no ${event || 'event'}`);
        assert.doesNotMatch(tagAt(index), /data-history-managed/, `${call.trim()} is data-history-managed`);
      }
    }
  }
  assert.ok(calls.length >= 30, 'the stack and layout controls call the transition');
  assert.match(html, /<button[^>]*data-history-managed[^>]*@click="featurePlacementActions\.resolveLayoutChange\('reset'\)"/);
});

test('the stack editors are created once, with the transition', () => {
  const sites = sources(`${WEB}/js`).flatMap((path) => [...readFileSync(path, 'utf8')
    .matchAll(/create(?:Circular|Linear)TrackSlotEditor\(\{([^}]*)\}\)/g)]
    .map((match) => [relative(`${WEB}/js`, path), match[1].replace(/\s+/g, ' ').trim()]));
  // Both editors take the Legend port their stack edits run in. The Circular
  // editor also takes the choice dialog and History's undoable step, which its
  // Hide GC content / skew question uses, and the reader of the record lengths
  // its estimates describe.
  assert.deepEqual(sites, [
    ['app/app-setup.js', 'state, changeTrackLayout, askChoice, runUndoable: history.runUndoable, estimateRecordLengths: readCircularEstimateLengths, changeTrackVisibility: retireShownLegendStyles'],
    ['app/app-setup.js', 'state, changeTrackLayout, changeTrackVisibility: retireShownLegendStyles']
  ]);
});

test('a stack edit runs in the Legend port once, and a deferred Reset apply runs in it on its own', () => {
  const log = [];
  const visibility = (change) => { log.push('enter'); const result = change(); log.push('leave'); return result; };
  // No lane is lost: the layout transition applies at once.
  const direct = trackLayoutWithVisibility((apply) => apply(), visibility);
  assert.equal(direct(() => { log.push('apply'); return 'done'; }), 'done');
  assert.deepEqual(log, ['enter', 'apply', 'leave']);
  // A lane is lost: the transition applies, restores, keeps the edit for its dialog, and returns false.
  log.length = 0;
  let pending = null;
  const asking = trackLayoutWithVisibility((apply) => { apply(); log.push('restore'); pending = apply; return false; }, visibility);
  assert.equal(asking(() => { log.push('apply'); }), false);
  assert.deepEqual(log, ['enter', 'apply', 'restore', 'leave']);
  // The Reset choice applies the kept edit later, inside the Legend port.
  log.length = 0;
  pending();
  assert.deepEqual(log, ['enter', 'apply', 'leave']);
});
