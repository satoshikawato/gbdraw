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

const WEB = 'gbdraw/web';
const MODES = ['circular', 'linear'];
const lane = (scope, side) => ({ scope, recordKey: `${scope}-record`, biologicalFeatureId: `${scope}-cds`,
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
  const steps = [];
  const placement = createFeaturePlacementActions({ state, getCommittedRequest: () => null, isCurrentFeature: () => true,
    history: { runUndoable: async (label, fn) => { steps.push(label); fn(); } } });
  const editor = createCircularTrackSlotEditor({ state, trackLayoutActions: placement.trackLayoutActions });
  editor.normalizeCircularTrackSlots();
  createLinearTrackSlotEditor({ state }).normalizeLinearTrackSlots();
  Object.assign(state.featurePlacementOverrides, { c: lane('circular', 'outward'), l: lane('linear', 'above') });
  const features = () => state.adv.circular_track_slots.findIndex((slot) => slot.renderer === 'features');
  const before = everything(state);
  const click = () => {
    const event = new Event('click');
    Object.defineProperty(event, 'currentTarget', { value: { tagName: 'BUTTON', textContent: '', title: 'Move outside Axis' } });
    return event;
  };
  assert.equal(editor.moveCircularTrackSlotOutside(features(), click()), false);
  assert.deepEqual({ ...placement.layoutChange },
    { open: true, count: 1, setting: 'Move outside Axis', value: '', scope: '' });
  assert.equal(everything(state), before);
  await placement.resolveLayoutChange('cancel');
  assert.equal(placement.layoutChange.open, false);
  assert.equal(everything(state), before);
  assert.deepEqual(steps, []);

  editor.moveCircularTrackSlotOutside(features(), click());
  await placement.resolveLayoutChange('reset');
  assert.equal(state.adv.circular_track_slots[features()].params.lane_direction, 'outside');
  // Only the lost row goes; the other mode's lane keeps its draft row (R2).
  assert.deepEqual(Object.keys(state.featurePlacementOverrides), ['l']);
  assert.deepEqual(steps, ['Change setting and reset Feature placements']);

  // The Circular panel's Separate Strands is the Linear predicate's input.
  const strands = { type: 'checkbox', checked: true, getAttribute: () => 'Separate Strands' };
  assert.equal(placement.changeLayoutSetting({ target: strands }, 'separate_strands'), false);
  assert.deepEqual({ ...placement.layoutChange },
    { open: true, count: 1, setting: 'Separate Strands', value: 'On', scope: 'Linear ' });
  assert.equal(strands.checked, false);
  assert.equal(state.form.separate_strands, false);
  await placement.resolveLayoutChange('cancel');
  assert.equal(state.form.separate_strands, false);
  assert.deepEqual(steps, ['Change setting and reset Feature placements']);

  // An edit that loses no lane applies now; the control's adapter records it.
  state.featurePlacementOverrides.l.placement = { kind: 'main' };
  assert.equal(placement.changeLayoutSetting({ target: { ...strands, checked: true } }, 'separate_strands'), undefined);
  assert.equal(state.form.separate_strands, true);
  assert.equal(placement.layoutChange.open, false);
  assert.equal(steps.length, 1);
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
const STACK = {
  circular: { normalize: 'normalizeCircularTrackSlots', offAxis: 'moveCircularTrackSlotOutside' },
  linear: { normalize: 'normalizeLinearTrackSlots', offAxis: 'moveLinearTrackSlotAbove' }
};

for (const mode of MODES) {
  test(`${mode} stack editor writes the feature-slot inputs only through the transition (R10)`, () => {
    const fixture = (offAxis) => {
      const state = draftState(mode);
      state.adv.linear_track_slots.push(
        { id: 'gc_content', renderer: 'dinucleotide_content', enabled: true, side: 'below', params: { nt: 'GC' } },
        { id: 'gap', renderer: 'spacer', enabled: true, side: 'below', height: '12px', params: {} });
      let raw = null;
      const windows = [];
      const editor = EDITORS[mode]({ state, trackLayoutActions: (actions) => {
        raw = actions;
        return Object.fromEntries(Object.entries(actions).map(([name, action]) => [name, (...args) => {
          const entry = project(state);
          try { return action(...args); } finally { windows.push([entry, project(state)]); }
        }]));
      } });
      editor[STACK[mode].normalize]();
      if (offAxis) raw[STACK[mode].offAxis](state.adv[`${mode}_track_slots`].indexOf(slotOf(state, mode)));
      // The simple controls now differ from the stack, so Reset changes it.
      Object.assign(state.form, { track_type: 'tuckin', linear_track_layout: 'below' });
      windows.length = 0;
      return { state, editor, raw, windows };
    };
    const transitional = Object.keys(fixture(false).raw);
    assert.ok(transitional.length >= 10, 'the stack edits run through trackLayoutActions');
    const changedBy = new Set();
    for (const offAxis of [false, true]) {
      const { editor: probe } = fixture(offAxis);
      for (const name of Object.keys(probe).filter((key) => typeof probe[key] === 'function')) {
        for (let call = 0; call < CALLS(fixture(offAxis).state, mode).length; call += 1) {
          const { state, editor, raw, windows } = fixture(offAxis);
          const args = CALLS(state, mode)[call];
          const lanes = laneSidesOf(state);
          const before = project(state);
          const whole = everything(state);
          const saved = saveTrackLayout(state);
          const label = `${mode} ${name} call ${call}${offAxis ? ' (feature row off the Axis)' : ''}`;
          try { (transitional.includes(name) ? raw : editor)[name](...args); } catch { /* an unfit call */ }
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
  'mode-profiles.js': { writeManagedState: 'other fields' },
  'services/config.js': { applyConfigData: 'restore', migratePersistedWebOptionValues: 'payload',
    reconcileDepthTrackStateAfterSessionFiles: 'restore', restoreStoredNonCanonicalConfig: 'payload' },
  'services/gallery-session-migration.js': { migratePersistedGalleryConfig: 'payload' },
  'services/reset.js': { resetSettings: 'restore' }
};
const FIELD = '(?:track_type|linear_track_layout|separate_strands|(?:circular|linear)_track_slots(?:_enabled|_axis_index)?)';
const WRITES = [
  new RegExp(`\\b(?:form|adv)\\.${FIELD}\\s*(?:=(?![=>])|\\+\\+|--)`, 'g'),
  /\badv\.(?:circular|linear)_track_slots(?:\[[^\]]*\])?\.(?:splice|push|pop|shift|unshift|sort|reverse|fill|copyWithin)\(/g,
  /\badv\.(?:circular|linear)_track_slots\[[^\]]+\]\s*=(?!=)/g,
  /\b(?:form|adv)\[[^\]]+\]\s*=(?!=)/g,
  /(?:replaceReactiveObject|replacePlainObject|safeDeepMerge|Object\.assign)\(\s*(?:state\.)?(?:form|adv)\b/g
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
    if (file.includes('generated') || /^app\/(?:circular|linear)-track-slots\.js$/.test(file)) continue;
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

test('templates reach the feature-slot inputs only through the transition (R10, R11)', () => {
  const html = readFileSync(`${WEB}/index.html`, 'utf8');
  const inputs = '(?:form\\.(?:track_type|linear_track_layout|separate_strands)|adv\\.(?:circular|linear)_track_slots\\w*'
    + '|entry\\.slot\\.(?:enabled|renderer|side|params\\.(?:lane_direction|lanes)))';
  for (const [, model] of html.matchAll(/v-model(?:\.\w+)*="([^"]*)"/g)) {
    assert.doesNotMatch(model.trim(), new RegExp(`^${inputs}$`), `v-model="${model}" writes a feature-slot input`);
  }
  const names = new Set(['featurePlacementActions.changeLayoutSetting']);
  for (const mode of MODES) {
    EDITORS[mode]({ state: draftState(mode), trackLayoutActions: (actions) => {
      Object.keys(actions).forEach((name) => names.add(name));
      return actions;
    } });
  }
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
  for (const match of html.matchAll(/@[\w.:-]+="([^"]*)"/g)) {
    const handler = match[1];
    assert.doesNotMatch(handler, new RegExp(`(?:^|[^.\\w])${inputs}\\s*=(?!=)`),
      `${handler} assigns a feature-slot input`);
    for (const name of names) {
      const escaped = name.replace(/[.$]/g, '\\$&');
      for (const [call, args] of handler.matchAll(new RegExp(`(?:^|[^.\\w])${escaped}\\(([^)]*)\\)`, 'g'))) {
        calls.push(name);
        // The control names itself for Cancel and focus; its History adapter
        // records an edit that loses no lane (R11).
        assert.match(args, /\$event/, `${call.trim()} passes no $event`);
        assert.doesNotMatch(tagAt(match.index), /data-history-managed/, `${call.trim()} is data-history-managed`);
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
  assert.deepEqual(sites, [['app/app-setup.js', 'state, trackLayoutActions'], ['app/app-setup.js', 'state, trackLayoutActions']]);
});
