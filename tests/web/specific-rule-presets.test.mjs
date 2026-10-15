import assert from 'node:assert/strict';
import { test } from 'node:test';
import { readFileSync } from 'node:fs';
import { dirname, join } from 'node:path';
import { fileURLToPath } from 'node:url';
import { parseSpecificRules } from '../../gbdraw/web/js/services/file-imports.js';
import { createFeatureRuleActions } from '../../gbdraw/web/js/app/feature-editor/rule-actions.js';
import { createRulePreparation } from '../../gbdraw/web/js/app/rule-matching.js';
import { diffLegendIntents } from '../../gbdraw/web/js/services/specific-color-rules.js';
import { evaluatePythonRules } from './helpers/python-rule-evaluator.mjs';
import { withDrawings } from './helpers/drawing-state.mjs';

// The preset schemes of the Colors panel: the tables shipped under `presets/`,
// the softened pair of each tool's table, and the settings a preset writes
// after its rules (`effects`), which the control says.
globalThis.window = { Vue: {
  ref: (value) => ({ value }), reactive: (value) => value,
  computed: (getter) => ({ get value() { return getter(); } }), nextTick: async () => {}
} };
const { state: appState } = await import('../../gbdraw/web/js/state.js');
const presets = appState.specificRulePresets;
const WEB = join(dirname(fileURLToPath(import.meta.url)), '..', '..', 'gbdraw', 'web');
const tableOf = (preset) => parseSpecificRules(readFileSync(join(WEB, preset.path), 'utf8')).rules;
const rowOf = ({ feat, qual, val, cap }) => [feat, qual, val, cap];

test('every preset scheme ships its table under presets/, and each table parses', () => {
  assert.deepEqual(presets.map(preset => preset.id), ['pharokka', 'pharokka_soft', 'bakta', 'bakta_soft']);
  for (const preset of presets) assert.ok(tableOf(preset).length > 0, preset.id);
});

test('a softened preset keeps the rows, captions and effects of its original in other colors', () => {
  for (const preset of presets.filter(entry => entry.id.endsWith('_soft'))) {
    const original = presets.find(entry => entry.id === preset.id.replace(/_soft$/, ''));
    const [before, after] = [tableOf(original), tableOf(preset)];
    assert.deepEqual(after.map(rowOf), before.map(rowOf), preset.id);
    assert.ok(after.every((rule, index) => rule.color !== before[index].color), `${preset.id} recolors every row`);
    assert.equal(preset.effects, original.effects, `${preset.id} writes what ${original.id} writes`);
    assert.ok(preset.label.startsWith(`${original.label.split(' ')[0]}, softened colors`), `${preset.id} is listed as a variant of ${original.id}`);
  }
});

// The apply path with the real rule owner: the rules of the fetched table
// replace the specific rules, then the preset's effects are written, as one
// transaction. The Bakta tables write the CDS default color and the Legend
// sizes; the Pharokka tables write no setting.
const setup = () => {
  const drawing = {
    manualSpecificRules: [], manualPriorityRules: [],
    extractedFeatures: { value: [{ type: 'CDS', svg_id: 'a', qualifiers: { note: ['COG:J'] } }] },
    featureColorOverrides: {}, featureOverrides: {}, results: { value: [{ name: 'figure', content: 'before' }] },
    svgResultIdentity: { value: 'before' }, fileLegendCaptions: { value: new Set() }, addedLegendCaptions: { value: new Set() },
    legendEntries: { value: [] }, files: { t_color: null }, legendColorOverrides: {}, labelReflowProcessing: { value: false },
    editableLabels: { value: [] }, newSpecRule: { feat: '', qual: '', val: '', color: '', cap: '' },
    currentColors: { value: { CDS: '#112233' } }, adv: { legend_box_size: null, legend_font_size: null },
    specificRulePresets: presets, selectedSpecificPreset: { value: '' }, specificRulePresetLoading: { value: false }
  };
  const state = withDrawings(drawing);
  const transactions = [];
  const transact = async (label, commit) => { await commit(); transactions.push(label); };
  const actions = createFeatureRuleActions({
    ref: value => ({ value }), computed: get => ({ get value() { return get(); } }), state,
    rulePreparation: createRulePreparation({ state, evaluate: evaluatePythonRules }),
    runUndoableCheckpoint: transact, runUndoable: transact,
    prepareFileLegendEntries: async (intents) => ({
      diff: diffLegendIntents(drawing.legendEntries.value, intents), isCurrent: () => true,
      apply: () => { drawing.legendEntries.value = intents; }
    }),
    projectPaletteAndRules: () => true, ports: { requestAutomaticRerender: () => true }, watch: () => () => {}
  });
  return { drawing, actions, transactions };
};
const servePresets = () => {
  const previous = globalThis.fetch;
  globalThis.fetch = async (path) => new Response(readFileSync(join(WEB, String(path))));
  return () => { globalThis.fetch = previous; };
};

for (const [id, rows, colors, adv] of [
  ['bakta', 26, { CDS: '#cccccc' }, { legend_box_size: 12, legend_font_size: 12 }],
  ['bakta_soft', 26, { CDS: '#cccccc' }, { legend_box_size: 12, legend_font_size: 12 }],
  ['pharokka', 16, { CDS: '#112233' }, { legend_box_size: null, legend_font_size: null }],
  ['pharokka_soft', 16, { CDS: '#112233' }, { legend_box_size: null, legend_font_size: null }]
]) {
  test(`applying the ${id} preset writes its ${rows} rules and its effects in one step`, async () => {
    const restore = servePresets();
    try {
      const { drawing, actions, transactions } = setup();
      drawing.selectedSpecificPreset.value = id;
      assert.equal(await actions.applySpecificRulePreset(), true);
      assert.equal(drawing.manualSpecificRules.length, rows);
      assert.deepEqual(drawing.manualSpecificRules.map(rowOf), tableOf(presets.find(preset => preset.id === id)).map(rowOf));
      assert.deepEqual(drawing.currentColors.value, colors);
      assert.deepEqual(drawing.adv, adv);
      assert.deepEqual(transactions, ['Apply specific color preset']);
    } finally {
      restore();
    }
  });
}
