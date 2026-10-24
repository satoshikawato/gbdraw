// The drawings of a project beyond the two mode drawings (PD-OI-090, B4).
// Reset resets each drawing by its own mode (R17-5): a drawing's reset leaves
// the other drawings and the one-per-Session values alone; the project part
// resets those values, as 0.14.0 Reset Settings did.
import assert from 'node:assert/strict';
import { File } from 'node:buffer';
import test from 'node:test';
import { installFakeSvgDom } from './fake-svg-dom.mjs';
import { modeDrawing } from './helpers/drawing-state.mjs';

globalThis.window = {
  Vue: {
    ref: (value) => ({ value }), shallowRef: (value) => ({ value }), reactive: (value) => value,
    isRef: (value) => Boolean(value && typeof value === 'object' && 'value' in value),
    computed: (getter) => ({ get value() { return getter(); } }), nextTick: async () => {},
    watch: () => () => {}, onMounted: () => {}, onUnmounted: () => {}, provide: () => {}
  },
  DOMPurify: { sanitize: (value) => value }
};
globalThis.document = { createElement: () => ({ click: () => {} }), querySelector: () => null };
installFakeSvgDom();
globalThis.File = File;

const { state, createDrawingState } = await import('../../gbdraw/web/js/state.js');
const { createAppSetup } = await import('../../gbdraw/web/js/app/app-setup.js');
const { bindDrawingInputFileTree } = await import('../../gbdraw/web/js/services/drawing/inputs.js');
const { resetProject, resetSettings } = await import('../../gbdraw/web/js/services/reset.js');
const { createDefaultLosatExecution } = await import('../../gbdraw/web/js/services/session-active-config-contract.js');
const { createDefaultEditorDraftState, createDefaultPaletteDraftState } = await import('../../gbdraw/web/js/state.js');

// The one-per-Session values (R17-5): Reset this drawing keeps them.
const sessionValues = () => JSON.parse(JSON.stringify({
  losatExecution: state.losatExecution,
  richFeaturePopup: state.richFeaturePopup.value,
  autoLabelReflowEnabled: state.autoLabelReflowEnabled.value,
  downloadDpi: state.downloadDpi.value,
  paletteInstantPreviewEnabled: state.paletteInstantPreviewEnabled.value
}));
const editSessionValues = () => {
  state.losatExecution.threadsPerJob = 2;
  state.richFeaturePopup.value = false;
  state.autoLabelReflowEnabled.value = true;
  state.downloadDpi.value = 72;
  state.paletteInstantPreviewEnabled.value = true;
};
const DEFAULT_SESSION_VALUES = {
  losatExecution: createDefaultLosatExecution(),
  richFeaturePopup: true,
  autoLabelReflowEnabled: false,
  downloadDpi: createDefaultEditorDraftState().downloadDpi,
  paletteInstantPreviewEnabled: createDefaultPaletteDraftState().paletteInstantPreviewEnabled
};

test('a drawing\'s reset leaves the other drawing and the Session values; the project reset restores those', () => {
  const circular = modeDrawing(state, 'circular');
  const linear = modeDrawing(state, 'linear');
  circular.form.plot_title = 'Circular title';
  linear.form.plot_title = 'Linear title';
  editSessionValues();
  const edited = sessionValues();
  assert.notDeepEqual(edited, JSON.parse(JSON.stringify(DEFAULT_SESSION_VALUES)));

  const artifact = resetSettings(state, circular);
  assert.equal(circular.form.plot_title, '');
  assert.equal(linear.form.plot_title, 'Linear title', 'the other drawing keeps its settings');
  assert.deepEqual(sessionValues(), edited, 'the one-per-Session values stay');
  assert.equal(artifact.appliedPaletteName, createDefaultPaletteDraftState().appliedPaletteName);

  resetProject(state);
  assert.deepEqual(sessionValues(), JSON.parse(JSON.stringify(DEFAULT_SESSION_VALUES)));
  assert.equal(linear.form.plot_title, 'Linear title', 'the project reset writes no drawing');
  linear.form.plot_title = '';
});

// The follow-ups run for every drawing by its mode, not for the drawing of a
// mode: a third drawing gets its Record Order rows and its comparison rings.
test('Reset all gives a second Circular drawing its Record Order rows and comparison rings', async () => {
  const app = createAppSetup();
  const history = window.__GBDRAW_HISTORY__;
  const second = createDrawingState('circular-2', 'circular');
  state.addDrawing(second, 'Circular 2', 2);
  history.addMember('circular-2');
  second.circularRecordList.value = [{ selector: '#1', record_id: 'NC_1' }, { selector: '#2', record_id: 'NC_2' }];
  second.adv.multi_record_positions.splice(0, Infinity, { selector: '#2', row: 3 });
  bindDrawingInputFileTree(state.inputFiles, second.slices.inputs, 'c_conservation_fastas', [new File(['>r\nACGT\n'], 'ring.fa')]);

  await app.resetSettings();
  assert.deepEqual(second.adv.multi_record_positions.map((position) => position.selector), ['#1', '#2'],
    'one Record Order row per loaded record, at its default');
  assert.equal(second.circularConservation.enabled, true, 'the kept ring file gives a ring');
  assert.equal(second.circularConservation.series.length, 1);
  assert.equal(second.files.c_conservation_fastas.length, 1, 'the ring file stays');
});
