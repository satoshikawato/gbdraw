// OV-1206: the Specific rules owner (`runSpecificRuleWrite`) is the one place
// that retries a rules import whose inputs changed. The file watcher does not
// retry it again, and an import the owner gave up on is reported with Retry,
// never dropped or unselected without a message.
import assert from 'node:assert/strict';
import { test } from 'node:test';

globalThis.window = {
  Vue: {
    ref: (value) => ({ value }), reactive: (value) => value,
    computed: (getter) => ({ get value() { return getter(); } }),
    nextTick: async () => {}
  },
  DOMPurify: { sanitize: (value) => value }
};
globalThis.document = {};

const { state } = await import('../../gbdraw/web/js/state.js');
const { createAuxiliaryFileImports, setupWatchers } = await import('../../gbdraw/web/js/app/watchers.js');
const { bindInputFile } = await import('../../gbdraw/web/js/services/drawing/inputs.js');
const { installShownDrawing } = await import('../../gbdraw/web/js/services/config.js');
// The shown drawing's card (R16-4): its binding and its one writer.
const card = (key) => state.activeDrawing().files[key];
const bindCard = (key, file) => bindInputFile(state.inputFiles, state.activeDrawing().slices.inputs, key, file);

// Any other port of the watcher owner: a function that does nothing.
const port = () => new Proxy({}, { get: () => () => {} });
let inputs = 0;
const owner = { calls: 0, outcome: /** @type {number | false} */ (false), removal: async () => true };
const imports = createAuxiliaryFileImports((target) => target, state.inputFiles);
const {
  importAuxiliaryTable, canRetryAuxiliaryImportFailure, retryAuxiliaryImportFailure, auxiliaryImportFailureOnArrival
} = setupWatchers({
  state,
  rulePreparation: /** @type {any} */ ({
    snapshot: () => ({ inputs, linearFiles: [], inputFiles: [] }),
    isCurrent: (before) => before.inputs === inputs
  }),
  trackAuxiliaryFileImport: imports.track,
  ref: (value) => ({ value }),
  computed: (getter) => ({ get value() { return getter(); } }),
  watch: () => () => {},
  nextTick: async () => {},
  onMounted: () => {},
  legendActions: port(),
  // The owner retries while another auxiliary file changes the rule inputs
  // during each of its attempts, then gives up (false, no error).
  featureActions: new Proxy({
    importSpecificRules: async () => {
      owner.calls += 1;
      if (owner.outcome === false) inputs += 3;
      return owner.outcome;
    },
    removeImportedSpecificRules: async (drawing) => owner.removal(drawing)
  }, { get: (target, key) => target[key] ?? (() => {}) }),
  legendLayout: port(),
  resultsManager: port(),
  runLabelReflow: () => {},
  refreshCircularRecordOrder: () => {},
  refreshLinearRecordSelectors: () => {},
  resetPreviewViewport: () => {},
  resetRightDrawer: () => {},
  closeLabelTextScopeDialog: () => {},
  closeHiddenLabelTextDialog: () => {},
  clearLabelBuildNotices: () => {},
  previewRuntime: port()
});
test('a rules import the owner gave up on keeps Retry through later input changes, and Retry applies it', async () => {
  const file = { name: 'rules.tsv', size: 1, text: async () => 'CDS\tproduct\tNADH\t#123456\tImported\n' };
  // The card holds the file's ref in the input file store.
  const ref = state.inputFiles.registerFile(file);
  state.errorLog.value = null;
  const result = await importAuxiliaryTable('t_color', file);

  assert.deepEqual({ result, calls: owner.calls, pending: imports.isPending() }, { result: false, calls: 1, pending: false });
  assert.equal(state.errorLog.value?.code, 'RULE_INPUTS_CHANGED');
  assert.equal(state.errorLog.value?.operation, 'evaluateRules');
  assert.equal(canRetryAuxiliaryImportFailure.value, true);
  // A rerender or another upload moves the inputs again: Retry runs on the new ones.
  inputs += 1;
  assert.equal(canRetryAuxiliaryImportFailure.value, true);

  owner.outcome = 1;
  assert.equal(await retryAuxiliaryImportFailure(), true);
  assert.deepEqual({ calls: owner.calls, selected: card('t_color'), error: state.errorLog.value },
    { calls: 2, selected: ref, error: null });
});

// REVIEW-1 N1: Remove on the card while the rule owner refuses to take the
// rows back (its matching failed, or a Legend row conflicts) leaves the rows,
// so the card comes back and the error is shown; nothing rejects unhandled.
test('a Specific rules removal the rule owner refuses brings the card back and reports why', async () => {
  const file = state.inputFiles.registerFile({ name: 'rules.tsv', size: 1, text: async () => '' });
  owner.removal = async () => { throw { code: 'LEGEND_ROW_CONFLICT' }; };
  state.errorLog.value = null;
  bindCard('t_color', file);
  assert.equal(await importAuxiliaryTable('t_color', null), false);
  assert.deepEqual({ selected: card('t_color'), operation: state.errorLog.value?.operation, pending: imports.isPending() },
    { selected: file, operation: 'evaluateRules', pending: false });
  owner.removal = async () => true;
});

// R16-4: a card changes only the shown drawing. Remove takes back the rows of
// the shown drawing alone, and the other drawing keeps its own card file.
test('Remove on a card takes back only the shown drawing\'s rows and binding', async () => {
  const file = state.inputFiles.registerFile({ name: 'rules.tsv', size: 1, text: async () => '' });
  const other = state.drawingList().find((drawing) => drawing !== state.activeDrawing());
  const shown = state.activeDrawing();
  bindInputFile(state.inputFiles, other.slices.inputs, 't_color', file);
  /** @type {unknown[]} */
  const removed = [];
  owner.removal = async (drawing) => {
    removed.push(drawing);
    return true;
  };
  bindCard('t_color', file);
  assert.equal(await importAuxiliaryTable('t_color', null), true);
  assert.deepEqual({ removed, selected: card('t_color'), other: other.files.t_color }, { removed: [shown], selected: null, other: file });
  bindInputFile(state.inputFiles, other.slices.inputs, 't_color', null);
  owner.removal = async () => true;
});

// The drawing switch as the composition root runs it (`showDrawingNow`): the
// arrival hook reads the departing drawing's message before the switch clears
// it (OV-372), and the arriving drawing's card failure is shown again.
const showDrawing = (drawing) => {
  const departing = state.activeDrawing();
  const arrivingCardFailure = auxiliaryImportFailureOnArrival(departing, drawing);
  installShownDrawing(drawing.id);
  state.errorLog.value = null;
  if (arrivingCardFailure !== null) state.errorLog.value = arrivingCardFailure;
};

// R16-4, D-15-126: a failed card import is the card state of its drawing. Its
// message and Retry show while that drawing is shown, on no other drawing, and
// come back with that drawing; Retry is offered only there.
test('a failed card import\'s message and Retry belong to the drawing whose card failed', async () => {
  const file = { name: 'rules.tsv', size: 1, text: async () => 'CDS\tproduct\tNADH\t#123456\tImported\n' };
  const failedIn = state.activeDrawing();
  bindCard('t_color', null);
  state.errorLog.value = null;
  owner.outcome = false;
  assert.equal(await importAuxiliaryTable('t_color', file), false);
  const failure = state.errorLog.value;
  assert.equal(failure?.code, 'RULE_INPUTS_CHANGED');
  assert.equal(canRetryAuxiliaryImportFailure.value, true);
  const other = state.drawingList().find((drawing) => drawing !== failedIn);
  showDrawing(other);
  assert.equal(state.errorLog.value, null, 'no message on another drawing');
  assert.equal(canRetryAuxiliaryImportFailure.value, false, 'another drawing shown: no Retry');
  assert.equal(await retryAuxiliaryImportFailure(), false);
  assert.equal(other.files.t_color, null, 'nothing was bound in the other drawing');
  showDrawing(failedIn);
  assert.equal(state.errorLog.value, failure, 'the message comes back with its drawing');
  assert.equal(canRetryAuxiliaryImportFailure.value, true);
  owner.outcome = 1;
  assert.equal(await retryAuxiliaryImportFailure(), true);
  assert.equal(state.inputFiles.fileOf(failedIn.files.t_color), file);
  assert.equal(state.errorLog.value, null);
  bindCard('t_color', null);
});

// A dismissed message ends the failure: it does not come back with its drawing.
test('a dismissed card failure does not come back when its drawing is shown again', async () => {
  const file = { name: 'rules.tsv', size: 1, text: async () => 'CDS\tproduct\tNADH\t#123456\tImported\n' };
  const failedIn = state.activeDrawing();
  bindCard('t_color', null);
  state.errorLog.value = null;
  owner.outcome = false;
  assert.equal(await importAuxiliaryTable('t_color', file), false);
  state.errorLog.value = null;
  const other = state.drawingList().find((drawing) => drawing !== failedIn);
  showDrawing(other);
  showDrawing(failedIn);
  assert.deepEqual({ error: state.errorLog.value, retry: canRetryAuxiliaryImportFailure.value }, { error: null, retry: false });
  owner.outcome = 1;
});

// B2b-1a-iii: a table card names its file by the ref of the input file store
// (`state.inputFiles`); the import reads the bytes through the store.
test('a table card holds the ref of the file it imported', async () => {
  const file = { name: 'priority.tsv', size: 1, text: async () => '' };
  state.errorLog.value = null;
  await importAuxiliaryTable('qualifier_priority', file);
  assert.equal(card('qualifier_priority'), state.inputFiles.registerFile(file));
  assert.equal(state.inputFiles.fileOf(card('qualifier_priority')), file);
  assert.equal(imports.wait(file), imports.wait(card('qualifier_priority')));
  bindCard('qualifier_priority', null);
});

// OV-1674: a failed import leaves the card unbound, so the store may release
// the file (History's budget keeps only named and bound files); Retry still
// imports the file the user chose.
test('Retry imports the chosen rules file after the store released it', async () => {
  const file = { name: 'rules.tsv', size: 1, text: async () => 'CDS\tproduct\tNADH\t#123456\tImported\n' };
  bindCard('t_color', null);
  state.errorLog.value = null;
  owner.outcome = false;
  assert.equal(await importAuxiliaryTable('t_color', file), false);
  assert.equal(card('t_color'), null);
  state.inputFiles.retainOnly(new Set());
  owner.outcome = 1;
  assert.equal(await retryAuxiliaryImportFailure(), true);
  assert.deepEqual({ file: state.inputFiles.fileOf(card('t_color')), error: state.errorLog.value },
    { file, error: null });
  bindCard('t_color', null);
});
