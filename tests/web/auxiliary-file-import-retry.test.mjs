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

// Any other port of the watcher owner: a function that does nothing.
const port = () => new Proxy({}, { get: () => () => {} });
let inputs = 0;
const owner = { calls: 0, outcome: /** @type {number | false} */ (false), removal: async () => true };
/** @type {{ source: Function, callback: Function }[]} */
const watchers = [];
const imports = createAuxiliaryFileImports((target) => target);
const { canRetryAuxiliaryImportFailure, retryAuxiliaryImportFailure } = setupWatchers({
  state,
  rulePreparation: /** @type {any} */ ({
    snapshot: () => ({ inputs, linearFiles: [], inputFiles: [] }),
    isCurrent: (before) => before.inputs === inputs
  }),
  trackAuxiliaryFileImport: imports.track,
  ref: (value) => ({ value }),
  computed: (getter) => ({ get value() { return getter(); } }),
  watch: (source, callback) => {
    watchers.push({ source: typeof source === 'function' ? source : () => source, callback });
    return () => {};
  },
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
    // The other drawing holds no rows of an earlier rules file (OV-1280).
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
const specificRulesWatcher = watchers.find(({ source }) => {
  state.files.t_color = 'probe';
  try { return source() === 'probe'; } finally { state.files.t_color = null; }
});

test('a rules import the owner gave up on keeps Retry through later input changes, and Retry applies it', async () => {
  const file = { name: 'rules.tsv', text: async () => 'CDS\tproduct\tNADH\t#123456\tImported\n' };
  state.errorLog.value = null;
  state.files.t_color = file;
  const result = await specificRulesWatcher.callback(file, null);

  assert.deepEqual({ result, calls: owner.calls, pending: imports.isPending() }, { result: false, calls: 1, pending: false });
  assert.equal(state.errorLog.value?.code, 'RULE_INPUTS_CHANGED');
  assert.equal(state.errorLog.value?.operation, 'evaluateRules');
  assert.equal(canRetryAuxiliaryImportFailure.value, true);
  // A rerender or another upload moves the inputs again: Retry runs on the new ones.
  inputs += 1;
  assert.equal(canRetryAuxiliaryImportFailure.value, true);

  owner.outcome = 1;
  assert.equal(await retryAuxiliaryImportFailure(), true);
  assert.deepEqual({ calls: owner.calls, selected: state.files.t_color, error: state.errorLog.value },
    { calls: 2, selected: file, error: null });
});

// REVIEW-1 N1: Remove on the card while the rule owner refuses to take the
// rows back (its matching failed, or a Legend row conflicts) leaves the rows,
// so the card comes back and the error is shown; nothing rejects unhandled.
test('a Specific rules removal the rule owner refuses brings the card back and reports why', async () => {
  const file = { name: 'rules.tsv', text: async () => '' };
  owner.removal = async () => { throw { code: 'LEGEND_ROW_CONFLICT' }; };
  state.errorLog.value = null;
  state.files.t_color = null;
  assert.equal(await specificRulesWatcher.callback(null, file), false);
  assert.deepEqual({ selected: state.files.t_color, operation: state.errorLog.value?.operation, pending: imports.isPending() },
    { selected: file, operation: 'evaluateRules', pending: false });
  owner.removal = async () => true;
});

// REVIEW-2 N2: the shown drawing's removal runs first; one it refuses leaves
// the drawing that is not shown unchanged, so the card comes back with no
// half-done removal.
test('a Specific rules removal the shown drawing refuses never changes the drawing not shown', async () => {
  const file = { name: 'rules.tsv', text: async () => '' };
  const shown = state.activeDrawing();
  /** @type {unknown[]} */
  const removed = [];
  owner.removal = async (drawing) => {
    removed.push(drawing);
    return drawing !== shown;
  };
  state.files.t_color = null;
  assert.equal(await specificRulesWatcher.callback(null, file), false);
  assert.deepEqual({ removed, selected: state.files.t_color }, { removed: [shown], selected: file });
  owner.removal = async () => true;
});
