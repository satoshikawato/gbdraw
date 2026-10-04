import assert from 'node:assert/strict';
import { mkdir, mkdtemp, readFile, writeFile } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { pathToFileURL } from 'node:url';
import { test } from 'node:test';

const repoRoot = process.cwd();
const sourcePath = join(repoRoot, 'gbdraw', 'web', 'js', 'app', 'history-inputs.js');
const indexPath = join(repoRoot, 'gbdraw', 'web', 'index.html');
const tempDir = await mkdtemp(join(tmpdir(), 'gbdraw-history-inputs-'));
await writeFile(join(tempDir, 'package.json'), '{"type":"module"}\n', 'utf8');
await mkdir(join(tempDir, 'app'), { recursive: true });
await writeFile(join(tempDir, 'app', 'history-inputs.js'), await readFile(sourcePath, 'utf8'), 'utf8');

const { isIgnoredTarget, setupHistoryInputs } = await import(pathToFileURL(join(tempDir, 'app', 'history-inputs.js')));
await mkdir(join(tempDir, 'services'), { recursive: true });
for (const name of ['history.js', 'history-files.js', 'runtime-test-hooks.js']) {
  await writeFile(
    join(tempDir, 'services', name),
    await readFile(join(repoRoot, 'gbdraw', 'web', 'js', 'services', name), 'utf8')
  );
}
const { createHistoryManager } = await import(pathToFileURL(join(tempDir, 'services', 'history.js')));

const createSelectHarness = async (attributes = {}) => {
  let intent = { labels_mode: 'out' };
  const history = createHistoryManager({
    buildIntent: () => ({ ...intent }),
    applyIntent: (restored) => { intent = { ...restored }; },
    buildCheckpoint: () => assert.fail('Select edits must use intent history'),
    applyCheckpoint: () => assert.fail('Select edits must restore intent history')
  });
  await history.initializeIntentBaseline();
  const listeners = new Map();
  const root = {
    addEventListener: (name, handler) => listeners.set(name, handler),
    removeEventListener: (name) => listeners.delete(name)
  };
  const select = {
    tagName: 'SELECT',
    closest: (selector) => {
      if (selector.startsWith('input,')) return select;
      return attributes.managed && selector.includes('data-history-managed') ? select : null;
    },
    ...attributes
  };
  const cleanup = setupHistoryInputs({ root, history, nextTick: () => Promise.resolve() });
  const settle = () => new Promise((resolve) => setImmediate(resolve));
  return {
    history, cleanup, settle,
    dispatch: (name) => listeners.get(name)?.({ target: select }),
    change: (value) => { intent.labels_mode = value; },
    value: () => intent.labels_mode
  };
};

for (const inputMethod of ['keyboard', 'pointer']) {
  test(`${inputMethod} select edits capture pre-change intent once and support Undo/Redo`, async () => {
    const h = await createSelectHarness();
    try {
      if (inputMethod === 'pointer') h.dispatch('pointerdown');
      h.dispatch('focusin');
      await h.settle();
      for (const value of ['both', 'none']) {
        h.dispatch('keydown');
        h.change(value);
        await h.settle();
        h.dispatch('change');
        await h.settle();
      }
      h.dispatch('focusout');
      await h.settle();
      assert.equal(h.history.getUndoCount(), 2);
      assert.equal(h.history.undoLabel(), 'Change setting');
      await h.history.undo();
      assert.equal(h.value(), 'both');
      await h.history.undo();
      assert.equal(h.value(), 'out');
      assert.equal(h.history.getUndoCount(), 0);
      await h.history.redo();
      assert.equal(h.value(), 'both');
      await h.history.redo();
      assert.equal(h.value(), 'none');
      assert.equal(h.history.getUndoCount(), 2);
    } finally {
      h.cleanup();
    }
  });
}

test('leaving an unchanged select closes its transaction without an Undo step', async () => {
  const h = await createSelectHarness();
  try {
    h.dispatch('focusin');
    await h.settle();
    h.dispatch('focusout');
    await h.settle();
    assert.equal(h.history.getUndoCount(), 0);
    // A later external change must not be swept into the abandoned focus transaction.
    h.change('both');
    const next = await h.history.begin('Later edit');
    assert.equal(next.before.labels_mode, 'both');
    h.history.cancel(next);
  } finally {
    h.cleanup();
  }
});

test('disabled and explicitly managed selects stay outside input-adapter history', async () => {
  for (const attributes of [
    { disabled: true },
    { managed: true }
  ]) {
    const h = await createSelectHarness(attributes);
    try {
      h.dispatch('focusin');
      h.dispatch('keydown');
      await h.settle();
      h.change('both');
      h.dispatch('change');
      h.dispatch('focusout');
      await h.settle();
      assert.equal(h.history.getUndoCount(), 0);
    } finally {
      h.cleanup();
    }
  }
});

const targetMatching = (attribute) => ({
  closest: (selector) => (String(selector).includes(attribute) ? {} : null)
});

assert.equal(isIgnoredTarget(null), false);
assert.equal(isIgnoredTarget({ closest: () => null }), false);
assert.equal(isIgnoredTarget(targetMatching('data-history-ignore')), true);
assert.equal(isIgnoredTarget(targetMatching('data-history-managed')), true);
assert.equal(isIgnoredTarget(targetMatching('data-history-scope="transient"')), true);

const indexHtml = await readFile(indexPath, 'utf8');
[
  '@click="resetSettings"',
  '@click="runAnalysis"',
  '@click="$refs.sessionInput.click()"',
  '@change="importSession"',
  // B22: the Add Seq step commits after the ring record label read.
  '@change="addCircularConservationComparisonFile"',
  // OV-09: a layout control and its placement dialog own one step (R10).
  `@change="featurePlacementActions.changeLayoutSetting($event, 'separate_strands')"`,
  `@click="featurePlacementActions.resolveLayoutChange('reset')"`
].forEach((handler) => {
  const escapedHandler = handler.replace(/[.*+?^${}()|[\]\\]/g, '\\$&');
  assert.match(
    indexHtml,
    new RegExp(`<(?:button|input)(?=[^>]*${escapedHandler})(?=[^>]*data-history-managed)[^>]*>`),
    `${handler} must bypass the generic input adapter and retain its explicit History boundary`
  );
});

console.log('history input tests passed');

test('adjacent numeric/unit focus transitions retain distinct operations instead of a closing transaction', async () => {
  let intent = { width: { value: '1.', unit: 'factor' } };
  const history = createHistoryManager({
    buildIntent: () => structuredClone(intent),
    applyIntent: restored => { intent = restored; },
    buildCheckpoint: () => assert.fail('Control edits use intent History'),
    applyCheckpoint: () => assert.fail('Control edits restore intent History')
  });
  await history.initializeIntentBaseline();
  const listeners = new Map();
  const root = {
    addEventListener: (name, handler) => listeners.set(name, handler),
    removeEventListener: name => listeners.delete(name)
  };
  const control = tagName => {
    const element = { tagName, type: 'text', closest: selector => selector.startsWith('input,') ? element : null };
    return element;
  };
  const select = control('SELECT');
  const text = control('INPUT');
  const dispatch = (name, target) => listeners.get(name)?.({ target });
  const settle = () => new Promise(resolve => setImmediate(resolve));
  const cleanup = setupHistoryInputs({ root, history, nextTick: () => Promise.resolve() });
  try {
    dispatch('focusin', select);
    await settle();
    intent.width.unit = 'px';
    dispatch('change', select);
    dispatch('focusout', select);
    dispatch('focusin', text);
    await settle();
    intent.width.value = '1e';
    dispatch('change', text);
    dispatch('focusout', text);
    dispatch('focusin', select);
    await settle();
    intent.width.unit = 'factor';
    dispatch('change', select);
    dispatch('focusout', select);
    await settle();
    assert.equal(history.getUndoCount(), 3);
    await history.undo();
    assert.deepEqual(intent.width, { value: '1e', unit: 'px' });
    await history.undo();
    assert.deepEqual(intent.width, { value: '1.', unit: 'px' });
    await history.undo();
    assert.deepEqual(intent.width, { value: '1.', unit: 'factor' });
    await history.redo();
    await history.redo();
    await history.redo();
    assert.deepEqual(intent.width, { value: '1e', unit: 'factor' });
  } finally {
    cleanup();
  }
});

// R11 boundary matrix (SE-02, SE-03, SE-04, N-18). A trusted browser event runs a
// microtask checkpoint between listeners; the harness reproduces that ordering.
await writeFile(
  join(tempDir, 'app', 'history-shortcuts.js'),
  await readFile(join(repoRoot, 'gbdraw', 'web', 'js', 'app', 'history-shortcuts.js'), 'utf8'),
  'utf8'
);
const { setupHistoryShortcuts } = await import(pathToFileURL(join(tempDir, 'app', 'history-shortcuts.js')));
const drainMicrotasks = async () => { for (let index = 0; index < 50; index += 1) await Promise.resolve(); };
const settleTasks = () => new Promise(resolve => setTimeout(resolve, 5));

const createBoundaryHarness = async () => {
  const state = { flag: false, choice: 'gb', mode: 'circular', prefix: '', labels: 'out', rings: 0 };
  const history = createHistoryManager({
    buildIntent: () => ({ ...state }),
    applyIntent: (restored) => { Object.assign(state, restored); },
    buildCheckpoint: () => assert.fail('Control edits use intent History'),
    applyCheckpoint: () => assert.fail('Control edits restore intent History')
  });
  await history.initializeIntentBaseline();
  const listeners = [];
  const root = {
    addEventListener: (type, handler, capture) => listeners.push({ type, handler, capture: Boolean(capture) }),
    removeEventListener: (type, handler, capture) => {
      const index = listeners.findIndex(entry => (
        entry.type === type && entry.handler === handler && entry.capture === Boolean(capture)
      ));
      if (index >= 0) listeners.splice(index, 1);
    }
  };
  const documentListeners = new Map();
  globalThis.document = {
    addEventListener: (type, handler) => documentListeners.set(type, handler),
    removeEventListener: (type) => documentListeners.delete(type)
  };
  const control = (tagName, type = '') => {
    const element = { tagName, type, isContentEditable: false };
    element.closest = (selector) => {
      if (selector.startsWith('input,')) return element;
      if (selector === 'button') return tagName === 'BUTTON' ? element : null;
      return null;
    };
    return element;
  };
  const controls = {
    checkbox: control('INPUT', 'checkbox'),
    radio: control('INPUT', 'radio'),
    select: control('SELECT', 'select-one'),
    text: control('INPUT', 'text'),
    button: control('BUTTON', 'button'),
    file: control('INPUT', 'file'),
    label: { tagName: 'LABEL', closest: () => null }
  };
  const cleanup = setupHistoryInputs({ root, history, nextTick: () => Promise.resolve() });
  setupHistoryShortcuts({ history, onMounted: callback => callback(), onUnmounted: () => {} });
  // Capture listeners, then the target-phase owner (for example v-model), then bubble.
  // A script-dispatched (untrusted) event runs no microtask checkpoint between listeners.
  const dispatch = async (type, target, { atTarget = null, untrusted = false, ...init } = {}) => {
    const event = { type, target, preventDefault: () => {}, ...init };
    const checkpoint = untrusted ? () => {} : drainMicrotasks;
    for (const entry of listeners.filter(item => item.type === type && item.capture)) {
      entry.handler(event);
      await checkpoint();
    }
    atTarget?.();
    await checkpoint();
    for (const entry of listeners.filter(item => item.type === type && !item.capture)) {
      entry.handler(event);
      await checkpoint();
    }
    if (type === 'keydown') documentListeners.get('keydown')?.(event);
    await drainMicrotasks();
  };
  const key = (target, keyName, modifiers = {}) => dispatch('keydown', target, {
    key: keyName, ctrlKey: false, metaKey: false, altKey: false, shiftKey: false, ...modifiers
  });
  const typeText = async (value) => {
    await dispatch('focusin', controls.text);
    await key(controls.text, value.at(-1));
    state.prefix = value;
  };
  return {
    state, history, controls, dispatch, key, typeText,
    counts: () => [history.getUndoCount(), history.getRedoCount()],
    cleanup: () => { cleanup(); delete globalThis.document; }
  };
};

const toggleThroughLabel = async (h) => {
  await h.dispatch('pointerdown', h.controls.label);
  await h.dispatch('click', h.controls.checkbox);
  await h.dispatch('change', h.controls.checkbox, { atTarget: () => { h.state.flag = !h.state.flag; } });
};
const toggleWithSpace = async (h) => {
  await h.dispatch('focusin', h.controls.checkbox);
  await h.key(h.controls.checkbox, ' ');
  await h.dispatch('click', h.controls.checkbox);
  await h.dispatch('change', h.controls.checkbox, { atTarget: () => { h.state.flag = !h.state.flag; } });
};
const moveRadioWithArrow = async (h) => {
  await h.dispatch('focusin', h.controls.radio);
  await h.key(h.controls.radio, 'ArrowRight');
  await h.dispatch('click', h.controls.radio);
  await h.dispatch('change', h.controls.radio, { atTarget: () => { h.state.choice = 'gff'; } });
};
const clickCheckboxFromText = async (h) => {
  await h.dispatch('pointerdown', h.controls.checkbox);
  await h.dispatch('change', h.controls.text);
  await h.dispatch('focusout', h.controls.text);
  await h.dispatch('focusin', h.controls.checkbox);
  await h.dispatch('click', h.controls.checkbox);
  await h.dispatch('change', h.controls.checkbox, { atTarget: () => { h.state.flag = !h.state.flag; } });
};
const clickButtonFromText = async (h) => {
  await h.dispatch('pointerdown', h.controls.button);
  await h.dispatch('change', h.controls.text);
  await h.dispatch('focusout', h.controls.text);
  await h.dispatch('click', h.controls.button, { atTarget: () => { h.state.mode = 'linear'; } });
  await settleTasks();
};

for (const [name, act, changed] of [
  ['checkbox label text click', toggleThroughLabel, { flag: true }],
  ['checkbox Space key', toggleWithSpace, { flag: true }],
  ['radio Arrow key', moveRadioWithArrow, { choice: 'gff' }]
]) {
  test(`SE-02: ${name} records exactly one Undo step`, async () => {
    const h = await createBoundaryHarness();
    try {
      const before = { ...h.state };
      await act(h);
      await settleTasks();
      assert.deepEqual(h.counts(), [1, 0]);
      assert.deepEqual({ ...h.state }, { ...before, ...changed });
      await h.history.undo();
      assert.deepEqual({ ...h.state }, before);
    } finally { h.cleanup(); }
  });
}

for (const [name, act, changed] of [
  ['checkbox click', clickCheckboxFromText, { flag: true }],
  ['mode button click', clickButtonFromText, { mode: 'linear' }]
]) {
  test(`SE-03: ${name} while a typed text field has focus records its own step`, async () => {
    const h = await createBoundaryHarness();
    try {
      const before = { ...h.state };
      await h.typeText('audit');
      await act(h);
      await settleTasks();
      assert.deepEqual(h.counts(), [2, 0]);
      await h.history.undo();
      assert.deepEqual({ ...h.state }, { ...before, prefix: 'audit' });
      await h.history.undo();
      assert.deepEqual({ ...h.state }, before);
      await h.history.redo();
      await h.history.redo();
      assert.deepEqual({ ...h.state }, { ...before, prefix: 'audit', ...changed });
    } finally { h.cleanup(); }
  });
}

// B21: a hidden file input, opened by a button or a label, mutates state in its
// own change handler. The step begins in the change capture phase.
const addRing = (h, init = {}) => h.dispatch('change', h.controls.file, {
  ...init, atTarget: () => { h.state.rings += 1; }
});
for (const [name, act] of [
  ['file chosen after its opening button', async (h) => {
    await h.dispatch('pointerdown', h.controls.button);
    await h.dispatch('click', h.controls.button);
    await settleTasks();
    await addRing(h);
  }],
  ['file set directly on the hidden input', (h) => addRing(h)],
  ['file set by a script-dispatched change', (h) => addRing(h, { untrusted: true })]
]) {
  test(`B21: ${name} records exactly one Undo step`, async () => {
    const h = await createBoundaryHarness();
    try {
      await act(h);
      await settleTasks();
      assert.deepEqual(h.counts(), [1, 0]);
      assert.equal(h.history.undoLabel(), 'Change uploaded file');
      await h.history.undo();
      assert.equal(h.state.rings, 0);
      await h.history.redo();
      assert.equal(h.state.rings, 1);
    } finally { h.cleanup(); }
  });
}

test('SE-04: shortcuts on a focused select undo and redo without recording the Undo as an edit', async () => {
  const h = await createBoundaryHarness();
  try {
    await h.dispatch('focusin', h.controls.select);
    await h.key(h.controls.select, 'ArrowDown');
    await h.dispatch('change', h.controls.select, { atTarget: () => { h.state.labels = 'both'; } });
    await h.dispatch('focusout', h.controls.select);
    await settleTasks();
    assert.deepEqual(h.counts(), [1, 0]);
    await h.dispatch('focusin', h.controls.select);
    await h.key(h.controls.select, 'z', { ctrlKey: true });
    await settleTasks();
    assert.equal(h.state.labels, 'out');
    await h.key(h.controls.select, 'z', { ctrlKey: true, shiftKey: true });
    await settleTasks();
    assert.equal(h.state.labels, 'both');
    await h.key(h.controls.select, 'y', { ctrlKey: true });
    await h.key(h.controls.select, 'z', { metaKey: true });
    await settleTasks();
    assert.equal(h.state.labels, 'out');
    await h.dispatch('focusout', h.controls.select);
    await settleTasks();
    assert.deepEqual(h.counts(), [0, 1]);
  } finally { h.cleanup(); }
});

test('text fields keep the browser undo for Ctrl+Z', async () => {
  const h = await createBoundaryHarness();
  try {
    await h.typeText('audit');
    await h.key(h.controls.text, 'z', { ctrlKey: true });
    await h.key(h.controls.text, 'y', { ctrlKey: true });
    await settleTasks();
    assert.equal(h.state.prefix, 'audit');
    assert.deepEqual(h.counts(), [0, 0]);
  } finally { h.cleanup(); }
});

test('N-18: a gesture owner settles the focused control before its own transaction', async () => {
  const h = await createBoundaryHarness();
  try {
    await h.typeText('audit');
    const gesture = await h.history.begin('Move legend', { source: 'legend-drag', owner: Symbol('legend-drag') });
    assert.deepEqual(h.counts(), [1, 0]);
    h.state.mode = 'moved';
    await h.history.commit(gesture);
    assert.deepEqual(h.counts(), [2, 0]);
    // An ownerless action joins the open transaction instead of splitting it.
    const owned = await h.history.begin('Owned edit', { owner: 'owner-a' });
    const joined = await h.history.begin('Joined edit');
    assert.equal(joined, owned);
    h.history.cancel(owned);
    await h.history.undo();
    assert.equal(h.state.mode, 'circular');
    assert.equal(h.state.prefix, 'audit');
  } finally { h.cleanup(); }
});
