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
  '@change="importSession"'
].forEach((handler) => {
  const escapedHandler = handler.replace(/[.*+?^${}()|[\]\\]/g, '\\$&');
  assert.match(
    indexHtml,
    new RegExp(`<(?:button|input)(?=[^>]*${escapedHandler})(?=[^>]*data-history-managed)[^>]*>`),
    `${handler} must bypass the generic input adapter and retain its explicit History boundary`
  );
});

console.log('history input tests passed');
