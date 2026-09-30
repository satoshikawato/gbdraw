import assert from 'node:assert/strict';
import test from 'node:test';

import { setupHistoryShortcuts } from '../../gbdraw/web/js/app/history-shortcuts.js';
import { assertKnownDefect } from './helpers/known-defect.mjs';

const installShortcuts = () => {
  const listeners = new Map();
  globalThis.document = {
    addEventListener: (type, listener) => listeners.set(type, listener),
    removeEventListener: (type) => listeners.delete(type)
  };
  const calls = [];
  setupHistoryShortcuts({
    history: { undo: () => calls.push('undo'), redo: () => calls.push('redo') },
    onMounted: (callback) => callback(),
    onUnmounted: () => {}
  });
  const press = ({ tagName, type = null, key = 'z', shiftKey = false }) => {
    const target = { tagName, type, isContentEditable: false };
    target.closest = () => target;
    listeners.get('keydown')({
      key, shiftKey, ctrlKey: true, metaKey: false, altKey: false, target, preventDefault: () => {}
    });
  };
  return { calls, press };
};

test('text entry keeps the browser undo while checkbox focus uses History', () => {
  const { calls, press } = installShortcuts();
  press({ tagName: 'INPUT', type: 'text' });
  press({ tagName: 'TEXTAREA' });
  assert.deepEqual(calls, []);
  press({ tagName: 'INPUT', type: 'checkbox' });
  press({ tagName: 'INPUT', type: 'checkbox', key: 'y' });
  assert.deepEqual(calls, ['undo', 'redo']);
});

test('Undo and Redo shortcuts work while a select has focus (SE-04 known defect)', async () => {
  const { calls, press } = installShortcuts();
  await assertKnownDefect('SE-04', () => {
    press({ tagName: 'SELECT' });
    press({ tagName: 'SELECT', shiftKey: true });
    press({ tagName: 'SELECT', key: 'y' });
    assert.deepEqual(calls, ['undo', 'redo', 'redo']);
  });
});
