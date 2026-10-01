import assert from 'node:assert/strict';
import test from 'node:test';

import { setupHistoryShortcuts } from '../../gbdraw/web/js/app/history-shortcuts.js';

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

// SE-04: a select has no native Undo, so History shortcuts apply to it.
test('Undo and Redo shortcuts work while a select has focus', () => {
  const { calls, press } = installShortcuts();
  press({ tagName: 'SELECT' });
  press({ tagName: 'SELECT', shiftKey: true });
  press({ tagName: 'SELECT', key: 'y' });
  assert.deepEqual(calls, ['undo', 'redo', 'redo']);
});
