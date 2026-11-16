// OV-3676: the Copy diagnostics outcome of an error card stays shown while the UI language may change, so the card
// keeps how to word it and shows it in the language shown now; the template shows `copyStatus` as it is.
import assert from 'node:assert/strict';
import test from 'node:test';
import { loadPageVue } from './helpers/vue-runtime.mjs';

const Vue = loadPageVue();
globalThis.window = { Vue, setTimeout, clearTimeout, innerWidth: 1440 };
const { OperationError } = await import('../../gbdraw/web/js/components.js');
const { SOURCE_LOCALE, setUiLocale } = await import('../../gbdraw/web/js/utils/i18n.js');
const { default: ja } = await import('../../gbdraw/web/js/utils/i18n-catalogs/ja.js');

/** The card's setup on an error, with its Details open, and a clipboard that writes or refuses. */
const card = (writeText) => {
  Object.defineProperty(globalThis, 'navigator', { configurable: true, value: { clipboard: writeText ? { writeText } : undefined } });
  const props = Vue.reactive({ error: { code: 'UNKNOWN' }, recovery: null, compact: false });
  const scope = Vue.effectScope();
  const shown = scope.run(() => OperationError.setup(props));
  shown.details.value = { open: true };
  return { shown, stop: () => scope.stop() };
};

for (const [name, writeText, key] of [
  ['copied', async () => {}, 'Diagnostics copied.'],
  ['refused', null, 'Clipboard unavailable. Select diagnostics and copy manually.']
]) {
  test(`the Copy outcome (${name}) shows in the language shown now`, async () => {
    const { shown, stop } = card(writeText);
    try {
      await shown.copyDiagnostics();
      assert.equal(shown.copyStatus.value, key);
      await setUiLocale('ja', async () => ({ default: ja }));
      assert.equal(shown.copyStatus.value, ja[key]);
    } finally {
      await setUiLocale(SOURCE_LOCALE);
      stop();
    }
  });
}
