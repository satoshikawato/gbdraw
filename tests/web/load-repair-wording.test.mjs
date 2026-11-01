// The Load repair dialog's sentences (utils/load-repairs.js, the one wording
// owner): each is one catalog key with placeholders, never a sentence built
// from a translated word (docs/internal/I18N_STYLE.md 1.3), and a list in a
// sentence is shown with the locale's list pattern (`listText`).
import assert from 'node:assert/strict';
import test from 'node:test';

import { SOURCE_LOCALE, setUiLocale } from '../../gbdraw/web/js/utils/i18n.js';
import { describeLoadRepairs, loadRepair } from '../../gbdraw/web/js/utils/load-repairs.js';

/** @param {import('../../gbdraw/web/js/utils/load-repairs.js').LoadRepair[]} repairs */
const items = (repairs) => describeLoadRepairs(repairs.map(loadRepair))?.items ?? [];

const PART = 'renderRequest.diagramOptions.colors.defaultColorsFile';
const REFERENCE = 'renderRequest.comparisons[0]';
const OTHER = 'otherModeResult.renderRequest';
const repairsWithoutName = () => [
  { code: 'other-result-set-dropped', path: OTHER, effect: 'dropped' },
  { code: 'request-part-dropped', path: PART, effect: 'defaulted' },
  { code: 'request-part-dropped', path: `${PART}.inner`, effect: 'dropped' },
  { code: 'request-reference-dropped', path: REFERENCE, effect: 'dropped' }
];

test('a repair without a name is one whole sentence; English is as before', () => {
  assert.deepEqual(items(repairsWithoutName()), [
    `Saved other diagram: its settings (${OTHER}) could not be read; that diagram and its Results are left out.`,
    `Settings: a part (${PART}) could not be read; its default is used.`,
    `Settings: a part (${PART}.inner) could not be read and was left out, with the settings inside it.`,
    `Settings: ${REFERENCE} names a file or record this Session does not hold; it is left out.`
  ]);
  assert.deepEqual(items([{ code: 'other-result-set-dropped', path: OTHER, effect: 'dropped', mode: 'linear' }]), [
    `Saved Linear diagram: its settings (${OTHER}) could not be read; that diagram and its Results are left out.`
  ]);
});

test('in Japanese, a repair sentence is translated whole, with no English word or fragment inside', async () => {
  try {
    await setUiLocale('ja');
    assert.deepEqual(items(repairsWithoutName()), [
      `保存されたもう一方の図：その設定(${OTHER})を読み取れなかったため、この図とその結果は除外されました。`,
      `設定：一部(${PART})を読み取れなかったため、既定値が使用されます。`,
      `設定：一部(${PART}.inner)を読み取れなかったため、その中の設定とともに除外されました。`,
      `設定：${REFERENCE}は、このセッションに存在しないファイルまたはレコードを指しているため、除外されます。`
    ]);
    assert.deepEqual(items([{ code: 'other-result-set-dropped', path: OTHER, effect: 'dropped', mode: 'circular' }]), [
      `保存された環状図：その設定(${OTHER})を読み取れなかったため、この図とその結果は除外されました。`
    ]);
  } finally {
    await setUiLocale(SOURCE_LOCALE);
  }
});

test('a Session part with no wording of its own names its path once', () => {
  assert.deepEqual(items([{ code: 'session-part-dropped', path: 'someSessionPart', effect: 'dropped' }]), [
    'Session file: someSessionPart could not be read and is left out.'
  ]);
});

test('the names in a repair sentence are a list: "a and b" in English, the locale pattern in Japanese', async () => {
  const fields = (count) => Array.from({ length: count }, (_, index) => (
    { code: 'session-field-dropped', path: `field${index + 1}`, effect: 'dropped' }
  ));
  assert.deepEqual(items(fields(2)), ['Session file: the fields field1 and field2 could not be read and were left out.']);
  assert.deepEqual(items(fields(3)), ['Session file: the fields field1, field2, and field3 could not be read and were left out.']);
  assert.deepEqual(items(fields(7)), [
    'Session file: the fields field1, field2, field3, field4, field5, and 2 more could not be read and were left out.'
  ]);
  try {
    await setUiLocale('ja');
    const [line] = items(fields(7));
    assert.ok(line.includes('field1、field2、field3、field4、field5、ほか2件'), line);
    assert.doesNotMatch(line, /, |\band\b/);
  } finally {
    await setUiLocale(SOURCE_LOCALE);
  }
});
