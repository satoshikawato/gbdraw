// A text, number, or color field that commits on change, blur, or Enter shows
// its model's value through `v-committed-value` (js/app/committed-value.js),
// never through `:value`: Vue patches a bound `value` on every render of the
// owning component and puts the model's value back over the text the user has
// typed but not committed (OV-3068). This guard reads every template (the root,
// the x-templates in index.html, and the inline templates of js/), and tests
// the directive's rule with a stand-in field.
import assert from 'node:assert/strict';
import { readFileSync, readdirSync } from 'node:fs';
import { join } from 'node:path';
import test from 'node:test';

import { committedText, committedValue, keepsTypedText } from '../../gbdraw/web/js/app/committed-value.js';

const WEB = join(process.cwd(), 'gbdraw', 'web');
const SOURCES = [
  ['index.html', readFileSync(join(WEB, 'index.html'), 'utf8')],
  ...readdirSync(join(WEB, 'js'), { recursive: true })
    .filter((path) => String(path).endsWith('.js'))
    .map((path) => [`js/${path}`, readFileSync(join(WEB, 'js', String(path)), 'utf8')])
];

// Input types whose `value` is not text the user types: a checkbox's or a
// radio's is its option value (its state is `checked`), and a range follows
// the pointer, as Vue's own v-model treats it.
const NOT_TYPED = new Set(['checkbox', 'radio', 'file', 'hidden', 'range', 'button', 'submit', 'reset', 'image']);
const BINDS_VALUE = /(?:^|\s)(?::|v-bind:)value=/;
const COMMITS = /(?:^|\s)(?:@|v-on:)(?:change|blur|keydown\.enter|keyup\.enter)\b/;
const DIRECTIVE = /(?:^|\s)v-committed-value=/;

/** The opening tags of `<input>` and `<textarea>` in `source`, skipping '>' inside quoted attribute values. */
const fieldTags = (source) => [...source.matchAll(/<(input|textarea)\b/g)].map((match) => {
  let quote = null;
  let index = match.index;
  for (; index < source.length; index += 1) {
    const char = source[index];
    if (quote) {
      if (char === quote) quote = null;
    } else if (char === '"' || char === "'") {
      quote = char;
    } else if (char === '>') {
      break;
    }
  }
  return { tag: match[1], text: source.slice(match.index, index + 1), line: source.slice(0, match.index).split('\n').length };
});

const isTyped = ({ tag, text }) => tag === 'textarea' || !NOT_TYPED.has(/\stype="([^"]*)"/.exec(text)?.[1] ?? 'text');

/** The fields of `sources` that break the rule, as `path:line`. */
const committedValueProblems = (sources) => sources.flatMap(([path, source]) => fieldTags(source).flatMap((field) => {
  const where = `${path}:${field.line}`;
  if (DIRECTIVE.test(field.text)) {
    return BINDS_VALUE.test(field.text) ? [`${where} binds both :value and v-committed-value`] : [];
  }
  return isTyped(field) && BINDS_VALUE.test(field.text) && COMMITS.test(field.text)
    ? [`${where} commits on change, blur, or Enter but shows :value; use v-committed-value`]
    : [];
}));

test('a field that commits on change, blur, or Enter shows its value through v-committed-value', () => {
  assert.deepEqual(committedValueProblems(SOURCES), []);
  // The directive is in use, so a template that lost it is not read as clean.
  assert.ok(SOURCES.some(([, source]) => DIRECTIVE.test(source)));
});

test('the guard finds :value on a field that commits, and only there', () => {
  const problems = (source) => committedValueProblems([['t.html', source]]);
  assert.equal(problems('<input :value="item.id" @change="rename($event.target.value)">').length, 1);
  assert.equal(problems('<input type="number"\n  :value="n"\n  @blur="set($event.target.value)">').length, 1);
  assert.equal(problems('<textarea :value="text" @change="read($event.target.value)"></textarea>').length, 1);
  assert.equal(problems('<input type="color" :value="c" @change="paint(c => c > 0)">').length, 1);
  assert.equal(problems('<input :value="a" @keydown.enter.prevent="commit()">').length, 1);
  assert.equal(problems('<input :value="a" v-committed-value="a" @change="commit()">').length, 1);
  assert.deepEqual(problems('<input v-committed-value="item.id" @change="rename($event.target.value)">'), []);
  assert.deepEqual(problems('<input :value="label" @input="setLabel($event.target.value)">'), []);
  assert.deepEqual(problems('<input type="radio" :value="choice" :checked="on" @change="pick(choice)">'), []);
  assert.deepEqual(problems('<select :value="mode" @change="setMode($event.target.value)"></select>'), []);
});

test('a render keeps the typed text only while the field has focus and the value is unchanged, or while composing', () => {
  assert.equal(keepsTypedText({ focused: true, composing: false, changed: false }), true);
  assert.equal(keepsTypedText({ focused: true, composing: false, changed: true }), false);
  assert.equal(keepsTypedText({ focused: false, composing: false, changed: false }), false);
  assert.equal(keepsTypedText({ focused: false, composing: true, changed: true }), true);
  assert.equal(committedText(null), '');
  assert.equal(committedText(undefined), '');
  assert.equal(committedText(12.5), '12.5');
});

/** A stand-in for an `<input>`: its value, focus, connection, and listeners. */
const standInField = () => {
  /** @type {Map<string, Set<() => void>>} */
  const listeners = new Map();
  const ownerDocument = { activeElement: /** @type {unknown} */ (null) };
  const field = {
    value: 'typed by the user',
    isConnected: true,
    ownerDocument,
    addEventListener: (type, listener) => { listeners.set(type, (listeners.get(type) ?? new Set()).add(listener)); },
    removeEventListener: (type, listener) => { listeners.get(type)?.delete(listener); },
    dispatch: (type) => { for (const listener of listeners.get(type) ?? []) listener(); },
    focus: () => { ownerDocument.activeElement = field; },
    blur: () => { ownerDocument.activeElement = null; field.dispatch('blur'); },
    listenerCount: () => [...listeners.values()].reduce((sum, set) => sum + set.size, 0)
  };
  return field;
};

test('the directive writes the model value at mount, on a render, after blur, and stops at unmount', async () => {
  globalThis.window = { Vue: { nextTick: (callback) => Promise.resolve().then(callback) } };
  const field = standInField();
  committedValue.beforeMount(field, { value: 'region_1' });
  assert.equal(field.value, 'region_1');

  field.focus();
  field.value = 'Repeat';
  committedValue.updated(field, { value: 'region_1' });
  assert.equal(field.value, 'Repeat', 'a render keeps the text typed into the focused field');

  field.dispatch('compositionstart');
  committedValue.updated(field, { value: 'other' });
  assert.equal(field.value, 'Repeat', 'a render leaves an open IME composition');
  field.dispatch('compositionend');

  committedValue.updated(field, { value: 'renamed' });
  assert.equal(field.value, 'renamed', 'a changed value replaces the text, also in the focused field');

  field.value = '  renamed  ';
  field.blur();
  assert.equal(field.value, '  renamed  ');
  await Promise.resolve();
  await Promise.resolve();
  assert.equal(field.value, 'renamed', 'after blur the field shows the value its commit left');

  field.value = 'typed again';
  committedValue.updated(field, { value: 'renamed' });
  assert.equal(field.value, 'renamed', 'a render writes the value into a field without focus');

  committedValue.updated(field, { value: null });
  assert.equal(field.value, '');

  committedValue.unmounted(field);
  assert.equal(field.listenerCount(), 0);
});
