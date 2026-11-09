// Every `<input>` and `<textarea>` the user types into shows its model's value
// through `v-committed-value` (js/app/committed-value.js), never through a value
// binding: Vue patches a bound `value` on every render of the owning component
// and puts the model's value back over the text the user has typed but not
// committed, or whose commit it refused (OV-3068). This guard reads every
// template (the root, the x-templates in index.html, and the inline templates
// of js/), and tests the directive's rule with a stand-in field.
import assert from 'node:assert/strict';
import { readFileSync, readdirSync } from 'node:fs';
import { join } from 'node:path';
import test from 'node:test';

import { committedText, committedValue, keepsTypedText, settleCommittedField } from '../../gbdraw/web/js/app/committed-value.js';

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
// A value binding in any of Vue's forms: `:value`, `v-bind:value`, `.value`
// (the `.prop` shorthand), with or without modifiers.
const BINDS_VALUE = /(?:^|\s)(?::|v-bind:|\.)value(?:\.[a-z]+)*=/;
// An object binding (`v-bind="{ value: x }"`, `v-bind="attrs"`) may set the value, and the guard cannot read it.
const BINDS_OBJECT = /(?:^|\s)v-bind="/;
const DIRECTIVE = /(?:^|\s)v-committed-value=/;
// A field with a static `readonly` takes no typing, so no text of the user's can be lost (the Copy diagnostics boxes).
const READONLY = /\sreadonly(?:[\s>=/]|$)/;

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
    return BINDS_VALUE.test(field.text) || BINDS_OBJECT.test(field.text) ? [`${where} binds a value besides v-committed-value`] : [];
  }
  if (!isTyped(field) || READONLY.test(field.text)) return [];
  if (BINDS_VALUE.test(field.text)) return [`${where} shows a bound value; use v-committed-value`];
  return BINDS_OBJECT.test(field.text) ? [`${where} binds an object the guard cannot read; bind the value with v-committed-value`] : [];
}));

test('a field the user types into shows its value through v-committed-value', () => {
  assert.deepEqual(committedValueProblems(SOURCES), []);
  // The directive is in use, so a template that lost it is not read as clean.
  assert.ok(SOURCES.some(([, source]) => DIRECTIVE.test(source)));
});

test('the guard finds every form of a value binding on a typed field, and only there', () => {
  const problems = (source) => committedValueProblems([['t.html', source]]);
  for (const source of [
    '<input :value="item.id" @change="rename($event.target.value)">',
    '<input type="number"\n  :value="n"\n  @blur="set($event.target.value)">',
    '<textarea :value="text" @change="read($event.target.value)"></textarea>',
    '<input type="color" :value="c" @change="paint(c => c > 0)">',
    '<input :value="a" @keydown.enter.prevent="commit()">',
    // A field whose model follows each keystroke loses text too: a refused keystroke is put back by the next render.
    '<input type="number" :value="row" @input="setRow($event.target.value)">',
    '<input v-bind:value="a">',
    '<input :value.prop="a" @focusout="commit()">',
    '<input .value="a" v-on="{ change: commit }">',
    '<input v-bind="{ value: a }" @keydown="commitOnEnter">',
    '<input v-bind="attrs">',
    '<input :value="a" v-committed-value="a" @change="commit()">',
    '<input v-committed-value="a" v-bind="{ value: a }">'
  ]) assert.equal(problems(source).length, 1, source);
  for (const source of [
    '<input v-committed-value="item.id" @change="rename($event.target.value)">',
    '<input type="radio" :value="choice" :checked="on" @change="pick(choice)">',
    '<input type="checkbox" :value="choice" v-model="picked">',
    '<select :value="mode" @change="setMode($event.target.value)"></select>',
    '<textarea :value="diagnostics" readonly></textarea>',
    '<input :aria-valuenow="a" :data-value="a">'
  ]) assert.deepEqual(problems(source), [], source);
});

test('a render keeps the typed text while the field has focus or a commit runs and the value is unchanged, or while composing', () => {
  const keeps = (focused, committing, composing, changed) => keepsTypedText({ focused, committing, composing, changed });
  assert.equal(keeps(true, false, false, false), true);
  assert.equal(keeps(false, true, false, false), true, 'a running commit keeps the text after blur');
  assert.equal(keeps(true, true, false, true), false, 'a changed value replaces the text');
  assert.equal(keeps(false, false, false, false), false);
  assert.equal(keeps(false, false, true, true), true);
  assert.equal(committedText(null), '');
  assert.equal(committedText(undefined), '');
  assert.equal(committedText(12.5), '12.5');
});

/** A stand-in for an `<input>` or a `<textarea>`: its value, focus, connection, and listeners. */
const standInField = (tagName = 'INPUT') => {
  /** @type {Map<string, Set<(event: object) => void>>} */
  const listeners = new Map();
  const ownerDocument = { activeElement: /** @type {unknown} */ (null) };
  const field = {
    tagName,
    value: 'typed by the user',
    isConnected: true,
    ownerDocument,
    addEventListener: (type, listener) => { listeners.set(type, (listeners.get(type) ?? new Set()).add(listener)); },
    removeEventListener: (type, listener) => { listeners.get(type)?.delete(listener); },
    dispatch: (type, event = {}) => { for (const listener of listeners.get(type) ?? []) listener({ type, ...event }); },
    focus: () => { ownerDocument.activeElement = field; },
    blur: () => { ownerDocument.activeElement = null; field.dispatch('blur'); },
    listenerCount: () => [...listeners.values()].reduce((sum, set) => sum + set.size, 0)
  };
  return field;
};
/** Vue's next tick, and a wait for every queued microtask. */
const useVue = () => { globalThis.window = { Vue: { nextTick: (callback) => Promise.resolve().then(callback) } }; };
const flush = () => new Promise((resolve) => { setImmediate(resolve); });
/** A commit that ends when the test says so. */
const deferred = () => {
  let resolve;
  const promise = new Promise((done) => { resolve = done; });
  return { promise, resolve };
};

test('the directive writes the model value at mount, on a render, after a commit ends, and stops at unmount', async () => {
  useVue();
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
  field.dispatch('keydown', { key: 'Enter', isComposing: true });
  await flush();
  assert.equal(field.value, 'Repeat', 'the Enter that ends a composition commits nothing');
  field.dispatch('compositionend');

  committedValue.updated(field, { value: 'renamed' });
  assert.equal(field.value, 'renamed', 'a changed value replaces the text, also in the focused field');

  field.value = '  renamed  ';
  field.dispatch('change');
  assert.equal(field.value, '  renamed  ');
  await flush();
  assert.equal(field.value, 'renamed', 'after a refused commit the focused field shows the stored value');

  field.value = 'typed again';
  field.blur();
  committedValue.updated(field, { value: 'renamed' });
  assert.equal(field.value, 'renamed', 'a render writes the value into a field without focus');

  committedValue.updated(field, { value: null });
  assert.equal(field.value, '');

  committedValue.unmounted(field);
  assert.equal(field.listenerCount(), 0);
});

test('an asynchronous commit keeps the typed text until it settles, then shows what it stored or kept', async () => {
  useVue();
  for (const stored of ['Kinases', 'Kinase']) {
    const field = standInField();
    committedValue.beforeMount(field, { value: 'Kinase' });
    field.focus();
    field.value = 'Kinases';
    const commit = deferred();
    assert.equal(settleCommittedField(field, commit.promise), commit.promise);
    field.dispatch('change');
    field.blur();
    committedValue.updated(field, { value: 'Kinase' });
    await flush();
    assert.equal(field.value, 'Kinases', 'a render and the blur while the commit runs keep the typed text');
    committedValue.updated(field, { value: stored });
    commit.resolve();
    await flush();
    assert.equal(field.value, stored, 'the settled commit shows the stored value, or the kept one after a refusal');
  }
});

test('a commit by Enter shows the stored value with focus kept, unless the user typed on, and a textarea Enter commits nothing', async () => {
  useVue();
  const field = standInField();
  committedValue.beforeMount(field, { value: 'CDS' });
  field.focus();
  field.value = '';
  field.dispatch('keydown', { key: 'Enter', isComposing: false });
  await flush();
  assert.equal(field.value, 'CDS', 'a refused Enter commit shows the stored name while the field keeps focus');

  const commit = deferred();
  field.value = 'gene';
  settleCommittedField(field, commit.promise);
  field.dispatch('keydown', { key: 'Enter', isComposing: false });
  field.value = 'gene2';
  commit.resolve();
  await flush();
  assert.equal(field.value, 'gene2', 'text typed after the commit started stays');

  const area = standInField('TEXTAREA');
  committedValue.beforeMount(area, { value: '>q' });
  area.focus();
  area.value = '>q\nACGT';
  area.dispatch('keydown', { key: 'Enter', isComposing: false });
  await flush();
  assert.equal(area.value, '>q\nACGT');
});

test('a field the directive has not mounted fails loudly', () => {
  useVue();
  const field = standInField();
  assert.throws(() => settleCommittedField(field, Promise.resolve()), /not bound with v-committed-value/);
  assert.throws(() => settleCommittedField(null, Promise.resolve()), /without its field/);
  assert.throws(() => committedValue.updated(field, { value: 'x' }), /not bound with v-committed-value/);
  assert.throws(() => committedValue.unmounted(field), /not bound with v-committed-value/);
});
