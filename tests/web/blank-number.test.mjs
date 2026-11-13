// A number field bound with `v-model` writes null when it is blank, never ''
// (js/app/blank-number.js, OV-3351). Vue writes '' for a cleared
// `v-model.number` field and for a cleared `<input type="number" v-model>`,
// so every such field carries `v-blank-null` right after its `v-model`, whose
// listener runs first. This guard reads every template (the root, the
// x-templates in index.html, and the inline templates of js/), with no
// exception list, and tests the directive's rule with a stand-in field.
import assert from 'node:assert/strict';
import { readFileSync, readdirSync } from 'node:fs';
import { join } from 'node:path';
import test from 'node:test';

import { blankNull } from '../../gbdraw/web/js/app/blank-number.js';

const WEB = join(process.cwd(), 'gbdraw', 'web');
const SOURCES = [
  ['index.html', readFileSync(join(WEB, 'index.html'), 'utf8')],
  ...readdirSync(join(WEB, 'js'), { recursive: true })
    .filter((path) => String(path).endsWith('.js'))
    .map((path) => [`js/${path}`, readFileSync(join(WEB, 'js', String(path)), 'utf8')])
];

/** The opening tags of `<input>` in `source`, skipping '>' inside quoted attribute values. */
const inputTags = (source) => [...source.matchAll(/<input\b/g)].map((match) => {
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
  return { text: source.slice(match.index, index + 1), line: source.slice(0, match.index).split('\n').length };
});

const MODEL = /(?:^|\s)v-model(\.[a-z.]+)?=/;
const DIRECTIVE = /(?:^|\s)v-blank-null(?=[\s>/]|$)/;

/** The modifiers of the field's `v-model` (`.lazy.number` gives `['lazy', 'number']`). */
const modelModifiers = (text) => (MODEL.exec(text)?.[1] ?? '').split('.').filter(Boolean);

/** Whether Vue writes a number, or '' when blank, for the field's `v-model`. */
const writesNumber = (text) => MODEL.test(text) && (modelModifiers(text).includes('number') || /\stype="number"/.test(text));

/** The fields of `sources` that break the rule, as `path:line`. */
const blankNumberProblems = (sources) => sources.flatMap(([path, source]) => inputTags(source).flatMap(({ text, line }) => {
  const where = `${path}:${line}`;
  const directive = DIRECTIVE.exec(text);
  if (!writesNumber(text)) return directive ? [`${where} has v-blank-null without a number v-model`] : [];
  if (!directive) return [`${where} writes '' when blank; add v-blank-null after its v-model`];
  // `.lazy` writes on change, after the directive's input listener ran: the blank would stay ''.
  if (modelModifiers(text).includes('lazy')) return [`${where} binds v-model.lazy, which writes '' on change after v-blank-null`];
  return directive.index < MODEL.exec(text).index ? [`${where} names v-blank-null before its v-model`] : [];
}));

test('every number field bound with v-model writes null when blank', () => {
  assert.deepEqual(blankNumberProblems(SOURCES), []);
  // The fields are there, so a template that lost them is not read as clean.
  const count = SOURCES.reduce((sum, [, source]) => sum + inputTags(source).filter(({ text }) => writesNumber(text)).length, 0);
  assert.ok(count >= 83, `${count} number fields`);
});

test('the guard finds a number field without the directive, and only there', () => {
  const problems = (source) => blankNumberProblems([['t.html', source]]);
  for (const source of [
    '<input type="number" v-model.number="adv.window_size">',
    '<input v-model.number="adv.window_size" class="x">',
    '<input type="number" min="0" v-model="panel.identity">',
    '<input type="number"\n  v-model.lazy.number="a">',
    '<input type="number" v-model.lazy.number="a" v-blank-null>',
    '<input type="number" v-model.lazy="a" v-blank-null>',
    '<input type="number" v-blank-null v-model.number="a">',
    '<input type="text" v-model="title" v-blank-null>'
  ]) assert.equal(problems(source).length, 1, source);
  for (const source of [
    '<input type="number" v-model.number="adv.window_size" v-blank-null>',
    '<input type="number" min="0" v-model="panel.identity" v-blank-null class="x">',
    '<input type="text" v-model="title">',
    '<input type="number" v-committed-value="a" @change="set($event.target.value)">',
    '<select v-model.number="slot.track_index"></select>'
  ]) assert.deepEqual(problems(source), [], source);
});

test('the directive writes null after v-model wrote the blank, through the assignment of the latest render', () => {
  const writes = [];
  /** @type {(() => void)[]} */
  const listeners = [];
  const field = { value: '', addEventListener: (_type, listener) => listeners.push(listener) };
  const vnode = (name) => ({ props: { 'onUpdate:modelValue': (value) => writes.push([name, value]) } });
  blankNull.created(field, null, vnode('first'));
  field.value = '12';
  listeners.forEach((listener) => listener());
  assert.deepEqual(writes, [], 'a number is v-model\'s to write');
  blankNull.beforeUpdate(field, null, vnode('second'));
  field.value = '';
  listeners.forEach((listener) => listener());
  assert.deepEqual(writes, [['second', null]]);
  blankNull.unmounted(field);
  listeners.forEach((listener) => listener());
  assert.deepEqual(writes, [['second', null]], 'an unmounted field writes nothing');
});
