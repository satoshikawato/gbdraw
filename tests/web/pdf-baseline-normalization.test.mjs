import assert from 'node:assert/strict';
import { test } from 'node:test';
import { normalizeTextBaselinesForPdf } from '../../gbdraw/web/js/services/export.js';

// The per-text algorithm the batched phases replace (read, write, read, write for each text).
const normalizeSequentially = (svg) => {
  svg.querySelectorAll('text').forEach((textEl) => {
    const baseline = textEl.getAttribute('dominant-baseline');
    if (!baseline || baseline === 'alphabetic' || baseline === 'auto') return;
    let before;
    try { before = textEl.getBBox(); } catch (error) { return; }
    textEl.setAttribute('dominant-baseline', 'alphabetic');
    let after;
    try { after = textEl.getBBox(); } catch (error) {
      textEl.setAttribute('dominant-baseline', baseline);
      return;
    }
    const dy = before.y - after.y;
    if (!Number.isFinite(dy) || Math.abs(dy) < 0.01) {
      textEl.removeAttribute('dominant-baseline');
      return;
    }
    const yAttr = textEl.getAttribute('y');
    if (yAttr) {
      const adjusted = yAttr.split(/[\s,]+/).map((v) => parseFloat(v)).map((v) => (Number.isFinite(v) ? v + dy : v));
      textEl.setAttribute('y', adjusted.join(' '));
    } else {
      const dyAttr = textEl.getAttribute('dy');
      const dyValue = dyAttr ? parseFloat(dyAttr) : 0;
      textEl.setAttribute('dy', Number.isFinite(dyValue) ? dyValue + dy : dy);
    }
    textEl.removeAttribute('dominant-baseline');
  });
};

const BASELINE_SHIFT = { middle: -4, central: -3, hanging: -8, alphabetic: 0 };

// A text whose box moves with its baseline; `layout` counts getBBox calls that follow a write.
const makeText = (attrs, { throwBefore = false, throwAfter = false } = {}, layout) => {
  const store = new Map(Object.entries(attrs));
  return {
    getAttribute: (name) => (store.has(name) ? store.get(name) : null),
    setAttribute: (name, value) => { store.set(name, String(value)); layout.dirty = true; },
    removeAttribute: (name) => { store.delete(name); layout.dirty = true; },
    getBBox() {
      const baseline = store.get('dominant-baseline') || 'alphabetic';
      if (baseline === 'alphabetic' ? throwAfter : throwBefore) throw new Error('not rendered');
      if (layout.dirty) { layout.forced += 1; layout.dirty = false; }
      const y = parseFloat(String(store.get('y') || '0').split(/[\s,]+/)[0]) + (parseFloat(store.get('dy')) || 0);
      return { y: y + BASELINE_SHIFT[baseline] };
    },
    snapshot: () => Object.fromEntries([...store.entries()].sort())
  };
};

const build = (specs) => {
  const layout = { dirty: false, forced: 0 };
  const texts = specs.map(([attrs, flags]) => makeText(attrs, flags, layout));
  const svg = { querySelectorAll: () => ({ forEach: (fn) => texts.forEach(fn) }) };
  return { svg, texts, layout };
};

const MIX = [
  [{ 'dominant-baseline': 'middle', y: '10' }],
  [{ 'dominant-baseline': 'central', y: '20, 30  40' }],
  [{ 'dominant-baseline': 'hanging', dy: '2' }],
  [{ 'dominant-baseline': 'hanging' }],
  [{ 'dominant-baseline': 'middle', y: '5' }, { throwBefore: true }],
  [{ 'dominant-baseline': 'middle', y: '6' }, { throwAfter: true }],
  [{ 'dominant-baseline': 'alphabetic', y: '7' }],
  [{ 'dominant-baseline': 'auto', y: '8' }],
  [{ y: '9' }],
  [{ 'dominant-baseline': 'central', y: 'abc 3' }]
];

test('batched baseline normalization produces the per-text sequential attributes', () => {
  const batched = build(MIX);
  const sequential = build(MIX);
  normalizeTextBaselinesForPdf(batched.svg);
  normalizeSequentially(sequential.svg);
  assert.deepEqual(batched.texts.map((t) => t.snapshot()), sequential.texts.map((t) => t.snapshot()));
  assert.equal(batched.texts[0].snapshot().y, '6');
  assert.equal(batched.texts[3].snapshot().dy, '-8');
  assert.equal(batched.texts[4].snapshot()['dominant-baseline'], 'middle', 'a text that cannot be measured before keeps its baseline');
  assert.equal(batched.texts[5].snapshot()['dominant-baseline'], 'middle', 'a text that cannot be measured after gets its baseline back');
});

test('baseline normalization forces at most two layouts however many texts there are', () => {
  const texts = Array.from({ length: 50 }, (_, i) => [{ 'dominant-baseline': i % 2 ? 'middle' : 'hanging', y: String(i) }]);
  const batched = build(texts);
  const sequential = build(texts);
  normalizeTextBaselinesForPdf(batched.svg);
  normalizeSequentially(sequential.svg);
  assert.ok(batched.layout.forced <= 2, `forced layouts: ${batched.layout.forced}`);
  assert.ok(sequential.layout.forced >= 99, 'the per-text order forces two layouts per text');
  assert.deepEqual(batched.texts.map((t) => t.snapshot()), sequential.texts.map((t) => t.snapshot()));
});

test('a text that cannot be measured after does not add a layout', () => {
  const texts = [
    ...Array.from({ length: 10 }, (_, i) => [{ 'dominant-baseline': 'middle', y: String(i) }]),
    [{ 'dominant-baseline': 'middle', y: '1' }, { throwAfter: true }]
  ];
  const batched = build(texts);
  normalizeTextBaselinesForPdf(batched.svg);
  assert.ok(batched.layout.forced <= 2, `forced layouts: ${batched.layout.forced}`);
});
