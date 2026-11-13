import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import test from 'node:test';
import { gunzipSync } from 'node:zlib';

import { LIVE_EDIT_DOMAINS } from '../../gbdraw/web/js/app/candidate-render.js';
import { admitFeatureCatalog } from '../../gbdraw/web/js/services/feature-catalog.js';
import { normalizeLogicalResults } from '../../gbdraw/web/js/services/result-normalization.js';
import {
  createEmptySvgMutationPlan, createSavedResultPlan, reconcileMountedResult
} from '../../gbdraw/web/js/services/svg-result-ingestion.js';
import { getLegendEntrySwatch as legendSwatch, pythonLegendKey } from '../../gbdraw/web/js/services/legend-svg.js';
import { sessionSourceFactsOf } from '../../gbdraw/web/js/services/session-drawings-upgrade.js';

// A Legend row's key is read in one place, `pythonLegendKey`: the key record a
// rename keeps, else the row's `data-legend-key`. The Result executor
// (services/svg-result-ingestion.js) reads every row that way, an editor row
// included, which is the same key only while no editor row carries a key
// record. These cases pin that invariant and the places the executor reads the
// key: its row index, the placement of an added row, and the reconcile that
// removes an editor row no add names. The helpers below are copies of the
// smallest set svg-result-ingestion.test.mjs uses.

// Session Load (services/config.js) reads the editor state through Vue state.
globalThis.window = { Vue: {
  ref: (value) => ({ value }), reactive: (value) => value, toRaw: (value) => value,
  computed: (getter) => ({ get value() { return getter(); } }), nextTick: async () => {}
} };
const { savedResultEdits } = await import('../../gbdraw/web/js/services/config.js');

class FakeElement {
  constructor(tagName, attributes = {}, children = []) {
    this.tagName = tagName;
    this.localName = tagName;
    this.attributes = new Map(Object.entries(attributes));
    this.children = [];
    this.parentElement = null;
    this.textContent = '';
    this.style = { removeProperty() {} };
    children.forEach((child) => this.appendChild(child));
  }

  get id() { return this.getAttribute('id') || ''; }
  getAttribute(name) { return this.attributes.has(name) ? this.attributes.get(name) : null; }
  setAttribute(name, value) { this.attributes.set(name, String(value)); }
  removeAttribute(name) { this.attributes.delete(name); }
  hasAttribute(name) { return this.attributes.has(name); }
  appendChild(child) { child.remove?.(); child.parentElement = this; this.children.push(child); return child; }
  insertBefore(child, reference) {
    child.remove?.();
    child.parentElement = this;
    const at = this.children.indexOf(reference);
    this.children.splice(at < 0 ? this.children.length : at, 0, child);
    return child;
  }
  remove() {
    if (!this.parentElement) return;
    this.parentElement.children = this.parentElement.children.filter((child) => child !== this);
    this.parentElement = null;
  }
  cloneNode(deep = false) {
    const clone = new FakeElement(this.tagName, Object.fromEntries(this.attributes));
    clone.textContent = this.textContent;
    if (deep) this.children.forEach((child) => clone.appendChild(child.cloneNode(true)));
    return clone;
  }
  getElementById(id) {
    return this.walk().find((element) => element.id === id) || null;
  }
  walk() {
    return [this, ...this.children.flatMap((child) => child.walk())];
  }
  matches(selector) {
    if (/^\[[\w-]+\](?:,\s*\[[\w-]+\])*$/.test(selector)) {
      return selector.split(',').some((part) => this.hasAttribute(part.trim().slice(1, -1)));
    }
    if (selector.startsWith('.')) return false;
    if (selector === 'path') return this.tagName === 'path';
    if (selector === 'text') return this.tagName === 'text';
    if (selector === 'textPath') return this.tagName === 'textPath';
    if (selector === 'g[data-legend-key]') {
      return this.tagName === 'g' && this.hasAttribute('data-legend-key');
    }
    if (selector === 'g[data-legend-key][display="none"]') {
      return this.tagName === 'g' && this.hasAttribute('data-legend-key') && this.getAttribute('display') === 'none';
    }
    if (selector === 'text[data-label-feature-id]') {
      return this.tagName === 'text' && this.hasAttribute('data-label-feature-id');
    }
    if (selector.includes('data-gbdraw-feature-id') || selector.includes('id^="f"')) {
      return ['path', 'polygon', 'rect'].includes(this.tagName)
        && (this.hasAttribute('data-gbdraw-feature-id') || this.id.startsWith('f'));
    }
    if (selector.startsWith('#')) return this.id === selector.slice(1);
    return false;
  }
  querySelectorAll(selector) {
    return this.walk().slice(1).filter((element) => element.matches(selector));
  }
  querySelector(selector) { return this.querySelectorAll(selector)[0] || null; }
}

const serializeNode = (node) => {
  const attributes = Array.from(node.attributes.entries())
    .map(([name, value]) => ` ${name}="${value}"`).join('');
  const content = `${node.textContent || ''}${node.children.map(serializeNode).join('')}`;
  return `<${node.tagName}${attributes}>${content}</${node.tagName}>`;
};

globalThis.XMLSerializer = class {
  serializeToString(node) { return serializeNode(node); }
};

const parseSvgText = (content) => {
  const decode = (value) => value.replace(/&quot;/g, '"').replace(/&apos;/g, "'")
    .replace(/&lt;/g, '<').replace(/&gt;/g, '>').replace(/&amp;/g, '&');
  const root = new FakeElement('#document');
  let parent = root;
  const tokens = /<!--[\s\S]*?-->|<!\[CDATA\[[\s\S]*?\]\]>|<[?!][^>]*>|<\/([\w:-]+)\s*>|<([\w:-]+)((?:\s+[\w:-]+="[^"]*")*)\s*(\/?)>/g;
  for (const [token, closing, name, attributes, selfClosing] of content.matchAll(tokens)) {
    if (closing) { parent = parent.parentElement; continue; }
    if (!name) continue;
    const element = parent.appendChild(new FakeElement(name, Object.fromEntries(
      [...attributes.matchAll(/([\w:-]+)="([^"]*)"/g)].map(([, key, value]) => [key, decode(value)])
    )));
    if (!selfClosing && !token.endsWith('/>')) parent = element;
  }
  return root.children[0];
};

const legendRow = (caption, y, attributes = {}) => {
  const row = new FakeElement('g', { 'data-legend-key': caption, ...attributes });
  row.appendChild(new FakeElement('path', {
    fill: '#aaaaaa', stroke: attributes['data-test-stroke'] || '#000000', 'stroke-width': '0.5', transform: `translate(0, ${y})`
  }));
  const label = new FakeElement('text', { transform: `translate(22, ${y})` });
  label.textContent = caption;
  row.appendChild(label);
  return row;
};

const legendSvg = (rows = [legendRow('CDS', 7), legendRow('tRNA', 31), legendRow('rRNA', 55)]) => {
  const root = new FakeElement('svg', { xmlns: 'http://www.w3.org/2000/svg' });
  const legend = root.appendChild(new FakeElement('g', { id: 'legend' }));
  const featureLegend = legend.appendChild(new FakeElement('g', { id: 'feature_legend' }));
  rows.forEach((row) => featureLegend.appendChild(row));
  return root;
};

const legendOperations = (operations = {}) => ({ ...createEmptySvgMutationPlan(1).operationsByResult[0], ...operations });

const legendRows = (svg) => svg.querySelectorAll('g[data-legend-key]').map((row) => [
  row.getAttribute('data-legend-key'), row.querySelector('text').textContent, row.getAttribute('display')
]);

const sessionFixture = (name) => JSON.parse(gunzipSync(readFileSync(new URL(
  `../fixtures/sessions/${name}`, import.meta.url
))).toString('utf8'));

const loadSavedResults = (session, {
  mode = 'circular', results = normalizeLogicalResults(session.results),
  featureCatalog = session.editorState.featureCatalog, patch = () => {}
} = {}) => {
  const admission = admitFeatureCatalog(featureCatalog, results, { adopt: true, mode });
  const plan = createSavedResultPlan(results, admission, savedResultEdits(session, mode, session.editorState, session.version, sessionSourceFactsOf(session.version, false)));
  return results.map((result, index) => {
    const svg = parseSvgText(result.content);
    patch(svg, index);
    (plan.operationsByResult[index].callerTransforms || []).forEach((transform) => transform(svg));
    return svg;
  });
};

const LEGEND_STRUCTURE = LIVE_EDIT_DOMAINS.legendStructure;
const KEY_RECORD = 'data-gbdraw-base-data-legend-key';
const rowsOf = (svg) => svg.querySelectorAll('g[data-legend-key]');
const isEditor = (row) => row.getAttribute('data-legend-owner') === 'direct-editor';
const add = (caption, extra = {}) => ({ caption, color: '#556677', xPos: null, yPos: null, ...extra });
const swatchFills = (svg) => rowsOf(svg).map((row) => [pythonLegendKey(row), isEditor(row), legendSwatch(row).getAttribute('fill')]);
const editorRecords = (svg) => rowsOf(svg).filter(isEditor).map((row) => row.getAttribute(KEY_RECORD));

// The invariant: an editor row is added under its own key and no operation of
// the executor records a key on it, so `pythonLegendKey` is its key.
test('no editor row carries a key record through add, rename, delete, reconcile and Load', () => {
  const svg = legendSvg();
  const full = legendOperations({
    legendAdds: [add('Mine'), add('Other', { name: 'Shown name' })],
    legendRenames: [{ from: 'CDS', to: 'Genes', allowMissing: false }],
    legendDeletes: [{ caption: 'tRNA', allowMissing: false }]
  });
  reconcileMountedResult(svg, full);
  assert.deepEqual(rowsOf(svg).map((row) => [row.getAttribute('data-legend-key'), row.getAttribute(KEY_RECORD)]), [
    ['Genes', 'CDS'], ['tRNA', null], ['rRNA', null], ['Mine', null], ['Other', null]
  ]);
  assert.deepEqual(rowsOf(svg).map(pythonLegendKey), ['CDS', 'tRNA', 'rRNA', 'Mine', 'Other']);
  assert.deepEqual(editorRecords(svg), [null, null]);

  // A second pass of the same edits changes nothing and records nothing more.
  reconcileMountedResult(svg, full);
  assert.deepEqual(editorRecords(svg), [null, null]);

  // A reconcile without the adds removes the editor rows and leaves Python's rows as renamed and deleted.
  reconcileMountedResult(svg, legendOperations({ ...full, legendAdds: [] }), { domains: LEGEND_STRUCTURE });
  assert.deepEqual(rowsOf(svg).map(pythonLegendKey), ['CDS', 'tRNA', 'rRNA']);
  assert.deepEqual(rowsOf(svg).filter(isEditor), []);

  // Load: a saved Result with an editor row, read again by the Load plan and then reconciled with and without its add.
  const session = sessionFixture('forced-label-underlay-legend-rows.v44.gbdraw-session.json.gz');
  const [loaded] = loadSavedResults(session);
  assert.deepEqual(rowsOf(loaded).filter(isEditor).map((row) => row.getAttribute('data-legend-key')), ['Added']);
  assert.deepEqual(editorRecords(loaded), [null]);
  reconcileMountedResult(loaded, legendOperations({ legendAdds: [add('Added')] }));
  assert.deepEqual(editorRecords(loaded), [null]);
  reconcileMountedResult(loaded, legendOperations(), { domains: LEGEND_STRUCTURE });
  assert.deepEqual(rowsOf(loaded).filter(isEditor), []);
});

// An editor row whose key is the name a renamed row of Python's now shows: each
// row is found by its own key, once, and the reconcile without the add removes
// the editor row only.
test('an editor row keyed like a renamed row of Python\'s is indexed once and removed alone', () => {
  const svg = legendSvg();
  const operations = legendOperations({
    legendAdds: [add('Mine'), add('Last')],
    legendRenames: [{ from: 'CDS', to: 'Mine', allowMissing: false }],
    // Each key addresses its own row: Python's by the key it drew, the editor row's by its own.
    legendFills: [{ caption: 'CDS', color: '#112233' }, { caption: 'Mine', color: '#445566' }]
  });
  reconcileMountedResult(svg, operations);
  assert.deepEqual(rowsOf(svg).map((row) => [row.getAttribute('data-legend-key'), pythonLegendKey(row), isEditor(row)]), [
    ['Mine', 'CDS', false], ['tRNA', 'tRNA', false], ['rRNA', 'rRNA', false], ['Mine', 'Mine', true], ['Last', 'Last', true]
  ]);
  assert.deepEqual(swatchFills(svg), [
    ['CDS', false, '#112233'], ['tRNA', false, '#aaaaaa'], ['rRNA', false, '#aaaaaa'], ['Mine', true, '#445566'], ['Last', true, '#556677']
  ]);
  // Placement is the order of the adds, after Python's rows, and stable.
  reconcileMountedResult(svg, operations);
  assert.deepEqual(rowsOf(svg).map((row) => [pythonLegendKey(row), isEditor(row)]), [
    ['CDS', false], ['tRNA', false], ['rRNA', false], ['Mine', true], ['Last', true]
  ]);

  // Without its add the editor row goes and Python's renamed row stays; the other editor row is kept.
  reconcileMountedResult(svg, legendOperations({ ...operations, legendAdds: [add('Last')] }), { domains: LEGEND_STRUCTURE });
  assert.deepEqual(rowsOf(svg).map((row) => [row.getAttribute('data-legend-key'), pythonLegendKey(row), isEditor(row)]), [
    ['Mine', 'CDS', false], ['tRNA', 'tRNA', false], ['rRNA', 'rRNA', false], ['Last', 'Last', true]
  ]);
});

// The placement of a new editor row: before the next editor row an add after it
// names, found by that row's key.
test('a new editor row is placed before the next editor row added after it', () => {
  const svg = legendSvg();
  reconcileMountedResult(svg, legendOperations({ legendAdds: [add('B'), add('C')] }));
  assert.deepEqual(rowsOf(svg).map(pythonLegendKey), ['CDS', 'tRNA', 'rRNA', 'B', 'C']);
  reconcileMountedResult(svg, legendOperations({ legendAdds: [add('A'), add('B'), add('C')] }));
  assert.deepEqual(rowsOf(svg).map(pythonLegendKey), ['CDS', 'tRNA', 'rRNA', 'A', 'B', 'C']);
  assert.deepEqual(editorRecords(svg), [null, null, null]);
});

// Load names a color-rule row an older build saved (D-19) by the key it reads
// from the row; a row with a key record is a renamed row of Python's and is
// skipped, whatever its shown key is.
test('Load names a saved rule row by its key and skips a row that carries a key record', () => {
  const session = sessionFixture('rule-captions-with-color-id.v46.gbdraw-session.json.gz');
  const results = normalizeLogicalResults(session.results);
  const admission = admitFeatureCatalog(session.editorState.featureCatalog, results, { adopt: true, mode: 'circular' });
  const edits = savedResultEdits(session, 'circular', session.editorState, session.version, sessionSourceFactsOf(session.version, false), session.modes.circular.config.rules);
  const named = (patch) => {
    const svg = parseSvgText(results[0].content);
    patch(svg);
    createSavedResultPlan(results, admission, edits).operationsByResult[0].callerTransforms.forEach((transform) => transform(svg));
    return rowsOf(svg).filter((row) => row.getAttribute('data-legend-key').startsWith('Complex IV'))
      .map((row) => [row.getAttribute('data-legend-key'), row.getAttribute('data-legend-name')]);
  };
  assert.deepEqual(named(() => {}), [['Complex IV [#ffe700]', 'Complex IV'], ['Complex IV [#1f77b4]', 'Complex IV']]);
  assert.deepEqual(named((svg) => rowsOf(svg).find((row) => row.getAttribute('data-legend-key') === 'Complex IV [#ffe700]')
    .setAttribute(KEY_RECORD, 'Complex IV [#ffe700]')), [['Complex IV [#ffe700]', null], ['Complex IV [#1f77b4]', 'Complex IV']]);
});

// The executor reads an editor row's key where it reads every row's: through
// `pythonLegendKey`. No operation writes a record on an editor row, so this
// state is built by hand (the row's key record is `Old`, its attribute `New`):
// the index, the placement of a new row, and the reconcile that removes an
// editor row no add names each read the one key.
test('the executor reads an editor row\'s key through the one owner accessor', () => {
  const build = () => {
    const svg = legendSvg();
    const editor = legendRow('New', 79, { 'data-legend-owner': 'direct-editor', [KEY_RECORD]: 'Old' });
    svg.querySelector('#feature_legend').appendChild(editor);
    return svg;
  };
  assert.deepEqual(rowsOf(build()).map(pythonLegendKey), ['CDS', 'tRNA', 'rRNA', 'Old']);

  // The index finds the row by its key.
  const filled = build();
  reconcileMountedResult(filled, legendOperations({ legendFills: [{ caption: 'Old', color: '#123456' }] }), { domains: ['legendFills'] });
  assert.deepEqual(swatchFills(filled).at(-1), ['Old', true, '#123456']);

  // A new editor row goes before the next editor row an add after it names, by that row's key.
  const placed = build();
  reconcileMountedResult(placed, legendOperations({ legendAdds: [add('Z'), add('Old')] }), { domains: LEGEND_STRUCTURE });
  assert.deepEqual(rowsOf(placed).map(pythonLegendKey), ['CDS', 'tRNA', 'rRNA', 'Z', 'Old']);

  // The reconcile keeps the row an add names by its key and removes it otherwise.
  const kept = build();
  reconcileMountedResult(kept, legendOperations({ legendAdds: [add('Old')] }), { domains: LEGEND_STRUCTURE });
  assert.deepEqual(rowsOf(kept).map(pythonLegendKey), ['CDS', 'tRNA', 'rRNA', 'Old']);
  const removed = build();
  // An add of `New` is no add of `Old`: the row is removed and `New` is added as a row of its own.
  reconcileMountedResult(removed, legendOperations({ legendAdds: [add('New')] }), { domains: LEGEND_STRUCTURE });
  assert.deepEqual(rowsOf(removed).map(pythonLegendKey), ['CDS', 'tRNA', 'rRNA', 'New']);
  assert.deepEqual(editorRecords(removed), [null]);
});
