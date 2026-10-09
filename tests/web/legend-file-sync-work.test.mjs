// S4 (perf 0.14.x, counts not timings), on the edit port (U3a A2a): the
// Legend rows the specific-color rules draw are computed on the listed rows, so
// their preparation copies no node of the mounted Result, and their `apply`
// writes the intent only: it neither lays out nor serializes the Result. The
// rule commit's one show (`projectPaletteAndRules` -> `showEditorIntent` ->
// preview-runtime `applyEditorOperations`) lays out a changed Legend and
// serializes the Result once (live-generate-parity-legend.playwright.spec.js,
// "a rule commit ... serializes the Result once").
import assert from 'node:assert/strict';
import { test } from 'node:test';
import { createLegendEntryActions } from '../../gbdraw/web/js/app/legend/entry-actions.js';
import { SPECIFIC_COLOR_FILE_OWNER } from '../../gbdraw/web/js/services/legend-svg.js';
import { withDrawings } from './helpers/drawing-state.mjs';

const ref = (value) => ({ value });
let cloned = 0;

// The DOM reads and writes of the Legend owner, on plain objects.
class Node {
  constructor(tagName, attributes = {}, children = [], text = '') {
    this.tagName = tagName;
    this.attributes = new Map(Object.entries(attributes));
    this.children = [];
    this.parentElement = null;
    this.ownText = text;
    this.style = {};
    children.forEach((child) => this.appendChild(child));
  }

  get id() { return this.getAttribute('id') || ''; }
  get textContent() { return this.ownText + this.children.map((child) => child.textContent).join(''); }
  set textContent(value) { this.children = []; this.ownText = String(value); }
  getAttribute(name) { return this.attributes.has(name) ? this.attributes.get(name) : null; }
  hasAttribute(name) { return this.attributes.has(name); }
  setAttribute(name, value) { this.attributes.set(name, String(value)); }
  removeAttribute(name) { this.attributes.delete(name); }
  appendChild(child) {
    child.remove();
    child.parentElement = this;
    this.children.push(child);
    return child;
  }

  remove() {
    if (!this.parentElement) return;
    this.parentElement.children = this.parentElement.children.filter((child) => child !== this);
    this.parentElement = null;
  }

  replaceWith(other) {
    const parent = this.parentElement;
    other.remove();
    parent.children[parent.children.indexOf(this)] = other;
    other.parentElement = parent;
    this.parentElement = null;
  }

  cloneNode(deep = false) {
    cloned += 1;
    const copy = new Node(this.tagName, Object.fromEntries(this.attributes), [], this.ownText);
    if (deep) this.children.forEach((child) => copy.appendChild(child.cloneNode(true)));
    return copy;
  }

  isEqualNode(other) {
    return this.tagName === other.tagName && this.ownText === other.ownText
      && this.attributes.size === other.attributes.size
      && [...this.attributes].every(([name, value]) => other.getAttribute(name) === value)
      && this.children.length === other.children.length
      && this.children.every((child, index) => child.isEqualNode(other.children[index]));
  }

  matches(selector) {
    const [, tag, id, attribute] = selector.match(/^([a-z]*)(?:#([\w-]+))?(?:\[([\w-]+)\])?$/) || [];
    return (!tag || this.tagName === tag) && (!id || this.id === id) && (!attribute || this.hasAttribute(attribute));
  }

  querySelectorAll(selector) {
    return this.children.flatMap((child) => [...(child.matches(selector) ? [child] : []), ...child.querySelectorAll(selector)]);
  }

  querySelector(selector) { return this.querySelectorAll(selector)[0] || null; }
  getElementById(id) { return this.querySelector(`#${id}`); }
}

const row = (caption, color, owner = null, y = 0) => new Node('g', {
  'data-legend-key': caption, transform: `translate(0,${y})`, ...(owner ? { 'data-legend-owner': owner } : {})
}, [new Node('path', { fill: color }), new Node('text', {}, [], caption)]);

const FEATURES = 200;
const setup = () => {
  const legend = new Node('g', { id: 'legend' }, [new Node('g', { id: 'feature_legend' }, [
    row('CDS', '#111111', null, 0), row('Rule A', '#222222', SPECIFIC_COLOR_FILE_OWNER, 20)
  ])]);
  const features = Array.from({ length: FEATURES }, (_, index) => new Node('path', { id: `f${index}`, fill: '#111111' }));
  const svg = new Node('svg', { 'data-gbdraw-composition-schema': '1' }, [new Node('g', { id: 'features' }, features), legend]);
  const state = withDrawings({
    results: ref([{ content: '<svg/>' }]), svgContainer: ref({ querySelector: (selector) => (selector === 'svg' ? svg : null) }),
    originalLegendOrder: ref(['CDS', 'Rule A']), originalLegendColors: ref({ CDS: '#111111' }),
    newLegendCaption: ref(''), newLegendColor: ref('#000000'), originalSvgStroke: ref({ color: null, width: null }),
    legendEntries: ref([
      { caption: 'CDS', originalCaption: 'CDS', color: '#111111', featureIds: [] },
      { caption: 'Rule A', originalCaption: 'Rule A', color: '#222222', featureIds: [] }
    ]), dormantLegendEntries: ref([]), deletedLegendEntries: ref([]),
    legendColorOverrides: {}, legendStrokeOverrides: {}
  });
  const work = { commits: 0, layouts: 0 };
  const commitActiveResultEdit = () => { work.commits += 1; return true; };
  const actions = createLegendEntryActions({ state, commitActiveResultEdit });
  // The layout owner (legend-layout/reposition-actions.js `refreshLegendGeometry`)
  // commits what it laid out unless asked not to.
  actions.setLegendGeometryChangedHandler(({ commit = true } = {}) => {
    work.layouts += 1;
    if (commit) commitActiveResultEdit();
  });
  return { svg, state, actions, work };
};
const intent = (caption, color) => ({ caption, color });

const listed = (state) => state.activeDrawing().legendEntries.value.map(({ caption, color }) => [caption, color]);
const fills = (svg) => svg.getElementById('legend').querySelectorAll('g[data-legend-key]')
  .map((entry) => [entry.getAttribute('data-legend-key'), entry.querySelector('path').getAttribute('fill')]);

test('a rule commit that leaves the Legend as it is copies no node and serializes nothing', async () => {
  const { svg, state, actions, work } = setup();
  cloned = 0;
  const legend = await actions.prepareFileLegendEntries([intent('Rule A', '#222222')], {
    previousFileIntents: [intent('Rule A', '#222222')]
  });
  assert.deepEqual([legend.diff.add, legend.diff.update, legend.diff.remove], [[], [], []]);
  assert.equal(cloned, 0, `${cloned} nodes copied`);
  const before = svg.getElementById('legend');
  assert.equal(legend.apply(), null, 'no Legend row to show');
  assert.deepEqual(work, { commits: 0, layouts: 0 });
  assert.equal(svg.getElementById('legend'), before);
  assert.deepEqual(listed(state), [['CDS', '#111111'], ['Rule A', '#222222']]);
});

test('a rule commit that changes a Legend row writes its intent and leaves the show to the one compile', async () => {
  const { svg, state, actions, work } = setup();
  cloned = 0;
  const legend = await actions.prepareFileLegendEntries([intent('Rule A', '#333333')], {
    previousFileIntents: [intent('Rule A', '#222222')]
  });
  assert.equal(legend.diff.update.length, 1);
  assert.equal(cloned, 0, `${cloned} nodes copied`);
  // A color change is a `legendFills` operation of the rule commit's one show.
  assert.equal(legend.apply(), null);
  assert.deepEqual(work, { commits: 0, layouts: 0 });
  assert.deepEqual(listed(state), [['CDS', '#111111'], ['Rule A', '#333333']]);
  assert.deepEqual(fills(svg), [['CDS', '#111111'], ['Rule A', '#222222']], 'apply writes no node');
  assert.equal(svg.querySelectorAll('path').length, FEATURES + 2);
});

test('a rule commit that adds a Legend row returns it for the one show and serializes nothing itself', async () => {
  const { svg, state, actions, work } = setup();
  cloned = 0;
  const legend = await actions.prepareFileLegendEntries([intent('Rule A', '#222222'), intent('Rule B', '#444444')], {
    previousFileIntents: [intent('Rule A', '#222222')]
  });
  assert.equal(cloned, 0, `${cloned} nodes copied`);
  assert.deepEqual(legend.apply(), { add: [{ caption: 'Rule B', color: '#444444' }], retire: [] });
  assert.deepEqual(work, { commits: 0, layouts: 0 });
  assert.deepEqual(listed(state), [['CDS', '#111111'], ['Rule A', '#222222'], ['Rule B', '#444444']]);
  assert.deepEqual(fills(svg), [['CDS', '#111111'], ['Rule A', '#222222']], 'apply writes no node');
});
