// A click on a drawn element of the Result opens the settings section that
// holds its settings: which element leads to which section, in every form its
// group takes; the pane's way to open a section by key; and the export, which
// keeps none of the preview's markers.
import assert from 'node:assert/strict';
import test from 'node:test';

import {
  DEPTH_GROUPS, DINUCLEOTIDE_CONTENT_GROUPS, DINUCLEOTIDE_SKEW_GROUPS, LENGTH_BAR_ID, PLOT_TITLE_ID, TICK_GROUPS,
  findTrackGroups, trackGroupSelector
} from '../../gbdraw/web/js/services/drawn-groups.js';
import {
  DRAWN_ELEMENT_SELECTOR, drawnElementSection, settingsSectionElement
} from '../../gbdraw/web/js/services/settings-sections.js';
import { stripTransientPreviewState } from '../../gbdraw/web/js/services/svg-serialization.js';

// One compound selector: an optional tag, an optional #id, an optional .class,
// and attribute tests (`[name]`, `[name="value"]`, `[name^="value"]`, `[name*="value"]`).
const matchesCompound = (element, compound) => {
  const [, tag, id, className, tests] = compound.match(/^([a-z]+)?(?:#([\w-]+))?(?:\.([\w-]+))?((?:\[[^\]]+\])*)$/) ?? [];
  if (!tag && !id && !className && !tests) throw new Error(`unsupported selector: ${compound}`);
  if (tag && element.tag !== tag) return false;
  if (id && element.getAttribute('id') !== id) return false;
  if (className && !String(element.getAttribute('class') ?? '').split(/\s+/).includes(className)) return false;
  return Array.from((tests || '').matchAll(/\[([\w-]+)(?:([\^*])?="([^"]*)")?\]/g)).every(([, name, operator, value]) => {
    const actual = element.getAttribute(name);
    if (actual === null) return false;
    if (value === undefined) return true;
    return operator === '^' ? actual.startsWith(value) : operator === '*' ? actual.includes(value) : actual === value;
  });
};

class FakeElement {
  constructor(tag, attributes = {}, children = []) {
    this.tag = tag;
    this.attributes = { ...attributes };
    this.parent = null;
    this.children = children;
    children.forEach((child) => { child.parent = this; });
  }

  getAttribute(name) { return this.attributes[name] ?? null; }

  setAttribute(name, value) { this.attributes[name] = String(value); }

  removeAttribute(name) { delete this.attributes[name]; }

  matches(selector) { return selector.split(',').some((part) => matchesCompound(this, part.trim())); }

  closest(selector) {
    for (let node = this; node; node = node.parent) if (node.matches(selector)) return node;
    return null;
  }

  contains(other) {
    for (let node = other; node; node = node.parent) if (node === this) return true;
    return false;
  }

  querySelectorAll(selector) {
    const found = [];
    const visit = (node) => node.children.forEach((child) => {
      if (child.matches(selector)) found.push(child);
      visit(child);
    });
    visit(this);
    return found;
  }

  querySelector(selector) { return this.querySelectorAll(selector)[0] ?? null; }
}

const g = (attributes, ...children) => new FakeElement('g', attributes, children);
const ink = (tag = 'path') => new FakeElement(tag);
const svgOf = (...children) => new FakeElement('svg', {}, children);
const slot = (id, renderer) => ({ id, 'data-gbdraw-slot-id': id, 'data-gbdraw-slot-renderer': renderer });

const withElement = (run) => {
  const previous = globalThis.Element;
  globalThis.Element = FakeElement;
  try { return run(); } finally { globalThis.Element = previous; }
};

// The section a click on `target` opens, with the group found.
const sectionOf = (target) => withElement(() => drawnElementSection(target));

test('the shared owner of the group patterns gives each kind one selector', () => {
  assert.equal(trackGroupSelector(TICK_GROUPS), 'g[data-gbdraw-slot-renderer="ticks"], g[id="tick"], g[id^="tick_"]');
  assert.equal(trackGroupSelector(DINUCLEOTIDE_CONTENT_GROUPS),
    'g[data-gbdraw-slot-renderer="dinucleotide_content"], g[id="gc_content"], g[id^="gc_content_"]');
  assert.equal(trackGroupSelector(DINUCLEOTIDE_SKEW_GROUPS),
    'g[data-gbdraw-slot-renderer="dinucleotide_skew"], g[id="skew"], g[id^="skew_"], g[id="gc_skew"], g[id^="gc_skew_"]');
  assert.deepEqual([LENGTH_BAR_ID, PLOT_TITLE_ID], ['length_bar', 'plot_title']);
  for (const selector of [TICK_GROUPS, DINUCLEOTIDE_CONTENT_GROUPS, DINUCLEOTIDE_SKEW_GROUPS].map(trackGroupSelector)
    .concat(`#${LENGTH_BAR_ID}`, `#${PLOT_TITLE_ID}`)) {
    assert.ok(DRAWN_ELEMENT_SELECTOR.includes(selector), selector);
  }
});

test('the groups of a kind are the marked ones, or without marks those found by id', () => {
  const marked = g(slot('track_slot_8bf4cb58aa6767fe_gc_skew', 'dinucleotide_skew'));
  const markedDefault = g({ id: 'skew', 'data-gbdraw-slot-renderer': 'dinucleotide_skew' });
  const child = g({ id: 'skew_axis' });
  const depth = g(slot('track_slot_depth', 'depth'));
  assert.deepEqual(findTrackGroups(svgOf(marked, markedDefault, child, depth), DINUCLEOTIDE_SKEW_GROUPS), [marked, markedDefault],
    'the marked groups leave out a group found only by id');
  assert.deepEqual(findTrackGroups(svgOf(marked, depth), DEPTH_GROUPS), [depth]);
  const stored = g({ id: 'gc_skew_record_2' });
  const storedCircular = g({ id: 'skew' });
  assert.deepEqual(findTrackGroups(svgOf(stored, storedCircular, g({ id: 'gc_content' })), DINUCLEOTIDE_SKEW_GROUPS),
    [storedCircular, stored]);
  assert.deepEqual(findTrackGroups(null, TICK_GROUPS), []);
});

test('Circular, the default layout: the ticks open Axis & Scale, a dinucleotide track its section', () => {
  const label = ink('textPath');
  const content = ink();
  const axis = ink('line');
  const skew = ink();
  const feature = ink();
  const svg = svgOf(
    g({ id: 'record' }, feature),
    g(slot('tick', 'ticks'), g({}, new FakeElement('text', {}, [label]))),
    g({ ...slot('gc_content', 'dinucleotide_content'), id: 'gc_content' }, content, g({ id: 'gc_content_axis' }, axis)),
    g({ ...slot('gc_skew', 'dinucleotide_skew'), id: 'skew' }, skew)
  );
  assert.equal(sectionOf(label), 'circular-axis');
  assert.equal(sectionOf(content), 'dinucleotide');
  assert.equal(sectionOf(axis), 'dinucleotide', 'a child group resolves to its track');
  assert.equal(sectionOf(skew), 'dinucleotide');
  assert.equal(sectionOf(feature), null);
  assert.equal(sectionOf(svg), null);
});

test('Circular, custom track slots and the grid: the slot renderer and the record ids', () => {
  const tick = ink('text');
  const atSkew = ink();
  const childInk = ink();
  const features = ink();
  const depth = ink();
  const custom = svgOf(
    g(slot('track_slot_099f388a5a722768_ticks', 'ticks'), tick),
    g(slot('track_slot_c235170817675430_a_skew_2', 'dinucleotide_skew'), atSkew),
    g(slot('track_slot_a9a17e0b3d691d4d_gc_content', 'dinucleotide_content'), g({ id: 'track_slot_child_0f0f0f0f0f0f0f0f' }, childInk)),
    g(slot('track_slot_3bf4b9268ef1dc03_features', 'features'), features),
    g(slot('track_slot_11aa22bb33cc44dd_depth', 'depth'), depth)
  );
  assert.equal(sectionOf(tick), 'circular-axis');
  assert.equal(sectionOf(atSkew), 'dinucleotide', 'another dinucleotide track');
  assert.equal(sectionOf(childInk), 'dinucleotide');
  assert.equal(sectionOf(features), null);
  assert.equal(sectionOf(depth), null);

  // A Result stored without slot marks is read by id.
  const gridTick = ink('text');
  const gridContent = ink();
  const gridSkew = ink();
  const grid = svgOf(g({ id: 'record_1' }, g({ id: 'tick_1' }, gridTick), g({ id: 'gc_content_1' }, gridContent), g({ id: 'skew_1' }, gridSkew)));
  assert.equal(sectionOf(gridTick), 'circular-axis');
  assert.equal(sectionOf(gridContent), 'dinucleotide');
  assert.equal(sectionOf(gridSkew), 'dinucleotide');
});

test('a record definition has no section, whatever its record is called', () => {
  // The definition group is named `<record id>_definition`: these ids start like a track's.
  for (const record of ['tick_chr1', 'skew_control', 'gc_content_demo', 'gc_skew_1']) {
    const text = ink('text');
    const row = g({ 'data-gbdraw-role': 'record-definition-row' }, text);
    const definition = g({ id: `${record}_definition`, 'data-gbdraw-role': 'record-definition', 'data-gbdraw-record-index': '0' }, row);
    svgOf(definition);
    assert.equal(sectionOf(text), null, record);
    assert.equal(sectionOf(definition), null, record);
    assert.equal(withElement(() => definition.matches(DRAWN_ELEMENT_SELECTOR)), true, 'the id alone reads like a track');
  }
});

test('a Result drawn by 0.13.0: the groups of a record have no section, whatever the record is called', () => {
  // That release named the group of a record's features `<record id>` and its definition
  // `<record id>_definition` (each with `_record_<n>` among several records), without a role.
  for (const record of ['tick_chr1', 'skew_control', 'gc_content_x', 'gc_skew_1']) {
    for (const suffix of ['', '_record_2']) {
      const axis = ink('line');
      const text = ink('text');
      const features = g({ id: `${record}${suffix}` }, axis, new FakeElement('path', { 'data-gbdraw-feature-id': 'f7cb30df9' }));
      const definition = g({ id: `${record}_definition${suffix}` }, text);
      svgOf(features, definition);
      for (const target of [axis, features, text, definition]) assert.equal(sectionOf(target), null, `${target.getAttribute('id') ?? target.tag} of ${record}${suffix}`);
      assert.equal(withElement(() => features.matches(DRAWN_ELEMENT_SELECTOR) && definition.matches(DRAWN_ELEMENT_SELECTOR)), true,
        'the ids alone read like a track');
    }
  }
  // Its tracks have the ids alone too, and open their section.
  const skew = ink();
  svgOf(g({ id: 'gc_skew' }, skew));
  assert.equal(sectionOf(skew), 'dinucleotide');
});

test('Linear: the scale opens Axis & Scale, the tracks of every record their section, the title Titles', () => {
  const bar = ink('line');
  const barText = ink('text');
  const title = ink('text');
  const content = ink();
  const skew = ink();
  const second = ink();
  const secondSkew = ink();
  const legendRow = ink();
  const comparison = ink('rect');
  const svg = svgOf(
    g({ id: 'length_bar', 'data-gbdraw-composition-role': 'primary' }, bar, barText),
    g({ id: 'plot_title', 'data-gbdraw-composition-role': 'title' }, title),
    g({ ...slot('gc_content', 'dinucleotide_content'), id: 'gc_content' }, content),
    g({ ...slot('gc_skew', 'dinucleotide_skew'), id: 'gc_skew' }, skew),
    g({ id: 'gc_content_record_2' }, second),
    g({ id: 'gc_skew_record_2' }, secondSkew),
    g({ id: 'legend' }, g({ 'data-legend-key': 'GC skew' }, legendRow), g({ 'data-gbdraw-role': 'comparison-legend' }, comparison))
  );
  assert.equal(sectionOf(bar), 'linear-axis');
  assert.equal(sectionOf(barText), 'linear-axis');
  assert.equal(sectionOf(title), 'titles');
  assert.equal(sectionOf(content), 'dinucleotide');
  assert.equal(sectionOf(skew), 'dinucleotide');
  assert.equal(sectionOf(second), 'dinucleotide');
  assert.equal(sectionOf(secondSkew), 'dinucleotide');
  assert.equal(sectionOf(legendRow), null, 'a Legend row is no drawn element with a section');
  assert.equal(sectionOf(comparison), null, 'the comparison legend');
  // An element of another document, and a target that is no element.
  assert.equal(sectionOf({}), null);
  assert.equal(sectionOf(null), null);
});

test('an export keeps no marker of a drawn element that opens a section', () => {
  const track = g({ ...slot('gc_content', 'dinucleotide_content'), 'data-gbdraw-preview-item': 'settings' });
  const row = g({ 'data-legend-key': 'tRNA', 'data-gbdraw-preview-item': '', role: 'button', tabindex: '0', 'aria-haspopup': 'dialog' });
  stripTransientPreviewState(svgOf(track, row));
  assert.deepEqual(track.attributes, slot('gc_content', 'dinucleotide_content'));
  assert.deepEqual(row.attributes, { 'data-legend-key': 'tRNA' });
});

// The pane's controller over a document that holds one closed section per key
// of the shown mode, 40 px below the top edge of the settings.
class FakeDetails {
  open = false;

  parentElement = null;
}
const createPane = async (mode, display = 'stacked') => {
  const calls = [];
  const sections = new Map();
  const sectionElement = (key) => {
    if (sections.has(key)) return sections.get(key);
    const details = new FakeDetails();
    const heading = {
      parentElement: details,
      focus: (options) => calls.push(['focus', key, options]),
      closest: (selector) => (selector === 'details' ? details : null)
    };
    const section = {
      dataset: { settingsSection: key },
      details,
      closest: (selector) => (selector === '[data-settings-section]' ? section : area),
      querySelector: () => heading,
      getBoundingClientRect: () => ({ top: 140, bottom: 400, height: 260 })
    };
    sections.set(key, section);
    return section;
  };
  const shown = { circular: ['titles', 'circular-axis', 'dinucleotide'], linear: ['titles', 'linear-axis', 'dinucleotide'] };
  const area = {
    contains: (node) => node instanceof FakeDetails,
    querySelector: (selector) => {
      const key = selector.match(/^\[data-settings-section="([^"]+)"\]$/)?.[1];
      return key && shown[mode].includes(key) ? sectionElement(key) : null;
    },
    getBoundingClientRect: () => ({ top: 100, bottom: 900 }),
    scrollBy: (options) => calls.push(['scroll', options.top])
  };
  globalThis.HTMLDetailsElement = FakeDetails;
  globalThis.getComputedStyle = () => ({ overflowY: 'auto' });
  globalThis.document = { querySelector: (selector) => (selector === '.settings-scroll' ? area : null), scrollingElement: null };
  globalThis.window = {
    localStorage: null, addEventListener() {}, removeEventListener() {},
    Vue: {
      computed: (read) => ({ get value() { return read(); } }), nextTick: async () => {},
      onMounted() {}, onUnmounted() {}, reactive: (value) => value, ref: (value) => ({ value }), watch() {}
    }
  };
  const state = { mode: { value: mode }, settingsPaneDisplay: { value: display }, settingsPaneDocked: { value: display === 'rail' ? 'layout' : null } };
  const { createSettingsPaneController } = await import('../../gbdraw/web/js/app/settings-pane.js');
  return { pane: createSettingsPaneController({ state }), state, calls, sections };
};

test('the pane opens the section of a key in the shown mode, at the top, with the focus on its heading', async () => {
  const { pane, calls, sections } = await createPane('circular');
  await pane.revealSettingsSection('dinucleotide');
  assert.deepEqual(calls, [['scroll', 40], ['focus', 'dinucleotide', { preventScroll: true }]]);
  assert.equal(sections.get('dinucleotide').details.open, true);
  assert.deepEqual(pane.settingsReveal.last, { key: 'dinucleotide', top: 140 });
});

test('a key of the other mode opens the section of the same title in the shown mode', async () => {
  // The Linear scale's section, asked for while Circular shows, and the reverse.
  const circular = await createPane('circular');
  await circular.pane.revealSettingsSection('linear-axis');
  assert.deepEqual(circular.calls.at(-1), ['focus', 'circular-axis', { preventScroll: true }]);
  assert.equal(circular.pane.settingsReveal.last.key, 'circular-axis');
  const linear = await createPane('linear', 'rail');
  await linear.pane.revealSettingsSection('circular-axis');
  assert.deepEqual(linear.calls.at(-1), ['focus', 'linear-axis', { preventScroll: true }]);
  assert.equal(linear.state.settingsPaneDocked.value, 'linear-axis', 'the rail docks the section');
  // A key without a section in the shown mode opens nothing.
  const before = linear.calls.length;
  await linear.pane.revealSettingsSection('pairwise-comparisons');
  assert.equal(linear.calls.length, before);
  assert.equal(linear.pane.settingsReveal.last.key, 'linear-axis');
});

test('the catalog finds the element of a section in the settings', () => {
  const found = {};
  assert.equal(settingsSectionElement({ querySelector: (selector) => (selector === '[data-settings-section="titles"]' ? found : null) }, 'titles'), found);
});
