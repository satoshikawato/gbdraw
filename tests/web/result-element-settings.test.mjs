// A click on a drawn element of the Result opens the settings section that
// holds its settings in a popup beside the Result: which element leads to
// which section, in every form its group takes; where the pane's owner shows
// the section (the popup, or the dock that already shows it) and how the pane
// takes it back; and the export, which keeps none of the preview's markers.
import assert from 'node:assert/strict';
import { readdirSync, readFileSync } from 'node:fs';
import { join } from 'node:path';
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
// of the shown mode, 40 px below the top edge of the settings. The popup host
// finds the same sections: a section is one element, in the pane or the popup.
class FakeDetails {
  open = false;

  parentElement = null;
}
// The module reads Vue once: its watchers and the window's listeners are those of the pane made last.
let watchers = [];
const Vue = {
  computed: (read) => ({ get value() { return read(); } }), nextTick: async (run) => run?.(),
  onMounted: (run) => run(), onUnmounted() {}, ref: (value) => ({ value }), watch: (source, run) => watchers.push({ source, run })
};
const createPane = async (mode, display = 'stacked') => {
  const calls = [];
  const listeners = new Map();
  watchers = [];
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
      scrollHeight: 500,
      clientHeight: 250,
      closest: (selector) => (selector === '[data-settings-section]' ? section : area),
      querySelector: () => heading,
      getBoundingClientRect: () => ({ top: 140, bottom: 400, height: 260 })
    };
    sections.set(key, section);
    return section;
  };
  const shown = { circular: ['titles', 'circular-axis', 'dinucleotide'], linear: ['titles', 'linear-axis', 'dinucleotide'] };
  const sectionOfSelector = (selector) => {
    const key = selector.match(/^\[data-settings-section="([^"]+)"\]$/)?.[1];
    return key && shown[mode].includes(key) ? sectionElement(key) : null;
  };
  const area = {
    contains: (node) => node instanceof FakeDetails,
    querySelector: sectionOfSelector,
    getBoundingClientRect: () => ({ top: 100, bottom: 900 }),
    scrollBy: (options) => calls.push(['scroll', options.top])
  };
  const popup = { offsetWidth: 320, offsetHeight: 300, contains: () => false, querySelector: sectionOfSelector };
  globalThis.HTMLDetailsElement = FakeDetails;
  globalThis.getComputedStyle = () => ({ overflowY: 'auto' });
  globalThis.document = {
    querySelector: (selector) => (selector === '.settings-scroll' ? area : null),
    getElementById: (id) => (id === 'settings-popup' ? popup : null),
    scrollingElement: null, activeElement: null
  };
  globalThis.window = {
    localStorage: null, innerWidth: 1600, innerHeight: 1000, Vue,
    addEventListener: (type, listener) => listeners.set(type, listener), removeEventListener() {}
  };
  const state = {
    mode: { value: mode }, settingsPaneDisplay: { value: display }, settingsPaneDocked: { value: display === 'rail' ? 'layout' : null },
    settingsPopup: { value: null }, svgContent: { value: '<svg/>' }, resultPanelTab: { value: 'preview' }
  };
  const { createSettingsPaneController } = await import('../../gbdraw/web/js/app/settings-pane.js');
  const pane = createSettingsPaneController({ state });
  // The window at a new size, as the browser tells the page.
  const resize = (width, height) => {
    Object.assign(globalThis.window, { innerWidth: width, innerHeight: height });
    listeners.get('resize')();
  };
  return { pane, state, calls, sections, resize, watchers };
};
// The Result's viewport, and a click on an item 20 px wide at (700, 500).
const RESULT = { top: 100, bottom: 900, left: 400, right: 1500 };
const CLICK = { x: 700, y: 500, avoid: { left: 690, right: 710, top: 490, bottom: 510 } };

test('the pane opens the section of a key in the shown mode, at the top, with the focus on its heading', async () => {
  const { pane, calls, sections } = await createPane('circular');
  await pane.revealSettingsSection('dinucleotide');
  assert.deepEqual(calls, [['scroll', 40], ['focus', 'dinucleotide', { preventScroll: true }]]);
  assert.equal(sections.get('dinucleotide').details.open, true);
  // A key of the other mode opens the section of the same title in the shown mode.
  await pane.revealSettingsSection('linear-axis');
  assert.deepEqual(calls.at(-1), ['focus', 'circular-axis', { preventScroll: true }]);
  const linear = await createPane('linear', 'rail');
  await linear.pane.revealSettingsSection('circular-axis');
  assert.deepEqual(linear.calls.at(-1), ['focus', 'linear-axis', { preventScroll: true }]);
  assert.equal(linear.state.settingsPaneDocked.value, 'linear-axis', 'the rail docks the section');
  // A key without a section in the shown mode opens nothing.
  const before = linear.calls.length;
  await linear.pane.revealSettingsSection('pairwise-comparisons');
  assert.equal(linear.calls.length, before);
});

test('a clicked element\'s section is shown in the popup beside the clicked point, inside the Result', async () => {
  const { pane, state, calls } = await createPane('circular');
  assert.equal(pane.settingsPane.value.floating, null);
  assert.equal(pane.showSettingsBesideResult('dinucleotide', CLICK, RESULT), true);
  // Right of the clicked item, as high as the section's content (300 + 500 - 250) allows.
  assert.deepEqual(state.settingsPopup.value, { section: 'dinucleotide', x: 722, y: 225, maxHeight: 663 });
  assert.equal(pane.settingsPane.value.floating, 'dinucleotide');
  assert.deepEqual(calls, [['focus', 'dinucleotide', { preventScroll: true }]], 'the pane does not scroll');
  // Another element: the popup shows its section, of the shown mode, in place of the first.
  assert.equal(pane.showSettingsBesideResult('linear-axis', { ...CLICK, x: 1400, avoid: { left: 1390, right: 1410, top: 490, bottom: 510 } }, RESULT), true);
  assert.deepEqual(state.settingsPopup.value, { section: 'circular-axis', x: 1058, y: 225, maxHeight: 663 }, 'left of an item near the right edge');
  // An item as wide as the Result (a ring) leaves no room beside it: beside the point.
  pane.showSettingsBesideResult('dinucleotide', { ...CLICK, avoid: { left: 420, right: 1480, top: 150, bottom: 850 } }, RESULT);
  assert.equal(state.settingsPopup.value.x, 712);
  assert.equal(pane.closeSettingsPopup(), true);
  assert.equal(pane.closeSettingsPopup(), false, 'nothing was open');
  assert.equal(state.settingsPopup.value, null);
});

// M-1: "over" is an overlap on both axes. A ring is as wide and as high as the Result; a title or a
// Linear track is wide and flat, and the popup has room below or above it.
test('the popup stands below or above a wide, flat item, and beside the point when neither side holds it', async () => {
  const { pane, state } = await createPane('linear');
  const title = { left: 420, right: 1480, top: 120, bottom: 150 };
  pane.showSettingsBesideResult('titles', { x: 950, y: 135, avoid: title }, RESULT);
  assert.deepEqual(state.settingsPopup.value, { section: 'titles', x: 790, y: 162, maxHeight: 726 }, 'below the title, clear of it');
  pane.showSettingsBesideResult('titles', { x: 950, y: 865, avoid: { ...title, top: 850, bottom: 880 } }, RESULT);
  assert.deepEqual(state.settingsPopup.value, { section: 'titles', x: 790, y: 288, maxHeight: 550 }, 'above it, ending 12 px over it');
  // A track in the middle of the Result: neither side holds two thirds of it.
  pane.showSettingsBesideResult('dinucleotide', { x: 950, y: 495, avoid: { ...title, top: 480, bottom: 510 } }, RESULT);
  assert.deepEqual([state.settingsPopup.value.x, state.settingsPopup.value.y], [962, 220], 'beside the point');
});

// S-1: the popup belongs to the Result that is shown.
test('the popup closes when the Result panel shows another tab, and when no Result is shown', async () => {
  for (const hide of [(state) => { state.resultPanelTab.value = 'run-info'; }, (state) => { state.svgContent.value = null; }]) {
    const pane = await createPane('circular');
    const shown = pane.watchers.find(({ source }) => source() === true);
    pane.pane.showSettingsBesideResult('titles', CLICK, RESULT);
    shown.run(shown.source());
    assert.equal(pane.state.settingsPopup.value.section, 'titles', 'a Result that stays shown keeps it');
    hide(pane.state);
    shown.run(shown.source());
    assert.equal(pane.state.settingsPopup.value, null);
  }
});

// S-2: the popup is 320 x 300 here; a drag keeps it 12 px inside the window.
test('a window that changed size moves the popup the least distance back inside it, and does not close it', async () => {
  const { pane, state, resize } = await createPane('circular');
  resize(900, 500);
  assert.equal(state.settingsPopup.value, null, 'no popup, nothing to move');
  pane.showSettingsBesideResult('dinucleotide', CLICK, RESULT);
  const placed = { section: 'dinucleotide', x: 722, y: 225, maxHeight: 663 };
  resize(1600, 1000);
  assert.deepEqual(state.settingsPopup.value, placed, 'inside the window: where it was, as high as it was');
  resize(900, 500);
  assert.deepEqual(state.settingsPopup.value, { section: 'dinucleotide', x: 568, y: 188, maxHeight: 0 });
  resize(1600, 1000);
  assert.deepEqual(state.settingsPopup.value, { section: 'dinucleotide', x: 568, y: 188, maxHeight: 0 }, 'a larger window moves nothing back');
  // A drag moves it through the same transition.
  pane.moveSettingsPopup({ x: 40, y: 60, maxHeight: 0 });
  assert.deepEqual(state.settingsPopup.value, { section: 'dinucleotide', x: 40, y: 60, maxHeight: 0 });
});

// The popup's place has one writer: the transitions of services/settings-pane-state.js.
test('no module but the pane-state service writes the popup', () => {
  const root = join(process.cwd(), 'gbdraw', 'web', 'js');
  const modules = readdirSync(root, { recursive: true }).filter((file) => file.endsWith('.js') && !file.startsWith('vendor'));
  const writers = modules.filter((file) => /settingsPopup\.value(?:\.\w+)?\s*=[^=]|Object\.assign\(\s*(?:state\.)?settingsPopup/.test(readFileSync(join(root, file), 'utf8')));
  assert.deepEqual(writers.map((file) => file.replaceAll('\\', '/')), ['services/settings-pane-state.js']);
  // A popup of the Result is dragged through a port that moves it; the drag holds no placement to write.
  const setup = readFileSync(join(root, 'app', 'app-setup.js'), 'utf8');
  assert.match(setup, /moveTo: settingsPaneActions\.moveSettingsPopup/);
  // Only the four sections a drawn element opens can be shown beside the Result: the port takes no other key.
  assert.match(readFileSync(join(root, 'app', 'settings-pane.js'), 'utf8'), /@typedef \{\(key: DrawnElementSectionKey, anchor: PopupAnchor, bounds: PopupBounds\) => boolean\} ShowSettingsBesideResult/);
});

test('the docked section takes the focus where it is, and a pane action takes a section back from the popup', async () => {
  const rail = await createPane('linear', 'rail');
  rail.pane.showSettingsPaneRail('linear-axis');
  assert.equal(rail.pane.showSettingsBesideResult('circular-axis', CLICK, RESULT), true);
  await new Promise((resolve) => { setImmediate(resolve); });
  assert.equal(rail.state.settingsPopup.value, null, 'no popup for the section the dock shows');
  assert.deepEqual(rail.calls.at(-1), ['focus', 'linear-axis', { preventScroll: true }]);
  // Another section floats beside the rail; opening it from the index docks it and closes the popup.
  rail.pane.showSettingsBesideResult('titles', CLICK, RESULT);
  assert.equal(rail.state.settingsPopup.value.section, 'titles');
  await rail.pane.revealSettingsSection('titles');
  assert.deepEqual([rail.state.settingsPopup.value, rail.state.settingsPaneDocked.value], [null, 'titles']);

  const stacked = await createPane('circular');
  stacked.pane.showSettingsBesideResult('titles', CLICK, RESULT);
  // A setting inside the floating section is brought into view in the pane.
  const titles = stacked.sections.get('titles');
  await stacked.pane.revealSetting(titles, { focus: titles.querySelector() });
  assert.equal(stacked.state.settingsPopup.value, null);
  assert.deepEqual(stacked.calls.at(-1), ['focus', 'titles', { preventScroll: true }]);
  // Opening another section in the pane leaves the popup open.
  stacked.pane.showSettingsBesideResult('titles', CLICK, RESULT);
  await stacked.pane.revealSettingsSection('dinucleotide');
  assert.equal(stacked.state.settingsPopup.value.section, 'titles');
});

// One section, one place: the popup shows the section's own element, moved
// there. A copy would give two controls one id and one History target.
test('each section a drawn element opens moves into the popup host of the app root', () => {
  const html = readFileSync(join(process.cwd(), 'gbdraw', 'web', 'index.html'), 'utf8');
  const moved = Array.from(html.matchAll(
    /<teleport (?:v-if="[^"]*" )?to="#settings-popup" defer :disabled="settingsPane\.floating !== '([^']+)'">\s*<div v-settings-section="'([^']+)'"/g
  ), ([, floating, section]) => [floating, section]);
  assert.deepEqual(moved, [['titles', 'titles'], ['linear-axis', 'linear-axis'], ['circular-axis', 'circular-axis'], ['dinucleotide', 'dinucleotide']]);
  assert.equal(html.match(/to="#settings-popup"/g).length, moved.length);
  assert.equal(html.match(/id="settings-popup"/g).length, 1);
  const app = html.slice(html.indexOf('<div id="app"'), html.indexOf('<footer data-app-footer'));
  assert.match(app, /<\/main>\s*(?:<!--[\s\S]*?-->\s*)?<div id="settings-popup" v-show="settingsPane\.floating"/, 'History listens on the app root');
});

test('the catalog finds the element of a section in the settings', () => {
  const found = {};
  assert.equal(settingsSectionElement({ querySelector: (selector) => (selector === '[data-settings-section="titles"]' ? found : null) }, 'titles'), found);
});
