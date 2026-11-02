// D-02 of the UI copy revision (PD-OI-037 rev 3): a timing note that applies
// on Generate, wholly or in part, shows only a clock icon; its text is the
// icon's tooltip and accessible description, and its accessible name is
// "Applies on Generate". A change that applies at once shows no mark, so no
// row renders "Live edit". The component is rendered with the page's Vue build
// into a plain element tree.
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import test from 'node:test';
import { loadPageVue } from './helpers/vue-runtime.mjs';

const Vue = loadPageVue();
globalThis.window = { Vue, setTimeout, clearTimeout, innerWidth: 1440 };
const { SettingNote, HelpTip } = await import('../../gbdraw/web/js/components.js');
const { SETTING_APPLICATION, SETTING_NOTE_APPLIES, settingApplicationNote } = await import('../../gbdraw/web/js/services/setting-application.js');
const { t, tm } = await import('../../gbdraw/web/js/utils/i18n.js');

const html = readFileSync(new URL('../../gbdraw/web/index.html', import.meta.url), 'utf8');
const helpTipTemplate = html.match(/<script type="text\/x-template" id="help-tip-template">([\s\S]*?)<\/script>/)[1];

// A minimal host: elements with attributes, children and text. HelpTip reads
// the name of its holder (`closest`) only when it has no `label`.
// `layout` gives every element the box a test sets (OV-2449).
const layout = { rect: { top: 400, bottom: 416, left: 100, width: 16 }, height: 40 };
const element = (tag) => ({
  tag, attributes: {}, children: [], parent: null, text: '', handlers: {}, style: {}, closest: () => null, matches: () => false,
  getBoundingClientRect: () => layout.rect, get offsetHeight() { return layout.height; }
});
const nodeOps = {
  createElement: (tag) => element(tag),
  createText: (text) => ({ ...element('#text'), text }),
  createComment: (text) => ({ ...element('#comment'), text }),
  setText: (node, text) => { node.text = text; },
  setElementText: (node, text) => { node.children = [{ ...element('#text'), text, parent: node }]; },
  insert: (child, parent, anchor) => {
    child.parent = parent;
    const at = anchor ? parent.children.indexOf(anchor) : -1;
    if (at < 0) parent.children.push(child); else parent.children.splice(at, 0, child);
  },
  remove: (child) => {
    const parent = child.parent;
    if (parent) parent.children.splice(parent.children.indexOf(child), 1);
  },
  parentNode: (node) => node.parent,
  nextSibling: (node) => (node.parent ? node.parent.children[node.parent.children.indexOf(node) + 1] ?? null : null),
  querySelector: () => body,
  patchProp: (node, key, previous, next) => {
    if (/^on[A-Z]/.test(key)) { node.handlers[key] = next; return; }
    if (key === 'style') { node.style = { ...next }; return; }
    if (next === null || next === undefined || next === false) delete node.attributes[key];
    else node.attributes[key] = String(next);
  }
};
let body = element('body');
const all = function* all(node) {
  yield node;
  for (const child of node.children) yield* all(child);
};
const textOf = (node) => [...all(node)].filter((inside) => inside.tag === '#text').map((inside) => inside.text).join('').replace(/\s+/g, ' ').trim();
const classes = (node) => (node.attributes.class || '').split(/\s+/);

/** Renders `<setting-note>` with the given props and committed-request predicate. */
const render = (props, applies = () => true) => {
  body = element('body');
  const root = element('div');
  const { createRenderer, h } = Vue;
  // The page's HelpTip reads its template from the document; here from index.html.
  const note = { ...SettingNote, components: { ...SettingNote.components, HelpTip: { ...HelpTip, template: helpTipTemplate } } };
  const app = createRenderer(nodeOps).createApp({ render: () => h(note, props) });
  app.config.globalProperties.t = t;
  app.config.globalProperties.tm = tm;
  app.provide(SETTING_NOTE_APPLIES, () => applies);
  app.mount(root);
  return { root, body, unmount: () => app.unmount() };
};

test('a note that applies on Generate is a clock icon named "Applies on Generate" with its text as the description', () => {
  for (const row of SETTING_APPLICATION.filter((candidate) => candidate.application === 'generate')) {
    const { root, body: page, unmount } = render({ setting: row.id });
    const nodes = [...all(root)];
    const wrapper = nodes.find((node) => node.attributes['data-applies-on-generate'] === row.id);
    assert.ok(wrapper, `${row.id}: the wrapper keeps data-applies-on-generate`);
    const button = [...all(wrapper)].find((node) => node.tag === 'button');
    assert.ok(button, `${row.id}: the clock is a button`);
    assert.equal(button.attributes['aria-label'], 'Applies on Generate', row.id);
    assert.ok([...all(button)].some((node) => node.tag === 'i' && classes(node).includes('ph-clock')), `${row.id}: a clock icon`);
    assert.equal(textOf(wrapper), '', `${row.id}: no visible text beside the icon`);
    const description = [...all(page)].find((node) => node.attributes.id === button.attributes['aria-describedby']);
    assert.ok(description, `${row.id}: the description exists`);
    assert.equal(textOf(description), row.note || 'Applies on Generate', `${row.id}: the full note is the description`);
    unmount();
  }
});

test('no live row renders "Live edit"; a live row with a Generate part shows the clock', () => {
  for (const row of SETTING_APPLICATION.filter((candidate) => candidate.application === 'live')) {
    const { root, body: page, unmount } = render({ setting: row.id });
    const shown = `${textOf(root)} ${textOf(page)}`;
    assert.doesNotMatch(shown, /Live edit/, row.id);
    const button = [...all(root)].find((node) => node.tag === 'button');
    // feature-visibility: leaving a feature out of protein search waits for Generate.
    assert.equal(Boolean(button), settingApplicationNote(row.id).display === 'clock', row.id);
    if (button) assert.equal(button.attributes['aria-label'], 'Applies on Generate', row.id);
    unmount();
  }
  // A Multi-Record Canvas change: the positions wait for Generate (partialNote).
  const partial = render({ setting: 'legend' }, (path) => path !== 'diagramOptions.output.legend');
  const button = [...all(partial.root)].find((node) => node.tag === 'button');
  assert.ok(button, 'a partial note shows the clock');
  const description = [...all(partial.body)].find((node) => node.attributes.id === button.attributes['aria-describedby']);
  assert.match(textOf(description), /the position applies on Generate after the Multi-Record Canvas change/);
  partial.unmount();
});

test('an export-only note stays a hint line, and a class that follows state shows the clock only while it waits', () => {
  const exportNote = render({ setting: 'font-embedding' });
  assert.equal(textOf(exportNote.root), 'Applies to SVG and Interactive SVG exports, not to Generate.');
  assert.ok(![...all(exportNote.root)].some((node) => node.tag === 'button'));
  exportNote.unmount();
  const waiting = render({ setting: 'colors', application: 'generate' });
  assert.ok([...all(waiting.root)].some((node) => node.attributes['data-applies-on-generate'] === 'colors'));
  waiting.unmount();
  const live = render({ setting: 'colors', application: 'live' });
  assert.equal(textOf(live.root), '');
  assert.ok(![...all(live.root)].some((node) => node.tag === 'button'));
  live.unmount();
});

test('a long tooltip that would run off the top of the window opens below its icon (OV-2449)', async () => {
  const open = async (rect, height) => {
    layout.rect = rect;
    layout.height = height;
    const shown = render({ setting: 'track-slots' });
    const button = [...all(shown.root)].find((node) => node.tag === 'button');
    button.handlers.onClick({ stopPropagation() {}, preventDefault() {} });
    for (let i = 0; i < 3; i += 1) await Vue.nextTick();
    const tooltip = [...all(shown.body)].find((node) => node.attributes.role === 'tooltip');
    assert.ok(tooltip, 'the tooltip opens');
    const style = { ...tooltip.style };
    shown.unmount();
    return style;
  };
  // 330 px of note above an icon 130 px from the top would start off screen.
  const below = await open({ top: 130, bottom: 146, left: 100, width: 16 }, 330);
  assert.equal(below.transform, 'translate(-50%, 0)');
  assert.equal(below.top, '158px');
  // With room above, it stays above.
  const above = await open({ top: 500, bottom: 516, left: 100, width: 16 }, 330);
  assert.equal(above.transform, 'translate(-50%, -100%)');
  assert.equal(above.top, '488px');
});
