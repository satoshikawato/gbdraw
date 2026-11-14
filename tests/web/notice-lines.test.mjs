// OV-3495, OV-3562: a preview notice above the figure lists its lines up to NOTICE_LINE_LIMIT; the rest opens on
// request, so a long list (56 record pairs of one LOSATP run, every skipped row of an annotation table, every
// comparison whose table IDs name no displayed record) never pushes the figure out of view. The three notices use
// one component, so the limit, the "and {count} more" line and the disclosure have one owner.
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import test from 'node:test';
import { loadPageVue } from './helpers/vue-runtime.mjs';

const Vue = loadPageVue();
globalThis.window = { Vue, setTimeout, clearTimeout, innerWidth: 1440 };
const { NOTICE_LINE_LIMIT, NoticeLines } = await import('../../gbdraw/web/js/components.js');

const element = (tag) => ({ tag, attributes: {}, children: [], parent: null, text: '' });
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
  remove: (child) => { if (child.parent) child.parent.children.splice(child.parent.children.indexOf(child), 1); },
  parentNode: (node) => node.parent,
  nextSibling: (node) => (node.parent ? node.parent.children[node.parent.children.indexOf(node) + 1] ?? null : null),
  patchProp: (node, key, previous, next) => {
    if (next === null || next === undefined || next === false) delete node.attributes[key];
    else node.attributes[key] = String(next);
  }
};
const all = function* all(node) {
  yield node;
  for (const child of node.children) yield* all(child);
};
const textOf = (node) => [...all(node)].filter((inside) => inside.tag === '#text').map((inside) => inside.text).join('').replace(/\s+/g, ' ').trim();

/** Renders `<notice-lines>` with a slot that shows each line as a `<p>`, as the notices do. */
const render = (lines) => {
  const root = element('div');
  const { createRenderer, h } = Vue;
  const app = createRenderer(nodeOps).createApp({
    render: () => h(NoticeLines, { lines }, { default: ({ line }) => h('p', { 'data-line': line }, line) })
  });
  app.mount(root);
  const details = [...all(root)].find((node) => node.tag === 'details') ?? null;
  const restLines = details ? [...all(details)].filter((node) => node.tag === 'p').map(textOf) : [];
  const shownLines = [...all(root)].filter((node) => node.tag === 'p').map(textOf).filter((line) => !restLines.includes(line));
  return { root, details, shownLines, restLines, unmount: () => app.unmount() };
};

test('a notice over the line limit shows its first lines and "and {count} more", which holds the rest', () => {
  assert.ok(NOTICE_LINE_LIMIT >= 2 && NOTICE_LINE_LIMIT <= 5, 'a few lines, not the list');
  const lines = Array.from({ length: 59 }, (_, index) => `Comparison ${index + 1}`);
  const shown = render(lines);
  assert.deepEqual(shown.shownLines, lines.slice(0, NOTICE_LINE_LIMIT - 1));
  assert.ok(shown.details, 'the rest is behind a disclosure');
  const summary = [...all(shown.details)].find((node) => node.tag === 'summary');
  assert.equal(textOf(summary), `and ${59 - NOTICE_LINE_LIMIT + 1} more`);
  assert.deepEqual(shown.restLines, lines.slice(NOTICE_LINE_LIMIT - 1));
  shown.unmount();
});

test('a notice that fits the line limit shows every line and no disclosure', () => {
  const lines = Array.from({ length: NOTICE_LINE_LIMIT }, (_, index) => `Row ${index + 1}`);
  const shown = render(lines);
  assert.deepEqual(shown.shownLines, lines);
  assert.equal(shown.details, null);
  shown.unmount();
});

// #2002 review: each notice is an atomic polite live region. Opening the disclosure adds the rest to the
// accessibility tree inside it, which would read the whole notice again; the rest is a region of its own with
// aria-live off, so opening it announces nothing. The summary stays in the notice's live region.
test('the rest of a notice is outside its live announcements', () => {
  const shown = render(['A', 'B', 'C', 'D', 'E']);
  const box = [...all(shown.details)].find((node) => node.attributes['data-notice-lines-box'] !== undefined);
  assert.equal(box.attributes['aria-live'], 'off');
  assert.equal(shown.details.attributes['aria-live'], undefined);
  const summary = [...all(shown.details)].find((node) => node.tag === 'summary');
  assert.equal(summary.attributes['aria-live'], undefined);
  assert.ok(![...all(summary)].some((node) => node.attributes['aria-live'] === 'off'));
  shown.unmount();
});

// The three notices of the preview render their lines through the one component: no notice keeps a list of its
// own, a second disclosure, or a second copy of the "and {count} more" line.
test('each preview notice renders its lines through notice-lines', () => {
  const html = readFileSync(new URL('../../gbdraw/web/index.html', import.meta.url), 'utf8');
  const block = (testId) => {
    const start = html.indexOf(`data-testid="${testId}"`);
    assert.ok(start > 0, testId);
    return html.slice(start, html.indexOf('data-testid=', start + 1));
  };
  for (const [testId, lines] of [
    ['annotation-resolution-notice', 'annotationWarnings'],
    ['comparison-id-notice', 'comparisonWarnings'],
    ['input-summary-notice', 'group.warnings']
  ]) {
    const markup = block(testId);
    assert.match(markup, new RegExp(`<notice-lines :lines="${lines.replace('.', '\\.')}"`), testId);
    assert.doesNotMatch(markup, /<details|moreLine\(|v-for="\(warning/, testId);
  }
  assert.doesNotMatch(html, /moreLine\(/, 'the page does not write "and {count} more" itself');
});
