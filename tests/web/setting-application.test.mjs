// OV-1035..1039 (R1, PD-OI-037 rev 2 / OIC-024): `SETTING_APPLICATION` in
// services/setting-application.js is the one owner of when a setting applies.
// index.html names a setting only by its key: an element with
// `data-setting="<key>"` is the scope of that row, a control belongs to the
// row of its nearest scope (itself or an ancestor), and a timing note is a
// `<setting-note setting="<key>">` that renders the row's note. No timing
// note is written as literal text.
//
// "Settings control": every input, select, textarea and button, and every
// component element except the presentation-only ones, inside an element
// marked `v-settings-section` (the settings pane's only containers; the
// catalog `SETTINGS_SECTIONS` owns their keys), unless it sits in a
// `data-history-ignore` or `data-history-scope="transient"` subtree, which
// the History input adapter (`isIgnoredTarget`) also leaves out because no
// Session setting changes there. Outside the settings pane the controls that
// write draft settings are the Result's Feature placement selects (popup and
// selection toolbar) and the record popup's definition fields; every other
// Result control is an editor edit (live), a review with its own Apply, or a
// dialog.
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { join } from 'node:path';
import test from 'node:test';

import {
  SETTING_APPLICATION, settingApplicationNote, settingApplicationSentence
} from '../../gbdraw/web/js/services/setting-application.js';
import { SETTINGS_SECTIONS } from '../../gbdraw/web/js/services/settings-sections.js';

const html = readFileSync(join(process.cwd(), 'gbdraw', 'web', 'index.html'), 'utf8');
const lineOf = (offset) => html.slice(0, offset).split('\n').length;

const VOID = new Set(['area', 'base', 'br', 'col', 'embed', 'hr', 'img', 'input', 'link', 'meta', 'source', 'track', 'wbr']);
const PRESENTATION_COMPONENTS = new Set(['help-tip', 'settings-section-title', 'setting-note', 'markup-preview']);

// The end of an opening tag, skipping '>' inside quoted attribute values.
const tagEnd = (from) => {
  let quote = null;
  for (let index = from; index < html.length; index += 1) {
    const char = html[index];
    if (quote) {
      if (char === quote) quote = null;
    } else if (char === '"' || char === "'") {
      quote = char;
    } else if (char === '>') {
      return index;
    }
  }
  throw new Error(`Unterminated tag at line ${lineOf(from)}`);
};

const parseAttributes = (text) => {
  const attributes = new Map();
  for (const match of text.matchAll(/([^\s=/"']+)(?:\s*=\s*("[^"]*"|'[^']*'|[^\s>"']+))?/g)) {
    const value = match[2] ?? '';
    attributes.set(match[1], value.replace(/^["']|["']$/g, ''));
  }
  return attributes;
};

// The element tree of index.html outside <script> and comments.
const parseTree = () => {
  const root = { tag: '#root', attributes: new Map(), parent: null, children: [], line: 0 };
  let current = root;
  let index = 0;
  while (index < html.length) {
    const open = html.indexOf('<', index);
    if (open < 0) break;
    if (html.startsWith('<!--', open)) {
      index = html.indexOf('-->', open) + 3;
      continue;
    }
    if (html.startsWith('</', open)) {
      const end = html.indexOf('>', open);
      const tag = html.slice(open + 2, end).trim().toLowerCase();
      let node = current;
      while (node && node.tag !== tag) node = node.parent;
      if (node) current = node.parent;
      index = end + 1;
      continue;
    }
    if (!/[a-zA-Z]/.test(html[open + 1] || '')) {
      index = open + 1;
      continue;
    }
    const end = tagEnd(open);
    const raw = html.slice(open + 1, end);
    const tag = raw.match(/^[^\s/>]+/)[0].toLowerCase();
    const element = {
      tag,
      attributes: parseAttributes(raw.slice(tag.length)),
      parent: current,
      children: [],
      line: lineOf(open)
    };
    current.children.push(element);
    index = end + 1;
    if (tag === 'script') {
      index = html.indexOf('</script>', index) + '</script>'.length;
      continue;
    }
    if (!VOID.has(tag) && !raw.endsWith('/')) current = element;
  }
  return root;
};

const tree = parseTree();
const descendants = function* descendants(element) {
  for (const child of element.children) {
    yield child;
    yield* descendants(child);
  }
};
const ancestorsOrSelf = function* ancestorsOrSelf(element) {
  for (let node = element; node; node = node.parent) yield node;
};
const allElements = [...descendants(tree)];
const scopeOf = (element) => [...ancestorsOrSelf(element)].find((node) => node.attributes.has('data-setting')) || null;
const enclosingScope = (scope) => (scope.parent ? scopeOf(scope.parent) : null);
const ROWS = new Map(SETTING_APPLICATION.map((row) => [row.id, row]));
const describe = (element) => `line ${element.line}: <${element.tag}${element.attributes.has('aria-label')
  ? ` aria-label="${element.attributes.get('aria-label')}"` : ''}>`;

const isControl = (element) => {
  if (['input', 'select', 'textarea', 'button'].includes(element.tag)) {
    return element.attributes.get('type') !== 'hidden';
  }
  return element.tag.includes('-') && !PRESENTATION_COMPONENTS.has(element.tag);
};
const isIgnored = (element) => [...ancestorsOrSelf(element)].some((node) => node.attributes.has('data-history-ignore')
  || node.attributes.get('data-history-scope') === 'transient');

const sections = allElements.filter((element) => element.attributes.has('v-settings-section'));
const settingsControls = sections.flatMap((section) => [...descendants(section)]
  .filter((element) => isControl(element) && !isIgnored(element)));
const RESULT_DRAFT_WRITERS = /featurePlacementActions\.setPlacement\(|setPreviewRecordField\(/;
const resultDraftControls = allElements.filter((element) => !sections.some((section) => [...ancestorsOrSelf(element)].includes(section))
  && RESULT_DRAFT_WRITERS.test(element.attributes.get('@change') || ''));
const scopes = allElements.filter((element) => element.attributes.has('data-setting'));
const notes = allElements.filter((element) => element.tag === 'setting-note');

test('the owner has one row per key with a known application', () => {
  const ids = SETTING_APPLICATION.map((row) => row.id);
  assert.equal(new Set(ids).size, ids.length, 'row keys are unique');
  for (const row of SETTING_APPLICATION) {
    assert.match(row.id, /^[a-z][a-z0-9-]*$/, `${row.id}: a stable key, not a caption`);
    assert.ok(['live', 'generate', 'apply', 'export'].includes(row.application), `${row.id}: ${row.application}`);
    if ('note' in row) assert.equal(row.application, 'generate', `${row.id}: only a Generate row qualifies its note`);
  }
  assert.ok(Object.isFrozen(SETTING_APPLICATION) && SETTING_APPLICATION.every(Object.isFrozen));
});

test('every settings section exists in index.html and has settings controls in a scope', () => {
  assert.deepEqual(
    sections.map((section) => section.attributes.get('v-settings-section').replace(/'/g, '')).sort(),
    SETTINGS_SECTIONS.map((section) => section.key).sort()
  );
  assert.ok(settingsControls.length > 250, `found ${settingsControls.length} settings controls`);
  assert.equal(resultDraftControls.length, 5, resultDraftControls.map(describe).join('\n'));
  const unscoped = [...settingsControls, ...resultDraftControls]
    .filter((element) => !scopeOf(element)).map(describe);
  assert.deepEqual(unscoped, [], 'every settings control belongs to a SETTING_APPLICATION row');
});

test('every scope and note names a row, and every row is used', () => {
  const named = [
    ...scopes.map((element) => [element.attributes.get('data-setting'), element]),
    ...notes.map((element) => [element.attributes.get('setting'), element])
  ];
  assert.deepEqual(named.filter(([id]) => !ROWS.has(id)).map(([id, element]) => `${id} at ${describe(element)}`), []);
  const sentences = [...html.matchAll(/settingApplicationSentence\('([^']+)'\)/g)].map((match) => match[1]);
  assert.deepEqual(sentences.filter((id) => !ROWS.has(id)), []);
  const used = new Set(scopes.map((element) => element.attributes.get('data-setting')));
  assert.deepEqual(SETTING_APPLICATION.filter((row) => !used.has(row.id)).map((row) => row.id), []);
  for (const note of notes) {
    const scope = scopeOf(note);
    assert.equal(scope?.attributes.get('data-setting'), note.attributes.get('setting'), `the note at ${describe(note)} sits in its own scope`);
  }
});

test('every scope that applies on Generate shows its note, and a nested scope that differs shows its own', () => {
  const missing = [];
  for (const scope of scopes) {
    const id = scope.attributes.get('data-setting');
    const row = ROWS.get(id);
    const parent = enclosingScope(scope);
    const differs = parent && ROWS.get(parent.attributes.get('data-setting'))?.application !== row.application;
    if (row.application !== 'generate' && !differs) continue;
    const own = [...descendants(scope)].some((element) => element.tag === 'setting-note'
      && element.attributes.get('setting') === id && scopeOf(element) === scope);
    if (!own) missing.push(`${id} at ${describe(scope)}`);
  }
  assert.deepEqual(missing, []);
  for (const control of resultDraftControls) {
    assert.equal(ROWS.get(scopeOf(control).attributes.get('data-setting')).application, 'generate', describe(control));
  }
});

test('the selection toolbar placement select has the note of the popup Feature placement (OV-1038)', () => {
  const placements = resultDraftControls.filter((element) => /setPlacement\(/.test(element.attributes.get('@change')));
  assert.equal(placements.length, 2);
  for (const control of placements) assert.equal(scopeOf(control).attributes.get('data-setting'), 'feature-placement', describe(control));
});

test('every Generate row renders the note, and live rows render Live edit', () => {
  for (const row of SETTING_APPLICATION) {
    const note = settingApplicationNote(row.id);
    if (row.application === 'generate') {
      assert.equal(note.text, row.note || 'Applies on Generate', row.id);
      assert.equal(note.application, 'generate');
      assert.equal(settingApplicationSentence(row.id), 'Changes apply on Generate.');
    } else if (row.application === 'live') {
      assert.equal(note.text, 'Live edit', row.id);
    } else {
      assert.equal(note.text, '', row.id);
    }
  }
  assert.throws(() => settingApplicationNote('no-such-setting'), /no-such-setting/);
  // The component (components.js `SettingNote`) renders exactly that text and,
  // for a Generate row, the marker attribute the browser tests select.
  const component = readFileSync(join(process.cwd(), 'gbdraw', 'web', 'js', 'components.js'), 'utf8')
    .match(/export const SettingNote = \{([\s\S]*?)\n\};/);
  assert.ok(component, 'components.js defines SettingNote');
  assert.match(component[1], /setup: \(props\) => \(\{ note: settingApplicationNote\(props\.setting\) \}\)/);
  assert.match(component[1], /:data-applies-on-generate="note\.application === 'generate' \? setting : null"/);
  assert.match(component[1], /\{\{ note\.text \}\}/);
});

test('the GC, skew and Depth toggles state today\'s split timing (OV-1036)', () => {
  const toggles = settingsControls.filter((element) => /show_gc|show_skew|show_depth|setCircularGcSuppressed|setCircularSkewSuppressed/
    .test(`${element.attributes.get('v-model') || ''} ${element.attributes.get('@change') || ''}`));
  assert.equal(toggles.length, 6, toggles.map(describe).join('\n'));
  for (const toggle of toggles) assert.equal(scopeOf(toggle).attributes.get('data-setting'), 'track-visibility', describe(toggle));
  assert.match(settingApplicationNote('track-visibility').text, /at once.*Generate/);
});

test('no timing note is written as literal text in index.html', () => {
  const page = html.replace(/<!--[\s\S]*?-->/g, '');
  const literal = [...page.matchAll(/appl(?:y|ies) on (?:the next )?Generate|data-applies-on-generate/gi)]
    .map((match) => `${match[0]} at line ${page.slice(0, match.index).split('\n').length}`);
  assert.deepEqual(literal, []);
});
