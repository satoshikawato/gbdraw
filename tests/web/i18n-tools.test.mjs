// The translation tools: the key owner (tools/i18n/source-keys.mjs), the
// catalog format (tools/i18n/catalogs.mjs), the sync tool (tools/i18n/sync.mjs)
// and the review page data (tools/i18n/review-page-data.mjs).
import assert from 'node:assert/strict';
import { createHash } from 'node:crypto';
import { cpSync, mkdtempSync, readFileSync, readdirSync, rmSync, writeFileSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import test from 'node:test';

import { CONTEXT_SEPARATOR, UI_LOCALES } from '../../gbdraw/web/js/utils/i18n.js';
import {
  CATALOG_LOCALES, CATALOG_ROOT, catalogProblems, formatCatalog, glossaryMatcher, parseCatalog, parseGlossary, readCatalog
} from '../../tools/i18n/catalogs.mjs';
import { reviewDraftText, reviewKeyId, reviewPageData } from '../../tools/i18n/review-page-data.mjs';
import { roleAt, sourceCalls, sourceKeys } from '../../tools/i18n/source-keys.mjs';
import { previousKeyOf, runSync } from '../../tools/i18n/sync.mjs';

const { keys } = await sourceKeys();

/** A copy of the catalogs that a test may change, and the sync tool run on it. */
const scratch = (t) => {
  const folder = mkdtempSync(join(tmpdir(), 'i18n-sync-'));
  t.after(() => rmSync(folder, { recursive: true, force: true }));
  const root = join(folder, 'catalogs');
  cpSync(CATALOG_ROOT, root, { recursive: true });
  // Start complete whatever dev's catalogs lack: a feature PR adds English only, and the batch before a release translates it.
  for (const locale of CATALOG_LOCALES) {
    const catalog = readCatalog(locale, root);
    for (const key of Object.keys(catalog.entries)) if (!keys.has(key)) delete catalog.entries[key];
    for (const { key, english } of keys.values()) catalog.entries[key] ??= english;
    writeFileSync(join(root, `${locale}.js`), formatCatalog(catalog));
  }
  const run = async (...argv) => {
    const lines = [];
    const code = await runSync(argv, { catalogRoot: root, print: (line) => lines.push(line) });
    return { code, lines };
  };
  const bytes = () => Object.fromEntries(readdirSync(root).sort().map((name) => [name, readFileSync(join(root, name), 'utf8')]));
  const edit = (locale, change) => {
    const catalog = readCatalog(locale, root);
    change(catalog.entries);
    writeFileSync(join(root, `${locale}.js`), formatCatalog(catalog));
  };
  const file = (name, records) => {
    const path = join(folder, name);
    writeFileSync(path, records.map((record) => `${JSON.stringify(record)}\n`).join(''));
    return path;
  };
  return { root, run, bytes, edit, file };
};

test('a key carries each call site with its role and part of the app', () => {
  const html = [
    '<div><button class="x" :title="t(\'Undo\')"><span>{{ t(\'Redo\') }}</span></button>',
    '<input :placeholder="t(\'Search\')"><label>{{ t(\'Name\') }}</label>',
    '<div role="dialog"><p>{{ t(\'Sure?\') }}</p></div><select><option>{{ tx(\'mode\', \'Linear\') }}</option></select>',
    '<help-tip :text="t(\'Help me\')"></help-tip><span :aria-label="t(\'Close\')"></span><p>{{ t(\'Plain\') }}</p></div>'
  ].join('\n');
  const roles = ['Undo', 'Redo', 'Search', 'Name', 'Sure?', 'Linear', 'Help me', 'Close', 'Plain']
    .map((text) => roleAt(html, html.indexOf(`'${text}'`) - (text === 'Linear' ? 12 : 2), true));
  assert.deepEqual(roles, ['title', 'button', 'placeholder', 'label', 'dialog', 'option', 'help', 'aria-label', 'other']);
  const module = 'const a = { label: msg(\'A\'), message: msg(\'B\'), other: msg(\'C\') };';
  assert.deepEqual(['A', 'B', 'C'].map((text) => roleAt(module, module.indexOf(`msg('${text}'`), false)), ['label', 'notice', 'other']);

  const heading = keys.get('Generating Diagram...');
  assert.deepEqual(heading.sites.map(({ file, role, area }) => [file, role, area]), [['index.html', 'heading', 'Page']]);
  assert.match(heading.sites[0].snippet, /t\('Generating Diagram\.\.\.'\)/);
  const plural = [...keys.values()].find(({ plural: pair }) => pair);
  assert.ok(keys.get(plural.plural.one).plural && keys.get(plural.plural.other).plural, 'both forms of a tn() pair name the pair');
  // A key that a t(), tx() or msg() call or a derived list also produces keeps one string, so it names no pair.
  const producers = new Map();
  for (const call of sourceCalls()) if (!('problem' in call)) producers.set(call.key, [...(producers.get(call.key) ?? []), call.name]);
  for (const entry of keys.values()) {
    if (!entry.plural) continue;
    const names = producers.get(entry.key) ?? [];
    assert.ok(names.length === entry.sites.length && names.every((name) => name === 'tn'), `${entry.key}: only tn() calls produce a plural key`);
  }
  const context = keys.get(`size scaling${CONTEXT_SEPARATOR}Linear`);
  assert.equal(context.context, 'size scaling');
  assert.equal(context.sites[0].area, 'Settings: Layout');
  assert.equal(keys.get(`member role${CONTEXT_SEPARATOR}anchor`).sites[0].file, 'js/app/pairwise-match-popup.js');
});

test('every catalog prints back to its own bytes; other text is refused', () => {
  for (const locale of CATALOG_LOCALES) {
    const text = readFileSync(join(CATALOG_ROOT, `${locale}.js`), 'utf8');
    assert.equal(formatCatalog(parseCatalog(text)), text, locale);
  }
  assert.throws(() => parseCatalog('export default {};\n'), /canonical/);
  assert.throws(() => parseCatalog('// x\nexport default Object.freeze({\n  "a": 1\n});\n'), /not a string/);
  assert.throws(() => parseCatalog('// x\nexport default Object.freeze({\n  "a": {"other": 1}\n});\n'), /not a string/);
  assert.equal(formatCatalog({ header: '// h\n', entries: { b: '2', a: '1', B: '3' } }), '// h\nexport default Object.freeze({\n  "B": "3",\n  "a": "1",\n  "b": "2"\n});\n');
  // Plural forms print on the key's row in CLDR category order, and read back to the same bytes.
  const plural = formatCatalog({ header: '// h\n', entries: { n: { other: 'o', few: 'f', many: 'm' } } });
  assert.equal(plural, '// h\nexport default Object.freeze({\n  "n": {"few": "f", "many": "m", "other": "o"}\n});\n');
  assert.equal(formatCatalog(parseCatalog(plural)), plural);
});

test('glossary hints: whole-word terms, longest first, names kept in every locale left out', () => {
  const header = `en\t${CATALOG_LOCALES.join('\t')}\tnote`;
  const row = (en, cell) => [en, ...CATALOG_LOCALES.map(() => cell), ''].join('\t');
  const glossary = parseGlossary([header, row('Record', 'R'), row('Record label', 'RL'), row('LOSAT', 'LOSAT'), row('bp', 'b')].join('\n'));
  const termsOf = glossaryMatcher(glossary);
  assert.deepEqual(termsOf('Show record labels from LOSAT'), [1, 0]);
  assert.deepEqual(termsOf('Recorded 100 bp'), []);
});

test('status counts missing and stale keys; --require-complete fails on either', async (t) => {
  const { run, edit } = scratch(t);
  assert.deepEqual(await run('status', '--require-complete'), {
    code: 0, lines: CATALOG_LOCALES.map((locale) => `${locale}: ${keys.size} of ${keys.size} translated, 0 missing, 0 stale`)
  });
  edit('fr', (entries) => {
    delete entries['Generating Diagram...'];
    entries['Old text'] = 'Ancien texte';
  });
  const { code, lines } = await run('status', '--locale', 'fr', '--require-complete');
  assert.equal(code, 1);
  assert.deepEqual(lines.slice(0, 3), [
    `fr: ${keys.size - 1} of ${keys.size} translated, 1 missing, 1 stale`, '  missing "Generating Diagram..."', '  stale "Old text"'
  ]);
  assert.equal((await run('status', '--locale', 'fr')).code, 0, 'without --require-complete status only reports');
  // A stale key alone: the PR-tier catalog check passes it, the release gate fails it.
  const stale = scratch(t);
  stale.edit('fr', (entries) => { entries['Old text'] = 'Ancien texte'; });
  assert.deepEqual(catalogProblems(keys, readCatalog('fr', stale.root).entries), []);
  const gate = await stale.run('status', '--locale', 'fr', '--require-complete');
  assert.equal(gate.code, 1);
  assert.ok(gate.lines.includes('  stale "Old text"'));
  assert.equal((await run('status', '--locale', 'de')).code, 2);
});

test('export lists the missing keys with what a drafter needs; --all lists every key', async (t) => {
  const { run, edit } = scratch(t);
  assert.deepEqual((await run('export', '--locale', 'ja')).lines, []);
  edit('ja', (entries) => {
    delete entries['Generating Diagram...'];
    entries['Generating diagram...'] = '旧訳';
  });
  const [line, ...more] = (await run('export', '--locale', 'ja')).lines;
  assert.deepEqual(more, []);
  const record = JSON.parse(line);
  assert.equal(record.key, 'Generating Diagram...');
  assert.deepEqual(record.sites.map(({ file, role }) => [file, role]), [['index.html', 'heading']]);
  assert.ok(record.glossary.some(({ en, term }) => en === 'Diagram' && term === '図'));
  assert.equal(record.translation, null);
  assert.deepEqual(record.previous, { key: 'Generating diagram...', translation: '旧訳', status: 'needs-check' });
  const all = (await run('export', '--locale', 'ja', '--all')).lines.map((text) => JSON.parse(text));
  assert.equal(all.length, keys.size);
  const counted = all.find(({ plural }) => plural);
  assert.ok(counted.key === counted.plural.one || counted.key === counted.plural.other);
  assert.equal(previousKeyOf('Save the session', ['Delete all records', 'Save this session']), 'Save this session');
  assert.equal(previousKeyOf('Save the session', ['Something else entirely']), null);
});

test('apply writes the canonical catalog, drops stale keys, and changes no byte when nothing changes', async (t) => {
  const { root, run, bytes, edit, file } = scratch(t);
  const before = bytes();
  const exported = (await run('export', '--locale', 'ko', '--all')).lines.map((text) => JSON.parse(text));
  assert.equal((await run('apply', '--locale', 'ko', file('all.jsonl', exported))).code, 0);
  assert.equal((await run('apply', '--locale', 'ko', file('none.jsonl', []))).code, 0);
  assert.deepEqual(bytes(), before, 'an export applied unchanged, or an empty file, leaves every catalog as it was');

  edit('ko', (entries) => {
    delete entries['Generating Diagram...'];
    entries['Old text'] = '옛 문장';
  });
  const { code, lines } = await run('apply', '--locale', 'ko', file('one.jsonl', [{ key: 'Generating Diagram...', translation: '다이어그램 생성 중…' }]));
  assert.deepEqual([code, lines], [0, ['ko: 1 lines applied (1 changed), 1 stale keys dropped']]);
  const text = bytes()['ko.js'];
  assert.equal(formatCatalog(parseCatalog(text)), text, 'canonical format');
  const after = parseCatalog(text).entries;
  assert.equal(after['Generating Diagram...'], '다이어그램 생성 중…');
  assert.ok(!Object.hasOwn(after, 'Old text'));
  assert.deepEqual(Object.keys(after), Object.keys(readCatalog('ko', root).entries));
  assert.equal(Object.keys(after).length, keys.size);
});

test('plural forms go through export and apply as an object; apply refuses them on a key that is no tn() form', async (t) => {
  const { run, bytes, file } = scratch(t);
  const counted = [...keys.values()].find(({ key, plural }) => plural && key === plural.other && /\{count\}/.test(key));
  const forms = { many: `${counted.key} (many)`, few: `${counted.key} (few)`, other: `${counted.key} (other)` };
  assert.deepEqual((await run('apply', '--locale', 'fr', file('plural.jsonl', [{ key: counted.key, translation: forms }]))).lines,
    ['fr: 1 lines applied (1 changed), 0 stale keys dropped']);
  assert.match(bytes()['fr.js'], new RegExp(`\\n  ${JSON.stringify(counted.key).replace(/[.*+?^${}()|[\]\\]/g, '\\$&')}: \\{"few": `));
  const exported = (await run('export', '--locale', 'fr', '--all')).lines.map((text) => JSON.parse(text));
  assert.deepEqual(exported.find(({ key }) => key === counted.key).translation, { few: forms.few, many: forms.many, other: forms.other });
  const before = bytes();
  assert.deepEqual((await run('apply', '--locale', 'fr', file('again.jsonl', exported))).lines,
    [`fr: ${keys.size} lines applied (0 changed), 0 stale keys dropped`], 'the export applied again changes nothing');
  assert.deepEqual(bytes(), before);

  const { code, lines } = await run('apply', '--locale', 'fr', file('refused.jsonl', [
    { key: 'Generating Diagram...', translation: { one: 'Génération…', other: 'Générations…' } },
    { key: counted.key, translation: { few: '{n}', other: forms.other } }
  ]));
  assert.equal(code, 1);
  assert.deepEqual(lines.map((line) => line.replace(/^refused .*?:(\d+): /, 'refused $1: ')), [
    'refused 1: "Generating Diagram...": plural forms are allowed only for a form of a tn() pair',
    `refused 2: ${JSON.stringify(counted.key)} few: placeholders {n} differ from the English {count}`,
    'fr: nothing written; fix the 2 refused line(s)'
  ]);
  assert.deepEqual(bytes(), before, 'nothing written');
});

test('apply refuses unknown keys, placeholder changes, repeated keys and bad lines, and then writes nothing', async (t) => {
  const { run, bytes, file } = scratch(t);
  const before = bytes();
  const counted = [...keys.values()].find(({ english }) => /\{count\}/.test(english));
  const path = file('bad.jsonl', [
    { key: 'Not a key of the sources', translation: 'x' },
    { key: counted.key, translation: '{n}' },
    { key: 'Generating Diagram...', translation: 'A' },
    { key: 'Generating Diagram...', translation: 'B' },
    { key: 'Undo', translation: 3 },
    { key: 'Cancel', translation: '  ' }
  ]);
  writeFileSync(path, `${readFileSync(path, 'utf8')}{oops\nnull\n[1]\n{"key": 5, "translation": "x"}\n{"key": "Cancel"}\n`);
  const { code, lines } = await run('apply', '--locale', 'es', path);
  assert.equal(code, 1);
  assert.deepEqual(lines.map((line) => line.replace(path, 'bad.jsonl')), [
    'refused bad.jsonl:1: "Not a key of the sources" is not a key of the sources',
    `refused bad.jsonl:2: ${JSON.stringify(counted.key)}: placeholders {n} differ from the English {count}`,
    'refused bad.jsonl:4: "Generating Diagram..." is given twice',
    'refused bad.jsonl:5: not an object with a string "key" and a string or object "translation"',
    'refused bad.jsonl:6: "Cancel": the translation is empty',
    'refused bad.jsonl:7: not a JSON line',
    'refused bad.jsonl:8: not an object with a string "key" and a string or object "translation"',
    'refused bad.jsonl:9: not an object with a string "key" and a string or object "translation"',
    'refused bad.jsonl:10: not an object with a string "key" and a string or object "translation"',
    'refused bad.jsonl:11: not an object with a string "key" and a string or object "translation"',
    'es: nothing written; fix the 10 refused line(s)'
  ]);
  assert.deepEqual(bytes(), before);
  assert.equal((await run('apply', '--locale', 'es')).code, 2, 'apply needs a file');
});

test('the review page data: ids are sha1(key)[:12], drafts are the catalogs, the screenshot map keeps only current keys', async () => {
  const kept = reviewKeyId('Cancel');
  const shots = { scale: 2, shots: [{ file: 'shots/s001.webp', w: 10, h: 10, state: 'x' }], keys: { [kept]: [[0, 1, 2, 3, 4]], ffffffffffff: [[0, 1, 1, 1, 1]] } };
  const files = await reviewPageData({ shots });
  assert.deepEqual([...files.keys()].sort(), [
    ...CATALOG_LOCALES.flatMap((locale) => [`draft-${locale}.json`, `notes-${locale}.json`]), 'glossary.json', 'keys.json', 'locales.json', 'shots.json', 'style.md'
  ].sort());
  assert.deepEqual(JSON.parse(files.get('locales.json')).map(({ tag }) => tag), [...CATALOG_LOCALES]);
  assert.equal(JSON.parse(files.get('locales.json'))[0].label, UI_LOCALES.find(({ tag }) => tag === CATALOG_LOCALES[0]).name);
  assert.equal(kept, createHash('sha1').update('Cancel').digest('hex').slice(0, 12));
  const rows = JSON.parse(files.get('keys.json'));
  assert.equal(rows.length, keys.size);
  const linear = rows.find(({ key }) => key === `size scaling${CONTEXT_SEPARATOR}Linear`);
  assert.deepEqual([linear.id, linear.en, linear.sense, linear.area, linear.kind], [reviewKeyId(linear.key), 'Linear', 'size scaling', 'Settings: Layout', 'option']);
  for (const row of rows) {
    assert.deepEqual(Object.keys(row), ['id', 'en', 'sense', 'key', 'area', 'where', 'kind', 'uses', 'ctx', 'g']);
    assert.doesNotMatch(row.where, /:\d+$/, 'no line number in the published data');
  }
  const glossary = JSON.parse(files.get('glossary.json'));
  assert.ok(rows.every(({ g }) => g.every((index) => glossary[index])));
  const ja = JSON.parse(files.get('draft-ja.json'));
  assert.equal(ja[reviewKeyId('Cancel')], readCatalog('ja').entries.Cancel);
  // The page shows one string per key: plural forms are joined in CLDR category order.
  assert.equal(reviewDraftText('{count} rows'), '{count} rows');
  assert.equal(reviewDraftText({ other: 'o {count}', many: 'm {count}', few: 'f {count}' }), 'few: f {count} | many: m {count} | other: o {count}');
  assert.deepEqual(JSON.parse(files.get('shots.json')).keys, { [kept]: [[0, 1, 2, 3, 4]] });
  assert.equal(files.get('style.md'), readFileSync(join(import.meta.dirname, '..', '..', 'docs', 'internal', 'I18N_STYLE.md'), 'utf8'));
});
