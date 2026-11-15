// OV-1280 (D-19): Remove on a file card leaves the drawing as if the file had
// never been loaded, keeping every hand edit made before or after it.
import test from 'node:test';
import assert from 'node:assert/strict';

const {
  isCardRow,
  importPriorityRows,
  markedRows,
  parseBlacklistWords,
  removeImportedPriorityRows,
  withImportedBlacklistWords,
  importedDefaultColors,
  withoutImportedBlacklistWords,
  withoutImportedDefaultColors,
  withoutImportedRows
} = await import('../../gbdraw/web/js/services/file-imports.js');
const {
  readColorRules,
  readFileImports,
  readQualifierPriorityRules,
  readWhitelistRules
} = await import('../../gbdraw/web/js/services/config-field-readers.js');
const {
  createDefaultAdv,
  createDefaultForm,
  validateCurrentWriterActiveConfig
} = await import('../../gbdraw/web/js/services/session-active-config-contract.js');
const { prepareSpecificColorImport } = await import('../../gbdraw/web/js/services/specific-color-rules.js');

const palette = { CDS: '#aaaaaa', tRNA: '#bbbbbb', rRNA: '#cccccc' };

test('the Default colors take back the user color or the palette color, and keep a later edit', () => {
  // CDS had a user color, tRNA the palette's; misc is not in the palette.
  const before = { ...palette, CDS: '#111111' };
  const userColorOf = (type) => (before[type] !== undefined && before[type] !== palette[type] ? before[type] : null);
  const fileColors = { CDS: '#ff0000', tRNA: '#00ff00', rRNA: '#0000ff', misc: '#123456' };
  const imported = importedDefaultColors(fileColors, userColorOf);
  const colors = { ...before, ...fileColors };
  assert.deepEqual(imported.CDS, { color: '#ff0000', handColor: '#111111' });
  assert.deepEqual(imported.tRNA, { color: '#00ff00', handColor: null });
  // A later edit of rRNA is the user's.
  const edited = { ...colors, rRNA: '#999999' };
  assert.deepEqual(withoutImportedDefaultColors(edited, imported, palette), { ...palette, CDS: '#111111', rRNA: '#999999' });
  // A color written in another case is the file's color still.
  assert.deepEqual(withoutImportedDefaultColors({ ...colors, CDS: '#FF0000' }, imported, palette).CDS, '#111111');
});

test('the Blacklist text takes back the text before, or the file words from an edited text', () => {
  for (const before of ['', 'hypothetical, unknown', 'putative,', '  spaced  ']) {
    const { text, imported } = withImportedBlacklistWords(before, ['alpha', 'beta']);
    assert.notEqual(text, before);
    assert.equal(withoutImportedBlacklistWords(text, imported), before);
  }
  const { text, imported } = withImportedBlacklistWords('hypothetical, unknown', ['alpha', 'beta']);
  // Words typed after the import stay; a removed or edited file word is the user's.
  assert.equal(withoutImportedBlacklistWords(`${text}, mine`, imported), 'hypothetical, unknown, mine');
  assert.equal(withoutImportedBlacklistWords(text.replace('beta', 'betas'), imported), 'hypothetical, unknown, betas');
  assert.equal(withoutImportedBlacklistWords('alpha, beta\nmine', imported), 'mine');
  // A file without words changes nothing and records nothing.
  assert.deepEqual(withImportedBlacklistWords('a', []), { text: 'a', imported: null });
});

// OV-2792: a Blacklist file and the Blacklist text are cut at the Arabic comma
// too, and Remove on the card cuts the text at the same separators.
test('a Blacklist file, and Remove on its card, cut at the Arabic comma as well', () => {
  assert.deepEqual(parseBlacklistWords('alpha،beta\ngamma, delta'), { words: ['alpha', 'beta', 'gamma', 'delta'], count: 4 });
  assert.deepEqual(parseBlacklistWords('a\r\nb,c'), { words: ['a', 'b', 'c'], count: 3 });
  const { text, imported } = withImportedBlacklistWords('فرضی، ناشناخته', ['alpha', 'beta']);
  assert.equal(text, 'فرضی، ناشناخته, alpha, beta');
  assert.equal(withoutImportedBlacklistWords(`${text}، mine`, imported), 'فرضی، ناشناخته، mine');
  assert.equal(withoutImportedBlacklistWords('alpha،beta،mine', imported), 'mine');
  // A text that ends with the Arabic comma takes no second separator.
  assert.equal(withImportedBlacklistWords('putative،', ['alpha']).text, 'putative،alpha');
  // CRLF is one line break: Remove leaves no lone CR.
  assert.equal(withoutImportedBlacklistWords('mine\r\nalpha\r\nbeta', imported), 'mine');
});

// A Session saved before OV-2792 stored a Blacklist file's words cut at ',' and line breaks only
// and trimmed by trim(): 'a،b' was one word, and 'c\u0085' kept its U+0085. Remove on the card
// reads those words as the Blacklist text is read now, so it takes them out of an edited text too.
test('Remove on the Blacklist card takes out the words an older Session stored unsplit', () => {
  const imported = { words: ['a،b', 'c\u0085'], before: 'mine' };
  // The text as that release wrote it, unchanged: the text before comes back.
  assert.equal(withoutImportedBlacklistWords('mine, a،b, c\u0085', imported), 'mine');
  // Edited after the import: the file's words go, the words typed since stay.
  assert.equal(withoutImportedBlacklistWords('mine, a،b, c\u0085, theirs', imported), 'mine, theirs');
  assert.equal(withoutImportedBlacklistWords('a، b\nc, theirs', imported), 'theirs');
});

test('the Qualifier priority rows take back the user order the file replaced', () => {
  const rows = [{ feat: 'CDS', order: 'gene' }, { feat: 'tRNA', order: 'product' }];
  importPriorityRows(rows, [{ feat: 'CDS', order: 'product,gene' }, { feat: 'rRNA', order: 'note' }, { feat: 'CDS', order: 'locus_tag' }]);
  assert.deepEqual(rows, [
    { feat: 'CDS', order: 'locus_tag', handOrder: 'gene', importedFrom: 'qualifier_priority' },
    { feat: 'tRNA', order: 'product' },
    { feat: 'rRNA', order: 'note', importedFrom: 'qualifier_priority' }
  ]);
  removeImportedPriorityRows(rows);
  assert.deepEqual(rows, [{ feat: 'CDS', order: 'gene' }, { feat: 'tRNA', order: 'product' }]);
});

test('the list rows a file added go; a hand row and an edited row stay', () => {
  const whitelist = [{ feat: 'CDS', qual: 'gene', key: 'mine' }, ...markedRows([{ feat: 'CDS', qual: 'gene', key: 'mine' }], 'whitelist')];
  assert.deepEqual(withoutImportedRows(whitelist, 'whitelist'), [{ feat: 'CDS', qual: 'gene', key: 'mine' }]);
  // D-ORIGIN-1: a Specific rules row of the file is marked `fromFile`, as on
  // main, so a Session up to 47 holds the mark already.
  assert.deepEqual(markedRows([{ val: 'a' }], 't_color'), [{ val: 'a', fromFile: true }]);
  assert.equal(isCardRow({ fromFile: true }, 't_color'), true);
  assert.equal(isCardRow({ importedFrom: 't_color' }, 't_color'), false);
  assert.equal(isCardRow({ fromFile: true }, 'whitelist'), false);
  const { nextRules } = prepareSpecificColorImport('CDS\tgene\talpha\t#112233\tAlpha\n',
    [{ feat: 'tRNA', qual: 'gene', val: 'hand', color: '#445566', cap: '' }]);
  assert.deepEqual(withoutImportedRows(nextRules, 't_color').map((rule) => rule.val), ['hand']);
  // Loading the card's next file replaces the rows its earlier file added.
  const replaced = prepareSpecificColorImport('CDS\tgene\tbeta\t#112233\tBeta\n', nextRules).nextRules;
  assert.deepEqual(replaced.map((rule) => [rule.val, rule.fromFile ?? null]), [['hand', null], ['beta', true]]);
});

test('Load keeps each mark and what the files set; it reads a foreign mark as a hand row', () => {
  assert.deepEqual(readQualifierPriorityRules([
    { feat: 'CDS', order: 'product', importedFrom: 'qualifier_priority', handOrder: 'gene' },
    { feat: 'tRNA', order: 'note', importedFrom: 'whitelist', handOrder: 'gene' }
  ]), {
    value: [{ feat: 'CDS', order: 'product', importedFrom: 'qualifier_priority', handOrder: 'gene' }, { feat: 'tRNA', order: 'note' }],
    lost: true,
    leaves: [{ at: '[1].importedFrom', stored: 'whitelist', used: undefined }, { at: '[1].handOrder', stored: 'gene', used: undefined }]
  });
  assert.deepEqual(readWhitelistRules([{ feat: 'CDS', qual: 'gene', key: 'a', importedFrom: 'whitelist' }]).value,
    [{ feat: 'CDS', qual: 'gene', key: 'a', importedFrom: 'whitelist' }]);
  assert.deepEqual(readWhitelistRules([{ feat: 'CDS', qual: 'gene', key: 'a', importedFrom: 't_color' }]).value,
    [{ feat: 'CDS', qual: 'gene', key: 'a' }]);
  assert.deepEqual(readColorRules([{ feat: 'CDS', qual: 'gene', val: 'a', color: '#112233', cap: '', fromFile: true }]).value,
    [{ feat: 'CDS', qual: 'gene', val: 'a', color: '#112233', cap: '', fromFile: true }]);
  const fileImports = { d_color: { CDS: { color: '#ff0000', handColor: null }, tRNA: { color: '#00ff00', handColor: '#111111' } },
    blacklist: { words: ['alpha'], before: 'unknown' } };
  assert.deepEqual(readFileImports(fileImports), { value: fileImports, lost: false });
  assert.deepEqual(readFileImports(undefined), { value: {}, lost: false });
  assert.deepEqual(readFileImports({ d_color: { CDS: { color: 1 } }, blacklist: { words: 'alpha' }, other: {} }), { value: {}, lost: true });
  // Save writes them in the current Session format.
  assert.doesNotThrow(() => validateCurrentWriterActiveConfig({
    mode: 'circular',
    storedConfig: {
      form: createDefaultForm(), adv: createDefaultAdv('circular'), fileImports,
      rules: [{ feat: 'CDS', qual: 'gene', val: 'a', color: '#112233', cap: '', fromFile: true }],
      qualifierPriorityRules: [{ feat: 'CDS', order: 'product', importedFrom: 'qualifier_priority', handOrder: 'gene' }],
      whitelist: [{ feat: 'CDS', qual: 'gene', key: 'a', importedFrom: 'whitelist' }]
    }
  }));
});
