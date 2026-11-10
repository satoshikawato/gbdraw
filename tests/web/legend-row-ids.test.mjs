// OV-3070: the Legend list keys its rows by an identity a rename keeps (`legendRowIds` in
// services/legend-draft-rows.js), not by their places, so a rename that moves a row keeps the row's element and the
// keyboard focus in it. The Legend list template keys its rows by these ids.
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import test from 'node:test';

const { legendRowIds: legendRowIdsOfKeys, rekeyLegendRowIds } = await import('../../gbdraw/web/js/services/legend-draft-rows.js');
const { legendEntryKey } = await import('../../gbdraw/web/js/services/drawing/paint.js');

/** The ids of a list's rows, read as the Legend list reads them: by the rows' keys. */
const legendRowIds = (rows) => legendRowIdsOfKeys(rows.map(legendEntryKey));

const row = (key, caption = key) => ({ caption, originalCaption: key, color: '#cccccc', xPos: 0, yPos: 0, featureIds: [] });
// A copy of a list as a Session Load or a History restore installs it: new row objects of the same keys.
const restored = (rows) => JSON.parse(JSON.stringify(rows));
/** The id of the row of `key` in `rows`. */
const idOf = (rows, key) => legendRowIds(rows)[rows.findIndex((entry) => entry.originalCaption === key)];

test('a rename that keeps the key, a move, and a Load keep each row\'s id', () => {
  const rows = [row('tRNA-1'), row('CDS-1'), row('rRNA-1')];
  const ids = legendRowIds(rows);
  assert.equal(new Set(ids).size, 3, 'each row has its own id');
  assert.deepEqual(legendRowIds([row('tRNA-1'), row('CDS-1', 'Coding'), row('rRNA-1')]), ids, 'a rename writes the caption only');
  assert.deepEqual(legendRowIds([rows[2], rows[0], rows[1]]), [ids[2], ids[0], ids[1]], 'a moved row keeps its id');
  assert.deepEqual(legendRowIds(restored(rows)), ids, 'a Load or a restore of the same keys');
});

test('a rule row\'s rename moves its id to the new key, and Undo and Redo find it again', () => {
  const before = [row('AT skew (+)-2'), row('CDS-2'), row('Cytochrome-2')];
  const renamed = idOf(before, 'Cytochrome-2');
  const others = [idOf(before, 'AT skew (+)-2'), idOf(before, 'CDS-2')];
  rekeyLegendRowIds(new Map([['Cytochrome-2', 'Cytochromes-2']]));
  // The rule commit lists the row under its new key, in Python's place for it, and adds a row.
  const after = [row('Cytochromes-2'), row('other proteins-2'), row('AT skew (+)-2'), row('CDS-2')];
  assert.equal(idOf(after, 'Cytochromes-2'), renamed, 'the renamed row keeps its id at its new place');
  assert.deepEqual([idOf(after, 'AT skew (+)-2'), idOf(after, 'CDS-2')], others);
  assert.ok(![renamed, ...others].includes(idOf(after, 'other proteins-2')), 'an added row has its own id');
  assert.equal(idOf(restored(before), 'Cytochrome-2'), renamed, 'Undo brings the old key back with the id');
  assert.equal(idOf(restored(after), 'Cytochromes-2'), renamed, 'Redo');
});

test('re-keys move all at once, and an old key listed beside its row\'s new key gets an id of its own', () => {
  const rows = [row('alpha-3'), row('beta-3')];
  const [alpha, beta] = legendRowIds(rows);
  rekeyLegendRowIds(new Map([['alpha-3', 'beta-3'], ['beta-3', 'alpha-3']]));
  assert.deepEqual(legendRowIds([row('beta-3'), row('alpha-3')]), [alpha, beta], 'a swap of two keys');
  rekeyLegendRowIds(new Map([['beta-3', 'gamma-3']]));
  // A new row takes the old key "beta-3" while the row renamed "gamma-3" is listed.
  const ids = legendRowIds([row('alpha-3'), row('beta-3'), row('gamma-3')]);
  assert.equal(ids[2], alpha, 'the re-keyed row keeps its id');
  assert.equal(new Set(ids).size, 3);
  assert.deepEqual(legendRowIds([row('gamma-3'), row('alpha-3'), row('beta-3')]), [ids[2], ids[0], ids[1]], 'and keeps it');
});

test('a merge leaves the other rows\' ids, and a row of a key listed twice still has an id of its own', () => {
  const rows = [row('ORF-4'), row('gene-4'), row('repeat-4')];
  const ids = legendRowIds(rows);
  assert.deepEqual(legendRowIds([rows[0], rows[2]]), [ids[0], ids[2]], 'a row merged into another leaves the list');
  assert.deepEqual(legendRowIds(rows), ids, 'and comes back with its id');
  const twice = legendRowIds([row('ORF-4'), row('ORF-4', 'ORF again')]);
  assert.equal(twice[0], ids[0]);
  assert.notEqual(twice[1], twice[0]);
});

test('the Legend list keys its rows by their ids, not by their places', () => {
  const html = readFileSync(new URL('../../gbdraw/web/index.html', import.meta.url), 'utf8');
  const list = html.match(/<div v-for="\(entry, idx\) in legendEntries"[^>]*data-legend-drawer-row/)?.[0];
  assert.ok(list, 'the Legend list row');
  assert.match(list, /:key="legendRowIds\[idx\]"/);
});
