// @ts-check
/**
 * What an exported Interactive SVG shows for each span of a match (OV-3542): Python's export finds the record a
 * span reads, or why it cannot read one, and writes it into the file; the file's viewer reads what Python wrote and
 * holds no rule of its own. Each case of tests/fixtures/embedded-match-span-cases.json goes through Python's export
 * (`embed_catalog_item`, which the CLI export and, part by part, the Web export run) and then through the viewer's own
 * functions, as the file runs them: what Python embeds must be the case's `embedded`, and what the viewer shows its
 * `shown` (the bases of a span, or the sentence of its reason).
 */
import assert from 'node:assert/strict';
import { execFileSync } from 'node:child_process';
import { readFileSync } from 'node:fs';
import test from 'node:test';

import * as schema from '../../gbdraw/web/js/contracts/schema.generated.js';
import * as slots from '../../gbdraw/web/js/services/drawn-attributes.generated.js';
import { STANDALONE_INTERACTIVE_SCRIPT } from '../../gbdraw/web/js/services/standalone-interactivity-assets.js';

const CASES = JSON.parse(readFileSync(new URL('../fixtures/embedded-match-span-cases.json', import.meta.url), 'utf8')).cases;

/** What Python's export embeds for every case, in one run. */
const pythonEmbeds = () => JSON.parse(execFileSync(process.env.PYTHON || 'python', ['-c', `
import json, sys
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from gbdraw.render.interactive_svg import CatalogSequenceRecords
from gbdraw.web_support.feature_catalog import embed_catalog_item

def record(row):
    return SeqRecord(Seq(row["sequence"]), id=row["id"], name=row["name"])

out = []
for case in json.load(sys.stdin):
    keys = [f"key-{index}" for index in range(len(case["records"]))]
    records = CatalogSequenceRecords(
        mode=case["mode"], records=tuple(map(record, case["records"])), unit_keys=tuple(keys),
        comparison_records=tuple(tuple(map(record, rows)) for rows in case["comparisonRecords"]),
    )
    item = embed_catalog_item({"recordKeys": keys, "biologicalFeatures": [], "comparisonMatches": case["matches"]}, records)
    out.append({key: item.get(key) for key in ("sequenceSources", "matchSpanReads")})
print(json.dumps(out))
`], { encoding: 'utf8', input: JSON.stringify(CASES) }));

/** The viewer script as an export writes it: Python's values in its slots. */
const SCRIPT = [
  [slots.STANDALONE_MATCH_FRAGMENT_OWN_ATTRIBUTES_SLOT, slots.MATCH_FRAGMENT_OWN_ATTRIBUTES],
  [slots.STANDALONE_MATCH_SPAN_REASON_TEXTS_SLOT, slots.MATCH_SPAN_REASON_TEXTS]
].reduce((text, [slot, value]) => (slot ? text.split(String(slot)).join(JSON.stringify(value)) : text),
  STANDALONE_INTERACTIVE_SCRIPT);

/** The source of each function the script defines, by name. */
const FUNCTIONS = new Map();
for (const match of SCRIPT.matchAll(/\n  function ([A-Za-z_$][\w$]*)\(/g)) {
  const start = /** @type {number} */ (match.index) + 1;
  let depth = 0;
  for (let index = SCRIPT.indexOf('{', start); index < SCRIPT.length; index += 1) {
    if (SCRIPT[index] === '{') depth += 1;
    if (SCRIPT[index] === '}') depth -= 1;
    if (depth === 0) {
      FUNCTIONS.set(match[1], SCRIPT.slice(start, index + 1));
      break;
    }
  }
}

/** A function of the script with every function of the script it calls. @param {string[]} names */
const withCallees = (names) => {
  const seen = new Set();
  const visit = (/** @type {string} */ name) => {
    if (seen.has(name) || !FUNCTIONS.has(name)) return;
    seen.add(name);
    for (const [, callee] of FUNCTIONS.get(name).matchAll(/\b([A-Za-z_$][\w$]*)\s*\(/g)) visit(callee);
  };
  names.forEach(visit);
  return [...seen].map((name) => FUNCTIONS.get(name)).join('\n');
};

/** The script's own top-level `var` statement for one name (a slot written in). @param {string} name */
const scriptVar = (name) => {
  const match = new RegExp(`\\n  var ${name} = [^\\n]*;`).exec(SCRIPT);
  return match ? match[0] : '';
};

/**
 * The spans the viewer shows for an embedded item, as the file's script computes them: its catalog expansion, then
 * the sequence bundle of each match.
 * @param {Record<string, unknown>} item
 */
const viewerSpans = (item) => new Function('item', `
  var metadata = null;
  var orthogroupsById = new Map();
  var featuresById = new Map();
  ${scriptVar('MATCH_SPAN_REASON_TEXTS')}
  ${withCallees(['expandCatalogItem', 'buildEmbeddedMatchBundle'])}
  var payload = expandCatalogItem(item, 'rich');
  var sequenceSources = payload.sequence_sources;
  var matchSpanReads = payload.match_span_reads || {};
  return payload.matches.map(function (match) {
    return { id: match.id, entries: buildEmbeddedMatchBundle(match).entries };
  });
`)(item);

test('Python\'s export embeds the source or the reason of each match span', () => {
  const embeds = pythonEmbeds();
  CASES.forEach((/** @type {any} */ testCase, /** @type {number} */ index) => {
    const embedded = { sequenceSources: [], matchSpanReads: null, ...Object.fromEntries(
      Object.entries(embeds[index]).filter(([, value]) => value !== null)) };
    assert.deepEqual(embedded, testCase.embedded, testCase.name);
  });
});

test('the viewer shows what Python embedded, with no lookup of its own', () => {
  const failures = CASES.flatMap((/** @type {any} */ testCase) => {
    const item = { recordKeys: testCase.records.map((/** @type {unknown} */ _row, /** @type {number} */ unit) => `key-${unit}`),
      biologicalFeatures: [], comparisonMatches: testCase.matches, ...testCase.embedded };
    const shown = Object.fromEntries(viewerSpans(item).map(({ id, entries }) => [id, Object.fromEntries(
      entries.map((/** @type {any} */ entry) => [entry.role, entry.available
        ? { sequence: entry.fasta.split('\n').slice(1).join('') }
        : { reason: entry.reason }]))]));
    return Object.entries(testCase.shown).flatMap(([id, roles]) => Object.entries(roles)
      .filter(([role, expected]) => JSON.stringify(shown[id]?.[role]) !== JSON.stringify(expected))
      .map(([role, expected]) => `${testCase.name}: ${id} ${role} shows ${JSON.stringify(shown[id]?.[role])}, not ${JSON.stringify(expected)}`));
  });
  assert.deepEqual(failures, []);
});
