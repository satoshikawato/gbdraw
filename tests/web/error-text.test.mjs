// User-facing error text in the UI language (D-17). The wording owner
// (`normalizeUserFacingError`, R6) returns the English `summary` and the parts
// it is composed of; `errorText` (utils/error-text.js) is the one display of
// those parts. These tests check that the English display equals the composed
// English, that every part text is a catalog key the catalog check knows, and
// that a locale translates each part, nested labels included, while values
// (identifiers, counts, Python's words) stay as written.
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import test from 'node:test';

import { diagnosticError, generationRecoveryGuidance, liveEditFailure, normalizeUserFacingError, operationErrorTitle, renderErrorText }
  from '../../gbdraw/web/js/utils/error-normalization.js';
import { errorText } from '../../gbdraw/web/js/utils/error-text.js';
import { SOURCE_LOCALE, setUiLocale, t, tm } from '../../gbdraw/web/js/utils/i18n.js';
import { errorTextKeys, errorVocabulary } from './helpers/error-text-keys.mjs';

const vocabulary = await errorVocabulary();
const keys = await errorTextKeys();

// Every code with every reason, and the locators, guidance and track room
// sentences that a context can add.
const corpus = () => {
  const models = [];
  for (const code of vocabulary.codes) {
    models.push(normalizeUserFacingError({ code }));
    for (const reason of vocabulary.reasons) models.push(normalizeUserFacingError({ code, context: { reason } }));
  }
  const contexts = [
    { field: 'window', inputOrdinal: 2, row: 3, column: 4, slotIndex: 0, seriesIndex: 1, configPath: 'objects.scale.interval' },
    { field: 'start', sessionTable: vocabulary.sessionTables[0], diagramMode: 'linear', queryRecordIndex: 0, subjectRecordIndex: 1 },
    { diagramMode: 'circular', featureCaption: 'dnaA (CDS)', featureId: 'f0000dead', featureType: 'CDS', featureStart: 1, featureEnd: 9, featureCount: 3 },
    { innerPx: 181, outerPx: 209, keptInsidePx: 34, keptInsideTracks: 'gc_content,gc_skew', slotTrack: 'ticks', reason: 'CANNOT_FIT' },
    { keptInsidePx: 12, slotTrack: 'gc_skew' },
    { slotTrack: 'ticks', keptInsidePx: 5, keptInsideTracks: 'gc_content' },
    { labelFitCounts: 'Overlapping label pairs: 2. Crossing leader pairs: 0.', exceptionType: 'ValueError', positionUnit: 'python-character', position: 4 },
    { columnCount: 5, recordCount: 2, reason: 'SESSION_TOO_NEW', sessionVersion: 99, latestSessionVersion: 47 },
    { reason: 'SESSION_VERSION_UNSUPPORTED', sessionVersion: 30, readSessionVersions: '27–33, 39–42, and 44' },
    { field: 'schema', reason: 'FIELDS' }
  ];
  for (const code of ['INPUT_INVALID', 'TRACK_LAYOUT', 'COMPARISON_INPUT', 'DECORATION_CONTINUITY', 'LABEL_LAYOUT']) {
    for (const context of contexts) models.push(normalizeUserFacingError({ code, context }));
  }
  for (const field of vocabulary.labelledFields) models.push(normalizeUserFacingError({ code: 'INPUT_INVALID', context: { field } }));
  models.push(normalizeUserFacingError({ code: 'INPUT_INVALID', operation: 'session-load', context: { field: 'schema', reason: 'JSON_FORMAT' } }));
  models.push(normalizeUserFacingError({ code: 'INPUT_INVALID', operation: 'align', context: { field: 'direction', reason: 'SOURCE_VIEW_CONFLICT' } }));
  return models;
};

/** Every text of a part and of the parts in its values. */
const partTexts = (part) => [
  ...(part.verbatim ? [] : [part.text]),
  ...Object.values(part.params || {}).flatMap((value) => (typeof value === 'object' ? partTexts(value) : []))
];

test('in English the display of the parts is the composed summary and note', () => {
  const models = corpus();
  assert.ok(models.length > 1000, `${models.length} models`);
  for (const model of models) {
    assert.equal(errorText(model.summaryParts), model.summary, model.code);
    const live = liveEditFailure(model);
    assert.equal(errorText(live.noteParts), live.note, model.code);
    assert.equal(renderErrorText(live.noteParts), live.note);
  }
});

test('every part text of the error text is a catalog key of the wording owner', () => {
  const unknown = new Set();
  for (const model of corpus()) {
    for (const parts of [model.summaryParts, liveEditFailure(model).noteParts]) {
      for (const text of parts.flatMap(partTexts)) if (!keys.has(text)) unknown.add(text);
    }
  }
  for (const operation of [...vocabulary.operations, 'unknown']) assert.ok(keys.has(operationErrorTitle(operation)), operation);
  for (const recovery of vocabulary.recoveries) assert.ok(keys.has(generationRecoveryGuidance(recovery)), recovery);
  assert.deepEqual([...unknown], []);
  console.log(`${keys.size} error text keys`);
});

test('a locale translates each part and the labels inside it; values stay as written', async () => {
  // A catalog that marks every key, so the test sees which text was looked up.
  const catalog = Object.freeze(Object.fromEntries([...keys].map((key) => [key, `«${key}»`])));
  try {
    await setUiLocale('ja', async () => ({ default: catalog }));
    const model = normalizeUserFacingError({ code: 'INPUT_INVALID',
      context: { field: 'window', reason: 'POSITIVE', row: 3, configPath: 'objects.scale.interval' } });
    assert.equal(model.summary, 'An input value is invalid. Line 3. Setting: objects.scale.interval. Field: Window.'
      + ' Use a finite value greater than zero.', 'the summary stays English in every locale');
    assert.equal(errorText(model.summaryParts), '«An input value is invalid.» «Line 3.» «Setting: objects.scale.interval.»'
      + ' «Field: «Window».» «Use a finite value greater than zero.»');
    const unlabelled = normalizeUserFacingError({ code: 'INPUT_INVALID', context: { field: 'record_selector' } });
    assert.equal(errorText(unlabelled.summaryParts), '«An input value is invalid.» «Field: record_selector.»',
      'a field without a label is an identifier');
    const room = normalizeUserFacingError({ code: 'TRACK_LAYOUT',
      context: { reason: 'CANNOT_FIT', keptInsidePx: 34, keptInsideTracks: 'gc_content,gc_skew', slotTrack: 'ticks',
        labelFitCounts: 'Overlapping label pairs: 2.' } });
    const remedies = ['Clear Show Coordinate Scale (with Custom Track Slots on, remove the Ticks row).',
      'Select Hide GC Content, or reduce GC Content Width (with Custom Track Slots on, the row width).',
      'Select Hide GC Skew, or reduce GC Skew Width (with Custom Track Slots on, the row width).'];
    assert.equal(room.summary, 'A circular track does not fit. Move the track, reduce widths, disable conflicting labels,'
      + ' or place it on the other side of the Axis. Below it, 34 px are kept for the tracks inside it, among them GC content'
      + ` and GC skew. To make room, do one of the following. ${remedies.join(' ')} Overlapping label pairs: 2.`);
    // A list of names is nested `{first} and {second}` parts. Each remedy is a whole sentence with
    // its control names inside it, so a translator never gets a fragment (review item 6).
    // Python's label fit counts are shown as sent.
    assert.equal(errorText(room.summaryParts), '«A circular track does not fit.» «Move the track, reduce widths, disable'
      + ' conflicting labels, or place it on the other side of the Axis.» «Below it, 34 px are kept for the tracks inside it,'
      + ` among them ««GC content» and «GC skew»».» «To make room, do one of the following.» ${remedies.map((text) => `«${text}»`).join(' ')}`
      + ' Overlapping label pairs: 2.');
    const table = normalizeUserFacingError({ code: 'TABLE_INVALID', context: { sessionTable: 'default-colors', row: 2 } });
    assert.match(errorText(table.summaryParts), /«Session table: «Default colors».»/);
    const live = liveEditFailure(normalizeUserFacingError({ code: 'WORKER_INIT' }));
    assert.equal(live.note, 'Live edit failed: direct edits already applied are kept; geometry may still need updating.'
      + ' The diagram runtime could not start. Save a Session and reload before retrying. Retry the live edit or use Generate.');
    assert.equal(errorText(live.noteParts), '«Live edit failed: direct edits already applied are kept; geometry may still need updating.»'
      + ' «The diagram runtime could not start. Save a Session and reload before retrying.» «Retry the live edit or use Generate.»');
    assert.equal(tm(operationErrorTitle('generate')), '«Generation Error»');
    assert.equal(tm(generationRecoveryGuidance('preserved')), '«Showing the last successful Result; your changes were not applied.»');
  } finally {
    await setUiLocale(SOURCE_LOCALE);
  }
});

test('a template shows an error model through errorText, never its English summary or note', () => {
  const html = readFileSync(new URL('../../gbdraw/web/index.html', import.meta.url), 'utf8');
  const bindings = [...html.matchAll(/\b[\w.?]*(?:error|Error|model)\??\.(?:summary|note)\b(?!Parts)/g)].map((match) => match[0]);
  assert.deepEqual(bindings, []);
  assert.match(html, /errorText\(labelReflowLastError\.noteParts\)/);
  assert.match(html, /errorText\(model\.summaryParts\)/);
});

// Review L4: i18n.js `t` / `tm` and the normalizer's `renderErrorText` fill `{name}` placeholders
// the same way, in English and in a locale that moves them, for the same vectors.
test('i18n substitution and the error text renderer fill the same placeholders the same way', async () => {
  const vectors = [
    ['Line {row}.', { row: 3 }], ['Pair: #{query} to #{subject}.', { query: 1, subject: 2 }], ['Sequence {n}.', { n: 0 }],
    ['Field: {field}.', {}], ['{a} and {a}', { a: 'x' }], ['Feature: {feature}.', { feature: '{row} $& $1' }],
    ['At {0} and { row } and {row}.', { row: 1 }], ['Value {toString}.', {}]
  ];
  const catalog = Object.freeze(Object.fromEntries(vectors.map(([text]) => [text, `[${[...text.matchAll(/\{[^}]*\}/g)].reverse().join('|')}]`])));
  for (const [text, params] of vectors) assert.equal(renderErrorText([{ text, params }]), t(text, params), text);
  try {
    await setUiLocale('ja', async () => ({ default: catalog }));
    for (const [text, params] of vectors) {
      assert.ok(tm(text, params).startsWith('['), 'the locale catalog is read');
      assert.equal(renderErrorText([{ text, params }], tm), tm(text, params), text);
    }
  } finally {
    await setUiLocale(SOURCE_LOCALE);
  }
});

// L2: the 1000-character cap keeps whole parts, so the English summary, the shown text, and the
// live-edit note are one cut and no sentence ends in the middle.
test('a summary above the cap keeps whole sentences and its display equals it', () => {
  const context = { reason: 'DERIVED_KEY_MISMATCH', field: 'window', inputOrdinal: 2, row: 3, column: 4, slotIndex: 0, seriesIndex: 1,
    diagramMode: 'circular', sessionTable: 'label-overrides', queryRecordIndex: 0, subjectRecordIndex: 1, featureCount: 3,
    featureCaption: 'dnaA chromosomal replication initiator protein DnaA, ATPase domain (CDS) [1]',
    innerPx: 181, outerPx: 209, keptInsidePx: 34, keptInsideTracks: 'gc_content,gc_skew', slotTrack: 'ticks',
    labelFitCounts: 'Overlapping label pairs: 12. Crossing leader pairs: 30. Labels outside: 4.',
    exceptionType: 'ValueError', positionUnit: 'python-character', position: 4, columnCount: 5, recordCount: 2 };
  const model = normalizeUserFacingError({ code: 'TRACK_LAYOUT', context });
  const dropped = ' Python exception: ValueError.';
  assert.ok(model.summary.length <= 1000 && model.summary.length + dropped.length > 1000, `${model.summary.length}`);
  assert.ok(model.summary.endsWith(' Loaded records: 2.'), model.summary.slice(-40));
  assert.equal(errorText(model.summaryParts), model.summary);
  const live = liveEditFailure(model);
  assert.equal(errorText(live.noteParts), live.note);
  const short = normalizeUserFacingError({ code: 'TRACK_LAYOUT', context }, { summaryLimit: 90 });
  assert.equal(short.summary, 'A circular track does not fit. Diagram: Circular. Session table: Label overrides.');
  assert.equal(errorText(short.summaryParts), short.summary);
  // A first sentence longer than the limit is cut as English text and shown as written.
  const cut = normalizeUserFacingError({ code: 'TRACK_LAYOUT' }, { summaryLimit: 8 });
  assert.deepEqual(cut.summaryParts, [{ text: 'A circul', verbatim: true }]);
  assert.equal(cut.summary, 'A circul');
});

// OV-1389: the overlay corrections name the placement options a track row offers, as the
// template labels them, and the LOSAT ring correction uses the UI's name of the rings.
test('the overlay and ring corrections name the options and rings the UI shows', () => {
  const html = readFileSync(new URL('../../gbdraw/web/index.html', import.meta.url), 'utf8');
  assert.match(html, /<option value="overlay">\{\{ t\('Overlay anchor'\) \}\}<\/option>/, 'Circular annotation placement');
  assert.match(html, /<option value="overlay">\{\{ t\('On axis'\) \}\}<\/option>/, 'Linear track placement');
  assert.match(html, /comparison rings/);
  const summary = (reason) => normalizeUserFacingError({ code: 'TRACK_INVALID', context: { reason } }).summary;
  for (const reason of ['OVERLAY_SIDE', 'OVERLAY_ANCHOR']) {
    assert.match(summary(reason), /Overlay anchor \(Circular\) or On axis \(Linear\)/, reason);
    assert.doesNotMatch(summary(reason), /Use Overlay|an Overlay annotation/, reason);
  }
  assert.match(summary('RING_LOSAT_PROGRAM'), /for comparison rings\.$/);
});
