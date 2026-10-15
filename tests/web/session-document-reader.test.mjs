// readSessionDocument (E2b-1): the band reader drops a content part it cannot
// read with a row whose path is the part's key, and R4 (`restoresNothing`)
// refuses a Session whose every held content part was dropped. E2b-2: the
// hand-written checks of the input files and the committed Result sets repair
// before the dialog what Save refuses.
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import test from 'node:test';
import { gunzipSync } from 'node:zlib';
import {
  readSessionAuthorityParts, readSessionDocument, restoresNothing, SESSION_TOP_LEVEL_AUTHORITY, validateCurrentComparisonAuthority,
  validateSessionAuthorityInventory
} from '../../gbdraw/web/js/services/session-authority.js';
import { validateWebFileBindings } from '../../gbdraw/web/js/services/session-resource-backing.js';

const fixture = (name) => {
  const bytes = readFileSync(new URL(`../fixtures/sessions/${name}`, import.meta.url));
  return JSON.parse(name.endsWith('.gz') ? gunzipSync(bytes) : bytes);
};
const V46 = 'legend-renamed-rows.v46.gbdraw-session.json.gz';
const V44 = 'lambda_basic_linear.v44-schema8.gbdraw-session.json.gz';
const rows = (document) => readSessionDocument(document).repairs.map(({ code, path }) => `${code}@${path}`);

test('a content part of another kind is dropped with one row at its key', () => {
  for (const [name, expected] of [
    [V46, {
      resources: 'session-part-dropped', results: 'session-part-dropped', editorState: 'session-part-dropped',
      modes: 'session-part-dropped', otherModeResult: 'session-part-dropped', orthogroupState: 'session-field-dropped',
      files: 'session-field-dropped', losatCache: 'cache-entry-dropped', losatDerivedCache: 'cache-entry-dropped',
      proteinIdentityManifest: 'cache-entry-dropped'
    }],
    [V44, {
      resources: 'session-part-dropped', results: 'session-part-dropped', editorState: 'session-part-dropped',
      features: 'session-field-dropped', orthogroupState: 'session-field-dropped', files: 'session-field-dropped',
      modes: 'session-field-dropped', otherModeResult: 'session-field-dropped', losatCache: 'cache-entry-dropped',
      losatDerivedCache: 'cache-entry-dropped', proteinIdentityManifest: 'cache-entry-dropped'
    }]
  ]) {
    const base = fixture(name);
    assert.deepEqual(rows(base), [], name);
    for (const [key, code] of Object.entries(expected)) {
      assert.notEqual(SESSION_TOP_LEVEL_AUTHORITY[key].part, null, key);
      // A binding that names a dropped resource is dropped with its own row (A6b).
      const own = rows({ ...structuredClone(base), [key]: true }).filter((row) => !row.startsWith('file-entry-dropped@webFiles.bindings.'));
      assert.deepEqual(own, [`${code}@${key}`], `${name} ${key}`);
    }
  }
});

test('R4: a Session whose every held content part was dropped restores nothing', () => {
  const stored = fixture(V46);
  const held = Object.keys(stored).filter((key) => SESSION_TOP_LEVEL_AUTHORITY[key]?.part
    && ![undefined, null].includes(stored[key]) && !(typeof stored[key] === 'object' && Object.keys(stored[key]).length === 0));
  assert.ok(held.includes('results') && held.includes('modes'), held.join());
  const dropped = (keys) => keys.map((path) => ({ code: 'session-part-dropped', path, effect: 'dropped' }));
  assert.equal(restoresNothing(stored, dropped(held)), true);
  assert.equal(restoresNothing(stored, dropped(held.filter((key) => key !== 'modes'))), false);
  // A row inside a part does not drop the part.
  assert.equal(restoresNothing(stored, dropped(held.map((key) => `${key}[0]`))), false);
  assert.equal(restoresNothing(stored, []), false);
});

test('R4 (E-OD21, OV-864): the row-less drops of a request Load cannot read count, as the settings-only reading drops them', () => {
  // The one producer of `request-dropped` (K22) drops the request at its part
  // key (or a place in it) and loads the file as settings only: the Results,
  // inputs, feature list and caches go with no row of their own.
  const cli = fixture('se08-main-linear-cli.v42.gbdraw-session.json.gz');
  const web = fixture(V44);
  assert.ok(!('config' in cli) && !('modes' in cli) && 'config' in web);
  for (const path of ['renderRequest', 'renderRequest.mode', 'renderRequest.layout']) {
    const k22 = [{ code: 'request-dropped', path, effect: 'dropped' }];
    // A CLI Session holds no settings: nothing is restored, as with a null request.
    assert.equal(restoresNothing(cli, k22), true, path);
    assert.equal(restoresNothing({ ...cli, renderRequest: null }, []), true);
    // A Web Session keeps its settings; dropped too, nothing is left.
    assert.equal(restoresNothing(web, k22), false, path);
    assert.equal(restoresNothing(web, [...k22, { code: 'session-part-dropped', path: 'config', effect: 'dropped' }]), true, path);
    assert.equal(restoresNothing({ ...structuredClone(fixture(V46)), modes: {} }, k22), true, path);
  }
  // Without that row, the same file restores its request.
  assert.equal(restoresNothing(cli, []), false);
});

test('R4 (OV-865): a part whose every held member was dropped, or whose members hold nothing, restores nothing', () => {
  const settingsOnly = { format: 'gbdraw-session', version: 46, renderRequest: null, title: 't' };
  const slices = (paths) => paths.map((path) => ({ code: 'session-part-dropped', path, effect: 'dropped' }));
  const both = { ...settingsOnly, modes: { circular: 5, linear: 5 } };
  assert.equal(restoresNothing(both, slices(['modes.circular', 'modes.linear'])), true);
  assert.equal(restoresNothing(both, slices(['modes.circular'])), false);
  assert.equal(restoresNothing({ ...settingsOnly, modes: { circular: 5 } }, slices(['modes.circular'])), true);
  assert.equal(restoresNothing({ ...settingsOnly, modes: { circular: null, linear: {} } }, []), true);
  // Every input file dropped one by one is the part dropped.
  const stored = fixture(V44_BOUND);
  const parts = Object.keys(stored).filter((key) => SESSION_TOP_LEVEL_AUTHORITY[key]?.part && key !== 'resources');
  const ids = Object.keys(stored.resources);
  const rest = slices(parts);
  assert.equal(restoresNothing(stored, [...rest, ...ids.map((id) => ({ code: 'file-entry-dropped', path: `resources.${id}`, effect: 'dropped' }))]), true);
  assert.equal(restoresNothing(stored, rest), false);
});

const V44_BOUND = 'inactive-class-m.v44.gbdraw-session.json.gz';
const V44_CIRCULAR = 'HmmtDNA_basic_circular.v44-schema8.gbdraw-session.json.gz';

test('the input files are read before the dialog: a resource or binding that cannot be read is dropped with a row (A8, A6b)', () => {
  const base = fixture(V44_BOUND);
  const [id] = Object.keys(base.resources);
  // A checksum or lastModified the app recomputes is dropped with no row.
  const recomputed = structuredClone(base);
  Object.assign(recomputed.resources[id], { checksum: 'abc', lastModified: -1 });
  const read = readSessionDocument(recomputed);
  assert.deepEqual(read.repairs, []);
  assert.equal(Object.hasOwn(read.value.resources[id], 'checksum'), false);
  assert.equal(Object.hasOwn(read.value.resources[id], 'lastModified'), false);
  assert.equal(recomputed.resources[id].checksum, 'abc', 'the stored document is not changed');
  // A resource that cannot be read: its row, and the row of the binding that named it.
  const unreadable = structuredClone(base);
  unreadable.resources[id].size = -1;
  unreadable.resources.unnamed = 5; // nothing names it: no row (E-OD06)
  const dropped = readSessionDocument(unreadable);
  assert.deepEqual(dropped.repairs.map(({ code, path }) => `${code}@${path}`),
    [`file-entry-dropped@resources.${id}`, 'file-entry-dropped@webFiles.bindings.c_gb']);
  assert.deepEqual(Object.keys(dropped.value.resources), []);
  assert.equal(dropped.value.webFiles.bindings.c_gb, null);
  // An inventory field no writer wrote, and a list element that names no resource.
  const listed = structuredClone(base);
  listed.webFiles.bindings.c_conservation_blasts = [{ resourceId: 'zzz', name: 'x', type: 'y', lastModified: 0 }];
  listed.webFiles.bindings.extra = 1;
  const read2 = readSessionDocument(listed);
  assert.deepEqual(read2.repairs.map(({ code, path }) => `${code}@${path}`),
    ['session-field-dropped@webFiles.bindings.extra', 'file-entry-dropped@webFiles.bindings.c_conservation_blasts[0]']);
  assert.deepEqual(read2.value.webFiles.bindings.c_conservation_blasts, []);
  assert.equal(Object.hasOwn(read2.value.webFiles.bindings, 'extra'), false);
  // Save refuses each of these as before (E-OD02).
  for (const document of [unreadable, listed]) {
    assert.throws(() => validateWebFileBindings(document.webFiles, document.resources));
  }
});

test('the legacy files of a Session 40+ are named with the inputs they held (A6c)', () => {
  const document = { ...fixture(V44), files: { c_gb: { name: 'phage.gb', data: 'x' }, linearSeqs: [{ gff: { name: 'phage.gff' } }] } };
  const [repair] = readSessionDocument(document).repairs;
  assert.deepEqual([repair.code, repair.path, repair.label], ['session-field-dropped', 'files', 'phage.gb, phage.gff']);
});

test('a committed Result set is read as far as it can be restored before the dialog (A7b, A7g, A9, A19, A20)', () => {
  const base = fixture(V44_CIRCULAR);
  const read = (document) => {
    const repairs = [];
    const value = readSessionAuthorityParts(document, document.version, {}, repairs);
    return { value, rows: repairs.map(({ code, path }) => `${code}@${path}`) };
  };
  assert.deepEqual(read(base).rows, []);
  const notSvg = structuredClone(base);
  notSvg.results[0].content = 'x';
  const dropped = read(notSvg);
  assert.deepEqual(dropped.rows, ['session-part-dropped@results']);
  assert.deepEqual([dropped.value.results, dropped.value.editorState.featureCatalog, Object.hasOwn(dropped.value, 'runMetadata')], [[], null, false]);
  const catalog = structuredClone(base);
  catalog.editorState.featureCatalog.schema = 1;
  catalog.runMetadata.annotationWarnings = 'x'; // derived again by the next Generate: no row
  const kept = read(catalog);
  assert.deepEqual(kept.rows, ['feature-recovery-failed@editorState.featureCatalog']);
  assert.equal(kept.value.results.length, base.results.length);
  assert.equal(Object.hasOwn(kept.value.runMetadata, 'annotationWarnings'), false);
  // The other mode's set: without its request it is dropped as a set; with a
  // fault of its own, as a part (the part key is the path, for R4).
  const linear = fixture(V44);
  const withOther = (other) => ({ ...structuredClone(base), version: 46, otherModeResult: other });
  const other = { renderRequest: linear.renderRequest, results: linear.results, editorState: { featureCatalog: null } };
  const otherRows = ({ rows }) => rows.filter((row) => row.includes('otherModeResult'));
  const noRequest = read(withOther({ ...other, renderRequest: undefined }));
  assert.deepEqual(otherRows(noRequest), ['other-result-set-dropped@otherModeResult.renderRequest']);
  assert.equal(Object.hasOwn(noRequest.value, 'otherModeResult'), false);
  // R4 counts the set dropped for its request as the part `otherModeResult`.
  const onlyOther = { format: 'gbdraw-session', version: 46, otherModeResult: other };
  assert.equal(restoresNothing(onlyOther, [{ code: 'other-result-set-dropped', path: 'otherModeResult.renderRequest', effect: 'dropped' }]), true);
  assert.deepEqual(otherRows(read(withOther(other))), ['session-part-dropped@otherModeResult']);
  // Save refuses each fault as before (E-OD02).
  for (const document of [notSvg, catalog]) {
    assert.throws(() => validateSessionAuthorityInventory(document, document.version));
  }
});

test('run notices that do not match their Results are dropped before the dialog with no row; Save refuses them (E-OD16)', () => {
  const base = fixture(V44_CIRCULAR);
  const read = (document) => {
    const repairs = [];
    const value = readSessionAuthorityParts(document, document.version, {}, repairs);
    return { value, rows: repairs.map(({ code, path }) => `${code}@${path}`) };
  };
  const notice = { recordKey: 'r', biologicalFeatureId: 'f', status: 'absent', kinds: ['placement'], resultIndex: 99 };
  for (const [field, value] of [['featureIdentityNotices', [notice]], ['comparisonWarnings', {}], ['annotationWarnings', 'x']]) {
    const stale = structuredClone(base);
    stale.runMetadata = { ...stale.runMetadata, [field]: value };
    // Derived run metadata: the next Generate reports it again.
    const kept = read(stale);
    assert.deepEqual(kept.rows, [], field);
    assert.equal(Object.hasOwn(kept.value.runMetadata, field), false, field);
    assert.equal(kept.value.results.length, base.results.length, field);
    assert.doesNotThrow(() => validateSessionAuthorityInventory(kept.value, kept.value.version, {}, { load: true }), field);
    // Save refuses it (E-OD02).
    assert.throws(() => validateSessionAuthorityInventory(stale, stale.version), field);
  }
  // The other mode's set: its stale notices are dropped and the set is kept.
  const gallery = (name) => JSON.parse(readFileSync(new URL(`../../gbdraw/web/gallery/sessions/${name}`, import.meta.url), 'utf8'));
  const circular = gallery('HmmtDNA_basic_circular.gbdraw-session.json');
  const linear = gallery('lambda_basic_linear.gbdraw-session.json');
  const other = {
    renderRequest: linear.renderRequest, results: linear.results,
    editorState: { featureCatalog: linear.editorState.featureCatalog, alignmentResetReceipt: null,
      legend: { originalOrder: [], originalColors: {} }, originalSvgStroke: { color: null, width: null } },
    ui: { selectedResultIndex: 0 }, runMetadata: { ...linear.runMetadata, featureIdentityNotices: [notice], comparisonWarnings: {} }
  };
  const both = { ...circular, otherModeResult: other };
  const kept = read(both);
  assert.deepEqual(kept.rows.filter((row) => row.includes('otherModeResult')), []);
  assert.ok(kept.value.otherModeResult);
  assert.equal(Object.hasOwn(kept.value.otherModeResult.runMetadata, 'featureIdentityNotices'), false);
  assert.equal(Object.hasOwn(kept.value.otherModeResult.runMetadata, 'comparisonWarnings'), false);
});

test('a Session 27-39 Result set Load could not admit is read again without its Results (OV-789, E-OD17)', () => {
  for (const name of ['scale-interval-negative-linear-cli.v30.gbdraw-session.json.gz', 'BGC0000708-BGC0000713.v39.gbdraw-session.json.gz']) {
    const session = JSON.parse(gunzipSync(readFileSync(new URL(`../fixtures/sessions/${name}`, import.meta.url))));
    const read = (uncomposed) => {
      const repairs = [];
      const value = readSessionAuthorityParts(session, session.version, {}, repairs, uncomposed);
      return { value, rows: repairs.map(({ code, path }) => `${code}@${path}`) };
    };
    assert.equal(read([]).value, session);
    const dropped = read(['']);
    assert.deepEqual(dropped.rows, ['session-part-dropped@results'], name);
    assert.deepEqual([dropped.value.results, Object.hasOwn(dropped.value, 'runMetadata'), Object.hasOwn(dropped.value.ui, 'selectedResultIndex')],
      [[], false, false], name);
    // The feature state recovered from the Results goes with them; the settings and edits stay.
    assert.deepEqual(['biologicalFeatures', 'extractedFeatures', 'featureCatalog', 'featureRecordIds', 'featureSelectorSafetyScope']
      .filter((field) => Object.hasOwn(dropped.value.features, field)), [], name);
    assert.deepEqual(Object.keys(dropped.value.features).filter((field) => !Object.hasOwn(session.features, field)), []);
    assert.equal(dropped.value.config, session.config);
  }
  const v39 = JSON.parse(gunzipSync(readFileSync(new URL('../fixtures/sessions/BGC0000708-BGC0000713.v39.gbdraw-session.json.gz', import.meta.url))));
  const kept = readSessionAuthorityParts(v39, 39, {}, [], ['']).features;
  assert.equal(kept.featureColorOverrides, v39.features.featureColorOverrides);
  assert.equal(kept.labelTextFeatureOverrides, v39.features.labelTextFeatureOverrides);
});

test('a Result set Load\'s composition transform could not compose is read again without its Results (A20, E-OD16)', () => {
  const gallery = (name) => JSON.parse(readFileSync(new URL(`../../gbdraw/web/gallery/sessions/${name}`, import.meta.url), 'utf8'));
  const circular = gallery('HmmtDNA_basic_circular.gbdraw-session.json');
  const linear = gallery('lambda_basic_linear.gbdraw-session.json');
  const other = {
    renderRequest: linear.renderRequest, results: linear.results,
    editorState: { featureCatalog: linear.editorState.featureCatalog, alignmentResetReceipt: null,
      legend: { originalOrder: [], originalColors: {} }, originalSvgStroke: { color: null, width: null } },
    ui: { selectedResultIndex: 0 }, runMetadata: linear.runMetadata
  };
  const both = { ...circular, otherModeResult: other };
  const read = (uncomposed) => {
    const repairs = [];
    const value = readSessionAuthorityParts(both, both.version, {}, repairs, uncomposed);
    return { value, rows: repairs.map(({ code, path }) => `${code}@${path}`) };
  };
  assert.deepEqual(read([]).rows, []);
  // The shown set: its Results and the members that name them.
  const shown = read(['']);
  assert.deepEqual(shown.rows, ['session-part-dropped@results']);
  assert.deepEqual([shown.value.results, shown.value.editorState.featureCatalog, Object.hasOwn(shown.value, 'runMetadata')], [[], null, false]);
  assert.ok(shown.value.otherModeResult);
  // The other mode's set needs a Result: it is dropped as a part.
  const dropped = read(['otherModeResult.']);
  assert.deepEqual(dropped.rows, ['session-part-dropped@otherModeResult']);
  assert.equal(Object.hasOwn(dropped.value, 'otherModeResult'), false);
  assert.equal(dropped.value.results.length, circular.results.length);
});

// E4d (review round 1, E-OD29). R1-M1: a comparison file binding the reader
// dropped leaves the edge with no file, and the load-mode check accepts it.
test('R1-M1: a comparison binding dropped at Load is not refused by the comparison check afterwards', () => {
  const session = fixture('q-frame-main-web-upload.v42.gbdraw-session.json.gz');
  session.webFiles.bindings.linearComparisons[0].file.resourceId = 'zz';
  const read = readSessionDocument(session);
  assert.deepEqual(read.repairs.map(({ code, path }) => `${code}@${path}`),
    ['file-entry-dropped@webFiles.bindings.linearComparisons[0].file']);
  assert.equal(read.value.webFiles.bindings.linearComparisons[0].file, null);
  assert.doesNotThrow(() => validateCurrentComparisonAuthority(read.value, { load: true }));
  // Save still requires each binding's file.
  assert.throws(() => validateCurrentComparisonAuthority(read.value), /requires a file resource binding/);
});

// R1-M2: Save refuses the committed artifacts of a settings-only Session as
// dev did (E-OD02); Load drops them with no row (E-OD05).
test('R1-M2: Save refuses a settings-only Session holding run metadata, older artifacts or cache entries of another kind', () => {
  const base = fixture('settings-only.v42.json.gz');
  assert.doesNotThrow(() => validateSessionAuthorityInventory(structuredClone(base), 42));
  for (const [field, value] of [['runMetadata', 'abc'], ['legacyArtifacts', [1]], ['losatCache', { entries: 'abc' }]]) {
    const session = { ...structuredClone(base), [field]: value };
    assert.throws(() => validateSessionAuthorityInventory(session, 42), /committed render artifacts/, field);
    const read = readSessionDocument(session);
    // A cache entry is dropped with a code the dialog does not list (E-OD05).
    assert.ok(read.repairs.every(({ code }) => code === 'cache-entry-dropped'), field);
    assert.notDeepEqual(read.value[field], value, field);
  }
});
