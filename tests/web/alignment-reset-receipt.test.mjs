import assert from 'node:assert/strict';
import test from 'node:test';
import { spawnSync } from 'node:child_process';
import { CANONICAL_REQUEST_SCHEMA } from '../../gbdraw/web/js/contracts/schema.generated.js';
const { buildSimilarityAlignmentResetReceipt: build, validateSimilarityAlignmentResetReceipt: admit,
  validateAlignmentResetReceiptShape: admitShape, readStoredSimilarityAlignmentResetReceipt: readStored
} = await import('../../gbdraw/web/js/services/session-active-config-contract.js');
const { promoteCanonicalRenderRequestToCurrent } = await import('../../gbdraw/web/js/services/session-request.js');
const { normalizeUserFacingError } = await import('../../gbdraw/web/js/utils/error-normalization.js');
// OV-130 (R6): a stale or malformed receipt is a diagnostic, not an unknown failure.
const stale = { code: 'ALIGNMENT_RESET_EVIDENCE' };
const source = Buffer.from('source bytes remain biological authority');
const resource = { name: 'record.gbk', kind: 'genbank', encoding: 'base64', data: source.toString('base64'), size: source.length };
const anchor = (recordKey) => ({ recordKey, biologicalFeatureId: `gene-${recordKey}`, sourceFeatureIndex: 1, stableFeatureSvgId: null });
const fixture = () => {
  const before = { resources: { source: resource }, renderRequest: {
    schema: CANONICAL_REQUEST_SCHEMA, mode: 'linear', records: ['a', 'b', 'c'].map(recordKey => ({
      recordKey, source: {kind: 'genbank', resourceId: 'source'}, cardinality: 'exactly_one',
      selector: { kind: 'recordId', id: recordKey }, windows: [],
      display: { isCircular: null, startCoordinate: null }, presentation: { reverseComplement: false }
    })), layout: { recordTranslations: ['a', 'b', 'c'].map(recordKey => ({recordKey, x: 7, y: -3})), similarityAlignment: null }
  }};
  const after = structuredClone(before);
  after.renderRequest.records[0].presentation.reverseComplement = true;
  after.renderRequest.records[1].presentation.reverseComplement = true;
  after.renderRequest.layout.recordTranslations[0].x = -147;
  after.renderRequest.layout.similarityAlignment = { schema: 2, groupId: 'og-1', reference: anchor('a'), records: [
    { recordKey: 'a', status: 'reference', rationale: 'reference', anchor: anchor('a') },
    { recordKey: 'b', status: 'aligned', rationale: 'user_selected', anchor: anchor('b') },
    { recordKey: 'c', status: 'skipped', rationale: 'skipped_by_user', anchor: null }
  ]};
  return { before, after };
};
test('receipt records actual absolute deltas including reference and only sparse reference x', async () => {
  const f = fixture(); const receipt = await build(f);
  assert.deepEqual(Object.keys(receipt).sort(), ['binding', 'directions', 'referenceDeltaX']);
  assert.deepEqual(receipt.directions, [{recordKey:'a', before:false, after:true}, {recordKey:'b', before:false, after:true}]);
  assert.deepEqual(receipt.referenceDeltaX, {recordKey:'a', deltaX:-154});
  assert.equal(await admit(receipt, f.after), receipt);
});
test('manual Reverse, cosmetic changes, translations and stable reorder retain the original evidence', async () => {
  const f = fixture(); const receipt = await build(f); const changed = structuredClone(f.after);
  changed.renderRequest.records[0].presentation.reverseComplement = false;
  changed.renderRequest.records[0].presentation.label = 'new label';
  // Both cardinalities select the same single source record; source bytes,
  // selector, record keys and exact active anchors retain the binding.
  changed.renderRequest.records[0].cardinality = 'all';
  changed.renderRequest.records.reverse(); changed.renderRequest.layout.similarityAlignment.records.reverse();
  changed.renderRequest.layout.recordTranslations.reverse();
  changed.renderRequest.layout.recordTranslations[0].y = 123;
  assert.equal(await admit(receipt, changed), receipt);
  assert.equal(receipt.directions[0].before, false);
});
test('resource ID remapping with unchanged bytes preserves the binding for Save/fresh Load', async () => {
  const f = fixture(); const receipt = await build(f);
  f.after.resources.renamed = f.after.resources.source; delete f.after.resources.source;
  f.after.renderRequest.records.forEach(record => { record.source.resourceId = 'renamed'; });
  assert.equal(await admit(receipt, f.after), receipt);
});
for (const [name, mutate] of [
  ['source bytes', f => { const bytes = Buffer.from('different source'); Object.assign(f.after.resources.source, {data:bytes.toString('base64'),size:bytes.length}); }],
  ['selector', f => { f.after.renderRequest.records[1].selector.id = 'another-record'; }],
  ['crop', f => { f.after.renderRequest.records[1].windows = [{start:2,end:10,reverseComplement:false,column:null}]; }],
  ['plan anchor', f => { f.after.renderRequest.layout.similarityAlignment.records[1].anchor.biologicalFeatureId = 'another-gene'; }],
  ['plan removal', f => { f.after.renderRequest.layout.similarityAlignment = null; }]
]) test(`reject changed ${name} without dropping evidence`, async () => {
  const f = fixture(); const receipt = await build(f); mutate(f);
  await assert.rejects(admit(receipt, f.after), stale);
});
test('missing historical evidence and empty modern delta are distinct', async () => {
  const f = fixture(); assert.equal(await admit(undefined, f.after), null);
  const same = structuredClone(f.after); const receipt = await build({before:same, after:f.after});
  assert.deepEqual(receipt.directions, []); assert.equal(receipt.referenceDeltaX, null);
  assert.equal(await admit(receipt, f.after), receipt);
});
test('malformed modern receipt never becomes absent historical evidence', async () => {
  const f = fixture(); const receipt = await build(f);
  const variants = [ {}, {...receipt, referenceDeltaX:undefined}, {...receipt, binding:'wrong'},
    {...receipt, directions:[...receipt.directions, receipt.directions[0]]},
    {...receipt, directions:[{recordKey:'c',before:false,after:true}]},
    {...receipt, directions:[{recordKey:'a',before:true,after:true}]},
    {...receipt, referenceDeltaX:{recordKey:'b',deltaX:1}},
    {...receipt, referenceDeltaX:{recordKey:'a',deltaX:Infinity}},
    {...receipt, mode:'left'} ];
  for (const invalid of variants) await assert.rejects(admit(invalid, f.after), stale);
});
test('a receipt whose binding is not a string fails the shape check and the public check alike', async () => {
  const f = fixture(); const receipt = await build(f);
  // A one-element array holding the hash coerces to it in the regex test, so only the type check rejects it.
  const arrayBinding = { ...receipt, binding: [receipt.binding] };
  assert.throws(() => admitShape(arrayBinding, f.after.renderRequest), stale);
  await assert.rejects(admit(arrayBinding, f.after), stale);
});
// OV-115 and OV-130: a Linear receipt checked against a Circular request (and any
// stale receipt) reaches Generate and Save as a specific diagnostic with a way out.
test('a stale receipt says what Reset alignment would restore and the ways out', async () => {
  const f = fixture(); const receipt = await build(f);
  const circular = { ...f.after.renderRequest, mode: 'circular' };
  assert.throws(() => admitShape(receipt, circular), stale);
  const changed = structuredClone(f.after); changed.renderRequest.records[1].selector.id = 'another-record';
  const failures = [];
  try { admitShape(receipt, circular); } catch (error) { failures.push(error); }
  await admit(receipt, changed).catch((error) => failures.push(error));
  assert.equal(failures.length, 2);
  for (const error of failures) {
    for (const [operation, stage] of [['generate', 'render'], ['session-save', 'unknown']]) {
      const model = normalizeUserFacingError(error, { operation, stage });
      assert.deepEqual([model.code, model.operation, model.stage, model.context, model.actions],
        ['ALIGNMENT_RESET_EVIDENCE', operation, 'request-validation', {}, ['edit-input']]);
      assert.equal(model.summary, 'What Reset alignment would restore no longer matches the current records'
        + ' or alignment plan. Undo the changes made after Align, or select Align\u2026 again.');
    }
  }
});
test('new Align replaces evidence relative to immediately preceding artifact', async () => {
  const a = fixture(); const first = await build(a); const b = structuredClone(a.after);
  b.renderRequest.records[1].presentation.reverseComplement = false;
  b.renderRequest.layout.recordTranslations[0].x += 33;
  const second = await build({before:a.after, after:b});
  assert.deepEqual(second.directions, [{recordKey:'b',before:true,after:false}]);
  assert.deepEqual(second.referenceDeltaX, {recordKey:'a',deltaX:33});
  assert.notDeepEqual(second, first);
});
test('Python and Web admit the same source-bound compact receipt and reject corrupted bindings', async () => {
  const f = fixture(); const receipt = await build(f);
  const session = {...f.after, editorState:{alignmentResetReceipt:receipt}};
  const script = `import json,sys; from gbdraw.session_io import _validate_alignment_reset_receipt, split_resource_bytes; document, sources = split_resource_bytes(json.load(sys.stdin)); _validate_alignment_reset_receipt(document, sources)`;
  const valid = spawnSync(process.env.PYTHON || 'python', ['-c', script], {input:JSON.stringify(session),encoding:'utf8'});
  if (valid.error) throw valid.error;
  assert.equal(valid.status, 0, valid.stderr);
  session.editorState.alignmentResetReceipt.binding = '0'.repeat(64);
  const invalid = spawnSync(process.env.PYTHON || 'python', ['-c', script], {input:JSON.stringify(session),encoding:'utf8'});
  if (invalid.error) throw invalid.error;
  assert.notEqual(invalid.status, 0); assert.match(invalid.stderr, /binding changed/);
});

test('Python and Web bind supplementary Unicode record keys in the same deterministic order', async () => {
  const f = fixture();
  const rename = value => {
    if (Array.isArray(value)) return value.forEach(rename);
    if (!value || typeof value !== 'object') return;
    if (value.recordKey) value.recordKey = ({a:'record-\u{10000}',b:'record-\uE000',c:'record-c'})[value.recordKey];
    Object.values(value).forEach(rename);
  };
  rename(f.before); rename(f.after);
  const receipt = await build(f);
  const script = 'import json,sys; from gbdraw.session_io import _validate_alignment_reset_receipt, split_resource_bytes; document, sources = split_resource_bytes(json.load(sys.stdin)); _validate_alignment_reset_receipt(document, sources)';
  const result = spawnSync(process.env.PYTHON || 'python', ['-c',script], {
    input:JSON.stringify({...f.after,editorState:{alignmentResetReceipt:receipt}}),encoding:'utf8'
  });
  if (result.error) throw result.error;
  assert.equal(result.status,0,result.stderr);
});

// OV-2772: v0.14.5 and dev (request schema 9) wrote a cropped Linear card's record selector in its region and bound
// the record's own selector, `null`; Load promotes that selector onto the record. A receipt binds the request as the
// file stored it, as Python's verifier reads it, and Load rebinds it to the promoted request. The bindings are the
// ones the v0.14.5 Web writes for these requests (`buildSimilarityAlignmentResetReceipt` at tag v0.14.5).
const storedBeforeWindows = (placement) => ({ schema: 9, mode: 'linear', records: ['a', 'b', 'c'].map(recordKey => ({
  recordKey, source: { kind: 'genbank', resourceId: 'source' }, cardinality: 'exactly_one',
  selector: recordKey === 'b' && placement === 'region' ? null : { kind: 'recordId', value: recordKey },
  region: recordKey === 'b'
    ? { selector: placement === 'region' ? { kind: 'recordId', value: 'b' } : null, start: 2, end: 10, reverseComplement: false }
    : null,
  display: { isCircular: null, startCoordinate: null }, presentation: { reverseComplement: false }
})), layout: { recordTranslations: ['a', 'b', 'c'].map(recordKey => ({ recordKey, x: 7, y: -3 })),
  similarityAlignment: fixture().after.renderRequest.layout.similarityAlignment } });
for (const [placement, binding] of [
  ['region', 'fe2664f3a1b0107bdf2e54ec4c0636a064039fb4484a758fa8109c0648964da4'],
  ['record', 'cf9c72f1cdbb85e4b89178634a5a22c54a7a6250e5357b985a457c2b3d8317b1']
]) test(`a receipt stored with request schema 9 binds the selector stored on the ${placement}`, async () => {
  const stored = storedBeforeWindows(placement);
  const canonical = { resources: { source: resource },
    renderRequest: promoteCanonicalRenderRequestToCurrent(structuredClone(stored)) };
  assert.deepEqual(canonical.renderRequest.records[1].selector, { kind: 'recordId', value: 'b' });
  const receipt = { binding, directions: [], referenceDeltaX: null };
  const read = await readStored(receipt, canonical, stored);
  assert.equal(read.stale, false);
  // Rebound to the promoted request: History, Align and Save read it as current.
  assert.equal(await admit(read.value, canonical), read.value);
  assert.deepEqual({ ...read.value, binding }, receipt);
  // A crop changed after the save still makes it stale.
  const changed = structuredClone(canonical);
  changed.renderRequest.records[1].windows[0].end = 11;
  assert.deepEqual(await readStored(receipt, changed, stored), { value: null, stale: true });
});
