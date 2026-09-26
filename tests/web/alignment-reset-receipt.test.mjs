import assert from 'node:assert/strict';
import test from 'node:test';
import { spawnSync } from 'node:child_process';
const { buildSimilarityAlignmentResetReceipt: build, validateSimilarityAlignmentResetReceipt: admit } = await import('../../gbdraw/web/js/services/session-active-config-contract.js');
const source = Buffer.from('source bytes remain biological authority');
const resource = { name: 'record.gbk', kind: 'genbank', encoding: 'base64', data: source.toString('base64'), size: source.length };
const anchor = (recordKey) => ({ recordKey, biologicalFeatureId: `gene-${recordKey}`, sourceFeatureIndex: 1, stableFeatureSvgId: null });
const fixture = () => {
  const before = { resources: { source: resource }, renderRequest: {
    schema: 8, mode: 'linear', records: ['a', 'b', 'c'].map(recordKey => ({
      recordKey, source: {kind: 'genbank', resourceId: 'source'}, cardinality: 'exactly_one',
      selector: { kind: 'recordId', id: recordKey }, region: null,
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
  ['crop', f => { f.after.renderRequest.records[1].region = {start:2,end:10,reverseComplement:false}; }],
  ['plan anchor', f => { f.after.renderRequest.layout.similarityAlignment.records[1].anchor.biologicalFeatureId = 'another-gene'; }],
  ['plan removal', f => { f.after.renderRequest.layout.similarityAlignment = null; }]
]) test(`reject changed ${name} without dropping evidence`, async () => {
  const f = fixture(); const receipt = await build(f); mutate(f);
  await assert.rejects(admit(receipt, f.after), /receipt/);
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
  for (const invalid of variants) await assert.rejects(admit(invalid, f.after), /receipt/);
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
  const script = `import json,sys; from gbdraw.session_io import _validate_alignment_reset_receipt; _validate_alignment_reset_receipt(json.load(sys.stdin))`;
  const valid = spawnSync(process.env.PYTHON || 'python', ['-c', script], {input:JSON.stringify(session),encoding:'utf8'});
  if (valid.error) throw valid.error;
  assert.equal(valid.status, 0, valid.stderr);
  session.editorState.alignmentResetReceipt.binding = '0'.repeat(64);
  const invalid = spawnSync(process.env.PYTHON || 'python', ['-c', script], {input:JSON.stringify(session),encoding:'utf8'});
  if (invalid.error) throw invalid.error;
  assert.notEqual(invalid.status, 0); assert.match(invalid.stderr, /binding changed/);
});
