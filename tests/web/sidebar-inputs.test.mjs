// Sidebar input cards: record readouts before the source records are read.
import assert from 'node:assert/strict';
import test from 'node:test';
import { cp, mkdtemp, writeFile } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { pathToFileURL } from 'node:url';

const tempRoot = await mkdtemp(join(tmpdir(), 'gbdraw-sidebar-inputs-'));
await cp(join(process.cwd(), 'gbdraw', 'web', 'js'), join(tempRoot, 'js'), { recursive: true });
await writeFile(join(tempRoot, 'package.json'), '{"type":"module"}', 'utf8');
const load = (path) => import(pathToFileURL(join(tempRoot, 'js', path)));

const { buildAnnotationRecordCatalog } = await load('app/annotations/record-catalog.js');
const { annotationRecordOptions, createAnnotationRecordSelector } = await load('app/annotations/record-selector.js');
const { coordinateTarget } = await load('app/annotations/target-actions.js');
const { createLinearRecordSelector } = await load('app/linear-record-selector.js');

const savedTarget = () => ({
  id: 'lsc',
  target: coordinateTarget({ recordId: 'NC_001879.2', start: 1, end: 86686 }),
  metadata: {}
});

// OV-373: a loaded Session reads no record bytes, so an annotation that names
// its record shows that record as not inspected, not as unavailable.
test('OV-373: a saved annotation target before Inspect is not called unavailable', () => {
  const pending = buildAnnotationRecordCatalog({
    mode: 'circular',
    circularSource: { sourceKey: 'circular', hasInput: true, status: 'deferred', records: [] }
  });
  const labels = annotationRecordOptions(pending, savedTarget()).map((option) => option.label);
  assert.equal(labels.some((label) => /unavailable/i.test(label)), false, labels.join(' | '));

  const selector = createAnnotationRecordSelector({ getCatalog: () => pending });
  assert.equal(selector.showsSelect(savedTarget()), false);
  assert.equal(selector.captionFor(savedTarget()), 'Record: NC_001879.2');

  // After Inspect the one record resolves; a single-record input needs no select.
  const ready = buildAnnotationRecordCatalog({
    mode: 'circular',
    circularSource: {
      sourceKey: 'circular', hasInput: true, status: 'ready',
      records: [{ selector: '#1', record_id: 'NC_001879.2', record_length: 155943 }]
    }
  });
  const readySelector = createAnnotationRecordSelector({ getCatalog: () => ready });
  assert.equal(readySelector.valueFor(savedTarget()), ready.records[0].key);
  assert.equal(readySelector.showsSelect(savedTarget()), false);
  assert.equal(readySelector.captionFor(savedTarget()), 'Record: NC_001879.2');

  // Several records keep the select.
  const multi = buildAnnotationRecordCatalog({
    mode: 'circular',
    circularSource: {
      sourceKey: 'circular', hasInput: true, status: 'ready',
      records: [{ selector: '#1', record_id: 'NC_001879.2' }, { selector: '#2', record_id: 'Other' }]
    }
  });
  assert.equal(createAnnotationRecordSelector({ getCatalog: () => multi }).showsSelect(savedTarget()), true);
});

// OV-379: before the records of a Linear File are read (a loaded Session
// defers discovery), the record select says "Not inspected", not "Loading".
test('OV-379: an unread Linear record select says Not inspected', async () => {
  const ref = (value) => ({ value });
  const state = {
    mode: ref('linear'),
    lInputType: ref('gb'),
    linearSeqs: [{ uid: 'row-a', gb: { name: 'a.gb' }, fasta: null, region_record_id: '' }]
  };
  let resolveRead;
  const controller = createLinearRecordSelector({
    state,
    reactive: (value) => value,
    recordReader: () => new Promise((resolve) => { resolveRead = resolve; }),
    logger: { warn: () => {} }
  });
  assert.deepEqual(controller.optionsFor(state.linearSeqs[0]), [
    { value: '', label: 'Not inspected', synthetic: false }
  ]);
  // A read in progress still says so.
  const refresh = controller.refresh();
  assert.equal(controller.optionsFor(state.linearSeqs[0])[0].label, 'Loading records...');
  resolveRead([{ selector: '#1', recordId: 'Rec', recordLength: 9 }]);
  await refresh;
  assert.equal(controller.statusFor(state.linearSeqs[0]), 'ready');
});
