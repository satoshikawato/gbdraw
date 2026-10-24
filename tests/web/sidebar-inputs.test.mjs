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
const { createAnnotationRecordSelector } = await load('app/annotations/record-selector.js');
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
  const selector = createAnnotationRecordSelector({ getCatalog: () => pending });
  assert.equal(selector.showsSelect(savedTarget()), false);
  assert.equal(selector.captionFor(savedTarget()), 'Record: NC_001879.2');

  // A catalog that failed to read its records says so beside the name.
  const failed = { status: 'error', records: [], issues: [], requiresSelection: false };
  assert.equal(createAnnotationRecordSelector({ getCatalog: () => failed }).captionFor(savedTarget()),
    'Record: NC_001879.2 (record list unavailable)');

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
  const { withModeDrawings } = await import(pathToFileURL(join(process.cwd(), 'tests/web/helpers/drawing-state.mjs')));
  const ref = (value) => ({ value });
  const rows = [{ uid: 'row-a', gb: { name: 'a.gb' }, fasta: null, region_record_id: '' }];
  const state = withModeDrawings({ mode: ref('linear'), lInputType: ref('gb') }, {
    circular: {}, linear: { id: 'linear', mode: 'linear', linearSeqs: rows }
  });
  let resolveRead;
  const controller = createLinearRecordSelector({
    state,
    reactive: (value) => value,
    recordReader: () => new Promise((resolve) => { resolveRead = resolve; }),
    logger: { warn: () => {} }
  });
  assert.deepEqual(controller.optionsFor(rows[0]), [
    { value: '', label: 'Not inspected', synthetic: false }
  ]);
  // A read in progress still says so.
  const refresh = controller.refresh();
  assert.equal(controller.optionsFor(rows[0])[0].label, 'Loading records...');
  resolveRead([{ selector: '#1', recordId: 'Rec', recordLength: 9 }]);
  await refresh;
  assert.equal(controller.statusFor(rows[0]), 'ready');
});

// Auto window and step values name their rule until the genome length is
// known, then show the value for the longest record the next Generate draws.
test('Auto window values follow the known Linear record lengths', async () => {
  const { createAutoValueDisplay } = await load('app/auto-value-display.js');
  const { withDrawings } = await import(pathToFileURL(join(process.cwd(), 'tests/web/helpers/drawing-state.mjs')));
  const linear = () => withDrawings({ mode: { value: 'linear' }, form: {}, adv: {} });
  const unknown = createAutoValueDisplay(linear());
  assert.equal(unknown.autoValueText('dinucleotideWindow'), 'Auto (by genome size)');
  assert.equal(unknown.autoValueText('depthWindow'), 'Auto (by genome size)');
  const known = createAutoValueDisplay(linear(), { linearRecordLengths: () => [40_000, 2_500_000] });
  assert.equal(known.autoValueText('dinucleotideWindow'), '10 kbp (auto)');
  assert.equal(known.autoValueText('dinucleotideStep'), '1 kbp (auto)');
  assert.equal(known.autoValueText('depthWindow'), '1 kbp (auto)');
});
