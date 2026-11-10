// A Linear Generate with an editable comparison never writes a card that draws a file of several records whole
// (OV-3272): its prepare step asks the record selector to split every automatic card into record cards before the
// request is built (app/app-setup.js `prepareGenerate`, app/record-discovery-requests.js
// `materializeAutomaticLinearRecords`). Python numbers a comparison's records over the records such a card draws,
// while the Web's unit plan counts the card as one, so a request that carried one beside a comparison could bind the
// comparison to the wrong record; Run Info refuses its Source recipe (tests/web/run-info.test.mjs).
import assert from 'node:assert/strict';
import { readFile } from 'node:fs/promises';
import { join } from 'node:path';
import test from 'node:test';

globalThis.window = {
  Vue: {
    computed: (definition) => ({ get value() { return definition(); } }), ref: (value) => ({ value }), watch: () => {}
  }
};
const { createRecordDiscoveryRuns } = await import('../../gbdraw/web/js/app/record-discovery-requests.js');

test('Generate prepares an editable Linear comparison by splitting every automatic card into record cards', async () => {
  const expansions = new Map([['waiting', false]]);
  /** @type {[string, boolean][][]} */
  const refreshed = [];
  const { materializeAutomaticLinearRecords } = createRecordDiscoveryRuns({
    pendingDiscoveryOf: () => ({ expansions, circularUpload: null }),
    requestRecordListAfterUpload: () => {},
    listedDrawingOf: () => null,
    rebindAnnotationRecords: () => {},
    refreshLinearRecords: async () => {
      refreshed.push([...expansions]);
      return { status: 'ok' };
    }
  });
  const drawing = {
    mode: 'linear',
    linearSeqs: [
      { uid: 'whole', region_record_id: '' }, { uid: 'picked', region_record_id: 'NC_1' },
      { uid: 'waiting', region_record_id: '' }
    ]
  };
  await materializeAutomaticLinearRecords(drawing);
  // The automatic card is split by the refresh that follows; a picked record stays, and a card that already waits for
  // its discovery keeps what its writer asked.
  assert.deepEqual(refreshed, [[['waiting', false], ['whole', true]]]);
  // A Circular drawing has no Linear cards to split.
  await materializeAutomaticLinearRecords({ mode: 'circular', linearSeqs: drawing.linearSeqs });
  assert.equal(refreshed.length, 1);
});

test('Generate runs that split before it builds the request of a drawing whose comparison is editable', async () => {
  const source = await readFile(join(process.cwd(), 'gbdraw', 'web', 'js', 'app', 'app-setup.js'), 'utf8');
  const prepare = source.slice(source.indexOf('prepareGenerate: async () => {'));
  assert.match(prepare.slice(0, 600), /if \(drawing\.importedComparisonIntent\.disposition === IMPORTED_COMPARISON_DISPOSITIONS\.EDITABLE\) \{\s*await materializeAutomaticLinearRecords\(drawing\);/);
});
