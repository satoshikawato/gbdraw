// A Linear Generate with an editable comparison never writes a card that draws a file of several records whole
// (OV-3272): its prepare step asks the record selector to split every automatic card into record cards before the
// request is built (app/record-discovery-requests.js `prepareLinearGenerate`,
// `materializeAutomaticLinearRecords`). Python numbers a comparison's records over the records such a card draws,
// while the Web's unit plan counts the card as one, so a request that carried one beside a comparison could bind the
// comparison to the wrong record; Run Info refuses its Source recipe (tests/web/run-info.test.mjs).
import assert from 'node:assert/strict';
import test from 'node:test';

globalThis.window = {
  Vue: {
    computed: (definition) => ({ get value() { return definition(); } }), ref: (value) => ({ value }), watch: () => {}
  }
};
const { createRecordDiscoveryRuns } = await import('../../gbdraw/web/js/app/record-discovery-requests.js');
const { IMPORTED_COMPARISON_DISPOSITIONS } = await import('../../gbdraw/web/js/services/imported-comparison-intent.js');

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

test('Generate splits the automatic cards only when the comparison is editable, and reads the cards again otherwise', async () => {
  /** @type {Map<string, Map<string, boolean>>} */
  const expansionsOf = new Map();
  /** @type {[string, [string, boolean][]][]} */
  const refreshed = [];
  const { prepareLinearGenerate } = createRecordDiscoveryRuns({
    pendingDiscoveryOf: (drawing) => {
      if (!expansionsOf.has(drawing.id)) expansionsOf.set(drawing.id, new Map());
      return { expansions: expansionsOf.get(drawing.id), circularUpload: null };
    },
    requestRecordListAfterUpload: () => {},
    listedDrawingOf: () => null,
    rebindAnnotationRecords: () => {},
    refreshLinearRecords: async (drawing) => {
      refreshed.push([drawing.id, [...expansionsOf.get(drawing.id) || []]]);
      return { status: 'ok' };
    }
  });
  const drawing = (id, disposition) => ({
    id, mode: 'linear', importedComparisonIntent: { disposition }, linearSeqs: [{ uid: 'whole', region_record_id: '' }]
  });
  await prepareLinearGenerate(drawing('editable', IMPORTED_COMPARISON_DISPOSITIONS.EDITABLE));
  const others = Object.values(IMPORTED_COMPARISON_DISPOSITIONS).filter(
    (disposition) => disposition !== IMPORTED_COMPARISON_DISPOSITIONS.EDITABLE
  );
  for (const disposition of others) await prepareLinearGenerate(drawing(disposition, disposition));
  assert.ok(others.length > 0);
  assert.deepEqual(refreshed, [
    ['editable', [['whole', true]]],
    ...others.map((disposition) => [disposition, []])
  ]);
});
