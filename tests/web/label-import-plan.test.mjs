import assert from 'node:assert/strict';
import test from 'node:test';

import {
  labelImportNotice,
  labelIntentSignature,
  labelTableExportLines,
  planLabelImport
} from '../../gbdraw/web/js/services/label-import-plan.js';
import { parseLabelOverrideTsv } from '../../gbdraw/web/js/services/label-override-table.js';
import { normalizeTsvCell } from '../../gbdraw/web/js/utils/tsv-cell.js';

// The pure middles of Load and Export Label TSV (app/feature-editor/label-table.js).
const tracked = (featureId, sourceText) => ({
  featureId,
  feature: { svg_id: featureId, type: 'CDS', record_key: 'record-1', biological_feature_id: `bio-${featureId}` },
  sourceText
});
const untracked = (sourceText) => ({ featureId: '', feature: null, sourceText });

test('a global label row edits a label by its source text, any other row only a label with a feature key', () => {
  const rows = parseLabelOverrideTsv('*\t*\tlabel\t^alpha$\tBULK\n*\tCDS\tproduct\t.*\tOWN\n');
  const labels = [tracked('fa', 'alpha'), untracked('alpha'), untracked('beta'), tracked('fc', ''), tracked('fd', 'delta')];
  const plan = planLabelImport(labels, rows, [0, 0, 1, 0, -1]);
  assert.deepEqual(plan.applicable.map(({ entry, matchedRow, tracked: hasKey }) => [entry.sourceText, matchedRow.labelText, hasKey]), [
    ['alpha', 'BULK', true],
    ['alpha', 'BULK', false]
  ]);
  assert.equal(plan.skippedNonTrackableCount, 2, 'the product row without a feature key and the global row without a source text');
});

test('each outcome of a load names its row count and why it applied or not (OV-131, OV-447, OV-449)', () => {
  const plan = (applicable, skippedNonTrackableCount = 0) => ({ plan: { applicable, skippedNonTrackableCount } });
  const edit = { entry: tracked('fa', 'alpha'), matchedRow: parseLabelOverrideTsv('*\t*\tlabel\t^alpha$\tA\n')[0], tracked: true };
  assert.deepEqual(labelImportNotice(1, { declined: 'label-intent' }), {
    text: 'Loaded 1 row. Not applied: label text, label visibility, or feature visibility changed while the table loaded. '
      + 'The changes were kept. Load the table again.',
    tone: 'warning'
  });
  assert.deepEqual(labelImportNotice(2, { declined: 'diagram' }), {
    text: 'Loaded 2 rows. Not applied: the displayed diagram changed while the table loaded. Load the table again.', tone: 'warning'
  });
  assert.match(labelImportNotice(1, { declined: 'rule-inputs' }).text, /^Loaded 1 row\. Not applied: colors, color rules, the Legend/);
  assert.deepEqual(labelImportNotice(1, { declined: 'no-diagram' }), {
    text: 'Loaded 1 row. No diagram is shown, so nothing was applied.', tone: 'warning'
  });
  const undecided = labelImportNotice(1, { undecided: 2 });
  assert.equal(undecided.text, 'Loaded 1 row. Not applied.');
  assert.match(undecided.items[0], /^location and record_location rows cannot be matched to 2 labels/);
  assert.deepEqual(labelImportNotice(1, plan([])), {
    text: 'Loaded 1 row. Not applied; the existing label edits were kept.', items: ['No row matched a label of the diagram.']
  });
  assert.deepEqual(labelImportNotice(1, plan([], 1)).items, ['The 1 matched label lacked a feature key.']);
  assert.deepEqual(labelImportNotice(1, plan([edit])), { text: 'Loaded 1 row. Applied to 1 label.', items: [] });
  assert.deepEqual(labelImportNotice(3, plan([edit, edit], 1)), {
    text: 'Loaded 3 rows. Applied to 2 labels.', items: ['1 match lacked a feature key and was not applied.']
  });
});

test('the label intent signature follows the per-feature and the bulk label edits', () => {
  const drawing = { featureOverrides: { a: { labelText: 'A' } }, labelTextBulkOverrides: {} };
  const before = labelIntentSignature(drawing);
  assert.equal(labelIntentSignature(structuredClone(drawing)), before);
  assert.notEqual(labelIntentSignature({ ...drawing, labelTextBulkOverrides: { alpha: 'B' } }), before);
});

test('Export Label TSV writes the saved label table, or else the bulk label edits as sorted global label rows', () => {
  const saved = [{ recordId: '*', featureType: 'CDS', qualifier: 'product', valueRegex: '^x$', labelText: 'X' }];
  assert.deepEqual(labelTableExportLines(saved, { alpha: 'A' }), ['*\tCDS\tproduct\t^x$\tX']);
  assert.deepEqual(labelTableExportLines([], { 'b.c': 'B\tC', alpha: 'A', '': 'dropped' }), [
    '*\t*\tlabel\t^alpha$\tA',
    `*\t*\tlabel\t^b\\.c$\t${normalizeTsvCell('B\tC')}`
  ]);
  assert.deepEqual(labelTableExportLines(null, {}), []);
});
