import assert from 'node:assert/strict';
import { test } from 'node:test';

import {
  anchorsAgree, classifyPlanAgainstCatalog, resetSourceMembers, restoredAxisTranslations, restoredPlacement
} from '../../gbdraw/web/js/services/alignment-plan-facts.js';

const anchor = (recordKey, biologicalFeatureId, sourceFeatureIndex = null, stableFeatureSvgId = null) => (
  { recordKey, biologicalFeatureId, sourceFeatureIndex, stableFeatureSvgId });
const plan = () => ({
  schema: 2, groupId: 'og-1', reference: anchor('a', 'ref', 1),
  records: [
    { recordKey: 'a', status: 'reference', rationale: 'reference', anchor: anchor('a', 'ref', 1) },
    { recordKey: 'b', status: 'aligned', rationale: 'user_selected', anchor: anchor('b', 'tb', 2) },
    { recordKey: 'c', status: 'aligned', rationale: 'user_selected', anchor: anchor('c', 'tc', 3) },
    { recordKey: 'd', status: 'skipped', rationale: 'skipped_by_user', anchor: null }
  ]
});
const feature = ({ recordKey, biologicalFeatureId, sourceFeatureIndex }) => ({ recordKey, biologicalFeatureId, sourceFeatureIndex });
const catalogOf = (...anchors) => anchors.map(feature);
const groupOf = (...anchors) => anchors;

test('anchors agree on record and feature unless an index or stable id both carry differs', () => {
  assert.equal(anchorsAgree(anchor('a', 'x', 1, 's'), anchor('a', 'x', null, null)), true);
  assert.equal(anchorsAgree(anchor('a', 'x', 1), anchor('a', 'x', 2)), false);
  assert.equal(anchorsAgree(anchor('a', 'x', 1, 's'), anchor('a', 'x', 1, 't')), false);
  assert.equal(anchorsAgree(anchor('a', 'x'), anchor('b', 'x')), false);
  assert.equal(anchorsAgree(null, anchor('a', 'x')), false);
});

test('the check before Generate classifies a plan by its catalog and its group', () => {
  const p = plan();
  const all = [p.reference, p.records[1].anchor, p.records[2].anchor];
  const classify = (catalogFeatures, groupAnchors) => classifyPlanAgainstCatalog({ plan: p, catalogFeatures, groupAnchors });
  assert.deepEqual(classify(catalogOf(...all), groupOf(...all)), { kind: 'resolve', goneRecordKeys: new Set() });
  // A Result without a catalog is not checked here.
  assert.deepEqual(classify([], groupOf(...all)), { kind: 'resolve', goneRecordKeys: new Set() });
  assert.deepEqual(classify(catalogOf(all[0], all[2]), groupOf(...all)), { kind: 'resolve', goneRecordKeys: new Set(['b']) });
  assert.deepEqual(classify(catalogOf(all[1], all[2]), groupOf(...all)), { kind: 'stale-reference',
    reason: 'The alignment reference feature is no longer found in its record.' });
  // Regrouped: the group no longer holds target c, every anchor still resolves.
  assert.deepEqual(classify(catalogOf(...all), groupOf(all[0], all[1])), { kind: 'current' });
  // Regrouped, and target b is gone from the catalog.
  assert.deepEqual(classify(catalogOf(all[0], all[2]), groupOf(all[0])), { kind: 'stale-reference',
    reason: 'An aligned target feature is no longer found in its record.' });
  // Two features match the reference: it no longer names one.
  assert.equal(classify(catalogOf(...all, all[0]), groupOf(...all)).kind, 'stale-reference');
  assert.deepEqual(classify(catalogOf(...all), groupOf(...all, all[0])), { kind: 'current' });
});

test('a Reset finds the catalog feature of each anchored record', () => {
  const p = plan();
  const catalogFeatures = catalogOf(p.reference, p.records[1].anchor, p.records[2].anchor).map((row) => ({ ...row, product: 'x' }));
  assert.deepEqual(resetSourceMembers(p, catalogFeatures).map(({ recordKey, product, stableFeatureSvgId }) => [recordKey, product, stableFeatureSvgId]),
    [['a', 'x', null], ['b', 'x', null], ['c', 'x', null]]);
  assert.throws(() => resetSourceMembers(p, catalogFeatures.slice(1)), /source anchor is unavailable or ambiguous/);
});

test('Reset restores the reference shift and, by scope, the directions before the latest Align', () => {
  const record = (recordKey, reverseComplement) => ({ recordKey, windows: [], presentation: { reverseComplement } });
  const records = [record('a', false), record('b', true)];
  const translations = [{ recordKey: 'a', x: 10, y: 0 }, { recordKey: 'b', x: 5, y: 40 }];
  const receipt = { binding: 'x', referenceDeltaX: { recordKey: 'a', deltaX: 4 },
    directions: [{ recordKey: 'b', before: false, after: true }] };
  const positions = restoredPlacement({ translations, records, receipt, scope: 'positions' });
  assert.deepEqual(positions, { translations: [{ recordKey: 'a', x: 6, y: 0 }, { recordKey: 'b', x: 5, y: 40 }],
    orientations: [{ recordKey: 'a', reverseComplement: false }, { recordKey: 'b', reverseComplement: true }],
    directionsChange: false });
  assert.equal(translations[0].x, 10, 'the committed translations are not changed');
  const both = restoredPlacement({ translations, records, receipt, scope: 'positions-and-directions' });
  assert.deepEqual(both.orientations, [{ recordKey: 'a', reverseComplement: false }, { recordKey: 'b', reverseComplement: false }]);
  assert.equal(both.directionsChange, true);
  assert.deepEqual(restoredPlacement({ translations, records, receipt: null, scope: 'positions' }).translations, translations);
  assert.throws(() => restoredPlacement({ translations: translations.slice(1), records, receipt, scope: 'positions' }),
    /reference placement is unavailable/);
});

test('Reset keeps each record axis where it was drawn once the helper projected the directions', () => {
  const translations = [{ recordKey: 'a', x: 6, y: 0 }, { recordKey: 'b', x: 5, y: 40 }];
  const vector = { a: false, b: false };
  const projection = { geometryOrientations: { a: false, b: false }, records: [
    { recordKey: 'a', beforeAxisY: 12, variants: [{ axisY: 12 }, { axisY: 20 }] },
    { recordKey: 'b', beforeAxisY: 30, variants: [{ axisY: 25 }, { axisY: 30 }] }
  ] };
  assert.deepEqual(restoredAxisTranslations({ translations, vector, projection }),
    [{ recordKey: 'a', x: 6, y: 0 }, { recordKey: 'b', x: 5, y: 45 }]);
  assert.throws(() => restoredAxisTranslations({ translations, vector, projection: null }), /geometry is not validated/);
  assert.throws(() => restoredAxisTranslations({ translations, vector: { a: false, b: true }, projection }),
    /geometry is not validated/);
});
