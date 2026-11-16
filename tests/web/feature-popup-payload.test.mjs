// The feature popup's payload (services/feature-popup.js) keeps the catalog keys of its texts, and the popup words them
// where it shows them, so a UI language change while the popup is open rewords them (OV-3676).
import assert from 'node:assert/strict';
import { test } from 'node:test';
import { buildFeaturePopup } from '../../gbdraw/web/js/services/feature-popup.js';
import { SOURCE_LOCALE, setUiLocale, tm } from '../../gbdraw/web/js/utils/i18n.js';

const feat = {
  id: 'rec-1\u0000kinase', recordKey: 'rec-1', biologicalFeatureId: 'kinase', record_id: 'REC', type: 'CDS',
  start: 0, end: 6, strand: 1, qualifiers: { product: ['kinase'] }
};
const view = {
  svgId: 'svg-kinase', legendName: 'CDS', color: '#000000', strokeColor: '#000000', strokeWidth: 0.5,
  featureVisibility: 'auto', identityEditable: true
};

test('a popup built in one language shows its Details rows and label reason in the language shown later', async () => {
  await setUiLocale('ja', async () => ({ default: {
    'Feature type': 'フィーチャータイプ', Length: '長さ',
    'No editable feature label for this feature in current diagram.': '現在の図には、このフィーチャーの編集可能なラベルがありません。'
  } }));
  let popup;
  try {
    popup = buildFeaturePopup(feat, view);
  } finally {
    await setUiLocale(SOURCE_LOCALE);
  }
  // The popup and Copy details word a row's label with `tm`.
  assert.deepEqual(popup.detailRows.map((row) => [row.key, tm(row.label), row.value]), [
    ['type', 'Feature type', 'CDS'],
    ['length', 'Length', '6 bp']
  ]);
  assert.equal(tm(popup.labelUnavailableReason), 'No editable feature label for this feature in current diagram.');
});
