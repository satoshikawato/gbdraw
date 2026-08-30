import assert from 'node:assert/strict';

import {
  createLinearTypographyController,
  reconcileImportedLinearTypographyLink
} from '../../gbdraw/web/js/app/linear-typography.js';
import { createAutoValueDisplay } from '../../gbdraw/web/js/app/auto-value-display.js';

const autoValues = createAutoValueDisplay({
  mode: { value: 'linear' },
  form: { scale_style: 'ruler' },
  adv: {}
});
assert.equal(autoValues.autoValueText('scaleFontSize'), '24/16 px (s/l auto)');
assert.equal(autoValues.autoValueText('rulerLabelFontSize'), '20/12 px (s/l auto)');

const linked = { value: true };
const adv = { scale_font_size: null, ruler_label_font_size: null };
const typography = createLinearTypographyController({ adv, linked });

assert.equal(linked.value, true);
typography.setScaleFontSize('18');
assert.deepEqual(adv, { scale_font_size: 18, ruler_label_font_size: 18 });

assert.equal(typography.setLinked(false), true);
assert.deepEqual(adv, { scale_font_size: 18, ruler_label_font_size: 18 });
typography.setRulerLabelFontSize('12');
assert.deepEqual(adv, { scale_font_size: 18, ruler_label_font_size: 12 });

assert.equal(typography.setLinked(true), true);
assert.deepEqual(adv, { scale_font_size: 18, ruler_label_font_size: 18 });
assert.equal(typography.setLinked(true), false);
assert.deepEqual(adv, { scale_font_size: 18, ruler_label_font_size: 18 });
typography.setScaleFontSize('21');
assert.deepEqual(adv, { scale_font_size: 21, ruler_label_font_size: 21 });
assert.equal(typography.setRulerLabelFontSize('9'), false);
assert.deepEqual(adv, { scale_font_size: 21, ruler_label_font_size: 21 });

const importedUnequal = { scale_font_size: 17, ruler_label_font_size: 29 };
const importedUnequalLink = { value: true };
assert.equal(reconcileImportedLinearTypographyLink({
  adv: importedUnequal,
  linked: importedUnequalLink,
  ui: { linearTypographyLinked: true }
}), false);
assert.deepEqual(importedUnequal, { scale_font_size: 17, ruler_label_font_size: 29 });

const importedEqual = { scale_font_size: 14, ruler_label_font_size: 14 };
for (const [ui, expected] of [
  [{}, false],
  [{ linearTypographyLinked: false }, false],
  [{ linearTypographyLinked: true }, true]
]) {
  const importedEqualLink = { value: true };
  assert.equal(reconcileImportedLinearTypographyLink({
    adv: importedEqual,
    linked: importedEqualLink,
    ui
  }), expected);
  assert.deepEqual(importedEqual, { scale_font_size: 14, ruler_label_font_size: 14 });
}
