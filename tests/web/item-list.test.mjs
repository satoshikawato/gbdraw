// The one owner of the separator between shown items (utils/item-list.js): English and every left-to-right language
// keep the bytes of `join(', ')`; Persian and Arabic write the Arabic comma, each item isolated and the list read right
// to left.
import assert from 'node:assert/strict';
import test from 'node:test';

import { SOURCE_LOCALE, UI_LOCALES, setUiLocale } from '../../gbdraw/web/js/utils/i18n.js';
import { itemListText } from '../../gbdraw/web/js/utils/item-list.js';

const RLI = '⁧';
const FSI = '⁨';
const PDI = '⁩';

test('English shows the items as join(", ") does', () => {
  for (const items of [[], ['NC_000913'], ['NC_000913', 'NC_002695'], ['contig start', 'contig end', 'x']]) {
    assert.equal(itemListText(items), items.join(', '), JSON.stringify(items));
  }
});

test('every left-to-right UI language keeps ", "; Persian and Arabic write ، with isolated items read right to left', async () => {
  try {
    for (const { tag } of UI_LOCALES) {
      await setUiLocale(tag, async () => ({ default: {} }));
      const shown = itemListText(['NC_1', 'gene (x)', 'شروع کانتیگ']);
      if (tag === 'fa' || tag === 'ar') {
        assert.equal(shown, `${RLI}${FSI}NC_1${PDI}، ${FSI}gene (x)${PDI}، ${FSI}شروع کانتیگ${PDI}${PDI}`, tag);
        assert.equal(itemListText([]), '', tag);
        assert.equal(itemListText(['NC_1']), `${RLI}${FSI}NC_1${PDI}${PDI}`, tag);
      } else {
        assert.equal(shown, 'NC_1, gene (x), شروع کانتیگ', tag);
      }
    }
  } finally {
    await setUiLocale(SOURCE_LOCALE);
  }
});
