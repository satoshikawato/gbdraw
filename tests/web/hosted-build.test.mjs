// The hosted analytics flag: `tools/prepare_cloudflare_pages.py` writes it
// into the gbdraw.app bundle, `utils/hosted-build.js` reads it once at setup,
// and the About section shows its Google Analytics notice only when it is set.
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { join } from 'node:path';
import test from 'node:test';

import { readHostedBuild } from '../../gbdraw/web/js/utils/hosted-build.js';

const SELECTOR = 'meta[name="gbdraw-hosted-analytics"]';

/** A document whose only element is the flag, when `hosted` is true. */
const documentOf = (hosted) => ({
  querySelector: (selector) => (hosted && selector === SELECTOR ? {} : null)
});

test('the flag is read from the bundle the build wrote', () => {
  assert.deepEqual(readHostedBuild(documentOf(true)), { googleAnalytics: true });
  assert.deepEqual(readHostedBuild(documentOf(false)), { googleAnalytics: false });
});

test('the About section shows the analytics notice only when the flag is set', () => {
  const html = readFileSync(join(process.cwd(), 'gbdraw', 'web', 'index.html'), 'utf8');
  const title = html.indexOf("t('Hosted Site Analytics')");
  assert.ok(title > 0, 'index.html holds the notice');
  assert.equal(html.indexOf("t('Hosted Site Analytics')", title + 1), -1, 'one notice');
  const open = html.lastIndexOf('<div', title);
  assert.match(html.slice(open, html.indexOf('>', open) + 1), /^<div v-if="hostedBuild\.googleAnalytics" /);
  assert.ok(open > html.indexOf("t('Privacy & Security')"), 'after the Privacy & Security box');
});
