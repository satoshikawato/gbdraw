// Loading indicator (Owner decision 2026-10-11): every loading place shows the
// gbdraw logo drawn in one stroke. The loader reveals the logo's own filled
// paths through two mask pens, so its finished frame is the logo itself. The
// boot splash holds a static copy of the Vue template because it renders
// before Vue mounts; this test keeps the two copies identical.
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { join } from 'node:path';
import test from 'node:test';

const webRoot = join(process.cwd(), 'gbdraw', 'web');
const html = readFileSync(join(webRoot, 'index.html'), 'utf8');
const logo = readFileSync(join(webRoot, 'assets', 'gbdraw-logo.svg'), 'utf8');
const outlinedTitle = readFileSync(join(webRoot, 'assets', 'gbdraw-logo-title-outlined.svg'), 'utf8');

const svgAt = (source, label) => {
  const start = source.indexOf('<svg class="gbdraw-loader"');
  assert.ok(start >= 0, `${label} has no loader SVG`);
  return source.slice(start, source.indexOf('</svg>', start) + '</svg>'.length);
};
const templateMatch = html.match(/<script type="text\/x-template" id="gbdraw-loader-template">([\s\S]*?)<\/script>/);
const bootMatch = html.match(/<div id="app-boot-splash"[\s\S]*?<\/noscript>/);

test('no loading place keeps a rotating spinner', () => {
  for (const pattern of [/animate-spin/, /ph-spinner/, /ph-circle-notch/, /boot-spinner/, /app-boot-spin/]) {
    assert.doesNotMatch(html, pattern);
  }
});

test('the boot splash copy equals the Vue template apart from its fixed mask ids', () => {
  assert.ok(templateMatch, 'missing #gbdraw-loader-template');
  assert.ok(bootMatch, 'missing #app-boot-splash');
  const normalized = svgAt(templateMatch[1], 'the template')
    .replace(/:id="ids\.([ab])"/g, 'id="gbdraw-loader-boot-$1"')
    .replace(/:mask="`url\(#\$\{ids\.([ab])\}\)`"/g, 'mask="url(#gbdraw-loader-boot-$1)"')
    .replace(/\s+/g, ' ');
  assert.equal(svgAt(bootMatch[0], 'the boot splash').replace(/\s+/g, ' '), normalized);
});

test('the loader reveals the logo asset paths under the logo transform', () => {
  const svg = svgAt(templateMatch[1], 'the template');
  const transform = logo.match(/id="g1"\s+transform="([^"]+)"/)[1];
  const logoPaths = [...logo.matchAll(/\sd="(m 67\.732299[^"]+)"/g)].map((match) => match[1]).sort();
  assert.equal(logoPaths.length, 2);
  const revealed = [...svg.matchAll(/<g :mask="[^"]+"><path transform="([^"]+)" d="([^"]+)"\/><\/g>/g)];
  assert.deepEqual(revealed.map((match) => match[1]), [transform, transform]);
  assert.deepEqual(revealed.map((match) => match[2]).sort(), logoPaths);
  assert.match(svg, /aria-hidden="true"/);
});

// The boot splash shows the outlined title logo's word beside the loader, so its frame
// starts with the loader's viewBox, and its mark is the logo asset's paths.
test('the outlined title logo starts with the loader frame and draws the logo asset paths', () => {
  const loaderViewBox = svgAt(templateMatch[1], 'the template').match(/viewBox="([^"]+)"/)[1].split(' ').map(Number);
  const titleViewBox = outlinedTitle.match(/viewBox="([^"]+)"/)[1].split(' ').map(Number);
  assert.deepEqual([titleViewBox[0], titleViewBox[1], titleViewBox[3]], [loaderViewBox[0], loaderViewBox[1], loaderViewBox[3]]);
  const [, transform, d] = outlinedTitle.match(/<path transform="([^"]+)" d="([^"]+)"\/>/);
  assert.equal(transform, logo.match(/id="g1"\s+transform="([^"]+)"/)[1]);
  const logoPaths = [...logo.matchAll(/\sd="m (67\.732299[^"]+)"/g)].map((match) => `M ${match[1]}`);
  assert.equal(d, logoPaths.join(' '));
});

test('every loading place renders the loader', () => {
  const places = [
    /<gbdraw-loader class="w-16 h-16 mb-6"><\/gbdraw-loader>\s*<h3[^>]*>Generating Diagram/,
    /<gbdraw-loader v-if="generationCancelRequested"><\/gbdraw-loader><i v-else class="ph ph-x-circle"/,
    /<gbdraw-loader v-if="sessionSavePending"><\/gbdraw-loader><i v-else class="ph ph-archive" aria-hidden="true"><\/i> \{\{ t\('Save Session'\) \}\}/,
    /<gbdraw-loader v-if="exampleSessionPending"><\/gbdraw-loader><i v-else class="ph ph-play-circle" aria-hidden="true"><\/i> \{\{ t\('Load example'\) \}\}/,
    /<gbdraw-loader v-if="exampleSessionPending"><\/gbdraw-loader><i v-else class="ph ph-play-circle" aria-hidden="true"><\/i> Load an example/,
    /<gbdraw-loader class="gbdraw-loader-on-fill text-xl"><\/gbdraw-loader> <span class="generate-bar-long">\{\{ generationCancelRequested \? 'Canceling\.\.\.' : 'Processing\.\.\.' \}\}/,
    /<gbdraw-loader><\/gbdraw-loader> Preparing features\.\.\./,
    /<gbdraw-loader><\/gbdraw-loader>\s*Reflowing label placement\.\.\./
  ];
  assert.deepEqual(places.filter((place) => !place.test(html)).map(String), []);
  assert.equal(html.match(/<gbdraw-loader\b/g).length, places.length);
});
