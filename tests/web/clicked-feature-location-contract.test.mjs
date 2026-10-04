import assert from 'node:assert/strict';
import fs from 'node:fs';
import test from 'node:test';

const read = (path) => fs.readFileSync(new URL(`../../gbdraw/web/js/${path}`, import.meta.url), 'utf8');

test('clickedFeatureLocation shows only the location the popup payload carries', () => {
  const source = read('app/app-setup.js');
  const start = source.indexOf('const clickedFeatureLocation = computed(');
  assert.notEqual(start, -1);
  const body = source.slice(start, source.indexOf('\n\n', start));
  assert.match(body, /clickedFeature\.value\?\.location \|\| ''/);
  // No envelope fallback that rebuilds start+1..end from cf.feat.
  assert.doesNotMatch(body, /\.feat\b|\.start\b|\.end\b|\.strand\b|\+ 1/);
});

test('every clickedFeature payload writer sets location with the shared formatter', () => {
  const actions = read('app/feature-editor/svg-actions.js');
  const builder = actions.slice(
    actions.indexOf('const buildClickedFeaturePayload ='),
    actions.indexOf('clickedFeature.value = buildClickedFeaturePayload(')
  );
  assert.match(builder, /const locationText = formatFeatureLocation\(feat\);/);
  assert.match(builder, /\blocation: locationText,/);
  // The popup payload is assigned from the builder only: no other production
  // module writes a non-null object to clickedFeature.value.
  const writers = [];
  for (const path of [
    'app/app-setup.js', 'app/feature-editor/color-actions.js', 'app/feature-editor/label-actions.js',
    'app/feature-editor/svg-actions.js', 'app/feature-editor/visibility-actions.js',
    'app/feature-search/preview-actions.js', 'app/feature-selection.js', 'app/orthogroups.js',
    'app/ui.js', 'app/watchers.js', 'services/reset.js', 'services/history-snapshot.js'
  ]) {
    for (const match of read(path).matchAll(/clickedFeature\.value = ([^;]+);/g)) {
      if (match[1].trim() !== 'null') writers.push(`${path}: ${match[1].trim()}`);
    }
  }
  assert.deepEqual(writers, ['app/feature-editor/svg-actions.js: buildClickedFeaturePayload(feat, featureElement, renderedSvgId)']);
});
