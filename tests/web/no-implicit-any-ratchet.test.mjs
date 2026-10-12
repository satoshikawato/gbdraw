// R14 "Typed boundaries" (gbdraw/web/CLAUDE.md): the noImplicitAny ratchet.
// `noImplicitAny` is not yet on in the guard config. Until it is,
// tests/web/types/tsconfig.no-implicit-any.json is the guard config
// (tests/web/types/tsconfig.json) plus `noImplicitAny`, and
// NO_IMPLICIT_ANY_BASELINE records how many diagnostics one `tsc` run with it
// reports per module under gbdraw/web/js/ (paths relative to it). A module
// without an entry, including a new module, has none. The guard asserts:
//   1. the ratchet config is exactly the guard config plus `noImplicitAny`;
//   2. every diagnostic of that run is in a module under gbdraw/web/js/;
//   3. no module reports more diagnostics than its entry.
//
// Lowering or removing an entry is optional: the test names the entries that
// may be lowered. Adding or raising an entry is an authority-only change:
// NO_IMPLICIT_ANY_BASELINE is a registered design-rule allowlist
// (tools/web-design-rule-guards.json, R14, kind `count-map`). When the map is
// empty, `noImplicitAny` is turned on in the guard config and this ratchet and
// its config are deleted.
import assert from 'node:assert/strict';
import { spawnSync } from 'node:child_process';
import { readFileSync } from 'node:fs';
import { createRequire } from 'node:module';
import { dirname, join } from 'node:path';
import test from 'node:test';
import { fileURLToPath } from 'node:url';

import { compareCountBaseline, reportLowerable } from './helpers/count-baseline.mjs';

const NO_IMPLICIT_ANY_BASELINE = {
  'app/annotations.js': 33,
  'app/annotations/record-catalog.js': 31,
  'app/annotations/record-selector.js': 30,
  'app/annotations/table-codec.js': 25,
  'app/annotations/target-actions.js': 6,
  'app/annotations/validation.js': 9,
  'app/app-setup.js': 343,
  'app/auto-value-display.js': 46,
  'app/candidate-render.js': 33,
  'app/circular-track-slots.js': 188,
  'app/circular-track-slots/measure-input.js': 4,
  'app/comparison-ui.js': 34,
  'app/depth-tracks.js': 4,
  'app/feature-editor/color-actions.js': 128,
  'app/feature-editor/feature-edit-table.js': 10,
  'app/feature-editor/label-actions.js': 147,
  'app/feature-editor/pattern-drafts.js': 12,
  'app/feature-editor/placement-actions.js': 22,
  'app/feature-editor/rule-actions.js': 61,
  'app/feature-editor/svg-actions.js': 89,
  'app/feature-editor/visibility-actions.js': 69,
  'app/feature-search/preview-actions.js': 6,
  'app/feature-search/preview-svg.js': 15,
  'app/feature-search/search-core.js': 129,
  'app/feature-selection.js': 52,
  'app/history-inputs.js': 15,
  'app/history-shortcuts.js': 2,
  'app/legend-layout.js': 11,
  'app/legend-layout/composition-actions.js': 133,
  'app/legend-layout/decoration-continuity.js': 28,
  'app/legend-layout/diagram-drag.js': 31,
  'app/legend-layout/reposition-actions.js': 4,
  'app/legend.js': 2,
  'app/legend/drag-actions.js': 7,
  'app/legend/entry-actions.js': 16,
  'app/legend/layout-actions.js': 1,
  'app/legend/sort-actions.js': 6,
  'app/legend/stroke-actions.js': 2,
  'app/legend/track-data-styles.js': 3,
  'app/linear-record-selector.js': 31,
  'app/linear-track-slots.js': 123,
  'app/linear-typography.js': 4,
  'app/losat-settings.js': 8,
  'app/orthogroups.js': 91,
  'app/pairwise-match-popup.js': 254,
  'app/palettes.js': 3,
  'app/preview-runtime.js': 24,
  'app/record-display-options.js': 75,
  'app/record-display/feature-anchor.js': 37,
  'app/record-display/feature-record-rotation.js': 45,
  'app/record-groups.js': 10,
  'app/result-paint-record.js': 2,
  'app/results.js': 4,
  'app/right-drawer.js': 1,
  'app/rule-matching.js': 36,
  'app/run-analysis.js': 304,
  'app/run-info.js': 222,
  'app/similarity-alignment.js': 268,
  'app/svg-styles.js': 18,
  'app/track-slot-colors.js': 5,
  'app/ui.js': 26,
  'app/watchers.js': 22,
  'services/config.js': 345,
  'services/feature-catalog.js': 80,
  'services/losat-cache.js': 59,
  'services/runtime-capabilities.js': 10,
  'services/session-authority.js': 19,
  'services/session-request.js': 484
};

const REPOSITORY_ROOT = join(dirname(fileURLToPath(import.meta.url)), '..', '..');
const RATCHET_CONFIG = 'tests/web/types/tsconfig.no-implicit-any.json';
const EXPECTED_RATCHET_CONFIG = { extends: './tsconfig.json', compilerOptions: { noImplicitAny: true } };
const DIAGNOSTIC = /^gbdraw\/web\/js\/(.+?)\((\d+),(\d+)\): error (TS\d+): /;

const tsc = () => {
  const require = createRequire(import.meta.url);
  let packagePath;
  try {
    packagePath = require.resolve('typescript/package.json');
  } catch {
    assert.fail('typescript is not installed: run `npm ci` (R14 runs the pinned devDependency)');
  }
  const pinned = JSON.parse(readFileSync(join(REPOSITORY_ROOT, 'package.json'), 'utf8')).devDependencies.typescript;
  const installed = JSON.parse(readFileSync(packagePath, 'utf8')).version;
  assert.equal(installed, pinned, `typescript ${installed} is installed but package.json pins ${pinned}: run \`npm ci\``);
  const run = spawnSync(
    process.execPath,
    [join(dirname(packagePath), 'bin', 'tsc'), '-p', RATCHET_CONFIG, '--pretty', 'false'],
    { cwd: REPOSITORY_ROOT, encoding: 'utf8', maxBuffer: 1 << 28 }
  );
  assert.ifError(run.error);
  return `${run.stdout}\n${run.stderr}`;
};

test('R14: the noImplicitAny config is the guard config plus noImplicitAny', () => {
  const config = JSON.parse(readFileSync(join(REPOSITORY_ROOT, RATCHET_CONFIG), 'utf8'));
  assert.deepEqual(
    config,
    EXPECTED_RATCHET_CONFIG,
    `${RATCHET_CONFIG} changes more than noImplicitAny; a config change is a rule change (authority-only)`
  );
});

test('R14: no module reports more noImplicitAny diagnostics than its NO_IMPLICIT_ANY_BASELINE entry', (t) => {
  Object.entries(NO_IMPLICIT_ANY_BASELINE).forEach(([path, count]) => {
    assert.ok(Number.isInteger(count) && count > 0, `NO_IMPLICIT_ANY_BASELINE['${path}'] must be a positive integer`);
  });
  const perModule = {};
  const unparsed = [];
  let inDiagnostic = false;
  tsc().split(/\r?\n/).forEach((line) => {
    if (!line.trim()) return;
    const match = line.match(DIAGNOSTIC);
    if (match) {
      perModule[match[1]] = (perModule[match[1]] || 0) + 1;
      inDiagnostic = true;
    } else if (/^\s/.test(line) && inDiagnostic) {
      // A continuation line of the diagnostic above.
    } else unparsed.push(line);
  });
  assert.deepEqual(unparsed, [], `tsc -p ${RATCHET_CONFIG} printed output outside gbdraw/web/js/`);

  const { above, lowerable } = compareCountBaseline(perModule, NO_IMPLICIT_ANY_BASELINE);
  const total = Object.values(perModule).reduce((sum, count) => sum + count, 0);
  t.diagnostic(`noImplicitAny: ${total} diagnostics in ${Object.keys(perModule).length} modules`);
  reportLowerable(t, 'NO_IMPLICIT_ANY_BASELINE', lowerable);
  assert.deepEqual(
    above.map(({ key, observed, recorded }) => `${key}: ${observed} noImplicitAny diagnostics, above its entry ${recorded}`),
    [],
    'Declare the types of the new parameters and variables (for example `@param {string} name`), '
      + 'or type the collection a callback iterates; adding or raising an entry is an authority-only change. '
      + `Run \`node node_modules/typescript/bin/tsc -p ${RATCHET_CONFIG} --pretty false\` to list the diagnostics.`
  );
});
