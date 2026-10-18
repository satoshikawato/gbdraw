// R14 "Typed boundaries" (gbdraw/web/CLAUDE.md): the explicit-any ratchet.
// Explicit `any` in JSDoc is type debt. EXPLICIT_ANY_BASELINE records it per
// module under gbdraw/web/js/ (paths relative to it); a module without an
// entry, including a new module, has none. The guard asserts that no module
// holds more than its entry: a module above its entry fails (type the value,
// use `unknown`, or narrow). Lowering or removing an entry is optional; the
// test names the entries that may be lowered.
//
// What counts: each JSDoc type node whose type, as the TypeScript checker of
// the R14 guard config (tests/web/types/tsconfig.json) resolves it, is
// intrinsic `any`: `any`, `*`, and each use of a type that resolves to `any`,
// such as a typedef alias of `any`, so an alias does not hide it (with
// `noImplicitAny` on, JSDoc `Object` is the `Object` type, not `any`). A node
// inside a counted node is not counted again; `Record<string, any>` counts
// once, for its `any` argument. `unknown`, `object`, and `Function` do not
// count.
// Every `.js` module is counted, generated modules included. Workers are
// outside the guard config; the checker resolves them in its inferred
// project, where JSDoc `Object` is also the `Object` type.
//
// EXPLICIT_ANY_BASELINE is a registered design-rule allowlist
// (tools/web-design-rule-guards.json, R14, kind `count-map`): adding an entry
// or raising one is an authority-only change. Owner decisions R15-5
// (2026-10-09) and BA-D04 (2026-10-10).
import assert from 'node:assert/strict';
import { readFileSync, readdirSync } from 'node:fs';
import { createRequire } from 'node:module';
import { dirname, join, sep } from 'node:path';
import test from 'node:test';
import { fileURLToPath } from 'node:url';

import { compareCountBaseline, reportLowerable } from './helpers/count-baseline.mjs';

const EXPLICIT_ANY_BASELINE = {
  'app/annotations.js': 8,
  'app/annotations/record-catalog.js': 3,
  'app/app-setup.js': 14,
  'app/candidate-render.js': 1,
  'app/circular-track-slots.js': 4,
  'app/feature-editor.js': 21,
  'app/feature-editor/color-actions.js': 40,
  'app/feature-editor/feature-edit-table.js': 7,
  'app/feature-editor/label-actions.js': 16,
  'app/feature-editor/placement-actions.js': 10,
  'app/feature-editor/rule-actions.js': 28,
  'app/feature-editor/svg-actions.js': 16,
  'app/feature-editor/visibility-actions.js': 5,
  'app/feature-search/preview-actions.js': 8,
  'app/feature-search/preview-svg.js': 2,
  'app/history-inputs.js': 2,
  'app/legend-layout.js': 7,
  'app/legend-layout/canvas-actions.js': 1,
  'app/legend-layout/decoration-continuity.js': 8,
  'app/legend-layout/diagram-drag.js': 6,
  'app/legend-layout/reposition-actions.js': 1,
  'app/legend.js': 6,
  'app/legend/drag-actions.js': 7,
  'app/legend/entry-actions.js': 1,
  'app/legend/sort-actions.js': 1,
  'app/legend/stroke-actions.js': 1,
  'app/legend/track-data-styles.js': 2,
  'app/linear-track-slots.js': 4,
  'app/orthogroups.js': 4,
  'app/preview-runtime.js': 20,
  'app/record-display-options.js': 15,
  'app/record-display/feature-record-rotation.js': 14,
  'app/results.js': 1,
  'app/right-drawer.js': 2,
  'app/rule-matching.js': 15,
  'app/run-analysis.js': 100,
  'app/run-info.js': 5,
  'app/similarity-alignment.js': 35,
  'app/svg-styles.js': 6,
  'app/track-slot-edits.js': 4,
  'app/ui.js': 3,
  'app/watchers.js': 19,
  'services/annotation-state.js': 1,
  'services/artifact-slot.js': 7,
  'services/config.js': 78,
  'services/depth-track-state.js': 4,
  'services/diagram-generation.js': 1,
  'services/feature-catalog.js': 3,
  'services/feature-edit-migration.js': 13,
  'services/feature-metadata-extraction.js': 3,
  'services/feature-placement.js': 12,
  'services/feature-utils.js': 2,
  'services/feature-visibility.js': 4,
  'services/file-imports.js': 1,
  'services/gallery-session-migration.js': 12,
  'services/gallery-session-publication.js': 10,
  'services/history-snapshot.js': 35,
  'services/history.js': 9,
  'services/linear-comparisons.js': 2,
  'services/losat-cache.js': 8,
  'services/losat.js': 2,
  'services/match-sequences.js': 11,
  'services/mode-scoped-migration.js': 26,
  'services/result-normalization.js': 1,
  'services/right-drawer-state.js': 8,
  'services/rule-matchers.js': 10,
  'services/session-active-config-contract.js': 23,
  'services/session-authority.js': 14,
  'services/session-file.js': 1,
  'services/session-request.js': 32,
  'services/session-resources.js': 4,
  'services/specific-color-rules.js': 1,
  'services/svg-result-ingestion.js': 3,
  'state.js': 5,
  'utils/error-normalization.js': 4
};

const REPOSITORY_ROOT = join(dirname(fileURLToPath(import.meta.url)), '..', '..');
const WEB_JS = join(REPOSITORY_ROOT, 'gbdraw', 'web', 'js');
const TSCONFIG = join(REPOSITORY_ROOT, 'tests', 'web', 'types', 'tsconfig.json');

const typescriptApi = async () => {
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
  const [{ API }, { isTypeNode, SyntaxKind }] = await Promise.all([
    import(require.resolve('typescript/unstable/sync')),
    import(require.resolve('typescript/unstable/ast'))
  ]);
  return { API, isTypeNode, SyntaxKind };
};

// The explicit-any nodes of every counted module: Map(path -> [{ line, type }]).
const explicitAnyByModule = async () => {
  const { API, isTypeNode, SyntaxKind } = await typescriptApi();
  const paths = readdirSync(WEB_JS, { recursive: true })
    .map((path) => path.split(sep).join('/'))
    .filter((path) => path.endsWith('.js'))
    .sort();
  const workers = paths.filter((path) => path.startsWith('workers/'));
  const api = new API({ cwd: REPOSITORY_ROOT });
  try {
    const snapshot = api.updateSnapshot({
      openProjects: [TSCONFIG],
      openFiles: workers.map((path) => join(WEB_JS, path))
    });
    const result = new Map();
    paths.forEach((path) => {
      const file = join(WEB_JS, path);
      const project = snapshot.getDefaultProjectForFile(file);
      if (!path.startsWith('workers/')) {
        assert.equal(project?.configFileName, TSCONFIG, `${path} is not in the program of ${TSCONFIG}`);
      }
      const source = project?.program.getSourceFile(file);
      assert.ok(source, `the TypeScript API returned no source file for ${path}`);
      // JSDoc type nodes, each with the index of the nearest enclosing one. A
      // JSDoc block can hang on more than one node, so nodes are keyed by span.
      const nodes = [];
      const enclosing = [];
      const seen = new Set();
      const visit = (node, inJsdoc, up) => {
        let next = up;
        if (inJsdoc && isTypeNode(node)) {
          const key = `${node.kind}:${node.pos}:${node.end}`;
          if (seen.has(key)) return;
          seen.add(key);
          nodes.push(node);
          enclosing.push(up);
          next = nodes.length - 1;
        }
        (node.jsDoc || []).forEach((block) => visit(block, true, next));
        node.forEachChild((child) => visit(child, inJsdoc, next));
      };
      visit(source, false, -1);
      if (!nodes.length) return;
      const types = project.checker.getTypeAtLocation(nodes);
      const isAny = nodes.map((node, index) => node.kind === SyntaxKind.JSDocAllType
        || types[index]?.intrinsicName === 'any');
      const lines = [];
      nodes.forEach((node, index) => {
        if (!isAny[index]) return;
        for (let up = enclosing[index]; up >= 0; up = enclosing[up]) if (isAny[up]) return;
        const text = source.text.slice(node.pos, node.end);
        const start = node.pos + text.length - text.trimStart().length;
        lines.push({ line: source.text.slice(0, start).split('\n').length, type: text.trim().replace(/\s+/g, ' ') });
      });
      if (lines.length) result.set(path, lines);
    });
    return result;
  } finally {
    api.close();
  }
};

test('R14: no module holds more explicit any in JSDoc than its EXPLICIT_ANY_BASELINE entry', async (t) => {
  Object.entries(EXPLICIT_ANY_BASELINE).forEach(([path, count]) => {
    assert.ok(Number.isInteger(count) && count > 0, `EXPLICIT_ANY_BASELINE['${path}'] must be a positive integer`);
  });
  const observed = await explicitAnyByModule();
  const { above, lowerable } = compareCountBaseline(
    Object.fromEntries([...observed].map(([path, lines]) => [path, lines.length])),
    EXPLICIT_ANY_BASELINE
  );
  const failures = above.map(({ key: path, observed: count, recorded }) => {
    // Every line once: `39 ×4` for four on one line, `54 Alias` for a type
    // other than a bare `any`.
    const byLine = new Map();
    observed.get(path).forEach(({ line, type }) => byLine.set(line, [...(byLine.get(line) || []), type]));
    const listed = [...byLine].map(([line, types]) => {
      const named = [...new Set(types.filter((type) => type !== 'any'))];
      return `${line}${types.length > 1 ? ` ×${types.length}` : ''}${named.length ? ` ${named.join(' ')}` : ''}`;
    });
    return `${path}: ${count} explicit any, above its entry ${recorded}. `
      + 'Type the value, use `unknown`, or narrow it; raising an entry is an authority-only change.\n'
      + `      gbdraw/web/js/${path} lines ${listed.join(', ')}`;
  });
  const total = [...observed.values()].reduce((sum, lines) => sum + lines.length, 0);
  t.diagnostic(`explicit any: ${total} in ${observed.size} modules`);
  reportLowerable(t, 'EXPLICIT_ANY_BASELINE', lowerable);
  assert.deepEqual(failures, [], `gbdraw/web/js/ holds more explicit any than EXPLICIT_ANY_BASELINE:\n  ${failures.join('\n  ')}`);
});
