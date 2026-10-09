// R14 "Typed boundaries" (gbdraw/web/CLAUDE.md): the explicit-any ratchet.
// Explicit `any` in JSDoc is type debt. EXPLICIT_ANY_BASELINE records it per
// module under gbdraw/web/js/ (paths relative to it); a module without an
// entry, including a new module, has none. The guard asserts that each module
// holds exactly its entry: a module above its entry fails (type the value, use
// `unknown`, or narrow), and a module below it fails until the entry is
// lowered, so the map only shrinks.
//
// What counts: each JSDoc type node whose type, as the TypeScript checker of
// the R14 guard config (tests/web/types/tsconfig.json) resolves it, is
// intrinsic `any`: `any`, `*`, `Object` (which JSDoc reads as `any` while
// `noImplicitAny` is off), and each use of a type that resolves to `any`, such
// as a typedef alias of `any`, so an alias does not hide it. A node inside a
// counted node is not counted again; `Record<string, any>` counts once, for
// its `any` argument. `unknown`, `object`, and `Function` do not count.
// Generated modules (`*.generated.js`) are not counted. Workers are outside
// the guard config; the checker resolves them in its inferred project, whose
// `noImplicitAny` is off as well.
//
// EXPLICIT_ANY_BASELINE is a registered design-rule allowlist
// (tools/web-design-rule-guards.json, R14, kind `count-map`): a pull request
// lowers or removes the entries of the modules it types, and adding an entry
// or raising one is an authority-only change. Owner decision R15-5
// (2026-10-09); values at dev ffa5a3ba.
import assert from 'node:assert/strict';
import { readFileSync, readdirSync } from 'node:fs';
import { createRequire } from 'node:module';
import { dirname, join, sep } from 'node:path';
import test from 'node:test';
import { fileURLToPath } from 'node:url';

const EXPLICIT_ANY_BASELINE = {
  'app/annotations.js': 8,
  'app/annotations/record-catalog.js': 3,
  'app/annotations/target-actions.js': 7,
  'app/annotations/validation.js': 2,
  'app/app-setup.js': 21,
  'app/auto-value-display.js': 2,
  'app/candidate-render.js': 12,
  'app/circular-track-slots.js': 12,
  'app/comparison-ui.js': 2,
  'app/feature-editor.js': 24,
  'app/feature-editor/color-actions.js': 48,
  'app/feature-editor/feature-edit-table.js': 8,
  'app/feature-editor/label-actions.js': 16,
  'app/feature-editor/pattern-drafts.js': 3,
  'app/feature-editor/placement-actions.js': 10,
  'app/feature-editor/rule-actions.js': 31,
  'app/feature-editor/svg-actions.js': 22,
  'app/feature-editor/visibility-actions.js': 8,
  'app/feature-search/preview-actions.js': 8,
  'app/feature-search/preview-svg.js': 2,
  'app/feature-search/search-core.js': 10,
  'app/feature-selection.js': 1,
  'app/history-inputs.js': 2,
  'app/legend-layout.js': 8,
  'app/legend-layout/canvas-actions.js': 1,
  'app/legend-layout/composition-actions.js': 13,
  'app/legend-layout/decoration-continuity.js': 8,
  'app/legend-layout/diagram-drag.js': 6,
  'app/legend-layout/reposition-actions.js': 1,
  'app/legend.js': 7,
  'app/legend/drag-actions.js': 7,
  'app/legend/entry-actions.js': 13,
  'app/legend/sort-actions.js': 1,
  'app/legend/stroke-actions.js': 10,
  'app/legend/track-data-styles.js': 7,
  'app/linear-record-selector.js': 12,
  'app/linear-track-slots.js': 4,
  'app/linear-typography.js': 2,
  'app/losat-settings.js': 2,
  'app/orthogroups.js': 10,
  'app/pairwise-match-popup.js': 21,
  'app/palettes.js': 2,
  'app/preview-runtime.js': 21,
  'app/record-display-options.js': 16,
  'app/record-display/feature-record-rotation.js': 24,
  'app/results.js': 1,
  'app/right-drawer.js': 2,
  'app/rule-matching.js': 24,
  'app/run-analysis.js': 139,
  'app/run-info.js': 13,
  'app/similarity-alignment.js': 40,
  'app/svg-styles.js': 8,
  'app/track-slot-colors.js': 4,
  'app/track-slot-edits.js': 4,
  'app/ui.js': 3,
  'app/watchers.js': 22,
  'services/annotation-state.js': 1,
  'services/artifact-slot.js': 10,
  'services/bounded-json-transport.js': 2,
  'services/circular-track-measure.js': 1,
  'services/circular-track-slot-model.js': 6,
  'services/config.js': 119,
  'services/depth-track-state.js': 11,
  'services/diagram-generation.js': 11,
  'services/diagram-resource-staging.js': 2,
  'services/feature-catalog.js': 7,
  'services/feature-edit-migration.js': 13,
  'services/feature-metadata-extraction.js': 7,
  'services/feature-placement.js': 12,
  'services/feature-utils.js': 2,
  'services/feature-visibility.js': 15,
  'services/file-imports.js': 1,
  'services/gallery-session-migration.js': 13,
  'services/gallery-session-publication.js': 18,
  'services/history-snapshot.js': 67,
  'services/history.js': 13,
  'services/imported-comparison-intent.js': 4,
  'services/layout-preferences.js': 3,
  'services/legacy-similarity-alignment.js': 3,
  'services/linear-comparisons.js': 14,
  'services/linear-sources.js': 4,
  'services/linear-track-slot-model.js': 1,
  'services/losat-cache.js': 8,
  'services/losat-runtime.js': 1,
  'services/losat.js': 2,
  'services/main-session-comparison-frame.js': 2,
  'services/match-sequences.js': 12,
  'services/mode-scoped-migration.js': 26,
  'services/orthogroup-feature-metadata.js': 7,
  'services/record-options.js': 1,
  'services/reset.js': 5,
  'services/result-normalization.js': 1,
  'services/right-drawer-state.js': 8,
  'services/rule-matchers.js': 11,
  'services/session-active-config-contract.js': 24,
  'services/session-authority.js': 24,
  'services/session-feature-recovery.js': 5,
  'services/session-file.js': 1,
  'services/session-request.js': 72,
  'services/session-resources.js': 10,
  'services/specific-color-rules.js': 2,
  'services/standalone-interactivity.js': 2,
  'services/svg-result-ingestion.js': 10,
  'services/svg-serialization.js': 2,
  'services/track-slot-display.js': 1,
  'services/track-slot-validation.js': 1,
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
    .filter((path) => path.endsWith('.js') && !path.endsWith('.generated.js'))
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

test('R14: each module holds exactly its EXPLICIT_ANY_BASELINE count of explicit any in JSDoc', async (t) => {
  Object.entries(EXPLICIT_ANY_BASELINE).forEach(([path, count]) => {
    assert.ok(Number.isInteger(count) && count > 0, `EXPLICIT_ANY_BASELINE['${path}'] must be a positive integer`);
  });
  const observed = await explicitAnyByModule();
  const failures = [];
  [...new Set([...observed.keys(), ...Object.keys(EXPLICIT_ANY_BASELINE)])].sort().forEach((path) => {
    const lines = observed.get(path) || [];
    const recorded = EXPLICIT_ANY_BASELINE[path] || 0;
    if (lines.length > recorded) {
      // Every line once: `39 ×4` for four on one line, `54 Alias` for a type
      // other than a bare `any`.
      const byLine = new Map();
      lines.forEach(({ line, type }) => byLine.set(line, [...(byLine.get(line) || []), type]));
      const listed = [...byLine].map(([line, types]) => {
        const named = [...new Set(types.filter((type) => type !== 'any'))];
        return `${line}${types.length > 1 ? ` ×${types.length}` : ''}${named.length ? ` ${named.join(' ')}` : ''}`;
      });
      failures.push(`${path}: ${lines.length} explicit any, above its entry ${recorded}. `
        + 'Type the value, use `unknown`, or narrow it; raising an entry is an authority-only change.\n'
        + `      gbdraw/web/js/${path} lines ${listed.join(', ')}`);
    } else if (lines.length < recorded) {
      failures.push(lines.length
        ? `${path}: ${lines.length} explicit any, below its entry ${recorded}: lower its entry to ${lines.length}`
        : `${path}: no explicit any left: remove its entry (${recorded})`);
    }
  });
  const total = [...observed.values()].reduce((sum, lines) => sum + lines.length, 0);
  t.diagnostic(`explicit any: ${total} in ${observed.size} modules`);
  assert.deepEqual(failures, [], `EXPLICIT_ANY_BASELINE does not match gbdraw/web/js/:\n  ${failures.join('\n  ')}`);
});
