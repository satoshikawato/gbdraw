// R14 "Typed boundaries" (gbdraw/web/CLAUDE.md). Modules under
// gbdraw/web/js/ except workers/ are checked by `tsc --noEmit` with JSDoc
// types (tests/web/types/tsconfig.json); a checked module starts with the
// line `// @ts-check`. The guard asserts:
//   1. every module is checked or listed in UNCHECKED_MODULES, never both;
//   2. one `tsc` run reports no diagnostic and its program holds every module;
//   3. a checked module types every parameter of its exported create*/setup*
//      factories, and no JSDoc outside a composition root names another
//      owner's factory result (`ReturnType<typeof createX>`);
//   4. no `@ts-ignore`, `@ts-expect-error`, or `@ts-nocheck`;
//   5. type imports resolve to a module and follow the R13 layers.
//
// UNCHECKED_MODULES is a registered design-rule allowlist
// (tools/web-design-rule-guards.json, R14, kind `set`): a pull request that
// adds `// @ts-check` to a module removes it here in the same change, and an
// addition is an authority-only change. A new module is checked from its
// first commit. Paths are relative to gbdraw/web/js/.
import assert from 'node:assert/strict';
import { spawnSync } from 'node:child_process';
import { readFileSync, readdirSync } from 'node:fs';
import { createRequire } from 'node:module';
import { dirname, isAbsolute, join, posix, relative, sep } from 'node:path';
import test from 'node:test';
import { fileURLToPath } from 'node:url';

import { maskJavaScript } from '../../tools/web-change-source.mjs';
import { WEB_OWNER_GRAPH_DEFAULTS } from '../../tools/web-owner-graph-detectors.mjs';

const UNCHECKED_MODULES = new Set([
  'app.js',
  'app/annotations.js',
  'app/annotations/record-catalog.js',
  'app/annotations/record-selector.js',
  'app/annotations/state.js',
  'app/annotations/table-codec.js',
  'app/annotations/target-actions.js',
  'app/annotations/validation.js',
  'app/auto-value-display.js',
  'app/candidate-render.js',
  'app/color-utils.js',
  'app/comparison-ui.js',
  'app/conservation-series.js',
  'app/current-option-values.js',
  'app/definition-line-style-state.js',
  'app/depth-track-state.js',
  'app/depth-tracks.js',
  'app/feature-dom.js',
  'app/feature-metadata-extraction.js',
  'app/feature-selection.js',
  'app/feature-selector.js',
  'app/feature-sequence-fasta.js',
  'app/feature-utils.js',
  'app/feature-visibility.js',
  'app/file-imports.js',
  'app/genbank-header.js',
  'app/history-inputs.js',
  'app/history-shortcuts.js',
  'app/layout-preferences.js',
  'app/linear-comparisons.js',
  'app/linear-label-visibility.js',
  'app/linear-record-layout.js',
  'app/linear-record-selector.js',
  'app/linear-sources.js',
  'app/linear-typography.js',
  'app/losat-cache.js',
  'app/losat-normalization.js',
  'app/losat-settings.js',
  'app/match-sequences.js',
  'app/orthogroups.js',
  'app/pairwise-match-popup.js',
  'app/palettes.js',
  'app/plot-title-position.js',
  'app/preview-runtime.js',
  'app/python-helpers.js',
  'app/record-discovery.js',
  'app/record-groups.js',
  'app/record-options.js',
  'app/record-source-coordinates.js',
  'app/results.js',
  'app/right-drawer.js',
  'app/run-info.js',
  'app/session-feature-metadata.js',
  'app/specific-color-rules.js',
  'app/track-slot-colors.js',
  'app/track-slot-display.js',
  'app/track-slot-validation.js',
  'app/ui.js',
  'components.js',
  'config.js',
  'mode-profiles.generated.js',
  'mode-profiles.js',
  'services/error-normalization.js',
  'services/export.js',
  'services/history-snapshot.js',
  'services/history.js',
  'services/losat-thread-plan.js',
  'services/losat.js',
  'services/pdf-fonts.js',
  'services/pyodide-assets.js',
  'services/standalone-interactivity-assets.js',
  'services/svg-result-ingestion.js',
  'services/svg-result-normalization.js',
  'services/svg-serialization.js',
  'services/text-download.js',
  'state.js',
  'utils/clipboard.js',
  'utils/feature-rendering.js',
  'utils/optional-positive-number.js',
  'utils/png.js',
  'utils/tsv-cell.js',
  'utils/zip.js',
  'web-ux-profile.js'
]);

const REPOSITORY_ROOT = join(dirname(fileURLToPath(import.meta.url)), '..', '..');
const WEB_JS = join(REPOSITORY_ROOT, 'gbdraw', 'web', 'js');
const TSCONFIG = 'tests/web/types/tsconfig.json';
const COMPOSITION_ROOTS = new Set(WEB_OWNER_GRAPH_DEFAULTS.compositionRoots);
const PRAGMA_LINE = /^\/\/ @ts-check\r?\n/;
const FACTORY_NAME = /^(?:create|setup)(?![a-z])/;
const UNDECLARED_TYPES = new Set(['any', '*', 'object', 'Object']);

const toPosix = (path) => path.split(sep).join('/');
const isWorker = (path) => path.startsWith('workers/');
const SOURCES = new Map(
  readdirSync(WEB_JS, { recursive: true })
    .map(toPosix)
    .filter((path) => path.endsWith('.js'))
    .sort()
    .map((path) => [path, readFileSync(join(WEB_JS, path), 'utf8')])
);
const MODULES = [...SOURCES.keys()].filter((path) => !isWorker(path));
const CHECKED = MODULES.filter((path) => PRAGMA_LINE.test(SOURCES.get(path)));

// The source with everything outside comments blanked (offsets and line
// breaks kept), so tags are read from comments only.
const commentText = (source) => {
  const code = maskJavaScript(source, { strings: false });
  const chars = source.split('');
  for (let index = 0; index < chars.length; index += 1) {
    if (code[index] === chars[index] && chars[index] !== '\n') chars[index] = ' ';
  }
  return chars.join('');
};
const COMMENTS = new Map([...SOURCES].map(([path, source]) => [path, commentText(source)]));
const lineOf = (text, index) => text.slice(0, index).split('\n').length;

// Index of the bracket closing the one at `open` in masked code, or -1.
const closingBracket = (code, open) => {
  const closers = { '(': ')', '[': ']', '{': '}' };
  const stack = [];
  for (let index = open; index < code.length; index += 1) {
    const character = code[index];
    if (closers[character]) stack.push(closers[character]);
    else if (character === ')' || character === ']' || character === '}') {
      if (stack.pop() !== character) return -1;
      if (!stack.length) return index;
    }
  }
  return -1;
};
const countParameters = (list) => {
  let depth = 0;
  let count = 0;
  let pending = false;
  for (const character of list) {
    if ('([{'.includes(character)) depth += 1;
    else if (')]}'.includes(character)) depth -= 1;
    if (character === ',' && depth === 0) {
      if (pending) count += 1;
      pending = false;
    } else if (!/\s/.test(character)) pending = true;
  }
  return count + (pending ? 1 : 0);
};

// The parameter count of the function literal that starts at `index` in
// masked code, or null when the initializer is not a function literal.
const functionLiteralParameters = (code, index) => {
  const rest = code.slice(index).replace(/^async\s+/, '');
  if (/^(?!function\b)[A-Za-z_$][\w$]*\s*=>/.test(rest)) return 1;
  const functionKeyword = rest.match(/^function\s*\*?\s*(?:[A-Za-z_$][\w$]*)?\s*(?=\()/);
  if (!functionKeyword && rest[0] !== '(') return null;
  const open = code.length - rest.length + (functionKeyword ? functionKeyword[0].length : 0);
  const close = closingBracket(code, open);
  if (close < 0 || (!functionKeyword && !/^\s*=>/.test(code.slice(close + 1)))) return null;
  return countParameters(code.slice(open + 1, close));
};
// The parameter count of a function declaration whose list opens at `open`.
const declarationParameters = (code, open) => {
  const close = closingBracket(code, open);
  return close < 0 ? null : countParameters(code.slice(open + 1, close));
};

// Exported create*/setup* functions: `export function`, `export const x =
// <function>`, and local functions named in an `export { … }` list.
const exportedFactories = (source) => {
  const code = maskJavaScript(source);
  const factories = [];
  const add = (name, start, parameters) => {
    if (FACTORY_NAME.test(name) && parameters !== null) factories.push({ name, start, parameters });
  };
  for (const match of code.matchAll(/(?<=^|\n)[ \t]*export\s+(?:async\s+)?function\b\s*\*?\s*([A-Za-z_$][\w$]*)\s*(?=\()/g)) {
    add(match[1], match.index + match[0].search(/\S/), declarationParameters(code, match.index + match[0].length));
  }
  for (const match of code.matchAll(/(?<=^|\n)[ \t]*export\s+(?:const|let|var)\s+([A-Za-z_$][\w$]*)\s*=\s*/g)) {
    add(match[1], match.index + match[0].search(/\S/), functionLiteralParameters(code, match.index + match[0].length));
  }
  for (const match of code.matchAll(/(?<=^|\n)[ \t]*export\s*\{([^}]*)\}(?!\s*from\b)/g)) {
    match[1].split(',').map((item) => item.trim().split(/\s+as\s+/)).forEach(([local, exported = local]) => {
      if (!local || !FACTORY_NAME.test(exported)) return;
      const escaped = local.replace(/[.*+?^${}()|[\]\\]/g, '\\$&');
      const declaration = code.match(new RegExp(`(?<=^|\\n)[ \\t]*(?:(?:async\\s+)?function\\s*\\*?\\s*${escaped}\\s*(?=\\()|(?:const|let|var)\\s+${escaped}\\s*=\\s*)`));
      if (!declaration) return;
      const start = declaration.index + declaration[0].search(/\S/);
      const end = declaration.index + declaration[0].length;
      const parameters = /^\s*(?:const|let|var)\b/.test(declaration[0])
        ? functionLiteralParameters(code, end)
        : declarationParameters(code, end);
      add(exported, start, parameters);
    });
  }
  return factories;
};

// The `/** … */` block directly above `start`, or null.
const jsdocAbove = (source, start) => {
  const before = source.slice(0, start).trimEnd();
  if (!before.endsWith('*/')) return null;
  const open = before.lastIndexOf('/**');
  return open < 0 ? null : before.slice(open);
};
// The top-level `@param` tags of a JSDoc block, before any `@typedef`,
// `@callback`, or `@overload` tag: `{ type, name }`, `type` null when absent.
const jsdocParameters = (block) => {
  const text = block.slice(3, -2).replace(/\n[ \t]*\*(?!\/)/g, '\n');
  const cut = text.search(/@(?:typedef|callback|overload)\b/);
  const scoped = cut < 0 ? text : text.slice(0, cut);
  const parameters = [];
  for (const match of scoped.matchAll(/@param\b/g)) {
    let index = match.index + match[0].length;
    while (/\s/.test(scoped[index] || '')) index += 1;
    let type = null;
    if (scoped[index] === '{') {
      const close = closingBracket(scoped, index);
      if (close > index) {
        type = scoped.slice(index + 1, close).replace(/\s+/g, ' ').trim();
        index = close + 1;
      }
    }
    const name = scoped.slice(index).match(/^\s*\[?\s*([A-Za-z_$][\w$]*(?:\.[\w$]+|\[\])*)/)?.[1] || '';
    if (!/[.[]/.test(name)) parameters.push({ type, name });
  }
  return parameters;
};
const isUndeclaredType = (type) => !type
  || UNDECLARED_TYPES.has(type.replace(/^\.\.\./, '').replace(/^[?!]/, '').replace(/=$/, '').trim());

// Relative type-import specifiers in comments: `@import … from '…'` and
// `import('…')`.
const typeImports = (comments) => [
  ...comments.matchAll(/@import\s*(?:\{[^}]*\}|\*\s*as\s+[A-Za-z_$][\w$]*|[A-Za-z_$][\w$]*)\s*from\s*(['"])([^'"\n]+)\1/g),
  ...comments.matchAll(/\bimport\s*\(\s*(['"])([^'"\n]+)\1\s*\)/g)
].map((match) => ({ specifier: match[2], index: match.index }));
const isLowerLayer = (path) => path === 'state.js' || path.startsWith('services/') || path.startsWith('utils/');

test('R14: every module is checked or listed in UNCHECKED_MODULES, never both', () => {
  const failures = [];
  MODULES.forEach((path) => {
    const checked = PRAGMA_LINE.test(SOURCES.get(path));
    const misplaced = [...COMMENTS.get(path).matchAll(/@ts-check\b/g)].find(({ index }) => !checked || index !== 3);
    if (misplaced) {
      failures.push(`${path}:${lineOf(COMMENTS.get(path), misplaced.index)}: \`@ts-check\` goes on line 1 as \`// @ts-check\` and nowhere else`);
    }
    if (checked && UNCHECKED_MODULES.has(path)) failures.push(`${path} has \`// @ts-check\`: remove it from UNCHECKED_MODULES`);
    if (!checked && !UNCHECKED_MODULES.has(path)) {
      failures.push(`${path} is not checked: add \`// @ts-check\` as its first line (a new module is checked from its first commit)`);
    }
  });
  [...UNCHECKED_MODULES].filter((path) => !MODULES.includes(path)).forEach((path) => {
    failures.push(`UNCHECKED_MODULES names ${path}, which is not a module under gbdraw/web/js/ outside workers/: remove it from UNCHECKED_MODULES`);
  });
  assert.deepEqual(failures, []);
});

test('R14: tsc reports no diagnostic and its program holds every module', () => {
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
    [join(dirname(packagePath), 'bin', 'tsc'), '-p', TSCONFIG, '--listFiles', '--pretty', 'false'],
    { cwd: REPOSITORY_ROOT, encoding: 'utf8', maxBuffer: 1 << 26 }
  );
  assert.ifError(run.error);
  const lines = `${run.stdout}\n${run.stderr}`.split(/\r?\n/).filter((line) => line.trim());
  const diagnostics = lines.filter((line) => !isAbsolute(line) || /\): error TS\d+/.test(line));
  assert.deepEqual(diagnostics, [], `tsc -p ${TSCONFIG} reported diagnostics`);
  assert.equal(run.status, 0, `tsc exited with ${run.status}`);
  const program = new Set(lines.map((line) => toPosix(relative(WEB_JS, line))));
  const missing = MODULES.filter((path) => !program.has(path));
  assert.deepEqual(missing, [], `modules missing from the tsc program: check "include" in ${TSCONFIG}`);
});

test('R14: checked factories declare their parameters; whole-object port types stay in composition roots', () => {
  const failures = [];
  CHECKED.forEach((path) => {
    const source = SOURCES.get(path);
    exportedFactories(source).forEach(({ name, start, parameters }) => {
      if (parameters === 0) return;
      const where = `${path}:${lineOf(source, start)} ${name}`;
      const block = jsdocAbove(source, start);
      if (!block) {
        failures.push(`${where}: add a JSDoc block with one \`@param {T}\` per parameter (${parameters})`);
        return;
      }
      const tags = jsdocParameters(block);
      if (tags.length !== parameters) {
        failures.push(`${where}: ${parameters} parameter(s) but ${tags.length} top-level \`@param\` tag(s)`);
      }
      tags.filter(({ type }) => isUndeclaredType(type)).forEach(({ type, name: tag }) => {
        failures.push(`${where}: \`@param\` ${tag || '(unnamed)'} has type ${type === null ? '(none)' : `{${type}}`}; declare the port or option type`);
      });
    });
  });
  MODULES.filter((path) => !COMPOSITION_ROOTS.has(path)).forEach((path) => {
    const comments = COMMENTS.get(path);
    for (const match of comments.matchAll(/ReturnType\s*<\s*typeof\s+(?:import\s*\([^)]*\)\s*\.\s*)?create/g)) {
      failures.push(`${path}:${lineOf(comments, match.index)}: \`ReturnType<typeof create…>\` is a whole-object port (R13); declare the port function types`);
    }
  });
  assert.deepEqual(failures, []);
});

test('R14: no @ts-ignore, @ts-expect-error, or @ts-nocheck', () => {
  const failures = [];
  COMMENTS.forEach((comments, path) => {
    for (const match of comments.matchAll(/@ts-(?:ignore|expect-error|nocheck)\b/g)) {
      failures.push(`${path}:${lineOf(comments, match.index)}: ${match[0]}; fix the type, cast with JSDoc, or fix the code`);
    }
  });
  assert.deepEqual(failures, []);
});

test('R14: type imports resolve to a module and follow the R13 layers', () => {
  const failures = [];
  MODULES.forEach((path) => {
    const comments = COMMENTS.get(path);
    typeImports(comments).forEach(({ specifier, index }) => {
      const where = `${path}:${lineOf(comments, index)} '${specifier}'`;
      const target = /^\.\.?\//.test(specifier) ? posix.normalize(posix.join(posix.dirname(path), specifier)) : null;
      if (!target || !MODULES.includes(target)) {
        failures.push(`${where}: a type import names a module under gbdraw/web/js/ outside workers/ by relative path`);
        return;
      }
      if (isLowerLayer(path) && target.startsWith('app/')) failures.push(`${where}: state.js, services/, and utils/ do not import types from app/`);
      if (COMPOSITION_ROOTS.has(target) && !COMPOSITION_ROOTS.has(path)) {
        failures.push(`${where}: only a composition root imports types from a composition root`);
      }
    });
  });
  assert.deepEqual(failures, []);
});
