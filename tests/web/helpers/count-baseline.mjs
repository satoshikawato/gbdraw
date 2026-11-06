import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { createRequire } from 'node:module';

// A registered count baseline holds "observed <= recorded" per key
// (docs/internal/WEB_CHANGE_POLICY.md "Design-rule co-change"): a key without
// an entry counts as 0, lowering or removing an entry is optional, and adding
// or raising one is an authority-only change. `above` lists the keys that fail
// and `lowerable` the entries a pull request may lower.
export const compareCountBaseline = (observed, baseline) => {
  const rows = [...new Set([...Object.keys(observed), ...Object.keys(baseline)])].sort()
    .map((key) => ({ key, observed: observed[key] ?? 0, recorded: baseline[key] ?? 0 }));
  return {
    above: rows.filter((row) => row.observed > row.recorded),
    lowerable: rows.filter((row) => row.observed < row.recorded)
      .map((row) => `${row.key} ${row.recorded} -> ${row.observed}`)
  };
};

// A copy of the module at `source` (a file URL) with `probe` appended, loaded
// as a data URL, so `probe` can export the size of a private table. A data URL
// resolves no relative specifier, so the copy's relative imports name the
// source files.
export const probeModule = (source, probe) => import(`data:text/javascript,${encodeURIComponent(`${
  readFileSync(source, 'utf8').replace(/(\bfrom\s*)(['"])(\.\.?\/[^'"]+)\2/g,
    (_, from, quote, specifier) => `${from}${quote}${new URL(specifier, source).href}${quote}`)
}\n${probe}`)}`);

// One `t.diagnostic` line that names the entries that may be lowered, if any.
export const reportLowerable = (t, name, lowerable) => {
  if (lowerable.length) t.diagnostic(`${name} may be lowered: ${lowerable.join(', ')}`);
};

const require = createRequire(import.meta.url);

// The path of the package.json of the TypeScript that package.json pins as a
// devDependency, after checking that the installed one is that version. The R14
// guards and the Web function size guard read the Web modules with it.
export const pinnedTypeScript = () => {
  let packagePath;
  try {
    packagePath = require.resolve('typescript/package.json');
  } catch {
    assert.fail('typescript is not installed: run `npm ci` (the Web guards run the pinned devDependency)');
  }
  const pinned = JSON.parse(readFileSync(new URL('../../../package.json', import.meta.url), 'utf8')).devDependencies.typescript;
  const installed = JSON.parse(readFileSync(packagePath, 'utf8')).version;
  assert.equal(installed, pinned, `typescript ${installed} is installed but package.json pins ${pinned}: run \`npm ci\``);
  return packagePath;
};

// The pinned TypeScript's API: its synchronous client `API` and its syntax-tree
// module `ast` (`SyntaxKind`, `isTypeNode`, ...).
export const typescriptApi = async () => {
  pinnedTypeScript();
  const [{ API }, ast] = await Promise.all([
    import(require.resolve('typescript/unstable/sync')),
    import(require.resolve('typescript/unstable/ast'))
  ]);
  return { API, ast };
};
