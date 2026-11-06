import { readFileSync } from 'node:fs';

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
