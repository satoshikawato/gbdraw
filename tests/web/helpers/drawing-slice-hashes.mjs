// The sha256 of each slice capture of each drawing (`DRAWING_SLICES` of
// tests/web/helpers/drawing-state.mjs, as tests/web/per-mode-drawings.test.mjs
// hashes them in node), and the comparisons the random walk
// (tests/web/helpers/live-random-walk.mjs) and the Gallery round-trip audit
// (tools/audit/gallery-roundtrip.audit.spec.js) make of two such sets: the one
// reader of a drawing's state by slice for both.

/**
 * The slice hashes of every drawing by drawing ID, and the shown drawing.
 * `drawings` hashes each capture as it is. `renamed` hashes it as a check
 * across a Session Load must compare it: object keys in order (a Load builds
 * a map such as `featureOverrides` in another insertion order; the request
 * reads them in canonical order, `canonicalFeatureOverrides`); the file
 * source of each ref (`sourceId`, services/history-files.js, read only per
 * drawing by the Linear File cards' grouping, `groupLinearSourceRecords`)
 * renamed in order of appearance within the drawing, so it says which refs
 * share a source and not which file of the project the store named it after;
 * then the runtime identities renamed in order of appearance within the
 * capture: the file IDs (`file-<n>`) and the Linear record and Circular File
 * card UIDs (`linear-seq-…`, `circular-source-…`, services/drawing/inputs.js),
 * which a page makes anew when it reads a Session; and without the fields of
 * `omit`. Two states that differ only by those have equal `renamed` hashes; a
 * check within one page compares `drawings`.
 * @typedef {{ active: string, drawings: Record<string, Record<string, string>>, renamed: Record<string, Record<string, string>> }} SliceHashes
 */
/**
 * A drawing whose slices differ between two hash sets, in a phase of a check.
 * @typedef {{ phase: string, drawing: string, slices: string[] }} SliceChange
 */
/**
 * Fields a `renamed` hash leaves out, by slice name: each a dotted path into
 * the slice's capture, `*` standing for every key or index at its place.
 * @typedef {Readonly<Record<string, readonly string[]>>} SliceFieldsOff
 */

/**
 * @param {import('@playwright/test').Page} page
 * @param {{ omit?: SliceFieldsOff }} [options]
 * @returns {Promise<SliceHashes>}
 */
export const sliceHashes = (page, { omit = {} } = {}) => page.evaluate(async (omit) => {
  const { state } = await import('/gbdraw/web/js/state.js');
  const { DRAWING_SLICES } = await import('/tests/web/helpers/drawing-state.mjs');
  const hex = async (text) => [...new Uint8Array(await crypto.subtle.digest('SHA-256', new TextEncoder().encode(text)))]
    .map((byte) => byte.toString(16).padStart(2, '0')).join('');
  const RUNTIME_IDENTITY = /(?:linear-seq|circular-source)-[0-9a-f]+(?:-[0-9a-f]+)*|file-\d+/g;
  const isObject = (value) => value !== null && typeof value === 'object';
  // Keys in order, compared with their runtime identities masked, so the
  // order does not depend on the names a page gave.
  const masked = (key) => key.replace(RUNTIME_IDENTITY, 'runtime-identity');
  const byMaskedKey = (a, b) => (masked(a) < masked(b) ? -1 : masked(a) > masked(b) ? 1 : 0);
  const canonical = (value) => (Array.isArray(value) ? value.map(canonical)
    : isObject(value) ? Object.fromEntries(Object.keys(value).sort(byMaskedKey).map((key) => [key, canonical(value[key])])) : value);
  const withoutField = (value, [head, ...rest]) => {
    if (!isObject(value)) return value;
    const copy = Array.isArray(value) ? [...value] : { ...value };
    for (const key of head === '*' ? Object.keys(copy) : [head]) {
      if (!Object.hasOwn(copy, key)) continue;
      if (rest.length) copy[key] = withoutField(copy[key], rest);
      else if (Array.isArray(copy)) copy[key] = '(left out)';
      else delete copy[key];
    }
    return copy;
  };
  const renamedText = (text) => {
    const names = new Map();
    return text.replace(RUNTIME_IDENTITY, (identity) => {
      if (!names.has(identity)) names.set(identity, `runtime-identity-${names.size + 1}`);
      return names.get(identity);
    });
  };
  const drawings = {};
  const renamed = {};
  for (const drawing of state.drawingList()) {
    drawings[drawing.id] = {};
    renamed[drawing.id] = {};
    const sources = new Map();
    const sourceName = (id) => {
      if (!sources.has(id)) sources.set(id, `source-${sources.size + 1}`);
      return sources.get(id);
    };
    const sourcesRenamed = (value) => (Array.isArray(value) ? value.map(sourcesRenamed)
      : isObject(value) ? Object.fromEntries(Object.entries(value)
        .map(([key, field]) => [key, key === 'sourceId' && typeof field === 'string' ? sourceName(field) : sourcesRenamed(field)])) : value);
    for (const [name, slice] of Object.entries(DRAWING_SLICES)) {
      const captured = slice.capture(drawing.slices[name]);
      drawings[drawing.id][name] = await hex(JSON.stringify(captured));
      const kept = (omit[name] || []).reduce((value, path) => withoutField(value, path.split('.')), captured);
      renamed[drawing.id][name] = await hex(renamedText(JSON.stringify(sourcesRenamed(canonical(kept)))));
    }
  }
  return { active: state.activeDrawing().id, drawings, renamed };
}, omit);

/**
 * The slices of `drawingId` whose captures differ between two hash sets, by
 * their `drawings` hashes or, with `by: 'renamed'`, their `renamed` hashes.
 * @param {SliceHashes} before @param {SliceHashes} after @param {string} drawingId
 * @param {{ by?: 'drawings' | 'renamed' }} [options]
 */
export const changedSlices = (before, after, drawingId, { by = 'drawings' } = {}) => Object.keys(before[by][drawingId])
  .filter((name) => before[by][drawingId][name] !== after[by][drawingId][name]);

/**
 * What a phase that must give back every drawing did not give back: the
 * changed slices of each drawing, a drawing it lost or added, and the shown
 * drawing when it differs. `by: 'renamed'` compares the captures with their
 * runtime identities renamed, as a check across a Session Load must.
 * @param {string} phase
 * @param {SliceHashes} before
 * @param {SliceHashes} after
 * @param {{ by?: 'drawings' | 'renamed' }} [options]
 * @returns {SliceChange[]}
 */
export const restoredChanges = (phase, before, after, { by = 'drawings' } = {}) => [
  ...[...new Set([...Object.keys(before.drawings), ...Object.keys(after.drawings)])].map((drawing) => ({
    phase,
    drawing,
    slices: !after.drawings[drawing] ? ['(the drawing is gone)'] : !before.drawings[drawing] ? ['(a new drawing)'] : changedSlices(before, after, drawing, { by })
  })),
  { phase, drawing: '(shown drawing)', slices: before.active === after.active ? [] : [`${before.active} shown before, ${after.active} after`] }
].filter(({ slices }) => slices.length);
