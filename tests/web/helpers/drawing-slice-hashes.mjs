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
 * reads them in canonical order, `canonicalFeatureOverrides`), sorted by key
 * and then by value, both with their runtime identities masked, so the order
 * does not depend on the names a page gave (two entries tied on both read as
 * the same text in either order once renamed); the file
 * source of each ref (`sourceId`, services/history-files.js, read only per
 * drawing by the Linear File cards' grouping, `groupLinearSourceRecords`)
 * renamed in order of appearance within the drawing, so it says which refs
 * share a source and not which file of the project the store named it after;
 * then the runtime identities renamed in order of appearance within the
 * capture: the file IDs (`file-<n>`) and the Linear record and Circular File
 * card UIDs (`linear-seq-…`, `circular-source-…`, services/drawing/inputs.js),
 * which a page makes anew when it reads a Session; and without the fields
 * left out (`leftOut`). Two states that differ only by those have equal
 * `renamed` hashes; a check within one page compares `drawings`.
 * `leftOut` names, per drawing and slice, the fields its `renamed` hash left
 * out, each a list of keys (runtime identities masked); a second call given
 * it leaves the same fields out, whatever they hold there.
 * @typedef {{ active: string, drawings: Record<string, Record<string, string>>, renamed: Record<string, Record<string, string>>,
 *   leftOut: FieldsLeftOut }} SliceHashes
 */
/** @typedef {Record<string, Record<string, string[][]>>} FieldsLeftOut */
/**
 * A drawing whose slices differ between two hash sets, in a phase of a check.
 * @typedef {{ phase: string, drawing: string, slices: string[] }} SliceChange
 */
/**
 * Fields a `renamed` hash leaves out, by slice name: each a dotted `path` into
 * the slice's capture, `*` standing for every key or index at its place, left
 * out where it holds one of `when` (a missing field holds null), or always
 * when `when` is absent.
 * @typedef {{ path: string, when?: readonly (string | number | boolean | null)[] }} FieldOff
 * @typedef {Readonly<Record<string, readonly FieldOff[]>>} SliceFieldsOff
 */

/**
 * The slice hashes of every drawing. `omit` names the fields to leave out
 * by what they hold here; `leftOut` (a `leftOut` an earlier call returned)
 * names the fields to leave out as that call did, so a check compares
 * two states without the fields where the first held a defect's value.
 * @param {import('@playwright/test').Page} page
 * @param {{ omit?: SliceFieldsOff, leftOut?: FieldsLeftOut }} [options]
 * @returns {Promise<SliceHashes>}
 */
export const sliceHashes = (page, { omit = {}, leftOut = null } = {}) => page.evaluate(async ({ omit, leftOut }) => {
  const { state } = await import('/gbdraw/web/js/state.js');
  const { DRAWING_SLICES } = await import('/tests/web/helpers/drawing-state.mjs');
  const hex = async (text) => [...new Uint8Array(await crypto.subtle.digest('SHA-256', new TextEncoder().encode(text)))]
    .map((byte) => byte.toString(16).padStart(2, '0')).join('');
  const RUNTIME_IDENTITY = /(?:linear-seq|circular-source)-[0-9a-f]+(?:-[0-9a-f]+)*|file-\d+/g;
  const isObject = (value) => value !== null && typeof value === 'object';
  const masked = (key) => String(key).replace(RUNTIME_IDENTITY, 'runtime-identity');
  const compare = (a, b) => (a < b ? -1 : a > b ? 1 : 0);
  // Keys in order by key, then by value, each with its runtime identities
  // masked, so the order does not depend on the names a page gave.
  const canonical = (value) => {
    if (Array.isArray(value)) return value.map(canonical);
    if (!isObject(value)) return value;
    const entries = Object.keys(value).map((key) => [key, canonical(value[key])])
      .map(([key, field]) => ({ key, field, order: [masked(key), masked(JSON.stringify(field) ?? '')] }));
    entries.sort((a, b) => compare(a.order[0], b.order[0]) || compare(a.order[1], b.order[1]));
    return Object.fromEntries(entries.map(({ key, field }) => [key, field]));
  };
  // The keys at a pattern's place (`*` each key or index), by masked key.
  const keysAt = (value, segment) => (segment === '*' ? Object.keys(value) : Object.keys(value).filter((key) => masked(key) === segment));
  // The fields a pattern names where they hold one of `when` (or always),
  // each as its masked keys.
  const fieldsAt = (value, [head, ...rest], when, at = []) => {
    if (!isObject(value)) return [];
    const keys = keysAt(value, head);
    if (rest.length) return keys.flatMap((key) => fieldsAt(value[key], rest, when, [...at, masked(key)]));
    const present = keys.length ? keys : head === '*' ? [] : [head];
    return present.filter((key) => !when || when.includes(Object.hasOwn(value, key) ? value[key] : null)).map((key) => [...at, masked(key)]);
  };
  const withoutField = (value, [head, ...rest]) => {
    if (!isObject(value)) return value;
    const copy = Array.isArray(value) ? [...value] : { ...value };
    for (const key of keysAt(copy, head)) {
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
  const fieldsLeftOut = {};
  for (const drawing of state.drawingList()) {
    drawings[drawing.id] = {};
    renamed[drawing.id] = {};
    fieldsLeftOut[drawing.id] = {};
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
      const fields = leftOut ? leftOut[drawing.id]?.[name] || []
        : (omit[name] || []).flatMap(({ path, when }) => fieldsAt(captured, path.split('.'), when));
      fieldsLeftOut[drawing.id][name] = fields;
      const kept = fields.reduce((value, field) => withoutField(value, field), captured);
      renamed[drawing.id][name] = await hex(renamedText(JSON.stringify(sourcesRenamed(canonical(kept)))));
    }
  }
  return { active: state.activeDrawing().id, drawings, renamed, leftOut: fieldsLeftOut };
}, { omit, leftOut });

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
