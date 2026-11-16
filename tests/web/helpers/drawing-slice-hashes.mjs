// The sha256 of each slice capture of each drawing (`DRAWING_SLICES` of
// tests/web/helpers/drawing-state.mjs, as tests/web/per-mode-drawings.test.mjs
// hashes them in node), and the comparisons the random walk
// (tests/web/helpers/live-random-walk.mjs) and the Gallery round-trip audit
// (tools/audit/gallery-roundtrip.audit.spec.js) make of two such sets: the one
// reader of a drawing's state by slice for both.

/**
 * The slice hashes of every drawing by drawing ID, and the shown drawing.
 * `renamed` hashes each capture with its runtime identities renamed in the
 * order they first appear in it: the file IDs (`file-<n>`,
 * services/history-files.js) and the Linear record and Circular File card UIDs
 * (`linear-seq-…`, `circular-source-…`, services/drawing/inputs.js), which a
 * page makes anew when it reads a Session. Two states that differ only by
 * those names have equal `renamed` hashes, so a check across a Session Load
 * compares them; a check within one page compares `drawings`.
 * @typedef {{ active: string, drawings: Record<string, Record<string, string>>, renamed: Record<string, Record<string, string>> }} SliceHashes
 */
/**
 * A drawing whose slices differ between two hash sets, in a phase of a check.
 * @typedef {{ phase: string, drawing: string, slices: string[] }} SliceChange
 */

/** @param {import('@playwright/test').Page} page @returns {Promise<SliceHashes>} */
export const sliceHashes = (page) => page.evaluate(async () => {
  const { state } = await import('/gbdraw/web/js/state.js');
  const { DRAWING_SLICES } = await import('/tests/web/helpers/drawing-state.mjs');
  const hex = async (text) => [...new Uint8Array(await crypto.subtle.digest('SHA-256', new TextEncoder().encode(text)))]
    .map((byte) => byte.toString(16).padStart(2, '0')).join('');
  const RUNTIME_IDENTITY = /(?:linear-seq|circular-source)-[0-9a-f]+(?:-[0-9a-f]+)*|file-\d+/g;
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
    for (const [name, slice] of Object.entries(DRAWING_SLICES)) {
      const text = JSON.stringify(slice.capture(drawing.slices[name]));
      drawings[drawing.id][name] = await hex(text);
      renamed[drawing.id][name] = await hex(renamedText(text));
    }
  }
  return { active: state.activeDrawing().id, drawings, renamed };
});

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
