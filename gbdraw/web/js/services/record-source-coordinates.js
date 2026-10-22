// @ts-check
// Input-file coordinates of displayed Linear records (PD-OI-076). One
// dependency-free formula for the match popup, FASTA headers and the
// Interactive SVG payload.

/**
 * The input-file span of a displayed record: `start`..`end` in input-file
 * coordinates, read in the direction `step` (1 or -1).
 * @typedef {{ start: number, end: number, step: 1 | -1 }} RecordSourceSpan
 */

/**
 * An element of the displayed SVG; an SVG element names its root as `ownerSVGElement`.
 * @typedef {Element & { ownerSVGElement?: SVGSVGElement | null }} RecordSourceElement
 */

const RECORD_SOURCE_ATTRIBUTES = ['start', 'end', 'step']
  .map((name) => `data-gbdraw-record-source-${name}`);

/**
 * The marked record groups of one SVG root by record index, read with one query
 * the first time the root is asked (OV-931). The renderer writes these
 * attributes once per record group and nothing edits them on a mounted root.
 * @type {WeakMap<Element, Map<string, Element>>}
 */
const markedRecordGroupsByRoot = new WeakMap();

/**
 * @param {Element} root
 * @returns {Map<string, Element>}
 */
const markedRecordGroups = (root) => {
  let groups = markedRecordGroupsByRoot.get(root);
  if (groups) return groups;
  groups = new Map();
  for (const group of root.querySelectorAll(`[data-gbdraw-record-index][${RECORD_SOURCE_ATTRIBUTES[2]}]`)) {
    const index = group.getAttribute('data-gbdraw-record-index') || '';
    // The first group in document order answers, as `querySelector` did.
    if (!groups.has(index)) groups.set(index, group);
  }
  markedRecordGroupsByRoot.set(root, groups);
  return groups;
};

/**
 * Read the input-file span of one displayed Linear record (PD-OI-076).
 *
 * The renderer marks only cropped or reverse-complemented records; any other
 * record shows its input file one to one and returns null.
 * @param {RecordSourceElement | null | undefined} element
 * @param {unknown} recordIndex
 * @returns {RecordSourceSpan | null}
 */
export const readRecordSourceSpan = (element, recordIndex) => {
  const root = element?.ownerSVGElement || element?.closest?.('svg') || null;
  const index = Number(recordIndex);
  if (!root?.querySelectorAll || !Number.isInteger(index) || index < 0) return null;
  const group = markedRecordGroups(root).get(String(index));
  if (!group) return null;
  const [start, end, step] = RECORD_SOURCE_ATTRIBUTES.map((name) => Number(group.getAttribute(name)));
  if (!Number.isInteger(start) || !Number.isInteger(end) || start < 1 || end < start) return null;
  return step === 1 || step === -1 ? { start, end, step } : null;
};

/**
 * Map record-local match coordinates to source coordinates, the one formula
 * for popups, FASTA headers and the Interactive SVG. `table` is the search
 * frame of comparison tables (the cropped record in the source strand) and is
 * returned only when it differs from the source coordinates.
 * @param {RecordSourceSpan | null | undefined} sourceSpan
 * @param {unknown} startRaw
 * @param {unknown} endRaw
 * @returns {{ start: number, end: number, table: { start: number, end: number } | null } | null}
 */
export const recordSourceInterval = (sourceSpan, startRaw, endRaw) => {
  const local = [Number(startRaw), Number(endRaw)];
  if (!sourceSpan || !local.every((value) => Number.isInteger(value) && value >= 1)) return null;
  const toSource = (/** @type {number} */ value) => (sourceSpan.step === 1
    ? sourceSpan.start + value - 1
    : sourceSpan.end - value + 1);
  const [start, end] = local.map(toSource);
  const table = { start: start - sourceSpan.start + 1, end: end - sourceSpan.start + 1 };
  return {
    start,
    end,
    table: table.start === start && table.end === end ? null : table
  };
};
