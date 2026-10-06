// @ts-check
// Input-file coordinates of displayed Linear records (PD-OI-076). One
// dependency-free formula for the match popup, FASTA headers and the
// Interactive SVG payload.

const RECORD_SOURCE_ATTRIBUTES = ['start', 'end', 'step']
  .map((name) => `data-gbdraw-record-source-${name}`);

/**
 * Read the input-file span of one displayed Linear record (PD-OI-076).
 *
 * The renderer marks only cropped or reverse-complemented records; any other
 * record shows its input file one to one and returns null.
 */
export const readRecordSourceSpan = (element, recordIndex) => {
  const root = element?.ownerSVGElement || element?.closest?.('svg') || null;
  const index = Number(recordIndex);
  if (!root?.querySelector || !Number.isInteger(index) || index < 0) return null;
  const group = root.querySelector(
    `[data-gbdraw-record-index="${index}"][${RECORD_SOURCE_ATTRIBUTES[2]}]`
  );
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
 */
export const recordSourceInterval = (sourceSpan, startRaw, endRaw) => {
  const local = [Number(startRaw), Number(endRaw)];
  if (!sourceSpan || !local.every((value) => Number.isInteger(value) && value >= 1)) return null;
  const toSource = (value) => (sourceSpan.step === 1
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
