// @ts-check
import { DEPTH_FILE_ENCODING, DEPTH_FILE_SCHEMA } from '../contracts/schema.generated.js';

const DEPTH_COLUMNS = ['reference_name', 'position', 'depth'];
const RUN_START = 0;
const RUN_STEP = 1;
const RUN_COUNT = 2;
const RUN_DEPTHS = 3;

/**
 * The encoded depth table: runs of `[start, step, count, depths]` per reference. Every field is validated on decode.
 * @typedef {{
 *   schema?: number,
 *   header?: unknown,
 *   lineEnding?: string,
 *   finalNewline?: boolean,
 *   records?: { id?: unknown, runs?: unknown }[]
 * }} DepthFilePayload
 */

/** @typedef {DepthFilePayload & { schema: 1, records: NonNullable<DepthFilePayload['records']> }} ValidDepthFilePayload */

/** @type {(payload: DepthFilePayload | null | undefined) => asserts payload is ValidDepthFilePayload} */
const assertDepthPayload = (payload) => {
  if (!payload || payload.schema !== DEPTH_FILE_SCHEMA || !Array.isArray(payload.records)) {
    throw new Error('Invalid encoded depth file.');
  }
};

/** @param {unknown} header */
const decodeHeader = (header) => {
  if (header === null || header === undefined) return null;
  if (!Array.isArray(header) || header.length !== DEPTH_COLUMNS.length) {
    throw new Error('Invalid encoded depth file header.');
  }
  return header.map((value) => String(value ?? '')).join('\t');
};

/**
 * @param {string} referenceName
 * @param {unknown} run
 * @param {string[]} lines
 */
const decodeRun = (referenceName, run, lines) => {
  if (!Array.isArray(run) || run.length !== 4 || !Array.isArray(run[RUN_DEPTHS])) {
    throw new Error('Invalid encoded depth run.');
  }

  const start = Number(run[RUN_START]);
  const step = Number(run[RUN_STEP]);
  const count = Number(run[RUN_COUNT]);
  const depths = run[RUN_DEPTHS];
  if (
    !Number.isSafeInteger(start) ||
    !Number.isSafeInteger(step) ||
    !Number.isSafeInteger(count) ||
    start <= 0 ||
    step <= 0 ||
    count <= 0 ||
    depths.length !== count
  ) {
    throw new Error('Invalid encoded depth coordinates.');
  }

  for (let idx = 0; idx < count; idx += 1) {
    const position = start + step * idx;
    if (!Number.isSafeInteger(position)) {
      throw new Error('Invalid encoded depth coordinate overflow.');
    }
    lines.push(`${referenceName}\t${position}\t${String(depths[idx] ?? '')}`);
  }
};

/** @param {DepthFilePayload | null | undefined} payload */
export const decodeDepthText = (payload) => {
  assertDepthPayload(payload);
  const lineEnding = payload.lineEnding === '\r\n' ? '\r\n' : '\n';
  /** @type {string[]} */
  const lines = [];
  const headerLine = decodeHeader(payload.header);
  if (headerLine !== null) lines.push(headerLine);

  payload.records.forEach((record) => {
    if (!record || !Array.isArray(record.runs)) {
      throw new Error('Invalid encoded depth record.');
    }
    const referenceName = String(record.id ?? '');
    record.runs.forEach((run) => decodeRun(referenceName, run, lines));
  });

  if (lines.length === 0) return '';
  const body = lines.join(lineEnding);
  return payload.finalNewline === false ? body : `${body}${lineEnding}`;
};

/** @param {{ encoding?: unknown, data?: { schema?: unknown } | null } | null | undefined} entry */
export const isEncodedDepthFileEntry = (entry) =>
  entry?.encoding === DEPTH_FILE_ENCODING && entry.data?.schema === DEPTH_FILE_SCHEMA;
