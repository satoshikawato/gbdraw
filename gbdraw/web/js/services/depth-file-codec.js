const DEPTH_FILE_SCHEMA = 1;
const DEPTH_FILE_ENCODING = 'gbdraw-depth-table-v1';

const DEPTH_COLUMNS = ['reference_name', 'position', 'depth'];
const RUN_START = 0;
const RUN_STEP = 1;
const RUN_COUNT = 2;
const RUN_DEPTHS = 3;

const assertDepthPayload = (payload) => {
  if (!payload || payload.schema !== DEPTH_FILE_SCHEMA || !Array.isArray(payload.records)) {
    throw new Error('Invalid encoded depth file.');
  }
};

const decodeHeader = (header) => {
  if (header === null || header === undefined) return null;
  if (!Array.isArray(header) || header.length !== DEPTH_COLUMNS.length) {
    throw new Error('Invalid encoded depth file header.');
  }
  return header.map((value) => String(value ?? '')).join('\t');
};

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

export const decodeDepthText = (payload) => {
  assertDepthPayload(payload);
  const lineEnding = payload.lineEnding === '\r\n' ? '\r\n' : '\n';
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

export const isEncodedDepthFileEntry = (entry) =>
  Boolean(entry) && entry.encoding === DEPTH_FILE_ENCODING && entry.data?.schema === DEPTH_FILE_SCHEMA;
