const DEPTH_FILE_SCHEMA = 1;
export const DEPTH_FILE_ENCODING = 'gbdraw-depth-table-v1';

const DEPTH_COLUMNS = ['reference_name', 'position', 'depth'];
const RUN_START = 0;
const RUN_STEP = 1;
const RUN_COUNT = 2;
const RUN_DEPTHS = 3;

const isIntegerText = (value) => /^[+-]?\d+$/.test(String(value || '').trim());

const parsePositiveSafeInteger = (value) => {
  const trimmed = String(value || '').trim();
  if (!isIntegerText(trimmed)) return null;
  const parsed = Number(trimmed);
  if (!Number.isSafeInteger(parsed) || parsed <= 0) return null;
  return parsed;
};

const isDepthText = (value) => {
  const trimmed = String(value || '').trim();
  if (!trimmed) return false;
  const parsed = Number(trimmed);
  return Number.isFinite(parsed) && parsed >= 0;
};

const hasDepthHeader = (fields) =>
  fields.length >= 3 &&
  (parsePositiveSafeInteger(fields[1]) === null || !isDepthText(fields[2]));

const makeRun = (position, depthValue) => [position, 1, 1, [depthValue]];

const appendDepthRow = (records, referenceName, position, depthValue) => {
  let record = records[records.length - 1];
  if (!record || record.id !== referenceName) {
    record = { id: referenceName, runs: [] };
    records.push(record);
  }

  let run = record.runs[record.runs.length - 1];
  if (!run) {
    record.runs.push(makeRun(position, depthValue));
    return;
  }

  const start = run[RUN_START];
  const count = run[RUN_COUNT];
  const step = run[RUN_STEP];

  if (count === 1) {
    const nextStep = position - start;
    if (nextStep > 0) {
      run[RUN_STEP] = nextStep;
      run[RUN_COUNT] = 2;
      run[RUN_DEPTHS].push(depthValue);
      return;
    }
    record.runs.push(makeRun(position, depthValue));
    return;
  }

  if (position === start + step * count) {
    run[RUN_COUNT] = count + 1;
    run[RUN_DEPTHS].push(depthValue);
    return;
  }

  record.runs.push(makeRun(position, depthValue));
};

export const encodeDepthText = (text) => {
  if (typeof text !== 'string' || !text) return null;
  const crlfCount = text.match(/\r\n/g)?.length || 0;
  const lfCount = text.match(/\n/g)?.length || 0;
  if (crlfCount > 0 && crlfCount !== lfCount) return null;

  const lineEnding = text.includes('\r\n') ? '\r\n' : '\n';
  const normalized = text.replace(/\r\n/g, '\n');
  if (normalized.includes('\r')) return null;
  const finalNewline = normalized.endsWith('\n');
  const lines = normalized.split('\n');
  if (finalNewline) lines.pop();
  if (lines.length === 0 || lines.some((line) => line === '')) return null;

  let lineIndex = 0;
  let header = null;
  const firstFields = lines[0].split('\t');
  if (firstFields.length !== DEPTH_COLUMNS.length) return null;

  if (hasDepthHeader(firstFields)) {
    header = firstFields;
    lineIndex = 1;
  }
  if (lineIndex >= lines.length) return null;

  const records = [];
  let rowCount = 0;
  for (; lineIndex < lines.length; lineIndex += 1) {
    const fields = lines[lineIndex].split('\t');
    if (fields.length !== DEPTH_COLUMNS.length) return null;

    const position = parsePositiveSafeInteger(fields[1]);
    const depthValue = String(fields[2] || '').trim();
    if (position === null || !isDepthText(depthValue)) return null;

    appendDepthRow(records, String(fields[0] ?? ''), position, depthValue);
    rowCount += 1;
  }

  if (rowCount === 0) return null;

  return {
    schema: DEPTH_FILE_SCHEMA,
    columns: DEPTH_COLUMNS,
    lineEnding,
    finalNewline,
    rowCount,
    header,
    records
  };
};

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
