const normalizeRecordLength = (value) => {
  const numeric = Number(value);
  return Number.isInteger(numeric) && numeric > 0 ? numeric : null;
};

export const normalizeSequenceRecords = (payload) => {
  if (payload?.error) throw new Error(String(payload.error));
  if (!Array.isArray(payload?.records)) throw new Error('Record list response is invalid.');

  const records = [];
  const seenSelectors = new Set();
  payload.records.forEach((entry, index) => {
    const selector = String(entry?.selector ?? `#${index + 1}`).trim();
    if (!selector || seenSelectors.has(selector)) return;
    seenSelectors.add(selector);
    records.push({
      selector,
      recordId: String(entry?.record_id ?? '').trim() || `Record_${index + 1}`,
      recordLength: normalizeRecordLength(entry?.record_length)
    });
  });

  if (records.length === 0) throw new Error('No records found.');
  return records;
};

export const discoverSequenceRecords = async ({
  file,
  format,
  pyodide,
  writeFileToFs,
  temporaryPath
}) => {
  if (!file) throw new Error('A sequence file is required.');
  if (!pyodide) throw new Error('Python environment is not ready.');
  if (typeof writeFileToFs !== 'function') throw new Error('File staging is unavailable.');
  if (!temporaryPath) throw new Error('A temporary path is required.');

  let listRecords = null;
  try {
    const staged = await writeFileToFs(file, temporaryPath);
    if (!staged) throw new Error('Could not stage the sequence file.');
    listRecords = pyodide.globals.get('list_sequence_records');
    if (typeof listRecords !== 'function') throw new Error('Record discovery helper is unavailable.');
    const rawPayload = listRecords(temporaryPath, format);
    const payload = JSON.parse(String(rawPayload || '{}'));
    return normalizeSequenceRecords(payload);
  } finally {
    listRecords?.destroy?.();
    try {
      pyodide.FS.unlink(temporaryPath);
    } catch (_error) {
      // The file may not have been staged when discovery fails early.
    }
  }
};
