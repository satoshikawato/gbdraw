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

const readRecordPayload = (pyodide, helperName, args) => {
  const helper = pyodide.globals.get(helperName);
  try {
    if (typeof helper !== 'function') throw new Error('Record discovery helper is unavailable.');
    const rawPayload = helper(...args);
    return normalizeSequenceRecords(JSON.parse(String(rawPayload || '{}')));
  } finally {
    helper?.destroy?.();
  }
};

const unlinkIfPresent = (pyodide, path) => {
  try {
    pyodide.FS.unlink(path);
  } catch (_error) {
    // The file may not have been staged when discovery fails early.
  }
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

  try {
    const staged = await writeFileToFs(file, temporaryPath);
    if (!staged) throw new Error('Could not stage the sequence file.');
    return readRecordPayload(pyodide, 'list_sequence_records', [temporaryPath, format]);
  } finally {
    unlinkIfPresent(pyodide, temporaryPath);
  }
};

export const discoverGffFastaRecords = async ({
  gffFile,
  fastaFile,
  pyodide,
  writeFileToFs,
  gffTemporaryPath,
  fastaTemporaryPath
}) => {
  if (!gffFile || !fastaFile) throw new Error('GFF3 and FASTA files are required.');
  if (!pyodide) throw new Error('Python environment is not ready.');
  if (typeof writeFileToFs !== 'function') throw new Error('File staging is unavailable.');
  if (!gffTemporaryPath || !fastaTemporaryPath) throw new Error('Temporary paths are required.');

  try {
    const gffStaged = await writeFileToFs(gffFile, gffTemporaryPath);
    if (!gffStaged) throw new Error('Could not stage the GFF3 file.');
    const fastaStaged = await writeFileToFs(fastaFile, fastaTemporaryPath);
    if (!fastaStaged) throw new Error('Could not stage the FASTA file.');
    return readRecordPayload(
      pyodide,
      'list_gff_fasta_records',
      [gffTemporaryPath, fastaTemporaryPath]
    );
  } finally {
    unlinkIfPresent(pyodide, gffTemporaryPath);
    unlinkIfPresent(pyodide, fastaTemporaryPath);
  }
};
