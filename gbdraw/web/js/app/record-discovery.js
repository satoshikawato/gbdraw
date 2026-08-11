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

const parseGenBankRecordText = (text) => {
  const records = String(text || '')
    .split(/^\/\/\s*$/m)
    .map((chunk, index) => {
      const locus = chunk.match(/^LOCUS\s+(\S+)(?:\s+(\d+)\s+(?:bp|aa)\b)?/m);
      if (!locus) return null;
      const accession = chunk.match(/^ACCESSION\s+(\S+)/m)?.[1];
      const version = chunk.match(/^VERSION\s+(\S+)/m)?.[1];
      return {
        selector: `#${index + 1}`,
        record_id: version || accession || locus[1],
        record_length: locus[2] ? Number(locus[2]) : null
      };
    })
    .filter(Boolean)
    .map((record, index) => ({ ...record, selector: `#${index + 1}` }));
  return normalizeSequenceRecords({ records });
};

const parseFastaRecordText = (text) => {
  const records = [];
  let current = null;
  String(text || '').split(/\r?\n/).forEach((line) => {
    if (line.startsWith('>')) {
      if (current) records.push(current);
      current = {
        selector: `#${records.length + 1}`,
        record_id: line.slice(1).trim().split(/\s+/)[0],
        record_length: 0
      };
    } else if (current) {
      current.record_length += line.replace(/\s+/g, '').length;
    }
  });
  if (current) records.push(current);
  return normalizeSequenceRecords({ records });
};

export const parseSequenceRecordText = (text, format) => {
  if (format === 'genbank') return parseGenBankRecordText(text);
  if (format === 'fasta') return parseFastaRecordText(text);
  throw new Error(`Unsupported format: ${String(format)}.`);
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
  readText = null,
  pyodide,
  writeFileToFs,
  temporaryPath
}) => {
  if (!file) throw new Error('A sequence file is required.');
  let textError = null;
  const readSourceText = typeof readText === 'function'
    ? () => readText(file)
    : (typeof file.text === 'function' ? () => file.text() : null);
  if (readSourceText) {
    try {
      return parseSequenceRecordText(await readSourceText(), format);
    } catch (error) {
      textError = error;
    }
  }
  if (!pyodide) throw textError || new Error('Python environment is not ready.');
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
  readText = null,
  pyodide,
  writeFileToFs,
  gffTemporaryPath,
  fastaTemporaryPath
}) => {
  if (!gffFile || !fastaFile) throw new Error('GFF3 and FASTA files are required.');
  let textError = null;
  const readSourceText = typeof readText === 'function'
    ? () => readText(fastaFile)
    : (typeof fastaFile.text === 'function' ? () => fastaFile.text() : null);
  if (readSourceText) {
    try {
      return parseSequenceRecordText(await readSourceText(), 'fasta');
    } catch (error) {
      textError = error;
    }
  }
  if (!pyodide) throw textError || new Error('Python environment is not ready.');
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
