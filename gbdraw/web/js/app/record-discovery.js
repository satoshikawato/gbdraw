import {
  DIAGRAM_HELPER_OPERATIONS,
  runDiagramHelperOperation
} from '../services/diagram-generation.js';
import { cloneFileBytesForTransfer } from '../services/file-content-cache.js';

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

export const discoverSequenceRecords = async ({
  file,
  format,
  readText = null,
  runHelperOperation = runDiagramHelperOperation
}) => {
  if (!file) throw new Error('A sequence file is required.');
  const readSourceText = typeof readText === 'function'
    ? () => readText(file)
    : (typeof file.text === 'function' ? () => file.text() : null);
  if (readSourceText) {
    try {
      return parseSequenceRecordText(await readSourceText(), format);
    } catch {
      // The packaged Worker parser handles variants beyond the lightweight text fast path.
    }
  }
  const response = await runHelperOperation(
    DIAGRAM_HELPER_OPERATIONS.LIST_SEQUENCE_RECORDS,
    {
      format,
      files: [{ role: 'source', bytes: await cloneFileBytesForTransfer(file) }]
    }
  );
  return normalizeSequenceRecords(response.result);
};

export const discoverGffFastaRecords = async ({
  gffFile,
  fastaFile,
  readText = null,
  runHelperOperation = runDiagramHelperOperation
}) => {
  if (!gffFile || !fastaFile) throw new Error('GFF3 and FASTA files are required.');
  const readSourceText = typeof readText === 'function'
    ? () => readText(fastaFile)
    : (typeof fastaFile.text === 'function' ? () => fastaFile.text() : null);
  if (readSourceText) {
    try {
      return parseSequenceRecordText(await readSourceText(), 'fasta');
    } catch {
      // Let the Worker validate the paired GFF3/FASTA record set together.
    }
  }
  const response = await runHelperOperation(
    DIAGRAM_HELPER_OPERATIONS.LIST_GFF_FASTA_RECORDS,
    {
      files: [
        { role: 'gff', bytes: await cloneFileBytesForTransfer(gffFile) },
        { role: 'fasta', bytes: await cloneFileBytesForTransfer(fastaFile) }
      ]
    }
  );
  return normalizeSequenceRecords(response.result);
};
