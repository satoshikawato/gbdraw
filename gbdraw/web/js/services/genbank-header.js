// @ts-check
// The single JavaScript reader of GenBank record headers. Record discovery,
// LOSAT FASTA extraction, and match sequences use it so their record IDs are
// the ones gbdraw.io.genome reads with Biopython.

// Biopython's GenBankScanner reads a header keyword from the first 12 columns
// and its value from the rest of the line; a value continues on lines that
// start with 12 spaces. A record starts only at a line that starts with
// "LOCUS" and seven spaces.
const GENBANK_INDENT = 12;
const GENBANK_SPACER = ' '.repeat(GENBANK_INDENT);
const GENBANK_RECORD_START = 'LOCUS       ';
const EMBL_RECORD_START = 'ID   ';
const GENBANK_HEADER_ENDS = new Set(['FEATURES', 'ORIGIN', 'CONTIG', 'BASE COUNT', 'WGS', 'TSA', 'TLS']);
const GENBANK_LINEAGE_ROOTS = new Set([
  'Bacteria.', 'Archaea.', 'Eukaryota.', 'Unclassified.', 'Viruses.',
  'cellular organisms.', 'other sequences.', 'unclassified sequences.'
]);
const GENBANK_LOCUS_PATTERN = /^LOCUS {7} *(\S+)(?:\s+(\d+)\s+(?:bp|aa)\b)?/;

// Python's str.strip and str.rstrip for the ORGANISM lines and the names the
// loader reports: String.trim also removes U+FEFF, which Python keeps. The
// U+FEFF stays in the organism and lineage, as Biopython reads them, so the
// name fast path declines on it (BROWSER_EQUAL_TEXT in record-discovery.js;
// OV-1833).
/** @param {string} text */
export const pythonStrip = (text) => text.replace(/^[^\S\uFEFF]+|[^\S\uFEFF]+$/g, '');
/** @param {string} text */
const pythonStripEnd = (text) => text.replace(/[^\S\uFEFF]+$/, '');

// One scan of a GenBank record header with the record ID and organism rules of
// Biopython's _FeatureConsumer, which gbdraw.io.genome uses. Returns null for a
// chunk without a LOCUS line and throws for a layout it cannot decide, so the
// caller declines to the Python reader.
/**
 * @typedef {{
 *   locus: string,
 *   accession: string,
 *   version: string,
 *   recordId: string,
 *   organism: string,
 *   lineage: string[],
 *   recordLength: number | null,
 *   topology: string
 * }} GenBankHeader
 */

/**
 * @param {unknown} chunk
 * @returns {GenBankHeader | null}
 */
export const scanGenBankHeader = (chunk) => {
  const lines = String(chunk || '').split(/\r?\n/);
  const locusIndex = lines.findIndex((line) => line.startsWith('LOCUS'));
  if (locusIndex < 0) return null;
  const locusLine = lines[locusIndex].trimEnd();
  const locusMatch = locusLine.match(GENBANK_LOCUS_PATTERN);
  const locusTokens = locusLine.slice(GENBANK_INDENT).trim().split(/\s+/).filter(Boolean);
  if (!locusMatch || (locusTokens.length > 1 && !locusMatch[2])) {
    throw new Error('This LOCUS line needs the Python reader.');
  }
  /** @type {string[]} */
  const header = [];
  for (const rawLine of lines.slice(locusIndex + 1)) {
    const line = pythonStripEnd(rawLine);
    if (line === '//' || GENBANK_HEADER_ENDS.has(line.slice(0, GENBANK_INDENT).trimEnd())) break;
    if (line) header.push(line);
  }
  /** @type {string[]} */
  const accessions = [];
  let recordId = '';
  let version = '';
  /** @type {number | null} */
  let versionSuffix = null;
  let organism = '';
  let lineage = '';
  /** @param {string} text */
  const addAccessions = (text) => {
    const tokens = text.replace(/;/g, ' ').split(/\s+/).filter(Boolean);
    tokens.forEach((token) => { if (!accessions.includes(token)) accessions.push(token); });
    if (!recordId && tokens.length > 0) recordId = accessions[0];
  };
  /** @param {number} index */
  const continuation = (index) => header[index + 1]?.startsWith(GENBANK_SPACER);
  for (let index = 0; index < header.length; index += 1) {
    const keyword = header[index].slice(0, GENBANK_INDENT).trim();
    let data = header[index].slice(GENBANK_INDENT).trim();
    if (keyword === 'ACCESSION') {
      while (continuation(index)) data += ` ${header[++index].slice(GENBANK_INDENT)}`;
      addAccessions(data);
    } else if (keyword === 'VERSION') {
      version = data.replace(/ {2,}/g, ' ').split(' GI:')[0];
      const parts = version.split('.');
      if (parts.length === 2 && /^\d+$/.test(parts[1])) {
        addAccessions(parts[0]);
        versionSuffix = Number(parts[1]);
      } else if (version) {
        recordId = version;
      }
    } else if (keyword === 'ORGANISM') {
      organism = pythonStrip(header[index].slice(GENBANK_INDENT));
      lineage = '';
      while (continuation(index)) {
        const line = header[++index];
        const value = pythonStrip(line.slice(GENBANK_INDENT));
        if (lineage || line.includes(';') || GENBANK_LINEAGE_ROOTS.has(value)) lineage += ` ${value}`;
        else if (value !== '.') organism += ` ${value}`;
      }
    }
  }
  const locus = locusMatch[1];
  if (!recordId) recordId = locus;
  else if (!recordId.includes('.') && versionSuffix !== null) recordId += `.${versionSuffix}`;
  return {
    locus,
    accession: accessions[0] || '',
    version,
    recordId,
    organism: pythonStrip(organism),
    lineage: lineageWordsOf(lineage),
    recordLength: locusMatch[2] ? Number(locusMatch[2]) : null,
    topology: locusLine.match(/^LOCUS\s+\S+\s+\d+\s+(?:bp|aa)\b.*\s(circular|linear)(?:\s|$)/)?.[1] || 'unknown'
  };
};

// The taxa of the ORGANISM lineage, highest first, as Biopython's
// _split_taxonomy gives annotations["taxonomy"] and
// infer_record_source_metadata reads it (stripped, empty taxa left out).
/** @param {string} lineage @returns {string[]} */
const lineageWordsOf = (lineage) => {
  const text = pythonStrip(lineage);
  if (!text || text === '.') return [];
  return (text.endsWith('.') ? text.slice(0, -1) : text).split(';').map(pythonStrip).filter(Boolean);
};

// The record ID the loader gives one GenBank record chunk, with the header
// values it is built from.
/**
 * @param {unknown} chunk
 * @returns {Pick<GenBankHeader, 'locus' | 'accession' | 'version' | 'recordId'> | null}
 */
export const genbankHeaderIds = (chunk) => {
  const header = scanGenBankHeader(chunk);
  if (!header) return null;
  const { locus, accession, version, recordId } = header;
  return { locus, accession, version, recordId };
};

// The number of records the loader reads from GenBank/DDBJ or EMBL text. As
// gbdraw.io.input_text.flat_file_format decides, a first non-blank line that
// starts with "ID" and three spaces makes the text EMBL, read one record per such
// line; other text is GenBank, one record per line that starts with "LOCUS" and
// seven spaces (GENBANK_RECORD_START).
/** @param {unknown} text */
export const countGenBankRecords = (text) => {
  const source = String(text || '');
  const recordStart = /^.*\S.*$/m.exec(source)?.[0].startsWith(EMBL_RECORD_START)
    ? EMBL_RECORD_START
    : GENBANK_RECORD_START;
  let count = 0;
  let lineStart = 0;
  while (lineStart >= 0 && lineStart < source.length) {
    if (source.startsWith(recordStart, lineStart)) count += 1;
    const newline = source.indexOf('\n', lineStart);
    lineStart = newline < 0 ? -1 : newline + 1;
  }
  return count;
};
