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
const GENBANK_HEADER_ENDS = new Set(['FEATURES', 'ORIGIN', 'CONTIG', 'BASE COUNT', 'WGS', 'TSA', 'TLS']);
const GENBANK_LINEAGE_ROOTS = new Set([
  'Bacteria.', 'Archaea.', 'Eukaryota.', 'Unclassified.', 'Viruses.',
  'cellular organisms.', 'other sequences.', 'unclassified sequences.'
]);
const GENBANK_LOCUS_PATTERN = /^LOCUS {7} *(\S+)(?:\s+(\d+)\s+(?:bp|aa)\b)?/;

// One scan of a GenBank record header with the record ID and organism rules of
// Biopython's _FeatureConsumer, which gbdraw.io.genome uses. Returns null for a
// chunk without a LOCUS line and throws for a layout it cannot decide, so the
// caller declines to the Python reader.
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
  const header = [];
  for (const rawLine of lines.slice(locusIndex + 1)) {
    const line = rawLine.trimEnd();
    if (line === '//' || GENBANK_HEADER_ENDS.has(line.slice(0, GENBANK_INDENT).trimEnd())) break;
    if (line) header.push(line);
  }
  const accessions = [];
  let recordId = '';
  let version = '';
  let versionSuffix = null;
  let organism = '';
  const addAccessions = (text) => {
    const tokens = text.replace(/;/g, ' ').split(/\s+/).filter(Boolean);
    tokens.forEach((token) => { if (!accessions.includes(token)) accessions.push(token); });
    if (!recordId && tokens.length > 0) recordId = accessions[0];
  };
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
      organism = data;
      let lineage = '';
      while (continuation(index)) {
        const line = header[++index];
        const value = line.slice(GENBANK_INDENT).trim();
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
    organism: organism.trim(),
    recordLength: locusMatch[2] ? Number(locusMatch[2]) : null,
    topology: locusLine.match(/^LOCUS\s+\S+\s+\d+\s+(?:bp|aa)\b.*\s(circular|linear)(?:\s|$)/)?.[1] || 'unknown'
  };
};

// The record ID the loader gives one GenBank record chunk, with the header
// values it is built from.
export const genbankHeaderIds = (chunk) => {
  const header = scanGenBankHeader(chunk);
  if (!header) return null;
  const { locus, accession, version, recordId } = header;
  return { locus, accession, version, recordId };
};
