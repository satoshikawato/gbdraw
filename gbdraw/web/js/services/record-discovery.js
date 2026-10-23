// @ts-check
import {
  DIAGRAM_HELPER_OPERATIONS,
  runDiagramHelperOperation
} from './diagram-generation.js';
import { scanGenBankHeader } from './genbank-header.js';
import { sha256Hex, textToBytes } from './byte-utils.js';
import {
  bytesToText,
  cloneFileBytesForTransfer,
  getSessionResourceSource,
  readFileBytes,
  readFileText
} from './file-content-cache.js';

/** @import { FileLike } from './session-resource-backing.js' */
/** @import { GenBankHeader } from './genbank-header.js' */
/** @import { UserFacingError } from '../utils/error-normalization.js' */

/**
 * The record catalog of the Circular source files: whose files it was read
 * from (`inputType`, `primaryFile`, `pairedFile`), and how reading went.
 * @typedef {object} CircularRecordDiscovery
 * @property {'idle' | 'loading' | 'deferred' | 'ready' | 'error'} status
 * @property {UserFacingError | ''} error The reader's error, or '' when there is none.
 * @property {'' | 'gb' | 'gff'} inputType
 * @property {FileLike | null} primaryFile
 * @property {FileLike | null} pairedFile
 * @property {{ selector: string, record_id: string, recordKey: string }[]} canonicalRecordIdentities
 */

/**
 * One record as a reader or the Worker reports it: Python names (`record_id`,
 * `record_length`, `topology`) or the names `normalizeSequenceRecords` returns.
 * @typedef {{
 *   selector?: unknown, record_id?: unknown, recordId?: unknown,
 *   record_length?: unknown, recordLength?: unknown,
 *   topology?: string, detectedTopology?: string,
 *   organism?: unknown, strain?: unknown,
 *   inferredDefinition?: unknown, inferred_definition?: unknown
 * }} SequenceRecordEntry
 */

/**
 * An unvalidated record-list response: `normalizeSequenceRecords` checks it.
 * @typedef {{ error?: unknown, records?: (SequenceRecordEntry | null | undefined)[] }} SequenceRecordsPayload
 */

/**
 * @typedef {{
 *   selector: string, recordId: string, recordLength: number | null, detectedTopology: string,
 *   organism?: string, strain?: string, inferredDefinition?: string
 * }} SequenceRecord
 */

/** @typedef {{ role: string, bytes: ArrayBuffer }} HelperFile */

/**
 * The Worker helper call: `runDiagramHelperOperation` answers it.
 * @template Result
 * @typedef {(
 *   operation: string, payload: { format?: string, files: HelperFile[] }
 * ) => Promise<{ result?: Result }>} HelperOperationRunner
 */

// A discovery error the reader reported for the source bytes is final for that
// exact source instance: reading the same files again returns it again, so a
// refresh or Generate reuses it. A Worker start-up, staging, transport, Cancel,
// or unclassified (UNKNOWN) failure is not a property of the bytes and reads again.
/** @param {UserFacingError | string | null | undefined} error The stored error, or '' when there is none. */
export const discoveryErrorIsFinal = (error) => typeof error === 'object' && error !== null
  && Boolean(error.code) && error.code !== 'UNKNOWN' && error.stage === 'helper';

/**
 * The members of the `state.js` state that `circularDiscoveryForInput` reads;
 * the record list passes through.
 * @template Entry
 * @typedef {{
 *   cInputType: { value: string },
 *   files: { c_gb: FileLike | null, c_gff: FileLike | null, c_fasta: FileLike | null },
 *   circularRecordDiscovery: CircularRecordDiscovery,
 *   circularRecordList: { value: Entry[] }
 * }} CircularDiscoveryState
 */

/**
 * A catalog is current only for the exact active source instance and input type.
 * @template Entry
 * @param {CircularDiscoveryState<Entry>} state
 */
export const circularDiscoveryForInput = (state) => {
  const inputType = state.cInputType.value;
  const primaryFile = inputType === 'gff' ? state.files.c_gff : state.files.c_gb;
  const pairedFile = inputType === 'gff' ? state.files.c_fasta : null;
  const hasInput = Boolean(primaryFile && (inputType !== 'gff' || pairedFile));
  const discovery = state.circularRecordDiscovery;
  const current = discovery.inputType === inputType
    && discovery.primaryFile === primaryFile && discovery.pairedFile === pairedFile;
  return {
    hasInput, current, primaryFile, pairedFile, inputType,
    status: hasInput ? (current ? discovery.status : 'deferred') : 'idle',
    error: current ? discovery.error : '',
    records: current && discovery.status === 'ready' ? state.circularRecordList.value : []
  };
};

/** @param {unknown} value */
const normalizeRecordLength = (value) => {
  const numeric = Number(value);
  return Number.isInteger(numeric) && numeric > 0 ? numeric : null;
};

const NON_ORGANISM_PATTERN = /^(?:synthetic construct|artificial sequence|unidentified(?: organism)?|unknown(?: organism)?|vector|(?:unidentified )?cloning vector|expression vector)$/i;

// This fast path mirrors gbdraw/core/record_metadata.py, which owns the same
// inference for the Worker. tests/web/record-metadata-inference.test.mjs and
// tests/test_record_metadata.py share one fixture table so the two cannot drift.
/** @param {{ organism?: unknown, strain?: unknown }} [fields] */
export const formatInferredOrganismStrain = ({ organism = '', strain = '' } = {}) => {
  const rawOrganism = String(organism || '').trim();
  const rawStrain = String(strain || '').trim();

  const isCandidatus = /^Candidatus\s+/i.test(rawOrganism);
  const nameWithoutCand = isCandidatus ? rawOrganism.replace(/^Candidatus\s+/i, '').trim() : rawOrganism;
  if (!nameWithoutCand || NON_ORGANISM_PATTERN.test(nameWithoutCand)) {
    return rawStrain;
  }

  const words = nameWithoutCand.split(/\s+/).filter(Boolean);

  let speciesPart = '';
  let rest = '';

  if (words.length >= 2) {
    const binomial = `${words[0]} ${words[1]}`;
    speciesPart = `${isCandidatus ? 'Candidatus ' : ''}<i>${binomial}</i>`;
    rest = words.slice(2).join(' ');
  } else if (words.length === 1) {
    speciesPart = `${isCandidatus ? 'Candidatus ' : ''}<i>${words[0]}</i>`;
  } else {
    speciesPart = rawOrganism;
  }

  if (rawStrain && !rest.toLowerCase().includes(rawStrain.toLowerCase())) {
    rest = rest ? `${rest} ${rawStrain}` : rawStrain;
  }

  return rest ? `${speciesPart} ${rest}`.trim() : speciesPart.trim();
};

/** @param {string} subject @returns {never} */
const declineToPythonReader = (subject) => {
  throw new Error(`${subject} needs the Python reader.`);
};

/**
 * @param {unknown} chunk
 * @param {GenBankHeader | null} [header]
 */
export const extractGenBankMetadata = (chunk, header = undefined) => {
  const text = String(chunk || '');
  const sourceMatch = text.match(/^ {5}source\s+[\s\S]*?(?=^ {5}[a-z]|\/\/)/mi);
  const sourceBlock = sourceMatch ? sourceMatch[0] : text;

  /** @param {string} key @param {string} block */
  const extractQualifier = (key, block) => {
    const regex = new RegExp(`/${key}="([^"]*(?:\\r?\\n {21}[^"]*)*)"`, 'i');
    const match = block.match(regex);
    if (!match) return '';
    return match[1].replace(/\r?\n\s+/g, ' ').trim();
  };

  // infer_record_source_metadata prefers the source /organism to ORGANISM.
  const organism = extractQualifier('organism', sourceBlock)
    || (header === undefined ? scanGenBankHeader(text) : header)?.organism || '';
  const strain = extractQualifier('strain', sourceBlock);
  const isolate = extractQualifier('isolate', sourceBlock);

  // infer_record_source_metadata reads /isolate before /strain.
  const effectiveStrain = isolate || strain || '';
  const inferredDefinition = formatInferredOrganismStrain({
    organism,
    strain: effectiveStrain
  });
  return { organism, strain: effectiveStrain, inferredDefinition };
};

/**
 * @param {SequenceRecordsPayload | null | undefined} payload
 * @returns {SequenceRecord[]}
 */
export const normalizeSequenceRecords = (payload) => {
  if (payload?.error) throw payload.error;
  if (!Array.isArray(payload?.records)) throw new Error('Record list response is invalid.');

  /** @type {SequenceRecord[]} */
  const records = [];
  const seenSelectors = new Set();
  payload.records.forEach((entry, index) => {
    const selector = String(entry?.selector ?? `#${index + 1}`).trim();
    if (!selector || seenSelectors.has(selector)) return;
    seenSelectors.add(selector);
    const topology = entry?.topology ?? entry?.detectedTopology;
    /** @type {SequenceRecord} */
    const record = {
      selector,
      recordId: String(entry?.record_id ?? entry?.recordId ?? '').trim() || `Record_${index + 1}`,
      recordLength: normalizeRecordLength(entry?.record_length ?? entry?.recordLength),
      detectedTopology: topology === 'circular' || topology === 'linear' ? topology : 'unknown'
    };
    const organism = String(entry?.organism ?? '').trim();
    const strain = String(entry?.strain ?? '').trim();
    const inferredDefinition = String(entry?.inferredDefinition ?? entry?.inferred_definition ?? '').trim();
    if (organism) record.organism = organism;
    if (strain) record.strain = strain;
    if (inferredDefinition) record.inferredDefinition = inferredDefinition;
    records.push(record);
  });

  if (records.length === 0) throw new Error('No records found.');
  return records;
};

/** @param {unknown} text */
const parseGenBankRecordText = (text) => {
  const records = String(text || '')
    .split(/^\/\/\s*$/m)
    .map((chunk) => {
      const header = scanGenBankHeader(chunk);
      if (!header) return null;
      const metadata = extractGenBankMetadata(chunk, header);
      return {
        record_id: header.recordId,
        record_length: header.recordLength,
        topology: header.topology,
        organism: metadata.organism,
        strain: metadata.strain,
        inferredDefinition: metadata.inferredDefinition
      };
    })
    .filter(Boolean)
    .map((record, index) => ({ ...record, selector: `#${index + 1}` }));
  return normalizeSequenceRecords({ records });
};

// Sequence entries with Biopython's FASTA rules; text before the first record
// or an entry without an ID is left to the Python reader.
/** @param {string[]} lines */
const fastaEntries = (lines) => {
  /** @type {{ id: string, length: number }[]} */
  const entries = [];
  /** @type {{ id: string, length: number } | null} */
  let current = null;
  lines.forEach((line) => {
    if (line.startsWith('>')) {
      const id = line.slice(1).trim().split(/\s+/)[0];
      if (!id) declineToPythonReader('A FASTA entry without an ID');
      current = { id, length: 0 };
      entries.push(current);
    } else if (current) {
      current.length += line.replace(/\s+/g, '').length;
    } else if (line.trim()) {
      declineToPythonReader('Text before the first FASTA entry');
    }
  });
  return entries;
};

/** @param {unknown} text */
const parseFastaRecordText = (text) => normalizeSequenceRecords({
  records: fastaEntries(String(text || '').split(/\r?\n/)).map((entry, index) => ({
    selector: `#${index + 1}`,
    record_id: entry.id,
    record_length: entry.length
  }))
});

// The record IDs BCBio builds from a GFF3 file: the seqid of every feature line
// before ##FASTA and every embedded FASTA entry. Lines BCBio would reject or
// match by NCBI-style IDs are left to the Python reader.
/** @param {unknown} text */
const gffRecordIds = (text) => {
  const lines = String(text || '').split(/\r?\n/);
  /** @type {Set<string>} */
  const ids = new Set();
  let index = 0;
  for (; index < lines.length; index += 1) {
    const line = lines[index].trim();
    if (line.startsWith('##')) {
      if (line.slice(2) === 'FASTA') {
        index += 1;
        break;
      }
      continue;
    }
    if (!line || line.startsWith('#')) continue;
    const columns = line.split('\t');
    const [seqid, , , start, end] = columns;
    if (columns.length < 8 || seqid === '.' || seqid.includes('|')
      || !/^\d+$/.test(start) || !/^\d+$/.test(end)) {
      declineToPythonReader('This GFF3 line');
    }
    ids.add(seqid);
  }
  fastaEntries(lines.slice(index)).forEach(({ id }) => {
    if (id.includes('|')) declineToPythonReader('An NCBI-style embedded FASTA ID');
    ids.add(id);
  });
  return ids;
};

// load_gff_fasta keeps the GFF3 records the FASTA file names, in FASTA order.
// A GFF3 record missing from the FASTA file is the loader's error to report.
/** @param {unknown} gffText @param {unknown} fastaText */
const parseGffFastaRecordText = (gffText, fastaText) => {
  const fasta = fastaEntries(String(fastaText || '').split(/\r?\n/));
  const fastaIds = new Set(fasta.map(({ id }) => id));
  if (fastaIds.size !== fasta.length) declineToPythonReader('A FASTA file with repeated IDs');
  const gffIds = gffRecordIds(gffText);
  gffIds.forEach((id) => {
    if (!fastaIds.has(id)) declineToPythonReader('A GFF3 record without a FASTA entry');
  });
  return normalizeSequenceRecords({
    records: fasta.filter(({ id }) => gffIds.has(id)).map((entry, index) => ({
      selector: `#${index + 1}`,
      record_id: entry.id,
      record_length: entry.length
    }))
  });
};

/** @param {unknown} text @param {unknown} format */
export const parseSequenceRecordText = (text, format) => {
  if (format === 'genbank') return parseGenBankRecordText(text);
  if (format === 'fasta') return parseFastaRecordText(text);
  throw new Error(`Unsupported format: ${String(format)}.`);
};

/** @param {Uint8Array} bytes */
const startsWithUtf8ByteOrderMark = (bytes) => (
  bytes[0] === 0xef && bytes[1] === 0xbb && bytes[2] === 0xbf
);

// The text the no-Worker fast path reads, or null to decline to the Worker.
// Biopython and BCBio read a leading UTF-8 byte order mark as content while the
// browser's text decoder drops it, so an upload that starts with one declines.
// A Session resource is a file Python already read for its Result.
/** @param {FileLike | null | undefined} file @returns {Promise<string | null>} */
const readFastPathText = async (file) => {
  if (getSessionResourceSource(file)) return readFileText(file);
  if (typeof file?.text !== 'function') return null;
  if (typeof file.arrayBuffer !== 'function') return file.text();
  const bytes = await readFileBytes(file);
  return startsWithUtf8ByteOrderMark(bytes) ? null : bytesToText(bytes);
};

// The fast path returns exactly the records the loader reads, or declines.
/**
 * @param {FileLike[]} files
 * @param {(...texts: string[]) => SequenceRecord[]} parse
 * @returns {Promise<SequenceRecord[] | null>}
 */
const fastPathRecords = async (files, parse) => {
  try {
    const texts = await Promise.all(files.map(readFastPathText));
    return texts.every((text) => typeof text === 'string') ? parse(...texts) : null;
  } catch {
    return null;
  }
};

// The one owner of what is read from a source File (OV-1587): the record lists
// of the Linear cards, the Circular refresh, a draft Save and the comparison
// frame, and the ring comparison sequence with its digest. Each read is kept
// once per File object (or per GFF3 and FASTA pair) and key: an uploaded File and
// a Session resource view never change their bytes. A failed read is not kept;
// each reader keeps an error it treats as final in its own state.
/**
 * @typedef {{ reads: Map<string, Promise<unknown>>, paired: WeakMap<FileLike, FileReads> }} FileReads
 */
/** @type {WeakMap<FileLike, FileReads>} */
const fileReads = new WeakMap();

/**
 * @param {WeakMap<FileLike, FileReads>} byFile
 * @param {FileLike} file
 * @returns {FileReads}
 */
const fileReadsOf = (byFile, file) => {
  let held = byFile.get(file);
  if (!held) byFile.set(file, held = { reads: new Map(), paired: new WeakMap() });
  return held;
};

/**
 * @template T
 * @param {readonly [FileLike] | readonly [FileLike, FileLike]} files The File, or a GFF3 File and its FASTA File.
 * @param {string} key What is read from them.
 * @param {() => Promise<T>} read
 * @returns {Promise<T>}
 */
const readOncePerFile = (files, key, read) => {
  const [file, pairedFile] = files;
  const owner = fileReadsOf(fileReads, file);
  const { reads } = pairedFile ? fileReadsOf(owner.paired, pairedFile) : owner;
  const held = /** @type {Promise<T> | undefined} */ (reads.get(key));
  if (held) return held;
  const started = read();
  reads.set(key, started);
  started.catch(() => { if (reads.get(key) === started) reads.delete(key); });
  return started;
};

// Each caller of a record list gets its own copy of the records.
/** @param {Promise<SequenceRecord[]>} read */
const copiedRecords = (read) => read.then((records) => records.map((record) => ({ ...record })));

/**
 * @param {{
 *   file: FileLike | null | undefined, format: string,
 *   runHelperOperation?: HelperOperationRunner<SequenceRecordsPayload>
 * }} options
 */
export const discoverSequenceRecords = async ({
  file,
  format,
  runHelperOperation = runDiagramHelperOperation
}) => {
  if (!file) throw new Error('A sequence file is required.');
  return copiedRecords(readOncePerFile([file], `records:${format}`, async () => {
    const records = await fastPathRecords([file], (text) => parseSequenceRecordText(text, format));
    if (records) return records;
    const response = await runHelperOperation(
      DIAGRAM_HELPER_OPERATIONS.LIST_SEQUENCE_RECORDS,
      {
        format,
        files: [{ role: 'source', bytes: await cloneFileBytesForTransfer(file) }]
      }
    );
    return normalizeSequenceRecords(response.result);
  }));
};

/**
 * A ring comparison file as the Python ring reader reads it (D12): the LOSAT
 * query FASTA, and the label the file names itself (GenBank / DDBJ: first
 * record's DEFINITION, then organism), or null (FASTA).
 * @typedef {{ fasta: string, recordLabel: string | null }} ComparisonSequence
 */
// Only the Python ring reader reads the file; there is no fast path. Add Seq and
// every Generate share one read per file (OV-1545), frozen, as every reader gets it.
/**
 * @param {{
 *   file: FileLike,
 *   runHelperOperation?: HelperOperationRunner<{ error?: unknown, fasta?: string, recordLabel?: string | null }>
 * }} options
 * @returns {Promise<ComparisonSequence>}
 */
export const readComparisonSequence = ({
  file,
  runHelperOperation = runDiagramHelperOperation
}) => readOncePerFile([file], 'comparison-sequence', async () => {
  const response = await runHelperOperation(
    DIAGRAM_HELPER_OPERATIONS.READ_COMPARISON_SEQUENCE,
    { files: [{ role: 'source', bytes: await cloneFileBytesForTransfer(file) }] }
  );
  if (response.result?.error) throw response.result.error;
  return Object.freeze({
    fasta: String(response.result?.fasta || ''), recordLabel: response.result?.recordLabel ?? null
  });
});

// The SHA-256 of a ring file's FASTA, the LOSAT cache key Generate builds from
// it (OV-1565): digested once per File, beside its FASTA. Only Generate asks.
/**
 * @param {Parameters<typeof readComparisonSequence>[0]} options
 * @returns {Promise<string>}
 */
export const readComparisonSequenceDigest = (options) => readOncePerFile(
  [options.file], 'comparison-sequence-sha256',
  async () => sha256Hex(textToBytes((await readComparisonSequence(options)).fasta))
);

/**
 * @param {{
 *   gffFile: FileLike | null | undefined, fastaFile: FileLike | null | undefined,
 *   runHelperOperation?: HelperOperationRunner<SequenceRecordsPayload>
 * }} options
 */
export const discoverGffFastaRecords = async ({
  gffFile,
  fastaFile,
  runHelperOperation = runDiagramHelperOperation
}) => {
  if (!gffFile || !fastaFile) throw new Error('GFF3 and FASTA files are required.');
  return copiedRecords(readOncePerFile([gffFile, fastaFile], 'records:gff-fasta', async () => {
    const records = await fastPathRecords([gffFile, fastaFile], parseGffFastaRecordText);
    if (records) return records;
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
  }));
};
