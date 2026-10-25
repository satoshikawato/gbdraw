// @ts-check
import {
  DIAGRAM_HELPER_OPERATIONS,
  runDiagramHelperOperation
} from './diagram-generation.js';
import { genbankHeaderIds, scanGenBankHeader } from './genbank-header.js';
import { ORGANISM_NAME_RULE } from '../contracts/schema.generated.js';
import { sha256Hex, textToBytes } from './byte-utils.js';
import { diagnosticError } from '../utils/error-normalization.js';
import {
  bytesToText,
  cloneFileBytesForTransfer,
  getSessionResourceSource,
  readFileBytes
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

/** @param {string} subject @returns {never} */
const declineToPythonReader = (subject) => {
  throw new Error(`${subject} needs the Python reader.`);
};

// The Auto organism name: this fast path runs the steps of split_organism_lines
// (gbdraw/core/record_metadata.py) on the word lists and patterns that
// tools/generate_web_contracts.py writes from it (ORGANISM_NAME_RULE), so the
// Linear Auto definition needs no Worker. tools/twins.json registers this twin;
// tests/web/record-metadata-inference.test.mjs and tests/test_record_metadata.py
// share one vector table so the two cannot drift. Python's string and pattern
// rules equal the browser's only for ASCII text (and the hybrid sign), so other
// text declines to the Worker (R4).
/** @typedef {{ source: string, flags: string }} NamePattern */
/** @param {NamePattern} pattern */
const searchPattern = ({ source, flags }) => new RegExp(source, flags);
/** @param {NamePattern} pattern Python's `match`: anchored at the start. */
const matchPattern = ({ source, flags }) => new RegExp(`^(?:${source})`, flags);
/** @param {NamePattern} pattern Python's `fullmatch`. */
const fullPattern = ({ source, flags }) => new RegExp(`^(?:${source})$`, flags);
const NAME = Object.freeze({
  nonOrganismNames: new Set(ORGANISM_NAME_RULE.nonOrganismNames),
  candidatusPrefix: matchPattern(ORGANISM_NAME_RULE.candidatusPrefix),
  nonOrganismStart: matchPattern(ORGANISM_NAME_RULE.nonOrganismStart),
  sequenceWords: new Set(ORGANISM_NAME_RULE.sequenceWords),
  virusWord: searchPattern(ORGANISM_NAME_RULE.virusWord),
  virusAcronym: fullPattern(ORGANISM_NAME_RULE.virusAcronym),
  ictvGenus: fullPattern(ORGANISM_NAME_RULE.ictvGenus),
  cellularLineages: new Set(ORGANISM_NAME_RULE.cellularLineages),
  virusLineage: ORGANISM_NAME_RULE.virusLineage,
  hybrid: searchPattern(ORGANISM_NAME_RULE.hybrid),
  genus: fullPattern(ORGANISM_NAME_RULE.genus),
  bracketedGenus: fullPattern(ORGANISM_NAME_RULE.bracketedGenus),
  aboveGenus: fullPattern(ORGANISM_NAME_RULE.aboveGenus),
  epithet: fullPattern(ORGANISM_NAME_RULE.epithet),
  latinRankWords: new Set(ORGANISM_NAME_RULE.latinRankWords),
  strainLineWords: new Set(ORGANISM_NAME_RULE.strainLineWords),
  notEpithets: new Set(ORGANISM_NAME_RULE.notEpithets),
  organismNoun: searchPattern(ORGANISM_NAME_RULE.organismNoun),
  statusSuffix: fullPattern(ORGANISM_NAME_RULE.statusSuffix),
  designationPart: searchPattern(ORGANISM_NAME_RULE.designationPart),
  lowercaseWord: fullPattern(ORGANISM_NAME_RULE.lowercaseWord),
  lowercaseStrain: fullPattern(ORGANISM_NAME_RULE.lowercaseStrain)
});
const BROWSER_EQUAL_TEXT = /^[\t\n\v\f\r\x20-\x7e×]*$/;

/** @param {string} text Python's `" ".join(text.split())`. */
const wordsOf = (text) => text.split(/\s+/).filter(Boolean);
/** @param {string} text */
const lettersOf = (text) => text.toLowerCase().replace(/[^0-9a-z]+/g, '');
/** @param {string} word @param {string} chars Python's `str.strip(chars)`. */
const stripChars = (word, chars) => {
  let start = 0;
  let end = word.length;
  while (start < end && chars.includes(word[start])) start += 1;
  while (end > start && chars.includes(word[end - 1])) end -= 1;
  return word.slice(start, end);
};

/**
 * The index of the first word where the name spells the strain, or null.
 * @param {string[]} tokens
 * @param {string} strain
 * @returns {number | null}
 */
const namesStrain = (tokens, strain) => {
  const target = lettersOf(strain);
  if (!target) return null;
  for (let index = 0; index < tokens.length; index += 1) {
    const parts = tokens[index].split(NAME.designationPart).filter(Boolean).map(lettersOf);
    if (parts.includes(target)) return index;
    let joined = '';
    for (const later of tokens.slice(index)) {
      joined += lettersOf(later);
      if (joined === target) return index;
      if (!target.startsWith(joined)) break;
    }
  }
  return null;
};

/**
 * The words after the taxon, then the strain unless the name holds it (OV-1830).
 * @param {string[]} rest
 * @param {string} strain
 * @param {string[]} nameTokens
 */
const strainLineOf = (rest, strain, nameTokens) => {
  const words = rest.join(' ');
  if (!lettersOf(strain) || namesStrain(nameTokens, strain) !== null) return words;
  const missing = strain.split(/\s*;\s*/).filter(Boolean).filter((part) => namesStrain(nameTokens, part) === null);
  return `${words} ${missing.join('; ')}`.trim();
};

/** @param {string} token */
const isEpithet = (token) => {
  const last = token.slice(token.lastIndexOf('-') + 1);
  return NAME.epithet.test(token) && !NAME.notEpithets.has(last) && !NAME.organismNoun.test(last);
};

/** @param {string[]} tokens @param {number} start */
const parentheticalEnd = (tokens, start) => {
  if (start >= tokens.length || !tokens[start].startsWith('(')) return start;
  let depth = 0;
  for (let end = start; end < tokens.length; end += 1) {
    depth += tokens[end].split('(').length - tokens[end].split(')').length;
    if (depth <= 0) return end + 1;
  }
  return start;
};

/** @param {string[]} words @returns {string | null} */
const nameSuffix = (words) => {
  const text = words.join(' ');
  if (NAME.statusSuffix.test(text)) return text;
  const inner = text.startsWith('(') && text.endsWith(')') ? wordsOf(text.slice(1, -1)) : [];
  if (inner.length >= 2 && NAME.genus.test(inner[0])) {
    if (isEpithet(inner[1])) return `(<i>${inner[0]} ${inner[1]}</i>${inner.slice(2).map((word) => ` ${word}`).join('')})`;
    if (inner[1] === 'sp.' || inner[1] === 'cf.') return `(<i>${inner[0]}</i> ${inner.slice(1).join(' ')})`;
  }
  return null;
};

/**
 * The taxon line's parts and the number of words they take, or null to keep the name whole.
 * @param {string[]} tokens
 * @param {number} strainAt
 * @returns {[string[], number] | null}
 */
const taxonWords = (tokens, strainAt) => {
  const bracketed = tokens[0].match(NAME.bracketedGenus);
  const genus = bracketed ? bracketed[1] : tokens[0];
  if (!(bracketed || NAME.genus.test(genus)) || NAME.aboveGenus.test(genus)) return null;
  const parts = [bracketed ? `[<i>${genus}</i>]` : `<i>${genus}</i>`];
  if (tokens.length === 1) return [parts, 1];
  if (tokens[1] === 'sp.') {
    const end = parentheticalEnd(tokens, 2);
    const homonym = end > 2 && tokens[2] === '(in:';
    return [[...parts, 'sp.', ...(homonym ? tokens.slice(2, end) : [])], homonym ? end : 2];
  }
  if (strainAt <= 1 || !isEpithet(tokens[1])) return null;
  parts.push(`<i>${tokens[1]}</i>`);
  let index = 2;
  while (index < strainAt && NAME.latinRankWords.has(tokens[index])) {
    const rank = tokens[index] === 'f.' && index + 1 < tokens.length && tokens[index + 1] === 'sp.' ? 'f. sp.' : tokens[index];
    const epithetAt = index + wordsOf(rank).length;
    if (epithetAt >= strainAt || !isEpithet(tokens[epithetAt])) return null;
    parts.push(`${rank} <i>${tokens[epithetAt]}</i>`);
    index = epithetAt + 1;
  }
  if (
    index === 2
    && index < strainAt
    && isEpithet(tokens[index])
    && (index + 1 === tokens.length || index + 1 === strainAt || !NAME.lowercaseWord.test(tokens[index + 1]))
  ) {
    parts.push(`<i>${tokens[index]}</i>`);
    index += 1;
  }
  for (;;) {
    const end = parentheticalEnd(tokens, index);
    const suffix = end > index ? nameSuffix(tokens.slice(index, end)) : null;
    if (!suffix) break;
    parts.push(suffix);
    index = end;
  }
  const rest = tokens.slice(index);
  if (rest.length > 0 && index < strainAt && /^[a-z]/.test(rest[0]) && !NAME.strainLineWords.has(rest[0])) return null;
  if (rest.some((word) => NAME.organismNoun.test(stripChars(word.toLowerCase(), '()\'",;').split('-').pop() || ''))) {
    return null;
  }
  return [parts, index];
};

/**
 * A record's Auto name as a taxon line and a strain line (text with `<i>`),
 * as split_organism_lines writes them.
 * @param {string} rawOrganism
 * @param {string} rawStrain
 * @param {string | null} lineageTop The first word of the record's lineage, or null.
 * @returns {{ taxon: string, strain: string }}
 */
const organismNameLines = (rawOrganism, rawStrain, lineageTop) => {
  const organism = wordsOf(rawOrganism).join(' ');
  const strain = wordsOf(rawStrain).join(' ');
  const prefix = organism.match(NAME.candidatusPrefix);
  const name = prefix ? organism.slice(prefix[0].length) : organism;
  const tokens = wordsOf(name);
  const whole = wordsOf(organism);
  // Bakta writes "." when it was given no organism.
  if (!/[A-Za-z]/.test(name) || NAME.nonOrganismNames.has(name.toLowerCase())) return { taxon: '', strain };
  const keepWhole = { taxon: organism, strain: strainLineOf([], strain, whole) };
  if (/[&<>]/.test(organism + strain)) return keepWhole;
  if (lineageTop === NAME.virusLineage || (
    !NAME.cellularLineages.has(lineageTop || '')
    && tokens.some((token) => NAME.virusWord.test(token) || NAME.virusAcronym.test(token))
  )) {
    if (!prefix && tokens.length === 2 && NAME.ictvGenus.test(tokens[0]) && NAME.epithet.test(tokens[1])) {
      return { taxon: `<i>${organism}</i>`, strain: keepWhole.strain };
    }
    return keepWhole;
  }
  if (
    NAME.nonOrganismStart.test(name)
    || NAME.hybrid.test(name)
    || tokens.some((token) => NAME.sequenceWords.has(stripChars(token.toLowerCase(), '()\'",;')))
  ) {
    return keepWhole;
  }
  const strainAt = namesStrain(tokens, strain);
  const noBoundary = strainAt === null || strainAt < 2 || NAME.lowercaseStrain.test(strain);
  const taxon = taxonWords(tokens, noBoundary ? tokens.length : strainAt);
  if (!taxon) return keepWhole;
  const [parts, used] = taxon;
  let line = parts.join(' ').replaceAll('</i> <i>', ' ');
  if (prefix) line = `<i>Ca.</i> ${line.replace(/<\/?i>/g, '')}`;
  return { taxon: line, strain: strainLineOf(tokens.slice(used), strain, whole) };
};

/**
 * The record's Auto name on one line, as format_inferred_definition writes it;
 * text the browser cannot read as Python does declines to the Worker.
 * @param {{ organism?: unknown, strain?: unknown, lineageTop?: string | null }} [fields]
 */
export const formatInferredOrganismStrain = ({ organism = '', strain = '', lineageTop = null } = {}) => {
  const rawOrganism = String(organism || '');
  const rawStrain = String(strain || '');
  if (!BROWSER_EQUAL_TEXT.test(rawOrganism) || !BROWSER_EQUAL_TEXT.test(rawStrain)) {
    declineToPythonReader('This organism name');
  }
  const lines = organismNameLines(rawOrganism, rawStrain, lineageTop);
  return [lines.taxon, lines.strain].filter(Boolean).join(' ');
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

  const scanned = header === undefined ? scanGenBankHeader(text) : header;
  // infer_record_source_metadata prefers the source /organism to ORGANISM.
  const organism = extractQualifier('organism', sourceBlock) || scanned?.organism || '';
  const strain = extractQualifier('strain', sourceBlock);
  const isolate = extractQualifier('isolate', sourceBlock);

  // infer_record_source_metadata reads /isolate before /strain.
  const effectiveStrain = isolate || strain || '';
  const inferredDefinition = formatInferredOrganismStrain({
    organism,
    strain: effectiveStrain,
    lineageTop: scanned?.lineageTop ?? null
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

// The text the browser readers (the record lists here, the LOSAT FASTA in
// run-analysis.js) read from a file's bytes, or null to decline to the Python
// reader. Biopython and BCBio read a leading UTF-8 byte order mark as content
// (Biopython then skips the first record) while the browser's text decoder
// drops it, so bytes that start with one decline.
/** @param {Uint8Array} bytes @returns {string | null} */
export const browserReaderText = (bytes) => (
  startsWithUtf8ByteOrderMark(bytes) ? null : bytesToText(bytes)
);

// The text the no-Worker fast path reads, or null to decline to the Worker. A
// Session resource view reads its bytes like a File, with the same rule (OV-1599).
/** @param {FileLike | null | undefined} file @returns {Promise<string | null>} */
const readFastPathText = async (file) => {
  if (getSessionResourceSource(file)) return browserReaderText(await readFileBytes(file));
  if (typeof file?.text !== 'function') return null;
  if (typeof file.arrayBuffer !== 'function') return file.text();
  return browserReaderText(await readFileBytes(file));
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
// frame, the ring comparison sequence with its digest, and what Generate reads
// for LOSAT (the parsed sequences, the Circular reference FASTA, and the
// Worker's FASTA and protein extractions). Each read is kept once per File
// object (or per GFF3 and FASTA pair) and key: an uploaded File and a Session
// resource view never change their bytes. A failed read is not kept; each
// reader keeps an error it treats as final in its own state.
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
 * @param {string} key What is read from them: its kind, then `:` and the settings of the reading.
 * @param {() => Promise<T>} read
 * @param {number} [limit] How many reads of this kind the File keeps, the oldest dropped first.
 * @returns {Promise<T>}
 */
export const readOncePerFile = (files, key, read, limit = Infinity) => {
  const [file, pairedFile] = files;
  const owner = fileReadsOf(fileReads, file);
  const { reads } = pairedFile ? fileReadsOf(owner.paired, pairedFile) : owner;
  const held = /** @type {Promise<T> | undefined} */ (reads.get(key));
  if (held) return held;
  if (Number.isFinite(limit)) {
    const kind = key.slice(0, key.indexOf(':') + 1);
    const sameKind = [...reads.keys()].filter((readKey) => readKey.startsWith(kind));
    if (sameKind.length >= limit) reads.delete(sameKind[0]);
  }
  const started = read();
  reads.set(key, started);
  started.catch(() => { if (reads.get(key) === started) reads.delete(key); });
  return started;
};

/**
 * A record of a sequence file, in file order, as the browser sequence readers (LOSAT, the restored match
 * sequences) parse it. `sequence` is null for a GenBank record without one (a CON record): it keeps its place,
 * so `#N` and the automatic first record name the records Python reads, and only the Python reader reports its
 * sequence. `aliases` are the names a match may give the record: its ID, then the GenBank VERSION, ACCESSION and
 * LOCUS names or the whole FASTA header.
 * @typedef {{ id: string, aliases: readonly string[], sequence: string | null }} ParsedSequenceRecord
 */
/** @param {unknown} text @returns {ParsedSequenceRecord[]} */
const parseFastaRecordsFast = (text) => {
  /** @type {ParsedSequenceRecord[]} */
  const records = [];
  let current = /** @type {{ id: string, header: string, parts: string[] } | null} */ (null);
  /** @param {{ id: string, header: string, parts: string[] }} entry */
  const push = ({ id, header, parts }) => records.push({
    id, aliases: [...new Set([id, header].filter(Boolean))], sequence: parts.join('').toUpperCase()
  });
  String(text || '').split(/\r?\n/).forEach((line) => {
    if (line.startsWith('>')) {
      if (current) push(current);
      const header = line.slice(1).trim();
      current = { id: header.split(/\s+/)[0] || `record_${records.length + 1}`, header, parts: [] };
      return;
    }
    if (current && line.trim()) current.parts.push(line.replace(/\s+/g, ''));
  });
  if (current) push(current);
  return records;
};
/** @param {unknown} text @returns {ParsedSequenceRecord[]} */
const parseGenbankRecordsFast = (text) => {
  /** @type {ParsedSequenceRecord[]} */
  const records = [];
  String(text || '').split(/^\/\/\s*$/m).forEach((chunk) => {
    const header = genbankHeaderIds(chunk);
    if (!header) return;
    const { recordId, version, accession, locus } = header;
    const originMatch = chunk.match(/\nORIGIN\b([\s\S]*)$/i);
    const sequence = originMatch ? originMatch[1].replace(/[^A-Za-z]/g, '').toUpperCase() : '';
    records.push({
      id: recordId,
      aliases: [...new Set([recordId, version, accession, locus].filter(Boolean))],
      sequence: sequence || null
    });
  });
  return records;
};
/**
 * @param {string} text
 * @param {string} fmt `genbank` or `fasta`.
 * @returns {ParsedSequenceRecord[]}
 */
export const parseSequenceRecordsFast = (text, fmt) => (
  fmt === 'genbank' ? parseGenbankRecordsFast(text) : parseFastaRecordsFast(text)
);
/**
 * A file's records, parsed once per File and format, frozen, or null when only the Python reader reads the
 * file (`browserReaderText`); `parsed` says whether this call read the file.
 * @param {FileLike} file
 * @param {string} fmt `genbank` or `fasta`.
 * @returns {Promise<{ records: readonly ParsedSequenceRecord[] | null, parsed: boolean }>}
 */
export const readParsedSequenceRecords = async (file, fmt) => {
  let parsed = false;
  const records = await readOncePerFile([file], `sequence-records:${fmt}`, async () => {
    parsed = true;
    // Decoded from the bytes Generate already read: a Session resource view would keep a text read beside them.
    const text = browserReaderText(await readFileBytes(file));
    return text === null ? null : Object.freeze(parseSequenceRecordsFast(text, fmt).map((record) => Object.freeze({
      ...record, aliases: Object.freeze(record.aliases)
    })));
  });
  return { records, parsed };
};

/**
 * One record of the Worker's `readAllSequenceRecords` reply, as `read_all_sequence_records` writes it: the ID, the
 * names a match may give the record (the ID first), and the sequence (null for a record without one).
 * @typedef {{ id: string, aliases: string[], sequence: string | null }} AllSequenceRecordReply
 */
/**
 * The Worker's `readAllSequenceRecords` reply: every record of the file, or the error the reader raised.
 * @typedef {{ records: AllSequenceRecordReply[] } | { error: unknown }} AllSequenceRecordsPayload
 */
/** @param {unknown} entry @returns {entry is AllSequenceRecordReply} */
const isAllSequenceRecordReply = (entry) => {
  if (typeof entry !== 'object' || entry === null) return false;
  const { id, aliases, sequence } = /** @type {Record<string, unknown>} */ (entry);
  return typeof id === 'string'
    && Array.isArray(aliases) && aliases.every((alias) => typeof alias === 'string')
    && (sequence === null || typeof sequence === 'string');
};
/**
 * The records of a `readAllSequenceRecords` reply, frozen. A reply that is not one fails as a helper protocol error,
 * as an invalid record list does; it is not kept.
 * @param {AllSequenceRecordsPayload | undefined} payload
 * @returns {readonly ParsedSequenceRecord[]}
 */
const allSequenceRecordsFromReply = (payload) => {
  if (payload && 'error' in payload) throw payload.error;
  if (!payload || !Array.isArray(payload.records) || !payload.records.every(isAllSequenceRecordReply)) {
    throw diagnosticError('HELPER_PROTOCOL', {}, {
      stage: 'helper', operation: DIAGRAM_HELPER_OPERATIONS.READ_ALL_SEQUENCE_RECORDS
    });
  }
  return Object.freeze(payload.records.map(({ id, aliases, sequence }) => Object.freeze({
    id, aliases: Object.freeze([...aliases]), sequence
  })));
};
/**
 * Every record of a sequence file, in file order, as the CLI loader reads it. The browser records
 * (`readParsedSequenceRecords`) serve when the browser reader and the record-list rules (`parseSequenceRecordText`)
 * both take the file; otherwise the Worker reads it (OV-1599): a file that starts with a UTF-8 byte order mark,
 * whose first record Python skips, or a FASTA entry without an ID, which Python names '' and the browser record_N.
 * Frozen, read once per File and format. The Circular conservation reference and the restored Linear and reference
 * match sequences read it.
 * @param {FileLike} file
 * @param {string} fmt `genbank` or `fasta`.
 * @param {HelperOperationRunner<AllSequenceRecordsPayload>} [runHelperOperation]
 * @returns {Promise<readonly ParsedSequenceRecord[]>}
 */
export const readAllSequenceRecords = (file, fmt, runHelperOperation = runDiagramHelperOperation) => (
  readOncePerFile([file], `all-sequence-records:${fmt}`, async () => {
    const { records } = await readParsedSequenceRecords(file, fmt);
    if (records && await fastPathRecords([file], (text) => parseSequenceRecordText(text, fmt))) return records;
    const { result } = await runHelperOperation(DIAGRAM_HELPER_OPERATIONS.READ_ALL_SEQUENCE_RECORDS, {
      format: fmt, files: [{ role: 'source', bytes: await cloneFileBytesForTransfer(file) }]
    });
    return allSequenceRecordsFromReply(result);
  })
);

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
 * query FASTA, the label the file names itself (GenBank / DDBJ: first
 * record's DEFINITION, then organism) or null (FASTA), and the ID of each
 * record of the FASTA, in its order ('' for a FASTA entry without an ID).
 * @typedef {{ fasta: string, recordLabel: string | null, recordIds: readonly string[] }} ComparisonSequence
 */
/**
 * The Worker's `readComparisonSequence` reply, as `read_comparison_sequence_json` writes it, or the error it raised.
 * @typedef {{ fasta: string, recordLabel: string | null, recordIds: string[] } | { error: unknown }} ComparisonSequencePayload
 */
/** @param {unknown} payload @returns {payload is { fasta: string, recordLabel: string | null, recordIds: string[] }} */
const isComparisonSequenceReply = (payload) => {
  if (typeof payload !== 'object' || payload === null) return false;
  const { fasta, recordLabel, recordIds } = /** @type {Record<string, unknown>} */ (payload);
  return typeof fasta === 'string' && (recordLabel === null || typeof recordLabel === 'string')
    && Array.isArray(recordIds) && recordIds.every((recordId) => typeof recordId === 'string');
};
// Only the Python ring reader reads the file; there is no fast path. Add Seq and
// every Generate share one read per file (OV-1545), frozen, as every reader gets it.
/**
 * @param {{
 *   file: FileLike,
 *   runHelperOperation?: HelperOperationRunner<ComparisonSequencePayload>
 * }} options
 * @returns {Promise<ComparisonSequence>}
 */
export const readComparisonSequence = ({
  file,
  runHelperOperation = runDiagramHelperOperation
}) => readOncePerFile([file], 'comparison-sequence', async () => {
  const { result } = await runHelperOperation(
    DIAGRAM_HELPER_OPERATIONS.READ_COMPARISON_SEQUENCE,
    { files: [{ role: 'source', bytes: await cloneFileBytesForTransfer(file) }] }
  );
  if (result && 'error' in result) throw result.error;
  if (!isComparisonSequenceReply(result)) {
    throw diagnosticError('HELPER_PROTOCOL', {}, {
      stage: 'helper', operation: DIAGRAM_HELPER_OPERATIONS.READ_COMPARISON_SEQUENCE
    });
  }
  const { fasta, recordLabel, recordIds } = result;
  return Object.freeze({ fasta, recordLabel, recordIds: Object.freeze([...recordIds]) });
});

/**
 * A ring comparison genome's records as the ring search read them, for the restored match sequences: the ID is
 * the record's ID, its only alias, and the sequence is upper case. The browser reads a FASTA file the browser
 * reader and the record-list rules take (no byte order mark, the first entry starts with '>', every entry has an
 * ID) and whose every entry has a sequence, as the comparison reader would; the comparison reader
 * (`readComparisonSequence`) reads every other file (a byte order mark, GenBank / DDBJ, an entry without an ID),
 * with the IDs its reply names. Frozen, once per File.
 * @param {FileLike} file
 * @param {HelperOperationRunner<ComparisonSequencePayload>} [runHelperOperation]
 * @returns {Promise<readonly ParsedSequenceRecord[]>}
 */
export const readComparisonSequenceRecords = (file, runHelperOperation = runDiagramHelperOperation) => (
  readOncePerFile([file], 'comparison-records', async () => {
    /** @param {string} id @param {string | null} sequence */
    const record = (id, sequence) => Object.freeze({ id, aliases: Object.freeze([id]), sequence });
    const { records } = await readParsedSequenceRecords(file, 'fasta');
    if (records && records.every(({ sequence }) => sequence)
      && await fastPathRecords([file], (text) => parseSequenceRecordText(text, 'fasta'))) {
      return Object.freeze(records.map(({ id, sequence }) => record(id, sequence)));
    }
    const { fasta, recordIds } = await readComparisonSequence({ file, runHelperOperation });
    const sequences = parseSequenceRecordsFast(fasta, 'fasta');
    if (sequences.length !== recordIds.length) {
      throw diagnosticError('HELPER_PROTOCOL', {}, {
        stage: 'helper', operation: DIAGRAM_HELPER_OPERATIONS.READ_COMPARISON_SEQUENCE
      });
    }
    return Object.freeze(sequences.map(({ sequence }, index) => record(recordIds[index], sequence)));
  })
);

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

/**
 * The selector a FASTA reader of a GFF3 and FASTA pair reads for a record selector. load_gff_fasta reads only the
 * records the GFF3 names, in FASTA order, while a FASTA reader reads every entry. So the automatic record (an empty
 * selector) and `#N` become the record ID the record list gives them (OV-1681). A record ID, a `#N` the list does
 * not hold, or a card without both files keeps its selector.
 * @param {{ gffFile: FileLike | null | undefined, fastaFile: FileLike | null | undefined, selector: string }} options
 * @returns {Promise<string>}
 */
export const gffFastaReaderSelector = async ({ gffFile, fastaFile, selector }) => {
  if ((selector && !selector.startsWith('#')) || !gffFile || !fastaFile) return selector;
  const records = await discoverGffFastaRecords({ gffFile, fastaFile });
  const record = selector ? records.find((entry) => entry.selector === selector) : records[0];
  return record?.recordId ?? selector;
};
