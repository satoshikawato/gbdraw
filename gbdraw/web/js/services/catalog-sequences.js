// @ts-check
// The one door to the sequences the feature catalog does not store (D5): a feature's nucleotide and amino acid
// sequences and warnings, the bases of match spans, and the catalog item an exported Interactive SVG embeds. Python
// reads them (gbdraw/web_support/catalog_sequences.py) from the request the shown Result was made from, as Copy FASTA
// of Select region reads its spans (`record-span-sequences.js`); the Web derives no sequence.
//
// The page keeps each feature's and each span's reply while the same Result is shown (the same committed request),
// so a second read waits for nothing. The replies live only here: never in the catalog, a History step, or a Session.
// A reply that arrives after the shown Result changed (Generate, Undo, Redo, a drawing switch, a Load) is discarded:
// the read resolves to `null` and its caller shows nothing of it.
/** @import { FeatureSequenceInput, FeatureSequencesOf, FeatureSequencesPayload, MatchSpanInput, MatchSpanSequence, MatchSequencesPayload, EmbedSequencesPayload } from './diagram-worker-protocol.js' */
/** @import { CommittedCanonicalSession } from './session-request.js' */
/** @import { FeatureLike } from './feature-utils.js' */
/** @import { FeatureCatalogItem } from './feature-catalog.js' */
/** @import { MatchSequenceRecords } from './pairwise-match-popup.js' */
import { diagnosticError } from '../utils/error-normalization.js';
import { buildFeatureSequenceFastas } from './feature-sequence-fasta.js';

/** @typedef {'readFeatureSequences' | 'readMatchSequences' | 'readEmbedSequences'} CatalogSequenceOperation */

/** @param {CatalogSequenceOperation} operation */
const failure = (operation) => /** @type {const} */ ({ stage: 'helper', operation });

/** @param {unknown} value @returns {value is Record<string, unknown>} */
const isRecord = (value) => Boolean(value) && typeof value === 'object' && !Array.isArray(value);

/** @param {unknown} value @returns {value is string[]} */
const isStringList = (value) => Array.isArray(value) && value.every((item) => typeof item === 'string');

/**
 * The `readFeatureSequences` reply, checked against the features asked: one per feature, in order, with the identity
 * asked and strings only.
 * @param {unknown} result
 * @param {readonly FeatureSequenceInput[]} features
 * @returns {readonly FeatureSequencesOf[]}
 */
export const readFeatureSequencesReply = (result, features) => {
  const read = isRecord(result) ? result.features : null;
  if (!Array.isArray(read) || read.length !== features.length || !read.every((item, index) => isRecord(item)
    && item.recordKey === features[index].recordKey && item.biologicalFeatureId === features[index].biologicalFeatureId
    && typeof item.nucleotideSequence === 'string' && typeof item.aminoAcidSequence === 'string'
    && isStringList(item.sequenceWarnings))) {
    throw diagnosticError('RESULT_INVALID', {}, failure('readFeatureSequences'));
  }
  return /** @type {readonly FeatureSequencesOf[]} */ (read);
};

const SPAN_REASONS = new Set([null, 'NEVER_SUPPLIED', 'RECORD_UNKNOWN', 'RECORD_AMBIGUOUS', 'OUT_OF_RANGE']);

/** @param {unknown} item @returns {item is MatchSpanSequence} */
const isSpanSequence = (item) => isRecord(item)
  && (item.status === 'ok' ? item.reason === null : item.status === 'unavailable' && item.reason !== null)
  && SPAN_REASONS.has(/** @type {null} */ (item.reason)) && (item.orientation === '+' || item.orientation === '-')
  && typeof item.recordId === 'string' && typeof item.sequence === 'string' && Number.isSafeInteger(item.length)
  && item.length === /** @type {string} */ (item.sequence).length;

/**
 * The `readMatchSequences` reply, checked against the spans asked: one per span, in order.
 * @param {unknown} result
 * @param {readonly MatchSpanInput[]} spans
 * @returns {readonly MatchSpanSequence[]}
 */
export const readMatchSequencesReply = (result, spans) => {
  const read = isRecord(result) ? result.spans : null;
  if (!Array.isArray(read) || read.length !== spans.length || !read.every(isSpanSequence)) {
    throw diagnosticError('RESULT_INVALID', {}, failure('readMatchSequences'));
  }
  return read;
};

/**
 * The `readEmbedSequences` reply: the catalog item asked, with the same record keys, in the export encoding.
 * @param {unknown} result
 * @param {Readonly<Record<string, unknown>>} item
 * @returns {Readonly<Record<string, unknown>>}
 */
export const readEmbedSequencesReply = (result, item) => {
  if (!isRecord(result) || !Array.isArray(result.features)
    || JSON.stringify(result.recordKeys) !== JSON.stringify(item.recordKeys)) {
    throw diagnosticError('RESULT_INVALID', {}, failure('readEmbedSequences'));
  }
  return result;
};

/** The feature span a read of many features sends at most, in bases: its reply stays far below the reply wall. */
export const FEATURE_READ_BASES = 8_000_000;
/** The characters the page keeps of replies for the shown Result; past it the oldest are dropped. */
export const KEPT_SEQUENCE_CHARACTERS = 64 * 1024 * 1024;

/** @param {FeatureSequenceInput} feature */
const featureKey = ({ recordKey, biologicalFeatureId }) => JSON.stringify([recordKey, biologicalFeatureId]);
/** @param {MatchSpanInput} span */
const spanKey = (span) => JSON.stringify(span.kind === 'record'
  ? [span.kind, span.recordKey, span.start, span.end]
  : [span.kind, span.sourceIndex, span.recordId, span.start, span.end]);

/**
 * The features of one read in chunks of at most `limit` bases each (a feature longer than the limit is a chunk of
 * its own), in the order asked.
 * @template {FeatureSequenceInput} F
 * @param {readonly F[]} features
 * @param {(feature: F) => number} basesOf
 * @param {number} limit
 * @returns {F[][]}
 */
export const featureReadChunks = (features, basesOf, limit) => {
  /** @type {F[][]} */
  const chunks = [];
  let bases = limit;
  features.forEach((feature) => {
    const size = Math.max(1, basesOf(feature));
    if (bases + size > limit) {
      chunks.push([]);
      bases = 0;
    }
    chunks[chunks.length - 1].push(feature);
    bases += size;
  });
  return chunks;
};

/**
 * The Worker helpers the door runs; each sends one payload and returns the helper's `result`.
 * @typedef {object} CatalogSequenceHelpers
 * @property {(payload: FeatureSequencesPayload) => Promise<{ result?: unknown }>} readFeatureSequences
 * @property {(payload: MatchSequencesPayload) => Promise<{ result?: unknown }>} readMatchSequences
 * @property {(payload: EmbedSequencesPayload) => Promise<{ result?: unknown }>} readEmbedSequences
 */

/**
 * A feature to read with the number of bases it spans (its catalog row's `start..end`), which only chunks the read.
 * @typedef {FeatureSequenceInput & { bases: number }} FeatureSequenceRead
 */

/**
 * The door the readers use (the Sequence tab, member Copy and Download, the sequence fields of Feature search, the
 * match popup's Copy and FASTA, Interactive SVG export).
 * @typedef {object} CatalogSequences
 * @property {(features: readonly FeatureSequenceRead[]) => Promise<ReadonlyMap<string, FeatureSequencesOf> | null>} readFeatures
 *   The features' sequences by `featureSequenceKey`, or null when the shown Result changed before they arrived.
 * @property {(spans: readonly MatchSpanInput[]) => Promise<readonly MatchSpanSequence[] | null>} readSpans
 *   The spans' bases in the order asked, or null as above.
 * @property {(item: Readonly<Record<string, unknown>>) => Promise<Readonly<Record<string, unknown>> | null>} readEmbedItem
 *   The shown Result's catalog item in the export encoding (not kept), or null as above.
 */

/** The key `readFeatures` returns a feature's sequences under. */
export const featureSequenceKey = featureKey;

/**
 * A catalog row's sequence identity, which its sequences are read and found by: the record key and the biological
 * feature ID of its catalog row, or null for a row without them.
 * @param {FeatureLike | null | undefined} feature
 * @returns {FeatureSequenceInput | null}
 */
export const featureSequenceIdentity = (feature) => {
  const recordKey = feature?.recordKey;
  const biologicalFeatureId = feature?.biologicalFeatureId;
  return typeof recordKey === 'string' && recordKey && typeof biologicalFeatureId === 'string' && biologicalFeatureId
    ? { recordKey, biologicalFeatureId } : null;
};

/**
 * The reads for catalog rows (searched features, Similarity group members, popup members): each identity once, with the bases of its catalog
 * row's `start..end`, which only chunk the read.
 * @param {readonly FeatureLike[]} features
 * @returns {FeatureSequenceRead[]}
 */
export const featureSequenceReads = (features) => {
  /** @type {Map<string, FeatureSequenceRead>} */
  const reads = new Map();
  features.forEach((feature) => {
    const identity = featureSequenceIdentity(feature);
    if (!identity || reads.has(featureSequenceKey(identity))) return;
    const start = Number(feature.start);
    const end = Number(feature.end);
    const bases = Number.isFinite(start) && Number.isFinite(end) ? Math.max(1, Math.abs(end - start)) : 1;
    reads.set(featureSequenceKey(identity), { ...identity, bases });
  });
  return [...reads.values()];
};


/**
 * The FASTA text of catalog rows for one sequence kind, read from Python in one read: each row's entry with its
 * header fields (`feature`) and the bases replied for its identity (`read`); a row without an identity, or whose reply
 * holds no sequence of the kind, is left out. '' when none has one; null when the door discarded the read.
 * @param {readonly { read: FeatureSequenceRead | null, feature: Record<string, unknown> }[]} rows
 * @param {'nt' | 'aa'} kind
 * @param {CatalogSequences['readFeatures']} readFeatures
 * @returns {Promise<string | null>}
 */
export const readRowsFasta = async (rows, kind, readFeatures) => {
  const asked = rows.flatMap(({ read }) => (read ? [read] : []));
  if (!asked.length) return '';
  const sequences = await readFeatures(asked);
  if (sequences === null) return null;
  return rows.map(({ read, feature }) => {
    const replied = read ? sequences.get(featureKey(read)) : undefined;
    if (!replied) return '';
    const fastas = buildFeatureSequenceFastas(feature, replied);
    const text = kind === 'aa' ? fastas.aminoAcidFasta : fastas.nucleotideFasta;
    return text ? `${text}\n` : '';
  }).join('');
};

/**
 * The records the shown Result's spans are read from, as a match popup names them to Python: the drawn units in
 * drawing order (a span's record index is a place in them), the units that draw a record ID (a Circular comparison
 * band names its reference by ID), and whether the committed request holds a comparison source's FASTA (a
 * never-supplied FASTA is known without reading). Null without a committed request (a Result-only Session 27-39),
 * whose spans are read after the next Generate.
 * @param {FeatureCatalogItem | null | undefined} item The shown Result's stored catalog item.
 * @param {CommittedCanonicalSession | null} committed
 * @returns {MatchSequenceRecords | null}
 */
export const matchSequenceRecordsOf = (item, committed) => {
  if (!committed || !item || !isStringList(item.recordKeys)) return null;
  const unitKeys = item.recordKeys;
  const windowOf = isRecord(item.windowOf) ? item.windowOf : {};
  /** @type {Map<string, Set<string>>} */
  const recordKeysById = new Map();
  (Array.isArray(item.biologicalFeatures) ? item.biologicalFeatures : []).forEach((row) => {
    const recordId = String(row?.record_id ?? row?.recordId ?? '').trim();
    if (!recordId || typeof row.recordKey !== 'string') return;
    if (!recordKeysById.has(recordId)) recordKeysById.set(recordId, new Set());
    /** @type {Set<string>} */ (recordKeysById.get(recordId)).add(row.recordKey);
  });
  const fastaFiles = committed.renderRequest.diagramOptions.conservationFastaFiles;
  return {
    unitKeys,
    unitKeysOfRecordId: (recordId) => {
      const keys = recordKeysById.get(recordId);
      return keys ? unitKeys.filter((unit) => keys.has(windowOf[unit] ?? unit)) : [];
    },
    comparisonSupplied: (sourceIndex) => Array.isArray(fastaFiles) && fastaFiles[sourceIndex] != null
  };
};

/**
 * The three reads as diagram Worker helper operations (the operation names are the Worker protocol's).
 * @param {(operation: string, payload: FeatureSequencesPayload | MatchSequencesPayload | EmbedSequencesPayload) => Promise<{ result?: unknown }>} run
 * @param {{ READ_FEATURE_SEQUENCES: string, READ_MATCH_SEQUENCES: string, READ_EMBED_SEQUENCES: string }} operations
 * @returns {CatalogSequenceHelpers}
 */
export const catalogSequenceHelpers = (run, operations) => ({
  readFeatureSequences: (payload) => run(operations.READ_FEATURE_SEQUENCES, payload),
  readMatchSequences: (payload) => run(operations.READ_MATCH_SEQUENCES, payload),
  readEmbedSequences: (payload) => run(operations.READ_EMBED_SEQUENCES, payload)
});

/**
 * @param {{ getCommittedSession: () => CommittedCanonicalSession | null, helpers: CatalogSequenceHelpers,
 *   keptCharacters?: number, featureReadBases?: number }} ports
 * @returns {CatalogSequences}
 */
export const createCatalogSequences = ({
  getCommittedSession, helpers, keptCharacters = KEPT_SEQUENCE_CHARACTERS, featureReadBases = FEATURE_READ_BASES
}) => {
  /**
   * The replies kept for one committed request; a read in flight is kept as its promise, so a second asker waits
   * for the same reply, and a failed read is forgotten, so Retry reads again.
   * @typedef {{ committed: CommittedCanonicalSession | null, features: Map<string, Promise<FeatureSequencesOf>>,
   *   spans: Map<string, Promise<MatchSpanSequence>>, sizes: Map<string, number>, characters: number }} Kept
   */
  /** @returns {Kept} */
  const fresh = (/** @type {CommittedCanonicalSession | null} */ committed) => ({
    committed, features: new Map(), spans: new Map(), sizes: new Map(), characters: 0
  });
  let kept = fresh(null);
  /** The replies of the shown Result: a Result of another request drops them. */
  const shown = () => {
    const committed = getCommittedSession();
    if (kept.committed !== committed) kept = fresh(committed);
    return kept;
  };
  /**
   * Keeps one settled reply's size and drops the oldest replies past the bound.
   * @param {Kept} owner @param {Map<string, Promise<unknown>>} map @param {string} key @param {number} size
   */
  const account = (owner, map, key, size) => {
    if (owner.committed !== kept.committed || map.get(key) === undefined) return;
    owner.sizes.set(`${map === owner.features ? 'f' : 's'}${key}`, size);
    owner.characters += size;
    for (const [sizeKey, sizeOf] of owner.sizes) {
      if (owner.characters <= keptCharacters) break;
      owner.sizes.delete(sizeKey);
      owner.characters -= sizeOf;
      (sizeKey[0] === 'f' ? owner.features : owner.spans).delete(sizeKey.slice(1));
    }
  };
  /**
   * The committed request a read sends. A shown Result without one (a Result-only Session 27-39) has no input to
   * read until the next Generate commits its request.
   * @param {CatalogSequenceOperation} operation
   */
  const committedOf = (operation) => {
    const owner = shown();
    if (!owner.committed) throw diagnosticError('SEQUENCE_REQUIRES_GENERATE', {}, failure(operation));
    return { owner, request: { canonicalRequest: owner.committed.renderRequest, resources: owner.committed.resources } };
  };

  /** @type {CatalogSequences['readFeatures']} */
  const readFeatures = async (features) => {
    const { owner, request } = committedOf('readFeatureSequences');
    const asked = [...new Map(features.map((feature) => [featureKey(feature), feature])).values()];
    const missing = asked.filter((feature) => !owner.features.has(featureKey(feature)));
    featureReadChunks(missing, (feature) => feature.bases, featureReadBases).forEach((chunk) => {
      const identities = chunk.map(({ recordKey, biologicalFeatureId }) => ({ recordKey, biologicalFeatureId }));
      const reply = helpers.readFeatureSequences({ ...request, features: identities })
        .then((result) => readFeatureSequencesReply(result?.result, identities));
      identities.forEach((identity, index) => {
        const key = featureKey(identity);
        const read = reply.then((sequences) => sequences[index]);
        owner.features.set(key, read);
        read.then((sequences) => account(owner, owner.features, key,
          sequences.nucleotideSequence.length + sequences.aminoAcidSequence.length), () => {
          if (owner.features.get(key) === read) owner.features.delete(key);
        });
      });
    });
    const read = await Promise.all(asked.map((feature) => /** @type {Promise<FeatureSequencesOf>} */ (
      owner.features.get(featureKey(feature)))));
    if (getCommittedSession() !== owner.committed) return null;
    return new Map(read.map((sequences) => [featureKey(sequences), sequences]));
  };

  /** @type {CatalogSequences['readSpans']} */
  const readSpans = async (spans) => {
    const { owner, request } = committedOf('readMatchSequences');
    const missing = [...new Map(spans.filter((span) => !owner.spans.has(spanKey(span)))
      .map((span) => [spanKey(span), span])).values()];
    if (missing.length) {
      const reply = helpers.readMatchSequences({ ...request, spans: missing })
        .then((result) => readMatchSequencesReply(result?.result, missing));
      missing.forEach((span, index) => {
        const key = spanKey(span);
        const read = reply.then((sequences) => sequences[index]);
        owner.spans.set(key, read);
        read.then((sequence) => account(owner, owner.spans, key, sequence.length), () => {
          if (owner.spans.get(key) === read) owner.spans.delete(key);
        });
      });
    }
    const read = await Promise.all(spans.map((span) => /** @type {Promise<MatchSpanSequence>} */ (
      owner.spans.get(spanKey(span)))));
    return getCommittedSession() === owner.committed ? read : null;
  };

  /** @type {CatalogSequences['readEmbedItem']} */
  const readEmbedItem = async (item) => {
    const { owner, request } = committedOf('readEmbedSequences');
    const reply = await helpers.readEmbedSequences({ ...request, item });
    const embedded = readEmbedSequencesReply(reply?.result, item);
    return getCommittedSession() === owner.committed ? embedded : null;
  };

  return { readFeatures, readSpans, readEmbedItem };
};
