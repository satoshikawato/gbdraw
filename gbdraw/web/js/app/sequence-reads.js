// @ts-check
// The sequence actions of the feature and match popups on a shown Result (D5): the feature popup's Sequence tab, the
// match popup's span Copy and FASTA, and the FASTA of a match popup's group members. The catalog stores no sequence,
// so each action reads its bases from Python through the catalog-sequence door when it is pressed (the tab when it is
// shown), shows "Reading the sequence…" while Python reads, and keeps the failure for the popup to show with Details,
// Copy diagnostics and Retry. A read the door discards (the shown Result changed while it read) writes nothing: not
// the popup, not the clipboard, not a file.
/** @import { CatalogSequences, SequenceReadState } from '../services/catalog-sequences.js' */
/** @import { FeaturePopup } from '../services/feature-popup.js' */
/** @import { MatchSequenceBundle, MatchSpanEntry, MemberRow } from '../services/pairwise-match-popup.js' */
/** @import { FeatureSequencesOf, MatchSpanInput } from '../services/diagram-worker-protocol.js' */
import { featureSequenceKey, featureSequenceReads, readRowsFasta } from '../services/catalog-sequences.js';
import { buildFeatureSequenceFastas } from '../services/feature-sequence-fasta.js';
import { matchSpanFasta, matchSpanUnavailableReason } from '../services/pairwise-match-popup.js';
import { diagnosticError, normalizeCaughtError } from '../utils/error-normalization.js';

/**
 * @typedef {object} SequenceReadPorts
 * @property {(popup: FeaturePopup) => boolean} isShown Whether the feature popup is the one the page shows.
 * @property {CatalogSequences['readFeatures']} readFeatures
 * @property {CatalogSequences['readSpans']} readSpans
 * @property {(pending: Promise<string>) => Promise<string>} copyPending Starts the clipboard write with the click.
 * @property {(filename: string, text: string, type?: string) => void} downloadText
 */

/**
 * Copy or Download of the nucleotide or amino acid FASTA of group members, with the read state of the popup that
 * asked (`memberAction` below).
 * @typedef {(state: SequenceReadState, rows: readonly Pick<MemberRow, 'read' | 'fastaFeature'>[], kind: 'nt' | 'aa',
 *   action: 'copy' | 'download', filename: string) => Promise<boolean>} MemberAction
 */

/** The read was discarded: nothing is written. */
const DISCARDED = Symbol('discarded');
/** The read holds no sequence: the popup says so. */
const EMPTY = Symbol('empty');
/** A span asked cannot be read: its row shows why, and nothing is written. */
const REFUSED = Symbol('refused');

/**
 * Whether every span of a Copy or FASTA can be read: the popup offers an action for spans only then, so a file or
 * the clipboard holds all the spans it names or none (both spans, as before D5).
 * @param {readonly MatchSpanEntry[]} entries
 */
export const spansReadable = (entries) => entries.length > 0 && entries.every((entry) => entry.available && entry.read);

/**
 * The action each read state shows: a later action on a state supersedes the earlier one, whose end then writes
 * nothing to it. @type {WeakMap<SequenceReadState, object>}
 */
const currentActions = new WeakMap();

/** @param {string | null} text */
const givenText = (text) => {
  if (text === null) throw DISCARDED;
  if (!text) throw EMPTY;
  return text;
};

/** @param {SequenceReadPorts} ports */
export const createSequenceReads = ({ isShown, readFeatures, readSpans, copyPending, downloadText }) => {
  /**
   * The Sequence tab of a feature popup: read once when the tab is shown, again on Retry. A read the door discards
   * (another Result is shown) is read again for the Result now shown, which reads the feature or says why it cannot,
   * while the popup is shown; a popup closed or showing another feature stops reading, and its tab reads when shown.
   * @param {FeaturePopup} popup
   */
  const showFeatureSequence = async (popup) => {
    const { sequence } = popup;
    if (sequence.status === 'reading' || sequence.status === 'ready') return;
    const read = featureSequenceReads([popup.feat])[0];
    Object.assign(sequence, { status: 'reading', error: null });
    try {
      /** @type {ReadonlyMap<string, FeatureSequencesOf> | null} */
      let sequences = null;
      // A feature without a catalog identity is one an older Session drew (from its SVG metadata): Python finds
      // features by identity, so its sequences are read after the next Generate.
      if (!read) throw diagnosticError('SEQUENCE_REQUIRES_GENERATE', {}, { operation: 'readFeatureSequences', stage: 'helper' });
      do sequences = await readFeatures([read]); while (sequences === null && isShown(popup));
      if (sequences === null) {
        Object.assign(sequence, { status: 'idle' });
        return;
      }
      const replied = sequences.get(featureSequenceKey(read));
      const { nucleotideFasta, aminoAcidFasta } = buildFeatureSequenceFastas(popup.feat, {
        nucleotideSequence: replied?.nucleotideSequence ?? '', aminoAcidSequence: replied?.aminoAcidSequence ?? ''
      });
      Object.assign(sequence, { status: 'ready', warnings: replied?.sequenceWarnings ?? [], nucleotideFasta, aminoAcidFasta });
    } catch (error) {
      Object.assign(sequence, { status: 'failed', error: normalizeCaughtError(error, { operation: 'readFeatureSequences', stage: 'helper' }) });
    }
  };

  /**
   * Copies or downloads the text one read makes, with the popup's read state: "Reading the sequence…" while Python
   * reads, then nothing (done or discarded), "No sequence available", the failure with Retry, or "Copy failed" with
   * Retry when the read gave the text and the browser refused the clipboard write.
   * @param {SequenceReadState} state
   * @param {() => Promise<string | null>} readText
   * @param {string | null} filename The file to download, or null to copy.
   * @param {'readMatchSequences' | 'readFeatureSequences'} operation
   * @returns {Promise<boolean>} Whether the text was copied or downloaded.
   */
  const runAction = async (state, readText, filename, operation) => {
    const action = {};
    currentActions.set(state, action);
    Object.assign(state, { status: 'reading', error: null, retry: null });
    try {
      let read = false;
      const pending = readText().then(givenText).then((text) => {
        read = true;
        return text;
      });
      if (filename === null) {
        try {
          await copyPending(pending);
        } catch (error) {
          if (!read) throw error;
          if (currentActions.get(state) === action) {
            Object.assign(state, { status: 'copy-failed', retry: () => runAction(state, readText, filename, operation) });
          }
          return false;
        }
      }
      else downloadText(filename, await pending, 'text/plain;charset=utf-8');
      if (currentActions.get(state) === action) Object.assign(state, { status: '' });
      return true;
    } catch (error) {
      if (currentActions.get(state) !== action) return false;
      if (error === DISCARDED || error === REFUSED) Object.assign(state, { status: '' });
      else if (error === EMPTY) Object.assign(state, { status: 'empty' });
      else {
        Object.assign(state, {
          status: 'failed',
          error: normalizeCaughtError(error, { operation, stage: 'helper' }),
          retry: () => runAction(state, readText, filename, operation)
        });
      }
      return false;
    }
  };

  /**
   * The FASTA text of match spans, read from Python: every span, or none. A span that cannot be read (known when
   * the popup opened, or from Python's reply) shows its reason in its row, and the text is refused (`REFUSED`), so no
   * file or clipboard holds part of what the action names. Null when the door discarded the read.
   * @param {readonly MatchSpanEntry[]} entries
   */
  const spansFasta = async (entries) => {
    if (!spansReadable(entries)) throw REFUSED;
    const replies = await readSpans(entries.map((entry) => /** @type {MatchSpanInput} */ (entry.read)));
    if (!replies) return null;
    entries.forEach((entry, index) => {
      const reply = replies[index];
      if (reply.status !== 'ok') {
        Object.assign(entry, { available: false, unavailableReason: matchSpanUnavailableReason(entry, reply) });
      }
    });
    if (!spansReadable(entries)) throw REFUSED;
    return entries.map((entry, index) => matchSpanFasta(entry, replies[index])).join('');
  };

  /**
   * Copy or FASTA of one span or both spans of a match popup.
   * @param {MatchSequenceBundle} bundle
   * @param {readonly MatchSpanEntry[]} entries
   * @param {'copy' | 'download'} action
   */
  const spanAction = (bundle, entries, action) => runAction(
    bundle.read, () => spansFasta(entries),
    action === 'copy' ? null : (entries.length > 1 ? bundle.combinedFilename : entries[0].filename), 'readMatchSequences'
  );

  /**
   * Copy or Download of the nucleotide or amino acid FASTA of one group member or of every member: the one owner of
   * a member read for the Similarity group popup of a match, the Similarity group in the feature popup, and the
   * Similarity groups drawer.
   * @type {MemberAction}
   */
  const memberAction = (state, rows, kind, action, filename) => runAction(
    state, () => readRowsFasta(rows.map((row) => ({ read: row.read, feature: row.fastaFeature })), kind, readFeatures),
    action === 'copy' ? null : filename, 'readFeatureSequences'
  );

  return { showFeatureSequence, spanAction, memberAction, spansReadable };
};

