// @ts-check
// The sequence actions of the feature and match popups on a shown Result (D5): the feature popup's Sequence tab, the
// match popup's span Copy and FASTA, and the FASTA of a match popup's group members. The catalog stores no sequence,
// so each action reads its bases from Python through the catalog-sequence door when it is pressed (the tab when it is
// shown), shows "Reading the sequence…" while Python reads, and keeps the failure for the popup to show with Details,
// Copy diagnostics and Retry. A read the door discards (the shown Result changed while it read) writes nothing: not
// the popup, not the clipboard, not a file.
/** @import { CatalogSequences } from '../services/catalog-sequences.js' */
/** @import { FeaturePopup } from '../services/feature-popup.js' */
/** @import { MatchSequenceBundle, MatchSpanEntry, MemberRow } from '../services/pairwise-match-popup.js' */
/** @import { FeatureSequencesOf, MatchSpanInput } from '../services/diagram-worker-protocol.js' */
import { featureSequenceKey, featureSequenceReads, readRowsFasta } from '../services/catalog-sequences.js';
import { buildFeatureSequenceFastas } from '../services/feature-sequence-fasta.js';
import { matchSpanFasta, matchSpanUnavailableReason } from '../services/pairwise-match-popup.js';
import { normalizeCaughtError } from '../utils/error-normalization.js';

/**
 * @typedef {object} SequenceReadPorts
 * @property {CatalogSequences['readFeatures']} readFeatures
 * @property {CatalogSequences['readSpans']} readSpans
 * @property {(pending: Promise<string>) => Promise<string>} copyPending Starts the clipboard write with the click.
 * @property {(filename: string, text: string, type?: string) => void} downloadText
 */

/** The read was discarded: nothing is written. */
const DISCARDED = Symbol('discarded');
/** The read holds no sequence: the popup says so. */
const EMPTY = Symbol('empty');

/** @param {string | null} text */
const givenText = (text) => {
  if (text === null) throw DISCARDED;
  if (!text) throw EMPTY;
  return text;
};

/** @param {SequenceReadPorts} ports */
export const createSequenceReads = ({ readFeatures, readSpans, copyPending, downloadText }) => {
  /**
   * The Sequence tab of a feature popup: read once when the tab is shown, again on Retry. A read the door discards
   * (another Result is shown) is read again for the Result now shown, which reads the feature or says why it cannot.
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
      while (read && sequences === null) sequences = await readFeatures([read]);
      const replied = read ? sequences?.get(featureSequenceKey(read)) : undefined;
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
   * reads, then nothing (done or discarded), "No sequence available", or the failure with Retry.
   * @param {MatchSequenceBundle['read']} state
   * @param {() => Promise<string | null>} readText
   * @param {string | null} filename The file to download, or null to copy.
   * @param {'readMatchSequences' | 'readFeatureSequences'} operation
   * @returns {Promise<boolean>} Whether the text was copied or downloaded.
   */
  const runAction = async (state, readText, filename, operation) => {
    Object.assign(state, { status: 'reading', error: null, retry: null });
    try {
      const pending = readText().then(givenText);
      if (filename === null) await copyPending(pending);
      else downloadText(filename, await pending, 'text/plain;charset=utf-8');
      Object.assign(state, { status: '' });
      return true;
    } catch (error) {
      if (error === DISCARDED) Object.assign(state, { status: '' });
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
   * The FASTA text of match spans, read from Python. A span Python cannot read shows its reason in its row and is
   * left out of the text. Null when the door discarded the read.
   * @param {readonly MatchSpanEntry[]} entries
   */
  const spansFasta = async (entries) => {
    const asked = entries.filter((entry) => entry.read);
    if (!asked.length) return '';
    const replies = await readSpans(asked.map((entry) => /** @type {MatchSpanInput} */ (entry.read)));
    if (!replies) return null;
    return asked.map((entry, index) => {
      const reply = replies[index];
      if (reply.status === 'ok') return matchSpanFasta(entry, reply);
      Object.assign(entry, { available: false, unavailableReason: matchSpanUnavailableReason(entry, reply) });
      return '';
    }).join('');
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
   * Copy or Download of the nucleotide or amino acid FASTA of one group member or of every member.
   * @param {MatchSequenceBundle['read']} state The popup's read state.
   * @param {readonly MemberRow[]} rows
   * @param {'nt' | 'aa'} kind
   * @param {'copy' | 'download'} action
   * @param {string} filename
   */
  const memberAction = (state, rows, kind, action, filename) => runAction(
    state, () => readRowsFasta(rows.map((row) => ({ read: row.read, feature: row.fastaFeature })), kind, readFeatures),
    action === 'copy' ? null : filename, 'readFeatureSequences'
  );

  return { showFeatureSequence, spanAction, memberAction };
};
