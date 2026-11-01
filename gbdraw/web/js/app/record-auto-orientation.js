// @ts-check
// "Orient records from comparisons" (Linear): Python decides, with the rule
// `--auto_orient` runs (gbdraw/web_support/auto_orientation.py), from the
// comparisons of the last Generate as the Session keeps them; the record
// display owner writes each changed record's Reverse complement as one History
// step. The next Generate draws it. The rule is not repeated here.
/** @import { AutoOrientationHelperPayload, AutoOrientationReply, RecordOrientationDecision } from '../services/diagram-worker-protocol.js' */
/** @import { CommittedCanonicalSession } from '../services/session-request.js' */
import { diagnosticError } from '../utils/error-normalization.js';
import { msg, t, tm } from '../utils/i18n.js';

const ORIENTATION_FAILURE = /** @type {const} */ ({ stage: 'helper', operation: 'decideAutoOrientation' });

/**
 * Read an `AutoOrientationReply` from the helper's result.
 * @param {unknown} value
 * @param {() => unknown} invalid Builds the caller's classified error for a malformed value.
 * @returns {AutoOrientationReply}
 */
const readAutoOrientationReply = (value, invalid) => {
  const records = value && typeof value === 'object' && !Array.isArray(value)
    ? /** @type {{ records?: unknown }} */ (value).records : null;
  if (!Array.isArray(records)) throw invalid();
  /** @type {RecordOrientationDecision[]} */
  const read = [];
  for (const entry of records) {
    const record = /** @type {{ recordKey?: unknown, reverseComplement?: unknown }} */ (
      entry && typeof entry === 'object' && !Array.isArray(entry) ? entry : {}
    );
    const { recordKey, reverseComplement } = record;
    if (typeof recordKey !== 'string' || !recordKey || typeof reverseComplement !== 'boolean'
      || read.some((known) => known.recordKey === recordKey)) throw invalid();
    read.push(Object.freeze({ recordKey, reverseComplement }));
  }
  return Object.freeze({ records: Object.freeze(read) });
};

/**
 * What the `writeDecidedOrientations` port did: `notLinear` when the shown
 * drawing is not Linear, `busy` (with the reason) when a Session operation
 * blocks edits, `stale` when the decisions do not name the committed records
 * or a record's orientation changed since the last Generate, `unchanged` when
 * no drawn orientation changes, otherwise the record IDs whose Reverse
 * complement it wrote.
 * @typedef {{ status: 'notLinear' } | { status: 'busy', reason: string } | { status: 'stale' } | { status: 'unchanged' }
 *   | { status: 'written', recordIds: string[] }} RecordOrientationWrite
 */

/**
 * @typedef {object} RecordAutoOrientationOptions
 * @property {<T>(getter: () => T) => { readonly value: T }} computed Vue `computed`
 * @property {<T>(value: T) => { value: T }} ref Vue `ref`
 * @property {() => CommittedCanonicalSession | null} getCommittedSession
 * @property {() => number} readCommittedRevision A reactive read that changes when Generate commits a Session.
 * @property {() => { reason: string } | null | undefined} busy A Session operation that blocks edits, with its reason.
 * @property {(payload: AutoOrientationHelperPayload) => Promise<{ result?: unknown }>} decideOrientation
 *   Runs the Worker helper `decideAutoOrientation`; the reply is read here.
 * @property {(decisions: AutoOrientationReply['records'], label: string) => Promise<RecordOrientationWrite>} writeDecidedOrientations
 * @property {(error: unknown) => string} errorMessage The user-facing text of a failure.
 */

/**
 * The committed Linear request has a nucleotide comparison for the rule to read.
 * @param {CommittedCanonicalSession | null} session
 */
const hasNucleotideComparison = (session) => (
  session?.renderRequest.mode === 'linear'
  && session.renderRequest.comparisons.some((comparison) => comparison.kind === 'nucleotideBlast')
);

/**
 * The status line after the writer answered.
 * @param {RecordOrientationWrite} written
 */
const writeStatusText = (written) => {
  switch (written.status) {
    case 'written':
      return t('Set Reverse complement from the comparisons: changed {records}. Generate to draw the new orientation.',
        { records: written.recordIds.join(', ') });
    case 'unchanged': return t('The comparisons give every record its current orientation.');
    case 'busy': return tm(written.reason);
    case 'notLinear': return t('The shown drawing changed while the records were being oriented. Orient the records again.');
    default:
      return t('The records or their orientation changed since the last Generate. Generate, then orient the records again.');
  }
};

/** @param {RecordAutoOrientationOptions} options */
export const createRecordAutoOrientation = ({
  computed, ref, getCommittedSession, readCommittedRevision, busy, decideOrientation, writeDecidedOrientations,
  errorMessage
}) => {
  const pending = ref(false);
  const status = ref('');
  const available = computed(() => {
    void readCommittedRevision();
    return hasNucleotideComparison(getCommittedSession());
  });
  const enabled = computed(() => available.value && !pending.value);

  const orient = async () => {
    const session = getCommittedSession();
    if (pending.value || !hasNucleotideComparison(session)) return;
    const blocked = busy();
    if (blocked) { status.value = tm(blocked.reason); return; }
    const committed = /** @type {CommittedCanonicalSession} */ (session);
    pending.value = true;
    status.value = t('Reading the comparisons of the last Generate...');
    try {
      const reply = readAutoOrientationReply(
        (await decideOrientation({
          canonicalRequest: committed.renderRequest, resources: committed.resources
        }))?.result,
        () => diagnosticError('RESULT_INVALID', {}, ORIENTATION_FAILURE)
      );
      // A Generate that committed meanwhile drew other records or comparisons.
      const written = getCommittedSession() === committed
        ? await writeDecidedOrientations(reply.records, msg('Orient records from comparisons'))
        : /** @type {RecordOrientationWrite} */ ({ status: 'stale' });
      status.value = writeStatusText(written);
    } catch (error) {
      status.value = errorMessage(error);
    } finally {
      pending.value = false;
    }
  };

  return { available, enabled, pending, status, orient };
};
