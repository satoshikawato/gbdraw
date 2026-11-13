// @ts-check
// The Review alignment panel and the Reset alignment dialog as display values: the "After Apply"
// strips (where each record, its features and its candidates lie once the chosen genes are lined
// up), the cards in their sections, the direction preview's sentences, and the Reset dialog's
// direction rows. Pure; the figure geometry is measured by the preview's owner
// (`app/feature-editor/svg-actions.js`), and the choices, the direction projection and the Reset
// receipt come from the alignment owner (`app/similarity-alignment.js`).
import { joinSentences, listText, msg, numberText, t, tm, tn, tx } from '../utils/i18n.js';
import { canonicalRecordReverseComplement } from './record-display-model.js';
/** @import { CanonicalRecordOrientation } from './record-display-model.js' */
/** @import { AlignmentResetReceipt } from './session-active-config-contract.js' */

/**
 * A horizontal span in figure units (SVG user units of the mounted Result,
 * which do not change with the preview zoom).
 * @typedef {{ left: number, right: number }} FigureSpan
 */

/**
 * One record as the mounted Result draws it now.
 * @typedef {object} AlignmentRecordGeometry
 * @property {string} recordKey
 * @property {FigureSpan} extent The record's axis and every drawn feature.
 * @property {FigureSpan[]} features Every drawn feature block.
 * @property {{ key: string, span: FigureSpan }[]} candidates The drawn candidates, by candidate key. A
 *   candidate drawn in several parts spans from its first block to its last.
 */

/**
 * @typedef {object} AlignmentFigureGeometry
 * @property {number} referenceCenter The reference feature's centre now.
 * @property {AlignmentRecordGeometry[]} records
 */

/**
 * What the draft says about one record.
 * @typedef {object} AlignmentStripRecord
 * @property {string} recordKey
 * @property {string} lineKey The key of the candidate that Apply puts on the reference line
 *   (the reference feature in its own record); '' when the record stays where it is.
 * @property {number | null} mirrorShift For a record Apply reverses: how far the centre of the
 *   `lineKey` anchor moves inside its record (Python's `afterCenterX - beforeCenterX`); null otherwise.
 */

/** A box on a strip, in strip pixels. @typedef {{ x: number, width: number }} StripBox */

/**
 * @typedef {object} AlignmentStripRow
 * @property {string} recordKey
 * @property {StripBox} baseline The record's extent.
 * @property {StripBox[]} features
 * @property {(StripBox & { key: string })[]} candidates
 */

/**
 * @typedef {object} AlignmentStripLayout
 * @property {number} width
 * @property {number} lineX The reference line, at the same x in every strip.
 * @property {number} scale Strip pixels per figure unit.
 * @property {AlignmentStripRow[]} rows In the order of `records`; a record without geometry has none.
 */

// The strips are drawn in a box this wide and stretched to the panel's width.
const STRIP_WIDTH = 350;
const MIN_FEATURE_WIDTH = 0.6;
const MIN_CANDIDATE_WIDTH = 4;

/** @param {FigureSpan} span */
const centreOf = (span) => (span.left + span.right) / 2;

/**
 * @param {FigureSpan} span
 * @param {number | null} axis
 * @returns {FigureSpan}
 */
const mirrored = (span, axis) => (
  axis === null ? span : { left: 2 * axis - span.right, right: 2 * axis - span.left }
);

/**
 * @param {{ geometry: AlignmentFigureGeometry | null, records: readonly AlignmentStripRecord[], width: number }} options
 * @returns {AlignmentStripLayout | null}
 */
export const layoutAlignmentStrips = ({ geometry, records, width }) => {
  if (!geometry || !(width > 0)) return null;
  const geometryByKey = new Map(geometry.records.map((record) => [record.recordKey, record]));
  const placed = records.flatMap(({ recordKey, lineKey, mirrorShift }) => {
    const drawn = geometryByKey.get(recordKey);
    if (!drawn) return [];
    const onLine = lineKey ? drawn.candidates.find(({ key }) => key === lineKey) : undefined;
    // Reverse complement mirrors a record about its axis centre. The projection gives how far
    // the anchor on the line moves, so the axis lies half that shift from the anchor; without a
    // drawn anchor it is the middle of the record's extent.
    const axis = mirrorShift === null ? null
      : onLine ? centreOf(onLine.span) + mirrorShift / 2 : centreOf(drawn.extent);
    const extent = mirrored(drawn.extent, axis);
    // A record with nothing on the line keeps its present offset from the reference feature.
    const lineAt = onLine ? centreOf(mirrored(onLine.span, axis)) : geometry.referenceCenter;
    return [{
      recordKey, extent, lineAt,
      features: drawn.features.map((span) => mirrored(span, axis)),
      candidates: drawn.candidates.map(({ key, span }) => ({ key, span: mirrored(span, axis) }))
    }];
  });
  const leftOfLine = Math.max(0, ...placed.map(({ extent, lineAt }) => lineAt - extent.left));
  const rightOfLine = Math.max(0, ...placed.map(({ extent, lineAt }) => extent.right - lineAt));
  if (!placed.length || !(leftOfLine + rightOfLine > 0)) return null;
  const scale = width / (leftOfLine + rightOfLine);
  const lineX = leftOfLine * scale;
  /**
   * @param {FigureSpan} span
   * @param {number} lineAt
   * @param {number} minWidth
   * @returns {StripBox}
   */
  const box = (span, lineAt, minWidth) => {
    const boxWidth = Math.max(minWidth, (span.right - span.left) * scale);
    return { x: lineX + (centreOf(span) - lineAt) * scale - boxWidth / 2, width: boxWidth };
  };
  return {
    width, lineX, scale,
    rows: placed.map(({ recordKey, extent, lineAt, features, candidates }) => ({
      recordKey,
      baseline: box(extent, lineAt, 0),
      features: features.map((span) => box(span, lineAt, MIN_FEATURE_WIDTH)),
      candidates: candidates.map(({ key, span }) => ({ key, ...box(span, lineAt, MIN_CANDIDATE_WIDTH) }))
    }))
  };
};

/** The key of the reference feature among a record's drawn candidates. */
export const ALIGNMENT_REFERENCE_KEY = 'reference';

/** @typedef {'reference' | 'choose' | 'several' | 'single' | 'none'} AlignmentReviewGroupId */

// The reference card comes first, without a heading; then the sections.
/** @type {readonly { id: AlignmentReviewGroupId, title: string, count: (records: number) => string }[]} */
const GROUPS = Object.freeze([
  { id: 'reference', title: '', count: () => '' },
  { id: 'choose', title: msg('Choose a gene'), count: (n) => tn(n, '{count} record awaiting a choice', '{count} records awaiting a choice') },
  { id: 'several', title: msg('Several candidates'), count: (n) => tn(n, '{count} record with several candidates', '{count} records with several candidates') },
  { id: 'single', title: msg('One candidate'), count: (n) => tn(n, '{count} record with one candidate', '{count} records with one candidate') },
  { id: 'none', title: msg('No usable candidate'), count: (n) => tn(n, '{count} record without a usable candidate', '{count} records without a usable candidate') }
]);

/**
 * A draft row, as far as the panel orders and marks it.
 * @typedef {object} ReviewRowLike
 * @property {string} recordKey
 * @property {readonly { key: string }[]} candidates
 * @property {{ kind: string, candidateKey?: string | null } | null} choice
 */

/**
 * One candidate of the figure overlay: a ring on the Result.
 * @typedef {object} AlignmentOverlayCandidate
 * @property {string} key
 * @property {Record<string, unknown>} anchor The feature identity the preview's owner finds the block by.
 * @property {string} name Shown in the tag above the ring.
 * @property {string} accessibleName The ring button's name.
 */

/**
 * One record of the figure overlay.
 * @typedef {object} AlignmentOverlayRecord
 * @property {string} recordKey
 * @property {string} chosenKey The candidate the draft lines up; '' when none.
 * @property {boolean} needsChoice
 * @property {boolean} skipped
 * @property {AlignmentOverlayCandidate[]} candidates
 */

/** A record needs a choice until the user selects a candidate or skips it. @param {{ choice: unknown }} row */
const needsChoice = ({ choice }) => choice === null;

/** @param {ReviewRowLike['choice']} choice */
const chosenKeyOf = (choice) => (choice?.kind === 'select' ? String(choice.candidateKey ?? '') : '');

/**
 * The draft rows as the figure overlay marks them.
 * @param {readonly { recordKey: string, recordLabel: string, choice: ReviewRowLike['choice'], candidates: readonly {
 *   key: string, anchor: Record<string, unknown>, name: string, label: string, coordinates: string }[] }[]} rows
 * @returns {AlignmentOverlayRecord[]}
 */
export const alignmentOverlayRecords = (rows) => rows.map((row) => ({
  recordKey: row.recordKey,
  chosenKey: chosenKeyOf(row.choice),
  needsChoice: needsChoice(row),
  skipped: row.choice?.kind === 'skip',
  candidates: row.candidates.map(({ key, anchor, name, label, coordinates }) => ({
    key, anchor, name, accessibleName: t('Select {label}, {coordinates} for {record}', { label, coordinates, record: row.recordLabel })
  }))
}));

/**
 * A record of the direction preview, in figure order.
 * @typedef {object} ReviewDirectionLike
 * @property {string} recordKey
 * @property {string} status
 * @property {boolean} reversed
 * @property {number | null} [beforeCenterX]
 * @property {number | null} [afterCenterX]
 */

/**
 * @template {ReviewRowLike} Row
 * @template {ReviewDirectionLike} Direction
 * @typedef {object} AlignmentReviewCard
 * @property {string} recordKey
 * @property {'reference' | 'several' | 'single' | 'none'} kind
 * @property {Direction} direction
 * @property {Row | null} row Null in the reference card.
 * @property {Row['candidates'][number][]} candidates Left to right as the figure draws them; undrawn ones last.
 * @property {boolean} needsChoice
 * @property {boolean} skipped
 * @property {string} chosenKey The candidate Apply puts on the reference line; '' when none.
 * @property {AlignmentStripRow | null} strip
 */

/**
 * @template {ReviewRowLike} Row
 * @template {ReviewDirectionLike} Direction
 * @typedef {object} AlignmentReviewView
 * @property {number} stripWidth
 * @property {number | null} lineX Null when the figure gives no strips.
 * @property {{ id: AlignmentReviewGroupId, title: string, cards: AlignmentReviewCard<Row, Direction>[] }[]} groups
 *   Only the groups that have a record, records in figure order.
 * @property {string} summary The record count and the count of each section, in words.
 */

/**
 * The keys of the records without a choice.
 * @param {readonly { recordKey: string, choice: unknown }[]} rows
 * @returns {string[]}
 */
export const recordsNeedingChoice = (rows) => rows.filter(needsChoice).map(({ recordKey }) => recordKey);

/**
 * @template {ReviewRowLike} Row
 * @template {ReviewDirectionLike} Direction
 * @param {object} options
 * @param {readonly Row[]} options.rows
 * @param {readonly Direction[]} options.directions
 * @param {AlignmentFigureGeometry | null} options.geometry
 * @param {readonly string[]} options.openedWithoutChoice The records without a choice when the review opened:
 *   they stay under "Choose a gene" and keep their radio rows while the user chooses.
 * @param {number} [options.width] The strip width.
 * @returns {AlignmentReviewView<Row, Direction>}
 */
export const buildAlignmentReviewView = ({ rows, directions, geometry, openedWithoutChoice, width = STRIP_WIDTH }) => {
  const choose = new Set(openedWithoutChoice);
  const rowByKey = new Map(rows.map((row) => [row.recordKey, row]));
  const drawnByKey = new Map((geometry?.records || []).map((record) => [record.recordKey, record]));
  /** @type {{ group: AlignmentReviewGroupId, card: AlignmentReviewCard<Row, Direction> }[]} */
  const placed = [];
  directions.forEach((direction) => {
    const { recordKey } = direction;
    const row = rowByKey.get(recordKey);
    if (direction.status === 'reference') {
      placed.push({ group: 'reference', card: { recordKey, direction, strip: null, kind: 'reference', row: null,
        candidates: [], needsChoice: false, skipped: false, chosenKey: ALIGNMENT_REFERENCE_KEY } });
      return;
    }
    if (!row) return;
    const drawnLeft = new Map(drawnByKey.get(recordKey)?.candidates.map(({ key, span }) => [key, span.left]));
    /** @type {Row['candidates'][number][]} */
    const candidates = [
      ...row.candidates.filter(({ key }) => drawnLeft.has(key))
        .sort((a, b) => Number(drawnLeft.get(a.key)) - Number(drawnLeft.get(b.key))),
      ...row.candidates.filter(({ key }) => !drawnLeft.has(key))
    ];
    const several = candidates.length > 1 || (candidates.length === 1 && choose.has(recordKey));
    placed.push({
      group: choose.has(recordKey) ? 'choose'
        : candidates.length > 1 ? 'several' : candidates.length ? 'single' : 'none',
      card: { recordKey, direction, strip: null, row, candidates,
        kind: several ? 'several' : candidates.length ? 'single' : 'none',
        needsChoice: needsChoice(row),
        skipped: row.choice?.kind === 'skip',
        chosenKey: chosenKeyOf(row.choice) }
    });
  });
  const layout = layoutAlignmentStrips({
    geometry, width,
    records: placed.map(({ card: { recordKey, chosenKey, direction } }) => ({
      recordKey, lineKey: chosenKey,
      mirrorShift: direction.reversed ? (Number(direction.afterCenterX) - Number(direction.beforeCenterX)) || 0 : null
    }))
  });
  const stripByKey = new Map(layout?.rows.map((row) => [row.recordKey, row]));
  placed.forEach(({ card }) => { card.strip = stripByKey.get(card.recordKey) ?? null; });
  const groups = GROUPS.map(({ id, title, count }) => ({
    id, title, count, cards: placed.filter(({ group }) => group === id).map(({ card }) => card)
  })).filter(({ cards }) => cards.length > 0);
  return {
    stripWidth: width,
    lineX: layout?.lineX ?? null,
    groups: groups.map(({ id, title, cards }) => ({ id, title: tm(title), cards })),
    // The list join stays English: it needs a locale list pattern.
    summary: tn(placed.length, '{count} record', '{count} records')
      + groups.filter(({ title }) => title).map(({ count, cards }) => `, ${count(cards.length)}`).join('') + '.'
  };
};

/** A feature's arrow as the figure draws it: right-facing, left-facing, or unknown. @typedef {-1 | 1 | null} DisplayedStrand */

/**
 * A record of the alignment owner's direction projection, as the direction preview reads it.
 * @typedef {object} ProjectedDirectionLike
 * @property {string} recordKey
 * @property {string} status
 * @property {string | null} exclusion Why Apply leaves the record's direction (or the record) as it is.
 * @property {DisplayedStrand} beforeArrow
 * @property {DisplayedStrand} afterArrow
 * @property {boolean} beforeReverseComplement
 * @property {boolean} afterReverseComplement
 */

/**
 * One record of the Review alignment direction preview, as the cards read it.
 * @typedef {object} DirectionPreviewRecord
 * @property {string} recordKey
 * @property {string} status
 * @property {string | null} exclusion
 * @property {DisplayedStrand} afterArrow
 * @property {boolean} beforeReverseComplement
 * @property {boolean} afterReverseComplement
 * @property {boolean} reversed
 * @property {string} shortLabel
 * @property {string} label
 * @property {string} description
 */

/**
 * The names of a record of the draft: its label, its accession, its definition.
 * @typedef {{ recordLabel?: string, recordName?: string, recordDetail?: string }} DraftRecordNames
 */

/** @param {unknown} arrow */
const arrowLabel = (arrow) => arrow === 1 ? t('right-facing') : arrow === -1 ? t('left-facing') : t('unknown');

/** @type {Readonly<Record<string, string>>} */
const directionExclusionLabels = Object.freeze({
  selection_required: msg('Select an anchor or Skip this record first.'),
  skipped_by_user: msg('Skipped by user.'),
  skipped_no_candidate: msg('No usable candidate in this record.'),
  skipped_unmappable: msg('Candidate center is outside the displayed crop.'),
  unusable_anchor: msg('The selected anchor cannot be displayed.'),
  unknown_strand: msg('The selected anchor has no known direction; it is never guessed.')
});

/** @param {string} code */
const exclusionReason = (code) => tm(directionExclusionLabels[code] || code);

/**
 * One record of the direction preview: its names and its before/after statement.
 * @template {ProjectedDirectionLike} R
 * @param {R} record
 * @param {DraftRecordNames & { recordLabel: string }} reference The exact reference's record.
 * @param {readonly (DraftRecordNames & { recordKey: string })[]} rows
 */
const directionRecordView = (record, reference, rows) => {
  const row = rows.find((entry) => entry.recordKey === record.recordKey);
  const isReference = record.status === 'reference';
  const label = isReference ? reference.recordLabel : row?.recordLabel || record.recordKey;
  const beforeDirection = arrowLabel(record.beforeArrow);
  const afterDirection = arrowLabel(record.afterArrow);
  // An anchor of unknown strand still moves its record; only the direction stays (PD-OI-027).
  const exclusionLabel = record.exclusion
    ? (record.exclusion === 'unknown_strand'
      ? t('Direction unchanged: {reason}', { reason: exclusionReason(record.exclusion) })
      : t('Unchanged: {reason}', { reason: exclusionReason(record.exclusion) }))
    : '';
  return { ...record, label, beforeDirection, afterDirection, exclusionLabel,
    shortLabel: (isReference ? reference.recordDetail || reference.recordName : row?.recordDetail || row?.recordName) || label,
    reversed: record.beforeReverseComplement !== record.afterReverseComplement,
    // The accessible before/after statement each card and its Custom select describe.
    description: (isReference
      ? t('{label} · exact reference: Current {before}; after Align {after}.', { label, before: beforeDirection, after: afterDirection })
      : t('{label}: Current {before}; after Align {after}.', { label, before: beforeDirection, after: afterDirection }))
      + (exclusionLabel ? ` ${exclusionLabel}` : '') };
};

/**
 * What Apply does to the records' directions, in sentences.
 * @param {readonly DirectionPreviewRecord[]} records
 * @param {{ mode: string }} intent The Review panel's direction choice.
 */
const directionEffectText = (records, intent) => {
  const eligible = records.filter(record => !record.exclusion);
  const reversed = records.filter(record => record.reversed);
  // A record still waiting for Select or Skip is counted by the footer, not here.
  const excluded = records.filter(record => record.exclusion
    && !['selection_required', 'unknown_strand'].includes(record.exclusion)).length;
  const unknownStrand = records.filter(record => record.exclusion === 'unknown_strand').length;
  const reversedText = tn(reversed.length, '{count} record ({names})', '{count} records ({names})',
    { names: listText(reversed.map(({ shortLabel }) => shortLabel)) });
  const sameAfter = eligible.length > 0 && eligible.every(({ afterArrow }) => afterArrow === eligible[0].afterArrow)
    ? (eligible[0].afterArrow === 1 ? 'right' : eligible[0].afterArrow === -1 ? 'left' : '') : '';
  let effect;
  if (intent.mode === 'custom') {
    effect = reversed.length ? t('Reverses {records}.', { records: reversedText }) : t('No record is reversed yet. Choose a direction in each card.');
  } else if (reversed.length) {
    // The direction is part of each whole sentence: a bare "right" or "left" is no key.
    effect = intent.mode === 'left'
      ? tn(eligible.length, 'Reverses {records} so the chosen feature points left.',
        'Reverses {records} so the {count} chosen features point left.', { records: reversedText })
      : tn(eligible.length, 'Reverses {records} so the chosen feature points right.',
        'Reverses {records} so the {count} chosen features point right.', { records: reversedText });
  } else if (sameAfter) {
    effect = joinSentences(t('No record is reversed.'), sameAfter === 'left'
      ? tn(eligible.length, 'The chosen feature already points left.', 'The {count} chosen features already point left.')
      : tn(eligible.length, 'The chosen feature already points right.', 'The {count} chosen features already point right.'));
  } else {
    effect = joinSentences(t('No record is reversed.'),
      tn(eligible.length, '{count} chosen feature keeps its current direction.', '{count} chosen features keep their current directions.'));
  }
  // Verb agreement lives in the two forms of one key.
  if (excluded) {
    effect = joinSentences(effect, tn(excluded, '{count} record unchanged.',
      '{count} records unchanged.'));
  }
  if (unknownStrand) {
    effect = joinSentences(effect, tn(unknownStrand, '{count} record keeps its direction.',
      '{count} records keep their direction.'));
  }
  return effect;
};

/**
 * The note on the reference record's left edge: the reference feature's centre stays fixed (D-11).
 * @param {string} name The reference feature's name.
 * @param {number} deltaX How far the reference record moves, in figure units at 100% zoom.
 */
const referenceShiftText = (name, deltaX) => {
  const shift = Math.round(deltaX);
  if (!shift) return '';
  return shift > 0
    ? t('{name} stays where it is; the record\'s left edge moves {px} px (at 100% zoom) to the right.', { name, px: numberText(Math.abs(shift)) })
    : t('{name} stays where it is; the record\'s left edge moves {px} px (at 100% zoom) to the left.', { name, px: numberText(Math.abs(shift)) });
};

/**
 * The direction preview of the Review panel: the projection with each record's names and
 * statement, the effect of Apply in sentences, and the reference record's shift.
 * @template {{ records: readonly ProjectedDirectionLike[], reference: { deltaX: number } }} P
 * @param {object} options
 * @param {P} options.projected The alignment owner's projection of the draft.
 * @param {DraftRecordNames & { recordLabel: string, name: string }} options.reference The draft's exact reference.
 * @param {readonly (DraftRecordNames & { recordKey: string })[]} options.rows The draft's rows.
 * @param {{ mode: string }} options.intent The Review panel's direction choice.
 */
export const describeDirectionPreview = ({ projected, reference, rows, intent }) => {
  /** @type {DirectionPreviewRecord[]} */
  const records = projected.records.map((record) => directionRecordView(record, reference, rows));
  return { ...projected, effect: directionEffectText(records, intent),
    referenceShift: referenceShiftText(reference.name, projected.reference.deltaX), records };
};

/** @param {boolean} reversed Whether a record reads opposite to its source file. */
const recordDirectionWord = (reversed) => reversed ? tx('record direction', 'reverse') : tx('record direction', 'forward');

/**
 * A record's name for one-line summaries: its accession, else its definition.
 * @param {{ name?: string, detail?: string } | undefined} names
 */
export const shortRecordName = (names) => names?.detail || names?.name || '';

/**
 * The Reset alignment dialog's direction option: the records the latest Align reversed, each
 * with its direction now and the one Reset restores, and why the option is unavailable.
 * @param {object} options
 * @param {AlignmentResetReceipt | null} options.receipt
 * @param {{ records: readonly (CanonicalRecordOrientation & { recordKey: string })[] } | null} options.request
 *   The committed render request.
 * @param {ReadonlyMap<string, { label: string, name: string, detail: string }>} options.names Record names by record key.
 */
export const resetPreviewView = ({ receipt, request, names }) => {
  const targets = (receipt?.directions || []).map((delta) => {
    // The receipt is read against the committed request, so each record it names is there.
    const record = /** @type {CanonicalRecordOrientation} */ (
      request?.records.find(({ recordKey }) => recordKey === delta.recordKey));
    const recordName = names.get(delta.recordKey);
    const current = canonicalRecordReverseComplement(record);
    return { recordKey: delta.recordKey, label: recordName?.label || t('Unavailable record'),
      name: shortRecordName(recordName) || t('Unavailable record'),
      current, restored: delta.before,
      currentWord: recordDirectionWord(current), restoredWord: recordDirectionWord(delta.before),
      laterManualEdit: current !== delta.after };
  });
  const count = targets.length;
  return { targets,
    directionsLabel: count ? tn(count, 'Positions and directions ({count} record)', 'Positions and directions ({count} records)')
      : t('Positions and directions'),
    disabledReason: !receipt
      ? t('Not available: this Session does not record which directions the latest Align changed.')
      : !count ? t('Not available: the latest Align did not change any record direction.') : '' };
};
