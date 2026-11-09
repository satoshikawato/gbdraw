// @ts-check
// The feature popup's payload (`state.clickedFeature`): what the popup shows
// of one feature and the edits it holds while open. The SVG owner
// (app/feature-editor/svg-actions.js) opens the popup with the values it reads
// off the Result and the drawing; this module builds the payload from them.
/** @import { FeatureLike } from './feature-utils.js' */
import {
  formatFeatureLength,
  formatFeatureLocation,
  getFeatureCaption,
  getQualifierFirstValue,
  normalizeStringArray,
  resolveFeatureProteinId
} from './feature-utils.js';
import { t } from '../utils/i18n.js';

/**
 * A qualifier of the Qualifiers tab: its name and values.
 * @typedef {{ key: string, values: string[] }} FeaturePopupQualifierRow
 */
/**
 * A row of the Details tab.
 * @typedef {{ key: string, label: string, value: string }} FeaturePopupDetailRow
 */
/**
 * The feature popup's tabs (the rich popup shows all four).
 * @typedef {'edit' | 'details' | 'qualifiers' | 'sequence'} FeaturePopupTab
 */

/**
 * The Sequence tab's read (D5): nothing read yet (`idle`), Python reading it (`reading`), read (`ready`, with the
 * FASTA texts and warnings), or failed (`failed`, with the error the tab shows with Details and Retry). The catalog
 * stores no sequence; the popup opens without one and reads it when the tab is shown.
 * @typedef {object} FeaturePopupSequence
 * @property {'idle' | 'reading' | 'ready' | 'failed'} status
 * @property {readonly string[]} warnings
 * @property {string} nucleotideFasta
 * @property {string} aminoAcidFasta
 * @property {unknown} error The caught failure (`normalizeCaughtError` reads it), or null.
 */

/** @returns {FeaturePopupSequence} */
export const unreadFeatureSequence = () => ({ status: 'idle', warnings: [], nucleotideFasta: '', aminoAcidFasta: '', error: null });

/**
 * The open feature popup. The feature's own fields (`feat` and the values
 * copied from it without conversion) are as the feature row carries them;
 * the label fields are filled by the label owner once the popup opens.
 * @typedef {object} FeaturePopup
 * @property {unknown} id
 * @property {string} svg_id The SVG ID the displayed Result draws it with ('' when it does not).
 * @property {string} rendered_svg_id
 * @property {unknown} stable_svg_id
 * @property {string} label The feature's caption.
 * @property {string} location
 * @property {unknown[]} locationParts
 * @property {string} color
 * @property {FeatureLike} feat
 * @property {FeaturePopupTab} activeTab
 * @property {string} recordId
 * @property {string} gene
 * @property {string} product
 * @property {number | null} recordIdx
 * @property {string} featureType
 * @property {number | null} start
 * @property {number | null} end
 * @property {string} strand
 * @property {object} qualifiers
 * @property {FeaturePopupQualifierRow[]} qualifierRows
 * @property {FeaturePopupSequence} sequence
 * @property {FeaturePopupDetailRow[]} detailRows
 * @property {string} legendName The Legend name the popup's field holds.
 * @property {string} appliedLegendName The Legend name the drawing applies.
 * @property {string} strokeColor
 * @property {number | ''} strokeWidth '' after a stroke reset when neither the Legend row nor the drawn feature has a width.
 * @property {string | null} strokeWidthDraft The stroke width typed or stepped into the popup's field and not committed yet; null: none, the field shows `strokeWidth`.
 * @property {string} labelKey
 * @property {string} labelText
 * @property {string} labelSourceText
 * @property {string} labelVisibility
 * @property {string} featureVisibility
 * @property {boolean} identityEditable Per-feature edits name the feature by its source identity (design Q4).
 * @property {unknown} proteinId
 * @property {unknown} sourceProteinId
 * @property {unknown} orthogroupId
 * @property {unknown} orthogroupMemberCount
 * @property {unknown} orthogroupRecordCoverage
 * @property {boolean} orthogroupRepresentative
 * @property {unknown} orthogroupMember
 * @property {boolean} hasEditableLabel
 * @property {string} labelUnavailableReason
 */

/**
 * What the SVG owner reads for the popup off the Result and the drawing.
 * @typedef {object} FeaturePopupView
 * @property {string} svgId The SVG ID the displayed Result draws the feature with ('' when it does not).
 * @property {string} legendName The feature's Legend name.
 * @property {string} color The drawn fill.
 * @property {string} strokeColor
 * @property {number} strokeWidth
 * @property {string} featureVisibility The feature's visibility override.
 * @property {boolean} identityEditable
 */

/**
 * @param {unknown} qualifiers
 * @returns {FeaturePopupQualifierRow[]}
 */
const normalizeQualifierRows = (qualifiers) => {
  if (!qualifiers || typeof qualifiers !== 'object' || Array.isArray(qualifiers)) return [];
  return Object.entries(qualifiers)
    .map(([key, value]) => {
      const values = normalizeStringArray(value);
      return { key: String(key || ''), values };
    })
    .filter((row) => row.key && row.values.length > 0)
    .sort((left, right) => left.key.localeCompare(right.key));
};

// Details lists what the popup header and the Edit tab do not show:
// the header names the feature, its product, record, and location, and the
// Similarity group section its group.
/**
 * @param {FeatureLike} feat
 * @returns {FeaturePopupDetailRow[]}
 */
const buildDetailRows = (feat) => [
  { key: 'type', label: t('Feature type'), value: feat.type },
  { key: 'length', label: t('Length'), value: formatFeatureLength(feat) },
  { key: 'protein_id', label: t('Protein ID'), value: resolveFeatureProteinId(feat, /** @type {Record<string, unknown> | null} */ (feat?.orthogroupMember || feat?.orthogroup_member || null)) }
]
  .map((row) => ({ ...row, value: row.value === null || row.value === undefined ? '' : String(row.value) }))
  .filter((row) => row.value !== '');

/**
 * The payload of the feature popup opened on `feat`.
 * @param {FeatureLike} feat
 * @param {FeaturePopupView} view
 * @returns {FeaturePopup}
 */
export const buildFeaturePopup = (feat, view) => {
  const defaultLabel = getFeatureCaption(feat);
  const locationText = formatFeatureLocation(feat);
  const locationParts = Array.isArray(feat.location_parts) ? feat.location_parts : [];
  const qualifierRows = normalizeQualifierRows(feat.qualifiers);

  return {
    id: feat.id,
    svg_id: view.svgId,
    rendered_svg_id: view.svgId,
    stable_svg_id: feat.stable_svg_id || feat.stable_feature_id || feat.svg_id || '',
    label: defaultLabel,
    location: locationText,
    locationParts,
    color: view.color,
    feat,
    activeTab: 'edit',
    recordId: String(feat.record_id || ''),
    gene: getQualifierFirstValue(feat, 'gene'),
    product: getQualifierFirstValue(feat, 'product'),
    recordIdx: Number.isInteger(Number(feat.record_idx)) ? Number(feat.record_idx) : null,
    featureType: String(feat.type || ''),
    start: Number.isFinite(Number(feat.start)) ? Number(feat.start) : null,
    end: Number.isFinite(Number(feat.end)) ? Number(feat.end) : null,
    strand: String(feat.strand || ''),
    qualifiers: feat.qualifiers && typeof feat.qualifiers === 'object' ? feat.qualifiers : {},
    qualifierRows,
    sequence: unreadFeatureSequence(),
    detailRows: buildDetailRows(feat),
    legendName: view.legendName,
    appliedLegendName: view.legendName,
    strokeColor: view.strokeColor,
    strokeWidth: view.strokeWidth,
    strokeWidthDraft: null,
    labelKey: '',
    labelText: '',
    labelSourceText: '',
    labelVisibility: 'default',
    featureVisibility: view.featureVisibility,
    identityEditable: view.identityEditable,
    proteinId: feat.proteinId || feat.protein_id || '',
    sourceProteinId: feat.sourceProteinId || feat.source_protein_id || '',
    orthogroupId: feat.orthogroupId || '',
    orthogroupMemberCount: feat.orthogroupMemberCount || 0,
    orthogroupRecordCoverage: feat.orthogroupRecordCoverage || 0,
    orthogroupRepresentative: Boolean(feat.orthogroupRepresentative),
    orthogroupMember: feat.orthogroupMember || null,
    hasEditableLabel: false,
    labelUnavailableReason: t('No editable feature label for this feature in current diagram.')
  };
};
