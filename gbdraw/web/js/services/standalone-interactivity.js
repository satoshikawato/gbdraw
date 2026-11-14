// @ts-check
/** @import { FeatureCatalogItem, FeatureCatalogRow } from './feature-catalog.js' */
/** @import { EmbeddedSequencesJson } from './catalog-sequences.js' */
/** @import { AliasSource } from './feature-identity.js' */
import { normalizeStringArray } from './feature-utils.js';
import { ownEntry } from '../utils/own-entry.js';
import { INTERACTIVE_EMBED_SCHEMA } from '../contracts/schema.generated.js';
import { STANDALONE_INTERACTIVE_SCRIPT, STANDALONE_INTERACTIVE_STYLE } from './standalone-interactivity-assets.js';
import { ensureSvgDefs } from './svg-serialization.js';
import { ANNOTATION_SELECTOR } from './drawn-groups.js';
import {
  MATCH_FRAGMENT_OWN_ATTRIBUTES, MATCH_SPAN_REASON_TEXTS, STANDALONE_MATCH_FRAGMENT_OWN_ATTRIBUTES_SLOT,
  STANDALONE_MATCH_SPAN_REASON_TEXTS_SLOT
} from './drawn-attributes.generated.js';

const SVG_NS = 'http://www.w3.org/2000/svg';
const FEATURE_ID_ATTRIBUTE = 'data-gbdraw-feature-id';
const RENDERED_FEATURE_ID_ATTRIBUTE = 'data-gbdraw-rendered-feature-id';
const FEATURE_SELECTOR = [
  `path[${FEATURE_ID_ATTRIBUTE}]`,
  `polygon[${FEATURE_ID_ATTRIBUTE}]`,
  `rect[${FEATURE_ID_ATTRIBUTE}]`,
  'path[id^="f"]',
  'polygon[id^="f"]',
  'rect[id^="f"]'
].join(', ');
const STANDALONE_MATCH_SELECTOR = [
  'path[data-gbdraw-match-id]',
  'path[data-match-kind]',
  'path[data-pairwise-match-style]'
].join(', ');
const INTERACTIVE_METADATA_ID = 'gbdraw-interactive-feature-metadata';
const INTERACTIVE_STYLE_ID = 'gbdraw-interactive-feature-style';
const INTERACTIVE_SCRIPT_ID = 'gbdraw-interactive-feature-script';
const INTERACTIVE_GLOW_FILTER_ID = 'gbdraw-interactive-feature-glow';
const INTERACTIVE_MATCH_GLOW_FILTER_ID = 'gbdraw-interactive-feature-match-glow';
const FEATURE_PART_SUFFIX_RE = /__part\d+$/;

/**
 * A label the Features list edits; its fields are only read through `String()`.
 * @typedef {object} StandaloneEditableLabel
 * @property {unknown} [featureId]
 * @property {unknown} [text]
 * @property {unknown} [sourceText]
 */

/**
 * What the caller of the Interactive SVG export hands over: the export
 * (`result-actions.js`) gives the embed item it read, and the capture
 * (`captureSvgExport` in `svg-serialization.js`) the popup mode and the
 * drawing's edits. Every field is optional and is normalized by
 * `normalizeStandaloneContext`.
 * @typedef {object} StandaloneInteractivityOptions
 * @property {string} [popupMode] `simple` or `rich`; anything else is `rich`.
 * @property {FeatureCatalogItem | null} [catalogItem] The shown Result's stored catalog item, which the file
 *   embeds with `embeddedSequences`; without both the export writes nothing.
 * @property {EmbeddedSequencesJson | null} [embeddedSequences] The item's sequences in the export encoding, as
 *   Python's encoder writes them (`embed_catalog_part`, which the CLI export writes too), read through the Worker
 *   (`CatalogSequences['readEmbedJson']`, D5): the biological rows with their sequences (inline, or as
 *   `sequenceSourceIndex` and the `translationFromAminoAcidSequence` mark) and the whole-record sources the item's
 *   matches read. The page derives none of these and writes them into the file as they are.
 * @property {StandaloneEditableLabel[]} [editableLabels]
 * @property {Record<string, unknown>} [labelTextFeatureOverrides] Label text by rendered Feature ID.
 * @property {Record<string, unknown>} [labelTextBulkOverrides] Label text by source text.
 * @property {Record<string, unknown>} [orthogroupNameOverrides]
 * @property {Record<string, unknown>} [orthogroupDescriptionOverrides]
 */

/**
 * The options after `normalizeStandaloneContext`: every field present.
 * @typedef {object} StandaloneContext
 * @property {'simple' | 'rich'} popupMode
 * @property {FeatureCatalogItem | null} catalogItem
 * @property {EmbeddedSequencesJson | null} embeddedSequences
 * @property {StandaloneEditableLabel[]} editableLabels
 * @property {Record<string, unknown>} labelTextFeatureOverrides
 * @property {Record<string, unknown>} labelTextBulkOverrides
 * @property {Record<string, unknown>} orthogroupNameOverrides
 * @property {Record<string, unknown>} orthogroupDescriptionOverrides
 */

/** @typedef {{ valid: boolean, value: string }} StandaloneTextAlias */
/** @typedef {{ x: number, y: number, width: number, height: number }} StandaloneViewBox */

/** @param {unknown} value */
const normalizeFeatureElementId = (value) =>
  String(value || '').trim().replace(FEATURE_PART_SUFFIX_RE, '');

/** @param {Element | null | undefined} element */
const getElementFeatureId = (element) =>
  normalizeFeatureElementId(
    element?.getAttribute?.(RENDERED_FEATURE_ID_ATTRIBUTE) ||
    element?.getAttribute?.(FEATURE_ID_ATTRIBUTE) ||
    element?.getAttribute?.('id') ||
    ''
  );

/**
 * @param {AliasSource | null | undefined} feature
 * @param {string} key
 * @returns {string}
 */
const firstQualifierValue = (feature, key) => {
  // `typeof` rules out the other values; an array reads as an empty map.
  const qualifiers = /** @type {Record<string, unknown>} */ (
    feature?.qualifiers && typeof feature.qualifiers === 'object'
      ? feature.qualifiers
      : {}
  );
  const values = normalizeStringArray(qualifiers[key]);
  return values.find((value) => value.trim()) || '';
};

/** @param {AliasSource | null | undefined} feature */
const getStandaloneFeatureLabel = (feature) => {
  const candidates = [
    feature?.label,
    feature?.gene,
    feature?.locus_tag,
    firstQualifierValue(feature, 'gene'),
    firstQualifierValue(feature, 'locus_tag'),
    firstQualifierValue(feature, 'product'),
    feature?.product,
    feature?.type,
    feature?.svg_id
  ];
  for (const candidate of candidates) {
    const label = String(candidate || '').trim();
    if (label) return label;
  }
  return 'Feature';
};

/** @param {unknown} value */
const normalizeFeatureIdKey = (value) => String(value || '').trim().toLowerCase();

/**
 * @param {AliasSource | null | undefined} payload
 * @param {readonly string[]} keys
 * @returns {StandaloneTextAlias}
 */
const standaloneConsistentTextAlias = (payload, keys) => {
  const values = new Set(keys
    .filter((key) => Object.prototype.hasOwnProperty.call(payload || {}, key))
    .map((key) => String(payload?.[key] ?? '').trim())
    .filter(Boolean));
  return {
    valid: values.size <= 1,
    value: values.size === 1 ? values.values().next().value : ''
  };
};

/**
 * @param {StandaloneInteractivityOptions} [options]
 * @returns {StandaloneContext}
 */
const normalizeStandaloneContext = (options = {}) => ({
  popupMode: normalizeStandalonePopupMode(options.popupMode),
  catalogItem: options.catalogItem && typeof options.catalogItem === 'object' ? options.catalogItem : null,
  embeddedSequences: options.embeddedSequences && typeof options.embeddedSequences === 'object'
    ? options.embeddedSequences : null,
  editableLabels: Array.isArray(options.editableLabels) ? options.editableLabels : [],
  labelTextFeatureOverrides:
    options.labelTextFeatureOverrides && typeof options.labelTextFeatureOverrides === 'object'
      ? options.labelTextFeatureOverrides
      : {},
  labelTextBulkOverrides:
    options.labelTextBulkOverrides && typeof options.labelTextBulkOverrides === 'object'
      ? options.labelTextBulkOverrides
      : {},
  orthogroupNameOverrides:
    options.orthogroupNameOverrides && typeof options.orthogroupNameOverrides === 'object'
      ? options.orthogroupNameOverrides
      : {},
  orthogroupDescriptionOverrides:
    options.orthogroupDescriptionOverrides && typeof options.orthogroupDescriptionOverrides === 'object'
      ? options.orthogroupDescriptionOverrides
      : {}
});

/**
 * @param {FeatureCatalogItem} item
 * @returns {FeatureCatalogItem}
 */
const cloneStandaloneCatalogItem = (item) => JSON.parse(JSON.stringify(item));

/**
 * @param {unknown} recordKey
 * @param {unknown} biologicalFeatureId
 * @returns {string}
 */
const standaloneCatalogFeatureKey = (recordKey, biologicalFeatureId) => (
  `${String(recordKey || '').trim()}\u0000${String(biologicalFeatureId || '').trim()}`
);

/**
 * @param {AliasSource | null | undefined} feature
 * @param {StandaloneContext} context
 * @returns {StandaloneEditableLabel | null}
 */
const getEditableLabelEntryForStandaloneFeature = (feature, context) => {
  const entries = Array.isArray(context.editableLabels) ? context.editableLabels : [];
  const candidates = [
    normalizeFeatureIdKey(feature?.svg_id),
    normalizeFeatureIdKey(feature?.id)
  ].filter(Boolean);
  if (candidates.length === 0) return null;
  return entries.find((entry) => candidates.includes(normalizeFeatureIdKey(entry?.featureId))) || null;
};

/**
 * @param {AliasSource | null | undefined} feature
 * @param {StandaloneContext} context
 * @returns {string}
 */
const getLabelTextFeatureOverride = (feature, context) => {
  const overrides = context.labelTextFeatureOverrides;
  const candidates = [
    String(feature?.svg_id || '').trim(),
    String(feature?.id || '').trim()
  ].filter(Boolean);
  for (const candidate of candidates) {
    if (Object.prototype.hasOwnProperty.call(overrides, candidate)) {
      const text = String(overrides[candidate] ?? '').trim();
      if (text) return text;
    }
  }
  return '';
};

/**
 * @param {AliasSource | null | undefined} feature
 * @param {string} fallbackLabel
 * @param {StandaloneContext} context
 * @returns {string}
 */
const getStandaloneDisplayLabel = (feature, fallbackLabel, context) => {
  const editableEntry = getEditableLabelEntryForStandaloneFeature(feature, context);
  const editableText = String(editableEntry?.text || '').trim();
  if (editableText) return editableText;

  const directOverride = getLabelTextFeatureOverride(feature, context);
  if (directOverride) return directOverride;

  const sourceCandidates = [
    editableEntry?.sourceText,
    fallbackLabel,
    feature?.product,
    feature?.gene,
    feature?.locus_tag
  ].map((value) => String(value || '').trim()).filter(Boolean);
  for (const sourceText of sourceCandidates) {
    const bulkText = String(ownEntry(context.labelTextBulkOverrides, sourceText) ?? '').trim();
    if (bulkText) return bulkText;
  }

  return fallbackLabel;
};

/**
 * @param {FeatureCatalogItem} sourceItem
 * @param {StandaloneContext} context
 * @returns {FeatureCatalogItem | null}
 */
const catalogItemWithStandaloneOverrides = (sourceItem, context) => {
  const item = cloneStandaloneCatalogItem(sourceItem);
  const biologicalByKey = new Map();
  for (const feature of (Array.isArray(item.biologicalFeatures) ? item.biologicalFeatures : [])) {
    const recordKey = String(feature?.recordKey || '').trim();
    const biologicalFeatureId = String(feature?.biologicalFeatureId || '').trim();
    if (!recordKey || !biologicalFeatureId) return null;
    const key = standaloneCatalogFeatureKey(recordKey, biologicalFeatureId);
    if (biologicalByKey.has(key)) return null;
    biologicalByKey.set(key, feature);
  }
  const renderedKeys = new Set();
  const renderedIds = new Set();
  // A null entry makes the check below return null, which drops this local clone.
  item.features = /** @type {FeatureCatalogRow[]} */ ((Array.isArray(item.features) ? item.features : []).map((reference) => {
    const recordKey = String(reference?.recordKey || '').trim();
    const biologicalFeatureId = String(reference?.biologicalFeatureId || '').trim();
    if (!recordKey || !biologicalFeatureId) return null;
    const key = standaloneCatalogFeatureKey(recordKey, biologicalFeatureId);
    const renderedId = String(reference?.svgId || '').trim();
    if (!key || !renderedId || renderedKeys.has(key) || renderedIds.has(renderedId)) return null;
    const biological = biologicalByKey.get(key);
    if (!biological) return null;
    renderedKeys.add(key);
    renderedIds.add(renderedId);
    const feature = {
      ...biological,
      id: key,
      svg_id: renderedId
    };
    const fallbackLabel = getStandaloneFeatureLabel(feature);
    const displayLabel = getStandaloneDisplayLabel(feature, fallbackLabel, context);
    return displayLabel && displayLabel !== fallbackLabel
      ? { ...reference, displayLabel }
      : reference;
  }));
  if (item.features.some((reference) => !reference)) return null;
  const groupIds = new Set();
  // A null entry makes the check below return null, which drops this local clone.
  item.orthogroups = /** @type {FeatureCatalogRow[]} */ ((Array.isArray(item.orthogroups) ? item.orthogroups : [])
    .map((group) => {
      const identity = standaloneConsistentTextAlias(
        group,
        ['id', 'orthogroupId', 'orthogroup_id']
      );
      const id = identity.valid ? identity.value : '';
      if (!id || groupIds.has(id)) return null;
      groupIds.add(id);
      const displayName = id
        ? String(ownEntry(context.orthogroupNameOverrides, id) || '').trim()
        : '';
      const description = id
        ? String(ownEntry(context.orthogroupDescriptionOverrides, id) || '').trim()
        : '';
      return {
        ...group,
        ...(displayName ? { display_name: displayName } : {}),
        ...(description ? { description } : {})
      };
    }));
  if (item.orthogroups.some((group) => !group)) return null;
  return item;
};

/** @param {unknown} popupMode */
const normalizeStandalonePopupMode = (popupMode) => (
  popupMode === 'simple' ? 'simple' : 'rich'
);

/**
 * @param {Element | null | undefined} element
 * @param {string} name
 */
const standaloneAttr = (element, name) =>
  String(element?.getAttribute?.(name) || '').trim();

/**
 * @param {SVGDefsElement} defs
 * @param {{ id: string, color: string, opacity: number, blurStdDeviation: number, slope: number, extent: number }} options
 */
const appendStandaloneFeatureGlowFilter = (
  defs,
  {
    id,
    color,
    opacity,
    blurStdDeviation,
    slope,
    extent
  }
) => {
  const existing = defs.querySelector(`#${CSS.escape(id)}`);
  if (existing?.parentNode) {
    existing.parentNode.removeChild(existing);
  }

  const filter = document.createElementNS(SVG_NS, 'filter');
  filter.setAttribute('id', id);
  filter.setAttribute('x', `-${extent}%`);
  filter.setAttribute('y', `-${extent}%`);
  filter.setAttribute('width', `${100 + (extent * 2)}%`);
  filter.setAttribute('height', `${100 + (extent * 2)}%`);
  filter.setAttribute('color-interpolation-filters', 'sRGB');

  const componentTransfer = document.createElementNS(SVG_NS, 'feComponentTransfer');
  componentTransfer.setAttribute('in', 'SourceGraphic');
  componentTransfer.setAttribute('result', 'gbdrawBrightenedFeature');
  ['R', 'G', 'B'].forEach((channel) => {
    const func = document.createElementNS(SVG_NS, `feFunc${channel}`);
    func.setAttribute('type', 'linear');
    func.setAttribute('slope', String(slope));
    componentTransfer.appendChild(func);
  });

  const blur = document.createElementNS(SVG_NS, 'feGaussianBlur');
  blur.setAttribute('in', 'SourceAlpha');
  blur.setAttribute('stdDeviation', String(blurStdDeviation));
  blur.setAttribute('result', 'gbdrawFeatureGlowBlur');

  const flood = document.createElementNS(SVG_NS, 'feFlood');
  flood.setAttribute('flood-color', color);
  flood.setAttribute('flood-opacity', String(opacity));
  flood.setAttribute('result', 'gbdrawFeatureGlowColor');

  const composite = document.createElementNS(SVG_NS, 'feComposite');
  composite.setAttribute('in', 'gbdrawFeatureGlowColor');
  composite.setAttribute('in2', 'gbdrawFeatureGlowBlur');
  composite.setAttribute('operator', 'in');
  composite.setAttribute('result', 'gbdrawFeatureGlow');

  const merge = document.createElementNS(SVG_NS, 'feMerge');
  ['gbdrawFeatureGlow', 'gbdrawBrightenedFeature'].forEach((resultName) => {
    const mergeNode = document.createElementNS(SVG_NS, 'feMergeNode');
    mergeNode.setAttribute('in', resultName);
    merge.appendChild(mergeNode);
  });

  filter.appendChild(componentTransfer);
  filter.appendChild(blur);
  filter.appendChild(flood);
  filter.appendChild(composite);
  filter.appendChild(merge);
  defs.appendChild(filter);
};

/** @param {SVGSVGElement} svg */
const ensureStandaloneFeatureGlowFilter = (svg) => {
  const defs = ensureSvgDefs(svg);
  appendStandaloneFeatureGlowFilter(defs, {
    id: INTERACTIVE_GLOW_FILTER_ID,
    color: '#2563eb',
    opacity: 0.85,
    blurStdDeviation: 3,
    slope: 1.2,
    extent: 35
  });
  appendStandaloneFeatureGlowFilter(defs, {
    id: INTERACTIVE_MATCH_GLOW_FILTER_ID,
    color: '#fbbf24',
    opacity: 0.32,
    blurStdDeviation: 1.5,
    slope: 1.04,
    extent: 25
  });
};

/** @param {Element} svg */
const removeExistingStandaloneInteractivityAssets = (svg) => {
  [
    INTERACTIVE_METADATA_ID,
    INTERACTIVE_STYLE_ID,
    INTERACTIVE_SCRIPT_ID,
    INTERACTIVE_GLOW_FILTER_ID,
    INTERACTIVE_MATCH_GLOW_FILTER_ID,
    'gbdraw-viewport-controls',
    'gbdraw-feature-search-controls',
    'gbdraw-feature-popup',
    'gbdraw-feature-hover-popup'
  ].forEach((id) => {
    const element = svg.querySelector(`#${CSS.escape(id)}`);
    if (element?.parentNode) {
      element.parentNode.removeChild(element);
    }
  });
};

/**
 * @param {Element} element
 * @param {string} token
 */
const addClassToken = (element, token) => {
  const existing = String(element.getAttribute('class') || '').trim();
  const tokens = new Set(existing ? existing.split(/\s+/) : []);
  tokens.add(token);
  element.setAttribute('class', Array.from(tokens).join(' '));
};

/**
 * @param {Element} element
 * @param {string} token
 */
const removeClassToken = (element, token) => {
  const existing = String(element.getAttribute('class') || '').trim();
  const tokens = existing ? existing.split(/\s+/).filter((entry) => entry !== token) : [];
  if (tokens.length) element.setAttribute('class', tokens.join(' '));
  else element.removeAttribute('class');
};

/**
 * @param {unknown} value
 * @returns {StandaloneViewBox | null}
 */
const parseSvgViewBox = (value) => {
  const parts = String(value || '')
    .trim()
    .split(/[\s,]+/)
    .map((part) => Number(part));
  if (parts.length < 4 || parts.some((part) => !Number.isFinite(part))) return null;
  if (parts[2] <= 0 || parts[3] <= 0) return null;
  return {
    x: parts[0],
    y: parts[1],
    width: parts[2],
    height: parts[3]
  };
};

/** @param {StandaloneViewBox} rect */
const formatSvgViewBox = (rect) =>
  [rect.x, rect.y, rect.width, rect.height].map((value) => String(value)).join(' ');

/** @param {Element} svg */
const resolveStandaloneOriginalGeometry = (svg) => {
  const dataViewBox = parseSvgViewBox(svg.getAttribute('data-gbdraw-original-viewbox'));
  const currentViewBox = parseSvgViewBox(svg.getAttribute('viewBox'));
  const width = parseFloat(svg.getAttribute('width') ?? '');
  const height = parseFloat(svg.getAttribute('height') ?? '');
  const fallbackViewBox =
    Number.isFinite(width) && Number.isFinite(height) && width > 0 && height > 0
      ? { x: 0, y: 0, width, height }
      : null;
  const viewBox = dataViewBox || currentViewBox || fallbackViewBox || { x: 0, y: 0, width: 900, height: 650 };
  const originalWidth = svg.getAttribute('data-gbdraw-original-width') ||
    svg.getAttribute('width') ||
    `${viewBox.width}px`;
  const originalHeight = svg.getAttribute('data-gbdraw-original-height') ||
    svg.getAttribute('height') ||
    `${viewBox.height}px`;
  return {
    viewBox,
    originalWidth,
    originalHeight
  };
};

/** @param {SVGSVGElement} svg */
const applyStandaloneViewportRoot = (svg) => {
  svg.setAttribute('width', '100vw');
  svg.setAttribute('height', '100vh');
  svg.setAttribute('preserveAspectRatio', 'xMidYMid meet');
  svg.style.setProperty('width', '100vw');
  svg.style.setProperty('height', '100vh');
  svg.style.setProperty('display', 'block');
  svg.style.setProperty('background', '#ffffff');
};

/**
 * The JSON text of the item an exported Interactive SVG embeds: `JSON.stringify` of the stored item with the
 * drawing's edits, whose biological rows, whole-record sources and match span reads are the export encoding Python
 * wrote, taken as text (the item Python's `embed_catalog_item` writes: the rows in their place, the sources and then
 * the reads after every stored key).
 * @param {FeatureCatalogItem} item
 * @param {EmbeddedSequencesJson} sequences
 * @returns {string}
 */
export const embeddedCatalogItemJson = (item, sequences) => {
  const { sequenceSources, matchSpanReads } = sequences;
  /** @type {Record<string, unknown>} */
  const fields = {
    ...item, biologicalFeatures: null,
    ...(sequenceSources === null ? {} : { sequenceSources: null }), ...(matchSpanReads === null ? {} : { matchSpanReads: null })
  };
  const text = (/** @type {string} */ key, /** @type {unknown} */ value) => (key === 'biologicalFeatures'
    ? sequences.biologicalFeatures
    : key === 'sequenceSources' && sequenceSources !== null ? sequenceSources
      : key === 'matchSpanReads' && matchSpanReads !== null ? matchSpanReads : JSON.stringify(value));
  return `{${Object.entries(fields).filter(([, value]) => value !== undefined && typeof value !== 'function')
    .map(([key, value]) => `${JSON.stringify(key)}:${text(key, value)}`).join(',')}}`;
};

/**
 * @param {SVGSVGElement | null | undefined} svg
 * @param {StandaloneInteractivityOptions} [options]
 * @returns {boolean}
 */
export const enrichSvgWithStandaloneInteractivity = (svg, options = {}) => {
  if (!svg) return false;

  const context = normalizeStandaloneContext(options);
  // The file embeds the stored item with the drawing's edits and the sequences Python encoded (D5); without both it
  // writes nothing.
  const { embeddedSequences } = context;
  if (!context.catalogItem || !embeddedSequences) return false;
  const catalogItem = catalogItemWithStandaloneOverrides(context.catalogItem, context);
  if (!catalogItem) return false;
  const originalGeometry = resolveStandaloneOriginalGeometry(svg);
  const normalizedPopupMode = normalizeStandalonePopupMode(context.popupMode);
  removeExistingStandaloneInteractivityAssets(svg);
  removeClassToken(svg, 'gbdraw-feature-search-active');
  removeClassToken(svg, 'gbdraw-feature-search-updating');
  svg.querySelectorAll('.gbdraw-interactive-feature--dimmed, .gbdraw-interactive-feature--match, .gbdraw-interactive-feature--active-match')
    .forEach((element) => {
      removeClassToken(element, 'gbdraw-interactive-feature--dimmed');
      removeClassToken(element, 'gbdraw-interactive-feature--match');
      removeClassToken(element, 'gbdraw-interactive-feature--active-match');
    });
  svg.querySelectorAll('[data-gbdraw-interactive-annotation]').forEach((element) => {
    element.removeAttribute('data-gbdraw-interactive-annotation');
    removeClassToken(element, 'gbdraw-interactive-annotation');
  });
  const featureIds = new Set(catalogItem.features.map((feature) => feature.svgId));
  svg.querySelectorAll(FEATURE_SELECTOR).forEach((element) => {
    const id = getElementFeatureId(element);
    if (!featureIds.has(id)) return;
    element.setAttribute('data-gbdraw-interactive-feature', 'true');
    addClassToken(element, 'gbdraw-interactive-feature');
  });
  const matchIds = new Set(catalogItem.comparisonMatches.map((match) => match.id));
  svg.querySelectorAll(STANDALONE_MATCH_SELECTOR).forEach((element) => {
    if (!matchIds.has(standaloneAttr(element, 'data-gbdraw-match-id'))) return;
    element.setAttribute('data-gbdraw-interactive-match', 'true');
    addClassToken(element, 'gbdraw-interactive-pairwise-match');
  });
  const annotationDomIds = new Set(
    (Array.isArray(catalogItem.annotations) ? catalogItem.annotations : [])
      .map((annotation) => String(annotation?.dom_id || annotation?.domId || '').trim())
      .filter(Boolean)
  );
  svg.querySelectorAll(ANNOTATION_SELECTOR).forEach((element) => {
    if (!annotationDomIds.has(String(element.id || '').trim())) return;
    element.setAttribute('data-gbdraw-interactive-annotation', 'true');
    addClassToken(element, 'gbdraw-interactive-annotation');
  });

  const metadata = document.createElementNS(SVG_NS, 'metadata');
  metadata.setAttribute('id', INTERACTIVE_METADATA_ID);
  metadata.setAttribute('data-popup-mode', normalizedPopupMode);
  // The export profile (D5-02): the embedded catalog's schema, which the CLI writes too, not the stored one.
  metadata.setAttribute('data-schema', String(INTERACTIVE_EMBED_SCHEMA));
  metadata.setAttribute('data-result-index', String(catalogItem.resultIndex));
  metadata.setAttribute('data-result-name', String(catalogItem.resultName || ''));
  metadata.textContent = `{"schema":${INTERACTIVE_EMBED_SCHEMA},"items":[${embeddedCatalogItemJson(catalogItem, embeddedSequences)}]}`;

  const style = document.createElementNS(SVG_NS, 'style');
  style.setAttribute('id', INTERACTIVE_STYLE_ID);
  style.setAttribute('type', 'text/css');
  style.textContent = STANDALONE_INTERACTIVE_STYLE;

  const script = document.createElementNS(SVG_NS, 'script');
  script.setAttribute('id', INTERACTIVE_SCRIPT_ID);
  script.setAttribute('type', 'application/ecmascript');
  // The script checks the fragments of a match with Python's list of the attributes a fragment holds on its own, and
  // shows Python's sentence for each reason a match span cannot be read; Python writes both in the same way.
  script.textContent = [
    [STANDALONE_MATCH_FRAGMENT_OWN_ATTRIBUTES_SLOT, JSON.stringify(MATCH_FRAGMENT_OWN_ATTRIBUTES)],
    [STANDALONE_MATCH_SPAN_REASON_TEXTS_SLOT, JSON.stringify(MATCH_SPAN_REASON_TEXTS)]
  ].reduce((text, [slot, json]) => {
    const parts = text.split(slot);
    if (parts.length !== 2) throw new Error('Malformed standalone interactive SVG asset: STANDALONE_INTERACTIVE_SCRIPT');
    return parts.join(json);
  }, STANDALONE_INTERACTIVE_SCRIPT);

  ensureStandaloneFeatureGlowFilter(svg);
  svg.appendChild(metadata);
  svg.appendChild(style);
  svg.appendChild(script);
  svg.setAttribute('data-gbdraw-original-viewbox', formatSvgViewBox(originalGeometry.viewBox));
  svg.setAttribute('data-gbdraw-original-width', originalGeometry.originalWidth);
  svg.setAttribute('data-gbdraw-original-height', originalGeometry.originalHeight);
  svg.setAttribute('data-gbdraw-interactive-svg', 'true');
  applyStandaloneViewportRoot(svg);
  return true;
};
