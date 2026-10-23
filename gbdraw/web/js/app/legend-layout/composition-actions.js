// @ts-check
import {
  prependTranslate,
  readLeadingTranslate,
  replaceLeadingTranslate
} from '../../services/svg-transform.js';
import { planLegendComposition } from '../../services/legend-layout.js';
import { legendPositionOfMode, readLegendPositionOfMode } from '../../services/layout-preferences.js';
import { LEGEND_POSITIONS } from '../../contracts/schema.generated.js';
import { LENGTH_BAR_ID, PLOT_TITLE_ID } from '../../services/drawn-groups.js';
/** @import { LayoutBox, CompositionSpacingPx, CompositionOverlayPolicy } from '../../services/legend-layout.js' */
/** @import { LoadRepair } from '../../utils/load-repairs.js' */

/**
 * Bounds as the metadata and the editor hold them.
 * @typedef {{ x: number, y: number, width: number, height: number }} CompositionBounds
 */
/** @typedef {'primary' | 'legend' | 'title'} CompositionRole */
/**
 * @typedef {object} CompositionTargetBase
 * @property {number[]} automaticTranslation
 * @property {CompositionRole} role
 * @property {string} selector
 */
/** @typedef {CompositionTargetBase & { finalBounds: CompositionBounds }} CompositionPrimaryTarget */
/** @typedef {CompositionTargetBase & { localBounds: CompositionBounds }} CompositionLocalTarget */
/**
 * The inputs Python laid the Legend out with exist together or not at all.
 * @typedef {object} CompositionLegendReflow
 * @property {number} colorRectSize
 * @property {number} lineHeight
 * @property {number} textXOffset
 * @property {number} [dpi]
 * @property {string | null} [fontFile]
 * @property {number} [fontSize]
 * @property {LayoutBox} [primaryLocalBounds]
 * @property {number} [wrapWidth]
 */
/**
 * @typedef {object} CompositionMetadata
 * @property {CompositionLocalTarget | null} legend
 * @property {CompositionLegendReflow | null} legendReflow
 * @property {string} legendSide
 * @property {CompositionBounds[]} overlayObstacles
 * @property {CompositionOverlayPolicy} overlayPolicy
 * @property {CompositionPrimaryTarget} primary
 * @property {CompositionSpacingPx} spacing
 * @property {CompositionLocalTarget | null} title
 * @property {string} titleSide
 * @property {true} [legacyNormalized]
 */
/**
 * The SVG groups one role is bound to, with the metadata that names them.
 * @typedef {{ metadata: CompositionTargetBase | null, targets: SVGGraphicsElement[] }} CompositionRoleBinding
 */
/**
 * @typedef {object} CompositionBinding
 * @property {CompositionMetadata} metadata
 * @property {CompositionRoleBinding} primary
 * @property {CompositionRoleBinding} legend
 * @property {CompositionRoleBinding} title
 */
/** @typedef {[number, number]} CompositionDelta */
/**
 * @typedef {object} CompositionMatrix
 * @property {number} a
 * @property {number} b
 * @property {number} c
 * @property {number} d
 * @property {number} e
 * @property {number} f
 */
/**
 * The saved user offsets of an older Session, as the Session config passes
 * them; `primary` is one offset, or one per primary group.
 * @typedef {object} LegacyUserDeltas
 * @property {number[] | number[][] | null} [primary]
 * @property {number[] | null} [lengthBar]
 * @property {number[] | null} [legend]
 * @property {number[] | null} [title]
 */
/**
 * @typedef {object} CompositionEditOptions
 * @property {string} [legendSide]
 * @property {string} [titleSide]
 * @property {LayoutBox | null} [legendLocalBox]
 * @property {CompositionBounds | null} [legendLocalBounds]
 * @property {CompositionBounds | null} [titleLocalBounds]
 * @property {Partial<Record<string, unknown>> | null} [canvasPadding]
 */
/**
 * @typedef {object} LegacyTargets
 * @property {Element[]} primary
 * @property {Element | null} legend
 * @property {Element | null} title
 */
/**
 * One alignment-record decision of the active plan; only the fields read here.
 * @typedef {object} LegacyRecordDecision
 * @property {unknown} [recordKey]
 * @property {unknown} [status]
 * @property {{ stableFeatureSvgId?: unknown, biologicalFeatureId?: unknown } | null} [anchor]
 */

export const COMPOSITION_SCHEMA_VERSION = 1;
export const COMPOSITION_SCHEMA_ATTRIBUTE = 'data-gbdraw-composition-schema';
export const COMPOSITION_METADATA_ATTRIBUTE = 'data-gbdraw-composition';
export const COMPOSITION_ROLE_ATTRIBUTE = 'data-gbdraw-composition-role';

const ROLE_SELECTORS = Object.freeze({
  primary: '[data-gbdraw-composition-role="primary"]',
  legend: '[data-gbdraw-composition-role="legend"]',
  title: '[data-gbdraw-composition-role="title"]'
});
// The legend sides of every mode. A loaded Result's side is one of its own
// mode's (`normalizeStoredLegendSide`); the checks below hold no mode.
/** @type {Set<unknown>} */
const LEGEND_SIDES = new Set([...LEGEND_POSITIONS.circular, ...LEGEND_POSITIONS.linear]);
/** @type {Set<unknown>} */
const TITLE_SIDES = new Set(['top', 'bottom', 'center', 'none']);
/** @type {Set<unknown>} */
const OVERLAY_CANDIDATE_SCORE_NAMES = new Set([
  'totalAnchorDistance',
  'xAnchorDistance',
  'yAnchorDistance',
  'nearEdgeX',
  'nearEdgeY'
]);
/** @type {Set<unknown>} */
const OVERLAY_CANVAS_GROWTH_CANDIDATE_NAMES = new Set(['horizontal', 'vertical']);
/** @type {Set<unknown>} */
const OVERLAY_CANVAS_GROWTH_SCORE_NAMES = new Set([
  'addedArea',
  'addedExtent',
  'candidateOrder'
]);
// Schema-free SVGs predate Python-emitted spacing, legend metrics, and overlay
// policy. Keep this frozen snapshot inside the one explicit legacy adapter;
// current SVGs never consume it.
const LEGACY_COMPOSITION_SPACING_V0 = Object.freeze({
  dockGapPx: 24,
  edgePaddingPx: 16,
  overlayClearancePx: 8,
  stackGapPx: 20,
  titleGapPx: 20
});
const LEGEND_LAYOUT_INPUT_FIELDS = Object.freeze([
  'dpi',
  'fontFile',
  'fontSize',
  'primaryLocalBounds',
  'wrapWidth'
]);
const LEGACY_LEGEND_LINE_HEIGHT_RATIO_V0 = 24 / 14;
const LEGACY_LEGEND_TEXT_OFFSET_RATIO_V0 = 22 / 14;
const LEGACY_OVERLAY_POLICY_V0 = Object.freeze({
  candidateScoreOrder: Object.freeze([
    'totalAnchorDistance',
    'xAnchorDistance',
    'yAnchorDistance',
    'nearEdgeX',
    'nearEdgeY'
  ]),
  canvasGrowthCandidateOrder: Object.freeze(['horizontal', 'vertical']),
  canvasGrowthScoreOrder: Object.freeze(['addedArea', 'addedExtent', 'candidateOrder']),
  quadrantBoundaryRatio: 0.5
});

export class CompositionMetadataError extends Error {
  /**
   * @param {string} message
   * @param {string} [code]
   */
  constructor(message, code = 'INVALID_COMPOSITION_METADATA') {
    super(message);
    this.name = 'CompositionMetadataError';
    this.code = code;
  }
}

/** @type {(message: string, code?: string) => never} */
const fail = (message, code) => {
  throw new CompositionMetadataError(message, code);
};

/**
 * @param {unknown} value
 * @returns {value is Record<string, unknown>}
 */
const isPlainObject = (value) => (
  Boolean(value) && typeof value === 'object' && !Array.isArray(value)
);

/**
 * @param {unknown} value
 * @param {string} path
 */
const finiteNumber = (value, path) => {
  if (typeof value !== 'number' || !Number.isFinite(value)) {
    fail(`${path} must be a finite number.`);
  }
  return value;
};

/**
 * @param {unknown} value
 * @param {string} path
 */
const nonNegativeNumber = (value, path) => {
  const number = finiteNumber(value, path);
  if (number < 0) fail(`${path} must be non-negative.`);
  return number;
};

/**
 * @param {unknown} value
 * @param {string} path
 */
const positiveNumber = (value, path) => {
  const number = finiteNumber(value, path);
  if (number <= 0) fail(`${path} must be positive.`);
  return number;
};

/**
 * @param {unknown} value
 * @param {string} path
 * @param {{ positive?: boolean }} [options]
 * @returns {CompositionBounds}
 */
const validateBounds = (value, path, { positive = false } = {}) => {
  if (!isPlainObject(value)) fail(`${path} must be an object.`);
  const bounds = {
    x: finiteNumber(value.x, `${path}.x`),
    y: finiteNumber(value.y, `${path}.y`),
    width: nonNegativeNumber(value.width, `${path}.width`),
    height: nonNegativeNumber(value.height, `${path}.height`)
  };
  if (positive && (bounds.width <= 0 || bounds.height <= 0)) {
    fail(`${path} must have positive width and height.`);
  }
  return bounds;
};

/**
 * Python's Aabb field for field (min/max form, not x/width).
 * @param {unknown} value
 * @param {string} path
 * @returns {LayoutBox}
 */
const validateCornerBounds = (value, path) => {
  if (!isPlainObject(value)) fail(`${path} must be an object.`);
  const bounds = {
    minX: finiteNumber(value.minX, `${path}.minX`),
    minY: finiteNumber(value.minY, `${path}.minY`),
    maxX: finiteNumber(value.maxX, `${path}.maxX`),
    maxY: finiteNumber(value.maxY, `${path}.maxY`)
  };
  if (bounds.maxX < bounds.minX || bounds.maxY < bounds.minY) {
    fail(`${path} must not be inverted.`);
  }
  return bounds;
};

/**
 * @param {unknown} value
 * @param {string} path
 * @returns {CompositionDelta}
 */
const validateTranslation = (value, path) => {
  if (!Array.isArray(value) || value.length !== 2) {
    fail(`${path} must contain exactly two numbers.`);
  }
  return [finiteNumber(value[0], `${path}[0]`), finiteNumber(value[1], `${path}[1]`)];
};

/**
 * The result carries the one bounds field named by `boundsKey`; callers narrow
 * it to the primary or local target type.
 * @param {unknown} value
 * @param {CompositionRole} role
 * @param {'finalBounds' | 'localBounds'} boundsKey
 * @returns {CompositionTargetBase & { finalBounds?: CompositionBounds, localBounds?: CompositionBounds }}
 */
const validateTarget = (value, role, boundsKey) => {
  if (!isPlainObject(value)) fail(`composition.${role} must be an object.`);
  if (value.role !== role) fail(`composition.${role}.role must be "${role}".`);
  if (value.selector !== ROLE_SELECTORS[role]) {
    fail(`composition.${role}.selector is not the schema 1 ${role} selector.`);
  }
  return {
    automaticTranslation: validateTranslation(
      value.automaticTranslation,
      `composition.${role}.automaticTranslation`
    ),
    [boundsKey]: validateBounds(
      value[boundsKey],
      `composition.${role}.${boundsKey}`,
      { positive: role === 'primary' }
    ),
    role,
    selector: ROLE_SELECTORS[role]
  };
};

/**
 * @param {unknown} value
 * @param {'legend' | 'title'} role
 * @returns {CompositionLocalTarget | null}
 */
const validateOptionalTarget = (value, role) => (
  value === null
    ? null
    : /** @type {CompositionLocalTarget} */ (validateTarget(value, role, 'localBounds')) // carries localBounds
);

/**
 * @param {unknown} value
 * @returns {CompositionSpacingPx}
 */
const validateSpacing = (value) => {
  if (!isPlainObject(value)) fail('composition.spacing must be an object.');
  return {
    dockGapPx: nonNegativeNumber(value.dockGapPx, 'composition.spacing.dockGapPx'),
    edgePaddingPx: nonNegativeNumber(value.edgePaddingPx, 'composition.spacing.edgePaddingPx'),
    overlayClearancePx: nonNegativeNumber(
      value.overlayClearancePx,
      'composition.spacing.overlayClearancePx'
    ),
    stackGapPx: nonNegativeNumber(value.stackGapPx, 'composition.spacing.stackGapPx'),
    titleGapPx: nonNegativeNumber(value.titleGapPx, 'composition.spacing.titleGapPx')
  };
};

/**
 * @param {unknown} value
 * @param {string} path
 * @param {Set<unknown>} allowed
 * @returns {string[]}
 */
const validateOrder = (value, path, allowed) => {
  if (!Array.isArray(value)) fail(`${path} must be an array.`);
  if (
    value.length !== allowed.size ||
    value.some((name) => typeof name !== 'string' || !allowed.has(name)) ||
    new Set(value).size !== value.length
  ) {
    fail(`${path} must contain each supported policy name exactly once.`);
  }
  return [...value];
};

/**
 * @param {unknown} value
 * @returns {CompositionOverlayPolicy}
 */
const validateOverlayPolicy = (value) => {
  if (!isPlainObject(value)) fail('composition.overlayPolicy must be an object.');
  const requiredFields = [
    'candidateScoreOrder',
    'canvasGrowthCandidateOrder',
    'canvasGrowthScoreOrder',
    'quadrantBoundaryRatio'
  ];
  if (
    Object.keys(value).length !== requiredFields.length ||
    requiredFields.some((field) => !Object.prototype.hasOwnProperty.call(value, field))
  ) {
    fail('composition.overlayPolicy has an invalid field set.');
  }
  return {
    candidateScoreOrder: validateOrder(
      value.candidateScoreOrder,
      'composition.overlayPolicy.candidateScoreOrder',
      OVERLAY_CANDIDATE_SCORE_NAMES
    ),
    canvasGrowthCandidateOrder: validateOrder(
      value.canvasGrowthCandidateOrder,
      'composition.overlayPolicy.canvasGrowthCandidateOrder',
      OVERLAY_CANVAS_GROWTH_CANDIDATE_NAMES
    ),
    canvasGrowthScoreOrder: validateOrder(
      value.canvasGrowthScoreOrder,
      'composition.overlayPolicy.canvasGrowthScoreOrder',
      OVERLAY_CANVAS_GROWTH_SCORE_NAMES
    ),
    quadrantBoundaryRatio: (() => {
      const ratio = finiteNumber(
        value.quadrantBoundaryRatio,
        'composition.overlayPolicy.quadrantBoundaryRatio'
      );
      if (ratio <= 0 || ratio >= 1) {
        fail('composition.overlayPolicy.quadrantBoundaryRatio must be between zero and one.');
      }
      return ratio;
    })()
  };
};

/**
 * @param {unknown} value
 * @param {{ required: boolean }} options
 * @returns {CompositionLegendReflow | null}
 */
const validateLegendReflow = (value, { required }) => {
  if (value === undefined) {
    fail(
      required
        ? 'composition.legendReflow is required when a legend is present.'
        : 'composition.legendReflow must be null when no legend is present.'
    );
  }
  if (value === null) {
    if (required) fail('composition.legendReflow is required when a legend is present.');
    return null;
  }
  if (!required) fail('composition.legendReflow requires a legend target.');
  if (!isPlainObject(value)) fail('composition.legendReflow must be an object.');
  const reflow = {
    colorRectSize: positiveNumber(
      value.colorRectSize,
      'composition.legendReflow.colorRectSize'
    ),
    lineHeight: positiveNumber(value.lineHeight, 'composition.legendReflow.lineHeight'),
    textXOffset: positiveNumber(
      value.textXOffset,
      'composition.legendReflow.textXOffset'
    )
  };
  // The inputs Python laid the Legend out with. SVGs written before Python
  // emitted them carry none; parsing keeps them so a composition edit writes
  // them back unchanged.
  const present = LEGEND_LAYOUT_INPUT_FIELDS.filter((field) => (
    Object.prototype.hasOwnProperty.call(value, field)
  ));
  if (present.length === 0) return reflow;
  if (present.length !== LEGEND_LAYOUT_INPUT_FIELDS.length) {
    fail('composition.legendReflow has an incomplete set of Python layout inputs.');
  }
  if (!Number.isInteger(value.dpi) || /** @type {number} */ (value.dpi) <= 0) { // an integer, so a number
    fail('composition.legendReflow.dpi must be a positive integer.');
  }
  if (value.fontFile !== null && (typeof value.fontFile !== 'string' || !value.fontFile.trim())) {
    fail('composition.legendReflow.fontFile must be a non-empty string or null.');
  }
  return {
    ...reflow,
    dpi: /** @type {number} */ (value.dpi), // checked just above
    fontFile: value.fontFile,
    fontSize: positiveNumber(value.fontSize, 'composition.legendReflow.fontSize'),
    primaryLocalBounds: validateCornerBounds(
      value.primaryLocalBounds,
      'composition.legendReflow.primaryLocalBounds'
    ),
    wrapWidth: nonNegativeNumber(value.wrapWidth, 'composition.legendReflow.wrapWidth')
  };
};

/**
 * @param {Element} svg
 * @returns {CompositionMetadata}
 */
export const parseCompositionMetadata = (svg) => {
  if (!svg?.getAttribute) fail('A root SVG element is required.', 'INVALID_SVG');
  const schema = svg.getAttribute(COMPOSITION_SCHEMA_ATTRIBUTE);
  const raw = svg.getAttribute(COMPOSITION_METADATA_ATTRIBUTE);
  if (schema === null && raw === null) {
    fail(
      'This SVG has no gbdraw composition metadata. Regenerate it, or import it through the legacy session adapter.',
      'MISSING_COMPOSITION_METADATA'
    );
  }
  if (schema === null || raw === null) {
    fail('The SVG has incomplete gbdraw composition metadata. Regenerate the diagram.');
  }
  if (String(schema) !== String(COMPOSITION_SCHEMA_VERSION)) {
    fail(`Unsupported gbdraw composition schema ${JSON.stringify(schema)}. Regenerate the diagram.`);
  }

  /** @type {unknown} */
  let source;
  try {
    source = JSON.parse(raw);
  } catch (error) {
    fail(`The SVG composition metadata is not valid JSON: ${error.message}`);
  }
  if (!isPlainObject(source)) fail('The SVG composition metadata must be an object.');
  if (!LEGEND_SIDES.has(source.legendSide)) fail('composition.legendSide is invalid.');
  if (!TITLE_SIDES.has(source.titleSide)) fail('composition.titleSide is invalid.');
  if (!Array.isArray(source.overlayObstacles)) {
    fail('composition.overlayObstacles must be an array.');
  }
  const legend = validateOptionalTarget(source.legend, 'legend');
  /** @type {CompositionMetadata} */
  const metadata = {
    legend,
    legendReflow: validateLegendReflow(source.legendReflow, {
      required: legend !== null
    }),
    legendSide: /** @type {string} */ (source.legendSide), // a member of LEGEND_SIDES
    overlayObstacles: source.overlayObstacles.map((bounds, index) => (
      validateBounds(bounds, `composition.overlayObstacles[${index}]`)
    )),
    overlayPolicy: validateOverlayPolicy(source.overlayPolicy),
    primary: /** @type {CompositionPrimaryTarget} */ (validateTarget(source.primary, 'primary', 'finalBounds')), // carries finalBounds
    spacing: validateSpacing(source.spacing),
    title: validateOptionalTarget(source.title, 'title'),
    titleSide: /** @type {string} */ (source.titleSide) // a member of TITLE_SIDES
  };
  if (source.legacyNormalized === true) metadata.legacyNormalized = true;
  return metadata;
};

/**
 * Gives a loaded Result's stored composition metadata a legend side of the
 * Result's mode: a side outside the mode's positions reads as the mode's
 * default legend position, listed in `load.repairs` at `load.path` for the
 * Load dialog, and case and outer spaces read as the position (OV-672,
 * BA-D26). Metadata that is not a JSON object stays for
 * `parseCompositionMetadata` to refuse. Returns whether it rewrote the side.
 * @param {Element} svg The root SVG element of a loaded Result.
 * @param {'circular' | 'linear'} mode The Result's mode.
 * @param {{ repairs: LoadRepair[] | null, path: string }} [load]
 * @returns {boolean}
 */
export const normalizeStoredLegendSide = (svg, mode, load = { repairs: null, path: '' }) => {
  const raw = svg.getAttribute(COMPOSITION_METADATA_ATTRIBUTE);
  /** @type {unknown} */
  let source = null;
  try {
    source = raw === null ? null : JSON.parse(raw);
  } catch {
    return false;
  }
  if (!isPlainObject(source)) return false;
  const { value: legendSide, repairs } = readLegendPositionOfMode(source.legendSide, mode, load.path, 'legend-side-defaulted');
  if (source.legendSide === legendSide) return false;
  load.repairs?.push(...repairs);
  svg.setAttribute(COMPOSITION_METADATA_ATTRIBUTE, JSON.stringify({ ...source, legendSide }));
  return true;
};

const STORED_METADATA_ATTRIBUTE = new RegExp(`\\s${COMPOSITION_METADATA_ATTRIBUTE}=(["'])(.*?)\\1`, 's');
/** @type {Record<string, string>} */
const XML_ENTITIES = { quot: '"', apos: "'", amp: '&', lt: '<', gt: '>' };
/** @param {string} text An attribute value as written in XML. */
const decodeAttributeText = (text) => text.replace(
  /&(?:#x([0-9a-f]+)|#([0-9]+)|(quot|apos|amp|lt|gt));/gi,
  (_entity, hex, decimal, name) => (hex !== undefined
    ? String.fromCodePoint(parseInt(hex, 16))
    : decimal !== undefined ? String.fromCodePoint(Number(decimal)) : XML_ENTITIES[name.toLowerCase()])
);

/**
 * Whether the root tag of a Result's text stores composition metadata whose
 * legend side is, as written, one of `mode`'s positions: a Result for which it
 * is false and that has metadata is parsed at Load for
 * `normalizeStoredLegendSide`.
 * @param {string} rootTag The `<svg ...>` tag of a Result's text.
 * @param {'circular' | 'linear'} mode The Result's mode.
 * @returns {boolean}
 */
export const storedLegendSideIsOfMode = (rootTag, mode) => {
  const value = STORED_METADATA_ATTRIBUTE.exec(rootTag)?.[2];
  if (value === undefined) return false;
  try {
    /** @type {unknown} */
    const source = JSON.parse(decodeAttributeText(value));
    return isPlainObject(source) && typeof source.legendSide === 'string'
      && legendPositionOfMode(source.legendSide, mode) === source.legendSide;
  } catch {
    return false;
  }
};

/**
 * @param {Element} svg
 * @param {string} selector
 * @returns {SVGGraphicsElement[]}
 */
const targetsFor = (svg, selector) => (
  // The role selectors match the composition <g> groups.
  /** @type {SVGGraphicsElement[]} */ (Array.from(svg.querySelectorAll?.(selector) || []))
);

/**
 * @param {Element} svg
 * @param {CompositionMetadata} metadata
 * @param {CompositionRole} role
 * @param {{ multiple?: boolean }} [options]
 * @returns {CompositionRoleBinding}
 */
const bindRole = (svg, metadata, role, { multiple = false } = {}) => {
  const targetMetadata = metadata[role];
  const targets = targetsFor(svg, ROLE_SELECTORS[role]);
  if (role === 'primary') {
    if (targets.length === 0) fail('Composition metadata does not bind any primary SVG group.');
  } else if (targetMetadata === null) {
    if (targets.length !== 0) fail(`Composition metadata has a stale ${role} role marker.`);
    return { metadata: null, targets: [] };
  } else if (targets.length !== 1) {
    fail(`Composition metadata must bind exactly one ${role} SVG group.`);
  }
  if (!multiple && targets.length > 1) fail(`Composition metadata binds duplicate ${role} groups.`);
  return { metadata: targetMetadata, targets };
};

/**
 * @param {Element} svg
 * @returns {CompositionBinding}
 */
export const bindCompositionMetadata = (svg) => {
  const metadata = parseCompositionMetadata(svg);
  return {
    metadata,
    primary: bindRole(svg, metadata, 'primary', { multiple: true }),
    legend: bindRole(svg, metadata, 'legend'),
    title: bindRole(svg, metadata, 'title')
  };
};

/**
 * @param {CompositionBounds} bounds
 * @param {number} dx
 * @param {number} dy
 * @returns {CompositionBounds}
 */
const translated = (bounds, dx, dy) => ({
  x: bounds.x + dx,
  y: bounds.y + dy,
  width: bounds.width,
  height: bounds.height
});
/** @param {CompositionBounds} bounds */
const maxX = (bounds) => bounds.x + bounds.width;
/** @param {CompositionBounds} bounds */
const maxY = (bounds) => bounds.y + bounds.height;
/**
 * @param {CompositionBounds[]} boundsList
 * @returns {CompositionBounds | null}
 */
const unionBounds = (boundsList) => {
  if (!boundsList.length) return null;
  const minX = Math.min(...boundsList.map((bounds) => bounds.x));
  const minY = Math.min(...boundsList.map((bounds) => bounds.y));
  const right = Math.max(...boundsList.map(maxX));
  const bottom = Math.max(...boundsList.map(maxY));
  return { x: minX, y: minY, width: right - minX, height: bottom - minY };
};

// The metadata and the editor hold bounds as {x, y, width, height}; Python's
// planner holds min/max boxes. The planner is the Python-order port in
// services/legend-layout.js (`planLegendComposition`), so a replan equals
// `plan_composition` given the same boxes; the payloads go back as
// `_bounds_payload` writes them.
/** @param {{ x: number, y: number, width: number, height: number }} bounds */
const toBox = (bounds) => ({
  minX: bounds.x, minY: bounds.y, maxX: bounds.x + bounds.width, maxY: bounds.y + bounds.height
});
/** @param {import('../../services/legend-layout.js').LayoutBox} box */
const toPayloadBounds = (box) => ({
  x: box.minX, y: box.minY, width: box.maxX - box.minX, height: box.maxY - box.minY
});

/**
 * Replan the composition. A `legendBox` (min/max, as the Legend layout port
 * returns it) takes precedence over `legendBounds`, and a `primaryBox` over
 * `primaryBounds`, so the planner reads those bounds exactly as Python holds
 * them.
 * @param {{
 *   primaryBounds: unknown, legendBounds?: unknown, titleBounds?: unknown,
 *   legendBox?: import('../../services/legend-layout.js').LayoutBox | null,
 *   primaryBox?: import('../../services/legend-layout.js').LayoutBox | null,
 *   legendSide?: string, titleSide?: string, spacing: unknown, overlayPolicy: unknown, overlayObstacles?: unknown[]
 * }} request
 */
export const planComposition = ({
  primaryBounds,
  primaryBox = null,
  legendBounds = null,
  titleBounds = null,
  legendBox = null,
  legendSide = 'none',
  titleSide = 'none',
  spacing,
  overlayPolicy,
  overlayObstacles = []
}) => {
  if (!LEGEND_SIDES.has(legendSide)) fail(`Unknown legend side ${JSON.stringify(legendSide)}.`);
  if (!TITLE_SIDES.has(titleSide)) fail(`Unknown title side ${JSON.stringify(titleSide)}.`);
  if (primaryBox) validateBounds(toPayloadBounds(primaryBox), 'primaryBox', { positive: true });
  const primary = primaryBox || toBox(validateBounds(primaryBounds, 'primaryBounds', { positive: true }));
  const resolvedSpacing = validateSpacing(spacing);
  const resolvedOverlayPolicy = validateOverlayPolicy(overlayPolicy);
  const obstacles = overlayObstacles.map((bounds, index) => validateBounds(bounds, `overlayObstacles[${index}]`));
  /** @type {import('../../services/legend-layout.js').LayoutBox | null} */
  let legend = null;
  if (legendSide !== 'none' && legendBox) {
    validateBounds(toPayloadBounds(legendBox), 'legendBox');
    legend = legendBox;
  } else if (legendSide !== 'none' && legendBounds) {
    legend = toBox(validateBounds(legendBounds, 'legendBounds'));
  }
  const title = titleBounds && titleSide !== 'none' ? validateBounds(titleBounds, 'titleBounds') : null;
  const plan = planLegendComposition({
    primary,
    legend,
    title: title ? toBox(title) : null,
    legendSide,
    titleSide,
    overlayObstacles: obstacles.map(toBox),
    spacing: resolvedSpacing,
    overlayPolicy: resolvedOverlayPolicy
  });
  const placements = Object.fromEntries(plan.placements.map((placement) => [
    placement.role,
    { automaticTranslation: placement.translation, finalBounds: toPayloadBounds(placement.finalBounds) }
  ]));
  return {
    width: plan.canvas.maxX - plan.canvas.minX,
    height: plan.canvas.maxY - plan.canvas.minY,
    placements,
    overlayObstacles: plan.overlayObstacles.map(toPayloadBounds),
    overlayPolicy: resolvedOverlayPolicy,
    spacing: resolvedSpacing
  };
};

/** @param {CompositionMetadata} metadata */
const localPrimaryBounds = (metadata) => translated(
  metadata.primary.finalBounds,
  -metadata.primary.automaticTranslation[0],
  -metadata.primary.automaticTranslation[1]
);

// The primary box Python planned with (`legendReflow.primaryLocalBounds`,
// recorded since Z1) while it is still the primary's box: rebuilt from x and
// width, a max edge can differ from Python's by one ulp. Null for a Result
// written before Python recorded it, or whose primary has changed since.
/**
 * @param {CompositionMetadata} metadata
 * @returns {LayoutBox | null}
 */
const recordedPrimaryBox = (metadata) => {
  const recorded = metadata.legendReflow?.primaryLocalBounds;
  if (!recorded) return null;
  const derived = toBox(localPrimaryBounds(metadata));
  const keys = /** @type {const} */ (['minX', 'minY', 'maxX', 'maxY']);
  return keys.every((key) => Math.abs(recorded[key] - derived[key]) <= 1e-6) ? recorded : null;
};

/**
 * @param {CompositionMetadata} metadata
 * @param {{
 *   legendSide?: string, titleSide?: string,
 *   legendLocalBounds?: CompositionBounds | null, legendLocalBox?: LayoutBox | null,
 *   titleLocalBounds?: CompositionBounds | null
 * }} [options]
 */
export const replanCompositionMetadata = (
  metadata,
  {
    legendSide = metadata.legendSide,
    titleSide = metadata.titleSide,
    legendLocalBounds = metadata.legend?.localBounds || null,
    legendLocalBox = null,
    titleLocalBounds = metadata.title?.localBounds || null
  } = {}
) => {
  if (!LEGEND_SIDES.has(legendSide)) fail(`Unknown legend side ${JSON.stringify(legendSide)}.`);
  if (!TITLE_SIDES.has(titleSide)) fail(`Unknown title side ${JSON.stringify(titleSide)}.`);
  const primaryAuto = metadata.primary.automaticTranslation;
  const sourceObstacles = metadata.overlayObstacles.map((bounds) => (
    translated(bounds, -primaryAuto[0], -primaryAuto[1])
  ));
  return planComposition({
    primaryBounds: localPrimaryBounds(metadata),
    primaryBox: recordedPrimaryBox(metadata),
    legendBounds: legendLocalBounds,
    legendBox: legendLocalBox,
    titleBounds: titleLocalBounds,
    legendSide,
    titleSide,
    spacing: metadata.spacing,
    overlayPolicy: metadata.overlayPolicy,
    overlayObstacles: sourceObstacles
  });
};

/**
 * @param {Element} target
 * @param {string} role
 */
const targetLeadingTranslate = (target, role) => {
  const transform = target.getAttribute?.('transform') || '';
  const leading = readLeadingTranslate(transform);
  if (!leading.found) {
    fail(`The ${role} composition target has no leading automatic translate. Regenerate the diagram.`);
  }
  return leading;
};

/**
 * @param {CompositionBinding} binding
 * @returns {{ primary: number[][], legend: number[] | null, title: number[] | null }}
 */
const userDeltas = (binding) => ({
  primary: binding.primary.targets.map((target) => {
    const current = targetLeadingTranslate(target, 'primary');
    return [
      current.x - binding.metadata.primary.automaticTranslation[0],
      current.y - binding.metadata.primary.automaticTranslation[1]
    ];
  }),
  legend: binding.legend.metadata
    ? (() => {
        const current = targetLeadingTranslate(binding.legend.targets[0], 'legend');
        return [
          current.x - binding.legend.metadata.automaticTranslation[0],
          current.y - binding.legend.metadata.automaticTranslation[1]
        ];
      })()
    : null,
  title: binding.title.metadata
    ? (() => {
        const current = targetLeadingTranslate(binding.title.targets[0], 'title');
        return [
          current.x - binding.title.metadata.automaticTranslation[0],
          current.y - binding.title.metadata.automaticTranslation[1]
        ];
      })()
    : null
});

/**
 * @param {CompositionMatrix} matrix
 * @param {number} x
 * @param {number} y
 */
const transformPoint = (matrix, x, y) => ({
  x: matrix.a * x + matrix.c * y + matrix.e,
  y: matrix.b * x + matrix.d * y + matrix.f
});

/**
 * @param {SVGGraphicsElement} target
 * @param {CompositionMatrix} matrix
 * @returns {CompositionMatrix}
 */
const rootRelativeMatrix = (target, matrix) => {
  const rootMatrix = target.ownerSVGElement?.getCTM?.();
  if (!rootMatrix || !(/** @type {const} */ (['a', 'b', 'c', 'd', 'e', 'f'])).every(
    (key) => Number.isFinite(rootMatrix[key])
  )) return matrix;
  const determinant = rootMatrix.a * rootMatrix.d - rootMatrix.b * rootMatrix.c;
  if (!Number.isFinite(determinant) || Math.abs(determinant) < Number.EPSILON) return matrix;
  const inverse = {
    a: rootMatrix.d / determinant,
    b: -rootMatrix.b / determinant,
    c: -rootMatrix.c / determinant,
    d: rootMatrix.a / determinant,
    e: (rootMatrix.c * rootMatrix.f - rootMatrix.d * rootMatrix.e) / determinant,
    f: (rootMatrix.b * rootMatrix.e - rootMatrix.a * rootMatrix.f) / determinant
  };
  return {
    a: inverse.a * matrix.a + inverse.c * matrix.b,
    b: inverse.b * matrix.a + inverse.d * matrix.b,
    c: inverse.a * matrix.c + inverse.c * matrix.d,
    d: inverse.b * matrix.c + inverse.d * matrix.d,
    e: inverse.a * matrix.e + inverse.c * matrix.f + inverse.e,
    f: inverse.b * matrix.e + inverse.d * matrix.f + inverse.f
  };
};

/**
 * @param {SVGGraphicsElement} target
 * @returns {CompositionBounds}
 */
const targetFinalBounds = (target) => {
  const bbox = target.getBBox();
  const sourceMatrix = target.getCTM?.();
  const matrix = sourceMatrix ? rootRelativeMatrix(target, sourceMatrix) : null;
  if (matrix && (/** @type {const} */ (['a', 'b', 'c', 'd', 'e', 'f'])).every((key) => Number.isFinite(matrix[key]))) {
    const corners = [
      transformPoint(matrix, bbox.x, bbox.y),
      transformPoint(matrix, bbox.x + bbox.width, bbox.y),
      transformPoint(matrix, bbox.x, bbox.y + bbox.height),
      transformPoint(matrix, bbox.x + bbox.width, bbox.y + bbox.height)
    ];
    const x = Math.min(...corners.map((point) => point.x));
    const y = Math.min(...corners.map((point) => point.y));
    const right = Math.max(...corners.map((point) => point.x));
    const bottom = Math.max(...corners.map((point) => point.y));
    return { x, y, width: right - x, height: bottom - y };
  }
  const leading = readLeadingTranslate(target.getAttribute?.('transform') || '');
  return translated(
    { x: bbox.x, y: bbox.y, width: bbox.width, height: bbox.height },
    leading.x,
    leading.y
  );
};

/**
 * @param {SVGGraphicsElement} target
 * @returns {CompositionBounds}
 */
export const measureCompositionTargetLocalBounds = (target) => {
  const leading = targetLeadingTranslate(target, target.getAttribute?.(COMPOSITION_ROLE_ATTRIBUTE) || 'unknown');
  return translated(targetFinalBounds(target), -leading.x, -leading.y);
};

/**
 * @param {Element} target
 * @param {number[]} automaticTranslation
 * @param {number[]} delta
 */
const setLeading = (target, automaticTranslation, delta) => {
  target.setAttribute(
    'transform',
    replaceLeadingTranslate(
      target.getAttribute('transform'),
      automaticTranslation[0] + delta[0],
      automaticTranslation[1] + delta[1]
    )
  );
};

/**
 * @param {CompositionRole} role
 * @param {number[]} automaticTranslation
 * @param {'finalBounds' | 'localBounds'} boundsKey
 * @param {CompositionBounds} bounds
 * @returns {CompositionTargetBase & { finalBounds?: CompositionBounds, localBounds?: CompositionBounds }}
 */
const targetPayload = (role, automaticTranslation, boundsKey, bounds) => ({
  automaticTranslation,
  [boundsKey]: bounds,
  role,
  selector: ROLE_SELECTORS[role]
});

/** @param {number} value */
const wireNumber = (value) => Object.is(value, -0) ? 0 : value;
/** @param {number[]} value */
const wireTranslation = (value) => value.map(wireNumber);
/** @param {CompositionBounds} bounds */
const wireBounds = (bounds) => ({
  height: wireNumber(bounds.height),
  width: wireNumber(bounds.width),
  x: wireNumber(bounds.x),
  y: wireNumber(bounds.y)
});

const CANVAS_PADDING_SIDES = Object.freeze(['top', 'right', 'bottom', 'left']);

/**
 * Normalize editor canvas padding; invalid or negative sides are 0.
 * @param {Partial<Record<string, unknown>> | null | undefined} value
 * @returns {Record<string, number>}
 */
export const normalizeCanvasPadding = (value) => Object.fromEntries(
  CANVAS_PADDING_SIDES.map((side) => {
    const number = Number(value?.[side]);
    return [side, Number.isFinite(number) && number >= 0 ? number : 0];
  })
);

/** @param {Partial<Record<string, unknown>> | null | undefined} value */
export const hasCanvasPadding = (value) => (
  Object.values(normalizeCanvasPadding(value)).some((number) => number !== 0)
);

/**
 * @param {Element} svg
 * @param {string} name
 * @param {string} value
 */
const setRootAttribute = (svg, name, value) => {
  if (svg.getAttribute(name) === value) return false;
  svg.setAttribute(name, value);
  return true;
};

/**
 * The bounds, widened by the edge padding, of each Legend or plot title the
 * user moved: a moved item keeps its offset and the canvas grows to hold it
 * (Owner, 2026-10-23). An unmoved or hidden item adds nothing. The bounds are
 * Python's (`localBounds`) at the item's current translate, so a detached
 * candidate needs no browser measurement.
 * @param {CompositionBinding | null} binding Null for a root without composition metadata.
 * @returns {CompositionBounds[]}
 */
const movedDecorationBounds = (binding) => {
  if (!binding) return [];
  const { metadata } = binding;
  const edge = metadata.spacing.edgePaddingPx;
  /** @type {Array<{ role: 'legend' | 'title', side: string, target: CompositionLocalTarget | null, element: SVGGraphicsElement | undefined }>} */
  const items = [
    { role: 'legend', side: metadata.legendSide, target: metadata.legend, element: binding.legend.targets[0] },
    { role: 'title', side: metadata.titleSide, target: metadata.title, element: binding.title.targets[0] }
  ];
  return items.flatMap(({ role, side, target, element }) => {
    if (!target || !element || side === 'none') return [];
    const current = targetLeadingTranslate(element, role);
    if (current.x === target.automaticTranslation[0] && current.y === target.automaticTranslation[1]) return [];
    const bounds = target.localBounds;
    return [{
      x: bounds.x + current.x - edge,
      y: bounds.y + current.y - edge,
      width: bounds.width + (2 * edge),
      height: bounds.height + (2 * edge)
    }];
  });
};

/**
 * Size the canvas of one SVG root: the automatic canvas, grown so that each
 * moved Legend or plot title stays inside it with the edge padding (Owner,
 * 2026-10-23; never smaller than the automatic canvas), then the editor
 * canvas padding around it. The automatic viewBox and size are recorded once
 * on the root, so applying the same state again changes nothing, and a Result
 * that neither grows nor is padded stays byte-identical. Returns whether the
 * root changed.
 * @param {Element} svg
 * @param {Partial<Record<string, unknown>> | null | undefined} value The canvas padding.
 * @param {CompositionBinding | null} [binding] The root's binding when the caller holds it.
 * @returns {boolean}
 */
export const applyCanvasToSvg = (svg, value, binding = null) => {
  const padding = normalizeCanvasPadding(value);
  const padded = Object.values(padding).some((number) => number !== 0);
  const moved = movedDecorationBounds(binding
    || (svg.getAttribute(COMPOSITION_SCHEMA_ATTRIBUTE) === String(COMPOSITION_SCHEMA_VERSION) ? bindCompositionMetadata(svg) : null));
  const recordedBase = svg.getAttribute('data-original-view-box');
  if (!recordedBase && !padded && !moved.length) return false;
  // A missing viewBox gives NaN below and throws, so `base` is a string past that check.
  const base = /** @type {string} */ (recordedBase || svg.getAttribute('viewBox'));
  const [x, y, width, height] = String(base || '').trim().split(/[\s,]+/).map(Number);
  if (![x, y, width, height].every(Number.isFinite) || width <= 0 || height <= 0) {
    throw new Error('The SVG base viewBox is invalid. Regenerate the diagram.');
  }
  const grows = moved.some((bounds) => (
    bounds.x < x || bounds.y < y || maxX(bounds) > x + width || maxY(bounds) > y + height
  ));
  if (!recordedBase && !padded && !grows) return false;
  const baseWidth = svg.getAttribute('data-original-width')
    || String(Number.parseFloat(svg.getAttribute('width') ?? '') || width);
  const baseHeight = svg.getAttribute('data-original-height')
    || String(Number.parseFloat(svg.getAttribute('height') ?? '') || height);
  let changed = false;
  if (!recordedBase) {
    changed = setRootAttribute(svg, 'data-original-view-box', base) || changed;
    changed = setRootAttribute(svg, 'data-original-width', baseWidth) || changed;
    changed = setRootAttribute(svg, 'data-original-height', baseHeight) || changed;
  }
  if (!padded && !grows) {
    changed = setRootAttribute(svg, 'viewBox', base) || changed;
    changed = setRootAttribute(svg, 'width', `${baseWidth}px`) || changed;
    return setRootAttribute(svg, 'height', `${baseHeight}px`) || changed;
  }
  const automatic = { x, y, width, height };
  // `grows` holds only when a moved item exists, so the union has it.
  const canvas = grows ? /** @type {CompositionBounds} */ (unionBounds([automatic, ...moved])) : automatic;
  const outerWidth = canvas.width + padding.left + padding.right;
  const outerHeight = canvas.height + padding.top + padding.bottom;
  /** @param {string} size @param {number} extent @param {number} outerExtent */
  const scaled = (size, extent, outerExtent) => {
    const number = Number.parseFloat(size);
    return number === extent ? outerExtent : number * (outerExtent / extent);
  };
  changed = setRootAttribute(svg, 'viewBox', [
    canvas.x - padding.left, canvas.y - padding.top, outerWidth, outerHeight
  ].map(wireNumber).join(' ')) || changed;
  changed = setRootAttribute(svg, 'width', `${wireNumber(scaled(baseWidth, width, outerWidth))}px`) || changed;
  return setRootAttribute(svg, 'height', `${wireNumber(scaled(baseHeight, height, outerHeight))}px`) || changed;
};

/**
 * @param {Element} svg
 * @param {CompositionEditOptions} [options]
 */
export const applyCompositionEdit = (svg, options = {}) => {
  const binding = bindCompositionMetadata(svg);
  const metadata = binding.metadata;
  const deltas = userDeltas(binding);
  const legendSide = options.legendSide ?? metadata.legendSide;
  const titleSide = options.titleSide ?? metadata.titleSide;
  const legendTarget = binding.legend.targets[0] || null;
  const titleTarget = binding.title.targets[0] || null;
  // The Legend and title bounds are Python's (the metadata), or the Legend
  // layout port's after a Legend edit (`legendLocalBox`), never the browser's
  // text extents, so a composition edit docks them where Python would. Only a
  // legacy Session's normalized metadata is measured on the mounted SVG.
  /** @param {SVGGraphicsElement | null} target @param {string} side @param {CompositionBounds | null | undefined} recorded */
  const localBoundsOf = (target, side, recorded) => (
    metadata.legacyNormalized && target && side !== 'none'
      ? measureCompositionTargetLocalBounds(target)
      : recorded || null
  );
  const legendLocalBox = options.legendLocalBox ?? null;
  const legendLocalBounds = options.legendLocalBounds
    ?? (legendLocalBox ? toPayloadBounds(legendLocalBox) : localBoundsOf(legendTarget, legendSide, metadata.legend?.localBounds));
  const titleLocalBounds = options.titleLocalBounds
    ?? localBoundsOf(titleTarget, titleSide, metadata.title?.localBounds);
  // Legacy sessions are admitted while detached, where browser geometry APIs
  // return empty boxes. Their renderer-authored canvas geometry is sufficient
  // for a safe, layout-neutral import. The first explicit layout edit happens
  // on the already-mounted SVG, so replace that conservative primary envelope
  // with the authoritative painted bounds as part of the same user action.
  const planningMetadata = metadata.legacyNormalized
    ? {
        ...metadata,
        primary: {
          ...metadata.primary,
          // bindRole rejects a Result without primary groups, so the union exists.
          finalBounds: wireBounds(/** @type {CompositionBounds} */ (unionBounds(
            binding.primary.targets.map(measureCompositionTargetLocalBounds)
          )))
        }
      }
    : metadata;
  const plan = replanCompositionMetadata(planningMetadata, {
    legendSide,
    titleSide,
    legendLocalBounds,
    legendLocalBox,
    titleLocalBounds
  });
  const primaryPlacement = plan.placements.primary;
  binding.primary.targets.forEach((target, index) => {
    setLeading(target, primaryPlacement.automaticTranslation, deltas.primary[index]);
  });

  /** @type {ReturnType<typeof targetPayload> | null} */
  let legendPayload = null;
  if (legendTarget && metadata.legend) {
    const placement = plan.placements.legend;
    const automatic = placement
      ? placement.automaticTranslation
      : metadata.legend.automaticTranslation;
    if (placement) setLeading(legendTarget, automatic, deltas.legend || [0, 0]);
    if (legendSide === 'none') legendTarget.setAttribute('display', 'none');
    else legendTarget.removeAttribute('display');
    legendPayload = targetPayload(
      'legend',
      wireTranslation(automatic),
      'localBounds',
      wireBounds(/** @type {CompositionBounds} */ (legendLocalBounds)) // set for a bound Legend
    );
  }

  /** @type {ReturnType<typeof targetPayload> | null} */
  let titlePayload = null;
  if (titleTarget && metadata.title) {
    const placement = plan.placements.title;
    const automatic = placement
      ? placement.automaticTranslation
      : metadata.title.automaticTranslation;
    if (placement) setLeading(titleTarget, automatic, deltas.title || [0, 0]);
    if (titleSide === 'none') titleTarget.setAttribute('display', 'none');
    else titleTarget.removeAttribute('display');
    titlePayload = targetPayload(
      'title',
      wireTranslation(automatic),
      'localBounds',
      wireBounds(/** @type {CompositionBounds} */ (titleLocalBounds)) // set for a bound title
    );
  }

  const nextMetadata = {
    legend: legendPayload,
    legendReflow: metadata.legendReflow,
    legendSide,
    overlayObstacles: plan.overlayObstacles.map(wireBounds),
    overlayPolicy: plan.overlayPolicy,
    primary: targetPayload(
      'primary',
      wireTranslation(primaryPlacement.automaticTranslation),
      'finalBounds',
      wireBounds(primaryPlacement.finalBounds)
    ),
    spacing: { ...plan.spacing },
    title: titlePayload,
    titleSide
  };
  const baseWidth = wireNumber(plan.width);
  const baseHeight = wireNumber(plan.height);
  svg.setAttribute('data-original-view-box', `0 0 ${baseWidth} ${baseHeight}`);
  svg.setAttribute('data-original-width', String(baseWidth));
  svg.setAttribute('data-original-height', String(baseHeight));
  svg.setAttribute(COMPOSITION_SCHEMA_ATTRIBUTE, String(COMPOSITION_SCHEMA_VERSION));
  svg.setAttribute(COMPOSITION_METADATA_ATTRIBUTE, JSON.stringify(nextMetadata));
  applyCanvasToSvg(svg, options.canvasPadding);
  return bindCompositionMetadata(svg);
};

/** @param {Element} svg */
export const resetCompositionUserDeltas = (svg) => {
  const binding = bindCompositionMetadata(svg);
  binding.primary.targets.forEach((target) => {
    setLeading(target, binding.metadata.primary.automaticTranslation, [0, 0]);
  });
  if (binding.legend.metadata) {
    setLeading(
      binding.legend.targets[0],
      binding.legend.metadata.automaticTranslation,
      [0, 0]
    );
  }
  if (binding.title.metadata) {
    setLeading(binding.title.targets[0], binding.title.metadata.automaticTranslation, [0, 0]);
  }
  return binding;
};

/**
 * @param {Element} svg
 * @param {{ primary?: unknown[], legend?: unknown, title?: unknown }} [deltas]
 */
export const applyCompositionUserDeltas = (svg, deltas = {}) => {
  const binding = bindCompositionMetadata(svg);
  let changed = false;

  /** @param {Element | undefined} target @param {number[]} automaticTranslation @param {unknown} requested */
  const applyDelta = (target, automaticTranslation, requested) => {
    if (!target || !Array.isArray(requested) || requested.length !== 2) return;
    const x = Number(requested[0]);
    const y = Number(requested[1]);
    if (!Number.isFinite(x) || !Number.isFinite(y)) return;
    const current = targetLeadingTranslate(target, target.getAttribute?.(COMPOSITION_ROLE_ATTRIBUTE) || 'unknown');
    const nextX = automaticTranslation[0] + x;
    const nextY = automaticTranslation[1] + y;
    if (current.x === nextX && current.y === nextY) return;
    setLeading(target, automaticTranslation, [x, y]);
    changed = true;
  };

  binding.primary.targets.forEach((target, index) => {
    applyDelta(target, binding.metadata.primary.automaticTranslation, deltas.primary?.[index]);
  });
  if (binding.legend.metadata) {
    applyDelta(
      binding.legend.targets[0],
      binding.legend.metadata.automaticTranslation,
      deltas.legend
    );
  }
  if (binding.title.metadata) {
    applyDelta(
      binding.title.targets[0],
      binding.title.metadata.automaticTranslation,
      deltas.title
    );
  }

  return { binding, changed };
};

/** @param {Element | null | undefined} element */
const elementTag = (element) => String(element?.tagName || '').toLowerCase().replace(/^.*:/, '');
/**
 * @param {SVGSVGElement} svg
 * @returns {LegacyTargets}
 */
const discoverLegacyTargets = (svg) => {
  const legend = svg.getElementById?.('legend') || null;
  const title = svg.getElementById?.(PLOT_TITLE_ID) || null;
  const primary = Array.from(svg.children || []).filter((element) => (
    elementTag(element) === 'g' && element !== legend && element !== title
  ));
  return { primary, legend, title };
};

/**
 * @param {unknown} value
 * @returns {CompositionDelta}
 */
const legacyDelta = (value) => {
  if (!Array.isArray(value) || value.length !== 2) return [0, 0];
  const x = Number(value[0]);
  const y = Number(value[1]);
  return Number.isFinite(x) && Number.isFinite(y) ? [x, y] : [0, 0];
};

/**
 * @param {Element | null} target
 * @returns {CompositionLegendReflow | null}
 */
const legacyLegendReflow = (target) => {
  if (!target) return null;
  let colorRectSize = 14;
  const firstColorRect = target.querySelector?.(
    'path[fill]:not([fill="none"]):not([fill^="url("])'
  );
  const pathData = firstColorRect?.getAttribute?.('d') || '';
  const sizeMatch = pathData.match(/L\s+([\d.]+),/);
  if (sizeMatch) {
    const measured = Number(sizeMatch[1]);
    if (Number.isFinite(measured) && measured > 0) colorRectSize = measured;
  }
  return {
    colorRectSize,
    lineHeight: LEGACY_LEGEND_LINE_HEIGHT_RATIO_V0 * colorRectSize,
    textXOffset: LEGACY_LEGEND_TEXT_OFFSET_RATIO_V0 * colorRectSize
  };
};

/**
 * @param {Element} target
 * @param {unknown} deltaValue
 * @returns {CompositionDelta}
 */
const factorLegacyUserDelta = (target, deltaValue) => {
  const [dx, dy] = legacyDelta(deltaValue);
  const source = target.getAttribute?.('transform') || '';
  if (dx === 0 && dy === 0) {
    target.setAttribute('transform', prependTranslate(source, 0, 0));
    return [0, 0];
  }
  const leading = readLeadingTranslate(source);
  const baseline = leading.found
    ? replaceLeadingTranslate(source, leading.x - dx, leading.y - dy)
    : prependTranslate(source, -dx, -dy);
  target.setAttribute('transform', prependTranslate(baseline, dx, dy));
  return [dx, dy];
};

/**
 * @param {unknown} value
 * @param {string} label
 * @returns {CompositionBounds}
 */
const legacyViewBox = (value, label) => {
  const parts = String(value || '')
    .trim()
    .split(/[\s,]+/)
    .map(Number);
  if (
    parts.length !== 4
    || !parts.every(Number.isFinite)
    || parts[2] <= 0
    || parts[3] <= 0
  ) {
    fail(`The legacy SVG has no valid ${label} geometry. Regenerate the diagram.`);
  }
  return { x: parts[0], y: parts[1], width: parts[2], height: parts[3] };
};

/**
 * @param {Element | null} target
 * @param {unknown} deltaValue
 * @returns {{ x: number, y: number }}
 */
const legacyBaseTranslation = (target, deltaValue) => {
  const [dx, dy] = legacyDelta(deltaValue);
  const leading = readLeadingTranslate(target?.getAttribute?.('transform') || '');
  return leading.found
    ? { x: leading.x - dx, y: leading.y - dy }
    : { x: -dx, y: -dy };
};

/**
 * @param {CompositionBounds} rootBounds
 * @param {Element | null} target
 * @param {string} side
 * @param {unknown} deltaValue
 * @returns {CompositionBounds | null}
 */
const legacyTitleBounds = (rootBounds, target, side, deltaValue) => {
  if (!target) return null;
  const translation = legacyBaseTranslation(target, deltaValue);
  const edgeMargin = 24;
  let height = 0;
  if (side === 'bottom') {
    height = 2 * (rootBounds.y + rootBounds.height - edgeMargin - translation.y);
  } else if (side === 'top') {
    height = 2 * (translation.y - rootBounds.y - edgeMargin);
  }
  if (!Number.isFinite(height) || height <= 0 || height > rootBounds.height) {
    const fontSize = Number(target.querySelector?.('text')?.getAttribute?.('font-size'));
    height = Number.isFinite(fontSize) && fontSize > 0
      ? Math.min(rootBounds.height, fontSize * 1.25)
      : rootBounds.height;
  }
  return {
    x: rootBounds.x,
    y: translation.y - height / 2,
    width: rootBounds.width,
    height
  };
};

/**
 * @param {SVGSVGElement} svg
 * @param {LegacyTargets} targets
 * @param {{ legendSide: string, titleSide: string, userDeltas: LegacyUserDeltas | null }} options
 */
const legacyAdmissionBounds = (
  svg,
  targets,
  { legendSide, titleSide, userDeltas: savedUserDeltas }
) => {
  const rootBounds = legacyViewBox(svg.getAttribute('viewBox'), 'root viewBox');
  const horizontalValue = svg.getAttribute('data-horizontal-viewbox');
  const verticalValue = svg.getAttribute('data-vertical-viewbox');
  const horizontal = horizontalValue === null
    ? null
    : legacyViewBox(horizontalValue, 'horizontal viewBox');
  const vertical = verticalValue === null
    ? null
    : legacyViewBox(verticalValue, 'vertical viewBox');

  let primary = rootBounds;
  let legend = targets.legend ? rootBounds : null;
  if (horizontal && vertical) {
    // Pre-composition linear SVGs persist two renderer-authored canvases.
    // H excludes vertical-legend width; V excludes horizontal-legend height.
    // Their H.width x V.height envelope is therefore the orientation-neutral
    // canvas and needs no browser layout measurement.
    primary = {
      x: horizontal.x,
      y: vertical.y,
      width: horizontal.width,
      height: vertical.height
    };
    if (targets.legend && (legendSide === 'top' || legendSide === 'bottom')) {
      const translation = legacyBaseTranslation(targets.legend, savedUserDeltas?.legend);
      const width = 2 * (horizontal.x + horizontal.width / 2 - translation.x);
      const height = horizontal.height - vertical.height;
      if (Number.isFinite(width) && width > 0 && Number.isFinite(height) && height > 0) {
        legend = { x: translation.x, y: translation.y, width, height };
      }
    }
  }

  return {
    primary,
    legend,
    title: legacyTitleBounds(rootBounds, targets.title, titleSide, savedUserDeltas?.title)
  };
};

/**
 * @param {SVGSVGElement} svg The root SVG element of a legacy Result; checked below.
 * @param {{ legendSide?: string, titleSide?: string, userDeltas?: LegacyUserDeltas | null }} [options]
 */
export const normalizeLegacyComposition = (
  svg,
  { legendSide = 'none', titleSide = 'none', userDeltas: savedUserDeltas = null } = {}
) => {
  if (!svg?.getAttribute) fail('A root SVG element is required.', 'INVALID_SVG');
  if (svg.getAttribute(COMPOSITION_SCHEMA_ATTRIBUTE) !== null || svg.getAttribute(COMPOSITION_METADATA_ATTRIBUTE) !== null) {
    return bindCompositionMetadata(svg);
  }
  if (!LEGEND_SIDES.has(legendSide)) fail(`Unknown legacy legend side ${JSON.stringify(legendSide)}.`);
  if (!TITLE_SIDES.has(titleSide)) fail(`Unknown legacy title side ${JSON.stringify(titleSide)}.`);
  const targets = discoverLegacyTargets(svg);
  if (targets.primary.length === 0) fail('The legacy SVG has no top-level primary groups.');

  const admissionBounds = legacyAdmissionBounds(svg, targets, {
    legendSide,
    titleSide,
    userDeltas: savedUserDeltas
  });

  const primaryDeltaValues = Array.isArray(savedUserDeltas?.primary?.[0])
    ? savedUserDeltas.primary
    : targets.primary.map((target) => (
        target.getAttribute?.('id') === LENGTH_BAR_ID
          ? savedUserDeltas?.lengthBar ?? savedUserDeltas?.primary
          : savedUserDeltas?.primary
      ));
  targets.primary.forEach((target, index) => {
    factorLegacyUserDelta(target, primaryDeltaValues[index] || [0, 0]);
    target.setAttribute(COMPOSITION_ROLE_ATTRIBUTE, 'primary');
  });
  /** @param {Element | null} target @param {'legend' | 'title'} role */
  const decorateLegacyTarget = (target, role) => {
    if (!target) return null;
    factorLegacyUserDelta(target, savedUserDeltas?.[role]);
    target.setAttribute(COMPOSITION_ROLE_ATTRIBUTE, role);
    // The admission bounds of a role exist exactly when its target does.
    const localBounds = /** @type {CompositionBounds} */ (admissionBounds[role]);
    return targetPayload(role, [0, 0], 'localBounds', wireBounds(localBounds));
  };
  const metadata = {
    legend: decorateLegacyTarget(targets.legend, 'legend'),
    legendReflow: legacyLegendReflow(targets.legend),
    legendSide,
    overlayObstacles: [],
    overlayPolicy: {
      candidateScoreOrder: [...LEGACY_OVERLAY_POLICY_V0.candidateScoreOrder],
      canvasGrowthCandidateOrder: [...LEGACY_OVERLAY_POLICY_V0.canvasGrowthCandidateOrder],
      canvasGrowthScoreOrder: [...LEGACY_OVERLAY_POLICY_V0.canvasGrowthScoreOrder],
      quadrantBoundaryRatio: LEGACY_OVERLAY_POLICY_V0.quadrantBoundaryRatio
    },
    primary: targetPayload('primary', [0, 0], 'finalBounds', wireBounds(admissionBounds.primary)),
    spacing: { ...LEGACY_COMPOSITION_SPACING_V0 },
    title: decorateLegacyTarget(targets.title, 'title'),
    titleSide,
    legacyNormalized: true
  };
  svg.setAttribute(COMPOSITION_SCHEMA_ATTRIBUTE, String(COMPOSITION_SCHEMA_VERSION));
  svg.setAttribute(COMPOSITION_METADATA_ATTRIBUTE, JSON.stringify(metadata));
  return bindCompositionMetadata(svg);
};

/** @param {Element} svg */
export const compositionUserDeltas = (svg) => userDeltas(bindCompositionMetadata(svg));

/**
 * @param {Element} target
 * @returns {Set<string>}
 */
const stableFeatureIds = (target) => new Set(
  Array.from(target.querySelectorAll?.('[data-gbdraw-stable-feature-id]') || [])
    .map((feature) => String(
      feature.getAttribute?.('data-gbdraw-stable-feature-id') || ''
    ).trim())
    .filter(Boolean)
);

/**
 * @param {Element} target
 * @param {LegacyRecordDecision[]} decisions
 */
const legacyRecordKey = (target, decisions) => {
  const available = stableFeatureIds(target);
  if (available.size === 0) return '';
  const matches = decisions.filter((decision) => {
    const anchor = decision?.anchor;
    return anchor && [anchor.stableFeatureSvgId, anchor.biologicalFeatureId]
      .map((value) => String(value || '').trim())
      .filter(Boolean)
      .some((value) => available.has(value));
  });
  if (matches.length > 1) {
    fail('Legacy composition metadata binds one record group to multiple alignment records.');
  }
  return matches.length === 1 ? String(matches[0].recordKey || '') : '';
};

/**
 * @param {Element} target
 * @param {{ x: number, y: number }} fallbackValue
 */
const legacyRenderedTranslation = (target, fallbackValue) => {
  // Before record-keyed translation metadata, ordinary Linear output stored
  // the final X translation as the record group's second translate. Y was the
  // renderer's row position, so the keyed Y baseline remains authoritative.
  // Grid output already had data-record-key; it is intentionally rejected by
  // this legacy-only path because its second translate also contains grid X.
  const outer = readLeadingTranslate(target.getAttribute?.('transform') || '');
  const inner = outer.found ? readLeadingTranslate(outer.tail) : { found: false };
  if (!inner.found) {
    fail('Legacy record composition has no materializable record translation.');
  }
  // `inner` is the found result of readLeadingTranslate here.
  return { x: /** @type {{ x: number }} */ (inner).x, y: fallbackValue.y };
};

// This is the only bridge from SVG composition placement to the record-keyed
// translation owner. Current SVGs bind directly by recordKey. Released
// pre-recordKey SVGs bind once through the active plan's stable feature IDs;
// DOM order, row number, rendered record ID, and array position never become
// persisted identity.
/**
 * @param {Element} svg
 * @param {({ recordKey?: unknown, x?: unknown, y?: unknown } | null | undefined)[]} [baseTranslations]
 * @param {unknown[]} [recordKeys]
 * @param {{ records?: unknown } | null} [activePlan] The active plan; only its `records` are read.
 */
export const materializeRecordTranslations = (
  svg,
  baseTranslations = [],
  recordKeys = [],
  activePlan = null
) => {
  const binding = bindCompositionMetadata(svg);
  const deltas = userDeltas(binding).primary;
  const fallback = new Map(
    (Array.isArray(baseTranslations) ? baseTranslations : []).map((entry) => [
      String(entry?.recordKey || ''),
      { x: Number(entry?.x), y: Number(entry?.y) }
    ])
  );
  const decisions = Array.isArray(activePlan?.records) ? activePlan.records : [];
  const decisionByRecord = new Map(
    decisions.map((decision) => [String(decision?.recordKey || ''), decision])
  );
  const byRecord = new Map();
  binding.primary.targets.forEach((target, index) => {
    const explicitRecordKey = String(
      target.getAttribute?.('data-record-key') || ''
    ).trim();
    const recordKey = explicitRecordKey || legacyRecordKey(target, decisions);
    if (!recordKey) return;
    if (byRecord.has(recordKey)) {
      fail(`Composition metadata binds duplicate record key ${JSON.stringify(recordKey)}.`);
    }
    const fallbackValue = fallback.get(recordKey) || { x: 0, y: 0 };
    const rawRenderedX = target.getAttribute?.('data-record-translation-x');
    const rawRenderedY = target.getAttribute?.('data-record-translation-y');
    const renderedX = rawRenderedX === null ? Number.NaN : Number(rawRenderedX);
    const renderedY = rawRenderedY === null ? Number.NaN : Number(rawRenderedY);
    if ((rawRenderedX === null) !== (rawRenderedY === null)) {
      fail(`Composition metadata has an incomplete translation for record ${JSON.stringify(recordKey)}.`);
    }
    let rendered = { x: fallbackValue.x, y: fallbackValue.y };
    if (Number.isFinite(renderedX) && Number.isFinite(renderedY)) {
      rendered = { x: renderedX, y: renderedY };
    } else if (rawRenderedX !== null || rawRenderedY !== null) {
      fail(`Composition metadata has a non-finite translation for record ${JSON.stringify(recordKey)}.`);
    } else if (!explicitRecordKey && activePlan) {
      rendered = legacyRenderedTranslation(target, fallbackValue);
    } else if (explicitRecordKey && activePlan) {
      fail(`Composition metadata cannot materialize record ${JSON.stringify(recordKey)}. Regenerate the diagram.`);
    }
    const [deltaX, deltaY] = deltas[index] || [0, 0];
    byRecord.set(recordKey, {
      recordKey,
      x: rendered.x + deltaX,
      y: rendered.y + deltaY
    });
  });
  const orderedKeys = Array.isArray(recordKeys) && recordKeys.length > 0
    ? recordKeys.map((key) => String(key || ''))
    : [...byRecord.keys()];
  return orderedKeys.map((recordKey) => {
    const materialized = byRecord.get(recordKey) || fallback.get(recordKey);
    const decision = decisionByRecord.get(recordKey);
    if (
      activePlan
      && decision?.status !== 'skipped'
      && !byRecord.has(recordKey)
    ) {
      fail(`Composition metadata cannot bind alignment record ${JSON.stringify(recordKey)}.`);
    }
    if (!materialized || !Number.isFinite(materialized.x) || !Number.isFinite(materialized.y)) {
      fail(`Composition metadata has no translation for record ${JSON.stringify(recordKey)}.`);
    }
    return {
      recordKey,
      x: Number(materialized.x),
      y: Number(materialized.y)
    };
  });
};
