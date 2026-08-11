import {
  prependTranslate,
  readLeadingTranslate,
  replaceLeadingTranslate
} from './transform-utils.js';

export const COMPOSITION_SCHEMA_VERSION = 1;
export const COMPOSITION_SCHEMA_ATTRIBUTE = 'data-gbdraw-composition-schema';
export const COMPOSITION_METADATA_ATTRIBUTE = 'data-gbdraw-composition';
export const COMPOSITION_ROLE_ATTRIBUTE = 'data-gbdraw-composition-role';

const ROLE_SELECTORS = Object.freeze({
  primary: '[data-gbdraw-composition-role="primary"]',
  legend: '[data-gbdraw-composition-role="legend"]',
  title: '[data-gbdraw-composition-role="title"]'
});
const LEGEND_SIDES = new Set([
  'left', 'right', 'top', 'bottom',
  'upper_left', 'upper_right', 'lower_left', 'lower_right', 'none'
]);
const TITLE_SIDES = new Set(['top', 'bottom', 'center', 'none']);
const OVERLAY_CANDIDATE_SCORE_NAMES = new Set([
  'totalAnchorDistance',
  'xAnchorDistance',
  'yAnchorDistance',
  'nearEdgeX',
  'nearEdgeY'
]);
const OVERLAY_CANVAS_GROWTH_CANDIDATE_NAMES = new Set(['horizontal', 'vertical']);
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
  constructor(message, code = 'INVALID_COMPOSITION_METADATA') {
    super(message);
    this.name = 'CompositionMetadataError';
    this.code = code;
  }
}

const fail = (message, code) => {
  throw new CompositionMetadataError(message, code);
};

const isPlainObject = (value) => (
  Boolean(value) && typeof value === 'object' && !Array.isArray(value)
);

const finiteNumber = (value, path) => {
  if (typeof value !== 'number' || !Number.isFinite(value)) {
    fail(`${path} must be a finite number.`);
  }
  return value;
};

const nonNegativeNumber = (value, path) => {
  const number = finiteNumber(value, path);
  if (number < 0) fail(`${path} must be non-negative.`);
  return number;
};

const positiveNumber = (value, path) => {
  const number = finiteNumber(value, path);
  if (number <= 0) fail(`${path} must be positive.`);
  return number;
};

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

const validateTranslation = (value, path) => {
  if (!Array.isArray(value) || value.length !== 2) {
    fail(`${path} must contain exactly two numbers.`);
  }
  return [finiteNumber(value[0], `${path}[0]`), finiteNumber(value[1], `${path}[1]`)];
};

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

const validateOptionalTarget = (value, role) => (
  value === null ? null : validateTarget(value, role, 'localBounds')
);

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
  return {
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
};

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
  const metadata = {
    legend,
    legendReflow: validateLegendReflow(source.legendReflow, {
      required: legend !== null
    }),
    legendSide: source.legendSide,
    overlayObstacles: source.overlayObstacles.map((bounds, index) => (
      validateBounds(bounds, `composition.overlayObstacles[${index}]`)
    )),
    overlayPolicy: validateOverlayPolicy(source.overlayPolicy),
    primary: validateTarget(source.primary, 'primary', 'finalBounds'),
    spacing: validateSpacing(source.spacing),
    title: validateOptionalTarget(source.title, 'title'),
    titleSide: source.titleSide
  };
  if (source.legacyNormalized === true) metadata.legacyNormalized = true;
  return metadata;
};

const targetsFor = (svg, selector) => Array.from(svg.querySelectorAll?.(selector) || []);

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

export const bindCompositionMetadata = (svg) => {
  const metadata = parseCompositionMetadata(svg);
  return {
    metadata,
    primary: bindRole(svg, metadata, 'primary', { multiple: true }),
    legend: bindRole(svg, metadata, 'legend'),
    title: bindRole(svg, metadata, 'title')
  };
};

const translated = (bounds, dx, dy) => ({
  x: bounds.x + dx,
  y: bounds.y + dy,
  width: bounds.width,
  height: bounds.height
});
const maxX = (bounds) => bounds.x + bounds.width;
const maxY = (bounds) => bounds.y + bounds.height;
const centerX = (bounds) => bounds.x + bounds.width / 2;
const centerY = (bounds) => bounds.y + bounds.height / 2;
const unionBounds = (boundsList) => {
  if (!boundsList.length) return null;
  const minX = Math.min(...boundsList.map((bounds) => bounds.x));
  const minY = Math.min(...boundsList.map((bounds) => bounds.y));
  const right = Math.max(...boundsList.map(maxX));
  const bottom = Math.max(...boundsList.map(maxY));
  return { x: minX, y: minY, width: right - minX, height: bottom - minY };
};
const intersects = (left, right, clearance = 0) => !(
  maxX(left) + clearance <= right.x ||
  maxX(right) + clearance <= left.x ||
  maxY(left) + clearance <= right.y ||
  maxY(right) + clearance <= left.y
);
const alignMin = (bounds, x, y) => ({
  translation: [x - bounds.x, y - bounds.y],
  bounds: { x, y, width: bounds.width, height: bounds.height }
});

const dockLegend = (primary, legend, side, spacing) => {
  if (side === 'left') {
    return alignMin(
      legend,
      primary.x - spacing.dockGapPx - legend.width,
      centerY(primary) - legend.height / 2
    );
  }
  if (side === 'right') {
    return alignMin(legend, maxX(primary) + spacing.dockGapPx, centerY(primary) - legend.height / 2);
  }
  if (side === 'top') {
    return alignMin(
      legend,
      centerX(primary) - legend.width / 2,
      primary.y - spacing.dockGapPx - legend.height
    );
  }
  return alignMin(legend, centerX(primary) - legend.width / 2, maxY(primary) + spacing.dockGapPx);
};

const overlayAxisRange = (minimum, maximum, itemSize, nearMinimum, boundaryRatio) => {
  let low = minimum;
  let high = maximum - itemSize;
  const midpointStart = low + boundaryRatio * (high - low);
  if (nearMinimum) high = Math.min(high, midpointStart);
  else low = Math.max(low, midpointStart);
  return low > high ? null : [low, high];
};

const overlayAxisValues = (range, itemSize, obstacles, clearance, axis) => {
  const values = new Set(range);
  obstacles.forEach((obstacle) => {
    if (axis === 'x') {
      values.add(obstacle.x - clearance - itemSize);
      values.add(maxX(obstacle) + clearance);
    } else {
      values.add(obstacle.y - clearance - itemSize);
      values.add(maxY(obstacle) + clearance);
    }
  });
  return [...values].filter((value) => value >= range[0] && value <= range[1]).sort((a, b) => a - b);
};

const compareScore = (left, right) => {
  for (let index = 0; index < left.length; index += 1) {
    if (left[index] !== right[index]) return left[index] - right[index];
  }
  return 0;
};

const overlayLegend = (primary, legend, side, obstacles, spacing, policy) => {
  const left = side === 'upper_left' || side === 'lower_left';
  const upper = side === 'upper_left' || side === 'upper_right';
  const anchorX = left ? primary.x : maxX(primary) - legend.width;
  const anchorY = upper ? primary.y : maxY(primary) - legend.height;
  const anchor = alignMin(legend, anchorX, anchorY);
  const conflicts = (bounds) => obstacles
    .map((obstacle, index) => intersects(bounds, obstacle, spacing.overlayClearancePx) ? index : -1)
    .filter((index) => index >= 0);
  const initialConflicts = conflicts(anchor.bounds);
  const xRange = overlayAxisRange(
    primary.x,
    maxX(primary),
    legend.width,
    left,
    policy.quadrantBoundaryRatio
  );
  const yRange = overlayAxisRange(
    primary.y,
    maxY(primary),
    legend.height,
    upper,
    policy.quadrantBoundaryRatio
  );
  if (xRange && yRange) {
    if (initialConflicts.length === 0) return anchor;
    const candidates = [];
    overlayAxisValues(xRange, legend.width, obstacles, spacing.overlayClearancePx, 'x').forEach((x) => {
      overlayAxisValues(yRange, legend.height, obstacles, spacing.overlayClearancePx, 'y').forEach((y) => {
        candidates.push([x, y]);
      });
    });
    const candidateScore = ([x, y]) => {
      const metrics = {
        totalAnchorDistance: Math.abs(x - anchorX) + Math.abs(y - anchorY),
        xAnchorDistance: Math.abs(x - anchorX),
        yAnchorDistance: Math.abs(y - anchorY),
        nearEdgeX: left ? x : -x,
        nearEdgeY: upper ? y : -y
      };
      return policy.candidateScoreOrder.map((name) => metrics[name]);
    };
    candidates.sort((a, b) => compareScore(candidateScore(a), candidateScore(b)));
    for (const [x, y] of candidates) {
      const candidate = alignMin(legend, x, y);
      if (conflicts(candidate.bounds).length === 0) return candidate;
    }
  }
  const horizontal = alignMin(
    legend,
    left ? primary.x - spacing.overlayClearancePx - legend.width : maxX(primary) + spacing.overlayClearancePx,
    anchorY
  );
  const vertical = alignMin(
    legend,
    anchorX,
    upper ? primary.y - spacing.overlayClearancePx - legend.height : maxY(primary) + spacing.overlayClearancePx
  );
  const candidatesByName = { horizontal, vertical };
  const canvasGrowthCandidates = policy.canvasGrowthCandidateOrder.map(
    (name) => candidatesByName[name]
  );
  const growthKey = (candidate, index) => {
    const union = unionBounds([primary, candidate.bounds]);
    const metrics = {
      addedArea: union.width * union.height - primary.width * primary.height,
      addedExtent: union.width - primary.width + union.height - primary.height,
      candidateOrder: index
    };
    return policy.canvasGrowthScoreOrder.map((name) => metrics[name]);
  };
  return canvasGrowthCandidates
    .map((candidate, index) => ({ candidate, score: growthKey(candidate, index) }))
    .sort((a, b) => compareScore(a.score, b.score))[0]
    .candidate;
};

const placeTitle = (primary, title, side, spacing, legendPlacement, legendSide) => {
  if (side === 'center') {
    return alignMin(title, centerX(primary) - title.width / 2, centerY(primary) - title.height / 2);
  }
  const sameSide = Boolean(legendPlacement) && side === legendSide;
  let placement;
  if (side === 'top') {
    const targetBottom = sameSide
      ? legendPlacement.bounds.y - spacing.stackGapPx
      : primary.y - spacing.titleGapPx;
    placement = alignMin(title, centerX(primary) - title.width / 2, targetBottom - title.height);
  } else {
    const targetTop = sameSide
      ? maxY(legendPlacement.bounds) + spacing.stackGapPx
      : maxY(primary) + spacing.titleGapPx;
    placement = alignMin(title, centerX(primary) - title.width / 2, targetTop);
  }
  if (
    legendPlacement && (legendSide === 'left' || legendSide === 'right') &&
    intersects(placement.bounds, legendPlacement.bounds)
  ) {
    const y = side === 'top'
      ? legendPlacement.bounds.y - spacing.stackGapPx - title.height
      : maxY(legendPlacement.bounds) + spacing.stackGapPx;
    placement = alignMin(title, centerX(primary) - title.width / 2, y);
  }
  return placement;
};

export const planComposition = ({
  primaryBounds,
  legendBounds = null,
  titleBounds = null,
  legendSide = 'none',
  titleSide = 'none',
  spacing,
  overlayPolicy,
  overlayObstacles = []
}) => {
  if (!LEGEND_SIDES.has(legendSide)) fail(`Unknown legend side ${JSON.stringify(legendSide)}.`);
  if (!TITLE_SIDES.has(titleSide)) fail(`Unknown title side ${JSON.stringify(titleSide)}.`);
  const primary = validateBounds(primaryBounds, 'primaryBounds', { positive: true });
  const resolvedSpacing = validateSpacing(spacing);
  const resolvedOverlayPolicy = validateOverlayPolicy(overlayPolicy);
  const obstacles = overlayObstacles.map((bounds, index) => validateBounds(bounds, `overlayObstacles[${index}]`));
  const working = [{ role: 'primary', translation: [0, 0], bounds: primary }];
  let legendPlacement = null;
  if (legendBounds && legendSide !== 'none') {
    const legend = validateBounds(legendBounds, 'legendBounds');
    if (legend.width > 0 && legend.height > 0) {
      legendPlacement = ['left', 'right', 'top', 'bottom'].includes(legendSide)
        ? dockLegend(primary, legend, legendSide, resolvedSpacing)
        : overlayLegend(
          primary,
          legend,
          legendSide,
          obstacles,
          resolvedSpacing,
          resolvedOverlayPolicy
        );
      working.push({ role: 'legend', ...legendPlacement });
    }
  }
  if (titleBounds && titleSide !== 'none') {
    const title = validateBounds(titleBounds, 'titleBounds');
    if (title.width > 0 && title.height > 0) {
      working.push({
        role: 'title',
        ...placeTitle(primary, title, titleSide, resolvedSpacing, legendPlacement, legendSide)
      });
    }
  }
  const painted = unionBounds(working.map((placement) => placement.bounds));
  const outerX = resolvedSpacing.edgePaddingPx - painted.x;
  const outerY = resolvedSpacing.edgePaddingPx - painted.y;
  const placements = Object.fromEntries(working.map((placement) => [
    placement.role,
    {
      automaticTranslation: [
        placement.translation[0] + outerX,
        placement.translation[1] + outerY
      ],
      finalBounds: translated(placement.bounds, outerX, outerY)
    }
  ]));
  return {
    width: painted.width + resolvedSpacing.edgePaddingPx * 2,
    height: painted.height + resolvedSpacing.edgePaddingPx * 2,
    placements,
    overlayObstacles: obstacles.map((bounds) => translated(bounds, outerX, outerY)),
    overlayPolicy: resolvedOverlayPolicy,
    spacing: resolvedSpacing
  };
};

const localPrimaryBounds = (metadata) => translated(
  metadata.primary.finalBounds,
  -metadata.primary.automaticTranslation[0],
  -metadata.primary.automaticTranslation[1]
);

export const replanCompositionMetadata = (
  metadata,
  {
    legendSide = metadata.legendSide,
    titleSide = metadata.titleSide,
    legendLocalBounds = metadata.legend?.localBounds || null,
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
    legendBounds: legendLocalBounds,
    titleBounds: titleLocalBounds,
    legendSide,
    titleSide,
    spacing: metadata.spacing,
    overlayPolicy: metadata.overlayPolicy,
    overlayObstacles: sourceObstacles
  });
};

const targetLeadingTranslate = (target, role) => {
  const transform = target.getAttribute?.('transform') || '';
  const leading = readLeadingTranslate(transform);
  if (!leading.found) {
    fail(`The ${role} composition target has no leading automatic translate. Regenerate the diagram.`);
  }
  return leading;
};

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

const transformPoint = (matrix, x, y) => ({
  x: matrix.a * x + matrix.c * y + matrix.e,
  y: matrix.b * x + matrix.d * y + matrix.f
});

const targetFinalBounds = (target) => {
  const bbox = target.getBBox();
  const matrix = target.getCTM?.();
  if (matrix && ['a', 'b', 'c', 'd', 'e', 'f'].every((key) => Number.isFinite(matrix[key]))) {
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

export const measureCompositionTargetLocalBounds = (target) => {
  const leading = targetLeadingTranslate(target, target.getAttribute?.(COMPOSITION_ROLE_ATTRIBUTE) || 'unknown');
  return translated(targetFinalBounds(target), -leading.x, -leading.y);
};

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

const targetPayload = (role, automaticTranslation, boundsKey, bounds) => ({
  automaticTranslation,
  [boundsKey]: bounds,
  role,
  selector: ROLE_SELECTORS[role]
});

const wireNumber = (value) => Object.is(value, -0) ? 0 : value;
const wireTranslation = (value) => value.map(wireNumber);
const wireBounds = (bounds) => ({
  height: wireNumber(bounds.height),
  width: wireNumber(bounds.width),
  x: wireNumber(bounds.x),
  y: wireNumber(bounds.y)
});

const editorPadding = (value) => {
  if (!isPlainObject(value)) return null;
  const padding = Object.fromEntries(['top', 'right', 'bottom', 'left'].map((side) => {
    const number = Number(value[side]);
    return [side, Number.isFinite(number) && number >= 0 ? number : 0];
  }));
  return Object.values(padding).some((number) => number !== 0) ? padding : null;
};

export const applyCompositionEdit = (svg, options = {}) => {
  const binding = bindCompositionMetadata(svg);
  const metadata = binding.metadata;
  const deltas = userDeltas(binding);
  const legendSide = options.legendSide ?? metadata.legendSide;
  const titleSide = options.titleSide ?? metadata.titleSide;
  const legendTarget = binding.legend.targets[0] || null;
  const titleTarget = binding.title.targets[0] || null;
  const legendLocalBounds = options.legendLocalBounds ?? (
    legendTarget && legendSide !== 'none'
      ? measureCompositionTargetLocalBounds(legendTarget)
      : metadata.legend?.localBounds || null
  );
  const titleLocalBounds = options.titleLocalBounds ?? (
    titleTarget && titleSide !== 'none'
      ? measureCompositionTargetLocalBounds(titleTarget)
      : metadata.title?.localBounds || null
  );
  const plan = replanCompositionMetadata(metadata, {
    legendSide,
    titleSide,
    legendLocalBounds,
    titleLocalBounds
  });
  const primaryPlacement = plan.placements.primary;
  binding.primary.targets.forEach((target, index) => {
    setLeading(target, primaryPlacement.automaticTranslation, deltas.primary[index]);
  });

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
      wireBounds(legendLocalBounds)
    );
  }

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
      wireBounds(titleLocalBounds)
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
  if (metadata.legacyNormalized) nextMetadata.legacyNormalized = true;
  const baseWidth = wireNumber(plan.width);
  const baseHeight = wireNumber(plan.height);
  if (svg.dataset) {
    svg.dataset.originalViewBox = `0 0 ${baseWidth} ${baseHeight}`;
    svg.dataset.originalWidth = String(baseWidth);
    svg.dataset.originalHeight = String(baseHeight);
  }
  const padding = editorPadding(options.canvasPadding);
  if (padding) {
    const width = baseWidth + padding.left + padding.right;
    const height = baseHeight + padding.top + padding.bottom;
    svg.setAttribute('viewBox', `${-padding.left} ${-padding.top} ${wireNumber(width)} ${wireNumber(height)}`);
    svg.setAttribute('width', `${wireNumber(width)}px`);
    svg.setAttribute('height', `${wireNumber(height)}px`);
  } else {
    svg.setAttribute('viewBox', `0 0 ${baseWidth} ${baseHeight}`);
    svg.setAttribute('width', `${baseWidth}px`);
    svg.setAttribute('height', `${baseHeight}px`);
  }
  svg.setAttribute(COMPOSITION_SCHEMA_ATTRIBUTE, String(COMPOSITION_SCHEMA_VERSION));
  svg.setAttribute(COMPOSITION_METADATA_ATTRIBUTE, JSON.stringify(nextMetadata));
  return bindCompositionMetadata(svg);
};

export const reconcileCompositionTitle = (
  svg,
  titleTarget,
  titleSide = 'none',
  { canvasPadding = null } = {}
) => {
  const metadata = parseCompositionMetadata(svg);
  if (!TITLE_SIDES.has(titleSide)) fail(`Unknown title side ${JSON.stringify(titleSide)}.`);

  targetsFor(svg, ROLE_SELECTORS.title).forEach((target) => {
    if (target !== titleTarget) target.removeAttribute(COMPOSITION_ROLE_ATTRIBUTE);
  });

  let titlePayload = null;
  if (titleTarget) {
    if (titleSide !== 'none') titleTarget.removeAttribute('display');
    let leading = readLeadingTranslate(titleTarget.getAttribute?.('transform') || '');
    if (!leading.found) {
      titleTarget.setAttribute(
        'transform',
        prependTranslate(titleTarget.getAttribute?.('transform'), 0, 0)
      );
      leading = readLeadingTranslate(titleTarget.getAttribute('transform'));
    }
    titleTarget.setAttribute(COMPOSITION_ROLE_ATTRIBUTE, 'title');
    titlePayload = metadata.title || targetPayload(
      'title',
      wireTranslation([leading.x, leading.y]),
      'localBounds',
      wireBounds(measureCompositionTargetLocalBounds(titleTarget))
    );
  }

  const nextMetadata = {
    ...metadata,
    title: titlePayload,
    titleSide: titleTarget ? titleSide : 'none'
  };
  svg.setAttribute(COMPOSITION_METADATA_ATTRIBUTE, JSON.stringify(nextMetadata));
  return applyCompositionEdit(svg, {
    titleSide: nextMetadata.titleSide,
    titleLocalBounds: titleTarget ? measureCompositionTargetLocalBounds(titleTarget) : null,
    canvasPadding
  });
};

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

const elementTag = (element) => String(element?.tagName || '').toLowerCase().replace(/^.*:/, '');
const discoverLegacyTargets = (svg) => {
  const legend = svg.getElementById?.('legend') || null;
  const title = svg.getElementById?.('plot_title') || null;
  const primary = Array.from(svg.children || []).filter((element) => (
    elementTag(element) === 'g' && element !== legend && element !== title
  ));
  return { primary, legend, title };
};

const legacyDelta = (value) => {
  if (!Array.isArray(value) || value.length !== 2) return [0, 0];
  const x = Number(value[0]);
  const y = Number(value[1]);
  return Number.isFinite(x) && Number.isFinite(y) ? [x, y] : [0, 0];
};

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

  const primaryDeltaValues = Array.isArray(savedUserDeltas?.primary?.[0])
    ? savedUserDeltas.primary
    : targets.primary.map((target) => (
        target.getAttribute?.('id') === 'length_bar'
          ? savedUserDeltas?.lengthBar ?? savedUserDeltas?.primary
          : savedUserDeltas?.primary
      ));
  targets.primary.forEach((target, index) => {
    factorLegacyUserDelta(target, primaryDeltaValues[index] || [0, 0]);
    target.setAttribute(COMPOSITION_ROLE_ATTRIBUTE, 'primary');
  });
  const primaryFinalBounds = unionBounds(
    targets.primary.map(measureCompositionTargetLocalBounds)
  );
  const decorateLegacyTarget = (target, role) => {
    if (!target) return null;
    factorLegacyUserDelta(target, savedUserDeltas?.[role]);
    target.setAttribute(COMPOSITION_ROLE_ATTRIBUTE, role);
    const localBounds = measureCompositionTargetLocalBounds(target);
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
    primary: targetPayload('primary', [0, 0], 'finalBounds', wireBounds(primaryFinalBounds)),
    spacing: { ...LEGACY_COMPOSITION_SPACING_V0 },
    title: decorateLegacyTarget(targets.title, 'title'),
    titleSide,
    legacyNormalized: true
  };
  svg.setAttribute(COMPOSITION_SCHEMA_ATTRIBUTE, String(COMPOSITION_SCHEMA_VERSION));
  svg.setAttribute(COMPOSITION_METADATA_ATTRIBUTE, JSON.stringify(metadata));
  return bindCompositionMetadata(svg);
};

export const compositionUserDeltas = (svg) => userDeltas(bindCompositionMetadata(svg));
