// @ts-check
/** @import { DrawingState } from '../state.js' */
import { ruleLegendCaptions } from '../services/specific-color-rules.js';
import {
  appliedFeatureColors,
  estimateColorFactor,
  interpolateColor,
  resolveCollinearMatchColor,
  resolvePairwiseLegendGradientColorKeys
} from '../utils/color-utils.js';
import { PAIRWISE_LEGEND_SELECTOR } from '../services/legend-svg.js';
import { getGroupsByBaseIds } from '../services/svg-result-normalization.js';
import { resolveTrackSlotSkewColorValue } from './track-slot-colors.js';

export { getGroupsByBaseIds } from '../services/svg-result-normalization.js';

const normalizeComparableColor = (value) => String(value || '').trim().toLowerCase();

const setColorAttributeIfChanged = (element, attribute, value) => {
  if (normalizeComparableColor(element.getAttribute(attribute)) === normalizeComparableColor(value)) {
    return false;
  }
  element.setAttribute(attribute, value);
  return true;
};

// Whether two palettes give every key the same color (the palette watcher's test).
export const paletteColorsEqual = (left, right) => {
  const keys = new Set([...Object.keys(left || {}), ...Object.keys(right || {})]);
  return Array.from(keys).every(
    (key) => normalizeComparableColor(left?.[key]) === normalizeComparableColor(right?.[key])
  );
};

const paletteColorKeysEqual = (left, right, keys) => keys.every(
  (key) => normalizeComparableColor(left?.[key]) === normalizeComparableColor(right?.[key])
);

// The palette color of a Legend row, by its key (caption).
/** @type {Record<string, string>} */
const keyToColorKey = {
  CDS: 'CDS',
  'D-loop': 'D-loop',
  repeat_region: 'repeat_region',
  tmRNA: 'tmRNA',
  tRNA: 'tRNA',
  rRNA: 'rRNA',
  ncRNA: 'ncRNA',
  misc_feature: 'misc_feature',
  mobile_element: 'mobile_element',
  'GC content': 'gc_content',
  'GC skew (+)': 'skew_high',
  'GC skew (-)': 'skew_low'
};
/**
 * @param {string} legendKey
 * @param {Record<string, string>} palette
 * @returns {string | null}
 */
const resolveOtherLegendColor = (legendKey, palette) => {
  if (!legendKey) return null;
  const lowerKey = legendKey.toLowerCase();
  if (lowerKey === 'other proteins') return palette.CDS || null;
  if (!lowerKey.startsWith('other ')) return null;
  let raw = legendKey.slice(6).trim();
  if (!raw) return null;
  if (raw.toLowerCase() === 'proteins') return palette.CDS || null;
  if (raw.endsWith('s')) raw = raw.slice(0, -1);
  return palette[raw] || null;
};
/**
 * @param {string} legendKey
 * @param {Record<string, string>} palette
 * @returns {string | null}
 */
const paletteLegendColor = (legendKey, palette) => {
  if (!legendKey) return null;
  const colorKey = keyToColorKey[legendKey];
  if (colorKey && palette[colorKey]) return palette[colorKey];
  if (palette[legendKey]) return palette[legendKey];
  return resolveOtherLegendColor(legendKey, palette);
};

/**
 * @typedef {object} SvgStylesOptions
 * @property {Record<string, any>} state App state (state.js; not yet typed).
 * @property {(source: any, callback: (...args: any[]) => void, options?: Record<string, any>) => any} watch Vue `watch`
 * @property {(callback?: () => void) => Promise<void>} nextTick Vue `nextTick`
 * @property {((reason: string) => boolean) | null} [commitActiveResultEdit]
 *   The preview owner's commit of an edit to the displayed Result (R1, R13).
 * @property {(options?: { recolor?: Record<string, any>, prepareRules?: boolean }) => boolean | Promise<boolean>} projectPaletteAndRules
 *   The root's projection of the palette and the rules (R3), which prepares the
 *   rule matches, paints the tracks through this owner, and shows the fills.
 */

/** @param {SvgStylesOptions} options */
export const createSvgStyles = ({
  state,
  watch,
  nextTick,
  // R13: the preview owner's commit of an edit to the displayed Result.
  commitActiveResultEdit = null,
  // R13: the composition root's projection of the palette and the rules (R3),
  // which prepares the rule matches and paints the tracks through this owner.
  projectPaletteAndRules
}) => {
  const {
    svgContent,
    appliedPaletteColors,
    pairwiseMatchFactors,
    svgContainer,
    mode
  } = state;

  const updatePairwiseLegendGradientStops = (pairwiseLegend, colors) => {
    let updated = false;
    pairwiseLegend.querySelectorAll('linearGradient').forEach((gradient) => {
      const legendKey = gradient.closest('g[data-legend-key]')?.getAttribute('data-legend-key') || '';
      const { minKey, maxKey } = resolvePairwiseLegendGradientColorKeys(legendKey);
      const minColor = colors[minKey];
      const maxColor = colors[maxKey];
      if (!minColor || !maxColor) return;
      const stops = gradient.querySelectorAll('stop');
      if (stops.length >= 2) {
        updated = setColorAttributeIfChanged(stops[0], 'stop-color', minColor) || updated;
        updated = setColorAttributeIfChanged(stops[1], 'stop-color', maxColor) || updated;
      }
    });
    return updated;
  };

  // A Legend row with a Legend color, or one a Specific color rule draws (also
  // a rule captioned like a palette key), keeps its color; the palette colors
  // the others.
  /** @param {DrawingState} drawing */
  const legendRowKeepsColor = (drawing) => {
    const legendCaption = ruleLegendCaptions({
      rules: drawing.manualSpecificRules, legendEntries: drawing.legendEntries?.value || [],
      originalLegendOrder: state.originalLegendOrder?.value || []
    });
    const ruleRows = new Set(drawing.manualSpecificRules.map(legendCaption));
    /** @param {string} caption */
    return (caption) => Boolean(drawing.legendColorOverrides[caption]) || ruleRows.has(caption);
  };
  // OV-276: the color the palette gives each Legend panel row of the active
  // drawing, as its repaint colors the swatch. The panel reads it; nothing
  // writes it into the rows, which are inputs of the rule preparation.
  const paletteLegendRowColors = () => {
    const drawing = state.activeDrawing();
    const keepsColor = legendRowKeepsColor(drawing);
    const colors = appliedFeatureColors(state);
    /** @type {Map<string, string>} */
    const rowColors = new Map();
    (drawing.legendEntries?.value || []).forEach((entry) => {
      const caption = entry?.caption;
      const color = caption && !keepsColor(caption) ? paletteLegendColor(caption, colors) : null;
      if (color) rowColors.set(caption, color);
    });
    return rowColors;
  };

  // The palette styles the tracks of the displayed Result with the settings of
  // the drawing of that Result's mode: its skew slot colors (OV-108). Feature
  // and Legend row fills are editor operations (`projectPaletteAndRules`).
  const applyPaletteToSvg = ({
    recolorPairwise = false,
    recolorCollinear = false
  } = {}) => {
    const resultMode = (state.generatedMode?.value ?? mode.value) === 'linear' ? 'linear' : 'circular';
    const drawing = state.drawings[resultMode];
    const busy = state.sessionOperationAvailability?.();
    if (busy) return busy;
    if (!svgContent.value) return;
    if (!svgContainer.value) return;

    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return;

    const colors = appliedFeatureColors(state);
    let updatedCount = 0;

    const gcContentGroups = getGroupsByBaseIds(
      svg,
      ['gc_content'],
      ['dinucleotide_content']
    );
    if (gcContentGroups.length > 0 && colors.gc_content) {
      gcContentGroups.forEach((gcContentGroup) => {
        const gcPaths = gcContentGroup.querySelectorAll('path');
        gcPaths.forEach((path) => {
          if (setColorAttributeIfChanged(path, 'fill', colors.gc_content)) updatedCount++;
        });
      });
    }

    const skewGroups = getGroupsByBaseIds(
      svg,
      ['skew', 'gc_skew'],
      ['dinucleotide_skew']
    );
    if (skewGroups.length > 0) {
      skewGroups.forEach((skewGroup) => {
        const slotId = String(skewGroup.getAttribute('data-gbdraw-slot-id') || '').trim();
        const customSlotsEnabled = resultMode === 'circular'
          ? drawing.adv.circular_track_slots_enabled
          : drawing.adv.linear_track_slots_enabled;
        const slots = resultMode === 'circular'
          ? drawing.adv.circular_track_slots
          : drawing.adv.linear_track_slots;
        const slot = customSlotsEnabled && Array.isArray(slots)
          ? slots.find((candidate) => (
              candidate?.enabled !== false &&
              candidate?.renderer === 'dinucleotide_skew' &&
              String(candidate?.id || '').trim() === slotId
            ))
          : null;
        const positiveColor = resolveTrackSlotSkewColorValue(/** @type {any} */ ({
          slot,
          key: 'positive_color',
          currentColors: colors
        }));
        const negativeColor = resolveTrackSlotSkewColorValue(/** @type {any} */ ({
          slot,
          key: 'negative_color',
          currentColors: colors
        }));
        const skewPaths = skewGroup.querySelectorAll('path');
        let pathIndex = 0;
        skewPaths.forEach((path) => {
          const fill = path.getAttribute('fill');
          if (fill && fill !== 'white' && fill !== 'none') {
            if (pathIndex === 0) {
              if (setColorAttributeIfChanged(path, 'fill', positiveColor)) updatedCount++;
            } else if (pathIndex === 1) {
              if (setColorAttributeIfChanged(path, 'fill', negativeColor)) updatedCount++;
            }
            pathIndex++;
          }
        });
      });
    }

    if (
      (recolorPairwise || recolorCollinear)
      && colors.pairwise_match_min
      && colors.pairwise_match_max
    ) {
      const committedPairwiseFactors = pairwiseMatchFactors.value || {};
      let nextPairwiseFactors = committedPairwiseFactors;
      const retainPairwiseFactor = (pathKey, factor) => {
        if (committedPairwiseFactors[pathKey] === factor) return;
        if (nextPairwiseFactors === committedPairwiseFactors) {
          nextPairwiseFactors = { ...committedPairwiseFactors };
        }
        nextPairwiseFactors[pathKey] = factor;
      };
      let compIdx = 1;
      let compGroup = svg.getElementById(`comparison${compIdx}`);
      while (compGroup) {
        const matchPaths = compGroup.querySelectorAll('path');
        matchPaths.forEach((path, pathIdx) => {
          const pathKey = `comp${compIdx}_path${pathIdx}`;
          const currentFill = path.getAttribute('fill');
          if (currentFill) {
            const collinearityBlockId = path.getAttribute('data-collinearity-block-id') || '';
            const collinearityColorMode = path.getAttribute('data-collinearity-color-mode') || '';
            const metadataText = path.getAttribute('data-identity-factor');
            const metadataFactor = metadataText === null || metadataText === '' ? NaN : Number(metadataText);
            const collinearColor = resolveCollinearMatchColor({
              blockId: collinearityBlockId,
              colorMode: collinearityColorMode,
              orientation: path.getAttribute('data-collinearity-orientation') || '',
              identityFactor: Number.isFinite(metadataFactor) ? metadataFactor : null,
              colors
            });
            if (collinearColor) {
              if (
                recolorCollinear
                && setColorAttributeIfChanged(path, 'fill', collinearColor)
              ) {
                updatedCount++;
              }
              return;
            }
            if (collinearityBlockId && !collinearityColorMode) return;
            if (!recolorPairwise) return;

            let factor;
            if (Number.isFinite(metadataFactor)) {
              factor = metadataFactor;
              retainPairwiseFactor(pathKey, factor);
            } else if (committedPairwiseFactors[pathKey] !== undefined) {
              factor = committedPairwiseFactors[pathKey];
            } else {
              const origMin = window._origPairwiseMin || '#FFE7E7';
              const origMax = window._origPairwiseMax || '#FF7272';
              factor = estimateColorFactor(currentFill, origMin, origMax);
              retainPairwiseFactor(pathKey, factor);
            }
            const newColor = interpolateColor(colors.pairwise_match_min, colors.pairwise_match_max, factor);
            if (setColorAttributeIfChanged(path, 'fill', newColor)) updatedCount++;
          }
        });
        compIdx++;
        compGroup = svg.getElementById(`comparison${compIdx}`);
      }
      if (nextPairwiseFactors !== committedPairwiseFactors) {
        pairwiseMatchFactors.value = nextPairwiseFactors;
      }
    }

    if (colors.pairwise_match_min && colors.pairwise_match_max) {
      const allPairwiseLegends = svg.querySelectorAll(PAIRWISE_LEGEND_SELECTOR);
      allPairwiseLegends.forEach((pairwiseLegend) => {
        if (updatePairwiseLegendGradientStops(pairwiseLegend, colors)) updatedCount++;
      });
    }

    if (updatedCount > 0) commitActiveResultEdit?.('palette');
  };

  const applyTrackVisibility = () => {
    const drawing = state.activeDrawing();
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    if (!svgContent.value) return;
    if (!svgContainer.value) return;

    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return;

    let updated = false;

    const gcContentGroups = getGroupsByBaseIds(
      svg,
      ['gc_content'],
      ['dinucleotide_content']
    );
    if (gcContentGroups.length > 0) {
      const shouldHide = mode.value === 'circular' ? drawing.form.suppress_gc : !drawing.form.show_gc;
      gcContentGroups.forEach((gcContentGroup) => {
        const currentDisplay = gcContentGroup.getAttribute('display');
        if (shouldHide && currentDisplay !== 'none') {
          gcContentGroup.setAttribute('display', 'none');
          updated = true;
        } else if (!shouldHide && currentDisplay === 'none') {
          gcContentGroup.removeAttribute('display');
          updated = true;
        }
      });
    }

    const skewGroups = getGroupsByBaseIds(
      svg,
      ['skew', 'gc_skew'],
      ['dinucleotide_skew']
    );
    if (skewGroups.length > 0) {
      const shouldHide = mode.value === 'circular' ? drawing.form.suppress_skew : !drawing.form.show_skew;
      skewGroups.forEach((skewGroup) => {
        const currentDisplay = skewGroup.getAttribute('display');
        if (shouldHide && currentDisplay !== 'none') {
          skewGroup.setAttribute('display', 'none');
          updated = true;
        } else if (!shouldHide && currentDisplay === 'none') {
          skewGroup.removeAttribute('display');
          updated = true;
        }
      });
    }

    const depthGroups = getGroupsByBaseIds(svg, ['depth'], ['depth']);
    if (depthGroups.length > 0) {
      const shouldHide = !drawing.form.show_depth;
      depthGroups.forEach((depthGroup) => {
        const currentDisplay = depthGroup.getAttribute('display');
        if (shouldHide && currentDisplay !== 'none') {
          depthGroup.setAttribute('display', 'none');
          updated = true;
        } else if (!shouldHide && currentDisplay === 'none') {
          depthGroup.removeAttribute('display');
          updated = true;
        }
      });
    }

    if (updated) {
      commitActiveResultEdit?.('track-visibility');
      console.log('Track visibility updated');
    }
  };

  watch(
    appliedPaletteColors,
    (colors, previousColors) => {
      if (state.semanticFileWatchersSuppressed?.value || paletteColorsEqual(colors, previousColors)) return;
      // Inferred comparison factors are lossy, so only re-interpolate a family
      // when one of its palette endpoints actually changed.
      const recolorPairwise = !paletteColorKeysEqual(
        colors,
        previousColors,
        ['pairwise_match_min', 'pairwise_match_max']
      );
      const recolorCollinear = !paletteColorKeysEqual(
        colors,
        previousColors,
        [
          'collinear_block_plus_min',
          'collinear_block_plus',
          'collinear_block_minus_min',
          'collinear_block_minus'
        ]
      );
      nextTick(() => projectPaletteAndRules({ recolor: { recolorPairwise, recolorCollinear } }));
    },
    { deep: true }
  );

  // Each drawing's track toggles; a mode switch changes neither drawing, so it
  // applies nothing (R10).
  /** @type {DrawingState[]} */ ([state.drawings.circular, state.drawings.linear]).forEach((drawing) => watch(
    () => [drawing.form.suppress_gc, drawing.form.suppress_skew, drawing.form.show_gc, drawing.form.show_skew, drawing.form.show_depth],
    () => {
      if (drawing !== state.activeDrawing()) return;
      if (state.semanticFileWatchersSuppressed?.value || state.sessionOperationAvailability?.()) return;
      applyTrackVisibility();
    }
  ));

  return {
    applyPaletteToSvg,
    paletteLegendRowColors,
    applyTrackVisibility
  };
};
