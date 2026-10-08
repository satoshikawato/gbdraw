// @ts-check
/** @import { DrawingState } from '../state.js' */
import { ruleMatcher } from '../services/rule-matchers.js';
import {
  estimateColorFactor,
  interpolateColor,
  resolveCollinearMatchColor,
  resolvePairwiseLegendGradientColorKeys
} from '../utils/color-utils.js';
import {
  getFeatureElementIndex,
  getFeatureFillElements,
  getFeatureIdentity
} from './feature-editor/svg-actions.js';
import { isFeatureFillTarget } from '../services/feature-dom.js';
import { getAllFeatureLegendGroups, PAIRWISE_LEGEND_SELECTOR, parseTransformXY } from '../services/legend-svg.js';
import { getFeatureOverride } from '../services/feature-override-identity.js';
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

const paletteColorsEqual = (left, right) => {
  const keys = new Set([...Object.keys(left || {}), ...Object.keys(right || {})]);
  return Array.from(keys).every(
    (key) => normalizeComparableColor(left?.[key]) === normalizeComparableColor(right?.[key])
  );
};

const paletteColorKeysEqual = (left, right, keys) => keys.every(
  (key) => normalizeComparableColor(left?.[key]) === normalizeComparableColor(right?.[key])
);

/**
 * @typedef {object} SvgStylesOptions
 * @property {Record<string, any>} state App state (state.js; not yet typed).
 * @property {(source: any, callback: (...args: any[]) => void, options?: Record<string, any>) => any} watch Vue `watch`
 * @property {(callback?: () => void) => Promise<void>} nextTick Vue `nextTick`
 * @property {((reason: string) => boolean) | null} [commitActiveResultEdit]
 *   The preview owner's commit of an edit to the displayed Result (R1, R13).
 * @property {(options?: { recolor?: Record<string, any>, prepareRules?: boolean }) => boolean | Promise<boolean>} projectPaletteAndRules
 *   The root's projection of the palette and the rules (R3), which prepares the
 *   rule matches and applies both through this owner.
 */

/** @param {SvgStylesOptions} options */
export const createSvgStyles = ({
  state,
  watch,
  nextTick,
  // R13: the preview owner's commit of an edit to the displayed Result.
  commitActiveResultEdit = null,
  // R13: the composition root's projection of the palette and the rules (R3),
  // which prepares the rule matches and applies both through this owner.
  projectPaletteAndRules
}) => {
  const {
    svgContent,
    extractedFeatures,
    featuresBySvgId,
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

  // The palette styles the displayed Result with the settings of the drawing
  // of that Result's mode: its rules, edits and skew slot colors (OV-108).
  const applyPaletteToSvg = ({
    recolorPairwise = false,
    recolorCollinear = false
  } = {}) => {
    const resultMode = (state.generatedMode?.value ?? mode.value) === 'linear' ? 'linear' : 'circular';
    const drawing = state.drawings[resultMode];
    const busy = state.sessionOperationAvailability?.();
    if (busy) return busy;
    if (!svgContent.value || !extractedFeatures.value.length) return;
    if (!svgContainer.value) return;
    const ruleMatches = ruleMatcher(drawing.manualSpecificRules);
    if (!ruleMatches.ready(extractedFeatures.value)) return;

    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return;

    const colors = appliedPaletteColors.value;
    const featurePaths = Array.from(getFeatureElementIndex(svg).values()).flat();
    const featureLookup = featuresBySvgId?.value || new Map();
    let updatedCount = 0;

    featurePaths.forEach((path) => {
      if (!isFeatureFillTarget(path)) return;
      const svgId = getFeatureIdentity(path);
      const feat = featureLookup.get(svgId);
      if (!feat) return;

      const paletteColor = colors[feat.type] || colors.default;
      if (!paletteColor) return;
      // A declined live match keeps the color Generate drew (R4).
      if (ruleMatches.declined(feat)) return;

      const hasSpecificRule = ruleMatches.matchesAny(feat) === true;

      if (!hasSpecificRule && !getFeatureOverride(drawing.featureColorOverrides, feat)) {
        const currentFill = path.getAttribute('fill');
        if (currentFill !== paletteColor) {
          path.setAttribute('fill', paletteColor);
          updatedCount++;
        }
      }
    });

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

    const featureLegendGroups = getAllFeatureLegendGroups(svg);
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
    const resolveLegendColor = (legendKey, palette) => {
      if (!legendKey) return null;
      const colorKey = keyToColorKey[legendKey];
      if (colorKey && palette[colorKey]) return palette[colorKey];
      if (palette[legendKey]) return palette[legendKey];
      return resolveOtherLegendColor(legendKey, palette);
    };
    featureLegendGroups.forEach((featureLegendGroup) => {
      if (!featureLegendGroup) return;

      const entryGroups = featureLegendGroup.querySelectorAll('g[data-legend-key]');

      if (entryGroups.length > 0) {
        entryGroups.forEach((entryGroup) => {
          const legendKey = entryGroup.getAttribute('data-legend-key');
          if (!legendKey) return;
          if (drawing.legendColorOverrides[legendKey]) return;

          const newColor = resolveLegendColor(legendKey, colors);
          if (!newColor) return;

          const paths = entryGroup.querySelectorAll('path');
          for (const path of paths) {
            const fill = path.getAttribute('fill');
            if (fill && fill !== 'none' && !fill.startsWith('url(')) {
              if (setColorAttributeIfChanged(path, 'fill', newColor)) updatedCount++;
              break;
            }
          }
        });
      } else {
        const texts = featureLegendGroup.querySelectorAll('text');
        const allPaths = featureLegendGroup.querySelectorAll('path');
        texts.forEach((textEl) => {
          const textContent = textEl.textContent?.trim();
          if (!textContent) return;
          if (drawing.legendColorOverrides[textContent]) return;

          const newColor = resolveLegendColor(textContent, colors);
          if (!newColor) return;

          const textPos = parseTransformXY(textEl.getAttribute('transform'));
          let bestPath = null;
          let bestX = -Infinity;
          for (const path of allPaths) {
            const pathPos = parseTransformXY(path.getAttribute('transform'));
            const fill = path.getAttribute('fill');
            if (
              Math.abs(pathPos.y - textPos.y) < 2 &&
              pathPos.x < textPos.x &&
              fill &&
              fill !== 'none' &&
              !fill.startsWith('url(')
            ) {
              if (pathPos.x > bestX) {
                bestX = pathPos.x;
                bestPath = path;
              }
            }
          }
          if (bestPath) {
            if (setColorAttributeIfChanged(bestPath, 'fill', newColor)) updatedCount++;
          }
        });
      }
    });

    if (colors.pairwise_match_min && colors.pairwise_match_max) {
      const allPairwiseLegends = svg.querySelectorAll(PAIRWISE_LEGEND_SELECTOR);
      allPairwiseLegends.forEach((pairwiseLegend) => {
        if (updatePairwiseLegendGradientStops(pairwiseLegend, colors)) updatedCount++;
      });
    }

    if (updatedCount > 0) commitActiveResultEdit?.('palette');
  };

  const applySpecificRulesToSvg = () => {
    const drawing = state.activeDrawing();
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    if (!svgContent.value || !extractedFeatures.value.length) return;
    if (!drawing.manualSpecificRules.length) return;
    if (!svgContainer.value) return;
    const ruleMatches = ruleMatcher(drawing.manualSpecificRules);
    if (!ruleMatches.ready(extractedFeatures.value)) return;

    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return;

    const featureElementIndex = getFeatureElementIndex(svg);
    let updatedCount = 0;

    extractedFeatures.value.forEach((feat) => {
      if (!feat.svg_id) return;
      // A declined live match keeps the color Generate drew (R4).
      if (ruleMatches.declined(feat)) return;

      const matchingRule = ruleMatches.first(feat);

      const elements = getFeatureFillElements(svg, feat.svg_id, featureElementIndex);
      if (elements.length > 0) {
        const newColor = matchingRule
          ? matchingRule.color
          : appliedPaletteColors.value[feat.type] || appliedPaletteColors.value.default;
        elements.forEach((el) => {
          if (el.getAttribute('fill') !== newColor) {
            el.setAttribute('fill', newColor);
            updatedCount++;
          }
        });
      }
    });

    if (updatedCount > 0) {
      commitActiveResultEdit?.('specific-rules');
      console.log(`Applied specific rules: updated ${updatedCount} elements`);
    }
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
    applySpecificRulesToSvg,
    applyTrackVisibility
  };
};
