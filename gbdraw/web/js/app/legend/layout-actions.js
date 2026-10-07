// @ts-check
/** @import { LayoutBox, LegendFontMetrics, LegendLayoutRow, LegendEntryLayout } from '../../services/legend-layout.js' */
import {
  buildCircularLegendLayout,
  buildLinearLegendLayout,
  circularLegendLocalBounds,
  linearLegendLocalBounds,
  loadLegendFontMetrics,
  resolveBundledFontFace
} from '../../services/legend-layout.js';
import { getComparisonLegendGroup, getLegendEntrySwatch } from '../../services/legend-svg.js';
import { parseCompositionMetadata } from '../legend-layout/composition-actions.js';

// The dpi Python lays the Legend out at unless the Result says otherwise
// (gbdraw/data/config.toml `canvas.dpi`); Results written before Python
// recorded its Legend layout inputs (`legendReflow.dpi`) were laid out at it.
const DEFAULT_LEGEND_DPI = 96;

/**
 * The Legend layout inputs Python records in the composition metadata
 * (`legendReflow`), when the Result has them.
 * @typedef {object} RecordedLegendInputs
 * @property {number} colorRectSize
 * @property {string | null} [fontFile]
 * @property {number} [fontSize]
 * @property {number} [dpi]
 * @property {number} [wrapWidth]
 */

/** @param {Element | null | undefined} parent */
const entryGroupsOf = (parent) => Array.from(parent?.children || []).filter(
  (child) => child.localName === 'g' && child.hasAttribute('data-legend-key')
);

/** @param {Element} element @param {string} id */
const childById = (element, id) => element.querySelector(`#${CSS.escape(id)}`);

/**
 * @param {Element | null} element
 * @param {number} x
 * @param {number} y
 */
const setTranslate = (element, x, y) => {
  element?.setAttribute('transform', `translate(${x}, ${y})`);
};

/**
 * The solid rows of one feature Legend group, in document order (the
 * Legend's order), as Python's legend table holds them.
 * @param {Element} group
 * @returns {{ rows: LegendLayoutRow[], entries: Element[], label: Element | null }}
 */
const readSolidRows = (group) => {
  const entries = entryGroupsOf(group);
  const rows = entries.map((entry) => {
    const swatch = getLegendEntrySwatch(entry);
    return /** @type {LegendLayoutRow} */ ({
      key: String(entry.getAttribute('data-legend-key') || ''),
      type: 'solid',
      stroke: swatch?.getAttribute('stroke') || 'none',
      strokeWidth: Number(swatch?.getAttribute('stroke-width') || 0)
    });
  });
  const label = entries.map((entry) => entry.querySelector('text')).find(Boolean) || null;
  return { rows, entries, label };
};

/**
 * The gradient rows of a comparison or conservation Legend: their keys, bar
 * strokes, and the lowest identity its scale shows.
 * @param {Element | null} group
 * @returns {LegendLayoutRow[]}
 */
const readGradientRows = (group) => {
  if (!group) return [];
  const labels = Array.from(group.querySelectorAll('text'))
    .map((text) => String(text.textContent || '').trim());
  const minLabel = labels.find((label) => /^-?\d+(?:\.\d+)?%$/.test(label) && label !== '100%');
  const minValue = minLabel ? Number(minLabel.slice(0, -1)) : 0;
  return entryGroupsOf(group).map((entry) => {
    const bar = Array.from(entry.querySelectorAll('path'))
      .find((path) => String(path.getAttribute('fill') || '').startsWith('url(')) || entry.querySelector('path');
    return /** @type {LegendLayoutRow} */ ({
      key: String(entry.getAttribute('data-legend-key') || ''),
      type: 'gradient',
      stroke: bar?.getAttribute('stroke') || 'none',
      strokeWidth: Number(bar?.getAttribute('stroke-width') || 0),
      minValue
    });
  });
};

/**
 * Write one solid entry's positions as Python's Legend writers do: the swatch
 * and the caption each carry a translate; the entry group carries none.
 * @param {Element} entry
 * @param {LegendEntryLayout} layout
 */
const placeEntry = (entry, layout) => {
  entry.removeAttribute('transform');
  setTranslate(getLegendEntrySwatch(entry) || entry.querySelector('path'), layout.rectX, layout.rectY);
  setTranslate(entry.querySelector('text'), layout.textX, layout.textY);
};

/**
 * @param {Element[]} entries
 * @param {LegendEntryLayout[]} layouts
 */
const placeEntries = (entries, layouts) => {
  const byKey = new Map(layouts.map((layout) => [layout.key, layout]));
  entries.forEach((entry) => {
    const layout = byKey.get(String(entry.getAttribute('data-legend-key') || ''));
    if (layout) placeEntry(entry, layout);
  });
};

// Python's Legend, laid out in the Web (zero shift): the rows the Legend holds
// now, in their document order, measured and placed by the port of Python's
// layout (services/legend-layout.js) with the inputs Python recorded for the
// Result. An unedited Legend keeps every position; an edited one gets the
// positions Python gives the same rows. The layout owner
// (app/legend-layout/reposition-actions.js) docks it with the bounds returned.
export const createLegendLayoutActions = () => {
  /** @type {LegendFontMetrics | null} */
  let fontMetrics = null;

  /**
   * Load the bundled-font metrics the layout reads, once. The mount binder
   * awaits it for a Result with a Legend, so the Legend edits that follow,
   * which are synchronous, find them loaded.
   * @returns {Promise<LegendFontMetrics>}
   */
  const prepareLegendLayout = async () => {
    fontMetrics = await loadLegendFontMetrics();
    return fontMetrics;
  };
  /** Whether the metrics are loaded, so a caller need not wait for them. */
  const isLegendLayoutReady = () => fontMetrics !== null;

  /**
   * Lay out the Legend of `svg` for `side` (default: its composition side) as
   * Python does, and return its local bounds for the composition. Null when
   * there is nothing to lay out (no Legend, a side of `none`, no rows, an
   * unknown structure) or the font metrics are not loaded yet.
   * @param {Element | null | undefined} svg
   * @param {{ side?: string }} [options]
   * @returns {LayoutBox | null}
   */
  const layOutLegend = (svg, { side } = {}) => {
    const legendGroup = /** @type {SVGSVGElement | null | undefined} */ (svg)?.getElementById?.('legend');
    if (!svg || !legendGroup) return null;
    if (!fontMetrics) {
      void prepareLegendLayout().catch(() => {});
      return null;
    }
    const metrics = fontMetrics;
    const metadata = parseCompositionMetadata(svg);
    const recorded = /** @type {RecordedLegendInputs | null} */ (metadata.legendReflow);
    if (!recorded) {
      throw new Error('This diagram has no legend reflow metadata. Regenerate it before editing the legend.');
    }
    const legendSide = side || metadata.legendSide;
    if (!legendSide || legendSide === 'none') return null;
    const horizontal = childById(legendGroup, 'legend_horizontal');
    const vertical = childById(legendGroup, 'legend_vertical');
    const linear = Boolean(horizontal && vertical);

    const featureGroup = linear
      ? childById(/** @type {Element} */ (horizontal), 'feature_legend_h')
      : (childById(legendGroup, 'feature_legend') || legendGroup);
    if (!featureGroup) return null;
    const solid = readSolidRows(featureGroup);
    const gradientGroup = linear
      ? getComparisonLegendGroup(horizontal)
      : childById(legendGroup, 'conservation_identity_legend');
    const rows = [...solid.rows, ...readGradientRows(gradientGroup)];
    if (rows.length === 0) return null;

    const sampleLabel = solid.label || gradientGroup?.querySelector('text') || null;
    const options = {
      side: legendSide,
      // Before Python recorded its inputs, the wrap width is the primary
      // item's width: what Linear passes; Circular passed its content width.
      wrapWidth: recorded.wrapWidth ?? metadata.primary.finalBounds.width,
      fontFile: recorded.fontFile || resolveBundledFontFace(metrics, sampleLabel?.getAttribute('font-family') || ''),
      fontSize: recorded.fontSize ?? Number(sampleLabel?.getAttribute('font-size')),
      dpi: recorded.dpi ?? DEFAULT_LEGEND_DPI,
      colorRectSize: recorded.colorRectSize
    };
    if (!Number.isFinite(options.fontSize) || options.fontSize <= 0) return null;

    if (linear) {
      const layout = buildLinearLegendLayout(rows, options, metrics);
      /** @type {Array<['horizontal' | 'vertical', Element, string]>} */
      const orientations = [
        ['horizontal', /** @type {Element} */ (horizontal), 'h'],
        ['vertical', /** @type {Element} */ (vertical), 'v']
      ];
      orientations.forEach(([name, group, suffix]) => {
        const orientation = layout[name];
        const featureLegend = childById(group, `feature_legend_${suffix}`);
        if (featureLegend) {
          if (orientation.featureX || orientation.featureY) {
            setTranslate(featureLegend, orientation.featureX, orientation.featureY);
          } else {
            featureLegend.removeAttribute('transform');
          }
          placeEntries(entryGroupsOf(featureLegend), orientation.feature.entries);
        }
        const comparison = getComparisonLegendGroup(group);
        if (comparison && orientation.gradient) {
          comparison.setAttribute('transform', suffix === 'h'
            ? `translate(${orientation.gradientX}, 0)${orientation.gradientY ? ` translate(0, ${orientation.gradientY})` : ''}`
            : `translate(0, ${orientation.gradientY})`);
        }
      });
      return linearLegendLocalBounds(layout, options.colorRectSize);
    }

    const layout = buildCircularLegendLayout(rows, options, metrics);
    placeEntries(solid.entries, layout.entries);
    const rect = options.colorRectSize;
    const outline = Array.from(legendGroup.children)
      .find((child) => child.localName === 'path' && child.getAttribute('fill') === 'none');
    outline?.setAttribute('d', `M 0,${-0.5 * rect} L ${layout.width},${-0.5 * rect} `
      + `L ${layout.width},${layout.height - 0.5 * rect} L 0,${layout.height - 0.5 * rect} z`);
    if (gradientGroup && layout.gradient) setTranslate(gradientGroup, layout.gradientX, layout.gradientY);
    return circularLegendLocalBounds(layout, rect);
  };

  return { isLegendLayoutReady, layOutLegend, prepareLegendLayout };
};
