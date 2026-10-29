// @ts-check
/** @import { CircularGradientLayout, CircularRingLayout, LayoutBox, LegendFontMetrics, LegendLayoutRow, LegendEntryLayout, LinearGradientLayout, RegisteredLegendFace } from '../../services/legend-layout.js' */
import {
  buildCircularLegendLayout,
  buildLinearLegendLayout,
  circularLegendLocalBounds,
  linearLegendLocalBounds,
  legendMeasurementFace,
  loadLegendFontMetrics,
  withRegisteredFontFaces
} from '../../services/legend-layout.js';
import { drawnLegendRowName, getComparisonLegendGroup, getLegendEntrySwatch, legendRowShown, pythonLegendKey } from '../../services/legend-svg.js';
import { parseCompositionMetadata } from '../legend-layout/composition-actions.js';

/**
 * The Legend layout inputs Python records in the composition metadata
 * (`legendReflow`), when the Result has them.
 * @typedef {object} RecordedLegendInputs
 * @property {number} colorRectSize
 * @property {string | null} [fontFile]
 * @property {number} [fontSize]
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
 * Legend's order), as Python's legend table holds them: the rows shown.
 * @param {Element} group
 * @returns {{ rows: LegendLayoutRow[], entries: Element[], label: Element | null }}
 */
const readSolidRows = (group) => {
  const entries = entryGroupsOf(group).filter(legendRowShown);
  const rows = entries.map((entry) => {
    const swatch = getLegendEntrySwatch(entry);
    return /** @type {LegendLayoutRow} */ ({
      key: pythonLegendKey(entry),
      name: drawnLegendRowName(entry),
      type: 'solid',
      stroke: swatch?.getAttribute('stroke') || 'none',
      strokeWidth: Number(swatch?.getAttribute('stroke-width') || 0)
    });
  });
  const label = entries.map((entry) => entry.querySelector('text')).find(Boolean) || null;
  return { rows, entries, label };
};

/**
 * The bar of a gradient entry: the path Python fills with the gradient.
 * @param {Element} entry
 * @returns {Element | null}
 */
const gradientBar = (entry) => Array.from(entry.children)
  .find((child) => child.localName === 'path' && String(child.getAttribute('fill') || '').startsWith('url(')) || null;

/**
 * The rows of a comparison or conservation Legend block and the lowest
 * identity its scale shows: the similarity-ring block's ring rows (a row with
 * `data-track-index`; whether it shows a swatch), or the gradient rows of a
 * comparison and of a Result drawn before the ring block (their bar strokes).
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
    if (entry.hasAttribute('data-track-index')) {
      return /** @type {LegendLayoutRow} */ ({
        key: pythonLegendKey(entry),
        name: drawnLegendRowName(entry),
        type: 'ring',
        stroke: 'none',
        strokeWidth: 0,
        minValue,
        swatch: Array.from(entry.children).some((child) => child.localName === 'path')
      });
    }
    const bar = gradientBar(entry);
    return /** @type {LegendLayoutRow} */ ({
      key: pythonLegendKey(entry),
      name: drawnLegendRowName(entry),
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
    const layout = byKey.get(pythonLegendKey(entry));
    if (layout) placeEntry(entry, layout);
  });
};

/**
 * Where a gradient block's parts go, in the block's own coordinates: the
 * positions Python's Legend layout gives them (legend/{circular,linear}_layout.py),
 * which its Legend writers draw.
 * @typedef {object} GradientPlacement
 * @property {boolean} compact Several gradients: one label column, one scale below the bars.
 * @property {Array<{ key: string, titleX: number, titleY: number, barX: number, barY: number }>} entries
 *   Each entry's title (a compact block's label) and bar.
 * @property {number} minLabelX
 * @property {number} maxLabelX
 * @property {number} scaleLabelY
 */

/** @param {LinearGradientLayout} gradient @returns {GradientPlacement} */
const linearGradientPlacement = (gradient) => ({
  compact: gradient.compact,
  entries: gradient.entries,
  minLabelX: gradient.minLabelX,
  maxLabelX: gradient.maxLabelX,
  scaleLabelY: gradient.scaleLabelY
});

/**
 * A compact entry's label is its title, and its bar sits at the block's `barX`.
 * @param {CircularGradientLayout} gradient
 * @returns {GradientPlacement}
 */
const circularGradientPlacement = (gradient) => ({
  compact: gradient.compact,
  entries: gradient.compact
    ? gradient.compactEntries.map((entry) => ({
      key: entry.key, titleX: entry.labelX, titleY: entry.labelY, barX: gradient.barX, barY: entry.barY
    }))
    : gradient.singleEntries,
  minLabelX: gradient.minLabelX,
  maxLabelX: gradient.maxLabelX,
  scaleLabelY: gradient.scaleLabelY
});

/**
 * Write a gradient block's parts (OV-1741): each entry's title and bar, and
 * the scale labels, which Python writes after a compact block's entries, or
 * after the title and bar of a single entry, minimum first, then 100%.
 * @param {Element} group
 * @param {GradientPlacement} placement
 */
const placeGradient = (group, placement) => {
  const entries = entryGroupsOf(group);
  const byKey = new Map(placement.entries.map((layout) => [layout.key, layout]));
  entries.forEach((entry) => {
    const layout = byKey.get(pythonLegendKey(entry));
    if (!layout) return;
    setTranslate(Array.from(entry.children).find((child) => child.localName === 'text') || null, layout.titleX, layout.titleY);
    setTranslate(gradientBar(entry), layout.barX, layout.barY);
  });
  const scaleParent = placement.compact ? group : entries[0];
  const texts = Array.from(scaleParent?.children || []).filter((child) => child.localName === 'text');
  const [minLabel, maxLabel] = placement.compact ? texts : texts.slice(1);
  setTranslate(minLabel || null, placement.minLabelX, placement.scaleLabelY);
  setTranslate(maxLabel || null, placement.maxLabelX, placement.scaleLabelY);
};

/**
 * Write the similarity-ring block's parts where Python's layout puts them
 * (legend/circular_layout.py `_build_ring_layout`): the `outer` and `inner`
 * captions, each ring row's swatch and name, and the band's bar, `identity`
 * label and scale labels, in the order Python writes them.
 * @param {Element} group
 * @param {CircularRingLayout} rings
 */
const placeRings = (group, rings) => {
  const children = Array.from(group.children);
  const [outer, inner] = children.filter((child) => child.localName === 'text');
  setTranslate(outer || null, 0, rings.outerCaptionY);
  setTranslate(inner || null, 0, rings.innerCaptionY);
  const byKey = new Map(rings.rows.map((row) => [row.key, row]));
  entryGroupsOf(group).forEach((entry) => {
    const layout = byKey.get(pythonLegendKey(entry));
    if (!layout) return;
    const parts = Array.from(entry.children);
    setTranslate(parts.find((child) => child.localName === 'path') || null, layout.swatchX, layout.y);
    setTranslate(parts.find((child) => child.localName === 'text') || null, layout.textX, layout.y);
  });
  const band = children.find((child) => child.getAttribute('data-gbdraw-role') === 'conservation-identity-band');
  const bandParts = Array.from(band?.children || []);
  setTranslate(bandParts.find((child) => child.localName === 'path') || null, 0, rings.bandY);
  const [label, minLabel, maxLabel] = bandParts.filter((child) => child.localName === 'text');
  setTranslate(label || null, rings.bandLabelX, rings.bandY);
  setTranslate(minLabel || null, rings.minLabelX, rings.scaleLabelY);
  setTranslate(maxLabel || null, rings.maxLabelX, rings.scaleLabelY);
};

// Python's Legend, laid out in the Web (zero shift): the rows the Legend holds
// now, in their document order, measured and placed by the port of Python's
// layout (services/legend-layout.js) with the inputs Python recorded for the
// Result. An unedited Legend keeps every position; an edited one gets the
// positions Python gives the same rows. The layout owner
// (app/legend-layout/reposition-actions.js) docks it with the bounds returned.
/**
 * @typedef {object} LegendLayoutActionOptions
 * @property {() => readonly RegisteredLegendFace[]} readRegisteredFontFaces
 *   The faces of the fonts registered in the page (the font store's), so a Legend in such a font is
 *   measured with it.
 */

/** @param {LegendLayoutActionOptions} options */
export const createLegendLayoutActions = ({ readRegisteredFontFaces }) => {
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
    const metrics = withRegisteredFontFaces(fontMetrics, readRegisteredFontFaces());
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
      fontFile: legendMeasurementFace(metrics, recorded.fontFile, sampleLabel?.getAttribute('font-family') || ''),
      fontSize: recorded.fontSize ?? Number(sampleLabel?.getAttribute('font-size')),
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
          placeGradient(comparison, linearGradientPlacement(orientation.gradient));
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
    if (gradientGroup && (layout.rings || layout.gradient)) setTranslate(gradientGroup, layout.gradientX, layout.gradientY);
    if (gradientGroup && layout.rings) placeRings(gradientGroup, layout.rings);
    else if (gradientGroup && layout.gradient) placeGradient(gradientGroup, circularGradientPlacement(layout.gradient));
    return circularLegendLocalBounds(layout, rect);
  };

  return { isLegendLayoutReady, layOutLegend, prepareLegendLayout };
};
