// @ts-check
import { parseTransform } from './svg-transform.js';
import { isAutoFeatureUnderlay } from './feature-dom.js';
import { legendSwatch, sameLegendLine } from './legend-row-names.js';
/** @import { DrawnLegendName, LegendRowName, LegendSwatchIdentity } from './legend-row-names.js' */
import { LEGEND_ORDER_RECORD, pythonDrawnAttribute, resultBaseAttribute } from './result-paint-bases.js';
import { legendEntryKey } from './drawing/paint.js';

export { parseTransform };

/**
 * @param {Element | null | undefined} parent
 * @param {string} id
 * @returns {Element | null}
 */
export const getLegendChildById = (parent, id) => {
  if (!parent) return null;
  for (const child of parent.children || []) {
    if (child.id === id) return child;
    if (
      id === 'pairwise_legend' &&
      (
        child.getAttribute?.('data-gbdraw-role') === 'comparison-legend' ||
        /^pairwise_legend_(?:h|v)$/.test(String(child.id || ''))
      )
    ) {
      return child;
    }
  }
  return null;
};

export const PAIRWISE_LEGEND_SELECTOR =
  '[data-gbdraw-role="comparison-legend"][data-gbdraw-orientation="h"], [data-gbdraw-role="comparison-legend"][data-gbdraw-orientation="v"], #pairwise_legend, #pairwise_legend_h, #pairwise_legend_v';

export const COMPARISON_LEGEND_SELECTOR =
  `[data-gbdraw-role="comparison-legend"], ${PAIRWISE_LEGEND_SELECTOR}, #conservation_identity_legend`;

/** @param {Element | null | undefined} parent */
export const getComparisonLegendGroup = (parent) =>
  getLegendChildById(parent, 'pairwise_legend') ||
  getLegendChildById(parent, 'conservation_identity_legend');

/** @param {Element | null | undefined} el */
export const isInsideComparisonLegend = (el) =>
  Boolean(el?.closest?.(COMPARISON_LEGEND_SELECTOR));

/** @param {string | null | undefined} transform */
export const parseTransformXY = (transform) => {
  if (!transform) return { x: 0, y: 0 };
  const number = '([+-]?(?:\\d+(?:\\.\\d*)?|\\.\\d+)(?:[eE][+-]?\\d+)?)';
  const match = String(transform).match(
    new RegExp(`translate\\(\\s*${number}(?:\\s*,\\s*|\\s+)${number}\\s*\\)`)
  );
  return match ? { x: parseFloat(match[1]), y: parseFloat(match[2]) } : { x: 0, y: 0 };
};

/**
 * @param {SVGSVGElement | null | undefined} svg
 * @returns {Element[]}
 */
export const getAllFeatureLegendGroups = (svg) => {
  if (!svg) return [];

  const legendGroup = svg.getElementById('legend');
  if (!legendGroup) return [];

  const horizontalLegend = legendGroup.querySelector('#legend_horizontal');
  const verticalLegend = legendGroup.querySelector('#legend_vertical');

  /** @type {Element[]} */
  const groups = [];

  if (horizontalLegend && verticalLegend) {
    const hFeatureLegend = horizontalLegend.querySelector('#feature_legend_h');
    const vFeatureLegend = verticalLegend.querySelector('#feature_legend_v');
    if (hFeatureLegend) groups.push(hFeatureLegend);
    if (vFeatureLegend) groups.push(vFeatureLegend);
    return groups;
  }

  const featureLegendGroup = legendGroup.querySelector('#feature_legend');
  return featureLegendGroup ? [featureLegendGroup] : [legendGroup];
};

// The owner mark of a Legend row a rule commit shows before Python draws it
// (the executor's `ifAbsent` addition); the rules own it, not the editor.
export const SPECIFIC_COLOR_FILE_OWNER = 'specific-color-file';
// D-28: the marker of a Legend row the naming pass hides because an earlier
// row has its name and its swatch (svg-result-ingestion.js `nameLegendRows`);
// its value is the key of that earlier row, which the Legend draws for both.
export const LEGEND_MERGED_ROW = 'data-legend-merged';

// The one absence predicate of a Legend row: a row is absent when the Result
// lacks it or a Legend delete hid it (the executor keeps a deleted row of
// Python's, hidden, so a reconcile without the delete shows it again).
/** @param {Element | null | undefined} row */
export const legendRowShown = (row) => Boolean(row) && row?.getAttribute('display') !== 'none';

// D-35: the keys of the rows a Legend row of a Result stands for: its own key,
// then the keys of the rows the naming pass merged into it. A Legend color or
// stroke edit of the row edits each of them, so it stays one row.
/** @param {SVGSVGElement | null | undefined} svg @param {PythonLegendKey} key @returns {PythonLegendKey[]} */
export const legendRowKeys = (svg, key) => [...new Set([key, ...getAllFeatureLegendGroups(svg).flatMap((group) => (
  Array.from(group.querySelectorAll('g[data-legend-key]'))
    .filter((row) => row.getAttribute(LEGEND_MERGED_ROW) === key)
    .map((row) => /** @type {PythonLegendKey} */ (pythonLegendKey(row)))
))])];

// The keys of the feature Legend rows a Result shows (a renamed row's key,
// not the name it shows).
/** @param {SVGSVGElement | null | undefined} svg @returns {Set<string>} */
export const shownLegendKeys = (svg) => new Set(getAllFeatureLegendGroups(svg).flatMap((group) => (
  Array.from(group.querySelectorAll('g[data-legend-key]')).filter(legendRowShown).map(pythonLegendKey)
)));

// The order a feature Legend group recorded as Python's, else null.
/** @param {Element} group @returns {string[] | null} */
const recordedGroupOrder = (group) => {
  const recorded = group.getAttribute(LEGEND_ORDER_RECORD);
  if (!recorded) return null;
  try {
    const order = JSON.parse(recorded);
    return Array.isArray(order) ? order.map(String) : null;
  } catch {
    return null;
  }
};

// The rows of a key each feature Legend group of a Result holds, shown or hidden.
/** @param {Element} group @param {string} key */
const rowsOfKey = (group, key) => Array.from(group.querySelectorAll('g[data-legend-key]'))
  .filter((row) => pythonLegendKey(row) === key);

// Whether a Result shows Python's own row of a key: neither an editor row nor
// a row a rule commit showed before Python drew it.
/** @param {SVGSVGElement | null | undefined} svg @param {string} key */
export const shownPythonLegendRow = (svg, key) => getAllFeatureLegendGroups(svg).some((group) => (
  rowsOfKey(group, key).some((row) => legendRowShown(row) && !row.getAttribute('data-legend-owner'))
));

// Whether a Result shows the row of `nextKey` where Python's row of `key`
// was, as Python keeps a relabeled row in its place: the Result has no row of
// Python's `key` (Python drew the change), or its shown rows follow Python's
// order (recorded when the rows were first reordered, else the drawn order)
// with `nextKey` in the place of `key`. A commit's row appended before the
// rerender does not (U3a review M1).
/** @param {SVGSVGElement | null | undefined} svg @param {string} key @param {string} nextKey */
export const legendRowTakesPlace = (svg, key, nextKey) => getAllFeatureLegendGroups(svg).every((group) => {
  const rows = Array.from(group.querySelectorAll('g[data-legend-key]'))
    .filter((row) => row.getAttribute('data-legend-owner') !== 'direct-editor');
  const drawn = rows.map(pythonLegendKey);
  if (!drawn.includes(key)) return true;
  const relabeled = (/** @type {string[]} */ keys) => [...new Set(keys.map((each) => (each === key ? nextKey : each)))];
  const python = relabeled(recordedGroupOrder(group) ?? drawn);
  const shown = relabeled(rows.filter(legendRowShown).map(pythonLegendKey));
  const inOrder = python.filter((each) => shown.includes(each));
  return shown.filter((each) => python.includes(each)).every((each, at) => each === inOrder[at]);
});

// Python's key of a Legend row: the key it drew, which a rename keeps as a
// record; the row of an editor addition has its own key.
/** @param {Element} row */
export const pythonLegendKey = (row) => String(pythonDrawnAttribute(row, 'data-legend-key') ?? '').trim();

// Whether a Result's bytes hold Python's row of a Legend key (shown or hidden):
// a Session saved before the executor kept deleted rows lacks it (O-2).
/** @param {SVGSVGElement | null | undefined} svg @param {string} key */
export const drawsPythonLegendRow = (svg, key) => getAllFeatureLegendGroups(svg).some((group) => (
  Array.from(group.querySelectorAll('g[data-legend-key]'))
    .some((row) => row.getAttribute('data-legend-owner') !== 'direct-editor' && pythonLegendKey(row) === key)
));

// Python's order of a Result's feature Legend rows as their Python keys: the
// order recorded when the rows were first reordered, else null.
/** @param {SVGSVGElement | null | undefined} svg @returns {string[] | null} */
export const recordedLegendOrder = (svg) => {
  const group = getAllFeatureLegendGroups(svg).find((each) => each.getAttribute(LEGEND_ORDER_RECORD));
  return group ? recordedGroupOrder(group) : null;
};

// Whether a Result shows Legend structure edits (a renamed, deleted, or added
// row, or another order), read from its records, so a Result the renderer
// drew is laid out as Python lays the edited rows out.
/** @param {SVGSVGElement | null | undefined} svg */
export const legendStructureEdited = (svg) => getAllFeatureLegendGroups(svg).some((group) => (
  group.hasAttribute(LEGEND_ORDER_RECORD)
  || Boolean(group.querySelector([
    'g[data-legend-owner="direct-editor"]',
    `g[${resultBaseAttribute('data-legend-key')}]`,
    `g[data-legend-key][${resultBaseAttribute('display')}]`
  ].join(', ')))
));

/**
 * @param {SVGSVGElement | null | undefined} svg
 * @returns {Element | null}
 */
export const getVisibleFeatureLegendGroup = (svg) => {
  if (!svg) return null;

  const legendGroup = svg.getElementById('legend');
  if (!legendGroup) return null;

  // The two layout groups of a Legend are SVG groups, which carry `style`.
  const horizontalLegend = /** @type {SVGElement | null} */ (legendGroup.querySelector('#legend_horizontal'));
  const verticalLegend = /** @type {SVGElement | null} */ (legendGroup.querySelector('#legend_vertical'));

  if (horizontalLegend && verticalLegend) {
    const hVisible =
      horizontalLegend.style.display !== 'none' &&
      horizontalLegend.getAttribute('display') !== 'none';
    const targetLegend = hVisible ? horizontalLegend : verticalLegend;
    const featureGroup = hVisible
      ? targetLegend.querySelector('#feature_legend_h')
      : targetLegend.querySelector('#feature_legend_v');
    return featureGroup || targetLegend;
  }

  const featureLegendGroup = legendGroup.querySelector('#feature_legend');
  return featureLegendGroup || legendGroup;
};

/** @param {SVGSVGElement | null | undefined} svg */
export const isCurrentLegendHorizontal = (svg) => {
  const legendGroup = svg?.getElementById('legend');
  if (!legendGroup) return false;

  // The two layout groups of a Legend are SVG groups, which carry `style`.
  const horizontalLegend = /** @type {SVGElement | null} */ (legendGroup.querySelector('#legend_horizontal'));
  const verticalLegend = /** @type {SVGElement | null} */ (legendGroup.querySelector('#legend_vertical'));

  if (horizontalLegend && verticalLegend) {
    const hVisible =
      horizontalLegend.style.display !== 'none' &&
      horizontalLegend.getAttribute('display') !== 'none';
    return hVisible;
  }

  return false;
};

// The solid swatch of a Legend row: the path Python draws in the row's group,
// found by that structure, never by the fill it shows now (a row of color
// `none` has a swatch of fill `none`, OV-1385). A row Python draws with a
// paint server (a gradient bar, a hatch pattern) has none.
/**
 * @param {Element | null | undefined} entryGroup
 * @returns {Element | null}
 */
export const getLegendEntrySwatch = (entryGroup) => {
  const swatch = Array.from(entryGroup?.children || []).find((child) => child.localName === 'path') || null;
  return String(pythonDrawnAttribute(swatch, 'fill') ?? '').startsWith('url(') ? null : swatch;
};

// The stroke Python drew on the swatch of the Legend row of `key`, shown or
// merged into a shown row (D-35), which a row stroke edit keeps as
// `originalStroke*`.
/**
 * @param {SVGSVGElement | null | undefined} svg
 * @param {string} key
 * @returns {{ originalStrokeColor: string | null, originalStrokeWidth: number | null }}
 */
export const drawnLegendRowStroke = (svg, key) => {
  const swatch = getAllFeatureLegendGroups(svg)
    .map((group) => getLegendEntrySwatch(Array.from(group.querySelectorAll('g[data-legend-key]'))
      .find((entry) => pythonLegendKey(entry) === key && (legendRowShown(entry) || entry.hasAttribute(LEGEND_MERGED_ROW))) || null))
    .find(Boolean) || null;
  const width = Number.parseFloat(String(pythonDrawnAttribute(swatch, 'stroke-width') ?? ''));
  return {
    originalStrokeColor: pythonDrawnAttribute(swatch, 'stroke'),
    originalStrokeWidth: Number.isFinite(width) ? width : null
  };
};

// The stroke Python drew on the Result's first feature path that is not an
// automatic underlay: the Result's block-stroke default, which the Session
// format keeps as `originalSvgStroke`; null without one.
/**
 * @param {SVGSVGElement | null | undefined} svg
 * @returns {{ color: string | null, width: number | null } | null}
 */
export const drawnBlockStroke = (svg) => {
  const block = Array.from(svg?.querySelectorAll?.('path[id^="f"]') || [])
    .find((path) => !isAutoFeatureUnderlay(path));
  if (!block) return null;
  const width = Number.parseFloat(String(pythonDrawnAttribute(block, 'stroke-width') ?? ''));
  return { color: pythonDrawnAttribute(block, 'stroke'), width: Number.isFinite(width) ? width : null };
};

/** @param {unknown} value */
const paintKey = (value) => String(value ?? '').trim().toLowerCase();

/**
 * Whether a feature stroke edit (the feature popup) sets a stroke color or width.
 * @param {unknown} override
 */
export const setsFeatureStroke = (override) => {
  if (!override || typeof override !== 'object') return false;
  const { strokeColor, strokeWidth } = /** @type {Record<string, unknown>} */ (override);
  return Boolean(String(strokeColor ?? '').trim()) || (strokeWidth !== undefined && strokeWidth !== null && strokeWidth !== '');
};

// The identities and display values of a Legend row and a drawn feature,
// each a distinct name, so tsc flags a value of one passed or compared as
// another (D-15-22). Each is cast once, at the reader that produces it.
/** @typedef {string & { readonly __brand: 'PythonLegendKey' }} PythonLegendKey Python's key of a Legend row: the caption it drew, which a rename keeps as a record. */
/** @typedef {string & { readonly __brand: 'ShownLegendKey' }} ShownLegendKey The key a Result's Legend row shows (`data-legend-key`); a rename rewrites it. */
/** @typedef {string & { readonly __brand: 'SwatchColor' }} SwatchColor The fill of a Legend row's swatch. */
/** @typedef {string & { readonly __brand: 'RenderedFeatureId' }} RenderedFeatureId A drawn feature's rendered ID. */

/**
 * Python's row of a Result's Legend: its key and the color the renderer gives
 * the row's features, its swatch fill as Python drew it (a Legend color edit
 * keeps it as a record).
 * @typedef {object} PythonLegendRow
 * @property {PythonLegendKey} key
 * @property {SwatchColor | null} color Null when the row has no swatch fill.
 * @property {boolean} [ruleLine] False when Python drew the swatch with another line than a color-rule
 *   row's (the Result's block stroke): a rule named like the row then draws a row of its own (D-13).
 *   `pythonLegendRows` gives it.
 */


/**
 * A drawing's Legend records stored under a row's key, in the form Load reads
 * and a commit writes them. Each is optional.
 * @typedef {object} LegendKeyedRecords
 * @property {readonly { caption?: unknown, originalCaption?: unknown }[]} [entries]
 * @property {readonly { caption?: unknown, originalCaption?: unknown }[]} [deletedEntries]
 * @property {readonly { caption?: unknown, originalCaption?: unknown }[]} [dormantEntries]
 * @property {readonly string[]} [originalOrder] The generated inventory.
 * @property {Readonly<Record<string, unknown>>} [originalColors]
 * @property {Readonly<Record<string, unknown>>} [colorOverrides]
 * @property {Readonly<Record<string, unknown>>} [strokeOverrides]
 * @property {readonly string[]} [addedCaptions]
 * @property {readonly string[]} [strokeOptionsOpen] The rows whose Stroke options are open.
 */

// The one move of a Legend row's records from its old key to its new key (F2,
// OV-660), for the rows a re-key map names: a rule row whose key the
// allocation changes (a rule commit) or the reading of a saved caption
// changes (Load). Every collection is moved at once, so a swap of two keys
// moves both; a record the new key already has stays, and the old one is
// dropped. Returns new collections: none is written in place (History and the
// artifact slot hold them). `dropped` collects each dropped entry and edit
// (Load lists them).
/**
 * A Legend entry or edit the re-key dropped: the collection, its old key and the key it moved to.
 * @typedef {{ collection: 'entries' | 'deletedEntries' | 'dormantEntries' | 'colorOverrides' | 'strokeOverrides',
 *   from: string, to: string }} DroppedLegendRecord
 */
/**
 * @template {LegendKeyedRecords} Records
 * @param {Records} records
 * @param {ReadonlyMap<PythonLegendKey, PythonLegendKey>} rekeyed
 * @param {DroppedLegendRecord[] | null} [dropped]
 * @returns {Records}
 */
export const rekeyLegendRecords = (records, rekeyed, dropped = null) => {
  if (rekeyed.size === 0) return records;
  const lookup = /** @type {ReadonlyMap<string, string>} */ (rekeyed);
  /** @param {unknown} key */
  const moved = (key) => lookup.get(String(key ?? '').trim());
  /**
   * @template {{ caption?: unknown, originalCaption?: unknown }} Entry
   * @param {readonly Entry[]} list @param {DroppedLegendRecord['collection']} collection @returns {Entry[]}
   */
  const entries = (list, collection) => {
    /** @type {Set<string>} */
    const held = new Set(list.map((entry) => legendEntryKey(entry)).filter((key) => !lookup.has(key)));
    return list.flatMap((entry) => {
      const key = moved(legendEntryKey(entry));
      if (key && held.has(key)) {
        dropped?.push({ collection, from: legendEntryKey(entry), to: key });
        return [];
      }
      const caption = moved(entry.caption);
      const originalCaption = moved(entry.originalCaption);
      return [caption || originalCaption ? { ...entry, ...(caption ? { caption } : {}), ...(originalCaption ? { originalCaption } : {}) } : entry];
    });
  };
  /**
   * @param {Readonly<Record<string, unknown>>} record
   * @param {DroppedLegendRecord['collection'] | null} [collection] Set for an edit: a dropped one is collected.
   * @returns {Record<string, unknown>}
   */
  const keyed = (record, collection = null) => {
    const kept = Object.entries(record).filter(([key]) => !lookup.has(key));
    const result = Object.fromEntries(kept);
    Object.entries(record).forEach(([key, value]) => {
      const to = lookup.get(key);
      if (to === undefined) return;
      if (!Object.hasOwn(result, to)) result[to] = value;
      else if (collection) dropped?.push({ collection, from: key, to });
    });
    return result;
  };
  /** @param {readonly string[]} list @returns {string[]} */
  const keys = (list) => [...new Set(list.map((key) => lookup.get(key) ?? key))];
  /** @type {LegendKeyedRecords} */
  const next = { ...records };
  if (records.entries) next.entries = entries(records.entries, 'entries');
  if (records.deletedEntries) next.deletedEntries = entries(records.deletedEntries, 'deletedEntries');
  if (records.dormantEntries) next.dormantEntries = entries(records.dormantEntries, 'dormantEntries');
  if (records.originalOrder) next.originalOrder = keys(records.originalOrder);
  if (records.originalColors) next.originalColors = keyed(records.originalColors);
  if (records.colorOverrides) next.colorOverrides = keyed(records.colorOverrides, 'colorOverrides');
  if (records.strokeOverrides) next.strokeOverrides = keyed(records.strokeOverrides, 'strokeOverrides');
  if (records.addedCaptions) next.addedCaptions = keys(records.addedCaptions);
  if (records.strokeOptionsOpen) next.strokeOptionsOpen = keys(records.strokeOptionsOpen);
  return /** @type {Records} */ (next);
};

/** @param {Element} row @returns {PythonLegendRow} */
export const pythonLegendRow = (row) => {
  const color = String(pythonDrawnAttribute(getLegendEntrySwatch(row), 'fill') ?? '').trim();
  return {
    key: /** @type {PythonLegendKey} */ (pythonLegendKey(row)),
    color: color ? /** @type {SwatchColor} */ (color) : null
  };
};

// Python's rows of a mounted Result's feature Legend, by Python's key; a
// Python row wins over an editor row of the same key, as the executor finds them.
/** @param {SVGSVGElement | null | undefined} svg @returns {Map<PythonLegendKey, PythonLegendRow>} */
export const pythonLegendRows = (svg) => {
  /** @type {Map<PythonLegendKey, PythonLegendRow>} */
  const rows = new Map();
  const drawn = getAllFeatureLegendGroups(svg).flatMap((group) => Array.from(group.querySelectorAll('g[data-legend-key]')));
  const editor = (/** @type {Element} */ row) => row.getAttribute('data-legend-owner') === 'direct-editor';
  const block = drawnBlockStroke(svg);
  const ruleLine = { stroke: block?.color, width: block?.width };
  [...drawn.filter((row) => !editor(row)), ...drawn.filter(editor)].forEach((element) => {
    const row = pythonLegendRow(element);
    if (!row.key || rows.has(row.key)) return;
    const swatch = getLegendEntrySwatch(element);
    rows.set(row.key, {
      ...row,
      ruleLine: sameLegendLine({ stroke: pythonDrawnAttribute(swatch, 'stroke'), width: pythonDrawnAttribute(swatch, 'stroke-width') }, ruleLine)
    });
  });
  return rows;
};

// The three facts of a Legend row the naming rule reads and writes (D-09,
// legend-row-names.js), each read from the row at one place. Its name is what
// the user or Python wrote: the key a rename shows, else `data-legend-name`,
// which Python and a rule row's addition write where the name is not the key,
// else the key. Its swatch is what the row shows. Its drawn name is its text,
// which the executor writes (services/svg-result-ingestion.js).
/** @param {Element} row @returns {LegendRowName} */
export const legendRowName = (row) => /** @type {LegendRowName} */ (
  (row.hasAttribute(resultBaseAttribute('data-legend-key')) ? null : row.getAttribute('data-legend-name'))
    ?? String(row.getAttribute('data-legend-key') || '')
);
/** @param {Element} row @returns {LegendSwatchIdentity} */
export const legendRowSwatch = (row) => {
  const stops = Array.from(row.querySelectorAll('stop'));
  if (stops.length > 1) {
    return legendSwatch({ type: 'gradient', minColor: stops[0].getAttribute('stop-color'), maxColor: stops[stops.length - 1].getAttribute('stop-color') });
  }
  const swatch = getLegendEntrySwatch(row) || row.querySelector('path');
  return legendSwatch({ fill: swatch?.getAttribute('fill'), stroke: swatch?.getAttribute('stroke'), width: swatch?.getAttribute('stroke-width') });
};
/** @param {Element} row @returns {DrawnLegendName} */
export const drawnLegendRowName = (row) => /** @type {DrawnLegendName} */ (String(row.querySelector('text')?.textContent ?? ''));

// The Legends a Result names its rows in: one, or the horizontal and the
// vertical copy of a Result that holds both. The rows of one are named
// together, in document order, as Python names the rows of its Legend table.
/** @param {{ querySelector: (selector: string) => Element | null } | null | undefined} svg @returns {Element[]} */
export const namedLegends = (svg) => {
  const legend = svg?.querySelector('#legend');
  if (!legend) return [];
  const copies = [legend.querySelector('#legend_horizontal'), legend.querySelector('#legend_vertical')];
  return copies[0] && copies[1] ? [copies[0], copies[1]] : [legend];
};

/**
 * What the editor intent says of the features a stroked Legend row reaches in
 * one Result; Python's row gives the rest (`legendRowFeatureIds`).
 * @typedef {object} LegendRowReach
 * @property {readonly RenderedFeatureId[]} listedIds The features the row lists (`featureIds`, kept by older
 *   Sessions); when there are any, the row reaches only them.
 * @property {readonly RenderedFeatureId[]} namedIds The features a feature color edit names into the row.
 * @property {readonly RenderedFeatureId[]} ownStrokeIds The features with a stroke edit of their own.
 * @property {SwatchColor | null} draftColor The color the draft gives the row's features where the Result
 *   predates it (a live preview of the row's rule or palette color); null: Python's row color. Used only
 *   where the Result draws the row.
 */

// The features a stroke on a Legend row reaches in one Result (R3, PD-OI-066,
// OV-123, OV-288): the features the row lists, else the features a feature
// color edit names into the row and the features drawn in the color the
// renderer gives the row's features, Python's row color (R15-4), never the
// color the row's swatch shows. A row the Result does not draw reaches no
// feature by color, whatever color the draft gives its key, as Generate draws
// it (OV-1987). A feature with a stroke edit of its own keeps it: the
// explicit per-feature edit wins over the row. The executor reads it
// for every Result it styles (services/svg-result-ingestion.js), and the
// feature popup for the feature it shows.
/**
 * @param {LegendRowReach} reach
 * @param {SwatchColor | null} pythonColor Python's row color in this Result; null without the row.
 * @param {Iterable<readonly [RenderedFeatureId, string]>} drawnFills The fill each drawn feature shows.
 * @returns {RenderedFeatureId[]}
 */
export const legendRowFeatureIds = ({ listedIds, namedIds, ownStrokeIds, draftColor }, pythonColor, drawnFills) => {
  const listed = new Set(listedIds);
  const named = new Set(namedIds);
  const own = new Set(ownStrokeIds);
  const color = pythonColor === null ? null : paintKey(draftColor ?? pythonColor);
  /** @type {Set<RenderedFeatureId>} */
  const reached = new Set();
  for (const [id, fill] of drawnFills) {
    if (!id || own.has(id)) continue;
    if (listed.size > 0
      ? listed.has(id)
      : named.has(id) || (color && color !== 'none' && paintKey(fill) === color)) reached.add(id);
  }
  return [...reached];
};

/**
 * The anchor of a Legend entry: its text position, else its swatch position.
 * @param {Element | null | undefined} entryGroup
 */
export const legendEntryAnchor = (entryGroup) => {
  const groupOffset = parseTransformXY(entryGroup?.getAttribute?.('transform'));
  const target = entryGroup?.querySelector?.('text') || getLegendEntrySwatch(entryGroup);
  const targetOffset = parseTransformXY(target?.getAttribute?.('transform'));
  return { x: groupOffset.x + targetOffset.x, y: groupOffset.y + targetOffset.y };
};

/**
 * A row of a Result's shown feature Legend group, as the Result holds it.
 * @typedef {object} ResultLegendRow
 * @property {PythonLegendKey} key The row's key: Python's (the record a rename keeps), else an editor row's own.
 * @property {ShownLegendKey} shownKey The key the row shows.
 * @property {LegendRowName} name The row's name (`legendRowName`).
 * @property {DrawnLegendName} drawnName The name the Legend draws for it.
 * @property {LegendSwatchIdentity} swatch The swatch it shows (`legendRowSwatch`).
 * @property {PythonLegendKey | null} recordedKey Python's key that a rename keeps as a record; null without one.
 * @property {boolean} editor A row the Legend editor added.
 * @property {boolean} pythonShown Python drew it shown, merged or not: a row a rule commit retired is not.
 * @property {boolean} shown The row shows.
 * @property {boolean} merged The naming pass hides it: an earlier row has its name and its swatch (D-28).
 * @property {SwatchColor} color The fill its swatch shows.
 * @property {number} xPos
 * @property {number} yPos
 */

// The rows of a Result's shown feature Legend group in their visual order
// (row-major anchors); null without a Legend. The Legend entry owner lists
// and names them from the drawing's intent.
/** @param {SVGSVGElement | null | undefined} svg @returns {ResultLegendRow[] | null} */
export const resultLegendRows = (svg) => {
  const group = getVisibleFeatureLegendGroup(svg);
  if (!group) return null;
  return Array.from(group.querySelectorAll('g[data-legend-key]')).flatMap((row) => {
    const shownKey = /** @type {ShownLegendKey} */ (String(row.getAttribute('data-legend-key') || ''));
    if (!shownKey) return [];
    const anchor = legendEntryAnchor(row);
    const recordedKey = row.getAttribute(resultBaseAttribute('data-legend-key'));
    return [{
      key: /** @type {PythonLegendKey} */ (pythonLegendKey(row)),
      shownKey,
      name: legendRowName(row),
      swatch: legendRowSwatch(row),
      drawnName: drawnLegendRowName(row),
      recordedKey: recordedKey === null ? null : /** @type {PythonLegendKey} */ (recordedKey),
      editor: row.getAttribute('data-legend-owner') === 'direct-editor',
      // A merged row is Python's shown row the naming pass hides (D-35).
      pythonShown: row.hasAttribute(LEGEND_MERGED_ROW) || pythonDrawnAttribute(row, 'display') !== 'none',
      shown: legendRowShown(row),
      merged: row.hasAttribute(LEGEND_MERGED_ROW),
      color: /** @type {SwatchColor} */ (getLegendEntrySwatch(row)?.getAttribute('fill') || '#cccccc'),
      xPos: anchor.x,
      yPos: anchor.y
    }];
  }).sort((a, b) => {
    const yDelta = a.yPos - b.yPos;
    if (Math.abs(yDelta) >= 1) return yDelta;
    const xDelta = a.xPos - b.xPos;
    return Math.abs(xDelta) >= 1 ? xDelta : a.shownKey.localeCompare(b.shownKey, undefined, { sensitivity: 'base' });
  });
};

/**
 * What the editor shows for a Legend row instead of its key: the row's name,
 * and the name the Legend draws for it (its name while the row is hidden);
 * and what a typed name is matched with: whether the row shows, and its
 * swatch and fill.
 * @typedef {{ name: LegendRowName, drawnName: DrawnLegendName, shown: boolean, merged: boolean, swatch: LegendSwatchIdentity, color: SwatchColor }} ShownLegendRowNames
 */

// The names of the rows of a Result's shown feature Legend group, by each
// row's key, in legend order. The Legend tab, the popups, and the dialogs show
// a row by its drawn name, never by its key (PD-OI-042 scenario revision 2).
/** @param {SVGSVGElement | null | undefined} svg @returns {Map<PythonLegendKey, ShownLegendRowNames>} */
export const shownLegendRowNames = (svg) => new Map((resultLegendRows(svg) || []).map((row) => [
  row.key,
  {
    name: row.name,
    drawnName: row.shown && row.drawnName ? row.drawnName : /** @type {DrawnLegendName} */ (/** @type {string} */ (row.name)),
    shown: row.shown,
    merged: row.merged,
    swatch: row.swatch,
    color: row.color
  }
]));

/**
 * @param {Element} entryGroup
 * @param {number} xPos
 * @param {number} yPos
 */
export const moveLegendEntryToAnchor = (entryGroup, xPos, yPos) => {
  if (!Number.isFinite(xPos) || !Number.isFinite(yPos)) return false;
  const current = legendEntryAnchor(entryGroup);
  const deltaX = xPos - current.x;
  const deltaY = yPos - current.y;
  if (Math.abs(deltaX) < 1e-6 && Math.abs(deltaY) < 1e-6) return false;
  if (entryGroup.hasAttribute?.('transform')) {
    const groupOffset = parseTransformXY(entryGroup.getAttribute('transform'));
    entryGroup.setAttribute(
      'transform',
      `translate(${groupOffset.x + deltaX}, ${groupOffset.y + deltaY})`
    );
    return true;
  }
  const transformedChildren = Array.from(entryGroup.querySelectorAll?.('[transform]') || []);
  if (transformedChildren.length === 0) {
    entryGroup.setAttribute('transform', `translate(${deltaX}, ${deltaY})`);
    return true;
  }
  transformedChildren.forEach((node) => {
    const position = parseTransformXY(node.getAttribute('transform'));
    node.setAttribute('transform', `translate(${position.x + deltaX}, ${position.y + deltaY})`);
  });
  return true;
};

/**
 * Legend entries in the default order (D-08): the generated entries ranked by
 * the generated order (a renamed entry ranks by its generated caption),
 * followed by the other entries in their current order. Without a generated
 * order, null.
 * @template {{ caption?: unknown, originalCaption?: unknown }} T
 * @param {readonly T[] | null | undefined} entries
 * @param {readonly unknown[] | null | undefined} generatedOrder
 * @returns {T[] | null}
 */
export const defaultLegendEntryOrder = (entries, generatedOrder) => {
  /** @type {Map<string, number>} */
  const rank = new Map();
  (Array.isArray(generatedOrder) ? generatedOrder : []).forEach((caption) => {
    const key = String(caption ?? '').trim();
    if (key && !rank.has(key)) rank.set(key, rank.size);
  });
  if (rank.size === 0) return null;
  const list = Array.isArray(entries) ? entries : [];
  /** @param {T} entry */
  const generated = (entry) => String(entry?.originalCaption || entry?.caption || '').trim();
  return [
    ...list.filter((entry) => rank.has(generated(entry)))
      // Both ranks exist: the filter above keeps the ranked entries.
      .sort((left, right) => /** @type {number} */ (rank.get(generated(left))) - /** @type {number} */ (rank.get(generated(right)))),
    ...list.filter((entry) => !rank.has(generated(entry)))
  ];
};

/**
 * The row keys of a Result's default order (D-08, OV-47): its generated
 * inventory, followed by the keys of the other entries in their current
 * order. A Result applies it through `orderLegendEntries`, which skips the
 * keys it does not draw.
 * @param {readonly { caption?: unknown, originalCaption?: unknown }[]} entries
 * @param {readonly string[]} inventory
 * @returns {PythonLegendKey[]}
 */
export const defaultLegendKeyOrder = (entries, inventory) => {
  const generated = new Set(inventory);
  return /** @type {PythonLegendKey[]} */ ([...new Set([
    ...inventory,
    ...entries.map(legendEntryKey).filter((key) => !generated.has(key))
  ])].filter(Boolean));
};

/**
 * Whether Legend entries are in an edited order, not the default order (D-08).
 * Generate replays an edited order; the extraction of a Generate keeps the
 * generated order only while an edited order is replayed.
 * @template {{ caption?: unknown, originalCaption?: unknown }} T
 * @param {readonly T[]} entries
 * @param {readonly unknown[] | null | undefined} generatedOrder
 */
export const isLegendOrderEdited = (entries, generatedOrder) => {
  const defaultOrder = defaultLegendEntryOrder(entries, generatedOrder);
  return defaultOrder !== null && defaultOrder.some((entry, index) => entry !== entries[index]);
};

/**
 * @param {Element | null | undefined} targetGroup
 * @returns {Element[]}
 */
const legendEntryGroups = (targetGroup) => {
  const direct = Array.from(targetGroup?.children || []).filter(
    (child) => String(child.tagName || '').toLowerCase() === 'g' && child.hasAttribute?.('data-legend-key')
  );
  return direct.length > 0
    ? direct
    : Array.from(targetGroup?.querySelectorAll?.('g[data-legend-key]') || []);
};


// One Legend group's entries in slot order (row-major anchor order).
/** @param {Element | null | undefined} targetGroup */
const legendEntriesInSlotOrder = (targetGroup) => legendEntryGroups(targetGroup)
  .map((entryGroup, index) => ({ entryGroup, index, anchor: legendEntryAnchor(entryGroup) }))
  .sort((left, right) => {
    const yDelta = left.anchor.y - right.anchor.y;
    if (Math.abs(yDelta) >= 1) return yDelta;
    return left.anchor.x - right.anchor.x || left.index - right.index;
  });

// Whether the group's entries listed in the order already follow it.
/**
 * @param {readonly { entryGroup: Element }[]} positioned
 * @param {ReadonlyMap<string, number>} rank
 */
const followsOrder = (positioned, rank) => {
  const listed = positioned
    .map(({ entryGroup }) => rank.get(pythonLegendKey(entryGroup)))
    .filter((value) => value !== undefined);
  return listed.every((value, index) => index === 0 || listed[index - 1] < value);
};

/**
 * Put one Legend group's entries in the requested order of row keys (a
 * renamed row ranks by its key, as two rows may show one name). Entries keep
 * the group's current slots (row-major anchor order) and receive them in the
 * new order; keys missing from the order follow in their current order.
 * Live Sort and Move and the Generate replay of the order (D-08) share this.
 * With `keepFollowed`, a group whose listed entries already follow the order
 * keeps its order, so the entries the order does not list keep their places:
 * a displayed batch Result draws entries that the order, read from another
 * Result, does not list (D-07, B18).
 * Returns whether any entry moved or changed document order.
 * @param {Element} targetGroup
 * @param {readonly unknown[] | null | undefined} captionOrder
 * @param {{ keepFollowed?: boolean }} [options]
 */
export const orderLegendEntries = (targetGroup, captionOrder, { keepFollowed = false } = {}) => {
  const entryGroups = legendEntryGroups(targetGroup);
  if (entryGroups.length < 2) return false;
  /** @type {Map<string, number>} */
  const rank = new Map();
  (Array.isArray(captionOrder) ? captionOrder : []).forEach((caption) => {
    const key = String(caption ?? '').trim();
    if (key && !rank.has(key)) rank.set(key, rank.size);
  });
  const caption = pythonLegendKey;
  const positioned = legendEntriesInSlotOrder(targetGroup);
  if (keepFollowed && followsOrder(positioned, rank)) return false;
  const ordered = positioned
    .map((item, slot) => ({ ...item, slot }))
    .sort((left, right) => {
      // `rank.has` holds in the branch that reads `rank.get`.
      const leftRank = rank.has(caption(left.entryGroup)) ? /** @type {number} */ (rank.get(caption(left.entryGroup))) : Infinity;
      const rightRank = rank.has(caption(right.entryGroup)) ? /** @type {number} */ (rank.get(caption(right.entryGroup))) : Infinity;
      return leftRank === rightRank ? left.slot - right.slot : leftRank - rightRank;
    });
  let changed = false;
  ordered.forEach(({ entryGroup }, slot) => {
    const { x, y } = positioned[slot].anchor;
    changed = moveLegendEntryToAnchor(entryGroup, x, y) || changed;
  });
  ordered.forEach(({ entryGroup }, index) => {
    if (entryGroups[index] !== entryGroup) changed = true;
    targetGroup.appendChild(entryGroup);
  });
  return changed;
};
