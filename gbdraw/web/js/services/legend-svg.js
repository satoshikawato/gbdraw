// @ts-check
import { parseTransform } from './svg-transform.js';
import { isAutoFeatureUnderlay } from './feature-dom.js';
import { LEGEND_ORDER_RECORD, pythonDrawnAttribute, resultBaseAttribute } from './result-paint-bases.js';

export { parseTransform };

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

export const getComparisonLegendGroup = (parent) =>
  getLegendChildById(parent, 'pairwise_legend') ||
  getLegendChildById(parent, 'conservation_identity_legend');

export const isInsideComparisonLegend = (el) =>
  Boolean(el?.closest?.(COMPARISON_LEGEND_SELECTOR));

export const parseTransformXY = (transform) => {
  if (!transform) return { x: 0, y: 0 };
  const number = '([+-]?(?:\\d+(?:\\.\\d*)?|\\.\\d+)(?:[eE][+-]?\\d+)?)';
  const match = String(transform).match(
    new RegExp(`translate\\(\\s*${number}(?:\\s*,\\s*|\\s+)${number}\\s*\\)`)
  );
  return match ? { x: parseFloat(match[1]), y: parseFloat(match[2]) } : { x: 0, y: 0 };
};

export const getAllFeatureLegendGroups = (svg) => {
  if (!svg) return [];

  const legendGroup = svg.getElementById('legend');
  if (!legendGroup) return [];

  const horizontalLegend = legendGroup.querySelector('#legend_horizontal');
  const verticalLegend = legendGroup.querySelector('#legend_vertical');

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

// The one absence predicate of a Legend row: a row is absent when the Result
// lacks it or a Legend delete hid it (the executor keeps a deleted row of
// Python's, hidden, so a reconcile without the delete shows it again).
/** @param {Element | null | undefined} row */
export const legendRowShown = (row) => Boolean(row) && row?.getAttribute('display') !== 'none';

// The keys of the feature Legend rows a Result shows.
/** @param {Element | null | undefined} svg @returns {Set<string>} */
export const shownLegendKeys = (svg) => new Set(getAllFeatureLegendGroups(svg).flatMap((group) => (
  Array.from(group.querySelectorAll('g[data-legend-key]')).filter(legendRowShown)
    .map((row) => String(row.getAttribute('data-legend-key') || '').trim())
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
  .filter((row) => String(row.getAttribute('data-legend-key') || '').trim() === key);

// Whether a Result shows Python's own row of a key: neither an editor row nor
// a row a rule commit showed before Python drew it.
/** @param {Element | null | undefined} svg @param {string} key */
export const shownPythonLegendRow = (svg, key) => getAllFeatureLegendGroups(svg).some((group) => (
  rowsOfKey(group, key).some((row) => legendRowShown(row) && !row.getAttribute('data-legend-owner'))
));

// Whether a Result shows the row of `nextKey` where Python's row of `key`
// was, as Python keeps a relabeled row in its place: the Result has no row of
// Python's `key` (Python drew the change), or its shown rows follow Python's
// order (recorded when the rows were first reordered, else the drawn order)
// with `nextKey` in the place of `key`. A commit's row appended before the
// rerender does not (U3a review M1).
/** @param {Element | null | undefined} svg @param {string} key @param {string} nextKey */
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
/** @param {Element | null | undefined} svg @param {string} key */
export const drawsPythonLegendRow = (svg, key) => getAllFeatureLegendGroups(svg).some((group) => (
  Array.from(group.querySelectorAll('g[data-legend-key]'))
    .some((row) => row.getAttribute('data-legend-owner') !== 'direct-editor' && pythonLegendKey(row) === key)
));

// Python's order of a Result's feature Legend rows as their Python keys: the
// order recorded when the rows were first reordered, else null.
/** @param {Element | null | undefined} svg @returns {string[] | null} */
export const recordedLegendOrder = (svg) => {
  const group = getAllFeatureLegendGroups(svg).find((each) => each.getAttribute(LEGEND_ORDER_RECORD));
  return group ? recordedGroupOrder(group) : null;
};

// Whether a Result shows Legend structure edits (a renamed, deleted, or added
// row, or another order), read from its records, so a Result the renderer
// drew is laid out as Python lays the edited rows out.
/** @param {Element | null | undefined} svg */
export const legendStructureEdited = (svg) => getAllFeatureLegendGroups(svg).some((group) => (
  group.hasAttribute(LEGEND_ORDER_RECORD)
  || Boolean(group.querySelector([
    'g[data-legend-owner="direct-editor"]',
    `g[${resultBaseAttribute('data-legend-key')}]`,
    `g[data-legend-key][${resultBaseAttribute('display')}]`
  ].join(', ')))
));

export const getVisibleFeatureLegendGroup = (svg) => {
  if (!svg) return null;

  const legendGroup = svg.getElementById('legend');
  if (!legendGroup) return null;

  const horizontalLegend = legendGroup.querySelector('#legend_horizontal');
  const verticalLegend = legendGroup.querySelector('#legend_vertical');

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

export const isCurrentLegendHorizontal = (svg) => {
  const legendGroup = svg?.getElementById('legend');
  if (!legendGroup) return false;

  const horizontalLegend = legendGroup.querySelector('#legend_horizontal');
  const verticalLegend = legendGroup.querySelector('#legend_vertical');

  if (horizontalLegend && verticalLegend) {
    const hVisible =
      horizontalLegend.style.display !== 'none' &&
      horizontalLegend.getAttribute('display') !== 'none';
    return hVisible;
  }

  return false;
};

export const getLegendEntrySwatch = (entryGroup) => Array.from(
  entryGroup?.querySelectorAll?.('path') || []
).find((path) => {
  const fill = path.getAttribute('fill');
  return fill && fill !== 'none' && !fill.startsWith('url(');
}) || null;

// The stroke Python drew on the swatch of the Legend row `caption`, which a
// row stroke edit keeps as `originalStroke*`.
/**
 * @param {Element | null | undefined} svg
 * @param {string} caption
 * @returns {{ originalStrokeColor: string | null, originalStrokeWidth: number | null }}
 */
export const drawnLegendRowStroke = (svg, caption) => {
  const swatch = getAllFeatureLegendGroups(svg)
    .map((group) => getLegendEntrySwatch(Array.from(group.querySelectorAll('g[data-legend-key]'))
      .find((entry) => entry.getAttribute('data-legend-key') === caption && legendRowShown(entry)) || null))
    .find(Boolean) || null;
  const width = Number.parseFloat(String(pythonDrawnAttribute(swatch, 'stroke-width') ?? ''));
  return {
    originalStrokeColor: pythonDrawnAttribute(swatch, 'stroke'),
    originalStrokeWidth: Number.isFinite(width) ? width : null
  };
};

// The stroke Python drew on the Result's first feature path that is not an
// automatic underlay, which a Session keeps as `originalSvgStroke` for the
// readers of Sessions saved before the executor's records; null without one.
/**
 * @param {Element | null | undefined} svg
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
 */

// The key of a listed Legend row: its key record (`originalCaption`), else
// its caption (a row saved before U3a has no record).
/** @param {{ caption?: unknown, originalCaption?: unknown }} entry @returns {PythonLegendKey} */
export const legendEntryKey = (entry) => /** @type {PythonLegendKey} */ (String(entry.originalCaption || entry.caption || '').trim());

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
/** @param {Element | null | undefined} svg @returns {Map<PythonLegendKey, PythonLegendRow>} */
export const pythonLegendRows = (svg) => {
  /** @type {Map<PythonLegendKey, PythonLegendRow>} */
  const rows = new Map();
  const drawn = getAllFeatureLegendGroups(svg).flatMap((group) => Array.from(group.querySelectorAll('g[data-legend-key]')));
  const editor = (/** @type {Element} */ row) => row.getAttribute('data-legend-owner') === 'direct-editor';
  [...drawn.filter((row) => !editor(row)), ...drawn.filter(editor)].forEach((element) => {
    const row = pythonLegendRow(element);
    if (row.key && !rows.has(row.key)) rows.set(row.key, row);
  });
  return rows;
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
 *   predates it (a live preview of the row's rule or palette color); null: Python's row color.
 */

// The features a stroke on a Legend row reaches in one Result (R3, PD-OI-066,
// OV-123, OV-288): the features the row lists, else the features a feature
// color edit names into the row and the features drawn in the color the
// renderer gives the row's features, Python's row color (R15-4), never the
// color the row's swatch shows. A feature with a stroke edit of its own keeps
// it: the explicit per-feature edit wins over the row. The executor reads it
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
  const color = paintKey(draftColor ?? pythonColor);
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

/** The anchor of a Legend entry: its text position, else its swatch position. */
export const legendEntryAnchor = (entryGroup) => {
  const groupOffset = parseTransformXY(entryGroup?.getAttribute?.('transform'));
  const target = entryGroup?.querySelector?.('text') || getLegendEntrySwatch(entryGroup);
  const targetOffset = parseTransformXY(target?.getAttribute?.('transform'));
  return { x: groupOffset.x + targetOffset.x, y: groupOffset.y + targetOffset.y };
};

/**
 * A row of a Result's shown feature Legend group, as the Result holds it.
 * @typedef {object} ResultLegendRow
 * @property {ShownLegendKey} shownKey The key the row shows.
 * @property {PythonLegendKey | null} recordedKey Python's key that a rename keeps as a record; null without one.
 * @property {boolean} editor A row the Legend editor added.
 * @property {boolean} pythonShown Python drew it shown: a row a rule commit retired is not.
 * @property {boolean} shown The row shows.
 * @property {SwatchColor} color The fill its swatch shows.
 * @property {number} xPos
 * @property {number} yPos
 */

// The rows of a Result's shown feature Legend group in their visual order
// (row-major anchors); null without a Legend. The Legend entry owner lists
// and names them from the drawing's intent.
/** @param {Element | null | undefined} svg @returns {ResultLegendRow[] | null} */
export const resultLegendRows = (svg) => {
  const group = getVisibleFeatureLegendGroup(svg);
  if (!group) return null;
  return Array.from(group.querySelectorAll('g[data-legend-key]')).flatMap((row) => {
    const shownKey = /** @type {ShownLegendKey} */ (String(row.getAttribute('data-legend-key') || ''));
    if (!shownKey) return [];
    const anchor = legendEntryAnchor(row);
    const recordedKey = row.getAttribute(resultBaseAttribute('data-legend-key'));
    return [{
      shownKey,
      recordedKey: recordedKey === null ? null : /** @type {PythonLegendKey} */ (recordedKey),
      editor: row.getAttribute('data-legend-owner') === 'direct-editor',
      pythonShown: pythonDrawnAttribute(row, 'display') !== 'none',
      shown: legendRowShown(row),
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
 */
export const defaultLegendEntryOrder = (entries, generatedOrder) => {
  const rank = new Map();
  (Array.isArray(generatedOrder) ? generatedOrder : []).forEach((caption) => {
    const key = String(caption ?? '').trim();
    if (key && !rank.has(key)) rank.set(key, rank.size);
  });
  if (rank.size === 0) return null;
  const list = Array.isArray(entries) ? entries : [];
  const generated = (entry) => String(entry?.originalCaption || entry?.caption || '').trim();
  return [
    ...list.filter((entry) => rank.has(generated(entry)))
      .sort((left, right) => rank.get(generated(left)) - rank.get(generated(right))),
    ...list.filter((entry) => !rank.has(generated(entry)))
  ];
};

/**
 * The captions of a Result's default order (D-08, OV-47): its generated
 * inventory with each renamed entry's current caption, followed by the other
 * entries in their current order. A Result applies it through
 * `orderLegendEntries`, which skips the captions it does not draw.
 * @param {{ caption?: string, originalCaption?: string }[]} entries
 * @param {string[]} inventory
 * @returns {string[]}
 */
export const defaultLegendCaptionOrder = (entries, inventory) => {
  const generated = new Set(inventory);
  const generatedCaption = (entry) => String(entry?.originalCaption || entry?.caption || '').trim();
  const currentCaption = new Map(entries.map((entry) => [generatedCaption(entry), String(entry?.caption || '').trim()]));
  return [...new Set([
    ...inventory.map((caption) => currentCaption.get(caption) || caption),
    ...entries.filter((entry) => !generated.has(generatedCaption(entry))).map((entry) => String(entry?.caption || '').trim())
  ])].filter(Boolean);
};

/**
 * Whether Legend entries are in an edited order, not the default order (D-08).
 * Generate replays an edited order; the extraction of a Generate keeps the
 * generated order only while an edited order is replayed.
 */
export const isLegendOrderEdited = (entries, generatedOrder) => {
  const defaultOrder = defaultLegendEntryOrder(entries, generatedOrder);
  return defaultOrder !== null && defaultOrder.some((entry, index) => entry !== entries[index]);
};

const legendEntryGroups = (targetGroup) => {
  const direct = Array.from(targetGroup?.children || []).filter(
    (child) => String(child.tagName || '').toLowerCase() === 'g' && child.hasAttribute?.('data-legend-key')
  );
  return direct.length > 0
    ? direct
    : Array.from(targetGroup?.querySelectorAll?.('g[data-legend-key]') || []);
};

const legendEntryCaption = (entryGroup) => String(entryGroup.getAttribute('data-legend-key') || '').trim();

// One Legend group's entries in slot order (row-major anchor order).
const legendEntriesInSlotOrder = (targetGroup) => legendEntryGroups(targetGroup)
  .map((entryGroup, index) => ({ entryGroup, index, anchor: legendEntryAnchor(entryGroup) }))
  .sort((left, right) => {
    const yDelta = left.anchor.y - right.anchor.y;
    if (Math.abs(yDelta) >= 1) return yDelta;
    return left.anchor.x - right.anchor.x || left.index - right.index;
  });

// Whether the group's entries listed in the order already follow it.
const followsOrder = (positioned, rank) => {
  const listed = positioned
    .map(({ entryGroup }) => rank.get(legendEntryCaption(entryGroup)))
    .filter((value) => value !== undefined);
  return listed.every((value, index) => index === 0 || listed[index - 1] < value);
};

/**
 * Put one Legend group's entries in the requested caption order. Entries keep
 * the group's current slots (row-major anchor order) and receive them in the
 * new order; captions missing from the order follow in their current order.
 * Live Sort and Move and the Generate replay of the order (D-08) share this.
 * With `keepFollowed`, a group whose listed entries already follow the order
 * keeps its order, so the entries the order does not list keep their places:
 * a displayed batch Result draws entries that the order, read from another
 * Result, does not list (D-07, B18).
 * Returns whether any entry moved or changed document order.
 */
export const orderLegendEntries = (targetGroup, captionOrder, { keepFollowed = false } = {}) => {
  const entryGroups = legendEntryGroups(targetGroup);
  if (entryGroups.length < 2) return false;
  const rank = new Map();
  (Array.isArray(captionOrder) ? captionOrder : []).forEach((caption) => {
    const key = String(caption ?? '').trim();
    if (key && !rank.has(key)) rank.set(key, rank.size);
  });
  const caption = legendEntryCaption;
  const positioned = legendEntriesInSlotOrder(targetGroup);
  if (keepFollowed && followsOrder(positioned, rank)) return false;
  const ordered = positioned
    .map((item, slot) => ({ ...item, slot }))
    .sort((left, right) => {
      const leftRank = rank.has(caption(left.entryGroup)) ? rank.get(caption(left.entryGroup)) : Infinity;
      const rightRank = rank.has(caption(right.entryGroup)) ? rank.get(caption(right.entryGroup)) : Infinity;
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
