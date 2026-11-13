// @ts-check
// Hover on the similarity rings of the shown Result (PD-OI-110 scenario
// revision 3): over a ring, that ring keeps its strength and every other ring
// is drawn at 35%, its Legend row is bold, a radial line at the pointer's angle
// spans all the rings and reaches the outer edge of the features, a popup gives the
// position (bp) and every ring's identity there, outermost first; over a ring's
// Legend name, the same rings and row without the line and the popup. Leaving
// restores everything; a click pins nothing.
// On the Multi-Record Canvas the hovered series keeps its strength in every
// record (`data-source-index`).
//
// Nothing is written into the mounted root (R1; Export and Save read it): the
// dimming and the bold row are rules of a <style> element beside the Result,
// and the line and the popup are an overlay over it; both are in the page only
// while a hover shows, leave it when the hover ends or the Result leaves the
// page, and never print. Every value comes from what Python
// drew: the ring groups, their hits' drawn spans, identities and paint, their
// radii (`drawn-attributes.generated.js`), the record's axis
// (`services/region-axis.js`) for the position, and each ring's name as its
// Legend row draws it (`services/legend-svg.js`, else the name Python wrote on
// the ring).
/** @import { CircularAxisGeometry, RegionAxisRecord } from '../../services/region-axis.js' */
import { HIT_DISPLAY_SPAN_ATTRIBUTES, RING_LEGEND_NAME_ATTRIBUTE, RING_RADIUS_ATTRIBUTES } from '../../services/drawn-attributes.generated.js';
import { legendCaptionRuns } from '../../services/legend-layout.js';
import { drawnLegendRowName } from '../../services/legend-svg.js';
import { axisPosition, framePoint, regionAxisFrame, sourceBaseAt } from '../../services/region-axis.js';
import { hexToRgb, resolveColorToHex } from '../../utils/color-utils.js';
import { numberText, t } from '../../utils/i18n.js';

/** A similarity ring's group in a Result. */
export const RING_SELECTOR = 'g[data-gbdraw-slot-renderer="sequence_conservation"]';
/** A ring's row of the similarity-ring Legend block. */
export const RING_LEGEND_ROW_SELECTOR = '#conservation_identity_legend g[data-source-index]';
/** The opacity of the rings other than the hovered one (the receipt's 35%). */
export const DIMMED_RING_OPACITY = 0.35;
/** The surface the preview mounts the Result in; the hover's rules apply only there. */
const PREVIEW_SURFACE = '.gbdraw-preview-surface';

/**
 * One drawn ring of a record, as Python drew it.
 * @typedef {object} DrawnRing
 * @property {Element} element
 * @property {string} sourceIndex The ring's series (`data-source-index`), the same in every record.
 * @property {string} name The ring's name as its Legend row draws it (`<i>` runs included).
 * @property {number} inner
 * @property {number} outer
 */

/**
 * A ring's hit at a drawn position: its identity, and the paint it is drawn with there.
 * @typedef {{ identity: number, color: string }} RingHit
 */

/** @typedef {{ swatch: HTMLElement, value: HTMLElement }} PopupCells A popup row's swatch and identity cells. */

/**
 * What the pointer is over: a ring at a position of its record, or a ring's Legend row.
 * @typedef {{ kind: 'ring', sourceIndex: string, ring: Element } | { kind: 'legend', sourceIndex: string } | null} RingHoverTarget
 */

/**
 * What a hovered ring's popup reads: its record's axis and frame, and the record's rings outermost first.
 * @typedef {{ axis: CircularAxisGeometry, frame: SVGGraphicsElement, rings: DrawnRing[] }} RingPopupSource
 */

/** @param {string | null} value */
const numberOrNull = (value) => (value === null || value.trim() === '' ? null : (Number.isFinite(Number(value)) ? Number(value) : null));

/**
 * The CSS that shows a hover of the series `sourceIndex` on screen: every other ring at 35%, the series' Legend row
 * bold.
 * @param {string} sourceIndex
 */
export const ringHoverStyle = (sourceIndex) => {
  const series = `[data-source-index="${CSS.escape(sourceIndex)}"]`;
  return '@media screen {\n'
    + `${PREVIEW_SURFACE} ${RING_SELECTOR}:not(${series}) { opacity: ${DIMMED_RING_OPACITY}; }\n`
    + `${PREVIEW_SURFACE} ${RING_LEGEND_ROW_SELECTOR}${series} text { font-weight: 700; }\n`
    + '}\n';
};

/**
 * Each ring's name as the Legend draws it: its Legend row's text where the Result draws one (it follows a rename
 * made in the Web), else the name Python wrote on the ring (`LegendTable.ring_names`, written also when no Legend is
 * drawn); null for a ring with neither.
 * @param {ParentNode} root @param {Element[]} rings
 * @returns {(string | null)[]}
 */
const ringNames = (root, rings) => {
  const rows = new Map(Array.from(root.querySelectorAll(RING_LEGEND_ROW_SELECTOR))
    .map((row) => [row.getAttribute('data-track-index'), row]));
  return rings.map((ring) => {
    const row = rows.get(ring.getAttribute('data-track-index'));
    return row ? drawnLegendRowName(row) : ring.getAttribute(RING_LEGEND_NAME_ATTRIBUTE);
  });
};

/**
 * The rings a record draws, outermost first by the radii Python wrote; null for a Result drawn before Python wrote
 * them and the rings' names (gbdraw 0.14.5 or earlier).
 * @param {ParentNode} root The mounted root, which holds the Legend.
 * @param {ParentNode} record The record's group, or the root of a Result with one record.
 * @returns {DrawnRing[] | null}
 */
export const ringsOfRecord = (root, record) => {
  const elements = Array.from(record.querySelectorAll(RING_SELECTOR));
  const radii = elements.map((element) => [RING_RADIUS_ATTRIBUTES.inner, RING_RADIUS_ATTRIBUTES.outer]
    .map((name) => numberOrNull(element.getAttribute(name))));
  const names = ringNames(root, elements);
  if (radii.some((pair) => pair.includes(null)) || names.includes(null)) return null;
  return elements.map((element, index) => ({
    element,
    sourceIndex: element.getAttribute('data-source-index') ?? '',
    name: String(names[index]),
    inner: Number(radii[index][0]),
    outer: Number(radii[index][1])
  })).sort((a, b) => b.outer - a.outer);
};

/** @param {string} color @param {number} opacity The colour a paint shows at an opacity, as CSS. */
const paintShown = (color, opacity) => {
  const rgb = opacity < 1 ? hexToRgb(resolveColorToHex(color) || '') : null;
  return rgb ? `rgba(${rgb.r}, ${rgb.g}, ${rgb.b}, ${opacity})` : color;
};

/**
 * A ring's hit at a drawn position: of its hits whose drawn span holds it, the last in document order (Python
 * paints that one on top), with its identity and the colour it shows (its fill at its fill opacity); null when no
 * hit does or Python wrote no spans.
 * @param {Element} ring
 * @param {number} position Bases from display index 0.
 * @returns {RingHit | null}
 */
export const ringHitAt = (ring, position) => {
  /** @type {RingHit | null} */
  let top = null;
  for (const hit of ring.querySelectorAll(`path[${HIT_DISPLAY_SPAN_ATTRIBUTES.start}]`)) {
    const start = numberOrNull(hit.getAttribute(HIT_DISPLAY_SPAN_ATTRIBUTES.start));
    const end = numberOrNull(hit.getAttribute(HIT_DISPLAY_SPAN_ATTRIBUTES.end));
    const identity = numberOrNull(hit.getAttribute('data-identity'));
    if (start === null || end === null || identity === null || position < start || position >= end) continue;
    top = { identity, color: paintShown(hit.getAttribute('fill') ?? '', numberOrNull(hit.getAttribute('fill-opacity')) ?? 1) };
  }
  return top;
};

/**
 * An identity as the popup shows it: one decimal, grouped like the position, and never 100.0 below 100.
 * @param {number | null} identity
 */
export const identityText = (identity) => {
  if (identity === null) return t('No hit');
  const shown = identity < 100 ? Math.min(Math.round(identity * 10) / 10, 99.9) : identity;
  return `${numberText(shown, 1)}%`;
};

/**
 * The record a ring is drawn in: its index and group (the root for a Result whose groups name no record).
 * @param {Element} ring @param {Element} root
 */
const recordOfRing = (ring, root) => {
  const group = ring.closest('[data-gbdraw-record-index]');
  return group && root.contains(group)
    ? { index: Number(group.getAttribute('data-gbdraw-record-index')), group }
    : { index: 0, group: root };
};

/**
 * What the popup of a hovered ring reads, or null when the Result cannot show it: no axis for the ring's record (a
 * Result loaded from a Session the `gbdraw` command saved) or no radii (a Result drawn by gbdraw 0.14.5 or earlier).
 * @param {Element} root @param {Element} ring The ring's group, or an element in it.
 * @param {RegionAxisRecord[]} axisRecords
 * @returns {RingPopupSource | null}
 */
export const ringPopupSource = (root, ring, axisRecords) => {
  const group = ring.closest(RING_SELECTOR);
  if (!group || !root.contains(group)) return null;
  const record = recordOfRing(group, root);
  const axis = axisRecords.find((entry) => entry.recordIndex === record.index)?.axis;
  const frame = regionAxisFrame(root, record.index);
  const rings = axis?.kind === 'circular' && frame ? ringsOfRecord(root, record.group) : null;
  return rings && frame && axis?.kind === 'circular' ? { axis, frame, rings } : null;
};

/**
 * Whether a ring (or a hit of it) shows the ring popup, which then replaces the match hover summary of its hits.
 * @param {Element} root @param {Element} element @param {RegionAxisRecord[]} axisRecords
 */
export const ringPopupShows = (root, element, axisRecords) => ringPopupSource(root, element, axisRecords) !== null;

/**
 * Where the pointer is on a hovered ring's record: the source base (R5), the drawn position, and the radial line's
 * ends in overlay coordinates: from the smallest to the largest of the rings' radii and the outer edge of the feature
 * band, so the line spans all the rings and reaches that edge. Null when the frame is not on screen.
 * @param {{ source: RingPopupSource, clientX: number, clientY: number, host: DOMRect | { left: number, top: number } }} input
 */
export const ringPointAt = ({ source: { axis, frame, rings }, clientX, clientY, host }) => {
  const mapped = framePoint(frame, clientX, clientY);
  const matrix = frame.getScreenCTM();
  if (!mapped || !matrix) return null;
  const position = axisPosition(axis, mapped.point).position;
  const angle = Math.atan2(mapped.point.y, mapped.point.x);
  const radii = [...rings.flatMap((entry) => [entry.inner, entry.outer]), axis.bandOuterPx];
  const end = (/** @type {number} */ radius) => new DOMPoint(radius * Math.cos(angle), radius * Math.sin(angle)).matrixTransform(matrix);
  const [from, to] = [end(Math.min(...radii)), end(Math.max(...radii))];
  return {
    base: sourceBaseAt(axis.source, Math.floor(position)),
    position,
    line: { x1: from.x - host.left, y1: from.y - host.top, x2: to.x - host.left, y2: to.y - host.top }
  };
};

/**
 * What a pointer event is over: a ring, else the ring Legend row whose name's box holds the point (the gaps between
 * its letters count; its identity ramp does not).
 * @param {Element} root @param {EventTarget | null} target @param {number} clientX @param {number} clientY
 * @returns {RingHoverTarget}
 */
export const ringHoverTarget = (root, target, clientX, clientY) => {
  const ring = target instanceof Element ? target.closest(RING_SELECTOR) : null;
  if (ring && root.contains(ring)) return { kind: 'ring', sourceIndex: ring.getAttribute('data-source-index') ?? '', ring };
  for (const row of root.querySelectorAll(RING_LEGEND_ROW_SELECTOR)) {
    const box = row.querySelector('text')?.getBoundingClientRect();
    if (box && clientX >= box.left && clientX <= box.right && clientY >= box.top && clientY <= box.bottom) {
      return { kind: 'legend', sourceIndex: row.getAttribute('data-source-index') ?? '' };
    }
  }
  return null;
};

/**
 * The popup's rows, outermost ring first: built once for a record's rings and the hovered series, so a move rewrites
 * only the swatches and identities (the vendored Tailwind runtime rescans every added node, OV-933).
 * @param {HTMLElement} list @param {DrawnRing[]} rings @param {string} sourceIndex
 * @returns {PopupCells[]}
 */
const buildPopupRows = (list, rings, sourceIndex) => {
  const values = rings.map((ring) => {
    const row = document.createElement('div');
    row.className = 'gbdraw-ring-hover-row';
    if (ring.sourceIndex === sourceIndex) row.setAttribute('data-hovered', '');
    const swatch = document.createElement('span');
    swatch.className = 'gbdraw-ring-hover-swatch';
    const name = document.createElement('span');
    name.className = 'gbdraw-ring-hover-name truncate';
    name.append(...legendCaptionRuns(ring.name).map((run) => {
      if (!run.italic) return run.text;
      const italic = document.createElement('i');
      italic.textContent = run.text;
      return italic;
    }));
    const value = document.createElement('span');
    value.className = 'gbdraw-ring-hover-value';
    row.append(swatch, name, value);
    return { row, swatch, value };
  });
  list.replaceChildren(...values.map(({ row }) => row));
  return values.map(({ swatch, value }) => ({ swatch, value }));
};

/**
 * A row's swatch and identity: the colour drawn at the position and its identity, or an empty swatch and No hit.
 * @param {PopupCells} cells @param {RingHit | null} hit
 */
const writeRingHit = ({ swatch, value }, hit) => {
  swatch.style.background = hit?.color ?? 'transparent';
  value.textContent = identityText(hit?.identity ?? null);
  value.toggleAttribute('data-no-hit', hit === null);
};

/**
 * @typedef {object} RingHoverOverlay
 * @property {HTMLStyleElement} style The rules of the dimming and the bold Legend row.
 * @property {HTMLDivElement} element The line and the popup, shown only at a ring point.
 * @property {SVGLineElement} line
 * @property {HTMLDivElement} popup
 * @property {HTMLDivElement} heading
 * @property {HTMLDivElement} list
 */

/**
 * The overlay beside the Result: the <style> of the dimming, the radial line, and the popup.
 * It is in the page only while a hover shows (a page without a hover holds none of it).
 * @param {HTMLElement} host
 * @returns {RingHoverOverlay}
 */
const createRingHoverOverlay = (host) => {
  const style = document.createElement('style');
  style.setAttribute('data-gbdraw-ring-hover-style', '');
  const element = document.createElement('div');
  element.className = 'gbdraw-ring-hover-overlay';
  element.setAttribute('aria-hidden', 'true');
  const drawing = document.createElementNS('http://www.w3.org/2000/svg', 'svg');
  const line = document.createElementNS('http://www.w3.org/2000/svg', 'line');
  line.setAttribute('class', 'gbdraw-ring-hover-line');
  drawing.append(line);
  const popup = document.createElement('div');
  popup.className = 'gbdraw-ring-hover-popup';
  const heading = document.createElement('div');
  heading.className = 'gbdraw-ring-hover-position';
  const list = document.createElement('div');
  popup.append(heading, list);
  element.append(drawing, popup);
  element.hidden = true;
  host.append(style, element);
  return { style, element, line, popup, heading, list };
};

/**
 * @typedef {object} RingHoverPorts
 * @property {Element} svg The mounted root of the shown Result.
 * @property {() => RegionAxisRecord[]} axisRecords The records of the shown Result whose axis Python wrote.
 * @property {() => boolean} held Whether hovers wait: a popup is open, a drag runs, a review holds the popups.
 * @property {() => boolean} isPreviewTransformActive Whether the preview pans or zooms.
 */

/**
 * Binds the ring hover to a mounted root; returns its release.
 * @param {RingHoverPorts} ports
 * @returns {() => void}
 */
export const attachRingHover = ({ svg, axisRecords, held, isPreviewTransformActive }) => {
  const host = svg.parentElement?.parentElement?.parentElement;
  if (!(host instanceof HTMLElement)) return () => {};
  /** @type {RingHoverOverlay | null} */
  let overlay = null;
  /** @type {MouseEvent | null} */
  let pointer = null;
  let frame = 0;
  let shownStyle = '';
  /** @type {{ rings: DrawnRing[], sourceIndex: string, cells: PopupCells[] }} */
  let shownRows = { rings: [], sourceIndex: '', cells: [] };
  const clear = () => {
    overlay?.style.remove();
    overlay?.element.remove();
    overlay = null;
    shownStyle = '';
    shownRows = { rings: [], sourceIndex: '', cells: [] };
  };
  // The Result leaves the page (another Result, drawing or mode, a Session reset): the hover goes with it.
  const unmounted = new MutationObserver(() => {
    if (!svg.isConnected) clear();
  });
  [svg.parentElement, svg.parentElement?.parentElement].forEach((node) => node && unmounted.observe(node, { childList: true }));
  const update = () => {
    frame = 0;
    const event = pointer;
    const target = event && svg.isConnected && !held() && !isPreviewTransformActive()
      ? ringHoverTarget(svg, event.target, event.clientX, event.clientY) : null;
    if (!event || !target) return clear();
    const shown = overlay ??= createRingHoverOverlay(host);
    const css = ringHoverStyle(target.sourceIndex);
    if (css !== shownStyle) shown.style.textContent = shownStyle = css;
    const box = host.getBoundingClientRect();
    const source = target.kind === 'ring' ? ringPopupSource(svg, target.ring, axisRecords()) : null;
    const point = source ? ringPointAt({ source, clientX: event.clientX, clientY: event.clientY, host: box }) : null;
    shown.element.hidden = !point;
    if (!source || !point) return;
    const { rings } = source;
    Object.entries(point.line).forEach(([name, value]) => shown.line.setAttribute(name, value.toFixed(1)));
    shown.heading.textContent = t('{position} bp', { position: numberText(point.base) });
    if (shownRows.sourceIndex !== target.sourceIndex || shownRows.rings.length !== rings.length
      || shownRows.rings.some((ring, index) => ring.element !== rings[index].element || ring.name !== rings[index].name)) {
      shownRows = { rings, sourceIndex: target.sourceIndex, cells: buildPopupRows(shown.list, rings, target.sourceIndex) };
    }
    rings.forEach((ring, index) => writeRingHit(shownRows.cells[index], ringHitAt(ring.element, point.position)));
    const { offsetWidth: width, offsetHeight: height } = shown.popup;
    const x = event.clientX - box.left + 16;
    const y = event.clientY - box.top + 16;
    shown.popup.style.left = `${Math.max(4, x + width > box.width - 4 ? x - width - 32 : x)}px`;
    shown.popup.style.top = `${Math.max(4, Math.min(y, box.height - height - 4))}px`;
  };
  /** @param {Event} event */
  const move = (event) => {
    if (!(event instanceof MouseEvent)) return;
    pointer = event;
    if (!frame) frame = window.requestAnimationFrame(update);
  };
  const leave = () => {
    pointer = null;
    if (!frame) frame = window.requestAnimationFrame(update);
  };
  svg.addEventListener('mousemove', move);
  svg.addEventListener('mouseleave', leave);
  svg.addEventListener('wheel', leave, { passive: true });
  svg.addEventListener('mousedown', leave);
  return () => {
    svg.removeEventListener('mousemove', move);
    svg.removeEventListener('mouseleave', leave);
    svg.removeEventListener('wheel', leave);
    svg.removeEventListener('mousedown', leave);
    if (frame) window.cancelAnimationFrame(frame);
    unmounted.disconnect();
    clear();
  };
};
