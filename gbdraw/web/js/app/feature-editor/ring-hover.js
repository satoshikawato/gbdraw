// @ts-check
// Hover on the similarity rings of the shown Result (PD-OI-110 scenario
// revision 3): over a ring, that ring keeps its strength and every other ring
// is drawn at 35%, its Legend row is bold, a radial line marks the pointer's
// angle, and a popup gives the position (bp) and every ring's identity there,
// outermost first; over a ring's Legend name, the same rings and row without
// the line and the popup. Leaving restores everything; a click pins nothing.
// On the Multi-Record Canvas the hovered series keeps its strength in every
// record (`data-source-index`).
//
// Nothing is written into the mounted root (R1; Export and Save read it): the
// dimming and the bold row are rules of a <style> element beside the Result,
// and the line and the popup are an overlay over it. Every value comes from
// what Python drew: the ring groups, their hits' drawn spans, identities and
// fills, their radii (`drawn-attributes.generated.js`), and the record's axis
// (`services/region-axis.js`) for the position.
/** @import { CircularAxisGeometry, RegionAxisRecord } from '../../services/region-axis.js' */
import { HIT_DISPLAY_SPAN_ATTRIBUTES, RING_RADIUS_ATTRIBUTES } from '../../services/drawn-attributes.generated.js';
import { axisPosition, framePoint, regionAxisFrame, sourceBaseAt } from '../../services/region-axis.js';
import { formatGroupedNumber } from '../../services/record-options.js';
import { t } from '../../utils/i18n.js';

/** A similarity ring's group in a Result. */
export const RING_SELECTOR = 'g[data-gbdraw-slot-renderer="sequence_conservation"]';
/** A ring's row of the similarity-ring Legend block. */
export const RING_LEGEND_ROW_SELECTOR = '#conservation_identity_legend g[data-source-index]';
/** The opacity of the rings other than the hovered one (the receipt's 35%). */
export const DIMMED_RING_OPACITY = 0.35;
/** The surface the preview mounts the Result in; the hover's rules apply only there. */
const PREVIEW_SURFACE = '.gbdraw-preview-surface';

/**
 * One drawn ring of a record, as Python drew it. `inner`/`outer` are null on a Result drawn before Python wrote them.
 * @typedef {object} DrawnRing
 * @property {Element} element
 * @property {string} sourceIndex The ring's series (`data-source-index`), the same in every record.
 * @property {string} name The ring's name, as its Legend row shows it.
 * @property {number | null} inner
 * @property {number | null} outer
 */

/** @typedef {{ swatch: HTMLElement, value: HTMLElement }} PopupCells A popup row's swatch and identity cells. */

/**
 * What the pointer is over: a ring at a position of its record, or a ring's Legend row.
 * @typedef {{ kind: 'ring', sourceIndex: string, ring: Element } | { kind: 'legend', sourceIndex: string } | null} RingHoverTarget
 */

/** @param {string | null} value */
const numberOrNull = (value) => (value === null || value.trim() === '' ? null : (Number.isFinite(Number(value)) ? Number(value) : null));

/**
 * The CSS that shows a hover of the series `sourceIndex`: every other ring at 35%, the series' Legend row bold.
 * @param {string} sourceIndex
 */
export const ringHoverStyle = (sourceIndex) => {
  const series = `[data-source-index="${CSS.escape(sourceIndex)}"]`;
  return `${PREVIEW_SURFACE} ${RING_SELECTOR}:not(${series}) { opacity: ${DIMMED_RING_OPACITY}; }\n`
    + `${PREVIEW_SURFACE} ${RING_LEGEND_ROW_SELECTOR}${series} text { font-weight: 700; }\n`;
};

/**
 * The rings a record draws, outermost first (by the radii Python wrote; in document order when it wrote none).
 * @param {ParentNode} record The record's group, or the root of a Result with one record.
 * @returns {DrawnRing[]}
 */
export const ringsOfRecord = (record) => {
  const rings = Array.from(record.querySelectorAll(RING_SELECTOR)).map((element) => ({
    element,
    sourceIndex: element.getAttribute('data-source-index') ?? '',
    name: element.getAttribute('data-track-label') ?? '',
    inner: numberOrNull(element.getAttribute(RING_RADIUS_ATTRIBUTES.inner)),
    outer: numberOrNull(element.getAttribute(RING_RADIUS_ATTRIBUTES.outer))
  }));
  return rings.every((ring) => ring.outer !== null) ? rings.sort((a, b) => Number(b.outer) - Number(a.outer)) : rings;
};

/**
 * A ring's hit at a drawn position: of its hits whose drawn span holds it, the one of highest identity (the one drawn
 * on top), with that identity and the colour it is drawn in; null when no hit does or Python wrote no spans.
 * @param {Element} ring
 * @param {number} position Bases from display index 0.
 * @returns {{ identity: number, color: string } | null}
 */
export const ringHitAt = (ring, position) => {
  /** @type {{ identity: number, color: string } | null} */
  let best = null;
  for (const hit of ring.querySelectorAll(`path[${HIT_DISPLAY_SPAN_ATTRIBUTES.start}]`)) {
    const start = numberOrNull(hit.getAttribute(HIT_DISPLAY_SPAN_ATTRIBUTES.start));
    const end = numberOrNull(hit.getAttribute(HIT_DISPLAY_SPAN_ATTRIBUTES.end));
    const identity = numberOrNull(hit.getAttribute('data-identity'));
    if (start === null || end === null || identity === null || position < start || position >= end) continue;
    if (best === null || identity > best.identity) best = { identity, color: hit.getAttribute('fill') ?? '' };
  }
  return best;
};

/** @param {number | null} identity */
export const identityText = (identity) => (identity === null ? t('No hit') : `${identity.toFixed(1)}%`);

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
 * Where the pointer is on a hovered ring's record: the source base (R5), the drawn position, and the radial line's
 * ends in overlay coordinates. Null when the Result has no axis or no spans for the record (a Result drawn before
 * Python wrote them: the rings still dim, without line or popup).
 * @param {{ root: Element, ring: Element, rings: DrawnRing[], axisRecords: RegionAxisRecord[], clientX: number, clientY: number, host: DOMRect }} input
 */
export const ringPointAt = ({ root, ring, rings, axisRecords, clientX, clientY, host }) => {
  const record = recordOfRing(ring, root);
  const axis = axisRecords.find((entry) => entry.recordIndex === record.index)?.axis;
  const frame = regionAxisFrame(root, record.index);
  const mapped = frame && axis?.kind === 'circular' ? framePoint(frame, clientX, clientY) : null;
  const matrix = frame?.getScreenCTM();
  if (!mapped || !matrix || axis?.kind !== 'circular' || rings.some((entry) => entry.inner === null)) return null;
  const position = axisPosition(axis, mapped.point).position;
  const angle = Math.atan2(mapped.point.y, mapped.point.x);
  const radii = [...rings.flatMap((entry) => [Number(entry.inner), Number(entry.outer)]), /** @type {CircularAxisGeometry} */ (axis).bandOuterPx];
  const end = (/** @type {number} */ radius) => new DOMPoint(radius * Math.cos(angle), radius * Math.sin(angle)).matrixTransform(matrix);
  const [from, to] = [end(Math.min(...radii)), end(Math.max(...radii))];
  return {
    base: sourceBaseAt(axis.source, Math.floor(position)),
    position,
    line: { x1: from.x - host.left, y1: from.y - host.top, x2: to.x - host.left, y2: to.y - host.top }
  };
};

/**
 * What a pointer event is over: a ring, else the Legend row whose box holds the point.
 * @param {Element} root @param {EventTarget | null} target @param {number} clientX @param {number} clientY
 * @returns {RingHoverTarget}
 */
export const ringHoverTarget = (root, target, clientX, clientY) => {
  const ring = target instanceof Element ? target.closest(RING_SELECTOR) : null;
  if (ring && root.contains(ring)) return { kind: 'ring', sourceIndex: ring.getAttribute('data-source-index') ?? '', ring };
  for (const row of root.querySelectorAll(RING_LEGEND_ROW_SELECTOR)) {
    const box = row.getBoundingClientRect();
    if (clientX >= box.left && clientX <= box.right && clientY >= box.top && clientY <= box.bottom) {
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
    name.className = 'gbdraw-ring-hover-name';
    name.textContent = ring.name;
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
 * @param {PopupCells} cells @param {{ identity: number, color: string } | null} hit
 */
const writeRingHit = ({ swatch, value }, hit) => {
  swatch.style.background = hit?.color ?? 'transparent';
  value.textContent = identityText(hit?.identity ?? null);
  value.toggleAttribute('data-no-hit', hit === null);
};

/**
 * The overlay beside the Result: the <style> of the dimming, the radial line, and the popup.
 * @param {HTMLElement} host
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
  const overlay = createRingHoverOverlay(host);
  /** @type {MouseEvent | null} */
  let pointer = null;
  let frame = 0;
  let shownStyle = '';
  /** @type {{ rings: Element[], sourceIndex: string, cells: PopupCells[] }} */
  let shownRows = { rings: [], sourceIndex: '', cells: [] };
  /** @param {string} css */
  const showStyle = (css) => {
    if (css !== shownStyle) overlay.style.textContent = shownStyle = css;
  };
  const clear = () => {
    showStyle('');
    overlay.element.hidden = true;
  };
  const update = () => {
    frame = 0;
    const event = pointer;
    const target = event && !held() && !isPreviewTransformActive()
      ? ringHoverTarget(svg, event.target, event.clientX, event.clientY) : null;
    if (!event || !target) return clear();
    showStyle(ringHoverStyle(target.sourceIndex));
    const box = host.getBoundingClientRect();
    const rings = target.kind === 'ring' ? ringsOfRecord(recordOfRing(target.ring, svg).group) : [];
    const point = target.kind === 'ring'
      ? ringPointAt({ root: svg, ring: target.ring, rings, axisRecords: axisRecords(), clientX: event.clientX, clientY: event.clientY, host: box })
      : null;
    overlay.element.hidden = !point;
    if (!point) return;
    Object.entries(point.line).forEach(([name, value]) => overlay.line.setAttribute(name, value.toFixed(1)));
    overlay.heading.textContent = t('{position} bp', { position: formatGroupedNumber(point.base) });
    if (shownRows.sourceIndex !== target.sourceIndex || shownRows.rings.length !== rings.length
      || shownRows.rings.some((element, index) => element !== rings[index].element)) {
      shownRows = { rings: rings.map((ring) => ring.element), sourceIndex: target.sourceIndex, cells: buildPopupRows(overlay.list, rings, target.sourceIndex) };
    }
    rings.forEach((ring, index) => writeRingHit(shownRows.cells[index], ringHitAt(ring.element, point.position)));
    const { offsetWidth: width, offsetHeight: height } = overlay.popup;
    const x = event.clientX - box.left + 16;
    const y = event.clientY - box.top + 16;
    overlay.popup.style.left = `${Math.max(4, x + width > box.width - 4 ? x - width - 32 : x)}px`;
    overlay.popup.style.top = `${Math.max(4, Math.min(y, box.height - height - 4))}px`;
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
    overlay.style.remove();
    overlay.element.remove();
  };
};
