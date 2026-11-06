// @ts-check
// The outline of Sequence search hits (PD-OI-117): a row click or a chip of
// the hit table outlines hits on the drawing. The outline is a preview over
// the Result, drawn in the overlay's own SVG and never in the Result, so it
// is no History step, is not saved, and is not exported. It follows each hit
// as Python drew it: the marks of the drawn hit, found by its set id and
// annotation id, are copied with their screen matrix; the Web computes no
// lane or band.
import { ANNOTATION_ID_ATTRIBUTE, ANNOTATION_SET_ID_ATTRIBUTE, WRAP_ROW_ATTRIBUTE } from '../../services/drawn-attributes.generated.js';
import { shallowRef } from '../../utils/reactivity.js';

/** The hits to outline: the set they belong to and their annotation ids. @typedef {{ setId: string, annotationIds: readonly string[] }} OutlineTarget */
/** One outlined mark: its path in its own coordinates and the matrix that places it over the overlay. @typedef {{ d: string, transform: string }} OutlineShape */

/**
 * The marks of one drawn hit: its paths and rectangles, directly in the hit's group or in one of its wrapped rows. The
 * steps of a gradient fill sit in a group of their own and are left out, so the outline is the mark's edge.
 * @param {Element} hit
 * @returns {SVGGraphicsElement[]}
 */
const marksOf = (hit) => Array.from(hit.querySelectorAll('path, rect')).flatMap((mark) => {
  const parent = mark.parentElement;
  const own = parent === hit || (parent !== null && parent.hasAttribute(WRAP_ROW_ATTRIBUTE) && parent.parentElement === hit);
  return own && mark instanceof SVGGraphicsElement ? [mark] : [];
});

/** @param {Element} mark The path data of a drawn mark (a rectangle as its four sides). */
const pathOf = (mark) => {
  if (mark.tagName.toLowerCase() !== 'rect') return mark.getAttribute('d') ?? '';
  const [x, y, width, height] = ['x', 'y', 'width', 'height'].map((name) => Number(mark.getAttribute(name) ?? 0));
  return `M${x} ${y}H${x + width}V${y + height}H${x}Z`;
};

/**
 * The ports of the outline: the mounted root of the shown Result, and the hits to outline (the search panel's).
 * @typedef {{ shownRoot: () => SVGSVGElement | null, target: () => OutlineTarget | null }} SearchHitOutlinePorts
 */

/** @param {SearchHitOutlinePorts} ports */
export const createSearchHitOutline = ({ shownRoot, target }) => {
  /** @type {{ value: OutlineShape[] }} */
  const shapes = shallowRef([]);
  /** @type {Element | null} */
  let overlayElement = null;
  let frameRequest = 0;
  // The marks found for a target in a root, found again when either changes or a mark leaves the page (a redraw).
  /** @type {{ root: Element | null, key: string, marks: SVGGraphicsElement[] }} */
  let found = { root: null, key: '', marks: [] };

  /** @param {SVGSVGElement} root @param {OutlineTarget} wanted */
  const marksFor = (root, wanted) => {
    const key = JSON.stringify([wanted.setId, wanted.annotationIds]);
    if (found.root !== root || found.key !== key || found.marks.some((mark) => !mark.isConnected)) {
      const ids = new Set(wanted.annotationIds);
      const hits = Array.from(root.querySelectorAll(`[${ANNOTATION_SET_ID_ATTRIBUTE}][${ANNOTATION_ID_ATTRIBUTE}]`))
        .filter((hit) => hit.getAttribute(ANNOTATION_SET_ID_ATTRIBUTE) === wanted.setId && ids.has(hit.getAttribute(ANNOTATION_ID_ATTRIBUTE) ?? ''));
      found = { root, key, marks: hits.flatMap(marksOf) };
    }
    return found.marks;
  };

  // The marks move with a pan, a zoom, or a Layout edit move: their screen matrices are read again on each frame
  // while the overlay is mounted.
  const followFrame = () => {
    frameRequest = 0;
    if (!overlayElement) return;
    const root = shownRoot();
    const wanted = target();
    const box = overlayElement.getBoundingClientRect();
    const next = root && wanted ? marksFor(root, wanted).flatMap((mark) => {
      const matrix = mark.getScreenCTM();
      return matrix ? [{
        d: pathOf(mark),
        transform: `matrix(${matrix.a} ${matrix.b} ${matrix.c} ${matrix.d} ${matrix.e - box.left} ${matrix.f - box.top})`
      }] : [];
    }) : [];
    if (JSON.stringify(next) !== JSON.stringify(shapes.value)) shapes.value = next;
    frameRequest = window.requestAnimationFrame(followFrame);
  };

  /** The overlay SVG, bound by the template while there are hits to outline. @param {Element | null} element */
  const bindOverlay = (element) => {
    overlayElement = element;
    if (!element) shapes.value = [];
    if (element && !frameRequest) frameRequest = window.requestAnimationFrame(followFrame);
  };

  return { shapes, bindOverlay };
};
