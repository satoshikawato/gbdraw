// @ts-check
// Runtime-only provenance shared by the Worker response decoder and SVG admission.
// The token is non-enumerable and cannot survive a JSON/Session round trip.
const CURRENT_WORKER_GENERATION_RESPONSE = Symbol('gbdraw.currentWorkerGenerationResponse');

/** @internal Mint provenance only while decoding a response from the current Worker path. */
/**
 * @template T
 * @param {T} response
 * @returns {T}
 */
export const markCurrentWorkerGenerationResponse = (response) => {
  if (!response || typeof response !== 'object' || Array.isArray(response)) {
    throw new TypeError('Current Worker response provenance requires a response object.');
  }
  Object.defineProperty(response, CURRENT_WORKER_GENERATION_RESPONSE, {
    value: true,
    enumerable: false
  });
  return response;
};

/** Test runtime-only Worker provenance without exposing its private Symbol. */
/** @param {{ [CURRENT_WORKER_GENERATION_RESPONSE]?: unknown } | null | undefined} value */
export const isCurrentWorkerGenerationResponse = (value) => Boolean(
  value?.[CURRENT_WORKER_GENERATION_RESPONSE]
);

/**
 * A children text of a Generate reply that the page sanitized once for its
 * revision (the sha256 of the UTF-8 text Python wrote; `bytes` is that
 * length). Only DOMPurify output is kept, never the raw text (CW-04).
 * @typedef {object} HeldSegment
 * @property {string} revision
 * @property {number} bytes
 * @property {string} sanitized
 */

/**
 * A placed content group of a decoded Generate reply: `children` is the
 * text the reply sent; a piece without it names a segment the page holds.
 * @typedef {{ kind: 'content', start: string, revision: string, bytes: number, children?: string }} DecodedContentPiece
 */

/**
 * A Generate Result's document as decoded, before admission: the root start
 * tag, the top-level pieces in order, and the root end tag.
 * @typedef {object} DecodedSvgDocument
 * @property {string} head
 * @property {'</svg>'} tail
 * @property {readonly (Readonly<DecodedContentPiece> | Readonly<{ kind: 'raw', text: string }>)[]} pieces
 */

/**
 * The key of a Result's held segments, `Map<revision, HeldSegment>`: on a
 * decoded candidate, the segments of the request's snapshot its document
 * names; on a committed Result, one segment for each revision of its document.
 * An enumerable own property, so spreads of the Result carry it; JSON,
 * structuredClone and Session save skip it, and it dies with the object.
 */
export const SEGMENTS = Symbol('gbdraw.svgResultSegments');
