// Runtime-only provenance shared by the Worker response decoder and SVG admission.
// The token is non-enumerable and cannot survive a JSON/Session round trip.
const CURRENT_WORKER_GENERATION_RESPONSE = Symbol('gbdraw.currentWorkerGenerationResponse');

/** @internal Mint provenance only while decoding a response from the current Worker path. */
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
export const isCurrentWorkerGenerationResponse = (value) => Boolean(
  value?.[CURRENT_WORKER_GENERATION_RESPONSE]
);
