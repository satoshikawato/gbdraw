import { sanitizeSvgContent } from './svg-sanitization.js';
import { serializeCleanSvg } from './svg-serialization.js';
import { collectRenderedFeatureIdentitiesFromSvgRoot } from './session-feature-metadata.js';
import { normalizeSvgResultIds } from './svg-result-normalization.js';

// This symbol is intentionally module-private. JSON/session data cannot reproduce it,
// while ordinary in-memory Result object spreads retain the runtime ownership record.
const COMMITTED_SVG_RESULT = Symbol('gbdraw.committedSvgResult');

let nextResultIdentity = 1;

const committedState = (result) => (
  result && typeof result === 'object' && result[COMMITTED_SVG_RESULT]
    ? result[COMMITTED_SVG_RESULT]
    : null
);

const parseSanitizedSvg = (
  content,
  parser = globalThis.DOMParser || globalThis.window?.DOMParser
) => {
  if (typeof parser !== 'function') {
    throw new Error('SVG parsing is unavailable.');
  }
  const document = new parser().parseFromString(content, 'image/svg+xml');
  if (
    document?.querySelector?.('parsererror')
    || String(document?.documentElement?.localName || '').toLowerCase() !== 'svg'
  ) {
    throw new Error('The diagram engine returned malformed SVG content.');
  }
  return document.documentElement;
};

const markCommitted = (result, metadata) => {
  result[COMMITTED_SVG_RESULT] = {
    identity: nextResultIdentity++,
    mounted: false,
    mountedContent: result.content,
    metadata
  };
  return result;
};

const admitSvgRoot = (result, svg) => {
  normalizeSvgResultIds(svg);
  const metadata = Object.freeze({
    renderedFeatureIdentities: collectRenderedFeatureIdentitiesFromSvgRoot(svg)
  });
  return markCommitted({
    ...result,
    content: serializeCleanSvg(svg)
  }, metadata);
};

export const parseCommittedSvgResultRoot = (
  result,
  parser = globalThis.DOMParser || globalThis.window?.DOMParser
) => {
  if (!isCommittedSvgResult(result)) {
    throw new Error('Only an admitted SVG Result can be parsed for projection.');
  }
  return parseSanitizedSvg(getCommittedSvgContent(result), parser);
};

export const admitProjectedSvgResult = (result, svg) => {
  if (!isCommittedSvgResult(result)) {
    throw new Error('Only an admitted SVG Result can receive a projected root.');
  }
  if (!svg || String(svg.localName || '').toLowerCase() !== 'svg') {
    throw new Error('A projected SVG root is required.');
  }
  return admitSvgRoot(result, svg);
};

export const isCommittedSvgResult = (result) => Boolean(committedState(result));

export const isCommittedSvgResultMounted = (result) => Boolean(
  committedState(result)?.mounted
);

export const getCommittedSvgContent = (result) => {
  const runtime = committedState(result);
  if (!runtime) return null;
  return runtime.mounted ? runtime.mountedContent : result.content;
};

export const getCommittedSvgResultMetadata = (result) => (
  committedState(result)?.metadata || null
);

export const markCommittedSvgResultMounted = (result) => {
  const runtime = committedState(result);
  if (!runtime) return false;
  runtime.mounted = true;
  runtime.mountedContent = result.content;
  return true;
};

export const markCommittedSvgResultUnmounted = (result) => {
  const runtime = committedState(result);
  if (!runtime) return false;
  runtime.mounted = false;
  runtime.mountedContent = result.content;
  return true;
};

/**
 * Admit one untrusted SVG Result into the application.
 *
 * The raw string is sanitized before parsing. A transform may inspect or mutate the
 * detached sanitized root. The detached root is serialized exactly once so the
 * committed Result already uses the same canonical representation as later exports.
 * The Result receives runtime-only committed state only after the transaction succeeds.
 */
export const ingestSvgResult = (
  result,
  {
    sanitizer = globalThis.DOMPurify || globalThis.window?.DOMPurify,
    parser = globalThis.DOMParser || globalThis.window?.DOMParser,
    transformSvg = null,
    resultIndex = 0
  } = {}
) => {
  if (!result || typeof result !== 'object' || Array.isArray(result)) {
    throw new Error('The diagram engine returned an invalid SVG Result.');
  }

  const sanitized = sanitizeSvgContent(result.content, sanitizer);
  const svg = parseSanitizedSvg(sanitized, parser);
  if (typeof transformSvg === 'function') {
    transformSvg(svg, { result, resultIndex });
  }
  return admitSvgRoot(result, svg);
};

/**
 * Reproject an already admitted Result without running untrusted-input sanitization again.
 * The detached root and replacement admission record are published only on success.
 */
export const reprojectCommittedSvgResult = (
  result,
  {
    parser = globalThis.DOMParser || globalThis.window?.DOMParser,
    transformSvg = null,
    resultIndex = 0
  } = {}
) => {
  if (!isCommittedSvgResult(result)) {
    throw new Error('Only an admitted SVG Result can be reprojected.');
  }
  const svg = parseCommittedSvgResultRoot(result, parser);
  if (typeof transformSvg === 'function') {
    transformSvg(svg, { result, resultIndex });
  }
  return admitProjectedSvgResult(result, svg);
};

/**
 * Admit a candidate derived from the mounted SVG without parsing or sanitizing it.
 * The supplied root is cloned so preparation cannot mutate the live mount.
 */
export const projectMountedSvgResult = (
  result,
  mountedSvg,
  { transformSvg = null, resultIndex = 0 } = {}
) => {
  if (!isCommittedSvgResult(result)) {
    throw new Error('Only an admitted SVG Result can be projected from a mount.');
  }
  if (!mountedSvg?.cloneNode) {
    throw new Error('A mounted SVG root is required for projection.');
  }
  const svg = mountedSvg.cloneNode(true);
  if (typeof transformSvg === 'function') {
    transformSvg(svg, { result, resultIndex });
  }
  return admitProjectedSvgResult(result, svg);
};

export const ingestSvgResults = (results, options = {}) => {
  if (!Array.isArray(results)) {
    throw new Error('The diagram engine returned an invalid Result list.');
  }
  return results.map((result, resultIndex) => (
    isCommittedSvgResult(result)
      ? result
      : ingestSvgResult(result, { ...options, resultIndex })
  ));
};
