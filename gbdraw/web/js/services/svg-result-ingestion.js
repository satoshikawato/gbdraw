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
  normalizeSvgResultIds(svg);
  const metadata = Object.freeze({
    renderedFeatureIdentities: collectRenderedFeatureIdentitiesFromSvgRoot(svg)
  });
  const content = serializeCleanSvg(svg);

  return markCommitted({
    ...result,
    content
  }, metadata);
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
