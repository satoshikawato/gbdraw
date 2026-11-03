// @ts-check
// A file name the app derives from a record ID or selector (FL-10). It is the
// name a browser keeps when it downloads the file, so a Result's name and its
// downloaded file's name are the same: the characters Chromium's download-name
// rule replaces anywhere (`"~*/:<>?\|`, control and format characters) become
// `_`, and the dots and white space it strips from both ends go. A Windows
// device name (CON, NUL, COM1, ...) takes a leading `_`, so Python accepts the
// name as one portable output prefix (gbdraw.api.requests.RenderOutputRequest).

const ILLEGAL = /["~*/:<>?\\|\p{Cc}\p{Cf}]/gu;
const ENDS = /^[\s.]+|[\s.]+$/gu;
const DEVICE_NAMES = new Set([
  'CON', 'PRN', 'AUX', 'NUL', 'CONIN$', 'CONOUT$',
  ...[...'123456789¹²³'].flatMap((suffix) => [`COM${suffix}`, `LPT${suffix}`])
]);

// The output prefix of a Linear request without one; the Output Prefix field shows it as its placeholder.
// Python owns the value (`RenderOutputRequest.output_prefix`); this module is the Web side of the registered
// twin (`concern:output-prefix`).
export const DEFAULT_LINEAR_OUTPUT_PREFIX = 'out';

/**
 * @param {unknown} value
 * @param {string} fallback The name when nothing is left.
 * @returns {string}
 */
export const downloadSafeName = (value, fallback) => {
  const name = String(value ?? '').replace(ILLEGAL, '_').replace(ENDS, '');
  if (!name) return fallback;
  return DEVICE_NAMES.has(name.split('.', 1)[0].replace(/ +$/, '').toUpperCase()) ? `_${name}` : name;
};
