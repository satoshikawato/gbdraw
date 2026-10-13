// @ts-check
import { FONT_FAMILY_RULES } from '../utils/font-family-rules.generated.js';

// Python's font-family rules (gbdraw/core/fonts.py), ported for the decisions
// the browser makes before Python runs: which family of a CSS list a drawing
// uses and the nearest substitute when a named family is unavailable
// (`resolve_font_family_among`: the bundled fonts plus the fonts registered
// in this page; the browser has no fontconfig). The rule data is generated
// from Python (utils/font-family-rules.generated.js), and the vectors of
// tests/fixtures/font_family_resolution_vectors.json pin the evaluation to
// Python in both pytest and node.

/** A CSS `font-family` list as gbdraw stores it ("'Name', generic"). @typedef {string} FontFamilyList */
/** One unquoted family name of a list ("DejaVu Serif"). @typedef {string} FontFamilyName */
/** `serif`, `sans-serif`, or `monospace`, as Python's `font_generic_family` reports it. @typedef {string} GenericFontFamily */

/**
 * @typedef {object} FontFamilyResolution
 * @property {FontFamilyList} requested
 * @property {FontFamilyList} effective The list a drawing uses: only families that resolve.
 * @property {readonly FontFamilyName[]} missing The requested families that do not resolve, in order.
 */

const GENERIC_FAMILIES = new Set(FONT_FAMILY_RULES.genericFamilies);

/** Python's `_normalize_family_name`. @param {string} value @returns {FontFamilyName} */
export const normalizeFamilyName = (value) => value.trim().replace(/^["']+|["']+$/g, '').trim();

/** Python's `_comparable_name`. @param {string} value */
export const comparableFamilyName = (value) => normalizeFamilyName(value).toLowerCase().replaceAll(' ', '');

/** Python's `_is_generic`. @param {string} family */
export const isGenericFamily = (family) => GENERIC_FAMILIES.has(normalizeFamilyName(family).toLowerCase());

/** Python's `parse_font_family_list`. @param {unknown} fontFamily @returns {FontFamilyName[]} */
export const parseFontFamilyList = (fontFamily) => String(fontFamily ?? '')
  .split(',')
  .map(normalizeFamilyName)
  .filter(Boolean);

/** Python's `_css_family_token`. @param {string} family */
const cssFamilyToken = (family) => {
  const name = normalizeFamilyName(family);
  if (GENERIC_FAMILIES.has(name.toLowerCase())) return name.toLowerCase();
  const quote = name.includes("'") ? '"' : "'";
  return `${quote}${name}${quote}`;
};

/**
 * The bundled face prefix a family name is an alias of ("LiberationSerif"), or null.
 * @param {string} family
 * @returns {string | null}
 */
export const bundledFacePrefix = (family) => FONT_FAMILY_RULES.bundledFamilyAliases[comparableFamilyName(family)] ?? null;

/**
 * Python's `custom_font_families_among`: the named families of a list that
 * gbdraw does not bundle, or that are loaded in this page (a loaded "Times New
 * Roman" draws with its own file here and with Liberation elsewhere; D-10). A
 * Session stores only their names (D-02), so the note under Font and the Save
 * Session notice name them. The one rule of "custom font".
 * @param {unknown} fontFamily
 * @param {Iterable<FontFamilyName>} registeredFamilies The families loaded in this page.
 * @returns {readonly FontFamilyName[]}
 */
export const customFontFamilies = (fontFamily, registeredFamilies) => {
  const registered = new Set(Array.from(registeredFamilies, comparableFamilyName));
  return Object.freeze(parseFontFamilyList(fontFamily).filter((family) => (
    !isGenericFamily(family) && (bundledFacePrefix(family) === null || registered.has(comparableFamilyName(family)))
  )));
};

/** The bundled face prefix of a list that resolves to nothing. */
export const FALLBACK_FACE_PREFIX = FONT_FAMILY_RULES.fallbackFace;

/**
 * Python's `_nearest_generic`.
 * @param {readonly string[]} generics
 * @param {readonly string[]} missing
 * @returns {GenericFontFamily}
 */
const nearestGeneric = (generics, missing) => {
  for (const generic of generics) {
    if (Object.hasOwn(FONT_FAMILY_RULES.liberationFamilyByGeneric, generic.toLowerCase())) return generic.toLowerCase();
  }
  const names = missing.join(' ').toLowerCase();
  if (FONT_FAMILY_RULES.monospaceNameHints.some((hint) => names.includes(hint))) return 'monospace';
  if ((names.includes('serif') && !names.includes('sans'))
    || FONT_FAMILY_RULES.serifNameHints.some((hint) => names.includes(hint))) return 'serif';
  return 'sans-serif';
};

/**
 * Python's `_resolve_font_family_with`: keep the generics and the named
 * families `available` accepts; when none of the named ones is left, insert the
 * nearest Liberation family before the first generic.
 * @param {FontFamilyList} fontFamily
 * @param {(family: FontFamilyName) => boolean} available
 * @returns {FontFamilyResolution}
 */
export const resolveFontFamilyWith = (fontFamily, available) => {
  const requested = String(fontFamily ?? '');
  /** @type {FontFamilyName[]} */
  const kept = [];
  /** @type {FontFamilyName[]} */
  const missing = [];
  for (const family of parseFontFamilyList(requested)) {
    if (isGenericFamily(family) || available(family)) kept.push(family);
    else missing.push(family);
  }
  if (missing.length === 0) return Object.freeze({ requested, effective: requested, missing: Object.freeze([]) });
  const generics = kept.filter(isGenericFamily);
  if (generics.length === kept.length) {
    const generic = nearestGeneric(generics, missing);
    if (generics.length === 0) kept.push(generic);
    const position = kept.findIndex(isGenericFamily);
    kept.splice(position, 0, FONT_FAMILY_RULES.liberationFamilyByGeneric[generic]);
  }
  return Object.freeze({
    requested,
    effective: kept.map(cssFamilyToken).join(', '),
    missing: Object.freeze(missing)
  });
};

/**
 * Python's `resolve_font_family_among`: what the Web app resolves, with the
 * bundled fonts and the families registered in this page available.
 * @param {FontFamilyList} fontFamily
 * @param {Iterable<FontFamilyName>} registeredFamilies
 * @returns {FontFamilyResolution}
 */
export const resolveWebFontFamily = (fontFamily, registeredFamilies) => {
  const registered = new Set(Array.from(registeredFamilies, comparableFamilyName));
  return resolveFontFamilyWith(
    fontFamily,
    (family) => bundledFacePrefix(family) !== null || registered.has(comparableFamilyName(family))
  );
};
