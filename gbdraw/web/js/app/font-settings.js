// @ts-check
/** @import { WebFontEntry, WebFontSource, RegisterFontsInWorker } from '../services/font-store.js' */
/** @import { FontFamilyList } from '../services/font-families.js' */
import {
  addWebFontsToDocument,
  mergeWebFonts,
  registerWebFonts,
  webFontFamilies
} from '../services/font-store.js';
import { bundledFacePrefix, comparableFamilyName, parseFontFamilyList } from '../services/font-families.js';
import { normalizeUserFacingError } from '../utils/error-normalization.js';

const { ref, reactive, computed } = window.Vue;

// The Font section (custom fonts D-01, D-06): the drawing's text font
// (`adv.font_family`, Applies on Generate), the fonts loaded in this page from
// a font file (A) or from the installed fonts (B, Local Font Access), and the
// note for a font a Session cannot carry (D-02).

/** The select's value of the default font (Liberation Sans, `adv.font_family` unset). */
export const DEFAULT_FONT_CHOICE = 'default';
/** The select's value that opens the font file picker. */
export const LOAD_FONT_FILE_CHOICE = 'load-file';
/** The select's value that opens the installed-font picker. */
export const INSTALLED_FONT_CHOICE = 'installed';

/** @type {ReadonlyArray<{ family: string, fontFamily: FontFamilyList | null }>} */
const BUNDLED_CHOICES = Object.freeze([
  { family: 'Liberation Sans', fontFamily: null },
  { family: 'Liberation Serif', fontFamily: "'Liberation Serif', serif" },
  { family: 'Liberation Mono', fontFamily: "'Liberation Mono', monospace" }
]);

/**
 * One option of the Font select.
 * @typedef {object} FontChoice
 * @property {string} value
 * @property {string} label
 * @property {'bundled' | 'loaded' | 'unavailable' | 'action'} kind
 */

/**
 * A face Local Font Access lists (`queryLocalFonts()`).
 * @typedef {object} InstalledFontData
 * @property {string} family
 * @property {string} fullName
 * @property {string} style
 * @property {() => Promise<Blob>} blob
 */

/**
 * What the Font section reads and writes of the app state.
 * @typedef {object} FontSettingsState
 * @property {{ value: WebFontEntry[] }} loadedFonts
 * @property {() => { adv: Record<string, unknown> }} activeDrawing
 */

/**
 * @typedef {object} FontSettingsOptions
 * @property {FontSettingsState} state App state (state.js): `loadedFonts` and the active drawing's `adv`.
 * @property {(label: string, fn: () => unknown) => Promise<unknown>} runUndoable History's one step (R11).
 * @property {RegisterFontsInWorker} registerFontsInWorker The Worker's `registerFonts` helper.
 */

/** @param {FontSettingsOptions} options */
export const createFontSettings = ({ state, runUndoable, registerFontsInWorker }) => {
  const loadError = ref('');
  const loading = ref(false);
  /** @type {{ open: boolean, query: string, families: { family: string, faces: InstalledFontData[] }[], error: string }} */
  const installedPicker = reactive({ open: false, query: '', families: [], error: '' });
  const localFontsSupported = typeof /** @type {{ queryLocalFonts?: unknown }} */ (globalThis).queryLocalFonts === 'function';

  /** @returns {readonly WebFontEntry[]} */
  const loadedFonts = () => state.loadedFonts.value;
  /** @returns {FontFamilyList | null} */
  const draftFontFamily = () => {
    const value = state.activeDrawing().adv.font_family;
    return typeof value === 'string' && value.trim() ? value : null;
  };

  /** The select value of a stored list: its first named family. @param {FontFamilyList | null} fontFamily */
  const choiceValueOf = (fontFamily) => {
    if (fontFamily === null) return DEFAULT_FONT_CHOICE;
    const first = parseFontFamilyList(fontFamily)[0] || '';
    const bundled = BUNDLED_CHOICES.find((choice) => comparableFamilyName(choice.family) === comparableFamilyName(first));
    if (bundled) return bundled.fontFamily === null ? DEFAULT_FONT_CHOICE : bundled.fontFamily;
    const loaded = webFontFamilies(loadedFonts()).find((entry) => comparableFamilyName(entry.family) === comparableFamilyName(first));
    return loaded ? loaded.fontFamily : fontFamily;
  };

  const choices = computed(() => {
    /** @type {FontChoice[]} */
    const list = BUNDLED_CHOICES.map((choice) => ({
      value: choice.fontFamily === null ? DEFAULT_FONT_CHOICE : choice.fontFamily,
      label: choice.fontFamily === null ? `${choice.family} (default)` : choice.family,
      kind: /** @type {const} */ ('bundled')
    }));
    webFontFamilies(loadedFonts()).forEach((entry) => {
      list.push({ value: entry.fontFamily, label: entry.family, kind: 'loaded' });
    });
    const current = choiceValueOf(draftFontFamily());
    if (!list.some((choice) => choice.value === current)) {
      const first = parseFontFamilyList(current)[0] || current;
      list.push({ value: current, label: `${first} (not loaded)`, kind: 'unavailable' });
    }
    list.push({ value: LOAD_FONT_FILE_CHOICE, label: 'Load font file…', kind: 'action' });
    if (localFontsSupported) list.push({ value: INSTALLED_FONT_CHOICE, label: 'Choose installed font…', kind: 'action' });
    return list;
  });

  const selectedChoice = computed({
    get: () => choiceValueOf(draftFontFamily()),
    /** @param {string} value */
    set: (value) => {
      if (value === LOAD_FONT_FILE_CHOICE || value === INSTALLED_FONT_CHOICE) return;
      state.activeDrawing().adv.font_family = value === DEFAULT_FONT_CHOICE ? null : value;
    }
  });

  /**
   * Whether the drawing's font is one a Session cannot carry: its first named
   * family is not bundled (D-02 note).
   */
  const customFontSelected = computed(() => {
    const first = parseFontFamilyList(draftFontFamily() || '')[0];
    return Boolean(first) && bundledFacePrefix(first) === null;
  });

  /**
   * Register font sources and select the family of the first one, as one History step.
   * @param {readonly WebFontSource[]} sources
   */
  const registerAndSelect = async (sources) => {
    if (sources.length === 0) return;
    loadError.value = '';
    loading.value = true;
    try {
      const drawing = state.activeDrawing();
      await runUndoable('Change font', async () => {
        const added = await registerWebFonts(sources, registerFontsInWorker);
        await addWebFontsToDocument(added);
        state.loadedFonts.value = mergeWebFonts(loadedFonts(), added);
        drawing.adv.font_family = added[0].fontFamily;
      });
    } catch (error) {
      loadError.value = normalizeUserFacingError(error, { code: 'INPUT_INVALID', operation: 'registerFonts' })?.summary || 'The font could not be loaded.';
    } finally {
      loading.value = false;
    }
  };

  /**
   * Load font files (A): the `change` of the Font file input.
   * @param {Event} event
   */
  const loadFontFiles = async (event) => {
    const input = /** @type {HTMLInputElement} */ (event.target);
    const files = Array.from(input.files || []);
    input.value = '';
    /** @type {WebFontSource[]} */
    const sources = [];
    for (const file of files) {
      sources.push({ bytes: await file.arrayBuffer(), source: 'file', label: file.name });
    }
    await registerAndSelect(sources);
  };

  /** Open the installed-font picker (B); asks for the browser permission on this click. */
  const openInstalledFonts = async () => {
    installedPicker.error = '';
    installedPicker.query = '';
    const query = /** @type {{ queryLocalFonts?: () => Promise<InstalledFontData[]> }} */ (globalThis).queryLocalFonts;
    if (typeof query !== 'function') return;
    try {
      const faces = await query();
      /** @type {Map<string, InstalledFontData[]>} */
      const byFamily = new Map();
      faces.forEach((face) => byFamily.set(face.family, [...(byFamily.get(face.family) || []), face]));
      installedPicker.families = [...byFamily].map(([family, list]) => ({ family, faces: list }))
        .sort((left, right) => left.family.localeCompare(right.family));
      installedPicker.open = true;
    } catch (_error) {
      installedPicker.error = 'Installed fonts are not available: the browser did not allow access to local fonts.';
      installedPicker.open = true;
    }
  };

  const filteredInstalledFamilies = computed(() => {
    const query = installedPicker.query.trim().toLowerCase();
    return query
      ? installedPicker.families.filter((entry) => entry.family.toLowerCase().includes(query))
      : installedPicker.families;
  });

  /**
   * Register every face of an installed family (Python picks the face of each
   * weight and style) and select it.
   * @param {string} family
   */
  const chooseInstalledFamily = async (family) => {
    const entry = installedPicker.families.find((candidate) => candidate.family === family);
    if (!entry) return;
    installedPicker.open = false;
    // The regular face first, so the family's first registered face is upright.
    const faces = [...entry.faces].sort((left, right) => Number(right.style === 'Regular') - Number(left.style === 'Regular'));
    /** @type {WebFontSource[]} */
    const sources = [];
    for (const face of faces) {
      sources.push({ bytes: await (await face.blob()).arrayBuffer(), source: 'local', label: face.fullName });
    }
    await registerAndSelect(sources);
  };

  const closeInstalledFonts = () => {
    installedPicker.open = false;
  };

  return {
    choices,
    selectedChoice,
    customFontSelected,
    localFontsSupported,
    loading,
    loadError,
    loadFontFiles,
    installedPicker,
    filteredInstalledFamilies,
    openInstalledFonts,
    chooseInstalledFamily,
    closeInstalledFonts
  };
};
