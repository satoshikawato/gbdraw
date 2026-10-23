// @ts-check
/** @import { WebFontEntry, WebFontSource, RegisterFontsInWorker } from '../services/font-store.js' */
/** @import { FontFamilyList, FontFamilyName } from '../services/font-families.js' */
/** @import { AskChoice } from './app-messages.js' */
import { loadWebFonts, mergeWebFonts, standardInstalledFaces, webFontFamilies, webFontFamilyNames } from '../services/font-store.js';
import {
  bundledFacePrefix,
  comparableFamilyName,
  customFontFamilies,
  isGenericFamily,
  parseFontFamilyList,
  resolveWebFontFamily
} from '../services/font-families.js';
import { SESSION_FONT_ACTIONS, missingSessionFonts, sessionFontQuestion } from '../services/session-fonts.js';
import { diagnosticError, normalizeUserFacingError } from '../utils/error-normalization.js';
import { msg, t } from '../utils/i18n.js';
import { errorText } from '../utils/error-text.js';

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
/** The font files Load font file… offers. */
const FONT_FILE_ACCEPT = '.ttf,.otf,font/ttf,font/otf';

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
 * @property {AskChoice} askChoice The app's choice dialog (the Session font question).
 * @property {() => unknown} committedFontFamily The text font of the shown Result's committed
 *   request (`requestFontFamily`), which its exports carry (OV-742).
 */

/** @param {FontSettingsOptions} options */
export const createFontSettings = ({ state, runUndoable, registerFontsInWorker, askChoice, committedFontFamily }) => {
  const loadError = ref('');
  const loading = ref(false);
  /** @type {{ open: boolean, query: string, families: { family: string, faces: InstalledFontData[] }[], error: string }} */
  const installedPicker = reactive({ open: false, query: '', families: [], error: '' });
  /** @type {((faces: InstalledFontData[] | null) => void) | null} */
  let answerInstalledPicker = null;
  const localFontsSupported = typeof /** @type {{ queryLocalFonts?: unknown }} */ (globalThis).queryLocalFonts === 'function';
  /**
   * Embed font in exported SVG (D-05): whether SVG, Interactive SVG, and PNG
   * exports carry the fonts gbdraw does not bundle. On by default; page state
   * like the loaded fonts, so a Session does not keep it and changing it is
   * no History step (exports record none).
   */
  const embedInExports = ref(true);
  /**
   * The fonts the diagram Worker lost (OV-686): they stay listed (D-10 (6)),
   * and loading the same bytes again registers them again (OV-740).
   * @type {Set<string>}
   */
  const lostIdentities = new Set();

  /** @returns {readonly WebFontEntry[]} */
  const loadedFonts = () => state.loadedFonts.value;
  /** @returns {FontFamilyName[]} */
  const loadedFamilies = () => webFontFamilyNames(loadedFonts());
  /** @returns {FontFamilyList | null} */
  const draftFontFamily = () => {
    const value = state.activeDrawing().adv.font_family;
    return typeof value === 'string' && value.trim() ? value : null;
  };

  /**
   * The select value of a stored list: the entry of its first named family,
   * a font loaded in this page before a bundled font of the same alias name
   * (OV-685: the user's explicit font wins).
   * @param {FontFamilyList | null} fontFamily
   */
  const choiceValueOf = (fontFamily) => {
    if (fontFamily === null) return DEFAULT_FONT_CHOICE;
    const first = parseFontFamilyList(fontFamily)[0] || '';
    const loaded = webFontFamilies(loadedFonts()).find((entry) => comparableFamilyName(entry.family) === comparableFamilyName(first));
    if (loaded) return loaded.fontFamily;
    // A bundled alias ("Arial") is the bundled family Python draws it with.
    const prefix = bundledFacePrefix(first);
    const bundled = prefix === null ? undefined : BUNDLED_CHOICES.find((choice) => bundledFacePrefix(choice.family) === prefix);
    if (bundled) return bundled.fontFamily === null ? DEFAULT_FONT_CHOICE : bundled.fontFamily;
    return fontFamily;
  };

  /**
   * The option of a font this page cannot draw: the Session's name, and the
   * nearest font a Generate draws it with (Python's substitution).
   * @param {FontFamilyList} fontFamily
   */
  const unavailableLabel = (fontFamily) => {
    const resolution = resolveWebFontFamily(fontFamily, loadedFamilies());
    const requested = resolution.replaced[0] || parseFontFamilyList(fontFamily)[0] || fontFamily;
    const nearest = parseFontFamilyList(resolution.effective).find((family) => !isGenericFamily(family)) || resolution.effective;
    return t('{requested} — not available, drawn in {nearest}', { requested, nearest });
  };

  const choices = computed(() => {
    /** @type {FontChoice[]} */
    const list = BUNDLED_CHOICES.map((choice) => ({
      value: choice.fontFamily === null ? DEFAULT_FONT_CHOICE : choice.fontFamily,
      label: choice.fontFamily === null ? t('{family} (default)', { family: choice.family }) : choice.family,
      kind: /** @type {const} */ ('bundled')
    }));
    webFontFamilies(loadedFonts()).forEach((entry) => {
      list.push({ value: entry.fontFamily, label: entry.family, kind: 'loaded' });
    });
    const current = choiceValueOf(draftFontFamily());
    if (!list.some((choice) => choice.value === current)) {
      list.push({ value: current, label: unavailableLabel(current), kind: 'unavailable' });
    }
    list.push({ value: LOAD_FONT_FILE_CHOICE, label: t('Load font file…'), kind: 'action' });
    if (localFontsSupported) list.push({ value: INSTALLED_FONT_CHOICE, label: t('Choose installed font…'), kind: 'action' });
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
   * Whether the drawing's font, or the font of the shown Result, names a
   * family a Session cannot carry (`customFontFamilies`, the rule of the Save
   * Session warning too; D-02 note). The Result's font counts because its
   * exports carry it until the next Generate (OV-742).
   */
  const customFontSelected = computed(() => [draftFontFamily(), committedFontFamily()]
    .some((fontFamily) => customFontFamilies(fontFamily, loadedFamilies()).length > 0));

  /**
   * Mark fonts the diagram Worker could not register again (OV-686), so
   * loading their bytes again registers them (OV-740).
   * @param {readonly string[]} identities
   */
  const markFontsLost = (identities) => {
    identities.forEach((identity) => lostIdentities.add(identity));
  };

  /**
   * Load font sources into the page: the browser checks them, Python
   * registers them, the preview gets them, and the store keeps them.
   * @param {readonly WebFontSource[]} sources
   * @returns {Promise<WebFontEntry[]>} The entry of each source, in order.
   */
  const loadSources = async (sources) => {
    const known = loadedFonts().filter((entry) => !lostIdentities.has(entry.identity));
    const loaded = await loadWebFonts(sources, { known, registerInWorker: registerFontsInWorker });
    loaded.added.forEach((entry) => lostIdentities.delete(entry.identity));
    state.loadedFonts.value = mergeWebFonts(loadedFonts(), loaded.added);
    return loaded.chosen;
  };

  /**
   * Read font sources; a read the browser refuses is a font-file error.
   * @param {() => Promise<WebFontSource[]>} read
   */
  const readSources = async (read) => {
    try {
      return await read();
    } catch (_error) {
      throw diagnosticError('INPUT_INVALID', { field: 'font_family', reason: 'FONT_FILE' }, { stage: 'read', operation: 'registerFonts' });
    }
  };

  /** @param {unknown} error A caught failure: a diagnostic or a Worker error. */
  const loadErrorSummary = (error) => {
    const model = normalizeUserFacingError(/** @type {Record<string, unknown>} */ (error), { code: 'INPUT_INVALID', operation: 'registerFonts' });
    return (model && errorText(model.summaryParts)) || t('The font could not be loaded.');
  };

  /**
   * Read and load fonts, then select the family of the first one as one
   * History step. `loading` holds from the read to the end.
   * @param {() => Promise<WebFontSource[]>} read
   */
  const loadAndSelect = async (read) => {
    loadError.value = '';
    loading.value = true;
    try {
      const sources = await readSources(read);
      if (sources.length === 0) return;
      const drawing = state.activeDrawing();
      await runUndoable(msg('Change font'), async () => {
        const [first] = await loadSources(sources);
        drawing.adv.font_family = first.fontFamily;
      });
    } catch (error) {
      loadError.value = loadErrorSummary(error);
    } finally {
      loading.value = false;
    }
  };

  /**
   * Let the user choose font files (A); resolves to none when the picker is closed.
   * @returns {Promise<File[]>}
   */
  const chooseFontFiles = () => new Promise((resolve) => {
    const input = document.createElement('input');
    input.type = 'file';
    input.accept = FONT_FILE_ACCEPT;
    input.hidden = true;
    input.dataset.fontFileInput = '';
    input.dataset.historyManaged = '';
    const settle = () => {
      const files = Array.from(input.files || []);
      input.remove();
      resolve(files);
    };
    input.addEventListener('change', settle, { once: true });
    input.addEventListener('cancel', settle, { once: true });
    document.body.append(input);
    input.click();
  });

  /** @param {readonly File[]} files @returns {Promise<WebFontSource[]>} */
  const fileSources = async (files) => {
    /** @type {WebFontSource[]} */
    const sources = [];
    for (const file of files) sources.push({ bytes: await file.arrayBuffer(), source: 'file', label: file.name });
    return sources;
  };

  /** @param {readonly InstalledFontData[]} faces @returns {Promise<WebFontSource[]>} */
  const installedSources = async (faces) => {
    /** @type {WebFontSource[]} */
    const sources = [];
    for (const face of standardInstalledFaces(faces)) {
      sources.push({ bytes: await (await face.blob()).arrayBuffer(), source: 'local', label: face.fullName });
    }
    return sources;
  };

  /**
   * Open the installed-font picker (B; asks for the browser permission on
   * this click) and resolve to the faces of the family chosen, or null.
   * @returns {Promise<InstalledFontData[] | null>}
   */
  const pickInstalledFamily = async () => {
    answerInstalledPicker?.(null);
    installedPicker.error = '';
    installedPicker.query = '';
    const query = /** @type {{ queryLocalFonts?: () => Promise<InstalledFontData[]> }} */ (globalThis).queryLocalFonts;
    if (typeof query !== 'function') return null;
    try {
      const faces = await query();
      /** @type {Map<string, InstalledFontData[]>} */
      const byFamily = new Map();
      faces.forEach((face) => byFamily.set(face.family, [...(byFamily.get(face.family) || []), face]));
      installedPicker.families = [...byFamily].map(([family, list]) => ({ family, faces: list }))
        .sort((left, right) => left.family.localeCompare(right.family));
    } catch (_error) {
      installedPicker.families = [];
      installedPicker.error = t('Installed fonts are not available: the browser did not allow access to local fonts.');
    }
    installedPicker.open = true;
    return new Promise((resolve) => { answerInstalledPicker = resolve; });
  };

  /** @param {InstalledFontData[] | null} faces */
  const settleInstalledPicker = (faces) => {
    const answer = answerInstalledPicker;
    answerInstalledPicker = null;
    installedPicker.open = false;
    answer?.(faces);
  };

  const filteredInstalledFamilies = computed(() => {
    const query = installedPicker.query.trim().toLowerCase();
    return query
      ? installedPicker.families.filter((entry) => entry.family.toLowerCase().includes(query))
      : installedPicker.families;
  });

  /**
   * The picker's choice of a family: its Regular, Bold, Italic, and Bold
   * Italic faces are registered (Python picks the face of each weight and style).
   * @param {string} family
   */
  const chooseInstalledFamily = (family) => {
    const entry = installedPicker.families.find((candidate) => candidate.family === family);
    if (entry) settleInstalledPicker(entry.faces);
  };

  const closeInstalledFonts = () => settleInstalledPicker(null);

  /** The select's Load font file… (A). */
  const loadFontFile = () => loadAndSelect(async () => fileSources(await chooseFontFiles()));

  /** The select's Choose installed font… (B). */
  const openInstalledFonts = async () => {
    const faces = await pickInstalledFamily();
    if (faces) await loadAndSelect(() => installedSources(faces));
  };

  /**
   * Before a Session replaces the current one: when a font it names is
   * neither bundled nor loaded, ask once (receipt 〔読み込むとき〕). Load font
   * file… and Use installed font… load a font and ask again until the
   * Session's fonts resolve; Load with <nearest> keeps the names, so Generate
   * draws the nearest font until the font is loaded. Resolves false for Cancel.
   * The check runs in the browser and constructs no diagram Worker.
   * @param {readonly unknown[]} fontFamilies Each drawing's `adv.font_family` in the Session.
   * @returns {Promise<boolean>}
   */
  const confirmSessionFonts = async (fontFamilies) => {
    /** @type {FontFamilyName | null} */
    let loadedFamily = null;
    /** @type {string} */
    let problem = '';
    for (;;) {
      const missing = missingSessionFonts(fontFamilies, loadedFamilies());
      if (!missing) return true;
      const question = sessionFontQuestion(missing, { installedFontsSupported: localFontsSupported, loadedFamily });
      const choice = await askChoice(problem ? { ...question, items: [...question.items, problem] } : question);
      loadedFamily = null;
      problem = '';
      if (choice === SESSION_FONT_ACTIONS.SUBSTITUTE) return true;
      if (choice !== SESSION_FONT_ACTIONS.FILE && choice !== SESSION_FONT_ACTIONS.INSTALLED) return false;
      try {
        /** @type {WebFontSource[]} */
        let sources = [];
        if (choice === SESSION_FONT_ACTIONS.FILE) {
          const files = await chooseFontFiles();
          sources = await readSources(() => fileSources(files));
        } else {
          const faces = await pickInstalledFamily();
          if (faces) sources = await readSources(() => installedSources(faces));
        }
        if (sources.length > 0) loadedFamily = (await loadSources(sources))[0].family;
      } catch (error) {
        problem = loadErrorSummary(error);
      }
    }
  };

  return {
    choices,
    selectedChoice,
    customFontSelected,
    embedInExports,
    markFontsLost,
    localFontsSupported,
    loading,
    loadError,
    loadFontFile,
    confirmSessionFonts,
    installedPicker,
    filteredInstalledFamilies,
    openInstalledFonts,
    chooseInstalledFamily,
    closeInstalledFonts
  };
};
