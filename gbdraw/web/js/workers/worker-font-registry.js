// @ts-check
import { sha256Hex } from '../services/byte-utils.js';

// The fonts registered in a diagram Worker's Python (custom fonts D-01): their
// files stay in FONT_DIRECTORY for the Worker's life, and each identity is the
// SHA-256 of the bytes. The main thread replays every font of the page into a
// replacement Worker (init `fonts`), and a Generate names the fonts it needs.

const FONT_DIRECTORY = '/gbdraw-web-fonts';
const FONT_IDENTITY_RE = /^[0-9a-f]{64}$/;

/**
 * A font file the main thread sends: the identity of its bytes, its file name
 * or installed face name, and the bytes.
 * @typedef {object} WorkerFontFile
 * @property {string} identity
 * @property {string} label
 * @property {ArrayBuffer | ArrayBufferView} bytes
 */

/** One file Python registers. @typedef {{ identity: string, path: string }} PythonFontEntry */

/**
 * The part of Pyodide's file system the registry uses.
 * @typedef {object} WorkerFontFileSystem
 * @property {(path: string) => { exists: boolean }} analyzePath
 * @property {(path: string) => void} mkdir
 * @property {(path: string, bytes: Uint8Array) => void} writeFile
 * @property {(path: string) => void} unlink
 */

/**
 * A font the replay into a new Worker could not register.
 * @typedef {{ identity: string, label: string }} LostWorkerFont
 */

/** @param {unknown} value @returns {Uint8Array | null} */
const fontBytes = (value) => {
  if (value instanceof ArrayBuffer) return new Uint8Array(value);
  if (ArrayBuffer.isView(value)) return new Uint8Array(value.buffer, value.byteOffset, value.byteLength);
  return null;
};

/**
 * @param {object} options
 * @param {WorkerFontFileSystem} options.fs
 * @param {(entries: PythonFontEntry[]) => { fonts?: unknown }} options.registerInPython
 *   `register_web_fonts_json`: registers all of the files or none, throws on failure.
 */
export const createWorkerFontRegistry = ({ fs, registerInPython }) => {
  /** @type {Set<string>} */
  const registered = new Set();

  /**
   * Register font files; a call registers all of them or none.
   * @param {unknown} fonts
   * @returns {Promise<unknown>} Python's facts per font, in order.
   */
  const register = async (fonts) => {
    if (!Array.isArray(fonts)) throw { code: 'HELPER_PROTOCOL', stage: 'request-validation' };
    if (!fs.analyzePath(FONT_DIRECTORY).exists) fs.mkdir(FONT_DIRECTORY);
    /** @type {PythonFontEntry[]} */
    const entries = [];
    /** @type {string[]} */
    const written = [];
    try {
      for (const font of fonts) {
        const file = /** @type {Partial<WorkerFontFile> | null} */ (font);
        const identity = String(file?.identity || '');
        const bytes = fontBytes(file?.bytes);
        if (!FONT_IDENTITY_RE.test(identity) || !bytes || await sha256Hex(bytes) !== identity) {
          throw { code: 'RESOURCE_INVALID', stage: 'resource-staging', context: { field: 'font_family' } };
        }
        const path = `${FONT_DIRECTORY}/${identity}`;
        if (!registered.has(identity)) {
          fs.writeFile(path, bytes);
          written.push(path);
        }
        entries.push({ identity, path });
      }
      const result = registerInPython(entries);
      entries.forEach(({ identity }) => registered.add(identity));
      return result.fonts;
    } catch (error) {
      written.forEach((path) => {
        try {
          fs.unlink(path);
        } catch (_error) {
          // The file was not written.
        }
      });
      throw error;
    }
  };

  /**
   * Register the page's fonts in a new Worker one by one, so one bad font
   * does not drop the others; returns the fonts it could not register.
   * @param {unknown} fonts
   * @returns {Promise<LostWorkerFont[]>}
   */
  const replay = async (fonts) => {
    /** @type {LostWorkerFont[]} */
    const lost = [];
    for (const font of Array.isArray(fonts) ? fonts : []) {
      try {
        await register([font]);
      } catch (_error) {
        const file = /** @type {Partial<WorkerFontFile> | null} */ (font);
        lost.push({ identity: String(file?.identity || ''), label: String(file?.label || '') });
      }
    }
    return lost;
  };

  /**
   * A Generate draws only with the fonts it names; a font the Worker lost
   * fails the Generate instead of falling back.
   * @param {unknown} fontIdentities
   */
  const requireRegistered = (fontIdentities) => {
    const missing = (Array.isArray(fontIdentities) ? fontIdentities : [])
      .filter((identity) => !registered.has(String(identity)));
    if (missing.length > 0) {
      throw { code: 'INPUT_INVALID', stage: 'request-validation', context: { field: 'font_family', reason: 'FONT_NOT_LOADED' } };
    }
  };

  return { register, replay, requireRegistered };
};
