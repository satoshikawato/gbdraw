// @ts-check
/** @import { FileLike } from './session-resource-backing.js' */
/** A File, or any object the History keeps as one. */
/** @typedef {{ name: string, size: number, type?: string, lastModified?: number }} HistoryFile */
/** What `registerFile` accepts: a File, a Session resource view, or anything else (which is not registered). @typedef {{ name?: string, size?: number, type?: string, lastModified?: number } | null | undefined} InputFileCandidate */
/**
 * A file binding: the frozen name of one stored file (`InputFileStore.registerFile`).
 * It is not a `Blob`; its bytes are `InputFileStore.fileOf(ref)`. `sourceId`
 * names where the bytes come from: the refs of two Session resource views of
 * one stored resource share it (a Linear File card groups its records by it);
 * any other file's is its `fileId`.
 * @typedef {Readonly<{ fileId: string, sourceId: string, name: string, type: string, size: number, lastModified: number }>} InputFileRef
 */

/**
 * @param {unknown} value
 * @returns {value is HistoryFile}
 */
const isFileLike = (value) =>
  Boolean(value) &&
  typeof value === 'object' &&
  typeof /** @type {{ name?: unknown }} */ (value).name === 'string' &&
  typeof /** @type {{ size?: unknown }} */ (value).size === 'number' &&
  // A ref names a stored file; it is never stored as one.
  !Object.hasOwn(/** @type {object} */ (value), 'fileId');

/**
 * @param {unknown} value
 * @returns {value is InputFileRef}
 */
export const isInputFileRef = (value) => Boolean(value) && typeof value === 'object'
  && typeof /** @type {{ fileId?: unknown }} */ (value).fileId === 'string';

/**
 * @param {unknown} value
 * @param {Set<string>} [target]
 * @returns {Set<string>}
 */
export const collectHistoryFileIds = (value, target = new Set()) => {
  if (Array.isArray(value)) {
    value.forEach((entry) => collectHistoryFileIds(entry, target));
    return target;
  }
  if (!value || typeof value !== 'object') return target;
  const { fileId } = /** @type {{ fileId?: unknown }} */ (value);
  if (typeof fileId === 'string' && fileId) {
    target.add(fileId);
  }
  Object.values(value).forEach((entry) => collectHistoryFileIds(entry, target));
  return target;
};

/**
 * The project's input file store (`state.inputFiles`): one per project,
 * shared by History and the input bindings. It holds each File (or Session
 * resource view) once, deduped by object identity, and names it by one frozen
 * `InputFileRef`; an input field or a History step holds the ref, and a byte
 * reader resolves it with `fileOf`. A file is released only when no History
 * step and no binding names it (`retainOnly` takes the union).
 */
/**
 * @typedef {object} InputFileStore
 * @property {(file: InputFileCandidate) => InputFileRef | null} registerFile
 * @property {(ref: InputFileRef | null | undefined) => FileLike | null} fileOf
 * @property {(fileIds: Iterable<string> | null | undefined) => void} retainOnly
 * @property {(fileIds: Iterable<string> | null | undefined) => number} estimateBytes
 * @property {(fileId: string) => boolean} has
 * @property {typeof collectHistoryFileIds} collectIds
 */
/**
 * @param {{ sourceOf?: (file: HistoryFile) => object | null | undefined }} [options] `sourceOf`: the
 *   object a file's bytes come from when other files share them (a Session resource view's descriptor);
 *   refs of files with one source share `sourceId`. The project's store is given the Session's
 *   (`state.js`), so this module imports nothing.
 * @returns {InputFileStore}
 */
export const createInputFileStore = ({ sourceOf = () => null } = {}) => {
  /** @type {Map<string, { file: HistoryFile, ref: InputFileRef }>} */
  const files = new Map();
  /** @type {WeakMap<HistoryFile, string>} */
  const idsByFile = new WeakMap();
  /** @type {WeakMap<object, string>} The source of each stored resource's views, by its descriptor. */
  const sourceIdsByDescriptor = new WeakMap();
  let nextId = 1;

  /**
   * The one frozen ref of a file, registering the file on its first call.
   * @param {InputFileCandidate} file
   * @returns {InputFileRef | null}
   */
  const registerFile = (file) => {
    // A ref of this store stays bound; a ref it has released names nothing.
    if (isInputFileRef(file)) return files.get(file.fileId)?.ref || null;
    if (!isFileLike(file)) return null;
    const existingId = idsByFile.get(file);
    const existing = existingId ? files.get(existingId) : undefined;
    if (existing) return existing.ref;
    const fileId = `file-${nextId}`;
    nextId += 1;
    idsByFile.set(file, fileId);
    const descriptor = sourceOf(file);
    const sourceId = descriptor ? (sourceIdsByDescriptor.get(descriptor) || fileId) : fileId;
    if (descriptor && !sourceIdsByDescriptor.has(descriptor)) sourceIdsByDescriptor.set(descriptor, sourceId);
    /** @type {InputFileRef} */
    const ref = Object.freeze({
      fileId,
      sourceId,
      name: file.name || 'file',
      type: file.type || '',
      size: file.size || 0,
      lastModified: file.lastModified || 0
    });
    files.set(fileId, { file, ref });
    return ref;
  };

  /**
   * The bytes a ref names: the registered File, or null once released.
   * @param {InputFileRef | null | undefined} ref
   * @returns {FileLike | null}
   */
  const fileOf = (ref) => (ref ? files.get(ref.fileId)?.file || null : null);

  /** @param {Iterable<string> | null | undefined} fileIds */
  const retainOnly = (fileIds) => {
    const retained = fileIds instanceof Set ? fileIds : new Set(fileIds || []);
    Array.from(files.keys()).forEach((fileId) => {
      if (!retained.has(fileId)) files.delete(fileId);
    });
  };

  /** @param {Iterable<string> | null | undefined} fileIds */
  const estimateBytes = (fileIds) => {
    const uniqueIds = fileIds instanceof Set ? fileIds : new Set(fileIds || []);
    let total = 0;
    uniqueIds.forEach((fileId) => {
      total += files.get(fileId)?.file?.size || 0;
    });
    return total;
  };

  /** @param {string} fileId */
  const has = (fileId) => files.has(fileId);

  return {
    collectIds: collectHistoryFileIds,
    estimateBytes,
    fileOf,
    has,
    registerFile,
    retainOnly
  };
};

