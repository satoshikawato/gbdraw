// @ts-check
/** A File, or any object the History keeps as one. */
/** @typedef {{ name: string, size: number, type?: string, lastModified?: number }} HistoryFile */
/**
 * A file binding: the frozen name of one stored file (`InputFileStore.registerFile`).
 * It is not a `Blob`; its bytes are `InputFileStore.fileOf(ref)`.
 * @typedef {Readonly<{ fileId: string, name: string, type: string, size: number, lastModified: number }>} InputFileRef
 */
/** @typedef {InputFileRef} HistoryFileDescriptor */
/** @typedef {HistoryFile | null | undefined | HistoryFileValue[]} HistoryFileValue */
/** @typedef {HistoryFileDescriptor | null | HistoryFileDescriptorValue[]} HistoryFileDescriptorValue */

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
 * @property {(file: HistoryFile | null | undefined) => InputFileRef | null} registerFile
 * @property {(ref: InputFileRef | null | undefined) => HistoryFile | null} fileOf
 * @property {(value: HistoryFileValue) => HistoryFileDescriptorValue} describeValue
 * @property {(value: HistoryDescriptorTree) => HistoryFileValue} restoreValue
 * @property {(fileIds: Iterable<string> | null | undefined) => void} retainOnly
 * @property {(fileIds: Iterable<string> | null | undefined) => number} estimateBytes
 * @property {(fileId: string) => boolean} has
 * @property {typeof collectHistoryFileIds} collectIds
 */
/** @typedef {InputFileRef | null | undefined | HistoryDescriptorTree[]} HistoryDescriptorTree */
/** @returns {InputFileStore} */
export const createHistoryFileStore = () => {
  /** @type {Map<string, { file: HistoryFile, ref: InputFileRef }>} */
  const files = new Map();
  /** @type {WeakMap<HistoryFile, string>} */
  const idsByFile = new WeakMap();
  let nextId = 1;

  /**
   * The one frozen ref of a file, registering the file on its first call.
   * @param {HistoryFile | null | undefined} file
   * @returns {InputFileRef | null}
   */
  const registerFile = (file) => {
    if (!isFileLike(file)) return null;
    const existingId = idsByFile.get(file);
    const existing = existingId ? files.get(existingId) : undefined;
    if (existing) return existing.ref;
    const fileId = `file-${nextId}`;
    nextId += 1;
    idsByFile.set(file, fileId);
    /** @type {InputFileRef} */
    const ref = Object.freeze({
      fileId,
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
   * @returns {HistoryFile | null}
   */
  const fileOf = (ref) => (ref ? files.get(ref.fileId)?.file || null : null);

  /**
   * @param {HistoryFileValue} value
   * @returns {HistoryFileDescriptorValue}
   */
  const describeValue = (value) => {
    if (Array.isArray(value)) return value.map((entry) => describeValue(entry));
    return registerFile(value);
  };

  /**
   * @param {HistoryDescriptorTree} value
   * @returns {HistoryFileValue}
   */
  const restoreValue = (value) => {
    if (Array.isArray(value)) return value.map((entry) => restoreValue(entry));
    return fileOf(value);
  };

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
    describeValue,
    estimateBytes,
    fileOf,
    has,
    registerFile,
    restoreValue,
    retainOnly
  };
};

