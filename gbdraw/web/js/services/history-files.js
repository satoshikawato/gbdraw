// @ts-check
/** A File, or any object the History keeps as one. */
/** @typedef {{ name: string, size: number, type?: string, lastModified?: number }} HistoryFile */
/** @typedef {{ fileId: string, name: string, type: string, size: number, lastModified: number }} HistoryFileDescriptor */
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
  typeof /** @type {{ size?: unknown }} */ (value).size === 'number';

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

export const createHistoryFileStore = () => {
  /** @type {Map<string, { file: HistoryFile }>} */
  const files = new Map();
  /** @type {WeakMap<HistoryFile, string>} */
  const idsByFile = new WeakMap();
  let nextId = 1;

  /**
   * @param {HistoryFile | null | undefined} file
   * @returns {string | null}
   */
  const registerFile = (file) => {
    if (!isFileLike(file)) return null;
    const existingId = idsByFile.get(file);
    if (existingId && files.has(existingId)) return existingId;
    const fileId = `file-${nextId}`;
    nextId += 1;
    idsByFile.set(file, fileId);
    files.set(fileId, { file });
    return fileId;
  };

  /**
   * @param {HistoryFile | null | undefined} file
   * @returns {HistoryFileDescriptor | null}
   */
  const describeFile = (file) => {
    const fileId = registerFile(file);
    if (!fileId || !file) return null;
    return {
      fileId,
      name: file.name || 'file',
      type: file.type || '',
      size: file.size || 0,
      lastModified: file.lastModified || 0
    };
  };

  /**
   * @param {HistoryFileValue} value
   * @returns {HistoryFileDescriptorValue}
   */
  const describeValue = (value) => {
    if (Array.isArray(value)) return value.map((entry) => describeValue(entry));
    return describeFile(value);
  };

  /**
   * @param {{ fileId?: unknown } | null | undefined} descriptor
   * @returns {HistoryFile | null}
   */
  const restoreFile = (descriptor) => {
    const fileId = String(descriptor?.fileId || '');
    if (!fileId) return null;
    return files.get(fileId)?.file || null;
  };

  /**
   * @typedef {{ fileId?: unknown } | null | undefined | HistoryDescriptorTree[]} HistoryDescriptorTree
   */
  /**
   * @param {HistoryDescriptorTree} value
   * @returns {HistoryFileValue}
   */
  const restoreValue = (value) => {
    if (Array.isArray(value)) return value.map((entry) => restoreValue(entry));
    return restoreFile(value);
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
    describeFile,
    describeValue,
    estimateBytes,
    has,
    registerFile,
    restoreFile,
    restoreValue,
    retainOnly
  };
};
