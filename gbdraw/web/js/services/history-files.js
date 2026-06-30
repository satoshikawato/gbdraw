const isFileLike = (value) =>
  Boolean(value) &&
  typeof value === 'object' &&
  typeof value.name === 'string' &&
  typeof value.size === 'number';

export const collectHistoryFileIds = (value, target = new Set()) => {
  if (Array.isArray(value)) {
    value.forEach((entry) => collectHistoryFileIds(entry, target));
    return target;
  }
  if (!value || typeof value !== 'object') return target;
  if (typeof value.fileId === 'string' && value.fileId) {
    target.add(value.fileId);
  }
  Object.values(value).forEach((entry) => collectHistoryFileIds(entry, target));
  return target;
};

export const createHistoryFileStore = () => {
  const files = new Map();
  const idsByFile = new WeakMap();
  let nextId = 1;

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

  const describeFile = (file) => {
    const fileId = registerFile(file);
    if (!fileId) return null;
    return {
      fileId,
      name: file.name || 'file',
      type: file.type || '',
      size: file.size || 0,
      lastModified: file.lastModified || 0
    };
  };

  const describeValue = (value) => {
    if (Array.isArray(value)) return value.map((entry) => describeValue(entry));
    return describeFile(value);
  };

  const restoreFile = (descriptor) => {
    const fileId = String(descriptor?.fileId || '');
    if (!fileId) return null;
    return files.get(fileId)?.file || null;
  };

  const restoreValue = (value) => {
    if (Array.isArray(value)) return value.map((entry) => restoreValue(entry));
    return restoreFile(value);
  };

  const retainOnly = (fileIds) => {
    const retained = fileIds instanceof Set ? fileIds : new Set(fileIds || []);
    Array.from(files.keys()).forEach((fileId) => {
      if (!retained.has(fileId)) files.delete(fileId);
    });
  };

  const estimateBytes = (fileIds) => {
    const uniqueIds = fileIds instanceof Set ? fileIds : new Set(fileIds || []);
    let total = 0;
    uniqueIds.forEach((fileId) => {
      total += files.get(fileId)?.file?.size || 0;
    });
    return total;
  };

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
