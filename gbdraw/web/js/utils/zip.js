// @ts-check
// The ZIP format (PKWARE APPNOTE) as gbdraw writes it: CRC-32, and the local
// header, central directory record and end records of an entry. Two writers use
// them: `createZipBlob` (several output files in one download) and the Session
// file container (`services/session-container.js`). No import: tests load this
// module alone.

// CRC-32 (ISO-HDLC) by slicing-by-8: eight tables, eight bytes per step.
const CRC32_TABLES = (() => {
  const tables = new Uint32Array(256 * 8);
  for (let i = 0; i < 256; i += 1) {
    let value = i;
    for (let bit = 0; bit < 8; bit += 1) {
      value = (value & 1) ? (0xEDB88320 ^ (value >>> 1)) : (value >>> 1);
    }
    tables[i] = value >>> 0;
  }
  for (let i = 0; i < 256; i += 1) {
    for (let slice = 1; slice < 8; slice += 1) {
      const previous = tables[(slice - 1) * 256 + i];
      tables[slice * 256 + i] = (previous >>> 8) ^ tables[previous & 0xFF];
    }
  }
  return tables;
})();

/**
 * The CRC-32 of the bytes before and `bytes`, from the CRC-32 `crc` of the
 * bytes before (0 for none), so a stream is checked one chunk at a time.
 * @param {number} crc
 * @param {Uint8Array} bytes
 * @returns {number}
 */
export const crc32Update = (crc, bytes) => {
  const t = CRC32_TABLES;
  let value = ~crc;
  let index = 0;
  for (const last = bytes.length - 8; index <= last; index += 8) {
    const low = (value ^ (bytes[index] | (bytes[index + 1] << 8) | (bytes[index + 2] << 16) | (bytes[index + 3] << 24))) >>> 0;
    value = t[1792 + (low & 0xFF)] ^ t[1536 + ((low >>> 8) & 0xFF)] ^ t[1280 + ((low >>> 16) & 0xFF)] ^ t[1024 + (low >>> 24)]
      ^ t[768 + bytes[index + 4]] ^ t[512 + bytes[index + 5]] ^ t[256 + bytes[index + 6]] ^ t[bytes[index + 7]];
  }
  for (; index < bytes.length; index += 1) value = t[(value ^ bytes[index]) & 0xFF] ^ (value >>> 8);
  return ~value >>> 0;
};

const textEncoder = new TextEncoder();

/**
 * @param {unknown} data
 * @returns {Uint8Array}
 */
const toBytes = (data) => {
  if (data instanceof Uint8Array) return data;
  if (data instanceof ArrayBuffer) return new Uint8Array(data);
  if (ArrayBuffer.isView(data)) return new Uint8Array(data.buffer, data.byteOffset, data.byteLength);
  return textEncoder.encode(String(data ?? ''));
};

const dosTimestamp = (date = new Date()) => {
  const year = Math.max(1980, date.getFullYear());
  const time = (
    (date.getHours() << 11) |
    (date.getMinutes() << 5) |
    Math.floor(date.getSeconds() / 2)
  ) & 0xFFFF;
  const day = Math.max(1, date.getDate());
  const month = Math.max(1, date.getMonth() + 1);
  const packedDate = (
    ((year - 1980) << 9) |
    (month << 5) |
    day
  ) & 0xFFFF;
  return { time, date: packedDate };
};

/**
 * @param {unknown} name
 * @param {unknown} fallback
 */
const sanitizeZipPath = (name, fallback) => {
  const raw = String(name || fallback || 'file').replace(/\\/g, '/');
  const parts = raw
    .split('/')
    .map((part) => part.trim())
    .filter((part) => part && part !== '.' && part !== '..');
  return parts.join('/') || String(fallback || 'file');
};

/**
 * @param {unknown} name
 * @param {Set<string>} usedNames The paths already taken; the result is added.
 */
const uniqueZipPath = (name, usedNames) => {
  const normalized = sanitizeZipPath(name, `file-${usedNames.size + 1}`);
  if (!usedNames.has(normalized)) {
    usedNames.add(normalized);
    return normalized;
  }

  const slashIndex = normalized.lastIndexOf('/');
  const directory = slashIndex >= 0 ? `${normalized.slice(0, slashIndex + 1)}` : '';
  const basename = slashIndex >= 0 ? normalized.slice(slashIndex + 1) : normalized;
  const dotIndex = basename.lastIndexOf('.');
  const stem = dotIndex > 0 ? basename.slice(0, dotIndex) : basename;
  const extension = dotIndex > 0 ? basename.slice(dotIndex) : '';
  let suffix = 2;
  let candidate = `${directory}${stem}-${suffix}${extension}`;
  while (usedNames.has(candidate)) {
    suffix += 1;
    candidate = `${directory}${stem}-${suffix}${extension}`;
  }
  usedNames.add(candidate);
  return candidate;
};

// A size or offset above this is written in the ZIP64 extra field, and the
// 32-bit field holds 0xFFFFFFFF; Python's `zipfile` uses the same limit.
export const ZIP64_LIMIT = 0x7FFFFFFF;
const ZIP32_MARK = 0xFFFFFFFF;
const UTF8_NAME_FLAG = 0x0800;

/**
 * One entry as its local header and its central directory record describe
 * it. `zip64` writes both sizes in the ZIP64 extra field even when they are
 * small (a test of the ZIP64 reader); larger sizes are written there anyway.
 * @typedef {{
 *   name: Uint8Array, method: number, crc: number, compressedSize: number, size: number,
 *   offset: number, time: number, date: number, zip64?: boolean
 * }} ZipEntryFields
 */

/** @param {ZipEntryFields} entry */
const sizesInZip64 = (entry) => Boolean(entry.zip64) || entry.size > ZIP64_LIMIT || entry.compressedSize > ZIP64_LIMIT;

/**
 * A ZIP64 extended information extra field holding `values` in order.
 * @param {number[]} values
 */
const zip64Extra = (values) => {
  const extra = new Uint8Array(values.length ? 4 + 8 * values.length : 0);
  if (!values.length) return extra;
  const view = new DataView(extra.buffer);
  view.setUint16(0, 1, true);
  view.setUint16(2, 8 * values.length, true);
  values.forEach((value, index) => view.setBigUint64(4 + 8 * index, BigInt(value), true));
  return extra;
};

/**
 * @param {ZipEntryFields} entry
 * @returns {Uint8Array<ArrayBuffer>} The entry's local header; its data follows it.
 */
export const zipLocalHeader = (entry) => {
  const zip64 = sizesInZip64(entry);
  const extra = zip64Extra(zip64 ? [entry.size, entry.compressedSize] : []);
  const header = new Uint8Array(30 + entry.name.length + extra.length);
  const view = new DataView(header.buffer);
  view.setUint32(0, 0x04034B50, true);
  view.setUint16(4, zip64 ? 45 : 20, true);
  view.setUint16(6, UTF8_NAME_FLAG, true);
  view.setUint16(8, entry.method, true);
  view.setUint16(10, entry.time, true);
  view.setUint16(12, entry.date, true);
  view.setUint32(14, entry.crc >>> 0, true);
  view.setUint32(18, zip64 ? ZIP32_MARK : entry.compressedSize, true);
  view.setUint32(22, zip64 ? ZIP32_MARK : entry.size, true);
  view.setUint16(26, entry.name.length, true);
  view.setUint16(28, extra.length, true);
  header.set(entry.name, 30);
  header.set(extra, 30 + entry.name.length);
  return header;
};

/**
 * @param {ZipEntryFields} entry
 * @returns {Uint8Array<ArrayBuffer>} The entry's central directory record.
 */
export const zipCentralRecord = (entry) => {
  const zip64Sizes = sizesInZip64(entry);
  const zip64Offset = entry.offset > ZIP64_LIMIT;
  const extra = zip64Extra([...(zip64Sizes ? [entry.size, entry.compressedSize] : []), ...(zip64Offset ? [entry.offset] : [])]);
  const version = extra.length ? 45 : 20;
  const record = new Uint8Array(46 + entry.name.length + extra.length);
  const view = new DataView(record.buffer);
  view.setUint32(0, 0x02014B50, true);
  view.setUint16(4, version, true);
  view.setUint16(6, version, true);
  view.setUint16(8, UTF8_NAME_FLAG, true);
  view.setUint16(10, entry.method, true);
  view.setUint16(12, entry.time, true);
  view.setUint16(14, entry.date, true);
  view.setUint32(16, entry.crc >>> 0, true);
  view.setUint32(20, zip64Sizes ? ZIP32_MARK : entry.compressedSize, true);
  view.setUint32(24, zip64Sizes ? ZIP32_MARK : entry.size, true);
  view.setUint16(28, entry.name.length, true);
  view.setUint16(30, extra.length, true);
  view.setUint32(42, zip64Offset ? ZIP32_MARK : entry.offset, true);
  record.set(entry.name, 46);
  record.set(extra, 46 + entry.name.length);
  return record;
};

/**
 * The end of central directory record, after the ZIP64 end record and its
 * locator when the count, size or offset needs them (or `zip64` asks for them).
 * @param {{ count: number, directoryOffset: number, directorySize: number, zip64?: boolean }} directory
 * @returns {Uint8Array<ArrayBuffer>}
 */
export const zipEndRecords = ({ count, directoryOffset, directorySize, zip64 = false }) => {
  const large = zip64 || count > 0xFFFF || directoryOffset > ZIP64_LIMIT || directorySize > ZIP64_LIMIT;
  const records = new Uint8Array((large ? 76 : 0) + 22);
  const view = new DataView(records.buffer);
  if (large) {
    view.setUint32(0, 0x06064B50, true);
    view.setBigUint64(4, 44n, true);
    view.setUint16(12, 45, true);
    view.setUint16(14, 45, true);
    view.setBigUint64(24, BigInt(count), true);
    view.setBigUint64(32, BigInt(count), true);
    view.setBigUint64(40, BigInt(directorySize), true);
    view.setBigUint64(48, BigInt(directoryOffset), true);
    view.setUint32(56, 0x07064B50, true);
    view.setBigUint64(64, BigInt(directoryOffset + directorySize), true);
    view.setUint32(72, 1, true);
  }
  const end = large ? 76 : 0;
  view.setUint32(end, 0x06054B50, true);
  view.setUint16(end + 8, Math.min(count, 0xFFFF), true);
  view.setUint16(end + 10, Math.min(count, 0xFFFF), true);
  view.setUint32(end + 12, large ? ZIP32_MARK : directorySize, true);
  view.setUint32(end + 16, large ? ZIP32_MARK : directoryOffset, true);
  return records;
};

/**
 * @typedef {{ name?: string, data?: string | ArrayBuffer | ArrayBufferView | null }} ZipFileEntry
 */

/**
 * @param {ZipFileEntry[] | null | undefined} files
 * @returns {Blob}
 */
export const createZipBlob = (files) => {
  /** @type {Uint8Array[]} */
  const localParts = [];
  /** @type {Uint8Array[]} */
  const centralParts = [];
  /** @type {Set<string>} */
  const usedNames = new Set();
  const { time, date } = dosTimestamp();
  let offset = 0;

  (Array.isArray(files) ? files : []).forEach((file) => {
    const name = textEncoder.encode(uniqueZipPath(file?.name, usedNames));
    const dataBytes = toBytes(file?.data);
    /** @type {ZipEntryFields} */
    const entry = {
      name, method: 0, crc: crc32Update(0, dataBytes), compressedSize: dataBytes.length, size: dataBytes.length,
      offset, time, date
    };
    const localHeader = zipLocalHeader(entry);
    localParts.push(localHeader, dataBytes);
    centralParts.push(zipCentralRecord(entry));
    offset += localHeader.length + dataBytes.length;
  });

  const directorySize = centralParts.reduce((sum, part) => sum + part.length, 0);
  const end = zipEndRecords({ count: centralParts.length, directoryOffset: offset, directorySize });
  // Entry bytes never sit in shared memory, which a Blob part rejects.
  return new Blob(/** @type {BlobPart[]} */ ([...localParts, ...centralParts, end]), { type: 'application/zip' });
};
