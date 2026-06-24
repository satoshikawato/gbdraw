const CRC32_TABLE = (() => {
  const table = new Uint32Array(256);
  for (let i = 0; i < 256; i += 1) {
    let value = i;
    for (let bit = 0; bit < 8; bit += 1) {
      value = (value & 1) ? (0xEDB88320 ^ (value >>> 1)) : (value >>> 1);
    }
    table[i] = value >>> 0;
  }
  return table;
})();

const textEncoder = new TextEncoder();

const toBytes = (data) => {
  if (data instanceof Uint8Array) return data;
  if (data instanceof ArrayBuffer) return new Uint8Array(data);
  if (ArrayBuffer.isView(data)) return new Uint8Array(data.buffer, data.byteOffset, data.byteLength);
  return textEncoder.encode(String(data ?? ''));
};

const crc32 = (bytes) => {
  let crc = 0xFFFFFFFF;
  for (let i = 0; i < bytes.length; i += 1) {
    crc = CRC32_TABLE[(crc ^ bytes[i]) & 0xFF] ^ (crc >>> 8);
  }
  return (crc ^ 0xFFFFFFFF) >>> 0;
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

const sanitizeZipPath = (name, fallback) => {
  const raw = String(name || fallback || 'file').replace(/\\/g, '/');
  const parts = raw
    .split('/')
    .map((part) => part.trim())
    .filter((part) => part && part !== '.' && part !== '..');
  return parts.join('/') || String(fallback || 'file');
};

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

const writeUint16 = (view, offset, value) => view.setUint16(offset, value, true);
const writeUint32 = (view, offset, value) => view.setUint32(offset, value >>> 0, true);

export const createZipBlob = (files) => {
  const localParts = [];
  const centralParts = [];
  const usedNames = new Set();
  const { time, date } = dosTimestamp();
  let offset = 0;

  (Array.isArray(files) ? files : []).forEach((file) => {
    const zipPath = uniqueZipPath(file?.name, usedNames);
    const nameBytes = textEncoder.encode(zipPath);
    const dataBytes = toBytes(file?.data);
    const crc = crc32(dataBytes);
    const localHeader = new Uint8Array(30 + nameBytes.length);
    const localView = new DataView(localHeader.buffer);
    writeUint32(localView, 0, 0x04034B50);
    writeUint16(localView, 4, 20);
    writeUint16(localView, 6, 0x0800);
    writeUint16(localView, 8, 0);
    writeUint16(localView, 10, time);
    writeUint16(localView, 12, date);
    writeUint32(localView, 14, crc);
    writeUint32(localView, 18, dataBytes.length);
    writeUint32(localView, 22, dataBytes.length);
    writeUint16(localView, 26, nameBytes.length);
    writeUint16(localView, 28, 0);
    localHeader.set(nameBytes, 30);
    localParts.push(localHeader, dataBytes);

    const centralHeader = new Uint8Array(46 + nameBytes.length);
    const centralView = new DataView(centralHeader.buffer);
    writeUint32(centralView, 0, 0x02014B50);
    writeUint16(centralView, 4, 20);
    writeUint16(centralView, 6, 20);
    writeUint16(centralView, 8, 0x0800);
    writeUint16(centralView, 10, 0);
    writeUint16(centralView, 12, time);
    writeUint16(centralView, 14, date);
    writeUint32(centralView, 16, crc);
    writeUint32(centralView, 20, dataBytes.length);
    writeUint32(centralView, 24, dataBytes.length);
    writeUint16(centralView, 28, nameBytes.length);
    writeUint16(centralView, 30, 0);
    writeUint16(centralView, 32, 0);
    writeUint16(centralView, 34, 0);
    writeUint16(centralView, 36, 0);
    writeUint32(centralView, 38, 0);
    writeUint32(centralView, 42, offset);
    centralHeader.set(nameBytes, 46);
    centralParts.push(centralHeader);

    offset += localHeader.length + dataBytes.length;
  });

  const centralOffset = offset;
  const centralSize = centralParts.reduce((sum, part) => sum + part.length, 0);
  const entryCount = centralParts.length;
  const eocd = new Uint8Array(22);
  const eocdView = new DataView(eocd.buffer);
  writeUint32(eocdView, 0, 0x06054B50);
  writeUint16(eocdView, 4, 0);
  writeUint16(eocdView, 6, 0);
  writeUint16(eocdView, 8, entryCount);
  writeUint16(eocdView, 10, entryCount);
  writeUint32(eocdView, 12, centralSize);
  writeUint32(eocdView, 16, centralOffset);
  writeUint16(eocdView, 20, 0);

  return new Blob([...localParts, ...centralParts, eocd], { type: 'application/zip' });
};

export const downloadZipFile = (filename, files) => {
  const blob = createZipBlob(files);
  const url = URL.createObjectURL(blob);
  const link = document.createElement('a');
  link.href = url;
  link.download = filename || 'gbdraw-cli-files.zip';
  link.addEventListener('click', (event) => {
    event.stopPropagation();
  }, { once: true });
  link.click();
  URL.revokeObjectURL(url);
};
