// @ts-check
const BASE64_ALPHABET = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/';

/** @param {string} character */
const isBase64Whitespace = (character) => (
  character === ' '
  || character === '\t'
  || character === '\n'
  || character === '\f'
  || character === '\r'
);

/** @param {Uint8Array<ArrayBuffer> | ArrayLike<number> | ArrayBuffer | null | undefined} value */
export const asBytes = (value) => (
  value instanceof Uint8Array ? value : new Uint8Array(value || [])
);

const BASE64_CODES = Uint8Array.from(BASE64_ALPHABET, (character) => character.charCodeAt(0));
const PADDING_CODE = 0x3D;

// Older engines write the ASCII text into one buffer; a binary-string btoa
// path runs at about 16 MB/s and freezes the first Generate of a large File.
/** @param {Uint8Array} bytes */
const bytesToBase64Portable = (bytes) => {
  const length = bytes.length;
  const output = new Uint8Array(Math.ceil(length / 3) * 4);
  const whole = length - (length % 3);
  let out = 0;
  for (let index = 0; index < whole; index += 3) {
    const triplet = (bytes[index] << 16) | (bytes[index + 1] << 8) | bytes[index + 2];
    output[out] = BASE64_CODES[triplet >> 18];
    output[out + 1] = BASE64_CODES[(triplet >> 12) & 63];
    output[out + 2] = BASE64_CODES[(triplet >> 6) & 63];
    output[out + 3] = BASE64_CODES[triplet & 63];
    out += 4;
  }
  if (whole < length) {
    const second = whole + 1 < length ? bytes[whole + 1] : 0;
    const triplet = (bytes[whole] << 16) | (second << 8);
    output[out] = BASE64_CODES[triplet >> 18];
    output[out + 1] = BASE64_CODES[(triplet >> 12) & 63];
    output[out + 2] = whole + 1 < length ? BASE64_CODES[(triplet >> 6) & 63] : PADDING_CODE;
    output[out + 3] = PADDING_CODE;
  }
  return new TextDecoder().decode(output);
};

/**
 * The one base64 encoder for resources and text tables. Engines with
 * `Uint8Array.prototype.toBase64` encode natively; both paths give the same
 * padded standard-alphabet text.
 * @param {Uint8Array<ArrayBuffer> | ArrayLike<number> | ArrayBuffer | null | undefined} value
 */
export const bytesToBase64 = (value) => {
  const bytes = asBytes(value);
  return typeof bytes.toBase64 === 'function'
    ? bytes.toBase64()
    : bytesToBase64Portable(bytes);
};


/**
 * Decodes with `atob` in bounded chunks straight into the final buffer; it
 * accepts and refuses what one `atob` of the whole text does.
 * @param {string | null | undefined} value
 */
export const base64ToBytes = (value) => {
  const raw = String(value || '');
  const encoded = /[ \t\n\f\r]/.test(raw) ? raw.replace(/[ \t\n\f\r]/g, '') : raw;
  const firstPadding = encoded.indexOf('=');
  if (firstPadding >= 0 && firstPadding < encoded.length - 2) {
    // Padding inside the text: the whole-text `atob` throws its own error.
    atob(encoded);
    throw new Error('Invalid base64 padding.');
  }
  const padding = encoded.length % 4 === 0
    ? (encoded.endsWith('==') ? 2 : encoded.endsWith('=') ? 1 : 0)
    : 0;
  const bytes = new Uint8Array(Math.max(0, Math.floor(encoded.length * 3 / 4) - padding));
  const chunkSize = 0x10000; // A multiple of four keeps each nonfinal chunk unpadded.
  let offset = 0;
  for (let start = 0; start < encoded.length; start += chunkSize) {
    const binary = atob(encoded.slice(start, start + chunkSize));
    if (offset + binary.length > bytes.length) throw new Error('Invalid base64 length.');
    for (let index = 0; index < binary.length; index += 1) bytes[offset + index] = binary.charCodeAt(index);
    offset += binary.length;
  }
  if (offset !== bytes.length) throw new Error('Invalid base64 length.');
  return bytes;
};

// Modern engines decode large adopted resources directly into their final
// buffer in bounded tasks, leaving room for paint and input without a
// resource-sized binary-string copy.
/**
 * @param {string | null | undefined} value
 * @returns {Promise<Uint8Array>}
 */
export const base64ToBytesInTasks = async (value) => {
  const raw = String(value || '');
  // Older engines keep their existing decoder and error behavior.
  if (typeof Uint8Array.prototype.setFromBase64 !== 'function') return base64ToBytes(raw);
  const encoded = /[ \t\n\f\r]/.test(raw)
    ? raw.replace(/[ \t\n\f\r]/g, '')
    : raw;
  const firstPadding = encoded.indexOf('=');
  if (firstPadding >= 0 && firstPadding < encoded.length - 2) {
    throw new Error('Invalid base64 padding.');
  }
  const padding = encoded.endsWith('==') ? 2 : encoded.endsWith('=') ? 1 : 0;
  const bytes = new Uint8Array(Math.max(0, Math.floor(encoded.length * 3 / 4) - padding));
  const chunkSize = 0x10000; // A multiple of four keeps each nonfinal chunk unpadded.
  let offset = 0;
  let deadline = performance.now() + 16;
  for (let start = 0; start < encoded.length; start += chunkSize) {
    const part = encoded.slice(start, start + chunkSize);
    const { read, written } = bytes.subarray(offset).setFromBase64(part);
    if (read !== part.length) throw new Error('Invalid base64 length.');
    offset += written;
    if (performance.now() >= deadline) {
      await new Promise((resolve) => setTimeout(resolve, 0));
      deadline = performance.now() + 16;
    }
  }
  if (offset !== bytes.byteLength) throw new Error('Invalid base64 length.');
  return bytes;
};

/**
 * @param {string | null | undefined} value
 * @param {number} byteLength
 * @returns {number | null}
 */
export const base64DecodedLastByte = (value, byteLength) => {
  if (!Number.isSafeInteger(byteLength) || byteLength <= 0) return null;
  const encoded = String(value || '');
  /** @type {number[]} */
  const terminalSextets = [];
  // The terminal byte depends only on the final two non-padding sextets and size modulo 3.
  for (let index = encoded.length - 1; index >= 0 && terminalSextets.length < 2; index -= 1) {
    const character = encoded[index];
    if (character === '=' || isBase64Whitespace(character)) continue;
    const sextet = BASE64_ALPHABET.indexOf(character);
    if (sextet < 0) return null;
    terminalSextets.push(sextet);
  }
  if (terminalSextets.length < 2) return null;

  const finalSextet = terminalSextets[0];
  const penultimateSextet = terminalSextets[1];
  if (byteLength % 3 === 1) {
    return (penultimateSextet << 2) | (finalSextet >> 4);
  }
  if (byteLength % 3 === 2) {
    return ((penultimateSextet & 0x0F) << 4) | (finalSextet >> 2);
  }
  return ((penultimateSextet & 0x03) << 6) | finalSextet;
};

/** @param {string | number | boolean | null | undefined} value */
export const textToBytes = (value) => new TextEncoder().encode(String(value));

/**
 * @param {Uint8Array<ArrayBuffer> | ArrayLike<number> | ArrayBuffer | null | undefined} value
 * @param {TextDecoderOptions} [options]
 */
export const bytesToText = (value, options) => (
  new TextDecoder('utf-8', options).decode(asBytes(value))
);

/** @param {string | number | boolean | null | undefined} value */
export const textToBase64 = (value) => bytesToBase64(textToBytes(value));

/**
 * The SHA-256 digest of `bytes`; Web Crypto copies them before it hashes.
 * @param {Uint8Array} bytes
 * @returns {Promise<ArrayBuffer>}
 */
export const sha256Digest = async (bytes) => {
  if (!globalThis.crypto?.subtle) {
    throw new Error('SHA-256 file identity requires Web Crypto.');
  }
  return globalThis.crypto.subtle.digest('SHA-256', /** @type {Uint8Array<ArrayBuffer>} */ (bytes));
};

/**
 * @param {Uint8Array<ArrayBuffer> | ArrayLike<number> | ArrayBuffer | null | undefined} value
 * @returns {Promise<string>}
 */
export const sha256Hex = async (value) => {
  const digest = await sha256Digest(asBytes(value));
  return Array.from(new Uint8Array(digest))
    .map((byte) => byte.toString(16).padStart(2, '0'))
    .join('');
};

// A gzip stream starts with these two bytes (RFC 1952), whatever the file is named.
export const GZIP_SIGNATURE = Object.freeze([0x1f, 0x8b]);

/** @param {Uint8Array} bytes */
export const isGzipBytes = (bytes) => GZIP_SIGNATURE.every((byte, index) => bytes[index] === byte);

/**
 * Hands each chunk of `stream` to `onChunk` until the stream ends, and rejects
 * with `overCap()` once the chunks pass `maxBytes`, before that chunk is handed
 * on or another is read. The one capped reader of an expanded stream: a gzip
 * input file (`expandGzipBytes`) and a Session (services/session-file.js).
 * @param {ReadableStream<Uint8Array>} stream
 * @param {number} maxBytes
 * @param {(chunk: Uint8Array) => void | Promise<void>} onChunk Awaited before the next chunk is read.
 * @param {() => Error} overCap The error past the cap.
 * @returns {Promise<void>}
 */
export const readCappedStream = async (stream, maxBytes, onChunk, overCap) => {
  const reader = stream.getReader();
  let byteCount = 0;
  try {
    while (true) {
      const { value, done } = await reader.read();
      if (done) return;
      byteCount += value.byteLength;
      if (byteCount > maxBytes) throw overCap();
      await onChunk(value);
    }
  } catch (error) {
    await reader.cancel(error).catch(() => {});
    throw error;
  } finally {
    reader.releaseLock();
  }
};

/**
 * Whether `bytes` are a gzip file this browser expands: without
 * `DecompressionStream`, a gzip input stays as it is (Python reads a gzip input
 * file as well, under its own cap, so the Worker still draws it).
 * @param {Uint8Array} bytes
 */
export const expandsAsGzip = (bytes) => isGzipBytes(bytes) && typeof DecompressionStream === 'function';

/**
 * Hands the bytes a gzip file holds to `onPiece` as they expand, in pieces of
 * `pieceBytes` (the last one shorter; one empty piece for an empty file), each
 * in a buffer of its own that `onPiece` may transfer. Expanding stops with
 * `overCap()` past `maxBytes`. The caller checks `expandsAsGzip` first.
 * @param {Uint8Array<ArrayBuffer>} bytes
 * @param {number} maxBytes
 * @param {() => Error} overCap
 * @param {number} pieceBytes
 * @param {(piece: Uint8Array<ArrayBuffer>) => Promise<void>} onPiece Awaited before more is expanded.
 * @returns {Promise<void>}
 */
export const expandGzipPieces = async (bytes, maxBytes, overCap, pieceBytes, onPiece) => {
  let piece = new Uint8Array(pieceBytes);
  let filled = 0;
  let total = 0;
  await readCappedStream(new Blob([bytes]).stream().pipeThrough(new DecompressionStream('gzip')), maxBytes, async (chunk) => {
    for (let offset = 0; offset < chunk.byteLength;) {
      const count = Math.min(pieceBytes - filled, chunk.byteLength - offset);
      piece.set(chunk.subarray(offset, offset + count), filled);
      filled += count;
      offset += count;
      total += count;
      if (filled === pieceBytes) {
        await onPiece(piece);
        piece = new Uint8Array(pieceBytes);
        filled = 0;
      }
    }
  }, overCap);
  if (filled > 0 || total === 0) await onPiece(piece.slice(0, filled));
};

/**
 * `pieces` joined in order into one buffer.
 * @param {readonly Uint8Array[]} pieces
 * @returns {Uint8Array<ArrayBuffer>}
 */
export const joinBytes = (pieces) => {
  const joined = new Uint8Array(pieces.reduce((length, piece) => length + piece.byteLength, 0));
  let offset = 0;
  for (const piece of pieces) {
    joined.set(piece, offset);
    offset += piece.byteLength;
  }
  return joined;
};
