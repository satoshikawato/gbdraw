// @ts-check
const BASE64_ALPHABET = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/';

const isBase64Whitespace = (character) => (
  character === ' '
  || character === '\t'
  || character === '\n'
  || character === '\f'
  || character === '\r'
);

export const asBytes = (value) => (
  value instanceof Uint8Array ? value : new Uint8Array(value || [])
);

export const bytesToBase64 = (value) => {
  const bytes = asBytes(value);
  const chunks = [];
  // Complete triplets let encoded chunks concatenate without interior padding.
  const chunkSize = 0x6000;
  for (let index = 0; index < bytes.length; index += chunkSize) {
    chunks.push(btoa(String.fromCharCode(...bytes.subarray(index, index + chunkSize))));
  }
  return chunks.join('');
};

export const base64ToBytes = (value) => {
  const binary = atob(String(value || ''));
  return Uint8Array.from(binary, (character) => character.charCodeAt(0));
};

// Modern engines decode large adopted resources directly into their final
// buffer in bounded tasks, leaving room for paint and input without a
// resource-sized binary-string copy.
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

export const base64DecodedLastByte = (value, byteLength) => {
  if (!Number.isSafeInteger(byteLength) || byteLength <= 0) return null;
  const encoded = String(value || '');
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

export const textToBytes = (value) => new TextEncoder().encode(String(value));

export const bytesToText = (value, options) => (
  new TextDecoder('utf-8', options).decode(asBytes(value))
);

export const textToBase64 = (value) => bytesToBase64(textToBytes(value));

export const sha256Hex = async (value) => {
  if (!globalThis.crypto?.subtle) {
    throw new Error('SHA-256 file identity requires Web Crypto.');
  }
  const bytes = asBytes(value);
  const digest = await globalThis.crypto.subtle.digest('SHA-256', bytes);
  return Array.from(new Uint8Array(digest))
    .map((byte) => byte.toString(16).padStart(2, '0'))
    .join('');
};
