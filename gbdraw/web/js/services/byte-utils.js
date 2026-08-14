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
  let binary = '';
  const chunkSize = 0x8000;
  for (let index = 0; index < bytes.length; index += chunkSize) {
    binary += String.fromCharCode(...bytes.subarray(index, index + chunkSize));
  }
  return btoa(binary);
};

export const base64ToBytes = (value) => {
  const binary = atob(String(value || ''));
  return Uint8Array.from(binary, (character) => character.charCodeAt(0));
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
