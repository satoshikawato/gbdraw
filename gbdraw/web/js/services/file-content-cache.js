const fileByteCache = new WeakMap();

import {
  isSessionResourceFileView,
  readSessionResourceBytes,
  readSessionResourceText,
  sessionResourceSource
} from './session-resource-backing.js';

const asBytes = (value) => (
  value instanceof Uint8Array ? value : new Uint8Array(value || [])
);

const asCacheKey = (file) => (
  file && (typeof file === 'object' || typeof file === 'function') ? file : null
);

export const readFileBytes = async (file) => {
  const key = asCacheKey(file);
  const lazyBytes = key && isSessionResourceFileView(file)
    ? readSessionResourceBytes(file)
    : null;
  if (!key || (!lazyBytes && typeof file.arrayBuffer !== 'function')) {
    throw new TypeError('A File-like object with arrayBuffer() is required.');
  }
  let pending = fileByteCache.get(key);
  if (!pending) {
    pending = Promise.resolve()
      .then(() => lazyBytes || file.arrayBuffer())
      .then((buffer) => new Uint8Array(buffer));
    fileByteCache.set(key, pending);
    if (!lazyBytes) {
      pending.catch(() => {
        if (fileByteCache.get(key) === pending) fileByteCache.delete(key);
      });
    }
  }
  return pending;
};

export const readFileText = async (file) => {
  const lazyText = isSessionResourceFileView(file)
    ? readSessionResourceText(file)
    : null;
  if (lazyText) return lazyText;
  if (typeof file?.arrayBuffer === 'function') {
    return bytesToText(await readFileBytes(file));
  }
  if (typeof file?.text === 'function') return file.text();
  throw new TypeError('A File-like object with arrayBuffer() or text() is required.');
};

export const getSessionResourceSource = (file) => sessionResourceSource(file);

export const cloneFileBytesForTransfer = async (file) => {
  const bytes = await readFileBytes(file);
  return bytes.slice().buffer;
};

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
