// What the app keeps in browser storage: the display of the settings pane and
// the UI language (services/browser-preferences.js, PD-OI-093 and PD-OI-106).
// A storage that is absent, throws, or holds another value gives the default,
// a failed write is silent, and no other module touches browser storage.
import assert from 'node:assert/strict';
import { readFileSync, readdirSync } from 'node:fs';
import { join, relative } from 'node:path';
import test from 'node:test';

import {
  readSettingsPaneDisplay, readUiLanguage, writeSettingsPaneDisplay, writeUiLanguage
} from '../../gbdraw/web/js/services/browser-preferences.js';

const KEY = 'gbdraw.settings-pane-display';
const LANGUAGE_KEY = 'gbdraw.ui-language';
const storageWith = (entries = {}) => {
  const values = new Map(Object.entries(entries));
  return {
    values,
    getItem: (key) => (values.has(key) ? values.get(key) : null),
    setItem: (key, value) => { values.set(key, String(value)); }
  };
};
const failing = (message) => () => { throw new DOMException(message, 'SecurityError'); };

test('the stored display is read under one key', () => {
  assert.equal(readSettingsPaneDisplay(() => storageWith({ [KEY]: 'rail' })), 'rail');
  assert.equal(readSettingsPaneDisplay(() => storageWith({ [KEY]: 'stacked' })), 'stacked');
  assert.equal(readSettingsPaneDisplay(() => storageWith({ 'settings-pane-display': 'rail' })), 'stacked');
});

test('no stored value, another value, and a storage that throws give the stacked display', () => {
  assert.equal(readSettingsPaneDisplay(() => storageWith()), 'stacked');
  for (const value of ['', 'Rail', 'rail ', 'icons', 'null', '{"display":"rail"}']) {
    assert.equal(readSettingsPaneDisplay(() => storageWith({ [KEY]: value })), 'stacked', JSON.stringify(value));
  }
  // The storage object itself is refused (blocked site data), or its read is.
  assert.equal(readSettingsPaneDisplay(failing('The operation is insecure.')), 'stacked');
  assert.equal(readSettingsPaneDisplay(() => ({ getItem: failing('read refused') })), 'stacked');
  assert.equal(readSettingsPaneDisplay(() => undefined), 'stacked');
});

test('a write stores the display, and a write that throws changes nothing and reports nothing', () => {
  const storage = storageWith();
  writeSettingsPaneDisplay('rail', () => storage);
  assert.deepEqual([...storage.values], [[KEY, 'rail']]);
  writeSettingsPaneDisplay('stacked', () => storage);
  assert.deepEqual([...storage.values], [[KEY, 'stacked']]);
  assert.equal(writeSettingsPaneDisplay('rail', failing('The operation is insecure.')), undefined);
  assert.equal(writeSettingsPaneDisplay('rail', () => ({ setItem: failing('QuotaExceededError') })), undefined);
  assert.equal(writeSettingsPaneDisplay('rail', () => undefined), undefined);
});

test('the UI language is read under its own key, and anything else gives English', () => {
  for (const tag of ['en', 'ja', 'zh-Hans', 'zh-Hant']) {
    assert.equal(readUiLanguage(() => storageWith({ [LANGUAGE_KEY]: tag })), tag);
  }
  // A first visit has nothing stored: English, whatever the browser's language.
  assert.equal(readUiLanguage(() => storageWith()), 'en');
  assert.equal(readUiLanguage(() => storageWith({ [KEY]: 'ja' })), 'en', 'the display key is not the language');
  for (const value of ['', 'JA', 'ja-JP', 'zh', 'zh-CN', 'zh-TW', 'fr', ' ja', '{"tag":"ja"}']) {
    assert.equal(readUiLanguage(() => storageWith({ [LANGUAGE_KEY]: value })), 'en', JSON.stringify(value));
  }
  assert.equal(readUiLanguage(failing('The operation is insecure.')), 'en');
  assert.equal(readUiLanguage(() => ({ getItem: failing('read refused') })), 'en');
  assert.equal(readUiLanguage(() => undefined), 'en');
});

test('a written UI language is kept beside the display, and a write that throws is silent', () => {
  const storage = storageWith({ [KEY]: 'rail' });
  writeUiLanguage('zh-Hant', () => storage);
  assert.deepEqual([...storage.values], [[KEY, 'rail'], [LANGUAGE_KEY, 'zh-Hant']]);
  writeUiLanguage('en', () => storage);
  assert.equal(storage.values.get(LANGUAGE_KEY), 'en');
  assert.equal(readSettingsPaneDisplay(() => storage), 'rail');
  assert.equal(writeUiLanguage('ja', failing('The operation is insecure.')), undefined);
  assert.equal(writeUiLanguage('ja', () => ({ setItem: failing('QuotaExceededError') })), undefined);
  assert.equal(writeUiLanguage('ja', () => undefined), undefined);
});

// The app's own modules and its page; vendored libraries are not the app's code.
const WEB = join(process.cwd(), 'gbdraw', 'web');
const sources = (directory) => readdirSync(directory, { withFileTypes: true }).flatMap((entry) => {
  const path = join(directory, entry.name);
  if (entry.isDirectory()) return sources(path);
  return entry.name.endsWith('.js') ? [path] : [];
});

test('browser storage is touched in the preferences module only', () => {
  const OWNER = 'js/services/browser-preferences.js';
  const files = [...sources(join(WEB, 'js')), join(WEB, 'index.html')];
  assert.ok(files.length > 200, `${files.length} files scanned`);
  const users = files.filter((path) => /\b(?:localStorage|sessionStorage|indexedDB)\b|document\s*\.\s*cookie/
    .test(readFileSync(path, 'utf8'))).map((path) => relative(WEB, path).replaceAll('\\', '/'));
  assert.deepEqual(users, [OWNER]);
  // The owner names local storage and nothing else.
  const owner = readFileSync(join(WEB, OWNER), 'utf8');
  assert.equal(/\b(?:sessionStorage|indexedDB)\b|document\s*\.\s*cookie/.test(owner), false);
  assert.equal(owner.split(`'${KEY}'`).length - 1, 1, 'one display key');
  assert.equal(owner.split(`'${LANGUAGE_KEY}'`).length - 1, 1, 'one language key');
  assert.deepEqual(owner.match(/'gbdraw\.[^']*'/g), [`'${KEY}'`, `'${LANGUAGE_KEY}'`], 'no other key');
});
