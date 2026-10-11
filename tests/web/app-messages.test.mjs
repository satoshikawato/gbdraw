import assert from 'node:assert/strict';
import { mkdtemp, readFile, readdir, writeFile } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join, relative } from 'node:path';
import { pathToFileURL } from 'node:url';
import { test } from 'node:test';

const repoRoot = process.cwd();
const webJsRoot = join(repoRoot, 'gbdraw', 'web', 'js');

const listModules = async (dir) => (await Promise.all(
  (await readdir(dir, { withFileTypes: true })).map((entry) => {
    const path = join(dir, entry.name);
    if (entry.isDirectory()) return listModules(path);
    return entry.name.endsWith('.js') ? [path] : [];
  })
)).flat();

// A call of the browser's own dialog, with or without `window.`/`globalThis.`;
// method calls such as `dialog.confirm(` and names such as `askConfirm(` do not count.
const NATIVE_DIALOG_CALL = /(?<![\w$.])(?:(?:window|globalThis)\.)?(alert|confirm)\s*\(/;

test('UX-A07, UX-C12: the Web app asks and reports in its own dialog and notice, never alert() or confirm()', async () => {
  const offenders = [];
  for (const path of await listModules(webJsRoot)) {
    const lines = (await readFile(path, 'utf8')).split('\n');
    lines.forEach((line, index) => {
      const code = line.replace(/\/\/.*$/, '').replace(/^\s*\*.*$/, '');
      if (NATIVE_DIALOG_CALL.test(code)) offenders.push(`${relative(repoRoot, path)}:${index + 1}: ${line.trim()}`);
    });
  }
  assert.deepEqual(offenders, []);
});

const tempDir = await mkdtemp(join(tmpdir(), 'gbdraw-app-messages-'));
await writeFile(join(tempDir, 'package.json'), '{"type":"module"}\n', 'utf8');
await writeFile(join(tempDir, 'app-messages.js'), await readFile(join(webJsRoot, 'app', 'app-messages.js'), 'utf8'), 'utf8');
const { createChoicePrompt, createNoticeBanner } = await import(pathToFileURL(join(tempDir, 'app-messages.js')));
const reactive = (value) => value;

test('the choice prompt resolves to the chosen action, and a new question cancels the open one', async () => {
  const { prompt, ask, resolve } = createChoicePrompt({ reactive });
  const first = ask({ title: 'Reset all settings?', actions: [{ id: 'reset', label: 'Reset all settings', tone: 'danger' }], cancelLabel: 'Keep settings' });
  assert.equal(prompt.show, true);
  assert.equal(prompt.cancelLabel, 'Keep settings');
  const second = ask({ title: 'Add this color rule?', actions: [{ id: 'add', label: 'Add rule' }] });
  assert.equal(await first, 'cancel');
  assert.deepEqual(prompt.actions, [{ tone: 'secondary', id: 'add', label: 'Add rule' }]);
  assert.equal(prompt.cancelLabel, 'Cancel');
  resolve('add');
  assert.equal(await second, 'add');
  assert.equal(prompt.show, false);
});

test('the notice banner shows warnings as a list and is cleared by dismiss', () => {
  const { notice, show, dismiss } = createNoticeBanner({ reactive });
  show({ text: 'Session loaded: a.gbdraw-session.json', items: ['', 'One warning'] });
  assert.deepEqual({ ...notice }, { text: 'Session loaded: a.gbdraw-session.json', items: ['One warning'], tone: 'warning' });
  show({ text: 'Session saved: a.gbdraw-session.json.' });
  assert.equal(notice.tone, 'info');
  assert.deepEqual(notice.items, []);
  dismiss();
  assert.equal(notice.text, '');
});
