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

// A call of the browser's own dialog, or any use of it through `window.` or
// `globalThis.` (`confirmFn = globalThis.confirm` calls it later); method
// calls such as `dialog.confirm(` and names such as `askConfirm(` do not count.
const NATIVE_DIALOG_CALL = /(?<![\w$.])(?:(?:window|globalThis|self)\.(?:alert|confirm)\b|(?:alert|confirm)\s*\()/;

test('the Web app asks and reports in its own dialog and notice, never alert() or confirm()', async () => {
  const offenders = [];
  for (const path of await listModules(webJsRoot)) {
    const lines = (await readFile(path, 'utf8')).split('\n');
    lines.forEach((line, index) => {
      const code = line.replace(/\/\/.*$/, '').replace(/^\s*\*.*$/, '');
      if (NATIVE_DIALOG_CALL.test(code)) offenders.push(`${relative(repoRoot, path)}:${index + 1}: ${line.trim()}`);
    });
  }
  assert.deepEqual(offenders, []);
  // The pattern finds the indirect use that the base code had.
  for (const code of ['alert("x")', 'window.confirm("x")', 'const confirmFn = globalThis.confirm;', 'self.alert']) {
    assert.equal(NATIVE_DIALOG_CALL.test(code), true, code);
  }
  for (const code of ['dialog.confirm(x)', 'askConfirm(x)', 'confirmLargeDownload(x)', 'window.confirmed']) {
    assert.equal(NATIVE_DIALOG_CALL.test(code), false, code);
  }
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
  assert.deepEqual({ ...notice }, {
    text: 'Session loaded: a.gbdraw-session.json', items: ['One warning'], tone: 'warning',
    announcement: 'Session loaded: a.gbdraw-session.json One warning'
  });
  show({ text: 'Session saved: a.gbdraw-session.json.' });
  assert.equal(notice.tone, 'info');
  assert.deepEqual(notice.items, []);
  dismiss();
  assert.equal(notice.text, '');
  assert.equal(notice.announcement, '');
});

// A live region speaks only a change of its text, so a notice that repeats
// the shown one changes the announced text and keeps the visible one.
test('a repeated notice is announced again', () => {
  const { notice, show } = createNoticeBanner({ reactive });
  const spoken = [];
  for (let count = 0; count < 3; count += 1) {
    show({ text: 'Label table saved.' });
    spoken.push(notice.announcement);
  }
  assert.equal(notice.text, 'Label table saved.');
  assert.deepEqual(spoken.map((entry) => entry.trim()), Array(3).fill('Label table saved.'));
  assert.notEqual(spoken[0], spoken[1]);
  assert.notEqual(spoken[1], spoken[2]);
});
