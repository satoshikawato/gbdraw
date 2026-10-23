// R15 (D-15-60): until a project holds several drawings of one mode (B), a
// drawing is resolved by its mode only through `state.drawingOfMode`
// (state.js), and only the modules `ID_BY_MODE_OWNERS` lists name a drawing ID
// by a mode. B moves each call to `activeDrawing()` or a `DrawingId`; the cap
// may only fall.
//
// Both literals are registered design-rule allowlists
// (tools/web-design-rule-guards.json, R15): `DRAWING_OF_MODE_CAP` (kind
// `count`) and `ID_BY_MODE_OWNERS` (kind `set`). Raising the cap or adding an
// owner is an authority-only change. This test holds a count literal, so it
// imports only Node built-ins (WEB_CHANGE_POLICY.md "Design-rule co-change").
import assert from 'node:assert/strict';
import { readdirSync, readFileSync } from 'node:fs';
import test from 'node:test';

const JS = new URL('../../gbdraw/web/js/', import.meta.url);
const webModules = (url = JS, prefix = '') => readdirSync(url, { withFileTypes: true }).flatMap((entry) => {
  const path = `${prefix}${entry.name}`;
  if (entry.isDirectory()) return webModules(new URL(`${entry.name}/`, url), `${path}/`);
  return entry.name.endsWith('.js') ? [path] : [];
});

const DRAWING_OF_MODE_CAP = 115;
const ID_BY_MODE_OWNERS = Object.freeze(['state.js', 'services/drawing-view.js']);
test('a drawing is resolved by mode only through state.drawingOfMode, within a falling cap', () => {
  const sources = webModules().map((file) => [file, readFileSync(new URL(file, JS), 'utf8')]);
  const calls = sources.reduce((sum, [, text]) => sum + (text.match(/\bdrawingOfMode\(/g) || []).length, 0);
  assert.ok(calls <= DRAWING_OF_MODE_CAP,
    `${calls} drawingOfMode( calls, above the cap ${DRAWING_OF_MODE_CAP}: resolve the shown drawing or a drawing ID instead`);
  assert.equal(calls, DRAWING_OF_MODE_CAP, `the drawingOfMode( calls fell to ${calls}: lower DRAWING_OF_MODE_CAP to it`);
  const byMode = sources.filter(([file, text]) => !ID_BY_MODE_OWNERS.includes(file) && /\bdrawingIdOfMode\(/.test(text))
    .map(([file]) => file);
  assert.deepEqual(byMode, [], 'name a drawing of a mode through state.drawingOfMode');
});
