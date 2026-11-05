// The Web's Figure size (issue #1427): what a drawing's form holds and how a stored length reads into it
// (services/figure-size-form.js), the field (app/figure-size-field.js), the size line under it, the root size it
// reads (services/svg-dimensions.js), the words of a refused size and of the larger-figure remedy (one wording
// owner each, utils/error-normalization.js and services/label-fit.js). The request and Session channel is tested in
// session-request, mode-scoped-session and run-info.
import assert from 'node:assert/strict';
import { execFileSync } from 'node:child_process';
import test from 'node:test';

import {
  DEFAULT_FIGURE_LENGTH_UNIT, FIGURE_LENGTH_UNITS, figureLengthAmountText, figureLengthFormValue, figureLengthInputValue,
  figureGrowth, figureLengthUnitOf, requestFigureSizeLeaves, requestedLengthText, shownFigureSize, splitFigureLength
} from '../../gbdraw/web/js/services/figure-size-form.js';
import { createFigureSizeField } from '../../gbdraw/web/js/app/figure-size-field.js';
import { getSvgDimensions } from '../../gbdraw/web/js/services/svg-dimensions.js';
import { createDefaultForm } from '../../gbdraw/web/js/services/session-active-config-contract.js';
import { normalizeUserFacingError } from '../../gbdraw/web/js/utils/error-normalization.js';
import { labelCapacityRemedyText } from '../../gbdraw/web/js/services/label-fit.js';
import { defaultSettingsContainers, installStoredSettings } from '../../gbdraw/web/js/services/config-field-readers.js';
import { withDrawings } from './helpers/drawing-state.mjs';

test('the units are Python\'s, and an empty side shows px (the Owner, 2026-11-04)', () => {
  assert.deepEqual([...FIGURE_LENGTH_UNITS], ['mm', 'cm', 'in', 'px']);
  assert.equal(DEFAULT_FIGURE_LENGTH_UNIT, 'px');
  const form = createDefaultForm();
  assert.deepEqual([form.figure_width, form.figure_height], [null, null], 'a new drawing is Auto');
});

test('a stored length reads as the form holds it; any other text is Python\'s to read, else lost', () => {
  assert.deepEqual(splitFigureLength('180mm'), { amount: '180', unit: 'mm' });
  assert.deepEqual(splitFigureLength('1.5e2px'), { amount: '1.5e2', unit: 'px' });
  for (const value of ['180 mm', '180MM', 'mm', '180', 180, null, '1,5in']) assert.equal(splitFigureLength(value), null, String(value));
  for (const raw of [undefined, null, 'auto', 'AUTO']) assert.deepEqual(figureLengthFormValue('width', raw), { value: null, lost: false });
  assert.deepEqual(figureLengthFormValue('height', '5in'), { value: '5in', lost: false });
  const readings = new Map([[JSON.stringify(['canvas.figure_width', '180 MM']), '180mm']]);
  assert.deepEqual(figureLengthFormValue('width', '180 MM', readings), { value: '180mm', lost: false });
  assert.deepEqual(figureLengthFormValue('height', '180 MM', readings), { value: null, lost: true }, 'a reading is per leaf');
  assert.deepEqual(figureLengthFormValue('width', 'tall'), { value: null, lost: true });
});

test('a stored draft\'s side reads as the field holds it; one the field cannot show is Auto, with its row', () => {
  for (const mode of /** @type {const} */ (['circular', 'linear'])) {
    const read = (form) => {
      const targets = defaultSettingsContainers(mode);
      const { losses } = installStoredSettings(targets, { form }, mode);
      return { form: [targets.form.figure_width, targets.form.figure_height], losses: losses.filter(({ at }) => at.startsWith('form.figure_')) };
    };
    assert.deepEqual(read({ figure_width: '180mm', figure_height: null }), { form: ['180mm', null], losses: [] });
    assert.deepEqual(read({}), { form: [null, null], losses: [] }, 'a draft before the field is Auto');
    const lost = read({ figure_width: '180 MM', figure_height: 'auto' });
    assert.deepEqual(lost.form, [null, null]);
    assert.deepEqual(lost.losses.map(({ at, stored }) => [at, stored]), [['form.figure_width', '180 MM']]);
  }
});

test('a side writes its amount as typed and its unit; Python reads and checks it', () => {
  assert.equal(figureLengthInputValue('', 'mm'), null, 'empty is Auto');
  assert.equal(figureLengthInputValue('  ', 'px'), null);
  assert.equal(figureLengthInputValue(' 180 ', 'mm'), '180mm');
  assert.equal(figureLengthInputValue('0', 'in'), '0in', 'Python refuses it on Generate');
  assert.equal(figureLengthUnitOf(null), 'px', 'an Auto value shows px');
  assert.equal(figureLengthUnitOf('5in'), 'in');
  assert.equal(figureLengthAmountText('180mm'), '180');
  assert.equal(figureLengthAmountText(null), '');
});

test('the request leaves are read from the override, else from a CLI request\'s full config', () => {
  const request = (configOverrides, config = null) => ({ diagramOptions: { configOverrides, config } });
  assert.deepEqual(requestFigureSizeLeaves(request({ 'canvas.figure_width': '180mm' })), { width: '180mm', height: null });
  assert.deepEqual(requestFigureSizeLeaves(request({}, { canvas: { figure_width: 'auto', figure_height: '5in' } })),
    { width: null, height: '5in' });
  assert.deepEqual(requestFigureSizeLeaves(request(null)), { width: null, height: null });
});

test('the size line names the requested unit first, else px then mm', () => {
  assert.deepEqual(shownFigureSize({ widthPx: 1400, heightPx: 1450, requested: { width: null, height: null } }),
    { first: '1400 × 1450 px', second: '370.4 × 383.6 mm' });
  assert.deepEqual(shownFigureSize({ widthPx: 180 * 96 / 25.4, heightPx: 541.2, requested: { width: '180mm', height: null } }),
    { first: '180 × 143.2 mm', second: '680 × 541 px' });
  assert.deepEqual(shownFigureSize({ widthPx: 720, heightPx: 480, requested: { width: null, height: '5in' } }),
    { first: '7.5 × 5 in', second: '720 × 480 px' });
  assert.deepEqual(shownFigureSize({ widthPx: 1000, heightPx: 500, requested: { width: '1000px', height: null } }),
    { first: '1000 × 500 px', second: '264.6 × 132.3 mm' });
});

test('the root size is read in px: a root in a requested unit reads its px viewBox (OV-2511)', () => {
  const root = (attributes) => ({ getAttribute: (name) => attributes[name] ?? null });
  assert.deepEqual(getSvgDimensions(root({ width: '1400', height: '1450.5' })), { width: 1400, height: 1450.5 });
  assert.deepEqual(getSvgDimensions(root({ width: '1400px', height: '1450px' })), { width: 1400, height: 1450 });
  assert.deepEqual(getSvgDimensions(root({ width: '180mm', height: '26.5327mm', viewBox: '0.0 0.0 680.3149606299213 100.28125' })),
    { width: 680.3149606299213, height: 100.28125 });
  assert.equal(getSvgDimensions(root({ width: '180mm' })), null);
  assert.equal(getSvgDimensions(null), null);
});

test('the field writes each side, keeps its unit as view state, and shows a value from elsewhere in its unit', () => {
  const drawing = { id: 'circular', mode: 'circular', form: createDefaultForm() };
  const store = withDrawings(drawing);
  const watchers = [];
  let shown = null;
  const field = createFigureSizeField({ state: store, watch: (_source, callback) => watchers.push(callback), shownSize: () => shown });
  const changed = () => watchers.forEach((callback) => callback());
  assert.equal(field.figureWidthUnit.value, 'px');
  assert.equal(field.figureWidthAmount.value, '');
  assert.equal(field.figureSizeLine.value, '', 'no line while the drawing shows no Result');
  // A unit change with no amount writes nothing.
  field.changeFigureSize({ target: { value: 'mm' } }, 'width', 'unit');
  changed();
  assert.equal(drawing.form.figure_width, null);
  assert.equal(field.figureWidthUnit.value, 'mm');
  field.changeFigureSize({ target: { value: '180' } }, 'width', 'amount');
  changed();
  assert.equal(drawing.form.figure_width, '180mm');
  // A unit change keeps the typed amount.
  field.changeFigureSize({ target: { value: 'in' } }, 'width', 'unit');
  changed();
  assert.equal(drawing.form.figure_width, '180in');
  // An Undo or a Load writes a value: it shows in its own unit.
  drawing.form.figure_height = '120cm';
  changed();
  assert.deepEqual([field.figureHeightUnit.value, field.figureHeightAmount.value], ['cm', '120']);
  field.changeFigureSize({ target: { value: '' } }, 'width', 'amount');
  changed();
  assert.equal(drawing.form.figure_width, null);
  assert.equal(field.figureWidthUnit.value, 'in', 'Auto keeps the unit');
  shown = { first: '180 × 143.2 mm', second: '680 × 541 px', larger: null };
  assert.equal(field.figureSizeLine.value, 'Figure: 180 × 143.2 mm (680 × 541 px)');
  assert.equal(field.figureSizeLarger.value, '');
  // Section 6: a canvas grown past the requested box says so. Python fits the figure to Legend edits and Canvas
  // padding (PD-OI-114 rev 2), so when they still leave it larger, the fit was refused: the note above says why
  // (review S3, D-60).
  shown = { first: '190 × 143.2 mm', second: '718 × 541 px', larger: { sides: [{ side: 'width', length: '180 mm' }], causes: ['legend-edits', 'canvas-padding'] } };
  assert.equal(field.figureSizeLine.value, 'Figure: 190 × 143.2 mm (718 × 541 px)');
  assert.equal(field.figureSizeLarger.value, 'Larger than the requested 180 mm width: the edit could not be fitted (see the note above the preview).');
  shown = { ...shown, larger: { sides: [{ side: 'width', length: '180 mm' }], causes: ['canvas-padding'] } };
  assert.equal(field.figureSizeLarger.value, 'Larger than the requested 180 mm width: the edit could not be fitted (see the note above the preview).');
  shown = { ...shown, larger: { sides: [{ side: 'width', length: '180 mm' }], causes: ['legend-edits', 'moved-items'] } };
  assert.equal(field.figureSizeLarger.value, 'Larger than the requested 180 mm width: the edit could not be fitted (see the note above the preview).');
  // A move is never fitted: the sentence names it.
  shown = { ...shown, larger: { sides: [{ side: 'height', length: '120 mm' }], causes: ['moved-items'] } };
  assert.equal(field.figureSizeLarger.value, 'Larger than the requested 120 mm height after moves of the Legend, title, or scale, which the fit does not include.');
});

// Review finding 1: the unit selects belong to the one field every drawing shares. A value that arrives from outside
// the field (a drawing switch, Reset Settings, Session Load, Undo) shows in its own unit, and an Auto value in px (the
// receipt: an empty side starts at px); only a side the field itself cleared keeps the unit the user chose.
test('a side that becomes Auto from outside the field shows px; a side the field cleared keeps its unit', () => {
  const a = { id: 'circular', mode: 'circular', form: createDefaultForm() };
  const b = { id: 'circular-2', mode: 'circular', form: createDefaultForm() };
  let shownDrawing = a;
  const watchers = [];
  const field = createFigureSizeField({ state: { activeDrawing: () => shownDrawing }, watch: (_source, callback) => watchers.push(callback), shownSize: () => null });
  const changed = () => watchers.forEach((callback) => callback());
  const shown = () => [field.figureWidthUnit.value, field.figureWidthAmount.value];
  const type = (value, part = 'amount') => { field.changeFigureSize({ target: { value } }, 'width', part); changed(); };
  type('mm', 'unit');
  type('85');
  assert.equal(a.form.figure_width, '85mm');
  // A switch to drawing B, whose width is Auto: px, and typing 1200 writes 1200 px.
  shownDrawing = b;
  changed();
  assert.deepEqual(shown(), ['px', ''], 'a drawing switch to an Auto side');
  type('1200');
  assert.equal(b.form.figure_width, '1200px');
  shownDrawing = a;
  changed();
  assert.deepEqual(shown(), ['mm', '85'], 'back in A, its own value in its own unit');
  // Reset Settings (the form back to its defaults), and Undo or Load of an Auto value.
  Object.assign(a.form, createDefaultForm());
  changed();
  assert.deepEqual(shown(), ['px', ''], 'Reset Settings');
  a.form.figure_width = '20cm';
  changed();
  assert.deepEqual(shown(), ['cm', '20'], 'a Load or Undo of a value');
  a.form.figure_width = null;
  changed();
  assert.deepEqual(shown(), ['px', ''], 'a Load or Undo of Auto');
  // The field clears a side: it keeps the unit chosen, also when the same value comes back by itself.
  type('in', 'unit');
  type('7');
  type('');
  assert.equal(a.form.figure_width, null);
  assert.deepEqual(shown(), ['in', ''], 'a side the field cleared keeps its unit');
  changed();
  assert.deepEqual(shown(), ['in', '']);
  type('7');
  assert.equal(a.form.figure_width, '7in');
});

test('a shown Result is larger than requested only past the box, by its Web-only causes', () => {
  const requested = { widthPx: 680, heightPx: null };
  const box = { width: 680, height: 541 };
  assert.equal(figureGrowth({ shown: box, base: box, requested, padded: true, moved: true }), null, 'within the box');
  assert.equal(figureGrowth({ shown: { width: 680.4, height: 900 }, base: box, requested, padded: false, moved: false }), null,
    'an Auto side and half a px are not larger');
  assert.deepEqual(figureGrowth({ shown: { width: 720, height: 541 }, base: { width: 700, height: 541 }, requested, padded: true, moved: false }),
    { sides: ['width'], causes: ['legend-edits', 'canvas-padding'] });
  assert.deepEqual(figureGrowth({ shown: { width: 720, height: 600 }, base: box, requested: { widthPx: 680, heightPx: 541 }, padded: false, moved: true }),
    { sides: ['width', 'height'], causes: ['moved-items'] });
  assert.equal(requestedLengthText('180mm'), '180 mm');
  assert.equal(requestedLengthText('180 MM'), null);
});

// What Python sends for a refused size, through the Web error adapter (gbdraw/web_support/error_adapter.py).
const pythonRefusals = () => JSON.parse(execFileSync(process.env.PYTHON || 'python', ['-c', `
import json, warnings
from pathlib import Path
import gbdraw
from gbdraw.exceptions import GbdrawError
from gbdraw.web_support.error_adapter import serialize_web_error
warnings.simplefilter("ignore")
records = gbdraw.read_genbank(Path("tests/test_inputs/HmmtDNA.gbk"))
out = []
for draw, options in ((gbdraw.draw_linear, gbdraw.LinearOptions(figure_width="40mm")),
                      (gbdraw.draw_linear, gbdraw.LinearOptions(figure_width="180mm", figure_height="10mm")),
                      (gbdraw.draw_circular, gbdraw.CircularOptions(figure_width="100mm"))):
    try:
        draw(records, options=options)
    except GbdrawError as error:
        out.append(serialize_web_error(error, operation="generate", stage="render"))
print(json.dumps(out))
`], { encoding: 'utf8' }));

test('a refused size names the parts that do not shrink and the way out, with Python\'s numbers', () => {
  const [narrow, short, circular] = pythonRefusals().map(normalizeUserFacingError);
  const { minimum, primary, legend, padding } = narrow.context;
  assert.equal(narrow.summary, 'An input value is invalid. The figure cannot be 40 mm wide:'
    + ` the parts that do not shrink need ${minimum} (record names and labels past the sequence ends ${primary},`
    + ` legend ${legend}, padding ${padding}). Set Figure width to at least ${minimum}, place the Legend at the top or`
    + ' bottom, or leave it empty (Auto).');
  const { minimum: tall, primary: records } = short.context;
  assert.ok(short.summary.startsWith(`An input value is invalid. The figure cannot be 10 mm tall at 180 mm wide: it needs ${tall} (records ${records}, `), short.summary);
  assert.ok(short.summary.endsWith(`Set Figure height to at least ${tall}, or leave it empty (Auto).`), short.summary);
  assert.match(circular.summary, /\(map at its smallest [0-9.]+ mm, legend [0-9.]+ mm, padding [0-9.]+ mm\)/);
  // Without its lengths a refusal says what to change.
  const bare = normalizeUserFacingError({ code: 'INPUT_INVALID', operation: 'generate', stage: 'render',
    context: { reason: 'FIGURE_TOO_NARROW', configPath: 'canvas.figure_width', requested: '40mm <b>' } });
  assert.ok(bare.summary.endsWith('The figure cannot be this narrow. Set a larger Figure width, or leave it empty (Auto).'), bare.summary);
  assert.equal(bare.context.requested, undefined, 'only a bounded length is kept');
});

test('the larger-figure remedy names Python\'s estimate for each requested side', () => {
  assert.deepEqual(labelCapacityRemedyText('figure_size', null), { text: 'Make the figure larger (Figure size in Layout).', params: {} });
  assert.deepEqual(labelCapacityRemedyText('figure_size', { width: '240.3 mm', height: null }),
    { text: 'Make the figure larger: Figure width about {width} (Figure size in Layout).', params: { width: '240.3 mm' } });
  assert.deepEqual(labelCapacityRemedyText('figure_size', { width: '240.3 mm', height: '9.5 in' }).params,
    { width: '240.3 mm', height: '9.5 in' });
  assert.deepEqual(labelCapacityRemedyText('filter', { width: '240.3 mm', height: null }).params, {});
});

// D-24: the plain SVG of a Result drawn to a requested size states the root size Python's `export_root_size` gives
// for the shown root's px size, on the export copy only; a Result drawn at Auto (no port) exports its root as shown.
test('a sized Result exports its SVG root in the requested unit, and the page root stays px', async () => {
  const { downloadSVG } = await import('../../gbdraw/web/js/services/export.js');
  const svg = {
    attributes: new Map([['width', '680.3149606299213px'], ['height', '541px'], ['viewBox', '0 0 680.3149606299213 541']]),
    getAttribute(name) { return this.attributes.get(name) ?? null; },
    setAttribute(name, value) { this.attributes.set(name, String(value)); },
    removeAttribute(name) { this.attributes.delete(name); },
    querySelectorAll: () => [],
    classList: { toggle() {}, remove() {}, length: 0 }
  };
  const page = svg;
  const serialized = [];
  const saved = { document: globalThis.document, XMLSerializer: globalThis.XMLSerializer };
  globalThis.document = { createElement: (tag) => {
    if (tag === 'template') return { set innerHTML(_markup) {}, content: { querySelector: () => ({ ...svg, attributes: new Map(svg.attributes) }) } };
    throw new Error(`the export created a ${tag}`);
  } };
  globalThis.XMLSerializer = class { serializeToString(root) { serialized.push(Object.fromEntries(root.attributes)); return '<svg/>'; } };
  try {
    const asked = [];
    const sizedSvgRoot = async (size) => { asked.push(size); return { width: '180mm', height: '143.1396mm' }; };
    const snapshot = { name: 'sized.svg', svgContent: '<svg/>' };
    await assert.rejects(downloadSVG(snapshot, { sizedSvgRoot }), /created a a/);
    assert.deepEqual(asked, [{ width: 680.3149606299213, height: 541 }]);
    assert.equal(serialized.at(-1).width, '180mm');
    assert.equal(serialized.at(-1).height, '143.1396mm');
    assert.equal(serialized.at(-1).viewBox, '0 0 680.3149606299213 541');
    assert.equal(page.getAttribute('width'), '680.3149606299213px', 'the page root stays px');
    serialized.length = 0;
    await assert.rejects(downloadSVG(snapshot, { sizedSvgRoot: null }), /created a a/);
    assert.equal(serialized.at(-1).width, '680.3149606299213px', 'an Auto Result exports its root as shown');
  } finally {
    Object.assign(globalThis, saved);
  }
});
