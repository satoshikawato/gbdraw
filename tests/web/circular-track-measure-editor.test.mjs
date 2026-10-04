import assert from 'node:assert/strict';
import { readFile } from 'node:fs/promises';
import { test } from 'node:test';
import {
  readCircularMeasure, writeCircularMeasureValue,
  changeCircularMeasureUnit, projectCircularMeasureDraft
} from '../../gbdraw/web/js/app/circular-track-slots/measure-editor.js';
import {
  parseOptionalCircularScalar, validateCustomTrackPlan
} from '../../gbdraw/web/js/app/track-slot-validation.js';
import { buildCircularTrackSlotPayload } from '../../gbdraw/web/js/app/circular-track-slots.js';
import {
  createDefaultAdv, createDefaultForm, validateCurrentWriterActiveConfig
} from '../../gbdraw/web/js/services/session-active-config-contract.js';
import { projectCanonicalSessionRequest } from '../../gbdraw/web/js/services/session-request.js';

const fixtures = JSON.parse(await readFile(new URL(
  '../../docs/internal/issue-619-implementation-plan-20260927/SESSION_RESULTS/scalar-fixtures.json', import.meta.url
)));
const hydrate = value => {
  if (value && typeof value === 'object') {
    if ('$number' in value) return Number(value.$number);
    return Object.fromEntries(Object.entries(value).map(([key, item]) => [key, hydrate(item)]));
  }
  return value;
};
const slotFor = (scalar, field = 'width') => ({
  id: 'features', renderer: 'features', enabled: true, side: 'inside',
  width: null, radius: null, inner_gap_px: null, outer_gap_px: null,
  z: 0, params: {}, [field]: scalar
});
const validateRow = slot => validateCustomTrackPlan({
  mode: 'circular', slots: [slot], axisIndex: 0, trackType: 'tuckin',
  depthTrackCount: 0, annotationSetIds: [], conservationSeries: []
});
const configFor = slot => ({
  form: createDefaultForm(),
  adv: { ...createDefaultAdv(), circular_track_slots_enabled: true, circular_track_slots: [slot] }
});

for (const fixture of fixtures) {
  test(`S00 scientific fixture: ${fixture.name}`, () => {
    assert.equal(globalThis.window, undefined);
    assert.equal(globalThis.document, undefined);
    const scalar = hydrate(fixture.input);
    const before = structuredClone(scalar);
    const view = readCircularMeasure(scalar);
    const projected = projectCircularMeasureDraft(scalar);
    assert.deepEqual(scalar, before, 'reading/projection preserves raw number, text and unit');
    for (const field of ['width', 'radius']) {
      const slot = slotFor(scalar, field);
      const config = configFor(slot);
      if (fixture.valid) {
        assert.equal(view.error, null);
        assert.equal(view.isAuto, fixture.canonical === null);
        assert.deepEqual(parseOptionalCircularScalar(scalar), fixture.canonical);
        assert.equal(validateRow(slot).rowIssues.size, 0);
        assert.doesNotThrow(() => validateCurrentWriterActiveConfig({ mode: 'circular', storedConfig: config }));
        assert.deepEqual(buildCircularTrackSlotPayload(slot)[field], fixture.canonical);
        const draft = writeCircularMeasureValue(view.valueText, view.selectedUnit ?? 'factor');
        assert.deepEqual(parseOptionalCircularScalar(draft), fixture.canonical, 'view → edit → exact canonical meaning');
        assert.deepEqual(parseOptionalCircularScalar(projected), fixture.canonical, 'request projection preserves meaning');
        if (typeof scalar?.value === 'string') assert.equal(draft.value, scalar.value, 'numeric lexeme is not rounded');
      } else {
        assert.equal(view.isAuto, false);
        assert.match(view.error, /positive finite px or factor/);
        assert.throws(() => parseOptionalCircularScalar(scalar));
        assert.throws(() => buildCircularTrackSlotPayload(slot));
        assert.ok(validateRow(slot).rowIssues.get(0).some(issue => issue.code === 'geometry_invalid'));
        assert.throws(() => validateCurrentWriterActiveConfig({ mode: 'circular', storedConfig: config }));
        if (scalar && typeof scalar === 'object') assert.strictEqual(projected, scalar, 'invalid typed projection does not repair');
        if (typeof scalar === 'string') assert.equal(view.valueText, scalar, 'invalid text stays visible');
      }
    }
    assert.deepEqual(scalar, before);
  });
}

test('selector owns plain numeric text; only complete valid suffix input overrides it', () => {
  const cases = [
    ['1.5', 'px', '1.5', 'px', 1.5],
    ['1.5', 'factor', '1.5', 'factor', 1.5],
    ['  1.  ', 'px', '1.', 'px', 1],
    ['1e-3', 'factor', '1e-3', 'factor', 0.001],
    ['20px', 'factor', '20', 'px', 20],
    ['1. PX', 'factor', '1.', 'px', 1],
    ['65%', 'px', '0.65', 'factor', 0.65],
    ['0.12345678901234567%', 'px', '0.0012345678901234567', 'factor', 0.0012345678901234567]
  ];
  for (const [text, selected, lexeme, unit, value] of cases) {
    const draft = writeCircularMeasureValue(text, selected);
    assert.deepEqual(draft, { value: lexeme, unit });
    assert.deepEqual(parseOptionalCircularScalar(draft), { value, unit });
  }
  for (const text of ['bad', '1e', '-', '.', '20p', '20pxx', '1epx', '0px', '-1%', 'Infinity', 'Infinitypx', 'NaN', '1e309', '5e-324%']) {
    const draft = writeCircularMeasureValue(text, 'px');
    assert.deepEqual(draft, { value: text, unit: 'px' });
    assert.equal(readCircularMeasure(draft).isAuto, false);
    assert.throws(() => parseOptionalCircularScalar(draft));
  }
  assert.throws(() => writeCircularMeasureValue('1', '%'), /unit/);
});

test('unit edit preserves valid and invalid lexemes, is pure, and returns same scalar for a no-op', () => {
  for (const text of ['1.', '1e-3', '0.12345678901234567', 'bad', '1e', '0', '-1', 'Infinity', '']) {
    const scalar = Object.freeze({ value: text, unit: 'factor' });
    assert.strictEqual(changeCircularMeasureUnit(scalar, 'factor'), scalar);
    const changed = changeCircularMeasureUnit(scalar, 'px');
    assert.deepEqual(changed, { value: text, unit: 'px' });
    assert.deepEqual(changeCircularMeasureUnit(changed, 'factor'), scalar);
    assert.deepEqual(scalar, { value: text, unit: 'factor' });
  }
  assert.deepEqual(changeCircularMeasureUnit('1. px', 'factor'), { value: '1.', unit: 'factor' });
  assert.deepEqual(changeCircularMeasureUnit('65%', 'px'), { value: '0.65', unit: 'px' });
  assert.deepEqual(changeCircularMeasureUnit('1e', 'px'), { value: '1e', unit: 'px' });
  assert.deepEqual(changeCircularMeasureUnit({ value: '1e', unit: 'em' }, 'px'), { value: '1e', unit: 'px' });
  assert.equal(changeCircularMeasureUnit(Infinity, 'px').value, Infinity);
  assert.throws(() => parseOptionalCircularScalar(changeCircularMeasureUnit(true, 'px')));
});

test('Auto has no canonical unit; explicit clear returns null while typed empty stays invalid on read/unit edit', () => {
  for (const scalar of [null, undefined, '']) {
    assert.deepEqual(readCircularMeasure(scalar), { valueText: '', selectedUnit: null, isAuto: true, error: null });
    assert.equal(changeCircularMeasureUnit(scalar, 'px'), null);
  }
  for (const unit of ['px', 'factor']) {
    assert.equal(writeCircularMeasureValue('', unit), null);
    assert.equal(writeCircularMeasureValue('  ', unit), null);
    const empty = { value: '', unit };
    assert.equal(readCircularMeasure(empty).isAuto, false);
    assert.match(readCircularMeasure(empty).error, /positive/);
    assert.strictEqual(projectCircularMeasureDraft(empty), empty);
  }
});

const tobacco = JSON.parse(await readFile(new URL(
  '../../gbdraw/web/gallery/sessions/tobacco-chloroplast.gbdraw-session.json', import.meta.url
)));

test('tobacco Gallery pairs and precision survive read and real canonical request projection without mutation', () => {
  const expected = [
    ['plastome_regions', 'width', '20', 'px'], ['plastome_regions', 'radius', '0.65', 'factor'],
    ['gc_content', 'width', '0.08', 'factor'], ['gc_content', 'radius', '0.56', 'factor']
  ];
  const before = structuredClone(tobacco);
  for (const [id, field, valueText, selectedUnit] of expected) {
    const scalar = tobacco.config.adv.circular_track_slots.find(slot => slot.id === id)[field];
    assert.deepEqual(readCircularMeasure(scalar), { valueText, selectedUnit, isAuto: false, error: null });
    const draft = writeCircularMeasureValue(valueText, selectedUnit);
    assert.deepEqual(parseOptionalCircularScalar(draft), scalar);
  }
  const projected = projectCanonicalSessionRequest(tobacco);
  for (const [id, field] of expected) {
    const raw = tobacco.config.adv.circular_track_slots.find(slot => slot.id === id)[field];
    const draft = projected.config.adv.circular_track_slots.find(slot => slot.id === id)[field];
    assert.deepEqual(parseOptionalCircularScalar(draft), raw);
  }
  assert.deepEqual(tobacco, before, 'request/config/Result read is pure');
});

test('current and historical request projection keeps zero pixel gaps separate from positive Circular scalars', () => {
  const current = structuredClone(tobacco);
  const slots = current.renderRequest.diagramOptions.tracks.circularTrackSlots;
  const slot = slots.find(slot => slot.id === 'gc_content');
  slot.innerGapPx = 0;
  slot.outerGapPx = 10;
  assert.equal(projectCanonicalSessionRequest(current).config.adv.circular_track_slots.find(slot => slot.id === 'gc_content').inner_gap_px, '0');
  for (const spacing of [0, '0px', { value: 0, unit: 'px' }, { value: 5, unit: 'px' }]) {
    const historical = structuredClone(tobacco);
    historical.renderRequest.schema = 2;
    delete historical.renderRequest.diagramOptions.featureOverrides;
    historical.renderRequest.diagramOptions.output.outputPrefix = 'legacy';
    delete historical.renderRequest.diagramOptions.featurePlacements;
    historical.renderRequest.records.forEach(record => delete record.display);
    historical.renderRequest.diagramOptions.tracks.circularTrackSlots = [{
      id: 'legacy', renderer: 'dinucleotide_content', enabled: true, side: 'inside',
      width: '20px', radius: '65%', spacing, z: 0, params: { nt: 'GC' }
    }];
    const before = structuredClone(historical);
    const restored = projectCanonicalSessionRequest(historical).config.adv.circular_track_slots[0];
    const pixels = typeof spacing === 'object' ? spacing.value : spacing === '0px' ? 0 : spacing;
    assert.equal(restored.inner_gap_px, String(pixels));
    assert.equal(restored.outer_gap_px, String(pixels));
    assert.deepEqual(buildCircularTrackSlotPayload(restored).width, { value: 20, unit: 'px' });
    assert.deepEqual(buildCircularTrackSlotPayload(restored).radius, { value: 0.65, unit: 'factor' });
    assert.deepEqual(historical, before);
    // A factor or percent is not a physical pixel gap; do not invent an R conversion.
    for (const invalid of [{ value: 0.1, unit: 'factor' }, { value: 5, unit: '%' }, '5%']) {
      historical.renderRequest.diagramOptions.tracks.circularTrackSlots[0].spacing = invalid;
      assert.throws(() => projectCanonicalSessionRequest(historical), /pixels/);
    }
  }
});
