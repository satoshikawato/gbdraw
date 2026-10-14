// @ts-check
// The generic repair of a value against a contract description
// (utils/contract-repair.js): the six rules of the stored-value repair, in
// order, the identity and no-mutation guarantees, the root row, and each hook
// of the policy.
import assert from 'node:assert/strict';
import test from 'node:test';

import { repairByDescription } from '../../gbdraw/web/js/utils/contract-repair.js';
import { shapeViolations } from '../../gbdraw/web/js/utils/contract-shape.js';
import { CANONICAL_REQUEST_DESCRIPTION } from '../../gbdraw/web/js/contracts/request-shape.generated.js';

/** @import { ContractDescription, ContractShape } from '../../gbdraw/web/js/utils/contract-shape.js' */
/** @import { RepairPolicy } from '../../gbdraw/web/js/utils/contract-repair.js' */

/** @type {ContractShape} */
const SOURCE = {
  kind: 'object',
  fields: { kind: { kind: 'literal', values: ['genbank'] }, resourceId: { kind: 'string' } },
  optional: [],
  open: false
};
/** @type {ContractShape} */
const PAIR = {
  kind: 'object',
  fields: {
    kind: { kind: 'literal', values: ['pair'] },
    left: { kind: 'integer' },
    right: { kind: 'integer' },
    note: { kind: 'string' }
  },
  optional: ['note'],
  open: false
};
/** @type {ContractShape} */
const TABLE = {
  kind: 'object',
  fields: {
    kind: { kind: 'literal', values: ['table'] },
    resourceId: { kind: 'string' },
    settings: { kind: 'ref', name: 'Settings' },
    selectors: { kind: 'array', items: { kind: 'string', pattern: '\\S' }, minItems: 1 }
  },
  optional: ['settings', 'selectors'],
  open: false
};

/** @type {ContractDescription} */
const DESCRIPTION = {
  root: 'Request',
  types: {
    Request: {
      kind: 'object',
      fields: {
        schema: { kind: 'literal', values: [9] },
        mode: { kind: 'literal', values: ['linear', 'circular'] },
        records: { kind: 'array', items: { kind: 'ref', name: 'Record' } },
        comparisons: { kind: 'array', items: { kind: 'ref', name: 'Comparison' } },
        formats: { kind: 'array', items: { kind: 'literal', values: ['svg', 'png'] } },
        options: { kind: 'ref', name: 'Options' },
        layout: { kind: 'ref', name: 'Layout' },
        window: { kind: 'tuple', items: [{ kind: 'integer' }, { kind: 'integer' }] },
        tags: { kind: 'record', values: { kind: 'string' } }
      },
      optional: ['comparisons', 'formats', 'layout', 'window', 'tags'],
      open: false,
      defaults: { formats: ['svg'], options: { dinucleotide: 'GC' } }
    },
    Record: {
      kind: 'object',
      fields: {
        source: { kind: 'ref', name: 'Source' },
        display: { kind: 'ref', name: 'Display' }
      },
      optional: ['display'],
      open: false
    },
    Source: SOURCE,
    Display: {
      kind: 'object',
      fields: { label: { kind: 'string' }, strand: { kind: 'literal', values: ['+', '-'] } },
      optional: [],
      open: false
    },
    Comparison: { kind: 'union', options: [{ kind: 'ref', name: 'Pair' }, { kind: 'ref', name: 'Table' }] },
    Pair: PAIR,
    Table: TABLE,
    Settings: {
      kind: 'object',
      fields: {
        identity: { kind: 'number', minimum: 0, maximum: 1 },
        method: { kind: 'literal', values: ['rbh', 'none'], spellings: { aliases: { reciprocal_best: 'rbh' }, foldHyphens: true, blank: 'rbh' } }
      },
      optional: ['method'],
      open: false
    },
    Options: {
      kind: 'object',
      fields: {
        dinucleotide: { kind: 'string' },
        width: { kind: 'union', options: [{ kind: 'integer' }, { kind: 'null' }] },
        ringFiles: { kind: 'array', items: { kind: 'string' } },
        ringLabels: { kind: 'array', items: { kind: 'string' } },
        ringColors: { kind: 'array', items: { kind: 'string' } }
      },
      optional: ['dinucleotide', 'width', 'ringFiles', 'ringLabels', 'ringColors'],
      open: false,
      defaults: { dinucleotide: 'GC' }
    },
    Layout: {
      kind: 'object',
      fields: {
        alignment: { kind: 'union', options: [{ kind: 'ref', name: 'Alignment' }, { kind: 'null' }] },
        gap: { kind: 'number' }
      },
      optional: ['alignment'],
      open: false
    },
    Alignment: {
      kind: 'object',
      fields: {
        reference: { kind: 'object', fields: { recordKey: { kind: 'string' } }, optional: [], open: false },
        decisions: { kind: 'array', items: { kind: 'string' } }
      },
      optional: [],
      open: false
    }
  }
};

const ROOT = Object.freeze({ type: 'Request', path: 'renderRequest' });

/** @type {RepairPolicy<'unknown-resource'>} */
const POLICY = {
  positionalArrays: ['records'],
  partBoundary: (pattern) => ['layout.alignment', 'records[*].display', 'options', 'comparisons[*]'].includes(pattern)
};

/** @returns {Record<string, unknown>} */
const validRequest = () => ({
  schema: 9,
  mode: 'linear',
  records: [
    { source: { kind: 'genbank', resourceId: 'r1' }, display: { label: 'A', strand: '+' } },
    { source: { kind: 'genbank', resourceId: 'r2' } }
  ],
  comparisons: [
    { kind: 'pair', left: 0, right: 1 },
    { kind: 'table', resourceId: 't1', settings: { identity: 0.5, method: 'rbh' } }
  ],
  formats: ['svg', 'png'],
  options: { dinucleotide: 'AT', width: null, ringFiles: ['a', 'b'], ringLabels: ['A', 'B'], ringColors: ['red', 'blue'] },
  layout: { alignment: { reference: { recordKey: 'k1' }, decisions: ['keep'] }, gap: 2 },
  window: [1, 2],
  tags: { a: 'x' }
});

/**
 * Deeply frozen, so a mutation throws in strict mode.
 * @template T
 * @param {T} value
 * @returns {T}
 */
const deepFreeze = (value) => {
  if (value && typeof value === 'object') {
    Object.values(value).forEach(deepFreeze);
    Object.freeze(value);
  }
  return value;
};

/**
 * Repairs a frozen copy of `request` and checks that the result conforms.
 * @param {Record<string, unknown>} request
 * @param {RepairPolicy<'unknown-resource'>} [policy]
 */
const repair = (request, policy = POLICY) => {
  const input = deepFreeze(structuredClone(request));
  const before = JSON.stringify(input);
  const result = repairByDescription(input, ROOT, DESCRIPTION, policy);
  assert.equal(JSON.stringify(input), before, 'the input is never changed');
  if (result.value !== undefined) {
    assert.deepEqual(shapeViolations(result.value, { kind: 'ref', name: 'Request' }, DESCRIPTION.types, 'renderRequest'), []);
  }
  return result;
};

/** @param {unknown} value @returns {Record<string, any>} */
const obj = (value) => /** @type {Record<string, any>} */ (value);

test('a conforming value comes back as the same object, without rows', () => {
  const input = deepFreeze(validRequest());
  const result = repairByDescription(input, ROOT, DESCRIPTION, POLICY);
  assert.equal(result.value, input);
  assert.deepEqual(result.rows, []);
});

test('rule 1: a lossless reading replaces the value without a row, when it conforms', () => {
  const request = validRequest();
  obj(request.comparisons)[1].settings.method = 'R B H';
  /** @type {string[]} */
  const asked = [];
  const result = repair(request, {
    ...POLICY,
    lossless: (pattern, value) => {
      asked.push(pattern);
      return typeof value === 'string' && pattern === 'comparisons[*].settings.method'
        ? { value: value.replaceAll(' ', '').toLowerCase() } : null;
    }
  });
  assert.deepEqual(result.rows, []);
  assert.equal(obj(result.value).comparisons[1].settings.method, 'rbh');
  assert.ok(asked.includes('comparisons[*].settings.method'));
  // A reading that does not conform is not used; the later rules apply.
  const rejected = repair(request, { ...POLICY, lossless: () => ({ value: 'RBH!' }) });
  assert.deepEqual(rejected.rows, [{
    path: 'renderRequest.comparisons[1].settings.method', pattern: 'comparisons[*].settings.method',
    reason: 'not-in-set', rule: 4, effect: 'dropped'
  }]);
});

test('rule 2: an unknown key is dropped', () => {
  const request = { ...validRequest(), extra: 1 };
  const result = repair(request);
  assert.deepEqual(result.rows, [{
    path: 'renderRequest.extra', pattern: 'extra', reason: 'unknown-key', rule: 2, effect: 'dropped'
  }]);
  assert.equal(Object.hasOwn(obj(result.value), 'extra'), false);
  assert.deepEqual(Object.keys(obj(result.value)), Object.keys(validRequest()));
});

test('rule 3: a field with a default takes it, present or missing', () => {
  const request = validRequest();
  obj(request.options).dinucleotide = 7;
  const result = repair(request);
  assert.deepEqual(result.rows, [{
    path: 'renderRequest.options.dinucleotide', pattern: 'options.dinucleotide',
    reason: 'wrong-kind', rule: 3, effect: 'defaulted', used: 'GC'
  }]);
  assert.equal(obj(result.value).options.dinucleotide, 'GC');
  const defaulted = repair({ ...validRequest(), formats: 'svg' });
  assert.deepEqual(defaulted.rows.map(({ rule, effect, used }) => ({ rule, effect, used })), [
    { rule: 3, effect: 'defaulted', used: ['svg'] }
  ]);
  // The default is a copy: changing the result never changes the description.
  obj(defaulted.value).formats.push('png');
  assert.deepEqual(DESCRIPTION.types.Request.kind === 'object' && DESCRIPTION.types.Request.defaults?.formats, ['svg']);
});

test('rule 4: an optional field without a default is dropped', () => {
  const request = validRequest();
  obj(request.options).width = 'wide';
  const result = repair(request);
  assert.deepEqual(result.rows, [{
    path: 'renderRequest.options.width', pattern: 'options.width', reason: 'wrong-kind', rule: 4, effect: 'dropped'
  }]);
  assert.equal(Object.hasOwn(obj(result.value).options, 'width'), false);
  // Rule 4 is for the place of the violation only: an item of an optional
  // tuple is not optional, and the tuple is not a part.
  const tuple = validRequest();
  tuple.window = [1, 'x'];
  assert.deepEqual(repair(tuple).rows, [{
    path: 'renderRequest', pattern: '', reason: 'wrong-kind', rule: 6, effect: 'dropped'
  }]);
  assert.deepEqual(repair(tuple, { ...POLICY, partBoundary: (pattern) => pattern === 'window' }).rows, [{
    path: 'renderRequest.window', pattern: 'window', reason: 'wrong-kind', rule: 6, effect: 'dropped'
  }]);
});

test('rule 5: an element nothing names by position is dropped, the innermost one', () => {
  const request = validRequest();
  obj(request.comparisons)[0].left = 'x';
  obj(request.formats).splice(1, 0, 'gif');
  const result = repair(request);
  assert.deepEqual(result.rows, [
    { path: 'renderRequest.comparisons[0]', pattern: 'comparisons[*]', reason: 'wrong-kind', rule: 5, effect: 'dropped' },
    { path: 'renderRequest.formats[1]', pattern: 'formats[*]', reason: 'not-in-set', rule: 5, effect: 'dropped' }
  ]);
  assert.equal(obj(result.value).comparisons.length, 1);
  assert.equal(obj(result.value).comparisons[0].kind, 'table');
  assert.deepEqual(obj(result.value).formats, ['svg', 'png']);
  // A required field missing inside an element drops the element.
  const missing = validRequest();
  delete obj(missing.comparisons)[1].resourceId;
  assert.deepEqual(repair(missing).rows, [{
    path: 'renderRequest.comparisons[1]', pattern: 'comparisons[*]', reason: 'missing-key', rule: 5, effect: 'dropped'
  }]);
  // An entry of a record is optional: it is dropped by rule 4.
  assert.deepEqual(repair({ ...validRequest(), tags: { a: 'x', b: 2 } }).rows, [{
    path: 'renderRequest.tags.b', pattern: 'tags.*', reason: 'wrong-kind', rule: 4, effect: 'dropped'
  }]);
});

test('rule 6: the enclosing part is dropped, outward until one can be', () => {
  const request = validRequest();
  obj(request.records)[0].display.strand = '?';
  obj(request.layout).alignment.reference.recordKey = 3;
  const result = repair(request);
  assert.deepEqual(result.rows, [
    { path: 'renderRequest.records[0].display', pattern: 'records[*].display', reason: 'not-in-set', rule: 6, effect: 'dropped' },
    { path: 'renderRequest.layout.alignment', pattern: 'layout.alignment', reason: 'wrong-kind', rule: 6, effect: 'dropped' }
  ]);
  assert.equal(Object.hasOwn(obj(result.value).records[0], 'display'), false);
  assert.equal(Object.hasOwn(obj(result.value).layout, 'alignment'), false);
  // A field that is not a part is not dropped on its own: `layout` is not a
  // part and `gap` is required, so the walk reaches the root.
  const gapless = validRequest();
  obj(gapless.layout).gap = 'x';
  assert.deepEqual(repair(gapless).rows, [{
    path: 'renderRequest', pattern: '', reason: 'wrong-kind', rule: 6, effect: 'dropped'
  }]);
});

test('rule 6 at the root returns no value and one root row', () => {
  for (const input of [null, 'request', [], { ...validRequest(), records: 'x' }, { ...validRequest(), mode: 'radial' }]) {
    const result = repairByDescription(deepFreeze(input), ROOT, DESCRIPTION, POLICY);
    assert.equal(result.value, undefined);
    assert.equal(result.rows.length, 1);
    assert.deepEqual({ ...result.rows[0], reason: 'any' }, {
      path: 'renderRequest', pattern: '', reason: 'any', rule: 6, effect: 'dropped'
    });
  }
  // A required element of a positional array cannot go: the walk reaches the root.
  const request = validRequest();
  obj(request.records)[1].source.kind = 'gff';
  assert.deepEqual(repair(request).rows, [{
    path: 'renderRequest', pattern: '', reason: 'not-in-set', rule: 6, effect: 'dropped'
  }]);
});

test('rows inside a dropped part are not reported; the part row is', () => {
  const request = validRequest();
  obj(request.comparisons)[1].settings = { identity: 'x', extra: 1 };
  obj(request.comparisons)[1].extra = 1;
  const result = repair(request);
  assert.deepEqual(result.rows, [{
    path: 'renderRequest.comparisons[1]', pattern: 'comparisons[*]', reason: 'wrong-kind', rule: 5, effect: 'dropped'
  }]);
});

test('reference hook: a conforming value it rejects is repaired like a violation', () => {
  const request = validRequest();
  obj(request.comparisons)[0].right = 5;
  obj(request.layout).alignment.reference.recordKey = 'gone';
  /** @type {RepairPolicy<'unknown-resource'>} */
  const policy = {
    ...POLICY,
    reference: (pattern, value) => {
      if (pattern === 'comparisons[*].right' && value === 5) return 'unknown-resource';
      if (pattern === 'layout.alignment.reference.recordKey' && value === 'gone') return 'unknown-resource';
      return null;
    }
  };
  const result = repair(request, policy);
  assert.deepEqual(result.rows, [
    { path: 'renderRequest.comparisons[0]', pattern: 'comparisons[*]', reason: 'unknown-resource', rule: 5, effect: 'dropped' },
    { path: 'renderRequest.layout.alignment', pattern: 'layout.alignment', reason: 'unknown-resource', rule: 6, effect: 'dropped' }
  ]);
  // Without a rejected value the request is returned as it is.
  const clean = deepFreeze(validRequest());
  assert.equal(repairByDescription(clean, ROOT, DESCRIPTION, policy).value, clean);
});

test('parallel arrays: rule 5 drops the same column of every array of the group', () => {
  const request = validRequest();
  obj(request.options).ringLabels = ['A', 7];
  const result = repair(request, {
    ...POLICY,
    parallelArrays: [['options.ringFiles', 'options.ringLabels', 'options.ringColors']]
  });
  assert.deepEqual(result.rows, [{
    path: 'renderRequest.options.ringLabels[1]', pattern: 'options.ringLabels[*]', reason: 'wrong-kind', rule: 5, effect: 'dropped'
  }]);
  assert.deepEqual(obj(result.value).options.ringFiles, ['a']);
  assert.deepEqual(obj(result.value).options.ringLabels, ['A']);
  assert.deepEqual(obj(result.value).options.ringColors, ['red']);
  // Without the group only that array loses its element.
  assert.deepEqual(obj(repair(request).value).options.ringFiles, ['a', 'b']);
});

test('discriminator: the option its field names is repaired, never another one', () => {
  const policy = { ...POLICY, discriminator: (/** @type {string} */ pattern) => (pattern === 'comparisons[*]' ? 'kind' : null) };
  // A table with two stray keys has as many violations as a pair with a wrong
  // kind and a stray key; the first of the fewest (a pair) cannot be read.
  const request = validRequest();
  obj(request.comparisons)[0] = { kind: 'table', resourceId: 't2', left: 0, right: 1 };
  assert.deepEqual(repair(request).rows, [{
    path: 'renderRequest.comparisons[0]', pattern: 'comparisons[*]', reason: 'not-in-set', rule: 5, effect: 'dropped'
  }]);
  const read = repair(request, policy);
  assert.deepEqual(read.rows.map(({ path, rule }) => ({ path, rule })), [
    { path: 'renderRequest.comparisons[0].left', rule: 2 },
    { path: 'renderRequest.comparisons[0].right', rule: 2 }
  ]);
  assert.deepEqual(obj(read.value).comparisons[0], { kind: 'table', resourceId: 't2' });
  // A kind that names no option drops the element.
  const unnamed = validRequest();
  obj(unnamed.comparisons)[0].kind = 'triple';
  assert.deepEqual(repair(unnamed, policy).rows, [{
    path: 'renderRequest.comparisons[0]', pattern: 'comparisons[*]', reason: 'not-in-set', rule: 5, effect: 'dropped'
  }]);
});

test('rule 1: a spelling the description names reads as its name, without a row', () => {
  const request = validRequest();
  obj(request.comparisons)[1].settings.method = ' Reciprocal-Best ';
  const result = repair(request);
  assert.deepEqual(result.rows, []);
  assert.equal(obj(result.value).comparisons[1].settings.method, 'rbh');
  const blank = validRequest();
  obj(blank.comparisons)[1].settings.method = '';
  assert.equal(obj(repair(blank).value).comparisons[1].settings.method, 'rbh');
});

test('a value outside its bounds is repaired like any violation', () => {
  const request = validRequest();
  obj(request.comparisons)[1].settings.identity = 1.5;
  assert.deepEqual(repair(request).rows, [{
    path: 'renderRequest.comparisons[1]', pattern: 'comparisons[*]', reason: 'out-of-range', rule: 5, effect: 'dropped'
  }]);
});

test('too few items: an array that starts or ends too short goes, by the rule for its place', () => {
  const request = validRequest();
  obj(request.comparisons)[1].selectors = [];
  assert.deepEqual(repair(request).rows, [{
    path: 'renderRequest.comparisons[1].selectors', pattern: 'comparisons[*].selectors', reason: 'too-short', rule: 4, effect: 'dropped'
  }]);
  // Its only item is dropped by rule 5; the array is then too short itself,
  // and it is optional (only the array's row is reported).
  const emptied = validRequest();
  obj(emptied.comparisons)[1].selectors = [' '];
  assert.deepEqual(repair(emptied).rows, [{
    path: 'renderRequest.comparisons[1].selectors', pattern: 'comparisons[*].selectors', reason: 'too-short', rule: 4, effect: 'dropped'
  }]);
  const kept = validRequest();
  obj(kept.comparisons)[1].selectors = ['a', ' '];
  const result = repair(kept);
  assert.deepEqual(result.rows.map(({ path, rule }) => ({ path, rule })), [
    { path: 'renderRequest.comparisons[1].selectors[1]', rule: 5 }
  ]);
  assert.deepEqual(obj(result.value).comparisons[1].selectors, ['a']);
});

test('a required field with a default takes it, at its own violation and as a part', () => {
  const wrong = { ...validRequest(), options: 'none' };
  assert.deepEqual(repair(wrong).rows, [{
    path: 'renderRequest.options', pattern: 'options', reason: 'wrong-kind', rule: 3, effect: 'defaulted', used: { dinucleotide: 'GC' }
  }]);
  const { options: _options, ...missing } = validRequest();
  const read = repair(missing);
  assert.deepEqual(read.rows.map(({ path, reason, rule }) => ({ path, reason, rule })), [
    { path: 'renderRequest.options', reason: 'missing-key', rule: 3 }
  ]);
  assert.deepEqual(obj(read.value).options, { dinucleotide: 'GC' });
  // `options` is a part: a violation inside it that its own fields cannot
  // absorb puts the default of the whole part.
  const inside = validRequest();
  obj(inside.options).ringFiles = 'x';
  obj(inside.options).extra = 1;
  const part = repair(inside, { ...POLICY, partBoundary: (pattern) => pattern === 'options' });
  assert.deepEqual(part.rows.map(({ path, rule }) => ({ path, rule })), [
    { path: 'renderRequest.options.ringFiles', rule: 4 },
    { path: 'renderRequest.options.extra', rule: 2 }
  ]);
  const deep = validRequest();
  obj(deep.options).ringFiles = [1];
  const ringPolicy = { ...POLICY, positionalArrays: ['records', 'options.ringFiles'], partBoundary: (/** @type {string} */ pattern) => pattern === 'options' };
  assert.deepEqual(repair(deep, ringPolicy).rows, [{
    path: 'renderRequest.options', pattern: 'options', reason: 'wrong-kind', rule: 6, effect: 'defaulted', used: { dinucleotide: 'GC' }
  }]);
});

test('an annotation colour: a spelling of a name reads as it, another text takes the default', () => {
  // annotations.models._color, described as data: a hex of 3 or 6 digits or a
  // colour name in any case, with outer spaces; a fill is unset as null or blank.
  /** @param {Record<string, unknown>} changes */
  const style = (changes) => ({
    stroke: '#404040', strokeWidth: 1.5, strokeDasharray: [], lineCap: 'tick', fill: null, fillOpacity: 0.2,
    hatch: { angle: 45, spacing: 5, color: '#666666', width: 1, cross: false }, labelColor: '#202020',
    labelFontSize: null, labelOrientation: 'auto', labelPosition: 'center', labelOffset: 4, ...changes
  });
  const root = { type: 'CanonicalRegionAnnotationStyle', path: 'style' };
  /** @type {RepairPolicy} */
  const policy = { positionalArrays: [], partBoundary: () => false };
  /** @param {Record<string, unknown>} changes */
  const read = (changes) => repairByDescription(style(changes), root, CANONICAL_REQUEST_DESCRIPTION, policy);
  const kept = style({ stroke: ' #AbC ', fill: '  ', labelColor: 'navy' });
  assert.equal(repairByDescription(kept, root, CANONICAL_REQUEST_DESCRIPTION, policy).value, kept);
  const spelled = read({ stroke: ' Red ', labelColor: 'NAVY' });
  assert.deepEqual(spelled.rows, []);
  assert.deepEqual([obj(spelled.value).stroke, obj(spelled.value).labelColor], ['red', 'navy']);
  const unusable = read({ stroke: 'notacolor', fill: '#abcd', labelColor: 'rgb(1,2,3)', hatch: { ...style({}).hatch, color: 'none' } });
  assert.deepEqual(unusable.rows.map(({ path, reason, rule, used }) => ({ path, reason, rule, used })), [
    { path: 'style.stroke', reason: 'pattern', rule: 3, used: '#404040' },
    { path: 'style.fill', reason: 'pattern', rule: 3, used: null },
    { path: 'style.hatch.color', reason: 'pattern', rule: 3, used: '#666666' },
    { path: 'style.labelColor', reason: 'pattern', rule: 3, used: '#202020' }
  ]);
});

test('a union: an object is repaired as the object option, even with more violations than null has', () => {
  /** @type {ContractDescription} */
  const description = {
    root: 'Holder',
    types: {
      Holder: { kind: 'object', fields: { style: { kind: 'union', options: [{ kind: 'ref', name: 'Style' }, { kind: 'null' }] } }, optional: [], open: false },
      Style: {
        kind: 'object', fields: { a: { kind: 'integer' }, b: { kind: 'integer' } }, optional: [], open: false, defaults: { a: 1, b: 2 }
      }
    }
  };
  const result = repairByDescription({ style: { a: 'x', b: 'y' } }, { type: 'Holder', path: 'h' }, description,
    { positionalArrays: [], partBoundary: () => false });
  assert.deepEqual(result.value, { style: { a: 1, b: 2 } });
  assert.deepEqual(result.rows.map(({ path, rule }) => ({ path, rule })), [{ path: 'h.style.a', rule: 3 }, { path: 'h.style.b', rule: 3 }]);
});
