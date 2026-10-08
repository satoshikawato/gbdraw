import assert from 'node:assert/strict';
import { execFileSync } from 'node:child_process';
import { readFileSync } from 'node:fs';
import test from 'node:test';
import { normalizeCircularGeometryShortcuts } from '../../gbdraw/web/js/services/circular-track-slot-model.js';
import { diagnosticError, normalizeCaughtError, normalizeUserFacingError } from '../../gbdraw/web/js/utils/error-normalization.js';
import { deserializeWorkerError, normalizeGenerationResponse } from '../../gbdraw/web/js/services/diagram-generation.js';
globalThis.self = {};
const { serializeError, callJsonHelper, resolveGenerationCleanupOutcome } = await import('../../gbdraw/web/js/workers/diagram-generation-worker.js');

const sentinel = 'PRIVATE_SENTINEL_'.repeat(5000);
const roundtrip = (source) => normalizeUserFacingError(deserializeWorkerError(structuredClone(serializeError(source))));
const source = { code: 'REGEX_SYNTAX', operation: 'evaluateRules', stage: 'rule-validation',
  context: { position: 1, positionUnit: 'python-character', reason: 'UNTERMINATED_SET',
    row: 2, field: 'pattern', pattern: sentinel, sequence: sentinel, path: sentinel },
  cause: new Error(sentinel), message: sentinel, name: sentinel, stdout: sentinel, stderr: sentinel,
  traceback: sentinel, stack: sentinel, notes: [sentinel],
  details: [{ label: sentinel, text: sentinel }],
  secondary: [{ code: 'CLEANUP_FAILED', stage: 'cleanup', message: sentinel }] };
const normalized = normalizeUserFacingError(source);
assert.deepEqual(roundtrip(source), normalized);
assert.deepEqual(normalizeUserFacingError(normalized), normalized);
assert.match(normalized.summary, /Python character position 1/);
assert.equal(normalized.context.row, 2);
assert.deepEqual(normalized.actions, ['edit-pattern', 'retry']);
assert.doesNotMatch(JSON.stringify(normalized), /PRIVATE_SENTINEL/);
assert.doesNotMatch(JSON.stringify(serializeError(source)), /PRIVATE_SENTINEL/);
assert.doesNotMatch(JSON.stringify(deserializeWorkerError(source)), /PRIVATE_SENTINEL/);
assert(normalized.summary.length <= 1000);
assert(normalized.details.length <= 8);
assert(normalized.details.every(({ text }) => text.length <= 4000));
for (const value of [sentinel, new Error(sentinel), { name: 'PythonError', message: `Traceback\nValueError: ${sentinel}` },
  { ...source, code: sentinel, context: { field: sentinel, row: -1, position: 1.5, reason: sentinel } }]) {
  const result = roundtrip(value);
  assert.equal(result.code, 'UNKNOWN');
  assert.doesNotMatch(JSON.stringify(result), /PRIVATE_SENTINEL|Traceback|ValueError/);
}
const unknown = normalizeUserFacingError({ code: 'UNKNOWN' });
assert.equal(unknown.stage, 'unknown');
assert.equal(unknown.context.position, undefined);
assert.equal(normalizeUserFacingError({ code: 'REGEX_SYNTAX', context: { position: 2 } }).context.position, undefined);
assert.deepEqual(normalizeUserFacingError({ ...source, context: { field: 'pattern', row: '2', position: Infinity,
  positionUnit: 'utf16', reason: 'bogus' } }).context, { field: 'pattern' });
assert(normalizeUserFacingError(source, { summaryLimit: 8 }).summary.length <= 8);
assert(normalizeUserFacingError(source, { detailLimit: 8 }).details[0].text.length <= 8);
assert(normalizeUserFacingError(source, { summaryLimit: 5000, detailLimit: 9000 }).summary.length <= 1000);

// Diagnostic numeric boundaries must survive transport exactly, while larger,
// fractional, nonfinite and nonnumeric values are omitted rather than clamped.
for (const key of ['position', 'row', 'column', 'inputOrdinal', 'recordIndex', 'seriesIndex', 'slotIndex', 'recordCount', 'columnCount', 'codepoint', 'innerPx', 'outerPx']) {
  const maximum = key === 'codepoint' ? 0x10ffff : 10000000;
  const bounded = value => roundtrip({ ...source, context: { [key]: value, positionUnit: 'python-character' } });
  for (const value of [0, maximum]) assert.equal(bounded(value).context[key], value, key);
  for (const value of [-1, maximum + 1, 0.5, Infinity, NaN, String(maximum)]) {
    assert.equal(bounded(value).context[key], undefined, key);
  }
}
const manyCleanup = roundtrip({ ...source, secondary: Array.from({ length: 20 }, () => ({
  code: 'CLEANUP_FAILED', stage: 'cleanup', message: sentinel
})) });
assert.equal(manyCleanup.code, 'REGEX_SYNTAX');
assert.deepEqual(manyCleanup.secondary, Array.from({ length: 2 }, () => ({ code: 'CLEANUP_FAILED', stage: 'cleanup' })));
assert.doesNotMatch(JSON.stringify(manyCleanup), /PRIVATE_SENTINEL/);

// Real Python producers and exact embedded adapters, through Worker serializer/client.
const invoke = (payload) => JSON.parse(execFileSync(process.env.PYTHON || 'python',
  ['tests/web/helpers/structured-error-oracle.py'], { input: JSON.stringify(payload), encoding: 'utf8' }));
const pattern = '😀[';
const renderError = invoke({ render: true, pattern }).error;
assert.equal(renderError.code, 'REGEX_SYNTAX');
assert.equal(renderError.context.position, 1);
assert.equal(renderError.operation, 'generate');
assert.deepEqual(roundtrip(renderError).context, renderError.context);
assert.deepEqual(normalizeGenerationResponse({ error: roundtrip(renderError) }).results.error.context, renderError.context);
for (const kind of ['color', 'label']) {
  for (const features of [[], [{ type: 'CDS', qualifiers: { product: ['unrelated'] }, selector: {}, record: 'private' }]]) {
    const rule = kind === 'color' ? { feat: 'CDS', qual: 'product', val: pattern }
      : { recordId: '*', featureType: 'CDS', qualifier: 'product', valueRegex: pattern };
    const args = [JSON.stringify(features), JSON.stringify([rule]), kind];
    let destroys = 0;
    const proxy = (_helperName, ...arguments_) => JSON.stringify(invoke({ helper: _helperName, args: arguments_ }));
    proxy.destroy = () => { destroys += 1; throw new Error(sentinel); };
    assert.throws(() => callJsonHelper({ globals: { get: (name) => {
      assert.equal(name, 'call_web_json_helper'); return proxy;
    } } }, 'evaluate_rules_json', args), (error) => {
      const model = roundtrip(error);
      assert.equal(model.code, 'REGEX_SYNTAX');
      assert.equal(model.stage, 'rule-validation');
      assert.equal(model.operation, 'evaluateRules');
      assert.equal(model.context.position, 1);
      if (kind === 'label') assert.equal(model.context.row, 1);
      assert.deepEqual(model.secondary, [{ code: 'CLEANUP_FAILED', stage: 'cleanup' }]);
      assert.doesNotMatch(JSON.stringify(model), /PRIVATE_SENTINEL|😀|unrelated/);
      return true;
    });
    assert.equal(destroys, 1);
  }
}
const cleanupOnlyProxy = () => '{}';
cleanupOnlyProxy.destroy = () => { throw new Error(sentinel); };
assert.throws(() => callJsonHelper({ globals: { get: () => cleanupOnlyProxy } }, 'measure_legend_text_json', []),
  (error) => error.code === 'CLEANUP_FAILED' && error.stage === 'cleanup');
const primary = { ...source, secondary: [] };
assert.equal(resolveGenerationCleanupOutcome({ result: { error: primary }, destroyError: new Error(sentinel) }).error, primary);
assert.equal(roundtrip(primary).code, 'REGEX_SYNTAX');
assert.deepEqual(roundtrip(primary).secondary, [{ code: 'CLEANUP_FAILED', stage: 'cleanup' }]);

// Existing native validators feed finite corrections; changing a validator's
// failure template must be detected rather than silently becoming UNKNOWN.
const options = await import('../../gbdraw/web/js/services/current-option-values.js');
for (const [validator, input, field, reason] of [
  ['requireCurrentProteinBlastpMaxHits', 0, 'protein_blastp_max_hits', 'POSITIVE_INTEGER'],
  ['requireCurrentProteinBlastpCandidateLimit', -1, 'protein_blastp_candidate_limit', 'POSITIVE_OR_AUTO'],
  ['requireCurrentCollinearMinAnchors', 0, 'collinear_min_anchors', 'POSITIVE_INTEGER'],
  ['requireCurrentCollinearMaxUnitGap', -1, 'collinear_max_unit_gap', 'NONNEGATIVE_INTEGER'],
  ['requireCurrentCollinearMaxDiagonalDrift', -1, 'collinear_max_diagonal_drift', 'NONNEGATIVE_INTEGER'],
  ['requireCurrentCollinearMaxConflicts', -1, 'collinear_max_conflicts', 'NONNEGATIVE_INTEGER'],
  ['requireCurrentCollinearMaxParalogLinks', 0, 'collinear_max_paralog_links_per_orthogroup', 'POSITIVE_INTEGER'],
  ['requireCurrentCollinearInferOrthogroups', 'PRIVATE_VALUE', 'collinear_infer_orthogroups', 'BOOLEAN'],
  ['requireCurrentCircularMultiRecordSizeMode', 'PRIVATE_VALUE', 'circular_multi_record_size_mode', 'CIRCULAR_MULTI_RECORD_SIZE_MODE'],
  ['requireCurrentLinearTrackLayout', 'PRIVATE_VALUE', 'linear_track_layout', 'LINEAR_TRACK_LAYOUT'],
  ['requireCurrentLinearLabelPlacement', 'PRIVATE_VALUE', 'linear_label_placement', 'LINEAR_LABEL_PLACEMENT'],
  ['requireCurrentProteinBlastpMode', 'PRIVATE_VALUE', 'protein_blastp_mode', 'PROTEIN_BLASTP_MODE'],
  ['requireCurrentCollinearSearchScope', 'PRIVATE_VALUE', 'collinear_search_scope', 'COLLINEAR_SEARCH_SCOPE'],
  ['requireCurrentCollinearUnitMode', 'PRIVATE_VALUE', 'collinear_unit_mode', 'COLLINEAR_UNIT_MODE'],
  ['requireCurrentCollinearAnchorMode', 'PRIVATE_VALUE', 'collinear_anchor_mode', 'COLLINEAR_ANCHOR_MODE'],
  ['requireCurrentCollinearMergeOrientation', 'PRIVATE_VALUE', 'collinear_merge_orientation', 'COLLINEAR_MERGE_ORIENTATION'],
  ['requireCurrentCollinearColorMode', 'PRIVATE_VALUE', 'collinear_color_mode', 'COLLINEAR_COLOR_MODE'],
  ['requireCurrentOrthogroupMembershipMode', 'PRIVATE_VALUE', 'orthogroup_membership_mode', 'ORTHOGROUP_MEMBERSHIP_MODE']
]) assert.throws(() => options[validator](input), (error) => {
  const model = roundtrip(error);
  assert.equal(model.code, 'INPUT_INVALID', validator);
  assert.deepEqual(model.context, { field, reason }, validator);
  assert.equal(model.stage, 'request-validation');
  assert.doesNotMatch(JSON.stringify(model), /PRIVATE_/);
  assert.match(model.summary, /Use|Choose/);
  return true;
});
const { validateAnnotationRecordTargets } = await import('../../gbdraw/web/js/app/annotations/validation.js');
const annotationError = validateAnnotationRecordTargets([{ id: 'PRIVATE_SET', annotations: [{ id: 'PRIVATE_ID',
  target: { kind: 'coordinateSpan', start: 0, end: 1 } }] }], { records: [] });
const annotationModel = roundtrip(diagnosticError(annotationError.code, annotationError.context));
assert.equal(annotationModel.code, 'ANNOTATION_TARGET');
assert.equal(annotationModel.context.reason, 'POSITIVE_INTEGER');
assert.doesNotMatch(JSON.stringify(annotationModel), /PRIVATE_/);
const { runDiagramHelperOperation } = await import('../../gbdraw/web/js/services/diagram-generation.js');
await assert.rejects(runDiagramHelperOperation('PRIVATE_OPERATION'), error =>
  error.code === 'HELPER_PROTOCOL' && error.stage === 'request-validation' && !JSON.stringify(error).includes('PRIVATE_'));
await assert.rejects(runDiagramHelperOperation('evaluateRules', []), error =>
  error.code === 'HELPER_PROTOCOL' && error.operation === 'evaluateRules');

assert.equal(serializeError(null, { operation: 'generate', stage: 'render' }).code, 'UNKNOWN');
// Keep independent Python and JS finite domains aligned without importing the
// privileged Worker protocol into the public wording owner.
// X-01 vocabulary parity: every code, reason and context key a Python producer
// or the Python adapter emits is defined by the JS wording owner.
const pythonContract = JSON.parse(execFileSync(process.env.PYTHON || 'python', ['-c', `
import json
from gbdraw.web_support import error_adapter as a
exact = list(a._EXACT.values())
print(json.dumps(dict(
    operations=sorted(a.OPERATIONS), stages=sorted(a.STAGES), fields=sorted(a.FIELDS),
    codes=sorted(a.DIAGNOSTIC_CODES | {code for code, _ in exact} | {row[1] for row in a._TEMPLATES}),
    reasons=sorted(a.DIAGNOSTIC_REASONS | {ctx["reason"] for _, ctx in exact if "reason" in ctx}
        | {row[2] for row in a._TEMPLATES} | set(a._CONSTRAINTS.values())),
    contextKeys=sorted(a._DIAGNOSTIC_INTEGER_KEYS | {"configPath"}),
    exceptionTypes=sorted(a.EXCEPTION_TYPE_NAMES))))
`], { encoding: 'utf8' }));
for (const code of pythonContract.codes) assert.equal(roundtrip({ code }).code, code, code);
for (const reason of pythonContract.reasons) {
  const model = roundtrip({ code: 'INPUT_INVALID', context: { reason } });
  assert.equal(model.context.reason, reason, reason);
  assert.notEqual(model.summary, normalizeUserFacingError({ code: 'INPUT_INVALID' }).summary, reason);
}
for (const key of pythonContract.contextKeys) {
  const value = key === 'configPath' ? 'objects.scale.interval' : 3;
  assert.equal(roundtrip({ code: 'INPUT_INVALID', context: { [key]: value } }).context[key], value, key);
}
const { DIAGRAM_HELPER_OPERATION_NAMES } = await import('../../gbdraw/web/js/services/diagram-worker-protocol.js');
assert.deepEqual(pythonContract.operations, ['unknown', 'generate', 'align', 'feature-extraction', 'export-svg', 'export-png', 'export-pdf', ...DIAGRAM_HELPER_OPERATION_NAMES].sort());
for (const operation of pythonContract.operations) assert.equal(roundtrip({ code: 'UNKNOWN', operation }).operation, operation);
for (const stage of pythonContract.stages) assert.equal(roundtrip({ code: 'UNKNOWN', stage }).stage, stage);
for (const field of pythonContract.fields) assert.equal(roundtrip({ code: 'INPUT_INVALID', context: { field } }).context.field, field);
for (const exceptionType of pythonContract.exceptionTypes) {
  assert.equal(roundtrip({ code: 'RENDER_FAILED', context: { exceptionType } }).context.exceptionType, exceptionType, exceptionType);
}
for (const exceptionType of ['PRIVATE_Error', 'valueerror', 7]) {
  assert.equal(roundtrip({ code: 'RENDER_FAILED', context: { exceptionType } }).context.exceptionType, undefined);
}

const { validateCustomTrackPlan, assertValidCustomTrackPlan, CustomTrackPlanValidationError } =
  await import('../../gbdraw/web/js/services/track-slot-validation.js');
const trackPlan = validateCustomTrackPlan({ mode: 'linear', axisIndex: 0, annotationSetIds: [], slots: [
  { id: 'PRIVATE_FEATURES', renderer: 'features', enabled: true, side: 'overlay', params: {} },
  { id: 'PRIVATE_ANNOTATION', renderer: 'annotations', enabled: true, side: 'overlay',
    params: { set_id: 'PRIVATE_SET', anchor_slot: 'PRIVATE_FEATURES', layer: 'foreground' } }
] });
assert.throws(() => assertValidCustomTrackPlan(trackPlan), error => {
  const model = roundtrip(error);
  assert.equal(model.code, 'TRACK_INVALID');
  assert.deepEqual(model.context, { reason: 'ANNOTATION_SET', field: 'set_id', slotIndex: 1 });
  assert.match(model.summary, /Select an existing annotation set/);
  assert.doesNotMatch(JSON.stringify(model), /PRIVATE_/);
  return true;
});
assert.equal(roundtrip(new CustomTrackPlanValidationError([{ code: 'new-unmapped-issue', message: 'PRIVATE_MSG' }])).code, 'VALIDATION_UNCLASSIFIED');

const glyphModel = roundtrip(new Error('PDF fonts do not contain U+4E00. Use SVG to retain this text.'));
assert.equal(glyphModel.code, 'PDF_GLYPH');
assert.equal(glyphModel.stage, 'font-validation');
assert.deepEqual(glyphModel.context, { codepoint: 0x4e00 });
assert.match(glyphModel.summary, /Use SVG to retain this text/);
assert.deepEqual(glyphModel.actions, ['use-svg']);
assert.equal(roundtrip(new Error('The selected PNG DPI is invalid.')).code, 'PNG_DPI');

const alignCause=normalizeUserFacingError(source,{operation:'align',stage:'render'});
assert.equal(alignCause.code,'REGEX_SYNTAX');
assert.equal(alignCause.stage,'rule-validation');
assert.equal(alignCause.operation,'align');
assert.deepEqual(alignCause.context,normalized.context);
assert.deepEqual(alignCause.secondary,normalized.secondary);
assert.deepEqual(normalizeUserFacingError(alignCause),alignCause);
assert.equal((alignCause.summary.match(/Python regular expression is invalid/g)||[]).length,1);

// Native boundary diagnostics retain corrections while discarding private IDs.
for (const [message, code, context] of [
  ['Circular region End (60000) exceeds the selected record length (50466).', 'REGION_INVALID', { field: 'region', reason: 'RECORD_BOUNDS' }],
  ['Circular region requires both Start and End coordinates.', 'REGION_INVALID', { field: 'region', reason: 'BOTH_ENDPOINTS' }],
  ['Circular region Start and End must be positive integers.', 'REGION_INVALID', { field: 'region', reason: 'POSITIVE_INTEGER' }],
  ['Circular region Start must not exceed End. Use Reverse complement to change display orientation.', 'REGION_INVALID', { field: 'region', reason: 'ORDER' }],
  ['The preserved comparison is missing a required resource.', 'COMPARISON_INPUT', { field: 'comparison', reason: 'REQUIRED' }],
  ['Canonical renderRequest records are required.', 'INPUT_INVALID', { field: 'schema', reason: 'RECORDS_REQUIRED' }],
  ['Session resource PRIVATE_RESOURCE has an invalid declared byte size.', 'INPUT_INVALID', { field: 'schema', reason: 'RESOURCE_SIZE' }],
  ['Custom Track Slots use an obsolete schema. Recreate the slots with schema version 2.', 'INPUT_INVALID', { field: 'schema', reason: 'TRACK_SCHEMA' }],
  ['Session contains unclassified top-level field(s): PRIVATE_FIELD', 'INPUT_INVALID', { field: 'schema', reason: 'SESSION_FIELDS' }]
]) {
  const model = roundtrip(new Error(message));
  assert.equal(model.code, code);
  assert.deepEqual(model.context, context);
  assert.equal(model.stage, 'request-validation');
  assert.doesNotMatch(JSON.stringify(model), /PRIVATE_/);
}
assert.equal(roundtrip(new Error('Circular region End (PRIVATE_VALUE) exceeds the selected record length (50466).')).code, 'UNKNOWN');

// Use actual track validators so producer wording and the public projection
// cannot drift; every physical field keeps its own bound and unit guidance.
for (const [mode, field, reason] of [
  ['linear', 'height', 'PIXEL_POSITIVE'], ['linear', 'spacing', 'PIXEL_NONNEGATIVE'],
  ['circular', 'inner_gap_px', 'PIXEL_NONNEGATIVE'], ['circular', 'outer_gap_px', 'PIXEL_NONNEGATIVE'],
  ['circular', 'radius', 'POSITIVE_SCALAR'], ['circular', 'width', 'POSITIVE_SCALAR']
]) {
  const plan = validateCustomTrackPlan({ mode, axisIndex: 0, annotationSetIds: [], slots: [{
    id: 'PRIVATE_SLOT', renderer: 'features', enabled: true, side: 'overlay', params: {}, [field]: '-1px'
  }] });
  assert.throws(() => assertValidCustomTrackPlan(plan), error => {
    const model = roundtrip(error);
    assert.equal(model.code, 'TRACK_INVALID');
    assert.deepEqual(model.context, { field, reason, slotIndex: 0 });
    assert.doesNotMatch(JSON.stringify(model), /PRIVATE_/);
    assert.match(model.summary, /px/);
    return true;
  });
}

const arrows = await import('../../gbdraw/web/js/utils/feature-rendering.js');
for (const [validator, field, reason] of [
  ['normalizeArrowHeadLengthRatio', 'arrow_head_length_ratio', 'POSITIVE_OR_AUTO'],
  ['normalizeArrowShaftWidthRatio', 'arrow_shaft_width_ratio', 'POSITIVE_UNIT_INTERVAL']
]) assert.throws(() => arrows[validator](0), error => {
  const model = roundtrip(error);
  assert.equal(model.code, 'INPUT_INVALID');
  assert.deepEqual(model.context, { field, reason });
  assert.match(model.summary, /Arrow/);
  return true;
});

const busy = roundtrip(new Error('A diagram generation request is already running.'));
assert.equal(busy.code, 'GENERATION_BUSY');
assert.equal(busy.stage, 'render');
assert.match(busy.summary, /Wait for it to finish/);
assert.deepEqual(busy.actions, ['retry', 'save-session']);

assert.equal(roundtrip(Object.assign(new Error('This browser does not support Session import Workers.'),
  { code: 'SESSION_IMPORT_UNAVAILABLE', stage: 'transport' })).code, 'SESSION_IMPORT_UNAVAILABLE');

// X-01 (PD-OI-046, R6): a user-facing validation producer yields a recognized
// diagnostic, not the UNKNOWN or unclassified fallback.
test('Circular geometry shortcut rejection is a recognized diagnostic (X-01)', () => {
  assert.throws(() => normalizeCircularGeometryShortcuts({ featureWidth: 0 }), (error) => {
    const model = normalizeUserFacingError(error, { operation: 'generate', stage: 'request-validation' });
    assert.equal(model.code, 'INPUT_INVALID');
    assert.deepEqual(model.context, { field: 'feature_width_circular', reason: 'POSITIVE_OR_AUTO' });
    assert.match(model.summary, /Field: Feature Width\. Use Auto or a finite value greater than zero\./);
    return true;
  });
});

// N-12: an error whose guidance is to Generate again offers that action, and a
// failed Save never offers Save Session again.
test('the operation error panel offers Generate and no Save after a failed Save (N-12)', () => {
  const exportInput = normalizeUserFacingError({ code: 'EXPORT_INPUT' }, { operation: 'export-png' });
  assert.ok(exportInput.actions.includes('generate'));
  const indexHtml = readFileSync(new URL('../../gbdraw/web/index.html', import.meta.url), 'utf8');
  // E1: an error about the other mode's Result runs Generate in that mode.
  assert.match(indexHtml, /errorDisplay\.actions\.includes\('generate'\) && errorDisplay\.operation !== 'generate'" type="button" @click="generateFromError"/);
  const hiddenMode = normalizeUserFacingError(diagnosticError(
    'SESSION_SAVE_REQUIRES_GENERATE', { diagramMode: 'linear' }, { operation: 'session-save', stage: 'result-admission' }
  ));
  assert.equal(hiddenMode.context.diagramMode, 'linear');
  assert.match(hiddenMode.summary, / Diagram: Linear\.$/);
  assert.deepEqual(hiddenMode.actions, ['generate']);
  assert.match(indexHtml, /errorDisplay\.actions\.includes\('save-session'\) && errorDisplay\.operation !== 'session-save'/);
  const saveFailure = normalizeUserFacingError(new Error('PRIVATE'), { operation: 'session-save' });
  assert.equal(saveFailure.code, 'UNKNOWN');
  assert.equal(saveFailure.operation, 'session-save');
});

// R6: locators the producer keeps are shown in the summary, wording only here.
test('summary shows Sequence, Line, Track row, Depth series, band and setting locators', () => {
  for (const [source, pattern] of [
    [diagnosticError('INPUT_REQUIRED', { inputOrdinal: 2 }), /^Supply GenBank input or matching GFF3 and FASTA inputs\. Sequence 2\.$/],
    [diagnosticError('COMPARISON_INPUT', { inputOrdinal: 1, reason: 'REQUIRED' }), / Comparison sequence 1\. Supply the required value\.$/],
    [diagnosticError('DECORATION_CONTINUITY', { inputOrdinal: 3, field: 'scale', reason: 'DECORATION_TARGET' }), / Result 3\. Field: scale\./],
    [diagnosticError('TABLE_INVALID', { row: 3, field: 'color', reason: 'COLOR' }), /^The table is invalid\. Line 3\. Field: color\. Use none/],
    [{ code: 'COMPARISON_INPUT', context: { reason: 'NONNEGATIVE_INTEGER', row: 4, column: 7 } }, / Line 4\. Column 7\. Use an integer of zero or greater\.$/],
    [{ code: 'TRACK_LAYOUT', context: { reason: 'CANNOT_FIT', slotIndex: 1, innerPx: 181, outerPx: 209 } },
      /^A circular track does not fit\. Track row 2\. Move the track.*the other side of the Axis\. Available band: 181–209 px\.$/],
    [{ code: 'DEPTH_INVALID', context: { reason: 'DEPTH_VALUES', seriesIndex: 0 } }, / Depth series 1\. Use integer positions/],
    [{ code: 'INPUT_INVALID', context: { reason: 'INTEGER', configPath: 'objects.scale.interval' } }, / Setting: objects\.scale\.interval\. Use an integer\.$/]
  ]) assert.match(roundtrip(source).summary, pattern);
  for (const configPath of ['PRIVATE VALUE', 'a', 'Objects.scale', `objects.${'x'.repeat(80)}`, 7]) {
    assert.equal(roundtrip({ code: 'INPUT_INVALID', context: { configPath } }).context.configPath, undefined);
  }
});

// X-01: the producer contract message is a fixed identifier, never a document value.
test('a falsy caught value normalizes to the unknown failure, not to no error (OV-70)', () => {
  assert.equal(normalizeUserFacingError(undefined), null);
  for (const value of [undefined, null, 0, '', false]) {
    const model = normalizeCaughtError(value, { stage: 'request-validation' });
    assert.deepEqual([model.code, model.stage], ['UNKNOWN', 'request-validation']);
    assert.match(model.summary, /failed without recognized diagnostic information/);
  }
  const thrown = diagnosticError('REGION_INVALID', { inputOrdinal: 1, reason: 'SELECT_RECORD_FOR_REGION' });
  assert.deepEqual(normalizeCaughtError(thrown, { stage: 'request-validation' }), normalizeUserFacingError(thrown, { stage: 'request-validation' }));
});

test('diagnosticError carries a fixed message, stage, operation and bounded context', () => {
  const error = diagnosticError('REGION_INVALID', { inputOrdinal: 1, reason: 'SELECT_RECORD_FOR_REGION' });
  assert.equal(error.message, 'REGION_INVALID/SELECT_RECORD_FOR_REGION');
  assert.equal(error.stage, 'request-validation');
  assert.equal(Object.hasOwn(error, 'operation'), false);
  const save = diagnosticError('SESSION_SAVE_REQUIRES_GENERATE', {}, { operation: 'session-save', stage: 'result-admission' });
  const model = normalizeUserFacingError(save);
  assert.deepEqual([model.code, model.operation, model.stage, model.actions], ['SESSION_SAVE_REQUIRES_GENERATE', 'session-save', 'result-admission', ['generate']]);
  assert.match(model.summary, /Generate once, then Save Session\./);
  // Literal codes only: no producer interpolates a document value into the code.
  // Pass-through forms re-raise a code another producer already chose.
  const passThrough = new Set(['issue.code', 'annotationError.code', 'model.code', 'code']);
  const files = execFileSync('git', ['grep', '-l', 'diagnosticError(', '--', 'gbdraw/web/js'], { encoding: 'utf8' })
    .split('\n').filter(Boolean);
  let literal = 0;
  for (const file of files) {
    const text = readFileSync(file, 'utf8');
    for (const match of text.matchAll(/diagnosticError\(\s*([^,)]*)/g)) {
      if (/^'[A-Z][A-Z0-9_]*'$/.test(match[1])) { literal += 1; continue; }
      assert.ok(passThrough.has(match[1].trim()), `${file}: diagnosticError(${match[1]}`);
    }
  }
  assert.ok(literal >= 20, literal);
});

test('a Circular placement limited by the center reservation names that cause', () => {
  for (const [reason, pattern] of [
    ['CANNOT_FIT', /Move the track, reduce widths/],
    ['DEFINITION_RESERVED', /center definition text limits the inside tracks\. Shorten Species or Strain, reduce Default font size, set a smaller Center Reserved Radius/],
    ['CENTER_RESERVED', /The Center Reserved Radius limits the inside tracks\. Set a smaller Center Reserved Radius or place tracks outside\./]
  ]) {
    const layout = roundtrip({ code: 'TRACK_LAYOUT', operation: 'generate', stage: 'render', context: { reason } });
    assert.equal(layout.code, 'TRACK_LAYOUT');
    assert.equal(layout.context.reason, reason);
    assert.match(layout.summary, /^A circular track does not fit\./);
    assert.match(layout.summary, pattern);
    assert.deepEqual(layout.actions, ['edit-track', 'retry']);
  }
});

// OV-09 (R6): Python reports the request row of an unsupported Feature placement
// lane; the Web adds that feature's caption and the summary names it.
test('an unsupported Feature placement lane names the feature and the ways out', () => {
  for (const [reason, pattern] of [['SPLIT_LANES', /or use split feature lanes \(Track Preset Middle\)\.$/],
    ['OVERLAY_LANES', /or use Track Layout Features on axis with Separate Strands off\.$/]]) {
    const engine = roundtrip({ code: 'FEATURE_PLACEMENT', operation: 'generate', stage: 'render', context: { reason, placementIndex: 4 } });
    assert.deepEqual([engine.code, engine.context], ['FEATURE_PLACEMENT', { reason, placementIndex: 4 }]);
    const named = normalizeUserFacingError({ ...engine, context: { ...engine.context, featureCaption: ' NADH\ndehydrogenase ' } });
    assert.match(named.summary, /^A Feature placement uses a lane that the current feature track does not have\. Feature: NADH dehydrogenase\. Set that feature's Feature placement to Auto or Main, /);
    assert.match(named.summary, pattern);
    assert.deepEqual(normalizeUserFacingError(named), named);
    assert.deepEqual(named.actions, ['edit-track', 'retry']);
  }
  const long = normalizeUserFacingError({ code: 'FEATURE_PLACEMENT', context: { featureCaption: 'x'.repeat(200) } });
  assert.equal(long.context.featureCaption.length, 80);
  assert.equal(normalizeUserFacingError({ code: 'FEATURE_PLACEMENT', context: { featureCaption: 7 } }).context.featureCaption, undefined);
});

// B9 (P07: offer only working actions): an unclassified failure inside the
// engine's render or result stage repeats for the same inputs. It is reported as
// a render failure naming the Python exception class, with Save Session and no
// Retry, instead of the input-validation text.
test('an unclassified engine failure is a render failure with its exception class and no Retry (B9)', () => {
  const validation = normalizeUserFacingError({ code: 'VALIDATION_UNCLASSIFIED' }).summary;
  const cases = [['ValidationError', 'render'], ['ParseError', 'render'], ['ValueError', 'render'],
    ['IndexError', 'render'], ['ZeroDivisionError', 'render'], ['KeyError', 'result-admission'], ['MemoryError', 'render']];
  const models = invoke({ raise: cases }).map(({ error }) => roundtrip(error));
  const memory = models.pop();
  for (const [index, model] of models.entries()) {
    const [exceptionType, stage] = cases[index];
    assert.deepEqual([model.code, model.operation, model.stage], ['RENDER_FAILED', 'generate', stage], exceptionType);
    assert.deepEqual(model.context, { exceptionType });
    assert.notEqual(model.summary, validation);
    assert.doesNotMatch(model.summary, /Input validation failed/);
    assert.match(model.summary, new RegExp(`^The diagram engine failed while drawing this diagram\\. .* Python exception: ${exceptionType}\\.$`));
    assert.match(model.details[0].text, new RegExp(`\\nexceptionType: ${exceptionType}$`));
    assert.equal(model.actions.includes('retry'), false, exceptionType);
    assert.deepEqual(model.actions, ['save-session']);
    assert.doesNotMatch(JSON.stringify(model), /PRIVATE_/);
  }
  // Memory depends on the runtime, not the inputs: that failure keeps Retry.
  assert.deepEqual([memory.code, memory.actions], ['UNKNOWN', ['retry', 'save-session']]);
});

// B10: two user-fixable engine failures carry a producer diagnostic through the
// real render wrapper, so they are input errors with a working action, not a
// RENDER_FAILED that offers only Save Session.
test('a display start beyond the record and an unparsable GenBank file are input errors (B10)', () => {
  const [display, genbank] = [invoke({ displayStart: 99999 }).error, invoke({ genbankText: 'LOCUS       PRIVATE_BROKEN\n     CDS  PRIVATE\n' }).error]
    .map((error) => roundtrip(error));
  assert.deepEqual([display.code, display.operation, display.stage], ['INPUT_INVALID', 'generate', 'render']);
  assert.deepEqual(display.context, { field: 'start', reason: 'DISPLAY_START_BOUNDS' });
  assert.match(display.summary, /Use a display start between 1 and the record length\./);
  assert.deepEqual(display.actions, ['edit-input', 'retry']);
  assert.deepEqual([genbank.code, genbank.operation, genbank.stage], ['INPUT_UNREADABLE', 'generate', 'render']);
  assert.deepEqual(genbank.context, {});
  assert.deepEqual(genbank.actions, ['select-input', 'retry']);
  for (const model of [display, genbank]) {
    assert.doesNotMatch(JSON.stringify(model), /PRIVATE_|RENDER_FAILED/);
  }
});

// OV-130 (R6): a preserved setting of the other diagram mode fails Generate and
// Save Session. Python names the first setting path and the mode it belongs to,
// so the summary names both and the two ways out, instead of the unclassified
// input-validation text. Retry repeats the same request, so it is not offered.
test('a preserved setting of the other diagram mode names the setting and the ways out (OV-130)', () => {
  const ways = (mode) => ` It applies only to ${mode} diagrams. Reset it under Preserved session settings, or switch to ${mode}.`;
  const generate = roundtrip(invoke({ configOverrides: { 'objects.blast_match.curve_tension': 0.3 } }).error);
  assert.deepEqual([generate.code, generate.operation, generate.stage], ['MODE_SETTING', 'generate', 'request-validation']);
  assert.deepEqual(generate.context, { reason: 'LINEAR_SETTING', configPath: 'objects.blast_match.curve_tension' });
  assert.equal(generate.summary, 'A preserved session setting does not apply to this diagram mode.'
    + ` Setting: objects.blast_match.curve_tension.${ways('Linear')}`);
  assert.deepEqual(generate.actions, ['edit-input']);
  // Save Session validates the preserved settings through the helper channel.
  const save = roundtrip(invoke({ helper: 'validate_web_config_overrides_json', args: ['linear', 'null',
    JSON.stringify({ 'objects.ticks.tick_width': 4, 'canvas.circular.radius': 1.2 }), '[]', true] }).error);
  assert.deepEqual([save.code, save.operation, save.stage], ['MODE_SETTING', 'validateConfigOverrides', 'helper']);
  assert.deepEqual(save.context, { reason: 'CIRCULAR_SETTING', configPath: 'canvas.circular.radius' });
  assert.equal(save.summary, `A preserved session setting does not apply to this diagram mode. Setting: canvas.circular.radius.${ways('Circular')}`);
  for (const model of [generate, save]) assert.doesNotMatch(model.summary, /Input validation failed/);
});

// B11 (R6: an offered action must be one that can work): a live edit rerenders
// the committed Session with the current editor tables, so repeating the edit
// sends the same request. The note above the Result comes from the failure's
// normalized model and offers Retry only for a failure that changes no input.
test('a live edit failure note offers Retry only when the same request can succeed (B11)', async () => {
  const html = readFileSync(new URL('../../gbdraw/web/index.html', import.meta.url), 'utf8');
  const note = html.match(/<p v-if="labelReflowProcessing \|\| labelReflowLastError"[^>]*>([^<]*)<\/p>/)[1];
  assert.doesNotMatch(note, /Retry the live edit/, 'the failure note is not fixed text');
  assert.match(note, /: labelReflowLastError\.note \}\}/);
  const { liveEditFailure } = await import('../../gbdraw/web/js/utils/error-normalization.js');
  assert.equal(liveEditFailure(null), null);
  const prefix = /^Live edit failed: direct edits already applied are kept; geometry may still need updating\. /;
  const repeating = [
    roundtrip({ code: 'RENDER_FAILED', operation: 'generate', stage: 'render', context: { exceptionType: 'ValidationError' } }),
    roundtrip({ code: 'INPUT_INVALID', operation: 'generate', stage: 'render', context: { field: 'start', reason: 'DISPLAY_START_BOUNDS' } }),
    roundtrip({ code: 'TRACK_LAYOUT', operation: 'generate', stage: 'render', context: { reason: 'CANNOT_FIT' } }),
    normalizeUserFacingError({ code: 'DECORATION_CONTINUITY', stage: 'result-admission', context: { field: 'title' } })
  ];
  for (const error of repeating) {
    const model = liveEditFailure(error);
    assert.deepEqual([model.code, model.actions], [error.code, ['generate']], error.code);
    assert.match(model.note, prefix);
    assert.ok(model.note.includes(error.summary), error.code);
    assert.doesNotMatch(model.note, /Retry the live edit/, error.code);
    assert.match(model.note, / Change the edit, or change the settings and use Generate\.$/);
  }
  const transient = [
    normalizeUserFacingError(new Error('PRIVATE_ worker transport failure'), { operation: 'generate', stage: 'render' }),
    normalizeUserFacingError({ code: 'WORKER_INIT', stage: 'initialization' }),
    normalizeUserFacingError({ code: 'GENERATION_BUSY' })
  ];
  for (const error of transient) {
    const model = liveEditFailure(error);
    assert.deepEqual([model.code, model.actions], [error.code, ['retry', 'generate']], error.code);
    assert.match(model.note, prefix);
    assert.match(model.note, new RegExp(` Retry the live edit or use Generate\\.$`));
    assert.doesNotMatch(JSON.stringify(model), /PRIVATE_/);
  }
  assert.equal(transient[0].code, 'UNKNOWN');
  const generate = liveEditFailure(normalizeUserFacingError(diagnosticError('LIVE_EDIT_REQUIRES_GENERATE')));
  assert.deepEqual(generate.actions, ['generate']);
  assert.equal(generate.note, 'Live edit failed: direct edits already applied are kept; geometry may still need updating. Generate the diagram to update its label placement.');
});

// OV-05..OV-07, OV-13 (R6): a Label visibility On that the diagram does not
// draw (an underlay feature, Embedded Only, a duplicate feature) names the
// feature and the popup choice that ends it, instead of UNKNOWN. Two labels
// bound to one feature repeat for the same inputs, so they are a render
// failure that names the feature and offers Save Session.
test('a forced label that the diagram does not draw names the feature (OV-06)', async () => {
  const { requireUniqueEditableLabelBindings } = await import('../../gbdraw/web/js/app/feature-editor/label-actions.js');
  const { liveEditFailure } = await import('../../gbdraw/web/js/utils/error-normalization.js');
  const label = (featureId) => ({ getAttribute: (name) => (name === 'data-label-feature-id' ? featureId : null) });
  const features = [
    { svg_id: 'f59dc64fc', type: 'repeat_region', start: 1000, end: 1600, note: 'PRIVATE_NOTE' },
    { svg_id: 'f1aff4c2b__instance_5_ef2d127de37b942b', type: 'CDS', start: 3000, end: 3600, product: 'PRIVATE_PRODUCT' }
  ];
  const failure = (labels, required, options = {}) => {
    try {
      requireUniqueEditableLabelBindings(labels, required, { features, ...options });
    } catch (error) {
      return roundtrip(error);
    }
    return null;
  };
  const missing = failure([label('f9b3094c1')], ['f59dc64fc']);
  assert.deepEqual([missing.code, missing.operation, missing.stage, missing.actions],
    ['LABEL_NOT_DRAWN', 'generate', 'render', ['edit-input']]);
  assert.deepEqual(missing.context, { reason: 'FORCED_LABEL', featureId: 'f59dc64fc',
    featureType: 'repeat_region', featureStart: 1001, featureEnd: 1600 });
  assert.equal(missing.summary, 'Label visibility is On for a feature, but this diagram does not draw its label.'
    + ' Feature: repeat_region 1001..1600, ID f59dc64fc.'
    + " Open the feature's popup and set Label visibility to Default, or change the setting that prevents its label.");
  assert.match(missing.details[0].text, /\nfeatureId: f59dc64fc\n/);
  const live = liveEditFailure(missing);
  assert.deepEqual(live.actions, ['generate']);
  assert.ok(live.note.includes(missing.summary));

  const several = failure([], ['f1aff4c2b__instance_5_ef2d127de37b942b', 'f59dc64fc']);
  assert.equal(several.code, 'LABEL_NOT_DRAWN');
  assert.match(several.summary, / Feature: CDS 3001\.\.3600, ID f1aff4c2b__instance_5_ef2d127de37b942b\. Features affected: 2\. Open/);
  assert.equal(failure([], ['f59dc64fc'], { allowMissing: true }), null);
  assert.equal(failure([label('f59dc64fc')], ['f59dc64fc']), null);

  const unknownFeature = failure([], ['f0000dead']);
  assert.deepEqual(unknownFeature.context, { reason: 'FORCED_LABEL', featureId: 'f0000dead' });
  assert.match(unknownFeature.summary, / Feature: ID f0000dead\. Open/);

  for (const allowMissing of [false, true]) {
    const ambiguous = failure([label('f59dc64fc'), label('f59dc64fc')], ['f59dc64fc'], { allowMissing });
    assert.deepEqual([ambiguous.code, ambiguous.stage, ambiguous.actions], ['RENDER_FAILED', 'render', ['save-session']]);
    assert.deepEqual(ambiguous.context, { featureId: 'f59dc64fc', featureType: 'repeat_region',
      featureStart: 1001, featureEnd: 1600 });
    assert.match(ambiguous.summary, / Feature: repeat_region 1001\.\.1600, ID f59dc64fc\.$/);
  }
  for (const model of [missing, several, live]) assert.doesNotMatch(JSON.stringify(model), /PRIVATE_/);

  // The feature locator is bounded like every other context value.
  assert.deepEqual(roundtrip({ code: 'LABEL_NOT_DRAWN', context: {
    featureId: 'PRIVATE VALUE', featureType: `CDS${'x'.repeat(80)}`, featureStart: 0, featureEnd: 1.5, featureCount: -1
  } }).context, {});
  assert.deepEqual(roundtrip({ code: 'LABEL_NOT_DRAWN', context: { featureType: "5'UTR", featureStart: 1, featureEnd: 20000000 } }).context,
    { featureType: "5'UTR", featureStart: 1, featureEnd: 20000000 });
});

const producerSummary = (code, context = {}) => normalizeUserFacingError({ code, stage: 'helper', context }).summary;

test('a region wholly beyond the record keeps the region correction (CI-01)', () => {
  assert.equal(producerSummary('REGION_INVALID', { field: 'region', reason: 'RECORD_BOUNDS' }),
    'The region is invalid. Field: region. Keep the region within the record length.');
});

test('a BLAST table error names its locator and correction, not the sequence-file sentence (CI-04)', () => {
  assert.equal(producerSummary('COMPARISON_INPUT', { reason: 'FINITE', row: 2, column: 11 }),
    'The comparison input is invalid. Line 2. Column 11. Use a finite number.');
  assert.match(producerSummary('COMPARISON_INPUT', { reason: 'RECORD_ID', column: 1 }),
    /^The comparison input is invalid\. Column 1\. A table row names another displayed record in this column\. Swap the query and subject columns/);
  // A producer that names no reason keeps the comparison-source guidance.
  assert.equal(producerSummary('COMPARISON_INPUT'),
    'The comparison input is invalid. Supply a comparison sequence file (FASTA, GenBank, or DDBJ) or BLAST outfmt 6/7 as required.');
});

test('a comparison error names its record pair as the comparison panel does (CI-04)', () => {
  assert.equal(producerSummary('COMPARISON_INPUT', { reason: 'FINITE', row: 2, column: 11, queryRecordIndex: 0, subjectRecordIndex: 1 }),
    'The comparison input is invalid. Pair: #1 to #2. Line 2. Column 11. Use a finite number.');
  // The locator is bounded like the other indexes and needs both endpoints.
  assert.doesNotMatch(producerSummary('COMPARISON_INPUT', { reason: 'FINITE', queryRecordIndex: 0 }), /Pair/);
  assert.doesNotMatch(producerSummary('COMPARISON_INPUT', { reason: 'FINITE', queryRecordIndex: -1, subjectRecordIndex: 1.5 }), /Pair/);
});

test('an Upload pair without its BLAST TSV names the missing file (CI-06)', () => {
  assert.match(producerSummary('COMPARISON_INPUT', { reason: 'BLAST_TSV_REQUIRED' }), /Choose a BLAST TSV for this pair/);
  assert.match(producerSummary('COMPARISON_INPUT', { reason: 'PAIR_TOPOLOGY' }), /adjacent rows/);
});

test('an outfmt 7 Fields line that cannot be read names the line (CI-07d)', () => {
  assert.match(producerSummary('COMPARISON_INPUT', { reason: 'OUTFMT7_FIELDS', row: 2 }), /^The comparison input is invalid\. Line 2\. List all 12 standard/);
});

test('a GenBank slot file without records says what it looks like (UJ-07)', () => {
  assert.equal(producerSummary('NO_RECORDS', { reason: 'FASTA_IN_GENBANK' }),
    'No records were found. Choose input containing records. This file looks like FASTA. Use GFF3 + FASTA input, or a GenBank/DDBJ flat file.');
  assert.match(producerSummary('NO_RECORDS', { reason: 'EMPTY_FILE' }), / The file is empty\.$/);
  assert.match(producerSummary('NO_RECORDS', { reason: 'NOT_GENBANK' }), / The file is not a GenBank\/DDBJ flat file: it has no record header line\.$/);
});

// TK-06: a Depth series with no file in any record names the series and offers
// both recoveries, whether the request builder or the Generate check finds it.
test('a Depth series without a file says to attach a TSV or remove the series (TK-06)', () => {
  const expected = 'The depth input or settings are invalid. Depth series 2. Attach a Depth TSV to this series, or remove the series.';
  assert.equal(producerSummary('DEPTH_INVALID', { reason: 'DEPTH_SERIES_SOURCE', seriesIndex: 1 }), expected);
  for (const message of [
    'Depth series #2 (logical track index 1) has no TSV source in any record.',
    'Depth series #2 (logical track index 1) has no TSV source in any record. Add a TSV or remove the series.'
  ]) {
    const model = roundtrip(new Error(message));
    assert.deepEqual([model.code, model.context], ['DEPTH_INVALID', { reason: 'DEPTH_SERIES_SOURCE', seriesIndex: 1 }]);
    assert.equal(model.summary, expected);
  }
});
