// Copy FASTA of Select region: the page sends the committed request (the one the shown Result was made from) with
// its resources and the spans, through the Worker helper `readRecordSpanSequences`; Python slices the bases
// (tests/test_web_record_spans.py). The page reads the reply, and the clipboard write starts with the click.
import assert from 'node:assert/strict';
import test from 'node:test';

import {
  DIAGRAM_HELPER_OPERATIONS, DIAGRAM_HELPER_STAGING_OPERATIONS, stagedHelperReferences
} from '../../gbdraw/web/js/services/diagram-worker-protocol.js';
import { readRecordSpanSequencesReply, requestRecordSpanSequences } from '../../gbdraw/web/js/services/record-span-sequences.js';

const committed = { renderRequest: { mode: 'circular', records: [] }, resources: { 'resource-1': { kind: 'web-file' } } };
const spans = [{ recordKey: 'record-1', start: 16_001, end: 501 }];
const reply = {
  spans: [{ recordKey: 'record-1', recordId: 'NC_012920.1', start: 16_001, end: 501, length: 1_070, wrapsOrigin: true,
    recordLength: 16_569,
    header: 'NC_012920.1:16001-501', sequence: 'ACGT' }],
  fasta: '>NC_012920.1:16001-501\nACGT\n'
};

test('the helper stages the resources of the committed request it is sent', () => {
  const operation = DIAGRAM_HELPER_OPERATIONS.READ_RECORD_SPAN_SEQUENCES;
  assert.equal(operation, 'readRecordSpanSequences');
  assert.ok(DIAGRAM_HELPER_STAGING_OPERATIONS.includes(operation));
  const canonicalRequest = { records: [] };
  assert.equal(stagedHelperReferences(operation, { canonicalRequest }), canonicalRequest);
});

test('Copy FASTA sends the committed request, its resources and the span, and reads the reply', async () => {
  /** @type {unknown[]} */
  const sent = [];
  const read = await requestRecordSpanSequences(committed, spans, async (payload) => {
    sent.push(payload);
    return { result: reply };
  });
  assert.deepEqual(sent, [{ canonicalRequest: committed.renderRequest, resources: committed.resources, spans }]);
  assert.equal(read.fasta, reply.fasta);
});

test('a Result without a committed request, and a malformed reply, fail with a typed code', async () => {
  let ran = false;
  await assert.rejects(requestRecordSpanSequences(null, spans, async () => { ran = true; return {}; }),
    (error) => error.code === 'RESOURCE_INVALID' && error.operation === 'readRecordSpanSequences');
  assert.equal(ran, false, 'no Python starts');
  for (const result of [null, { spans: [], fasta: '' }, { ...reply, fasta: 1 }, { ...reply, spans: [{ ...reply.spans[0], recordKey: 'other' }] }]) {
    assert.throws(() => readRecordSpanSequencesReply(result, spans),
      (error) => error.code === 'RESULT_INVALID' && error.operation === 'readRecordSpanSequences');
  }
});

test('the clipboard write starts before the text is ready, and a failed read is the caller\'s to report', async () => {
  /** @type {unknown[]} */
  const written = [];
  const restore = { ClipboardItem: globalThis.ClipboardItem, isSecureContext: globalThis.isSecureContext, navigator: globalThis.navigator };
  globalThis.ClipboardItem = class { constructor(items) { this.items = items; } };
  Object.defineProperty(globalThis, 'isSecureContext', { value: true, configurable: true });
  Object.defineProperty(globalThis, 'navigator', { configurable: true, value: { clipboard: {
    write: async ([item]) => { written.push(await (await item.items['text/plain']).text()); }
  } } });
  try {
    const { copyPendingTextToClipboard } = await import('../../gbdraw/web/js/utils/clipboard.js');
    let resolve;
    const pending = new Promise((done) => { resolve = done; });
    const copied = copyPendingTextToClipboard(pending);
    resolve(reply.fasta);
    assert.equal(await copied, reply.fasta);
    assert.deepEqual(written, [reply.fasta]);
    const failure = Object.assign(new Error('read failed'), { code: 'RESOURCE_INVALID' });
    await assert.rejects(copyPendingTextToClipboard(Promise.reject(failure)), (error) => error === failure);
  } finally {
    globalThis.ClipboardItem = restore.ClipboardItem;
    Object.defineProperty(globalThis, 'isSecureContext', { value: restore.isSecureContext, configurable: true });
    Object.defineProperty(globalThis, 'navigator', { value: restore.navigator, configurable: true });
  }
});
