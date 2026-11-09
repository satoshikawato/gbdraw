// Text crosses between Python and JavaScript intact whatever its address in
// the Pyodide heap (OV-2985). After one large Generate the Worker's heap stays
// above 2 GiB, and Pyodide 0.29 read and wrote a string above 2^31 with a
// character beyond U+00FF at a wrong address: Python-to-JavaScript text came
// back as NUL characters or failed to convert, and JavaScript-to-Python text
// arrived as NUL characters. Each case runs a crossing the Worker uses, with
// every string above 2 GiB: a helper call (its JSON arguments and its JSON
// reply, built with ensure_ascii=False) and the Generate result's conversion.
import assert from 'node:assert/strict';
import { after, test } from 'node:test';
import { HIGH_ADDRESS, allocateAboveTwoGiB, loadVendoredPyodide } from './helpers/pyodide-high-heap.mjs';

globalThis.self ??= {};
const { callJsonHelper, convertGenerationResult } = await import('../../gbdraw/web/js/workers/diagram-generation-worker.js');

// Latin-1, two-byte and four-byte characters: Python stores each kind differently.
const TEXTS = Object.freeze({
  latin1: 'Homo sapiens µ é',
  ucs2: 'ヒト 1 番染色体 – µ',
  ucs4: 'ヒト 1 番染色体 – µ \u{1F9EC}'
});
// A Python string literal of `text` written in ASCII, so that the source itself crosses intact.
const pythonLiteral = (text) => `"${Array.from(text, (character) => {
  const point = character.codePointAt(0);
  return point < 128 ? character : `\\U${point.toString(16).padStart(8, '0')}`;
}).join('')}"`;
// Long enough that Python allocates each string in its own block above 2 GiB.
const REPEAT = 64;

const pyodide = await loadVendoredPyodide();
pyodide.runPython(`
import json

ADDRESSES = []

def call_web_json_helper(helper_name, request_json):
    # The shape of a helper in gbdraw/web/js/app/python-helpers.js: JSON in, JSON out.
    ADDRESSES.append(id(request_json))
    request = json.loads(request_json)
    reply = json.dumps({"helper": helper_name, "title": request["title"], "kind": request["kind"]}, ensure_ascii=False)
    ADDRESSES.append(id(reply))
    return reply

def generate_result(text):
    # The strings of a Generate result that travel by toJs, not by a reply buffer.
    head = '<svg xmlns="http://www.w3.org/2000/svg"><title>' + text + "</title>"
    raw = '<text class="definition">' + text + "</text>"
    start = '<g id="record-' + text + '">'
    name = "record " + text
    message = "The record " + text + " cannot be drawn."
    label_fit = json.dumps({"rows": [{"label": text}]}, ensure_ascii=False)
    ADDRESSES.extend(id(value) for value in (head, raw, start, name, message, label_fit))
    return {
        "results": [{"name": name, "document": {"head": head, "tail": "</svg>", "pieces": [
            {"kind": "raw", "text": raw}, {"kind": "content", "start": start}
        ]}}],
        "error": {"message": message},
        "labelFit": label_fit,
    }
`);
const takeAddresses = () => {
  const addresses = pyodide.runPython('[ADDRESSES.copy(), ADDRESSES.clear()][0]').toJs();
  assert.ok(addresses.length > 0);
  assert.deepEqual(addresses.filter((address) => address < HIGH_ADDRESS), [], 'a string sat below 2 GiB');
};

const release = allocateAboveTwoGiB(pyodide);
after(release);

for (const [kind, text] of Object.entries(TEXTS)) {
  const title = text.repeat(REPEAT);

  test(`a helper call carries ${kind} text both ways above 2 GiB`, () => {
    let reply;
    try {
      reply = callJsonHelper(pyodide, 'validate_config_overrides_json', [JSON.stringify({ title, kind })]);
    } catch (error) {
      assert.fail(`the helper call failed: ${error?.code ?? error}`);
    }
    takeAddresses();
    assert.deepEqual(reply, { helper: 'validate_config_overrides_json', title, kind });
  });

  test(`the Generate result carries ${kind} text above 2 GiB`, () => {
    const make = pyodide.runPython(`lambda: generate_result(${pythonLiteral(title)})`);
    const resultHandle = make();
    make.destroy();
    let result;
    try {
      result = convertGenerationResult(resultHandle);
    } finally {
      resultHandle.destroy();
    }
    takeAddresses();
    const [drawn] = result.results;
    assert.equal(drawn.name, `record ${title}`, 'Result name');
    assert.equal(drawn.document.head, `<svg xmlns="http://www.w3.org/2000/svg"><title>${title}</title>`, 'SVG head');
    assert.equal(drawn.document.pieces[0].text, `<text class="definition">${title}</text>`, 'raw piece');
    assert.equal(drawn.document.pieces[1].start, `<g id="record-${title}">`, 'start tag');
    assert.equal(result.error.message, `The record ${title} cannot be drawn.`, 'error payload');
    assert.deepEqual(JSON.parse(result.labelFit), { rows: [{ label: title }] }, 'label fit payload');
  });
}
