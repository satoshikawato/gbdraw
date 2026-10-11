// FL-10: the typed Output Prefix reaches Python as typed, which rejects a value
// that is not one portable file name and names the field; a name derived from
// a record ID follows the browser's download-name rule, so a Result's name and
// its downloaded file's name are the same.
import assert from 'node:assert/strict';
import { execFileSync } from 'node:child_process';
import { fileURLToPath } from 'node:url';
import { test } from 'node:test';

globalThis.window = {
  Vue: {
    ref: (value) => ({ value }), reactive: (value) => value,
    computed: (getter) => ({ get value() { return getter(); } }),
    nextTick: async () => {}
  },
  DOMPurify: { sanitize: (value) => value }
};

const { state, createDefaultAdv, createDefaultForm } = await import('../../gbdraw/web/js/state.js');
const { buildCanonicalRenderRequest } = await import('../../gbdraw/web/js/services/session-request.js');
const { downloadSafeName } = await import('../../gbdraw/web/js/utils/download-names.js');

const REPO = fileURLToPath(new URL('../../', import.meta.url));
const genbank = (recordId) => {
  const text = `LOCUS       FL10                       40 bp    DNA     linear   UNK 01-JAN-1980
DEFINITION  Output prefix fixture.
ACCESSION   FL10
VERSION     ${recordId}
KEYWORDS    .
SOURCE      .
  ORGANISM  .
            .
FEATURES             Location/Qualifiers
     CDS             1..30
                     /product="test protein"
ORIGIN
        1 atgcatgcat gcatgcatgc atgcatgcat gcatgcatgc
//
`;
  return { name: 'fl10.gb', type: 'text/plain', size: text.length, lastModified: 0, data: btoa(text) };
};
// The output prefix of a Circular request: the typed prefix, or the one record's ID.
const circularPrefix = ({ prefix = '', recordId = 'FL10.1' } = {}) => {
  const drawing = state.drawings.circular;
  Object.assign(drawing.form, createDefaultForm(), { prefix });
  Object.assign(drawing.adv, createDefaultAdv('circular'));
  state.mode.value = 'circular';
  state.cInputType.value = 'gb';
  state.circularRecordList.value = [{ record_id: recordId, record_length: 40 }];
  const canonical = buildCanonicalRenderRequest({
    state, drawing, filesData: { c_gb: genbank(recordId), linearSeqs: [] }, comparisonPlanSnapshot: null
  });
  return canonical.renderRequest.output.prefix;
};
// Python's check of an output prefix: null when accepted, else its diagnostic.
const pythonPrefixDiagnostic = (prefix) => JSON.parse(execFileSync('python', ['-c', `
import json, sys
from gbdraw.api.requests import RenderOutputRequest
from gbdraw.exceptions import ValidationError
try:
    RenderOutputRequest(output_prefix=sys.argv[1])
    print('null')
except ValidationError as error:
    print(json.dumps(getattr(error, 'diagnostic', None) or {'message': str(error)}))
`, prefix], { encoding: 'utf8', cwd: REPO, env: { ...process.env, PYTHONPATH: REPO } }));

test('a typed Output Prefix reaches Python as typed, which names the field when it rejects it', () => {
  assert.equal(circularPrefix({ prefix: '../../x' }), '../../x');
  assert.equal(circularPrefix({ prefix: '  report  ' }), 'report');
  assert.deepEqual(pythonPrefixDiagnostic('../../x'),
    { code: 'INPUT_INVALID', field: 'output_prefix', reason: 'FILENAME' });
  assert.equal(pythonPrefixDiagnostic('report'), null);
});

test('a name derived from a record ID follows the download-name rule and Python accepts it', () => {
  const cases = {
    'gi|1|ref|X': 'gi_1_ref_X', '../../x': '_.._x', '.hidden': 'hidden', 'a b.': 'a b',
    CON: '_CON', 'con.txt': '_con.txt', 'NC_001416.1': 'NC_001416.1', '#1': '#1', 'a\tb': 'a_b'
  };
  for (const [value, expected] of Object.entries(cases)) {
    assert.equal(downloadSafeName(value, 'out'), expected, value);
    assert.equal(pythonPrefixDiagnostic(expected), null, `Python accepts ${expected}`);
  }
  assert.equal(downloadSafeName('...', 'out'), 'out');
  assert.equal(circularPrefix({ recordId: 'gi|1|ref|X' }), 'gi_1_ref_X');
});
