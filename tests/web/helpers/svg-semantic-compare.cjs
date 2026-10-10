// "Generate equals the loaded preview": compares two SVG files with
// tests/utils/svg_compare.compare_svgs, ignoring the label-binding and
// record-placement attributes that a saved preview and a fresh Generate may
// number differently. Gallery publication parity
// (tests/web/contracts/gallery-publication-parity.serial.spec.js) and the
// promotion user journeys (tools/audit/user-journey.audit.spec.js) share it.
const { spawnSync } = require('node:child_process');
const { resolve } = require('node:path');

const IGNORED_ATTRIBUTES = [
  'data-label-feature-id',
  'data-gbdraw-label-binding-schema',
  'data-record-key',
  'data-record-translation-x',
  'data-record-translation-y'
];

const COMPARE_COMMAND = [
  'import json, sys',
  'from tests.utils.svg_compare import compare_svgs',
  'result = compare_svgs(sys.argv[1], sys.argv[2], ignored_attributes=set(json.loads(sys.argv[3])))',
  'print(result.message)',
  'print("\\n".join(result.differences))',
  'raise SystemExit(0 if result.equal else 1)'
].join(';');

// Exit status 0 means equal; `report` holds the comparison message and differences.
const compareSvgFiles = (left, right, {
  repoRoot = resolve(process.env.GBDRAW_REPO || process.cwd()),
  python = process.env.GBDRAW_PYTHON || 'python'
} = {}) => {
  const run = spawnSync(python, ['-c', COMPARE_COMMAND, left, right, JSON.stringify(IGNORED_ATTRIBUTES)], {
    cwd: repoRoot,
    encoding: 'utf8'
  });
  return { status: run.status, equal: run.status === 0, report: `${run.stdout || ''}\n${run.stderr || ''}`.trim() };
};

module.exports = { compareSvgFiles };
