// OV-671 and OV-672 (BA-D26): the side a loaded Result's Legend is laid out
// with is read through the positions of the Result's mode. A side outside
// them reads as the mode's default legend position; case and outer spaces
// read as the position. A Result saved without composition metadata (Session
// 27 to 30) whose Session stored no side keeps `none`.
import { installSessionImportWorker } from './helpers/session-import-node.mjs';
import assert from 'node:assert/strict';
import { execFileSync } from 'node:child_process';
import { readFileSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import test from 'node:test';
import { gunzipSync } from 'node:zlib';
import { installFakeSvgDom } from './fake-svg-dom.mjs';

globalThis.window = {
  Vue: {
    ref: (value) => ({ value }),
    reactive: (value) => value,
    computed: (getter) => ({ get value() { return getter(); } }),
    nextTick: async () => {}
  },
  DOMPurify: { sanitize: (value) => value }
};
globalThis.document = {};
installFakeSvgDom();
globalThis.alert = () => {};
installSessionImportWorker();

const { importSession, setUnmanagedConfigOverrideValidator } = await import('../../gbdraw/web/js/services/config.js');
// The Worker's typed override check, run by Python (as session-cli-compatibility does).
setUnmanagedConfigOverrideValidator((payload) => ({ result: JSON.parse(execFileSync('python', ['-c', `
import json, sys
from gbdraw.web_support.config_overrides import validate_web_config_overrides_json
p = json.load(sys.stdin)
print(validate_web_config_overrides_json(p['mode'], json.dumps(p['config']),
    json.dumps(p['configOverrides']), json.dumps(p['managedPaths']), p['requireUnmanagedOnly']))
`], { input: JSON.stringify(payload), encoding: 'utf8', cwd: fileURLToPath(new URL('../../', import.meta.url)) })) }));
const { transformLegacyResultSvg } = await import('../../gbdraw/web/js/app/app-setup.js');
const {
  applyCompositionEdit,
  COMPOSITION_METADATA_ATTRIBUTE
} = await import('../../gbdraw/web/js/app/legend-layout/composition-actions.js');

const REPO = new URL('../../', import.meta.url);
const fixture = (name) => {
  const bytes = readFileSync(new URL(`tests/fixtures/sessions/${name}`, REPO));
  return JSON.parse((bytes[0] === 0x1f && bytes[1] === 0x8b ? gunzipSync(bytes) : bytes).toString('utf8'));
};

/** @param {{ getAttribute: (name: string) => string | null }} svg */
const storedSide = (svg) => {
  const raw = svg.getAttribute(COMPOSITION_METADATA_ATTRIBUTE);
  return raw === null ? null : JSON.parse(raw).legendSide;
};

// Load through the composition root's port, recording each Result it transforms.
const load = async (session) => {
  /** @type {{ svg: any, side: unknown, data: any }[]} */
  const transformed = [];
  const port = Object.assign((svg, data) => {
    const changed = transformLegacyResultSvg(svg, data);
    transformed.push({ svg, side: storedSide(svg), data });
    return changed;
  }, { appliesToContent: (...args) => transformLegacyResultSvg.appliesToContent(...args) });
  const outcome = await importSession({
    target: { files: [new Blob([JSON.stringify(session)], { type: 'application/json' })], value: 'selected' }
  }, { transformLegacyResultSvg: port });
  return { outcome, transformed };
};

// A valid, a cased and spaced, an empty, and an outside value in each mode.
// A value that is not text is outside too; one the old `||` chain skipped
// (null, false) is no stored side.
const CASES = {
  linear: [['top', 'top'], [' Top ', 'top'], ['', 'none'], ['upper_left', 'bottom'], ['banana', 'bottom'],
    [5, 'bottom'], [false, 'none']],
  circular: [['upper_left', 'upper_left'], [' Upper_Left ', 'upper_left'], ['', 'none'], ['banana', 'left'],
    [true, 'left'], [null, 'none']]
};
// Sessions 27 to 30: no stored request, a Result without composition metadata.
const LEGACY = {
  linear: 'scale-interval-negative-linear-cli.v30.gbdraw-session.json.gz',
  circular: 'scale-interval-zero-circular-cli.v30.gbdraw-session.json.gz'
};

for (const [mode, cases] of Object.entries(CASES)) {
  for (const [stored, expected] of cases) {
    test(`OV-671: a ${mode} Session 30 Result with legend ${JSON.stringify(stored)} is laid out with ${expected}`, async () => {
      const session = fixture(LEGACY[mode]);
      assert.equal(session.version, 30);
      assert.equal(session.ui.mode, mode);
      session.config.form.legend = stored;
      const { outcome, transformed } = await load(session);
      assert.equal(outcome?.status, 'ok', JSON.stringify(outcome?.error));
      assert.equal(transformed.length, 1);
      assert.equal(transformed[0].side, expected);
    });
  }
}

// Session 40+: the Result carries composition metadata (written by Python).
const CURRENT = {
  linear: 'se08-main-linear-cli.v42.gbdraw-session.json.gz',
  circular: 'scale-interval-zero-circular-cli.v44.gbdraw-session.json.gz'
};
const METADATA_CASES = {
  linear: [['right', 'right'], [' Top ', 'top'], ['', 'bottom'], ['upper_left', 'bottom']],
  circular: [['upper_left', 'upper_left'], [' Upper_Left ', 'upper_left'], ['', 'left'], ['banana', 'left']]
};
const decodeAttribute = (value) => value.replace(/&quot;/g, '"').replace(/&lt;/g, '<').replace(/&gt;/g, '>')
  .replace(/&amp;/g, '&');
// The metadata as a hand-edited file may hold it. The fake DOM keeps attribute
// text as written, so the test writes it with single quotes and no entities.
const withStoredSide = (content, side) => content.replace(
  /\sdata-gbdraw-composition="([^"]*)"/,
  (_match, value) => {
    const metadata = JSON.parse(decodeAttribute(value));
    metadata.legendSide = side;
    return ` data-gbdraw-composition='${JSON.stringify(metadata)}'`;
  }
);

for (const [mode, cases] of Object.entries(METADATA_CASES)) {
  for (const [stored, expected] of cases) {
    test(`OV-672: a ${mode} Result whose composition metadata names ${JSON.stringify(stored)} loads with ${expected}`, async () => {
      const session = fixture(CURRENT[mode]);
      assert.equal(session.renderRequest.mode, mode);
      session.results[0].content = withStoredSide(session.results[0].content, stored);
      const { outcome, transformed } = await load(session);
      assert.equal(outcome?.status, 'ok', JSON.stringify(outcome?.error));
      if (stored === expected) {
        // A side of the mode's positions is admitted as written, unparsed.
        assert.equal(transformed.length, 0);
      } else {
        assert.equal(transformed.length, 1);
        assert.equal(transformed[0].side, expected);
      }
    });
  }
}

test('OV-671: a 0.13.0 Linear Result drawn with upper_left keeps its SVG at Load, and its first Legend edit docks the Legend at the bottom', async () => {
  const session = fixture(LEGACY.linear);
  session.config.form.legend = 'upper_left';
  const { svg } = (await load(session)).transformed[0];
  const legend = svg.getElementById('legend');
  assert.ok(legend);
  // Load only prepends a zero translation: the Legend stays where 0.13.0 drew it.
  assert.match(legend.getAttribute('transform'), /^translate\(0,0\) /);
  assert.equal(legend.getAttribute('display'), null);
  // The relayout a Legend edit runs (`refreshLegendGeometry`) at the stored side.
  applyCompositionEdit(svg, { legendSide: storedSide(svg), legendLocalBox: null });
  const metadata = JSON.parse(svg.getAttribute(COMPOSITION_METADATA_ATTRIBUTE));
  const primary = metadata.primary.finalBounds;
  const legendTop = metadata.legend.automaticTranslation[1] + metadata.legend.localBounds.y;
  assert.equal(metadata.legendSide, 'bottom');
  assert.ok(legendTop >= primary.y + primary.height, `${legendTop} below ${JSON.stringify(primary)}`);
  assert.equal(legend.getAttribute('display'), null);
});

test('OV-672: the other mode\'s Result set reads a stored side outside its mode as the mode\'s default', async () => {
  // tests/utils/two_mode_session.py: a Circular Session 46 with a Linear `otherModeResult`.
  const session = JSON.parse(execFileSync('python', ['-c',
    'import json; from tests.utils.two_mode_session import two_mode_session; print(json.dumps(two_mode_session()))'
  ], { encoding: 'utf8', cwd: fileURLToPath(REPO), maxBuffer: 1 << 30 }));
  assert.equal(session.otherModeResult.renderRequest.mode, 'linear');
  session.otherModeResult.results[0].content = withStoredSide(session.otherModeResult.results[0].content, 'upper_left');
  const { outcome, transformed } = await load(session);
  assert.equal(outcome?.status, 'ok', JSON.stringify(outcome?.error));
  assert.deepEqual(transformed.map(({ side, data }) => [data.mode, data.composition, side]), [['linear', null, 'bottom']]);
});
