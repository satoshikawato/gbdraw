// OV-671 and OV-672 (BA-D26): the side a loaded Result's Legend is laid out
// with is read through the positions of the Result's mode. A side outside
// them reads as the mode's default legend position, listed in the Load dialog
// (BA-D25, PLAN 12.7 (4)); case and outer spaces read as the position. A Result saved without composition metadata (Session
// 27 to 30) whose Session stored no side keeps `none`.
import { installSessionImportWorker } from './helpers/session-import-node.mjs';
import { evaluatePython, validateWebConfigOverrides } from './helpers/python-evaluator.mjs';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
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
setUnmanagedConfigOverrideValidator(async (payload) => ({ result: await validateWebConfigOverrides(payload) }));
const { transformLegacyResultSvg } = await import('../../gbdraw/web/js/app/app-setup.js');
const {
  applyCompositionEdit,
  COMPOSITION_METADATA_ATTRIBUTE
} = await import('../../gbdraw/web/js/app/legend-layout/composition-actions.js');
const { state } = await import('../../gbdraw/web/js/state.js');

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
  /** @type {string[][]} */
  const dialogs = [];
  const outcome = await importSession({
    target: { files: [new Blob([JSON.stringify(session)], { type: 'application/json' })], value: 'selected' }
  }, { transformLegacyResultSvg: port, confirmLoadRepairs: async ({ items }) => { dialogs.push(items); return true; } });
  // The Legend rows: a Session 30 also lists its feature recovery, which
  // needs the diagram Worker that this node test does not run (B50).
  const legendRows = dialogs.map((items) => items.filter((line) => line.startsWith('Legend'))).filter((items) => items.length > 0);
  return { outcome, transformed, dialogs: legendRows };
};

// A valid, a cased and spaced, an empty, and an outside value in each mode.
// A value that is not text is outside too; a falsy one the old `||` chain
// read as `none` (null, false, 0) and blank text are no stored side.
const CASES = {
  linear: [['top', 'top'], [' Top ', 'top'], ['', 'none'], ['  ', 'none'], ['upper_left', 'bottom'],
    ['banana', 'bottom'], [5, 'bottom'], [false, 'none']],
  circular: [['upper_left', 'upper_left'], [' Upper_Left ', 'upper_left'], ['', 'none'], ['banana', 'left'],
    [true, 'left'], [null, 'none'], [0, 'none']]
};
// BA-D26, BA-D25: the Load dialog lists a stored side outside the mode, once,
// before the Load; case, outer spaces and an unset side are no row.
const sideRows = (stored, expected, path) => (
  stored && !(typeof stored === 'string' && stored.trim() === '') && expected !== String(stored).trim().toLowerCase()
    ? [[`Legend of a saved Result: the side "${stored}" (${path}) is not one of its mode's positions; the Result is read with ${expected}.`]]
    : []
);

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
      const { outcome, transformed, dialogs } = await load(session);
      assert.equal(outcome?.status, 'ok', JSON.stringify(outcome?.error));
      assert.equal(transformed.length, 1);
      assert.equal(transformed[0].side, expected);
      assert.deepEqual(dialogs, sideRows(stored, expected, 'config.form.legend'));
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
      const { outcome, transformed, dialogs } = await load(session);
      assert.equal(outcome?.status, 'ok', JSON.stringify(outcome?.error));
      assert.deepEqual(dialogs, sideRows(stored, expected, 'results'));
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
  const session = await evaluatePython('from tests.utils.two_mode_session import two_mode_session\nresult = two_mode_session()');
  assert.equal(session.otherModeResult.renderRequest.mode, 'linear');
  session.otherModeResult.results[0].content = withStoredSide(session.otherModeResult.results[0].content, 'upper_left');
  const { outcome, transformed, dialogs } = await load(session);
  assert.equal(outcome?.status, 'ok', JSON.stringify(outcome?.error));
  assert.deepEqual(transformed.map(({ side, data }) => [data.mode, data.composition, side]), [['linear', null, 'bottom']]);
  assert.deepEqual(dialogs, sideRows('upper_left', 'bottom', 'otherModeResult.results'));
});

// OV-677 (BA-D31, BA-OD118): release 0.13.0 kept its Result across a mode
// switch and saved the shown mode in `ui.mode`; the Result's mode is the last
// Generate's (`cliInvocation.mode`). Load reads the Result as that mode, with
// the side it was drawn with (`ui.generatedLegendPosition`), and opens on it;
// the flat draft's settings stay with the shown mode.
const MODE_SWITCH = [
  // [fixture, Result mode, shown mode, drawn side, the shown mode's form legend]
  ['mode-switch-circular-result-linear-shown.v30.gbdraw-session.json.gz', 'circular', 'linear', 'left', 'bottom'],
  ['mode-switch-linear-result-circular-shown.v30.gbdraw-session.json.gz', 'linear', 'circular', 'top', 'left']
];
for (const [name, resultMode, shownMode, drawnSide, shownLegend] of MODE_SWITCH) {
  test(`OV-677: a 0.13.0 Session saved with ${shownMode} shown opens its ${resultMode} Result as ${resultMode}, laid out with ${drawnSide}`, async () => {
    const session = fixture(name);
    assert.deepEqual([session.version, session.ui.mode, session.cliInvocation.mode], [30, shownMode, resultMode]);
    const { outcome, transformed, dialogs } = await load(session);
    assert.equal(outcome?.status, 'ok', JSON.stringify(outcome?.error));
    assert.deepEqual(dialogs, []);
    assert.deepEqual(transformed.map(({ side, data }) => [data.mode, side]), [[resultMode, drawnSide]]);
    assert.deepEqual([state.mode.value, state.generatedMode.value, state.generatedLegendPosition.value], [resultMode, resultMode, drawnSide]);
    assert.equal(state.results.value.length, 1);
    // The flat draft's legend position is the shown mode's.
    assert.equal(state.drawings[shownMode].form.legend, shownLegend);
    // The first Legend edit keeps the Legend on the side 0.13.0 drew it.
    const { svg } = transformed[0];
    applyCompositionEdit(svg, { legendSide: storedSide(svg), legendLocalBox: null });
    assert.equal(JSON.parse(svg.getAttribute(COMPOSITION_METADATA_ATTRIBUTE)).legendSide, drawnSide);
  });
}

// The accepted residual of PD-OI-096: 0.13.0 wrote no `cliInvocation` for a
// Generate it could not replay (a Linear LOSAT run among them), so such a
// Session's Result still opens as the shown mode.
test('OV-677: a 0.13.0 Session saved after a mode switch without the last Generate\'s mode opens its Result as the shown mode', async () => {
  const session = fixture(MODE_SWITCH[0][0]);
  delete session.cliInvocation;
  const { outcome, transformed } = await load(session);
  assert.equal(outcome?.status, 'ok', JSON.stringify(outcome?.error));
  assert.deepEqual(transformed.map(({ data }) => data.mode), ['linear']);
  assert.equal(state.mode.value, 'linear');
});

// A blank drawn side is no stored side, read as `none` (PD-OI-096 〔値がないとき〕),
// as the same-mode path reads `config.form.legend`.
test('OV-677: a blank drawn side of a Result saved after a mode switch reads as none', async () => {
  const session = fixture(MODE_SWITCH[0][0]);
  session.ui.generatedLegendPosition = '  ';
  const { outcome, transformed, dialogs } = await load(session);
  assert.equal(outcome?.status, 'ok', JSON.stringify(outcome?.error));
  assert.deepEqual(dialogs, []);
  assert.deepEqual(transformed.map(({ side, data }) => [data.mode, side]), [['circular', 'none']]);
});

// OV-969 (BA-OD121): the title of such a Result was drawn with the Result's
// mode's title side, which 0.13.0 saved in `ui` (Circular:
// `generatedCircularPlotTitlePosition`, Linear: `linearPlotTitlePosition`);
// the stored form's side is the shown mode's. The fixtures draw no title, so
// each case adds one where 0.13.0 drew it, and the Circular case names `top`.
const TITLED_MODE_SWITCH = [
  // [fixture, Result mode, drawn title side, the shown mode's form title side, the title group]
  ['mode-switch-circular-result-linear-shown.v30.gbdraw-session.json.gz', 'circular', 'top', 'bottom',
    '<g id="plot_title" transform="translate(608,60)"><text x="0" y="0" font-size="32" text-anchor="middle">REGEX</text></g>'],
  ['mode-switch-linear-result-circular-shown.v30.gbdraw-session.json.gz', 'linear', 'bottom', 'none',
    '<g id="plot_title" transform="translate(1168,370)"><text x="0" y="0" font-size="32" text-anchor="middle">REGEX</text></g>']
];
for (const [name, resultMode, drawnTitle, shownTitle, titleGroup] of TITLED_MODE_SWITCH) {
  test(`OV-969: a 0.13.0 ${resultMode} Result saved after a mode switch is laid out with the title side ${drawnTitle} it was drawn with`, async () => {
    const session = fixture(name);
    session.results[0].content = session.results[0].content.replace(/<\/svg>\s*$/, `${titleGroup}</svg>`);
    if (resultMode === 'circular') session.ui.generatedCircularPlotTitlePosition = drawnTitle;
    assert.equal(resultMode === 'circular' ? session.ui.generatedCircularPlotTitlePosition : session.ui.linearPlotTitlePosition, drawnTitle);
    assert.equal(session.config.adv.plot_title_position, shownTitle);
    const { outcome, transformed } = await load(session);
    assert.equal(outcome?.status, 'ok', JSON.stringify(outcome?.error));
    assert.deepEqual(transformed.map(({ data }) => [data.mode, data.composition?.titleSide]), [[resultMode, drawnTitle]]);
    const { svg } = transformed[0];
    applyCompositionEdit(svg, { legendSide: storedSide(svg), legendLocalBox: null });
    assert.equal(JSON.parse(svg.getAttribute(COMPOSITION_METADATA_ATTRIBUTE)).titleSide, drawnTitle);
  });
}

