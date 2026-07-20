import assert from 'node:assert/strict';
import { mkdir, mkdtemp, readFile, writeFile } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { pathToFileURL } from 'node:url';

const repoRoot = process.cwd();
const tempDir = await mkdtemp(join(tmpdir(), 'gbdraw-orthogroup-stable-identity-'));
await writeFile(join(tempDir, 'package.json'), '{"type":"module"}\n', 'utf8');
await mkdir(join(tempDir, 'app'), { recursive: true });
await mkdir(join(tempDir, 'services'), { recursive: true });

const copyModule = async (sourceRelative, targetRelative) => {
  await writeFile(
    join(tempDir, targetRelative),
    await readFile(join(repoRoot, sourceRelative), 'utf8'),
    'utf8'
  );
};

await copyModule('gbdraw/web/js/app/feature-utils.js', 'app/feature-utils.js');
await copyModule('gbdraw/web/js/app/feature-sequence-fasta.js', 'app/feature-sequence-fasta.js');
await copyModule('gbdraw/web/js/app/losat-normalization.js', 'app/losat-normalization.js');
await copyModule('gbdraw/web/js/services/text-download.js', 'services/text-download.js');
await writeFile(join(tempDir, 'app', 'feature-editor-svg-actions.js'), `
  export const FEATURE_SELECTOR = '[data-gbdraw-feature-id]';
  export const getFeatureIdentity = (element) => String(
    element?.getAttribute?.('data-gbdraw-feature-id') || ''
  ).replace(/__part\\d+$/, '');
  export const getFeatureElements = (svg, featureId) => (
    svg?.featureElements?.get?.(String(featureId || '')) || []
  );
`, 'utf8');

const orthogroupSource = (await readFile(
  join(repoRoot, 'gbdraw/web/js/app/orthogroups.js'),
  'utf8'
)).replace(
  "from './feature-editor/svg-actions.js';",
  "from './feature-editor-svg-actions.js';"
);
await writeFile(join(tempDir, 'app', 'orthogroups.js'), orthogroupSource, 'utf8');

globalThis.window = {
  Vue: {
    computed: (getter) => ({ get value() { return getter(); } })
  }
};
let copiedText = '';
globalThis.navigator = {
  clipboard: {
    writeText: async (value) => { copiedText = String(value); }
  }
};

const { createOrthogroupEditor } = await import(
  pathToFileURL(join(tempDir, 'app', 'orthogroups.js'))
);

const ref = (value) => ({ value });
const visibleElement = {
  attrs: new Map([
    ['data-gbdraw-feature-id', 'stable-x_record_1'],
    ['data-gbdraw-stable-feature-id', 'stable-x'],
    ['data-gbdraw-record-index', '0'],
    ['stroke', '#111111'],
    ['stroke-width', '0.5']
  ]),
  getAttribute(name) { return this.attrs.get(name) ?? null; },
  hasAttribute(name) { return this.attrs.has(name); },
  setAttribute(name, value) { this.attrs.set(name, String(value)); },
  removeAttribute(name) { this.attrs.delete(name); }
};
const svg = {
  featureElements: new Map([['stable-x_record_1', [visibleElement]]]),
  querySelectorAll(selector) {
    if (selector === '[data-gbdraw-feature-id]') return [visibleElement];
    if (selector.includes('[data-og-original-stroke]')) {
      return visibleElement.hasAttribute('data-og-original-stroke') ? [visibleElement] : [];
    }
    return [];
  }
};

const group = {
  id: 'og_1',
  members: [
    { recordIndex: 0, recordId: 'record-a', featureSvgId: 'stable-x', proteinId: 'visible' },
    { recordIndex: 1, recordId: 'record-b', featureSvgId: 'stable-x', proteinId: 'hidden' }
  ]
};
const state = {
  orthogroups: ref([group]),
  orthogroupNameOverrides: {},
  orthogroupDescriptionOverrides: {},
  selectedOrthogroupId: ref('og_1'),
  orthogroupSearch: ref(''),
  orthogroupSortMode: ref('id'),
  selectedOrthogroupAlignmentFeature: ref(''),
  svgContainer: ref({ querySelector: () => svg }),
  showRightDrawer: ref(false),
  rightDrawerTab: ref('orthogroups'),
  showFeaturePanel: ref(false),
  showLegendPanel: ref(false),
  linearSeqs: [{ definition: 'record-a' }, { definition: 'record-b' }],
  extractedFeatures: ref([
    {
      fileIdx: 0,
      svg_id: 'stable-x_record_1',
      stable_svg_id: 'stable-x',
      rendered_svg_id: 'stable-x_record_1'
    }
  ]),
  biologicalFeatures: ref([
    { fileIdx: 0, svg_id: 'stable-x', stable_svg_id: 'stable-x', nucleotide_sequence: 'AAAA' },
    { fileIdx: 1, svg_id: 'stable-x', stable_svg_id: 'stable-x', nucleotide_sequence: 'CCCC' }
  ])
};

const editor = createOrthogroupEditor({ state, runAnalysis: null });
const members = editor.getEnrichedOrthogroupMembers(group);
assert.deepEqual(members.map((member) => member.nucleotideSequence), ['AAAA', 'CCCC']);
assert.equal(editor.getOrthogroupSequenceCount(group, 'nt'), 2);
await editor.copyOrthogroupSequences(group, 'nt');
assert.match(copiedText, /AAAA/);
assert.match(copiedText, /CCCC/);

editor.highlightOrthogroupById('og_1');
assert.equal(visibleElement.getAttribute('stroke'), '#2563eb');
assert.equal(visibleElement.getAttribute('stroke-width'), '2.4');
assert.equal(svg.featureElements.has('stable-x_record_2'), false);

state.orthogroups.value = [{
  id: 'og_hidden_only',
  members: [
    { recordIndex: 1, recordId: 'record-b', featureSvgId: 'stable-x', proteinId: 'hidden' }
  ]
}];
editor.highlightOrthogroupById('og_hidden_only');
assert.equal(visibleElement.getAttribute('stroke'), '#111111');
assert.equal(visibleElement.getAttribute('stroke-width'), '0.5');

console.log('orthogroup stable identity tests passed');
