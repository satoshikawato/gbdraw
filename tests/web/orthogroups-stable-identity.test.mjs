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
Object.defineProperty(globalThis, 'navigator', {
  configurable: true,
  value: {
    clipboard: {
      writeText: async (value) => { copiedText = String(value); }
    }
  }
});

const { createOrthogroupEditor } = await import(
  pathToFileURL(join(tempDir, 'app', 'orthogroups.js'))
);
const {
  isInternalProteinDisplayId,
  resolveDisplayProteinId
} = await import(pathToFileURL(join(tempDir, 'app', 'feature-utils.js')));

assert.equal(isInternalProteinDisplayId(`h_${'c'.repeat(26)}`), true);
assert.equal(isInternalProteinDisplayId(`f_${'d'.repeat(64)}`), true);
assert.equal(
  resolveDisplayProteinId(null, {
    proteinId: `h_${'c'.repeat(26)}`,
    label: 'display fallback'
  }),
  'display fallback'
);
assert.equal(
  resolveDisplayProteinId(null, { proteinId: `h_${'c'.repeat(26)}` }),
  ''
);
const runtimeHandle = `h_${'e'.repeat(26)}`;
const featureAnalysisId = `f_${'f'.repeat(64)}`;
const v35TransportId = `record@instance|alias~f_${'a'.repeat(64)}`;
assert.equal(
  resolveDisplayProteinId(
    {
      displayProteinId: runtimeHandle,
      sourceProteinId: featureAnalysisId,
      qualifiers: { protein_id: [v35TransportId] },
      locusTag: 'WP_SAFE.1'
    },
    { label: runtimeHandle }
  ),
  'WP_SAFE.1'
);
assert.equal(
  resolveDisplayProteinId(
    {
      displayProteinId: runtimeHandle,
      sourceProteinId: featureAnalysisId,
      qualifiers: { protein_id: [v35TransportId] }
    },
    { label: runtimeHandle }
  ),
  ''
);
assert.equal(
  resolveDisplayProteinId({
    qualifiers: { protein_id: [runtimeHandle, 'WP_ARRAY_FALLBACK.1'] }
  }),
  'WP_ARRAY_FALLBACK.1'
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
    {
      recordIndex: 0,
      recordId: 'record-a',
      featureSvgId: 'stable-x',
      proteinId: `h_${'a'.repeat(26)}`
    },
    {
      recordIndex: 1,
      recordId: 'record-b',
      featureSvgId: 'stable-x',
      proteinId: `h_${'b'.repeat(26)}`
    }
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
    {
      fileIdx: 0,
      record_id: 'record-a',
      svg_id: 'stable-x',
      stable_svg_id: 'stable-x',
      locus_tag: 'LOC_A',
      nucleotide_sequence: 'AAAA',
      amino_acid_sequence: 'MK'
    },
    {
      fileIdx: 1,
      record_id: 'record-b',
      svg_id: 'stable-x',
      stable_svg_id: 'stable-x',
      locus_tag: 'LOC_B',
      nucleotide_sequence: 'CCCC',
      amino_acid_sequence: 'MP'
    }
  ])
};

const editor = createOrthogroupEditor({ state, runAnalysis: null });
const members = editor.getEnrichedOrthogroupMembers(group);
assert.deepEqual(members.map((member) => member.nucleotideSequence), ['AAAA', 'CCCC']);
assert.deepEqual(members.map((member) => member.displayProteinId), ['LOC_A', 'LOC_B']);
assert.equal(editor.getOrthogroupSequenceCount(group, 'nt'), 2);
await editor.copyOrthogroupSequences(group, 'nt');
assert.match(copiedText, /AAAA/);
assert.match(copiedText, /CCCC/);
await editor.copyOrthogroupSequences(group, 'aa');
assert.match(copiedText, />LOC_A\b/);
assert.match(copiedText, />LOC_B\b/);
assert.doesNotMatch(copiedText, /h_[a-z2-7]{26}/);

const membersWithoutStableIds = editor.getEnrichedOrthogroupMembers({
  id: 'og_no_stable_ids',
  members: [
    {
      recordId: 'record-c',
      proteinId: runtimeHandle,
      sourceProteinId: 'WP_VISIBLE.1',
      label: featureAnalysisId
    },
    {
      recordId: 'record-d',
      proteinId: runtimeHandle,
      displayProteinId: featureAnalysisId,
      sourceProteinId: v35TransportId,
      label: runtimeHandle
    }
  ]
});
assert.deepEqual(
  membersWithoutStableIds.map((member) => member.displayProteinId),
  ['WP_VISIBLE.1', '']
);

const downloadedFilenames = [];
globalThis.document = {
  createElement: () => ({
    addEventListener() {},
    click() { downloadedFilenames.push(this.download); }
  }),
  body: {
    appendChild() {},
    removeChild() {}
  }
};
Object.defineProperty(globalThis.URL, 'createObjectURL', {
  configurable: true,
  value: () => 'blob:orthogroup-test'
});
Object.defineProperty(globalThis.URL, 'revokeObjectURL', {
  configurable: true,
  value: () => {}
});
editor.downloadOrthogroupMemberSequence(
  {
    proteinId: runtimeHandle,
    sourceProteinId: 'WP_DOWNLOAD.1',
    aminoAcidSequence: 'MK'
  },
  'aa',
  { id: 'og_download' }
);
editor.downloadOrthogroupMemberSequence(
  {
    proteinId: runtimeHandle,
    stableFeatureSvgId: featureAnalysisId,
    aminoAcidSequence: 'MK'
  },
  'aa',
  { id: 'og_download' }
);
assert.deepEqual(downloadedFilenames, [
  'og_download_WP_DOWNLOAD.1_aa.faa',
  'og_download_member_aa.faa'
]);
assert.doesNotMatch(downloadedFilenames.join('\n'), /(?:h_[a-z2-7]{26}|f_[0-9a-f]{64})/);

const indexSource = await readFile(join(repoRoot, 'gbdraw/web/index.html'), 'utf8');
assert.doesNotMatch(
  indexSource,
  /member\.displayProteinId\s*\|\|\s*member\.(?:sourceProteinId|locusTag|label)/
);
assert.doesNotMatch(
  indexSource,
  /currentMember\.displayProteinId\s*\|\|\s*clickedOrthogroupDetail\.currentMember\.sourceProteinId/
);

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
