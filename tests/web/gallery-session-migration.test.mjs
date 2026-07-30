import assert from 'node:assert/strict';
import { cp, mkdtemp, readFile, writeFile } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { pathToFileURL } from 'node:url';
import { gunzipSync } from 'node:zlib';

const repoRoot = process.cwd();
const sourceRoot = join(repoRoot, 'gbdraw', 'web', 'js');
const tempRoot = await mkdtemp(join(tmpdir(), 'gbdraw-gallery-session-migration-'));
await cp(sourceRoot, join(tempRoot, 'js'), { recursive: true });
await writeFile(join(tempRoot, 'package.json'), '{"type":"module"}', 'utf8');

const { promoteGallerySessionToCurrent } = await import(
  pathToFileURL(join(tempRoot, 'js', 'services', 'gallery-session-migration.js'))
);
const { projectCanonicalSessionRequest } = await import(
  pathToFileURL(join(tempRoot, 'js', 'services', 'session-request.js'))
);

const sessionRoot = join(repoRoot, 'gbdraw', 'web', 'gallery', 'sessions');
const loadSession = async (name) => {
  const raw = await readFile(join(sessionRoot, name));
  const decoded = raw[0] === 0x1f && raw[1] === 0x8b ? gunzipSync(raw) : raw;
  return JSON.parse(decoded.toString('utf8'));
};

const resourceText = (session, ref) => Buffer.from(
  session.resources[ref.resourceId].data,
  'base64'
).toString('utf8');

const syntheticResource = (kind, name, text) => ({
  kind,
  name,
  type: 'text/plain',
  size: Buffer.byteLength(text),
  lastModified: 0,
  encoding: 'base64',
  data: Buffer.from(text).toString('base64')
});

const syntheticCliSession = {
  version: 33,
  config: {
    cliOptions: {
      rawArgs: ['--record_label', 'CLI label', '--reverse_complement', '1']
    }
  },
  renderRequest: {
    schema: 2,
    mode: 'linear',
    records: [{
      recordKey: 'record-1',
      source: { kind: 'genbank', resourceId: 'record-source' },
      selector: null,
      region: null,
      presentation: null
    }],
    diagramOptions: {
      config: { form: { legend: 'right' }, marker: 'must-survive' },
      output: {
        outputPrefix: 'cli',
        legend: 'right',
        plotTitlePosition: 'bottom'
      },
      featureShapes: {}
    },
    layout: {},
    comparisons: [],
    output: { prefix: 'cli', formats: ['interactive_svg'], overwrite: true }
  },
  resources: {
    'record-source': syntheticResource('genbank', 'cli.gb', 'LOCUS synthetic')
  }
};
const originalCliRequest = structuredClone(syntheticCliSession.renderRequest);
const promotedSyntheticCli = promoteGallerySessionToCurrent(syntheticCliSession);
assert.equal(promotedSyntheticCli.renderRequest.schema, 5);
assert.equal(promotedSyntheticCli.renderRequest.grouping, 'single');
assert.deepEqual(promotedSyntheticCli.renderRequest.diagramOptions.output, {
  legend: 'right',
  plotTitlePosition: 'bottom'
});
assert.equal(promotedSyntheticCli.renderRequest.records[0].presentation.label, 'CLI label');
assert.equal(promotedSyntheticCli.renderRequest.records[0].presentation.reverseComplement, false);
assert.equal(promotedSyntheticCli.renderRequest.diagramOptions.featureShapes.repeat_region, 'underlay');
assert.deepEqual(
  promotedSyntheticCli.renderRequest.diagramOptions.config,
  originalCliRequest.diagramOptions.config
);
assert.deepEqual(syntheticCliSession.renderRequest, originalCliRequest);

const syntheticGuiSession = {
  version: 33,
  config: {
    form: {
      prefix: 'gui',
      legend: 'left',
      labels_mode: 'out'
    },
    adv: {
      def_font_size: 31,
      evalue: 0.01,
      min_bitscore: 50,
      identity: 30,
      alignment_length: 0,
      depth_tick_interval: 15
    },
    losat: {
      blastp: {
        collinearMaxGeneGap: 4
      }
    },
    palette: 'pending-palette',
    colors: { CDS: '#111111' }
  },
  ui: {
    appliedPaletteName: 'applied-palette',
    appliedPaletteColors: { CDS: '#abcdef' }
  },
  renderRequest: {
    schema: 2,
    mode: 'circular',
    records: [{
      recordKey: 'record-1',
      source: { kind: 'genbank', resourceId: 'record-source' },
      selector: null,
      region: null,
      presentation: {
        label: null,
        subtitle: null,
        reverseComplement: false,
        gridRow: null,
        gridColumn: null
      }
    }],
    diagramOptions: {
      configOverrides: { show_labels: false, circular_definition_font_size: 10 },
      tracks: { circularTrackSlots: null, circularTrackAxisIndex: null },
      output: { outputPrefix: 'old', legend: 'right', plotTitlePosition: 'none' },
      selectedFeaturesSet: ['CDS'],
      featureShapes: { repeat_region: 'rectangle' },
      dinucleotide: 'GC'
    },
    layout: {},
    comparisons: [],
    output: { prefix: 'old', formats: ['interactive_svg'], overwrite: true }
  },
  resources: {
    'record-source': syntheticResource('genbank', 'gui.gb', 'LOCUS synthetic'),
    'colors-default-colors-file': syntheticResource(
      'colors-default-colors-file',
      'stale-colors.tsv',
      'CDS\t#000000\n'
    )
  }
};
const promotedSyntheticGui = promoteGallerySessionToCurrent(syntheticGuiSession);
const syntheticGuiOptions = promotedSyntheticGui.renderRequest.diagramOptions;
assert.equal(promotedSyntheticGui.renderRequest.schema, 5);
assert.equal(promotedSyntheticGui.renderRequest.grouping, 'single');
assert.equal(promotedSyntheticGui.renderRequest.output.prefix, 'gui');
assert.equal(promotedSyntheticGui.renderRequest.output.overwrite, false);
assert.deepEqual(syntheticGuiOptions.output, {
  legend: 'left',
  plotTitlePosition: 'none'
});
assert.equal(syntheticGuiOptions.output.legend, 'left');
assert.equal(syntheticGuiOptions.configOverrides['labels.circular.scope'], 'outer');
assert.equal(
  syntheticGuiOptions.configOverrides['objects.definition.circular.font_size'],
  31
);
assert.equal(
  syntheticGuiOptions.configOverrides['objects.depth.large_tick_interval'],
  15
);
assert.equal(syntheticGuiOptions.colors.defaultColorsPalette, 'applied-palette');
assert.match(
  resourceText(promotedSyntheticGui, syntheticGuiOptions.colors.defaultColorsFile),
  /CDS\t#abcdef/
);
assert.doesNotMatch(
  resourceText(promotedSyntheticGui, syntheticGuiOptions.colors.defaultColorsFile),
  /#000000/
);
const projectedSyntheticGui = projectCanonicalSessionRequest({
  renderRequest: promotedSyntheticGui.renderRequest,
  resources: promotedSyntheticGui.resources,
  webFiles: promotedSyntheticGui.webFiles
});
assert.equal(projectedSyntheticGui.config.form.legend, 'left');
assert.equal(projectedSyntheticGui.config.adv.def_font_size, 31);
assert.equal(projectedSyntheticGui.config.colors.CDS, '#abcdef');
assert.equal(
  promotedSyntheticGui.config.losat.blastp.collinearMaxUnitGap,
  4
);
assert.equal(
  Object.hasOwn(
    promotedSyntheticGui.config.losat.blastp,
    'collinearMaxGeneGap'
  ),
  false
);
assert.equal(
  Object.hasOwn(promotedSyntheticGui.config.adv, 'depth_tick_interval'),
  false
);

const staleCurrentConfig = structuredClone(promotedSyntheticGui);
staleCurrentConfig.config.losat.blastp.collinearMaxGeneGap =
  staleCurrentConfig.config.losat.blastp.collinearMaxUnitGap;
delete staleCurrentConfig.config.losat.blastp.collinearMaxUnitGap;
const repairedCurrentConfig = promoteGallerySessionToCurrent(staleCurrentConfig);
assert.equal(repairedCurrentConfig.config.losat.blastp.collinearMaxUnitGap, 4);
assert.equal(
  Object.hasOwn(repairedCurrentConfig.config.losat.blastp, 'collinearMaxGeneGap'),
  false
);
const currentWithoutWebConfig = structuredClone(promotedSyntheticGui);
delete currentWithoutWebConfig.config;
assert.equal(
  promoteGallerySessionToCurrent(currentWithoutWebConfig).config,
  undefined
);

const hmmt = await loadSession('HmmtDNA_ATskew.gbdraw-session.json');
const promotedHmmt = promoteGallerySessionToCurrent(hmmt);
const hmmtOptions = promotedHmmt.renderRequest.diagramOptions;
assert.equal(promotedHmmt.renderRequest.schema, 5);
assert.equal(hmmtOptions.configOverrides['labels.circular.scope'], 'outer');
assert.equal(hmmtOptions.configOverrides['objects.definition.circular.font_size'], 28);
assert.equal(hmmtOptions.featureShapes.repeat_region, 'underlay');
assert.ok(hmmtOptions.tracks.circularTrackSlots.some((slot) => (
  slot.id === 'a_skew_2' &&
  slot.renderer === 'dinucleotide_skew' &&
  slot.params.nt === 'AT' &&
  slot.params.legend_label === 'AT skew'
)));
assert.equal(
  hmmtOptions.tracks.circularTrackSlots.find((slot) => slot.id === 'features')
    .params.lane_direction,
  'split'
);
assert.match(
  resourceText(promotedHmmt, hmmtOptions.colors.defaultColorsFile),
  /CDS\t#84b9ec/
);
assert.match(
  resourceText(promotedHmmt, hmmtOptions.qualifierPriorityFile),
  /CDS\tgene/
);

const bgc = await loadSession('BGC0000708-BGC0000713.gbdraw-session.json');
const promotedBgc = promoteGallerySessionToCurrent(bgc);
const bgcOptions = promotedBgc.renderRequest.diagramOptions;
assert.deepEqual(
  promotedBgc.renderRequest.records.map((record) => record.presentation.reverseComplement),
  [false, false, false, false, false]
);
assert.deepEqual(
  promotedBgc.renderRequest.records.map((record) => record.presentation.label),
  [
    '<i>Streptomyces lividus</i> CBS 844.73',
    '<i>Streptomyces fradiae</i> ATCC 10745',
    '<i>Streptomyces fradiae</i> MCIMB 8233',
    '<i>Streptomyces rimosus</i> subsp. <i>paromomycinus</i> NRRL 2455',
    '<i>Streptomyces ribosidificus</i> ATCC 21294'
  ]
);
assert.deepEqual(
  promotedBgc.renderRequest.records.map((record) => record.presentation.subtitle),
  [
    'Lividomycin biosynthetic gene cluster',
    'Neomycin biosynthetic gene cluster',
    'Neomycin biosynthetic gene cluster',
    'Paromomycin biosynthetic gene cluster',
    'Ribostamycin biosynthetic gene'
  ]
);
assert.equal(
  bgcOptions.plotTitle,
  'Aminoglycoside biosynthetic gene clusters from <i>Streptomyces</i> spp.'
);
assert.equal(
  bgcOptions.configOverrides[
    'objects.definition.linear.line_styles.name.font_weight'
  ],
  'bold'
);
assert.match(
  resourceText(promotedBgc, bgcOptions.colors.colorTableFile),
  /Core biosynthetic genes/
);
assert.deepEqual(
  promotedBgc.renderRequest.comparisons,
  bgc.renderRequest.comparisons
);

const wssv = await loadSession('WSSV_genome_comparison.gbdraw-session.json');
const wssvSchema2 = structuredClone(wssv);
wssvSchema2.version = 33;
wssvSchema2.renderRequest.schema = 2;
delete wssvSchema2.renderRequest.grouping;
wssvSchema2.renderRequest.output = Array.isArray(wssvSchema2.renderRequest.output)
  ? wssvSchema2.renderRequest.output[0]
  : wssvSchema2.renderRequest.output;
wssvSchema2.renderRequest.diagramOptions.output.outputPrefix =
  wssvSchema2.renderRequest.output.prefix;
const wssvWithoutCanonicalConservation = {
  ...wssvSchema2,
  renderRequest: {
    ...wssvSchema2.renderRequest,
    diagramOptions: {
      ...wssvSchema2.renderRequest.diagramOptions,
      conservationBlastFiles: []
    }
  }
};
const promotedWssv = promoteGallerySessionToCurrent(wssvWithoutCanonicalConservation);
const wssvOptions = promotedWssv.renderRequest.diagramOptions;
assert.equal(wssvOptions.conservationBlastFiles.length, 20);
assert.deepEqual(
  wssvOptions.conservationLabels,
  wssv.config.circularConservation.series.map((series) => series.label)
);
assert.deepEqual(
  wssvOptions.conservationColors,
  wssv.config.circularConservation.series.map((series) => series.color)
);
for (const ref of wssvOptions.conservationBlastFiles) {
  assert.ok(promotedWssv.resources[ref.resourceId]);
  assert.ok(resourceText(promotedWssv, ref).length > 0);
}
const invalidWssvCacheEntries = wssv.losatCache.entries.map((entry, index) => (
  index === 0 ? { ...entry, flow: 'linear-comparison' } : entry
));
assert.throws(
  () => promoteGallerySessionToCurrent({
    ...wssvWithoutCanonicalConservation,
    losatCache: { ...wssv.losatCache, entries: invalidWssvCacheEntries }
  }),
  /20 series but 19 reusable LOSAT result/
);

const majani = await loadSession('majanivirus_orthogroup.gbdraw-session.json.gz');
const promotedMajani = promoteGallerySessionToCurrent(majani);
assert.equal(promotedMajani.renderRequest.diagramOptions.output.legend, 'right');
assert.deepEqual(
  promotedMajani.renderRequest.records.map((record) => record.presentation.label),
  [
    'Marsupenaeus japonicus endogenous nimavirus',
    'Melicertus latisulcatus majanivirus',
    'Penaeus monodon majanivirus A',
    'Penaeus semisulcatus majanivirus',
    'Penaeus monodon majanivirus B',
    'Litopenaeus vannamei majanivirus',
    'Trachysalambria curvirostris majanivirus',
    'Metapenaeus ensis majanivirus',
    'Metapenaeus joyneri majanivirus'
  ]
);
assert.deepEqual(
  promotedMajani.renderRequest.comparisons,
  majani.renderRequest.comparisons
);
assert.deepEqual(
  promotedMajani.renderRequest.diagramOptions.config,
  majani.renderRequest.diagramOptions.config
);
assert.equal(
  promotedMajani.renderRequest.diagramOptions.featureShapes.repeat_region,
  'underlay'
);

const vibrio = await loadSession('vibrio-harveyi-group-collinear.gbdraw-session.json.gz');
const promotedVibrio = promoteGallerySessionToCurrent(vibrio);
const promotedVibrioComparison = promotedVibrio.renderRequest.comparisons.find(
  (comparison) => comparison.kind === 'generatedProteinComparison'
);
assert.equal(
  promotedVibrioComparison.settings.collinearityParams.parameters.maxUnitGap,
  2
);

for (const session of [promotedBgc, promotedMajani]) {
  for (const comparison of session.renderRequest.comparisons) {
    if (!comparison.resourceId) continue;
    assert.ok(
      session.resources[comparison.resourceId],
      `missing retained comparison resource ${comparison.resourceId}`
    );
  }
}

for (const [name, promoted] of [
  ['HmmtDNA_ATskew', promotedHmmt],
  ['BGC0000708-BGC0000713', promotedBgc],
  ['WSSV_genome_comparison', promotedWssv]
]) {
  for (const [resourceId, resource] of Object.entries(promoted.resources)) {
    assert.ok(
      !resource.name.startsWith(`${resourceId}-${resourceId}-`),
      `${name} resource ${resourceId} repeats its canonical prefix`
    );
  }
  const promotedAgain = promoteGallerySessionToCurrent(promoted);
  assert.deepEqual(
    promotedAgain.renderRequest,
    promoted.renderRequest,
    `${name} renderRequest changed on a second promotion`
  );
  assert.deepEqual(
    promotedAgain.resources,
    promoted.resources,
    `${name} resources changed on a second promotion`
  );
}

console.log('gallery session migration assertions passed');
