import assert from 'node:assert/strict';
import { test } from 'node:test';
import { readFile } from 'node:fs/promises';
import { gunzipSync } from 'node:zlib';
import { adoptCurrentSessionDocument, hasBiologicalSessionInputs, isSettingsOnlySessionDocument, readSessionDocument,
  validateSessionAuthorityInventory } from '../../gbdraw/web/js/services/session-authority.js';
import { buildSessionResources } from '../../gbdraw/web/js/services/session-resources.js';
import { adoptCurrentSessionResources } from '../../gbdraw/web/js/services/session-resource-backing.js';
import { projectSettingsOnlySession } from '../../gbdraw/web/js/services/session-request.js';
import { readFileBytes } from '../../gbdraw/web/js/services/file-content-cache.js';
import { createInputFileStore } from '../../gbdraw/web/js/services/history-files.js';
import { unreadFieldRepairs } from '../../gbdraw/web/js/services/config-field-readers.js';
import { createDefaultAdv, createDefaultForm, validateCurrentWriterActiveConfig } from '../../gbdraw/web/js/services/session-active-config-contract.js';

const settings = JSON.parse(gunzipSync(await readFile(new URL('../fixtures/sessions/settings-only.v42.json.gz', import.meta.url))));
const full = JSON.parse(await readFile(new URL('../fixtures/sessions/single.v41-bindings1.json', import.meta.url)));

test('settings-only document owns its full configuration and resources without a canonical owner', () => {
  const adopted = adoptCurrentSessionDocument(structuredClone(settings), 42);
  assert.equal(adopted.canonical, null);
  const projected = projectSettingsOnlySession(adopted.document, adoptCurrentSessionResources(settings.resources));
  assert.deepEqual(projected.config, settings.config);
  assert.equal(projected.config.form.labels_mode, 'both');
  assert.equal(projected.config.form.track_type, 'middle');
  assert.equal(projected.config.form.show_labels_linear, 'all');
  assert.equal(isSettingsOnlySessionDocument({ ...full, renderRequest: null }), false);
  assert.throws(() => adoptCurrentSessionDocument({ ...full, renderRequest: null }, 41), /canonical render request/);
});

test('auxiliary file bytes and metadata round trip through the existing resource owner', async () => {
  const bytes = new TextEncoder().encode('CDS\tgene,product,locus_tag\r\n');
  const file = Object.assign(new Blob([bytes], { type: 'text/tab-separated-values' }), {
    name: 'priority.tsv', lastModified: 123456
  });
  // The card names its file by ref of the input file store.
  const inputFiles = createInputFileStore();
  const resources = await buildSessionResources(
    { inputFiles, files: { qualifier_priority: inputFiles.registerFile(file) }, linearSeqs: [] }, null);
  const document = { ...structuredClone(settings), ...resources };
  assert.equal(adoptCurrentSessionDocument(document, 42).canonical, null);
  const projected = projectSettingsOnlySession(document, adoptCurrentSessionResources(document.resources));
  const restored = projected.files.qualifier_priority;
  assert.deepEqual(await readFileBytes(restored), bytes);
  assert.equal(restored.name, file.name);
  assert.equal(restored.type, file.type);
  assert.equal(restored.lastModified, file.lastModified);
  // Load binds the loaded file to the card by ref.
  const again = await buildSessionResources(
    { inputFiles, files: { qualifier_priority: inputFiles.registerFile(restored) }, linearSeqs: [] }, null);
  assert.deepEqual(again, resources);
});

test('inactive biological inputs cannot be classified or saved as settings-only', async () => {
  for (const inventory of [
    { c_gb: new Blob(['real source']) }, { c_gff: {} },
    { c_conservation_fastas: [null, {}] }, { c_conservation_sequence_sources: [{}] },
    { linearSeqs: [{ gb: {} }] }, { linearSeqs: [{ gff: {}, fasta: {} }] }
  ]) {
    assert.equal(hasBiologicalSessionInputs(inventory), true);
    await assert.rejects(buildSessionResources({ files: inventory, linearSeqs: inventory.linearSeqs }, null), /canonical render request/);
  }
  const inactive = structuredClone(full);
  inactive.ui = { mode: 'linear' };
  assert.ok(adoptCurrentSessionDocument(inactive, 41).canonical);
});

test('OV-566: a settings-only Session refuses a Linear row that is not an object, as Python does', () => {
  // `validateWebFileBindings` lets such a row through; Python's
  // `_validate_settings_only_session` refuses it with the same message.
  for (const rows of [[null], ['x'], [[]], [{ uid: 'a' }, null]]) {
    const document = structuredClone(settings);
    document.webFiles.bindings.linearSeqs = rows;
    assert.throws(() => adoptCurrentSessionDocument(document, 42), /biological sources/, JSON.stringify(rows));
  }
  const document = structuredClone(settings);
  document.webFiles.bindings.linearSeqs = [{ uid: 'a' }, { uid: 'b', gb: null }];
  assert.equal(adoptCurrentSessionDocument(document, 42).canonical, null);
});

// R6 (R1-M3): the binding check refuses with the typed Session-fields error.
const sessionFields = (error) => error?.code === 'INPUT_INVALID' && error.context?.reason === 'FIELDS';
const invalid = [
  ['missing request', d => { delete d.renderRequest; }, /canonical render request/],
  ['biological binding', d => { d.resources = full.resources; d.webFiles = full.webFiles; }, /biological sources/],
  ['dangling binding', d => { d.webFiles.bindings.whitelist = { resourceId: 'absent', name: 'list.tsv', type: '', lastModified: 0 }; }, sessionFields],
  ['binding schema', d => { d.webFiles.bindings.schema = 99; }, sessionFields],
  // Load read the committed Result set before (`readSessionAuthorityParts`); the settings-only check refuses it.
  ['saved Result', d => { d.results = [{ name: 'old.svg', content: '<svg/>' }]; }, /committed render artifacts/],
  ['committed provenance', d => { d.cliInvocation = {}; }, /committed render artifacts/],
  ['empty annotation warning artifact', d => { d.runMetadata = { annotationWarnings: [] }; }, /committed render artifacts/],
  ['track geometry artifact', d => { d.runMetadata = { trackSlotGeometry: {} }; }, /committed render artifacts/]
];
for (const [name, mutate, error] of invalid) {
  test(`settings-only admission rejects ${name}`, () => {
    const document = structuredClone(settings);
    mutate(document);
    assert.throws(() => adoptCurrentSessionDocument(document, 42), error);
  });
}

// E-OD02: Save refuses these; Load drops a resource no binding names with no
// row (E-OD06) and reads a setting container of another kind with a row (E1a).
for (const [name, mutate, error] of [
  ['unbound resource', d => { d.resources = full.resources; }, /unbound resource/],
  ['malformed config', d => { d.config.form = []; }, /active form/]
]) {
  test(`settings-only Save refuses ${name}, which Load reads`, () => {
    const document = structuredClone(settings);
    mutate(document);
    assert.throws(() => validateSessionAuthorityInventory(structuredClone(document), 42), error);
    assert.equal(adoptCurrentSessionDocument(document, 42).canonical, null);
  });
}

// E1a (E-OD02): Load reads an unknown setting of a settings-only Session with a
// row before the dialog; Save's check of the configuration it writes still refuses it.
test('settings-only Load lists an unknown setting, which Save refuses', () => {
  const document = structuredClone(settings);
  document.config.adv.unsupportedSetting = true;
  assert.equal(adoptCurrentSessionDocument(structuredClone(document), 42).canonical, null);
  assert.deepEqual(unreadFieldRepairs(document.config, 'config', 'circular').map(({ code, path }) => `${code}@${path}`),
    ['setting-replaced@config.adv.unsupportedSetting']);
  const written = { form: createDefaultForm(), adv: { ...createDefaultAdv('circular'), unsupportedSetting: true } };
  assert.throws(() => validateCurrentWriterActiveConfig({ mode: 'circular', storedConfig: written }), /unknown.*field/);
});

// OV-860 (E-OD06, E-OD02): Load reads a settings-only Session (a null request)
// without a Web input inventory as one that binds nothing, with no row for the
// inventory it never held; Save still refuses the inventory it lacks.
test('OV-860: settings-only Load reads a missing Web input inventory as the empty one; Save refuses it', () => {
  for (const mutate of [(d) => { delete d.webFiles.bindings; }, (d) => { delete d.webFiles; }, (d) => { delete d.webFiles.bindings.c_gff; }]) {
    const document = structuredClone(settings);
    mutate(document);
    assert.throws(() => validateSessionAuthorityInventory(structuredClone(document), 42), /explicit Web input inventory/);
    const read = readSessionDocument(document);
    assert.deepEqual(read.repairs, []);
    const { bindings } = read.value.webFiles;
    assert.deepEqual([bindings.c_gb, bindings.c_gff, bindings.c_fasta, Array.isArray(bindings.linearSeqs)], [null, null, null, true]);
    assert.equal(adoptCurrentSessionDocument(read.value, 42).canonical, null);
  }
});

test('released schema-1 and schema-2 full documents retain admission', () => {
  validateSessionAuthorityInventory(full, 41);
  const schema2 = structuredClone(full);
  schema2.webFiles.bindings.schema = 2;
  validateSessionAuthorityInventory(schema2, 41);
});
