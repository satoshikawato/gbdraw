import assert from 'node:assert/strict';
import { test } from 'node:test';
import { readFile } from 'node:fs/promises';
import { gunzipSync } from 'node:zlib';
import { adoptCurrentSessionDocument, hasBiologicalSessionInputs, isSettingsOnlySessionDocument,
  validateSessionAuthorityInventory } from '../../gbdraw/web/js/services/session-authority.js';
import { buildSessionResources } from '../../gbdraw/web/js/services/session-resources.js';
import { adoptCurrentSessionResources } from '../../gbdraw/web/js/services/session-resource-backing.js';
import { projectSettingsOnlySession } from '../../gbdraw/web/js/services/session-request.js';
import { readFileBytes } from '../../gbdraw/web/js/services/file-content-cache.js';

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
  const resources = await buildSessionResources({ files: { qualifier_priority: file }, linearSeqs: [] }, null);
  const document = { ...structuredClone(settings), ...resources };
  assert.equal(adoptCurrentSessionDocument(document, 42).canonical, null);
  const projected = projectSettingsOnlySession(document, adoptCurrentSessionResources(document.resources));
  const restored = projected.files.qualifier_priority;
  assert.deepEqual(await readFileBytes(restored), bytes);
  assert.equal(restored.name, file.name);
  assert.equal(restored.type, file.type);
  assert.equal(restored.lastModified, file.lastModified);
  const again = await buildSessionResources({ files: projected.files, linearSeqs: [] }, null);
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

const invalid = [
  ['missing request', d => { delete d.renderRequest; }, /canonical render request/],
  ['biological binding', d => { d.resources = full.resources; d.webFiles = full.webFiles; }, /biological sources/],
  ['dangling binding', d => { d.webFiles.bindings.whitelist = { resourceId: 'absent', name: 'list.tsv', type: '', lastModified: 0 }; }, /Missing canonical resource/],
  ['unbound resource', d => { d.resources = full.resources; }, /unbound resource/],
  ['malformed config', d => { d.config.form = []; }, /active form/],
  ['unknown config', d => { d.config.adv.unsupportedSetting = true; }, /unknown.*field/],
  ['binding schema', d => { d.webFiles.bindings.schema = 99; }, /binding schema/],
  ['saved Result', d => { d.results = [{ name: 'old.svg', content: '<svg/>' }]; }, /feature catalog/],
  ['committed provenance', d => { d.cliInvocation = {}; }, /committed render artifacts/]
];
for (const [name, mutate, error] of invalid) {
  test(`settings-only admission rejects ${name}`, () => {
    const document = structuredClone(settings);
    mutate(document);
    assert.throws(() => adoptCurrentSessionDocument(document, 42), error);
  });
}

test('released schema-1 and schema-2 full documents retain admission', () => {
  validateSessionAuthorityInventory(full, 41);
  const schema2 = structuredClone(full);
  schema2.webFiles.bindings.schema = 2;
  validateSessionAuthorityInventory(schema2, 41);
});
