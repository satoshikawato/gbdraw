// OV-2732: the resources of a Session file container are file-backed (their
// bytes are the container's entries, copied to the browser, PD-OI-121) and
// hold no base64 text. Load reads the editor tables a request's projection
// reads as text (`editorTableResourceIds`) once, from their copies, and passes
// them to the projection as `resourceTexts`, which reads a table's text there
// first. A table a descriptor holds as base64 text (an older JSON Session) is
// read by the projection itself, as before.
import assert from 'node:assert/strict';
import { readFile } from 'node:fs/promises';
import { join } from 'node:path';
import test from 'node:test';
import { readSavedSession } from './helpers/session-documents.mjs';
import { savedDrawingView } from './helpers/saved-drawing-view.mjs';
import { fileBackedDescriptor } from '../../gbdraw/web/js/services/session-resource-backing.js';
import * as sessionRequest from '../../gbdraw/web/js/services/session-request.js';

const { projectCanonicalSessionRequest, readStoredRenderRequest, committedFeatureVisibilityMatches } = sessionRequest;

// The table texts, each in the form its parser reads.
const TABLES = {
  'linear-color-table': ['colorTableFile', 'CDS\tproduct\tkinase\t#123456\n'],
  'linear-whitelist': ['labelWhitelistFile', 'CDS\tgene\talpha\n'],
  'linear-priority': ['qualifierPriorityFile', 'CDS\tgene,product\n'],
  'linear-visibility': ['featureVisibilityTableFile', '*\tCDS\tproduct\tmajor capsid\toff\n'],
  'linear-label-overrides': ['labelOverrideFile', '*\t*\tlabel\t^alpha$\tRenamed alpha\n']
};

/** The Python-written Linear drawing (no draft) of a Session 51, with every editor table, its resources as base64. */
const base64Canonical = async () => {
  const document = await readSavedSession(new Blob([await readFile(join(
    process.cwd(), 'tests', 'fixtures', 'sessions', 'drawings-circular-linear.v51.gbdraw'
  ))]));
  const view = savedDrawingView(document, 'linear');
  const resources = structuredClone(view.resources);
  const options = structuredClone(view.renderRequest.diagramOptions);
  assert.equal(options.colors.defaultColors.resourceId, 'linear-colors-default-colors');
  for (const [resourceId, [field, text]] of Object.entries(TABLES)) {
    const bytes = Buffer.from(text, 'utf8');
    resources[resourceId] = {
      kind: resourceId, name: `${resourceId}.tsv`, type: 'text/tab-separated-values', size: bytes.length,
      lastModified: 0, encoding: 'base64', data: bytes.toString('base64')
    };
    if (field === 'colorTableFile') options.colors[field] = { resourceId, representation: 'file' };
    else options[field] = { resourceId, representation: 'file' };
  }
  const read = readStoredRenderRequest({ ...view.renderRequest, diagramOptions: options }, {
    featureCatalog: null, legacyOrthogroupState: null, resourceIds: new Set(Object.keys(resources)), path: 'renderRequest'
  });
  assert.ok(!('kind' in read), JSON.stringify(read));
  return { renderRequest: read.value, resources, webFiles: view.webFiles || {} };
};

/** The same resources as a container Load holds them: each a file-backed descriptor without `encoding` or `data`. */
const fileBackedResources = (resources) => {
  const reads = [];
  const backed = Object.fromEntries(Object.entries(resources).map(([resourceId, descriptor]) => {
    const { encoding: _encoding, data, ...fields } = descriptor;
    const bytes = Buffer.from(data, 'base64');
    return [resourceId, fileBackedDescriptor(fields, async () => {
      reads.push(resourceId);
      return new Uint8Array(bytes);
    })];
  }));
  return { backed, reads };
};

/** What the projection read of the editor tables. */
const tableValues = (projection) => ({
  defaultColors: projection.config.colors,
  rules: projection.config.rules,
  whitelist: projection.config.whitelist,
  priority: projection.config.qualifierPriorityRules,
  visibility: projection.semanticFeatureState.featureVisibilityManualRules,
  labelOverrides: projection.semanticFeatureState.labelOverrideRows
});

const TABLE_IDS = ['linear-colors-default-colors', ...Object.keys(TABLES)].sort();

test('editorTableResourceIds names every editor table a request projection reads as text', async () => {
  const { renderRequest } = await base64Canonical();
  assert.equal(typeof sessionRequest.editorTableResourceIds, 'function');
  assert.deepEqual([...sessionRequest.editorTableResourceIds(renderRequest)].sort(), TABLE_IDS);
  assert.deepEqual(sessionRequest.editorTableResourceIds({ diagramOptions: {} }), []);
});

test('a container drawing without a draft projects its editor tables from the texts Load read (OV-2732)', async () => {
  const canonical = await base64Canonical();
  const expected = tableValues(projectCanonicalSessionRequest(canonical));
  assert.equal(expected.rules.length, 1);
  assert.equal(expected.whitelist.length, 1);
  assert.equal(expected.labelOverrides.length, 1);
  const { backed, reads } = fileBackedResources(canonical.resources);
  // The cause: the projection's own reader takes a descriptor's base64 text only.
  assert.throws(() => projectCanonicalSessionRequest({ ...canonical, resources: backed }), /has no text payload/);
  assert.equal(typeof sessionRequest.readEditorTableTexts, 'function');
  const resourceTexts = await sessionRequest.readEditorTableTexts([canonical.renderRequest], backed);
  // Each table is read once, and no input.
  assert.deepEqual([...reads].sort(), TABLE_IDS);
  assert.deepEqual([...resourceTexts.keys()].sort(), TABLE_IDS);
  assert.deepEqual(tableValues(projectCanonicalSessionRequest({ ...canonical, resources: backed, resourceTexts })), expected);
  // A base64 descriptor is read by the projection itself, as before.
  assert.equal((await sessionRequest.readEditorTableTexts([canonical.renderRequest], canonical.resources)).size, 0);
});

test('a table that is not UTF-8 fails as it does from base64 text', async () => {
  const canonical = await base64Canonical();
  const { backed } = fileBackedResources(canonical.resources);
  backed['linear-whitelist'] = fileBackedDescriptor(
    { kind: 'linear-whitelist', name: 'w.tsv', type: '', size: 2, lastModified: 0 }, async () => new Uint8Array([0xff, 0xfe])
  );
  await assert.rejects(sessionRequest.readEditorTableTexts([canonical.renderRequest], backed), TypeError);
});

// OV-2733: Generate reuses a committed request's LOSATP results only under the Feature visibility rules it was drawn
// with; after a container Load the committed table is file-backed, so the caller passes the texts it read.
test('committed Feature visibility rules of a container Load compare from the texts read (OV-2733)', async () => {
  const canonical = await base64Canonical();
  const rules = TABLES['linear-visibility'][1];
  assert.equal(committedFeatureVisibilityMatches(canonical, rules), true);
  const { backed } = fileBackedResources(canonical.resources);
  const committed = { ...canonical, resources: backed };
  const texts = await sessionRequest.readEditorTableTexts([committed.renderRequest], backed);
  assert.equal(committedFeatureVisibilityMatches(committed, rules, [], texts), true);
  assert.equal(committedFeatureVisibilityMatches(committed, '', [], texts), false);
});
