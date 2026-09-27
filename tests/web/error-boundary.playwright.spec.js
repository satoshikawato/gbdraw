const { test, expect } = require('@playwright/test');
const { openApp, getDiagramWorkerActivity } = require('./helpers/app-lifecycle.cjs');

test('real Python adapters retain causes through one lazy Worker and successful retry', async ({ page }) => {
  test.setTimeout(180000);
  const consoleMessages = [];
  page.on('console', message => consoleMessages.push(message.text()));
  page.on('pageerror', error => consoleMessages.push(error.message));
  await openApp(page);
  expect((await getDiagramWorkerActivity(page)).constructions).toBe(0);
  const observations = await page.evaluate(async () => {
    const service = await import('/gbdraw/web/js/services/diagram-generation.js');
    const { normalizeUserFacingError } = await import('/gbdraw/web/js/services/error-normalization.js');
    const pattern = '😀[PRIVATE_PATTERN_SENTINEL';
    const failures = [];
    for (const kind of ['color', 'label']) {
      for (const features of [[], [{ type: 'CDS', qualifiers: { product: ['unrelated'] }, selector: {}, record: 'PRIVATE_RECORD_SENTINEL' }]]) {
        const rules = kind === 'color' ? [{ feat: 'CDS', qual: 'product', val: pattern }]
          : [{ recordId: '*', featureType: 'CDS', qualifier: 'product', valueRegex: pattern }];
        try {
          await service.runDiagramHelperOperation(service.DIAGRAM_HELPER_OPERATIONS.EVALUATE_RULES, { kind, features, rules });
          throw new Error('Syntax failure was accepted');
        } catch (error) {
          const model = normalizeUserFacingError(error);
          failures.push({ kind, model, again: normalizeUserFacingError(model) });
        }
      }
    }
    const session = await (await fetch('/gbdraw/web/gallery/sessions/HmmtDNA_basic_circular.gbdraw-session.json')).json();
    const request = structuredClone(session.renderRequest);
    request.diagramOptions.colors.colorTable = { resourceId: 'error-colors', representation: 'canonicalTsv' };
    const resource = (value) => {
      const bytes = new TextEncoder().encode(`feature_type\tqualifier_key\tvalue\tcolor\tcaption\nCDS\tproduct\t${value}\t#ff0000\t\n`);
      return { kind: 'canonical-tsv', name: 'PRIVATE_FILE_SENTINEL.tsv', encoding: 'base64',
        type: 'text/tab-separated-values', size: bytes.length, data: btoa(String.fromCharCode(...bytes)) };
    };
    const failed = await service.runDiagramGeneration({ request,
      resources: { ...session.resources, 'error-colors': resource(pattern) } });
    const render = normalizeUserFacingError(failed.results.error);
    const successful = await service.runDiagramGeneration({ request,
      resources: { ...session.resources, 'error-colors': resource('(?i)NADH') } });
    const helperRetry = await service.runDiagramHelperOperation(service.DIAGRAM_HELPER_OPERATIONS.EVALUATE_RULES, {
      kind: 'color', features: [], rules: [{ feat: 'CDS', qual: 'product', val: '(?P<enzyme>NADH)' }]
    });
    return { failures, render, success: successful.results.length,
      hasCatalog: Boolean(successful.metadata.featureCatalog), helperRetry: helperRetry.result };
  });
  for (const { kind, model, again } of observations.failures) {
    expect(model.code).toBe('REGEX_SYNTAX');
    expect(model.operation).toBe('evaluateRules');
    expect(model.stage).toBe('rule-validation');
    expect(model.context.position).toBe(1);
    expect(model.context.positionUnit).toBe('python-character');
    if (kind === 'label') expect(model.context.row).toBe(1);
    expect(again).toEqual(model);
  }
  expect(observations.render.code).toBe('REGEX_SYNTAX');
  expect(observations.render.operation).toBe('generate');
  expect(observations.render.stage).toBe('rule-validation');
  expect(observations.render.context.position).toBe(1);
  expect(observations.success).toBeGreaterThan(0);
  expect(observations.hasCatalog).toBe(true);
  expect(observations.helperRetry.winners).toEqual([]);
  expect(JSON.stringify(observations)).not.toContain('PRIVATE_');
  expect(consoleMessages.join('\n')).not.toContain('PRIVATE_');
  const activity = await getDiagramWorkerActivity(page);
  expect(activity.constructions).toBe(1);
  expect(activity.instances[0].initializations).toBe(1);
});
