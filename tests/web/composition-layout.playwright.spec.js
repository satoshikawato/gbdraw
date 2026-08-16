const { test, expect } = require('@playwright/test');
const { join, resolve } = require('node:path');

const repoRoot = resolve(process.env.GBDRAW_REPO || process.cwd());
const frozenV39Session = join(
  repoRoot,
  'tests',
  'fixtures',
  'sessions',
  'BGC0000708-BGC0000713.v39.gbdraw-session.json.gz'
);

test('schema 1 composition survives browser edit, drag, history, reset, export, and sanitization boundaries', async ({ page }) => {
  const runtimeRequests = [];
  page.on('request', (request) => runtimeRequests.push(request.url()));
  await page.goto('/');
  const origin = new URL(page.url()).origin;
  await page.addScriptTag({ url: '/gbdraw/web/vendor/dompurify/purify.min.js' });

  const result = await page.evaluate(async ({ origin }) => {
    const composition = await import(`${origin}/gbdraw/web/js/app/legend-layout/composition-actions.js`);
    const sanitizer = await import(`${origin}/gbdraw/web/js/services/svg-sanitization.js`);
    const roleAttribute = composition.COMPOSITION_ROLE_ATTRIBUTE;
    const target = (role, automaticTranslation, boundsKey, bounds) => ({
      automaticTranslation,
      [boundsKey]: bounds,
      role,
      selector: `[${roleAttribute}="${role}"]`
    });
    const metadata = {
      legend: target('legend', [142, 91], 'localBounds', {
        x: -2, y: -5, width: 30, height: 20
      }),
      legendReflow: {
        colorRectSize: 14,
        lineHeight: 24,
        textXOffset: 22
      },
      legendSide: 'right',
      overlayObstacles: [],
      overlayPolicy: {
        candidateScoreOrder: [
          'totalAnchorDistance',
          'xAnchorDistance',
          'yAnchorDistance',
          'nearEdgeX',
          'nearEdgeY'
        ],
        canvasGrowthCandidateOrder: ['horizontal', 'vertical'],
        canvasGrowthScoreOrder: ['addedArea', 'addedExtent', 'candidateOrder'],
        quadrantBoundaryRatio: 0.5
      },
      primary: target('primary', [16, 56], 'finalBounds', {
        x: 16, y: 56, width: 100, height: 80
      }),
      spacing: {
        dockGapPx: 24,
        edgePaddingPx: 16,
        overlayClearancePx: 8,
        stackGapPx: 20,
        titleGapPx: 20
      },
      title: target('title', [41, 26], 'localBounds', {
        x: -5, y: -10, width: 60, height: 20
      }),
      titleSide: 'top'
    };
    document.body.innerHTML = [
      '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 186 152" width="186px" height="152px">',
      `<g id="plot" ${roleAttribute}="primary" transform="translate(20,53) scale(1) rotate(0)">`,
      '<rect x="0" y="0" width="100" height="80" fill="#abcdef"/></g>',
      `<g id="legend" ${roleAttribute}="legend" transform="translate(149,96) rotate(0) matrix(1 0 0 1 0 0)">`,
      '<rect id="legend-box" x="-2" y="-5" width="30" height="20" fill="#123456"/></g>',
      `<g id="plot_title" ${roleAttribute}="title" transform="translate(39,30) scale(1)">`,
      '<rect id="title-box" x="-5" y="-10" width="60" height="20" fill="#654321"/></g>',
      '</svg>'
    ].join('');
    let svg = document.querySelector('svg');
    svg.setAttribute(composition.COMPOSITION_SCHEMA_ATTRIBUTE, '1');
    svg.setAttribute(composition.COMPOSITION_METADATA_ATTRIBUTE, JSON.stringify(metadata));

    const beforeDeltas = composition.compositionUserDeltas(svg);
    document.getElementById('legend-box').setAttribute('width', '80');
    document.getElementById('legend-box').setAttribute('height', '25');
    const next = composition.applyCompositionEdit(svg, {
      legendSide: 'bottom',
      canvasPadding: { top: 3, right: 4, bottom: 5, left: 2 }
    });
    const afterDeltas = composition.compositionUserDeltas(svg);
    const afterTransforms = {
      primary: document.getElementById('plot').getAttribute('transform'),
      legend: document.getElementById('legend').getAttribute('transform'),
      title: document.getElementById('plot_title').getAttribute('transform')
    };
    const colorAfterLayout = document.getElementById('legend-box').getAttribute('fill');
    const automaticSnapshot = new XMLSerializer().serializeToString(svg);

    composition.resetCompositionUserDeltas(svg);
    const resetSnapshot = new XMLSerializer().serializeToString(svg);
    const resetDeltas = composition.compositionUserDeltas(svg);

    document.body.innerHTML = automaticSnapshot;
    svg = document.querySelector('svg');
    const undoDeltas = composition.compositionUserDeltas(svg);
    document.body.innerHTML = resetSnapshot;
    svg = document.querySelector('svg');
    const redoDeltas = composition.compositionUserDeltas(svg);

    const unsafeExport = resetSnapshot.replace(
      '</svg>',
      '<script>alert(1)</script><g onclick="alert(2)"></g></svg>'
    );
    const sanitized = sanitizer.sanitizeSvgContent(unsafeExport, window.DOMPurify);
    document.body.innerHTML = sanitized;
    const sanitizedSvg = document.querySelector('svg');
    const sanitizedBinding = composition.bindCompositionMetadata(sanitizedSvg);
    const sanitizedRoleCount = sanitizedSvg.querySelectorAll(`[${roleAttribute}]`).length;
    const automaticDocument = new DOMParser().parseFromString(automaticSnapshot, 'image/svg+xml');

    return {
      beforeDeltas,
      afterDeltas,
      resetDeltas,
      undoDeltas,
      redoDeltas,
      afterTransforms,
      colorAfterLayout,
      metadata: next.metadata,
      viewBox: automaticDocument.documentElement.getAttribute('viewBox'),
      originalViewBox: automaticDocument.documentElement.getAttribute('data-original-view-box'),
      sanitizedSchema: sanitizedSvg.getAttribute(composition.COMPOSITION_SCHEMA_ATTRIBUTE),
      sanitizedLegendSide: sanitizedBinding.metadata.legendSide,
      sanitizedRoleCount,
      sanitizerRemovedScript: !sanitizedSvg.querySelector('script'),
      sanitizerRemovedHandler: !sanitizedSvg.querySelector('[onclick]')
    };
  }, { origin });

  expect(result.beforeDeltas).toEqual({
    primary: [[4, -3]],
    legend: [7, 5],
    title: [-2, 4]
  });
  expect(result.afterDeltas).toEqual(result.beforeDeltas);
  expect(result.undoDeltas).toEqual(result.beforeDeltas);
  expect(result.resetDeltas).toEqual({ primary: [[0, 0]], legend: [0, 0], title: [0, 0] });
  expect(result.redoDeltas).toEqual(result.resetDeltas);
  expect(result.afterTransforms.primary).toMatch(/ scale\(1\) rotate\(0\)$/);
  expect(result.afterTransforms.legend).toMatch(/ rotate\(0\) matrix\(1 0 0 1 0 0\)$/);
  expect(result.afterTransforms.title).toMatch(/ scale\(1\)$/);
  expect(result.colorAfterLayout).toBe('#123456');
  expect(result.metadata.legendSide).toBe('bottom');
  expect(result.metadata.title.finalBounds).toBeUndefined();
  const primary = result.metadata.primary.finalBounds;
  const title = result.metadata.title.localBounds;
  const titleAutomatic = result.metadata.title.automaticTranslation;
  expect(title.x + titleAutomatic[0] + title.width / 2)
    .toBeCloseTo(primary.x + primary.width / 2, 8);
  expect(result.viewBox).toMatch(/^-2 -3 /);
  expect(result.originalViewBox).toMatch(/^0 0 /);
  expect(result.sanitizedSchema).toBe('1');
  expect(result.sanitizedLegendSide).toBe('bottom');
  expect(result.sanitizedRoleCount).toBe(3);
  expect(result.sanitizerRemovedScript).toBe(true);
  expect(result.sanitizerRemovedHandler).toBe(true);
  expect(runtimeRequests.every((url) => url.startsWith(origin))).toBe(true);
});

test('legacy normalization is explicit and malformed current metadata never falls back', async ({ page }) => {
  await page.goto('/');
  const origin = new URL(page.url()).origin;
  const result = await page.evaluate(async ({ origin }) => {
    const composition = await import(`${origin}/gbdraw/web/js/app/legend-layout/composition-actions.js`);
    const legacyDocument = new DOMParser().parseFromString([
      '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 200 120"',
      ' data-horizontal-viewbox="0 0 200 120" data-vertical-viewbox="0 0 160 100">',
      '<g id="record" transform="translate(10,20) scale(1)"><rect width="40" height="20"/></g>',
      '<g id="legend" transform="translate(75,92) rotate(0)"><rect width="60" height="20"/></g>',
      '</svg>'
    ].join(''), 'image/svg+xml');
    const legacy = legacyDocument.documentElement;
    const detachedNativeBounds = (() => {
      try {
        const bounds = legacy.getElementById('record').getBBox();
        return { x: bounds.x, y: bounds.y, width: bounds.width, height: bounds.height };
      } catch (error) {
        return { error: String(error?.message || error) };
      }
    })();
    const normalized = composition.normalizeLegacyComposition(legacy, {
      legendSide: 'bottom',
      userDeltas: { primary: [3, -4], legend: [5, 2] }
    });
    const normalizedResult = {
      legacyNormalized: normalized.metadata.legacyNormalized,
      admissionPrimary: normalized.metadata.primary.finalBounds,
      admissionLegend: normalized.metadata.legend.localBounds,
      deltas: composition.compositionUserDeltas(legacy),
      primaryTransform: legacy.getElementById('record').getAttribute('transform'),
      legendTransform: legacy.getElementById('legend').getAttribute('transform')
    };

    document.body.replaceChildren(document.adoptNode(legacy));
    const liveRecord = legacy.getElementById('record');
    const liveBox = liveRecord.getBBox();
    const liveMatrix = liveRecord.getCTM();
    const liveNativeGeometry = {
      bounds: { x: liveBox.x, y: liveBox.y, width: liveBox.width, height: liveBox.height },
      matrix: Object.fromEntries(['a', 'b', 'c', 'd', 'e', 'f'].map((key) => [key, liveMatrix[key]]))
    };
    const edited = composition.applyCompositionEdit(legacy, { legendSide: 'right' });
    const editedResult = {
      primary: edited.metadata.primary.finalBounds,
      legend: edited.metadata.legend.localBounds,
      legacyNormalized: edited.metadata.legacyNormalized,
      deltas: composition.compositionUserDeltas(legacy)
    };

    const malformed = document.createElementNS('http://www.w3.org/2000/svg', 'svg');
    malformed.setAttribute(composition.COMPOSITION_SCHEMA_ATTRIBUTE, '1');
    let error = null;
    try {
      composition.normalizeLegacyComposition(malformed);
    } catch (caught) {
      error = { name: caught.name, code: caught.code, message: caught.message };
    }
    return {
      ...normalizedResult,
      detachedNativeBounds,
      liveNativeGeometry,
      edited: editedResult,
      malformedError: error,
      malformedHasMetadata: malformed.hasAttribute(composition.COMPOSITION_METADATA_ATTRIBUTE)
    };
  }, { origin });

  expect(result.legacyNormalized).toBe(true);
  expect(result.detachedNativeBounds).toEqual({ x: 0, y: 0, width: 0, height: 0 });
  expect(result.liveNativeGeometry.bounds).toEqual({ x: 0, y: 0, width: 40, height: 20 });
  expect(result.admissionPrimary).toEqual({ x: 0, y: 0, width: 200, height: 100 });
  expect(result.admissionLegend).toEqual({ x: 70, y: 90, width: 60, height: 20 });
  expect(result.deltas).toEqual({ primary: [[3, -4]], legend: [5, 2], title: null });
  expect(result.primaryTransform).toBe('translate(3,-4) translate(7,24) scale(1)');
  expect(result.legendTransform).toBe('translate(5,2) translate(70,90) rotate(0)');
  expect(result.edited.primary).toMatchObject({ width: 40, height: 20 });
  expect(result.edited.legend.x).toBeCloseTo(70, 8);
  expect(result.edited.legend.y).toBeCloseTo(90, 8);
  expect(result.edited.legend.width).toBeCloseTo(60, 8);
  expect(result.edited.legend.height).toBeCloseTo(20, 8);
  expect(result.edited.legacyNormalized).toBeUndefined();
  expect(result.edited.deltas).toEqual(result.deltas);
  expect(result.malformedError).toMatchObject({
    name: 'CompositionMetadataError',
    code: 'INVALID_COMPOSITION_METADATA'
  });
  expect(result.malformedError.message).toMatch(/incomplete/i);
  expect(result.malformedHasMetadata).toBe(false);
});

test('detached Legend admission preserves renderer anchors and updates wrapped composition geometry', async ({ page }) => {
  await page.goto('/');
  const origin = new URL(page.url()).origin;
  await page.addScriptTag({ url: '/gbdraw/web/vendor/dompurify/purify.min.js' });

  const result = await page.evaluate(async ({ origin }) => {
    const { createEditorSvgProjection } = await import(
      `${origin}/gbdraw/web/js/app/editor-svg-projection.js`
    );
    const { ingestSvgResult } = await import(
      `${origin}/gbdraw/web/js/services/svg-result-ingestion.js`
    );
    const { prepareReflowResultCommit } = await import(
      `${origin}/gbdraw/web/js/app/candidate-render.js`
    );
    const role = 'data-gbdraw-composition-role';
    const target = (name, automaticTranslation, boundsKey, bounds) => ({
      automaticTranslation,
      [boundsKey]: bounds,
      role: name,
      selector: `[${role}="${name}"]`
    });
    const metadata = (legendBounds) => ({
      legend: target('legend', [16, 16], 'localBounds', legendBounds),
      legendReflow: { colorRectSize: 14, lineHeight: 24, textXOffset: 22 },
      legendSide: 'top',
      overlayObstacles: [],
      overlayPolicy: {
        candidateScoreOrder: [
          'totalAnchorDistance',
          'xAnchorDistance',
          'yAnchorDistance',
          'nearEdgeX',
          'nearEdgeY'
        ],
        canvasGrowthCandidateOrder: ['horizontal', 'vertical'],
        canvasGrowthScoreOrder: ['addedArea', 'addedExtent', 'candidateOrder'],
        quadrantBoundaryRatio: 0.5
      },
      primary: target('primary', [16, 78], 'finalBounds', {
        x: 16, y: 78, width: 300, height: 100
      }),
      spacing: {
        dockGapPx: 24,
        edgePaddingPx: 16,
        overlayClearancePx: 8,
        stackGapPx: 20,
        titleGapPx: 20
      },
      title: null,
      titleSide: 'none'
    });
    const entry = (caption, x, y) => [
      `<g data-legend-key="${caption}">`,
      `<path d="M 0,-7 L 14,-7 L 14,7 L 0,7 z" fill="#123456" stroke="#222" stroke-width="2" transform="translate(${x - 22},${y})"/>`,
      `<text dominant-baseline="central" font-family="sans-serif" font-size="14" text-anchor="start" transform="translate(${x},${y})">${caption}</text>`,
      '</g>'
    ].join('');
    const source = ({ entries, legendBounds }) => {
      const composition = JSON.stringify(metadata(legendBounds))
        .replaceAll('&', '&amp;')
        .replaceAll('"', '&quot;');
      return [
        `<svg xmlns="http://www.w3.org/2000/svg" data-gbdraw-composition-schema="1" data-gbdraw-composition="${composition}" viewBox="0 0 332 194" width="332px" height="194px">`,
        `<g ${role}="primary" id="plot" transform="translate(16,78)"><rect width="300" height="100" fill="#ddd"/></g>`,
        `<g ${role}="legend" id="legend" transform="translate(16,16)">`,
        '<path d="M -1,-1 L 261,-1 L 261,49 L -1,49 z" fill="none" stroke="none" stroke-width="0"/>',
        ...entries,
        '</g></svg>'
      ].join('');
    };
    const admit = (content, canonical) => {
      const projection = createEditorSvgProjection(canonical);
      return ingestSvgResult(
        { content, name: 'detached.svg', format: 'svg' },
        { transformSvg: (svg) => projection.project(svg).changed }
      );
    };
    const read = (content) => {
      const document = new DOMParser().parseFromString(content, 'image/svg+xml');
      const svg = document.documentElement;
      const entries = Object.fromEntries(Array.from(svg.querySelectorAll('g[data-legend-key]')).map(
        (group) => [
          group.getAttribute('data-legend-key'),
          group.querySelector('text').getAttribute('transform')
        ]
      ));
      return {
        entries,
        metadata: JSON.parse(svg.getAttribute('data-gbdraw-composition')),
        viewBox: svg.getAttribute('viewBox')
      };
    };
    const detachedProbe = new DOMParser().parseFromString(
      source({ entries: [entry('A', 22, 7)], legendBounds: { x: -1, y: -1, width: 70, height: 16 } }),
      'image/svg+xml'
    ).getElementById('legend').getBBox();

    const singleCanonical = {
      legendEntries: [
        { caption: 'A', color: '#123456' },
        { caption: 'B', color: '#123456' }
      ],
      originalLegendOrder: ['A'],
      addedLegendCaptions: ['B']
    };
    const single = read(admit(source({
      entries: [entry('A', 22, 7)],
      legendBounds: { x: -1, y: -1, width: 70, height: 16 }
    }), singleCanonical).content);

    const wrappedSource = source({
      entries: [entry('A', 60, 7), entry('B', 160, 7), entry('C', 135, 31)],
      legendBounds: { x: 37, y: -1, width: 168, height: 40 }
    });
    const wrappedCanonical = {
      legendEntries: [
        { caption: 'A', color: '#123456' },
        { caption: 'B', color: '#123456' },
        { caption: 'C', color: '#123456' },
        { caption: 'D', color: '#123456' },
        { caption: 'A deliberately wide added Legend caption', color: '#123456' }
      ],
      originalLegendOrder: ['A', 'B', 'C'],
      addedLegendCaptions: ['D', 'A deliberately wide added Legend caption']
    };
    const admittedWrappedResult = admit(wrappedSource, wrappedCanonical);
    const wrapped = read(admittedWrappedResult.content);
    const measurementContext = document.createElement('canvas').getContext('2d');
    measurementContext.font = 'normal normal 14px sans-serif';
    const expectedDAnchorX = 135 + 14 + measurementContext.measureText('C').width + 44;
    const repeated = prepareReflowResultCommit({
      results: [{
        content: admittedWrappedResult.content,
        name: 'detached.svg',
        format: 'svg'
      }],
      ...wrappedCanonical
    });
    const reflowed = read(repeated.results[0].content);
    return {
      detachedProbe: {
        x: detachedProbe.x,
        y: detachedProbe.y,
        width: detachedProbe.width,
        height: detachedProbe.height
      },
      single,
      expectedDAnchorX,
      wrapped,
      reflowed
    };
  }, { origin });

  expect(result.detachedProbe).toEqual({ x: 0, y: 0, width: 0, height: 0 });
  expect(result.single.entries.A).toBe('translate(22,7)');
  expect(result.single.entries.B).toMatch(/^translate\([^,]+, 7\)$/);

  expect(result.wrapped.entries.A).toBe('translate(60,7)');
  expect(result.wrapped.entries.B).toBe('translate(160,7)');
  expect(result.wrapped.entries.C).toBe('translate(135,31)');
  const position = (transform) => transform.match(/[+-]?[\d.]+/g).map(Number);
  const [, cY] = position(result.wrapped.entries.C);
  const [dX, dY] = position(result.wrapped.entries.D);
  const [, wideY] = position(
    result.wrapped.entries['A deliberately wide added Legend caption']
  );
  expect(dX).toBeCloseTo(result.expectedDAnchorX, 8);
  expect(dY).toBe(cY);
  expect(wideY).toBeGreaterThan(dY);
  expect(result.wrapped.metadata.legend.localBounds.height).toBeGreaterThan(40);
  expect(result.wrapped.viewBox).not.toBe('0 0 332 194');
  expect(result.reflowed.entries).toEqual(result.wrapped.entries);
  expect(result.reflowed.metadata).toEqual(result.wrapped.metadata);
  expect(result.reflowed.viewBox).toBe(result.wrapped.viewBox);
});

test('fill-only Legend replay preserves real renderer composition geometry', async ({ page }) => {
  await page.goto('/');
  const result = await page.evaluate(async () => {
    const source = await fetch('/examples/tutorial-9-arrow-geometry-linear.svg').then(
      (response) => response.text()
    );
    const svg = new DOMParser().parseFromString(source, 'image/svg+xml').documentElement;
    const entry = svg.querySelector('g[data-legend-key]');
    const caption = entry?.getAttribute('data-legend-key') || '';
    const swatch = Array.from(entry?.querySelectorAll('path') || []).find((path) => {
      const fill = path.getAttribute('fill');
      return fill && fill !== 'none' && !fill.startsWith('url(');
    });
    const geometry = () => ({
      composition: svg.getAttribute('data-gbdraw-composition'),
      height: svg.getAttribute('height'),
      originalViewBox: svg.getAttribute('data-original-view-box'),
      transforms: Array.from(svg.querySelectorAll('[data-gbdraw-composition-role]')).map(
        (element) => [
          element.getAttribute('data-gbdraw-composition-role'),
          element.getAttribute('transform')
        ]
      ),
      viewBox: svg.getAttribute('viewBox'),
      width: svg.getAttribute('width')
    });
    const beforeGeometry = geometry();
    const beforeFill = swatch?.getAttribute('fill') || '';
    const nextFill = beforeFill.toLowerCase() === '#abcdef' ? '#123456' : '#abcdef';
    const { createEditorSvgProjection } = await import(
      '/gbdraw/web/js/app/editor-svg-projection.js'
    );
    const projection = createEditorSvgProjection({
      legendColorOverrides: { [caption]: nextFill }
    });
    const projectionResult = projection.project(svg);
    return {
      afterFill: swatch?.getAttribute('fill') || '',
      afterGeometry: geometry(),
      beforeFill,
      beforeGeometry,
      caption,
      nextFill,
      projectionResult
    };
  });

  expect(result.caption).not.toBe('');
  expect(result.beforeFill).not.toBe(result.nextFill);
  expect(result.afterFill).toBe(result.nextFill);
  expect(result.projectionResult.changed).toBe(true);
  expect(result.afterGeometry).toEqual(result.beforeGeometry);
});

test('frozen v39 session is admitted once with usable composition geometry', async ({ page }) => {
  test.setTimeout(180000);
  await page.goto('/gbdraw/web/index.html', { waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.__GBDRAW_APP__);
  const input = page.locator(
    'input[type="file"][accept*="application/json"][accept*="application/gzip"]'
  );
  const dialogPromise = page.waitForEvent('dialog', { timeout: 120000 });
  await input.setInputFiles(frozenV39Session);
  const dialog = await dialogPromise;
  expect(dialog.message()).toBe('Session loaded successfully!');
  await dialog.accept();

  await page.waitForFunction(() => {
    const app = window.__GBDRAW_APP__;
    return Boolean(
      app?.results?.length === 1
      && document.querySelector('svg[data-gbdraw-composition-schema="1"]')
    );
  }, null, { timeout: 120000 });

  const admitted = await page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    const content = String(app.results[0].content || '');
    const resultDocument = new DOMParser().parseFromString(content, 'image/svg+xml');
    const resultSvg = resultDocument.documentElement;
    const metadata = JSON.parse(resultSvg.getAttribute('data-gbdraw-composition'));
    const liveSvg = document.querySelector('svg[data-gbdraw-composition-schema="1"]');
    const liveRecordBounds = liveSvg.querySelector('g[data-gbdraw-record-index="0"]').getBBox();
    return {
      mainRuntimeStatePresent: Object.prototype.hasOwnProperty.call(app, 'pyodideReady'),
      viewBox: resultSvg.getAttribute('viewBox'),
      metadata,
      liveRecordBounds: {
        width: liveRecordBounds.width,
        height: liveRecordBounds.height
      }
    };
  });

  expect(admitted.mainRuntimeStatePresent).toBe(false);
  expect(admitted.viewBox).toBe('0 0 2480.0 718.8395510426063');
  expect(admitted.metadata.legacyNormalized).toBe(true);
  expect(admitted.metadata.primary.finalBounds).toEqual({
    x: 0,
    y: 0,
    width: 2480,
    height: 661.849967709273
  });
  expect(admitted.metadata.legend.localBounds.x).toBeCloseTo(1063.1019345238096, 8);
  expect(admitted.metadata.legend.localBounds.y).toBeCloseTo(590.0791343759396, 8);
  expect(admitted.metadata.legend.localBounds.width).toBeCloseTo(353.79613095238096, 8);
  expect(admitted.metadata.legend.localBounds.height).toBeCloseTo(56.98958333333337, 8);
  expect(admitted.liveRecordBounds.width).toBeGreaterThan(0);
  expect(admitted.liveRecordBounds.height).toBeGreaterThan(0);

  await page.evaluate(() => {
    window.__GBDRAW_APP__.form.legend = 'right';
  });
  await page.waitForFunction(() => {
    const app = window.__GBDRAW_APP__;
    const svg = document.querySelector('svg[data-gbdraw-composition-schema="1"]');
    if (!svg || app?.form?.legend !== 'right') return false;
    const metadata = JSON.parse(svg.getAttribute('data-gbdraw-composition'));
    return metadata.legendSide === 'right' && metadata.legacyNormalized !== true;
  }, null, { timeout: 120000 });

  const edited = await page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    const liveSvg = document.querySelector('svg[data-gbdraw-composition-schema="1"]');
    const liveMetadata = JSON.parse(liveSvg.getAttribute('data-gbdraw-composition'));
    const resultDocument = new DOMParser().parseFromString(
      String(app.results[0].content || ''),
      'image/svg+xml'
    );
    const resultMetadata = JSON.parse(
      resultDocument.documentElement.getAttribute('data-gbdraw-composition')
    );
    return { liveMetadata, resultMetadata, errorLog: app.errorLog };
  });
  expect(edited.errorLog).toBeNull();
  expect(edited.liveMetadata.legendSide).toBe('right');
  expect(edited.liveMetadata.legacyNormalized).toBeUndefined();
  expect(edited.liveMetadata.primary.finalBounds.width).toBeCloseTo(2280, 3);
  expect(edited.liveMetadata.primary.finalBounds.height).toBeCloseTo(568.1312770843506, 3);
  expect(edited.resultMetadata).toEqual(edited.liveMetadata);
});
