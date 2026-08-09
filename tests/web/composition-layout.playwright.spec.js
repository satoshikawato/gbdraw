const { test, expect } = require('@playwright/test');

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
    document.body.innerHTML = [
      '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 200 100">',
      '<g id="record" transform="translate(10,20) scale(1)"><rect width="40" height="20"/></g>',
      '<g id="legend" transform="translate(100,5) rotate(0)"><rect width="30" height="15"/></g>',
      '</svg>'
    ].join('');
    const legacy = document.querySelector('svg');
    const normalized = composition.normalizeLegacyComposition(legacy, {
      legendSide: 'right',
      userDeltas: { primary: [3, -4], legend: [5, 2] }
    });
    const normalizedResult = {
      legacyNormalized: normalized.metadata.legacyNormalized,
      deltas: composition.compositionUserDeltas(legacy),
      primaryTransform: document.getElementById('record').getAttribute('transform'),
      legendTransform: document.getElementById('legend').getAttribute('transform')
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
      malformedError: error,
      malformedHasMetadata: malformed.hasAttribute(composition.COMPOSITION_METADATA_ATTRIBUTE)
    };
  }, { origin });

  expect(result.legacyNormalized).toBe(true);
  expect(result.deltas).toEqual({ primary: [[3, -4]], legend: [5, 2], title: null });
  expect(result.primaryTransform).toBe('translate(3,-4) translate(7,24) scale(1)');
  expect(result.legendTransform).toBe('translate(5,2) translate(95,3) rotate(0)');
  expect(result.malformedError).toMatchObject({
    name: 'CompositionMetadataError',
    code: 'INVALID_COMPOSITION_METADATA'
  });
  expect(result.malformedError.message).toMatch(/incomplete/i);
  expect(result.malformedHasMetadata).toBe(false);
});
