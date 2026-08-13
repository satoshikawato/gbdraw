import { expect, test } from '@playwright/test';

const projectionFixture = `
<svg xmlns="http://www.w3.org/2000/svg">
  <g data-gbdraw-slot-id="gc_content" data-gbdraw-slot-renderer="dinucleotide_content">
    <path data-testid="gc-path" fill="#aaaaaa" d="M0 0h1v1z"/>
  </g>
  <g data-gbdraw-slot-id="gc_skew" data-gbdraw-slot-renderer="dinucleotide_skew">
    <path data-testid="gc-skew-high" fill="#111111" d="M0 0h1v1z"/>
    <path data-testid="gc-skew-low" fill="#222222" d="M0 0h1v1z"/>
  </g>
  <g data-gbdraw-slot-id="at_skew" data-gbdraw-slot-renderer="dinucleotide_skew">
    <path data-testid="at-skew-high" fill="#111111" d="M0 0h1v1z"/>
    <path data-testid="at-skew-low" fill="#fedcba" d="M0 0h1v1z"/>
  </g>
  <g id="comparison1">
    <path data-gbdraw-pairwise-match-id="pair-1" data-identity-factor="0.25"
      fill="#000000" d="M0 0h1v1z"/>
    <path data-gbdraw-pairwise-match-id="block-1" data-identity-factor="0.5"
      data-collinearity-block-id="b1" data-collinearity-color-mode="orientation_identity"
      data-collinearity-orientation="plus" fill="#000000" d="M0 0h1v1z"/>
  </g>
  <g id="legend">
    <g id="feature_legend">
      <g data-legend-key="GC content"><path fill="#aaaaaa" d="M0 0h1v1z"/></g>
      <g data-legend-key="GC skew (+)"><path fill="#111111" d="M0 0h1v1z"/></g>
      <g data-legend-key="GC skew (-)"><path fill="#222222" d="M0 0h1v1z"/></g>
      <g data-legend-key="AT skew (+)"><path fill="#111111" d="M0 0h1v1z"/></g>
      <g data-legend-key="AT skew (-)"><path fill="#fedcba" d="M0 0h1v1z"/></g>
    </g>
    <g data-gbdraw-role="comparison-legend" data-gbdraw-orientation="v">
      <g data-legend-key="Pairwise match identity">
        <linearGradient><stop stop-color="#000000"/><stop stop-color="#ffffff"/></linearGradient>
      </g>
      <g data-legend-key="Collinear identity">
        <linearGradient><stop stop-color="#000000"/><stop stop-color="#ffffff"/></linearGradient>
      </g>
    </g>
  </g>
</svg>`;

const beforeColors = {
  gc_content: '#aaaaaa',
  skew_high: '#111111',
  skew_low: '#222222',
  pairwise_match_min: '#000000',
  pairwise_match_max: '#ffffff',
  collinear_block_plus_min: '#202020',
  collinear_block_plus: '#808080',
  collinear_block_minus_min: '#303030',
  collinear_block_minus: '#909090'
};

const afterColors = {
  ...beforeColors,
  gc_content: '#bbbbbb',
  skew_high: '#333333',
  skew_low: '#444444',
  pairwise_match_min: '#200000',
  pairwise_match_max: '#a00000',
  collinear_block_plus_min: '#002000',
  collinear_block_plus: '#00a000'
};

test('detached palette projection covers GC, skew, pairwise, collinear, and legends', async ({ page }) => {
  await page.goto('/gbdraw/web/index.html', { waitUntil: 'domcontentloaded' });
  const result = await page.evaluate(async ({ fixture, before, after }) => {
    const {
      buildNonFeaturePaletteProjectionContext,
      prepareNonFeaturePaletteProjection
    } = await import('/gbdraw/web/js/app/svg-styles.js');
    const source = new DOMParser().parseFromString(fixture, 'image/svg+xml').documentElement;
    const context = buildNonFeaturePaletteProjectionContext({
      defaultDinucleotide: 'GC',
      trackSlots: [
        { id: 'gc_skew', renderer: 'dinucleotide_skew', params: { nt: 'GC' } },
        {
          id: 'at_skew',
          renderer: 'dinucleotide_skew',
          params: { nt: 'AT', negative_color: '#fedcba' }
        },
        { id: 'gc_content', renderer: 'dinucleotide_content', params: { nt: 'GC' } }
      ]
    });
    const prepared = prepareNonFeaturePaletteProjection({
      svg: source,
      beforeColors: before,
      afterColors: after,
      context
    });
    const projected = prepared.svg;
    const fill = (root, selector) => root.querySelector(selector)?.getAttribute('fill') || '';
    const legendFill = (caption) => fill(
      projected,
      `g[data-legend-key="${caption}"] path`
    );
    const gradients = [...projected.querySelectorAll(
      '[data-gbdraw-role="comparison-legend"] linearGradient'
    )].map((gradient) => [...gradient.querySelectorAll('stop')]
      .map((stop) => stop.getAttribute('stop-color')));
    return {
      originalGc: fill(source, '[data-testid="gc-path"]'),
      originalPair: fill(source, '[data-gbdraw-pairwise-match-id="pair-1"]'),
      gc: fill(projected, '[data-testid="gc-path"]'),
      gcSkew: [
        fill(projected, '[data-testid="gc-skew-high"]'),
        fill(projected, '[data-testid="gc-skew-low"]')
      ],
      atSkew: [
        fill(projected, '[data-testid="at-skew-high"]'),
        fill(projected, '[data-testid="at-skew-low"]')
      ],
      pair: fill(projected, '[data-gbdraw-pairwise-match-id="pair-1"]'),
      collinear: fill(projected, '[data-gbdraw-pairwise-match-id="block-1"]'),
      legends: {
        gc: legendFill('GC content'),
        gcHigh: legendFill('GC skew (+)'),
        gcLow: legendFill('GC skew (-)'),
        atHigh: legendFill('AT skew (+)'),
        atLow: legendFill('AT skew (-)')
      },
      gradients,
      changed: prepared.changed,
      counters: prepared.counters,
      unsupported: prepared.unsupported
    };
  }, { fixture: projectionFixture, before: beforeColors, after: afterColors });

  expect(result.originalGc).toBe('#aaaaaa');
  expect(result.originalPair).toBe('#000000');
  expect(result.gc).toBe('#bbbbbb');
  expect(result.gcSkew).toEqual(['#333333', '#444444']);
  expect(result.atSkew).toEqual(['#333333', '#fedcba']);
  expect(result.pair).toBe('#400000');
  expect(result.collinear).toBe('#006000');
  expect(result.legends).toEqual({
    gc: '#bbbbbb',
    gcHigh: '#333333',
    gcLow: '#444444',
    atHigh: '#333333',
    atLow: '#fedcba'
  });
  expect(result.gradients).toEqual([
    ['#200000', '#a00000'],
    ['#002000', '#00a000']
  ]);
  expect(result.changed).toBe(true);
  expect(result.unsupported).toEqual([]);
  expect(result.counters).toMatchObject({
    gcPaths: 1,
    skewPaths: 3,
    comparisonPaths: 2,
    legendSwatches: 4,
    gradientStops: 4
  });
});

test('projection fails closed for a custom skew slot without ownership metadata', async ({ page }) => {
  await page.goto('/gbdraw/web/index.html', { waitUntil: 'domcontentloaded' });
  const result = await page.evaluate(async ({ fixture, before, after }) => {
    const { prepareNonFeaturePaletteProjection } = await import(
      '/gbdraw/web/js/app/svg-styles.js'
    );
    const source = new DOMParser().parseFromString(fixture, 'image/svg+xml').documentElement;
    try {
      prepareNonFeaturePaletteProjection({
        svg: source,
        beforeColors: before,
        afterColors: after
      });
      return { accepted: true };
    } catch (error) {
      return {
        accepted: false,
        code: error.code,
        details: error.details
      };
    }
  }, { fixture: projectionFixture, before: beforeColors, after: afterColors });

  expect(result.accepted).toBe(false);
  expect(result.code).toBe('NON_FEATURE_PALETTE_PROJECTION_UNSUPPORTED');
  expect(result.details).toContainEqual({
    kind: 'skew-slot-ownership',
    slotId: 'at_skew'
  });
});
