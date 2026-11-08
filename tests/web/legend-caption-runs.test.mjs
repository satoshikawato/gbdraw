// A Legend caption written with <i> draws its italic runs as tspans (OV-2760).
// The Web splits a caption as Python does (gbdraw/legend/metrics.py
// `legend_caption_runs`, tests/test_legend_caption_markup.py holds the same
// cases), and the naming pass writes and reads a row's drawn name through one
// writer and one reader, so a caption drawn as tspans reads back as its markup.
import assert from 'node:assert/strict';
import test from 'node:test';

import { legendCaptionRuns } from '../../gbdraw/web/js/services/legend-layout.js';
import { drawnLegendRowName, writeDrawnLegendRowName } from '../../gbdraw/web/js/services/legend-svg.js';
import { installFakeSvgDom } from './fake-svg-dom.mjs';

const runs = (...pairs) => pairs.map(([text, style]) => ({ text, italic: style === 'italic' }));

test('a caption splits into italic runs only inside the italic subset, as Python splits it', () => {
  const cases = [
    ['CDS', runs(['CDS', 'normal'])],
    ['', runs(['', 'normal'])],
    ['<i>Danio rerio</i> (zebrafish)', runs(['Danio rerio', 'italic'], [' (zebrafish)', 'normal'])],
    ['<i>E. coli</i> and <i>V. cholerae</i>', runs(['E. coli', 'italic'], [' and ', 'normal'], ['V. cholerae', 'italic'])],
    ['<i>Homo sapiens</i>', runs(['Homo sapiens', 'italic'])],
    ['<i>a\tb\nc</i> 🧬', runs(['a\tb\nc', 'italic'], [' 🧬', 'normal'])],
    ['<b>x</b>', runs(['<b>x</b>', 'normal'])],
    ['a & <i>b</i>', runs(['a & <i>b</i>', 'normal'])],
    ['<i>a &amp; b</i>', runs(['<i>a &amp; b</i>', 'normal'])],
    ['<i></i>x', runs(['<i></i>x', 'normal'])],
    ['<i><i>a</i></i>', runs(['<i><i>a</i></i>', 'normal'])],
    ['<I>a</I>', runs(['<I>a</I>', 'normal'])],
    ['<i>a', runs(['<i>a', 'normal'])],
    ['a</i>', runs(['a</i>', 'normal'])],
    ['1 < 2', runs(['1 < 2', 'normal'])],
    ['<i>a\r</i>', runs(['<i>a\r</i>', 'normal'])],
    ['<i>a\ud800</i>', runs(['<i>a\ud800</i>', 'normal'])]
  ];
  cases.forEach(([caption, expected]) => assert.deepEqual(legendCaptionRuns(caption), expected, JSON.stringify(caption)));
});

test('the drawn name a row is written with reads back from its tspans', () => {
  installFakeSvgDom();
  const svg = new DOMParser().parseFromString('<svg><g data-legend-key="k"><path/><text>CDS</text></g></svg>').documentElement;
  const row = svg.querySelector('g[data-legend-key]');
  const label = row.querySelector('text');
  assert.equal(drawnLegendRowName(row), 'CDS');
  ['<i>Danio rerio</i> (zebrafish)', '<i>Homo sapiens</i>', 'x <i>a</i><i>b</i> y', 'a & <i>b</i>', 'tRNA'].forEach((name) => {
    writeDrawnLegendRowName(label, name);
    assert.equal(drawnLegendRowName(row), name, name);
  });
  writeDrawnLegendRowName(label, '<i>Danio rerio</i> (zebrafish)');
  assert.deepEqual(label.children.map((span) => [span.localName, span.getAttribute('font-style'), span.textContent]), [
    ['tspan', 'italic', 'Danio rerio'], ['tspan', null, ' (zebrafish)']
  ]);
  writeDrawnLegendRowName(label, 'a & <i>b</i>');
  assert.equal(label.children.length, 0, 'a caption outside the subset is one text node');
});

test('a Python row drawn as italic runs reads as the markup Python drew', () => {
  installFakeSvgDom();
  const svg = new DOMParser().parseFromString(
    '<svg><g data-legend-key="&lt;i&gt;Danio rerio&lt;/i&gt; (zebrafish)"><path/><text><tspan font-style="italic">Danio rerio</tspan><tspan> (zebrafish)</tspan></text></g></svg>'
  ).documentElement;
  assert.equal(drawnLegendRowName(svg.querySelector('g[data-legend-key]')), '<i>Danio rerio</i> (zebrafish)');
});

test('the naming pass keeps a Python italic caption and draws a renamed caption in italics', async () => {
  installFakeSvgDom();
  const { createEmptySvgMutationPlan, reconcileMountedResult } = await import('../../gbdraw/web/js/services/svg-result-ingestion.js');
  const svg = new DOMParser().parseFromString(
    '<svg viewBox="0 0 100 100"><g id="legend"><g id="feature_legend">'
    + '<g data-legend-key="Danio"><path fill="#123456"/>'
    + '<text><tspan font-style="italic">Danio rerio</tspan><tspan> (zebrafish)</tspan></text></g>'
    + '<g data-legend-key="CDS"><path fill="#654321"/><text>CDS</text></g>'
    + '</g></g></svg>'
  ).documentElement;
  // The fake DOM does not decode attribute entities; Python writes the name where it is not the key.
  svg.querySelector('g[data-legend-key]').setAttribute('data-legend-name', '<i>Danio rerio</i> (zebrafish)');
  const empty = createEmptySvgMutationPlan(1).operationsByResult[0];
  const unchanged = reconcileMountedResult(svg, empty);
  assert.equal(unchanged.changed, false, 'a caption drawn as Python drew it is not rewritten');
  const [danio, cds] = svg.querySelectorAll('g[data-legend-key]');
  assert.equal(danio.querySelector('text').children.length, 2);
  reconcileMountedResult(svg, { ...empty, legendRenames: [{ from: 'CDS', to: '<i>Homo sapiens</i>' }] });
  assert.deepEqual(cds.querySelector('text').children.map((span) => [span.getAttribute('font-style'), span.textContent]), [
    ['italic', 'Homo sapiens']
  ]);
});
