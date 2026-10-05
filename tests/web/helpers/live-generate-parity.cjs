// Live-edit vs Generate parity (PD-OI-066 LIVE-EDIT-EQUALS-REGENERATION, R1,
// R3): a live edit shows on the displayed Result what the next Generate from
// the same draft draws. semanticSnapshot reads the displayed Result SVG into a
// state keyed by stable IDs, without geometry; expectLiveEqualsGenerate runs
// Generate and diffs the two states, minus the differences that
// tests/web/contracts/live-generate-parity-allowed.json allows.
const { expect, test } = require('@playwright/test');
const { generateAndWaitForResult } = require('./app-lifecycle.cjs');
const ALLOWED_DIFFERENCES = require('../contracts/live-generate-parity-allowed.json').entries;

// No Generate, label reflow, automatic rerender, rule matching, Session import,
// or History step is running, over two frames in a row.
const settleLive = async (page, { timeout = 120_000 } = {}) => {
  const idle = () => {
    const app = window.__GBDRAW_APP__;
    const history = window.__GBDRAW_HISTORY__;
    return !app.processing && !app.labelReflowProcessing && !app.sessionImportPending
      && !app.ruleMatchingPending && !history?.capturing?.value && !history?.restoring?.value
      && !history?.mutationPending?.();
  };
  for (let stableFrames = 0; stableFrames < 2;) {
    await page.waitForFunction(idle, null, { timeout });
    await page.evaluate(() => new Promise((resolve) => (
      requestAnimationFrame(() => requestAnimationFrame(resolve))
    )));
    stableFrames = await page.evaluate(idle) ? stableFrames + 1 : 0;
  }
};

// Shows the batch Result at `index` through the Result selector and waits until
// it is mounted.
const showResult = async (page, index) => {
  await page.locator('h2 select').selectOption({ index });
  await expect.poll(() => page.evaluate(async (wanted) => {
    const { isCommittedSvgResultMounted } = await import('/gbdraw/web/js/services/svg-result-ingestion.js');
    const app = window.__GBDRAW_APP__;
    return app.selectedResultIndex === wanted && isCommittedSvgResultMounted(app.results[wanted]);
  }, index)).toBe(true);
  await settleLive(page);
};

// The displayed Result as the reader sees it, keyed by stable IDs:
// - features[renderedId]: drawn parts and their fill and stroke colors;
// - labels[renderedId]: the text of each drawn label;
// - leaders[renderedId]: the drawn leader segments of the label;
// - legend: the drawn legend rows in reading order (caption, swatch colors);
// - texts: every other drawn text (ticks, definitions, titles).
// An element is drawn unless it or an ancestor has display none, or its
// nearest visibility setting is hidden or collapse; a hidden element and an
// absent one read the same. Positions, sizes, and paths are left out.
const semanticSnapshot = (page) => page.evaluate(() => {
  const app = window.__GBDRAW_APP__;
  const root = app.svgContainer?.querySelector('svg');
  if (!root) throw new Error('No displayed Result SVG');
  const property = (element, name) => {
    const inline = element.style?.getPropertyValue?.(name);
    return String(inline || element.getAttribute(name) || '').trim();
  };
  const unrendered = 'defs, clipPath, mask, pattern, marker, symbol, linearGradient, radialGradient';
  const drawn = (element) => {
    if (element.closest(unrendered)) return false;
    let visibility = '';
    for (let node = element; node && node.nodeType === 1; node = node.parentNode) {
      if (property(node, 'display') === 'none') return false;
      if (!visibility) visibility = property(node, 'visibility');
      if (node === root) break;
    }
    return !['hidden', 'collapse'].includes(visibility);
  };
  const canvas = document.createElement('canvas').getContext('2d');
  const color = (raw) => {
    const value = String(raw || '').trim();
    if (!value || ['none', 'transparent', 'currentcolor', 'inherit'].includes(value.toLowerCase())) {
      return value.toLowerCase();
    }
    const reference = value.match(/^url\(\s*['"]?#([^'")]+)['"]?\s*\)/);
    if (reference) {
      const target = root.querySelector(`#${CSS.escape(reference[1])}`);
      const stops = [...(target?.querySelectorAll('stop') || [])].map((stop) => color(property(stop, 'stop-color')));
      return `url(${stops.join(' ')})`;
    }
    canvas.fillStyle = '#010203';
    canvas.fillStyle = value;
    const normalized = String(canvas.fillStyle);
    return normalized === '#010203' && value.toLowerCase() !== '#010203' ? value.toLowerCase() : normalized;
  };
  // Fill and stroke inherit; the nearest ancestor-or-self setting applies.
  const paint = (element, name) => {
    for (let node = element; node && node.nodeType === 1; node = node.parentNode) {
      const value = property(node, name);
      if (value) return color(value);
      if (node === root) break;
    }
    return '';
  };
  const sortedUnique = (values) => [...new Set(values)].sort();
  const textOf = (element) => element.textContent.replace(/\s+/g, ' ').trim();

  const shapes = 'path, rect, polygon, polyline, circle, ellipse, line, use';
  const featureParts = new Map();
  for (const element of root.querySelectorAll(`:is(${shapes})[data-gbdraw-feature-id]`)) {
    if (!drawn(element)) continue;
    const id = element.getAttribute('data-gbdraw-rendered-feature-id') || element.getAttribute('data-gbdraw-feature-id');
    if (!featureParts.has(id)) featureParts.set(id, []);
    featureParts.get(id).push(element);
  }
  const features = Object.fromEntries([...featureParts].sort(([left], [right]) => left.localeCompare(right))
    .map(([id, elements]) => [id, {
      parts: sortedUnique(elements.map((element) => element.getAttribute('data-gbdraw-feature-part') || element.localName)),
      fill: sortedUnique(elements.map((element) => paint(element, 'fill'))),
      stroke: sortedUnique(elements.map((element) => paint(element, 'stroke')))
    }]));

  const labels = {};
  const leaders = {};
  const boundParts = [...root.querySelectorAll('[data-label-feature-id]')].filter(drawn);
  for (const element of boundParts) {
    const id = element.getAttribute('data-label-feature-id');
    if (element.localName === 'text') {
      labels[id] = { text: [...(labels[id]?.text || []), textOf(element)].sort() };
    } else if (element.matches(shapes)) {
      leaders[id] = { count: (leaders[id]?.count || 0) + 1, stroke: sortedUnique([...(leaders[id]?.stroke || []), paint(element, 'stroke')]) };
    }
  }

  // Reading order: rows top to bottom, then left to right, on the screen.
  const legendRows = [...root.querySelectorAll('g[data-legend-key]')].filter(drawn).map((entry) => {
    const text = [...entry.querySelectorAll('text')].find(drawn);
    const swatch = [...entry.querySelectorAll(shapes)].find((element) => drawn(element) && !['', 'none'].includes(paint(element, 'fill')));
    const box = (text || entry).getBoundingClientRect();
    return {
      key: entry.getAttribute('data-legend-key'),
      caption: text ? textOf(text) : '',
      fill: swatch ? paint(swatch, 'fill') : '',
      stroke: swatch ? paint(swatch, 'stroke') : '',
      x: box.x,
      y: box.y
    };
  });
  const legend = legendRows
    .sort((left, right) => (Math.abs(left.y - right.y) < 2 ? left.x - right.x : left.y - right.y))
    .map(({ x, y, ...row }) => row);

  const texts = [...root.querySelectorAll('text')]
    .filter((element) => drawn(element) && !element.closest('[data-label-feature-id], g[data-legend-key]'))
    .map(textOf).filter(Boolean).sort();

  return {
    context: {
      mode: app.mode,
      resultIndex: app.selectedResultIndex,
      resultCount: app.results.length,
      autoReflow: Boolean(app.autoLabelReflowEnabled)
    },
    features,
    labels,
    leaders,
    legend,
    texts
  };
});

// One difference per ID and attribute: { kind, id, attribute, live, generated }.
// kind is feature, label, leader, legend, or text; a missing item reads as
// attribute 'drawn' (true/false).
const diffSemanticSnapshots = (live, generated) => {
  const differences = [];
  const same = (left, right) => JSON.stringify(left) === JSON.stringify(right);
  for (const kind of ['features', 'labels', 'leaders']) {
    const singular = kind.slice(0, -1);
    const ids = [...new Set([...Object.keys(live[kind]), ...Object.keys(generated[kind])])].sort();
    for (const id of ids) {
      const left = live[kind][id];
      const right = generated[kind][id];
      if (!left || !right) {
        differences.push({ kind: singular, id, attribute: 'drawn', live: Boolean(left), generated: Boolean(right) });
        continue;
      }
      for (const attribute of [...new Set([...Object.keys(left), ...Object.keys(right)])]) {
        if (!same(left[attribute], right[attribute])) {
          differences.push({ kind: singular, id, attribute, live: left[attribute], generated: right[attribute] });
        }
      }
    }
  }
  const rows = (legend) => new Map(legend.map((row) => [row.caption, row]));
  const liveRows = rows(live.legend);
  const generatedRows = rows(generated.legend);
  for (const caption of [...new Set([...liveRows.keys(), ...generatedRows.keys()])]) {
    const left = liveRows.get(caption);
    const right = generatedRows.get(caption);
    if (!left || !right) {
      differences.push({ kind: 'legend', id: caption, attribute: 'drawn', live: Boolean(left), generated: Boolean(right) });
      continue;
    }
    for (const attribute of ['key', 'fill', 'stroke']) {
      if (left[attribute] !== right[attribute]) {
        differences.push({ kind: 'legend', id: caption, attribute, live: left[attribute], generated: right[attribute] });
      }
    }
  }
  const order = (legend) => legend.map((row) => row.caption).filter((caption) => liveRows.has(caption) && generatedRows.has(caption));
  if (!same(order(live.legend), order(generated.legend))) {
    differences.push({ kind: 'legend', id: '*', attribute: 'order', live: order(live.legend), generated: order(generated.legend) });
  }
  if (!same(live.texts, generated.texts)) {
    const missing = (left, right) => {
      const remaining = [...right];
      return left.filter((text) => {
        const index = remaining.indexOf(text);
        if (index < 0) return true;
        remaining.splice(index, 1);
        return false;
      });
    };
    differences.push({
      kind: 'text', id: '*', attribute: 'drawn', live: missing(live.texts, generated.texts), generated: missing(generated.texts, live.texts)
    });
  }
  return differences;
};

// The context an allowed-difference condition reads: the live Result's
// context, and whether the feature's label is drawn on both sides.
const differenceContext = (difference, live, generated) => ({
  ...live.context,
  labelDrawnInBoth: Boolean(live.labels[difference.id] && generated.labels[difference.id])
});

// The allowed-difference entry for a difference: its kind, one of its
// attributes, and every condition field equal to the difference's context.
const allowedBy = (difference, context, allowed) => allowed.find((entry) => (
  entry.kind === difference.kind
  && entry.attributes.includes(difference.attribute)
  && Object.entries(entry.condition || {}).every(([field, value]) => context[field] === value)
)) || null;

const describeDifference = ({ kind, id, attribute, live, generated }) => (
  `${kind} ${id} ${attribute}: live ${JSON.stringify(live)}, generated ${JSON.stringify(generated)}`
);

// Snapshots the live Result, runs Generate, snapshots again, and fails with one
// line per unexpected difference. `allowed` defaults to every entry of the
// allowed-difference table; `label` names the case in the message.
const expectLiveEqualsGenerate = async (page, { allowed = ALLOWED_DIFFERENCES, label = 'live edit' } = {}) => {
  await settleLive(page);
  const live = await semanticSnapshot(page);
  await generateAndWaitForResult(page);
  await settleLive(page);
  // Generate displays its first Result; compare the Result the edit showed.
  if (await page.evaluate(() => window.__GBDRAW_APP__.selectedResultIndex) !== live.context.resultIndex) {
    await showResult(page, live.context.resultIndex);
  }
  const generated = await semanticSnapshot(page);
  const differences = diffSemanticSnapshots(live, generated).map((difference) => ({
    ...difference,
    allowedBy: allowedBy(difference, differenceContext(difference, live, generated), allowed)?.id || null
  }));
  const unexpected = differences.filter((difference) => !difference.allowedBy);
  const tolerated = differences.filter((difference) => difference.allowedBy);
  // The report lists each allowed difference with the entry that allows it.
  for (const difference of tolerated) {
    test.info().annotations.push({
      type: 'live-generate-parity allowed',
      description: `${describeDifference(difference)} (${difference.allowedBy})`
    });
  }
  expect(unexpected.map(describeDifference),
    `${label}: the live Result differs from the Result Generate draws from the same draft `
    + `(${JSON.stringify(live.context)})`).toEqual([]);
  return { live, generated, differences, tolerated };
};

module.exports = {
  diffSemanticSnapshots,
  expectLiveEqualsGenerate,
  semanticSnapshot,
  settleLive,
  showResult
};
