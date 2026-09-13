import { createSvgStyles } from '../../../gbdraw/web/js/app/svg-styles.js';

const ref = value => ({ value });
globalThis.XMLSerializer = class { serializeToString(svg) { return svg.snapshot(); } };
export const fixture = (colors, rules) => {
  const attrs = { id: 'f1', 'data-gbdraw-feature-id': 'f1', 'data-gbdraw-feature-part': 'block', fill: '#000000' };
  const feature = { svg_id: 'f1', type: 'unlisted_type', start: 1, end: 10, qualifiers: { gene: ['sample'] } };
  const path = { getAttribute: key => attrs[key] ?? null, setAttribute: (key, value) => { attrs[key] = value; } };
  const depthAttrs = {};
  const depth = { getAttribute: key => depthAttrs[key] ?? null,
    setAttribute: (key, value) => { depthAttrs[key] = value; }, removeAttribute: key => { delete depthAttrs[key]; } };
  const rootAttrs = { xmlns: 'http://www.w3.org/2000/svg', 'xmlns:xlink': 'http://www.w3.org/1999/xlink' };
  const svg = { querySelectorAll: selector => selector.startsWith('g[') ? (selector.includes('depth') ? [depth] : []) : selector.includes('data-gbdraw-feature-id') ? [path] : [],
    getElementById: () => null, getAttribute: key => rootAttrs[key] ?? null,
    setAttribute: (key, value) => { rootAttrs[key] = value; }, removeAttribute: key => { delete rootAttrs[key]; }, cloneNode: () => svg,
    snapshot: () => JSON.stringify({ fill: attrs.fill, depth: depthAttrs }) };
  const state = { svgContent: ref('<svg/>'), extractedFeatures: ref([feature]), featuresBySvgId: ref(new Map([['f1', feature]])),
    appliedPaletteColors: ref(colors), manualSpecificRules: rules, featureColorOverrides: {}, legendColorOverrides: {},
    pairwiseMatchFactors: ref({}), results: ref([{ content: JSON.stringify({ fill: '#000000', depth: {} }) }]), selectedResultIndex: ref(0),
    skipCaptureBaseConfig: ref(false), svgContainer: ref({ querySelector: () => svg }), adv: {}, mode: ref('circular'), form: { show_depth: true } };
  const actions = createSvgStyles({ state, watch() {}, nextTick: fn => fn(), legendActions: { getAllFeatureLegendGroups: () => [] } });
  return { actions, state, attrs, depthAttrs };
};
