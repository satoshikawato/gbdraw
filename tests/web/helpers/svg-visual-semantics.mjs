// Promoted from SESSION 05A13-B audit-state.js. Same-browser numbers retain
// their full precision; CLI/Wasm tolerance belongs to the separate CLI oracle.
import { stripTransientPreviewState } from '../../../gbdraw/web/js/services/svg-serialization.js';

const properties = ['x','y','x1','y1','x2','y2','cx','cy','r','rx','ry','width','height','d','points','transform','viewBox','fill','fill-opacity','fill-rule','stroke','stroke-width','stroke-opacity','stroke-dasharray','stroke-linecap','stroke-linejoin','opacity','display','visibility','font-family','font-size','font-weight','font-style','text-anchor','dominant-baseline','alignment-baseline','baseline-shift','letter-spacing','clip-path','mask','filter','marker-start','marker-mid','marker-end','href','xlink:href','offset','stop-color','stop-opacity','gradientUnits','gradientTransform','spreadMethod','preserveAspectRatio'];
const numbers = new Set(['x','y','x1','y1','x2','y2','cx','cy','r','rx','ry','width','height','d','points','transform','viewBox','fill-opacity','stroke-width','stroke-opacity','stroke-dasharray','opacity','font-size','letter-spacing','offset','stop-opacity']);
const normalize = (key, value) => numbers.has(key)
  ? value.replace(/[-+]?(?:\d*\.\d+|\d+\.?)(?:[eE][-+]?\d+)?/g, number => String(Number(number))) : value;
const attributes = element => {
  const attrs = Object.fromEntries(properties.map(key => [key, element.getAttribute(key)])
    .filter(([, value]) => value !== null).map(([key, value]) => [key, normalize(key, value)]));
  // Inline CSS has precedence over presentation attributes, including on groups.
  for (const declaration of (element.getAttribute('style') || '').split(';').filter(value => value.trim())) {
    const colon = declaration.indexOf(':');
    const key = declaration.slice(0, colon).trim(), value = declaration.slice(colon + 1).trim();
    if (key !== 'cursor') attrs[key] = normalize(key, value);
  }
  // Diagram drag adds opacity:1 to groups; this is the SVG default, unlike
  // inherited visibility/paint or any real group transform.
  if (attrs.opacity === '1') delete attrs.opacity;
  return Object.fromEntries(Object.entries(attrs).sort(([a], [b]) => a.localeCompare(b)));
};

export const visualSemanticsFromRoot = root => {
  const nodes = [root, ...root.querySelectorAll('path,text,tspan,textPath,circle,ellipse,rect,line,polyline,polygon,image,use,linearGradient,radialGradient,stop,style')];
  return nodes.map(element => {
    const ancestors = [];
    for (let parent = element.parentElement; parent; parent = parent.parentElement) {
      const attrs = attributes(parent);
      if (Object.keys(attrs).length) ancestors.push(attrs);
      if (parent === root) break;
    }
    return {
      tag: element.localName,
      feature: element.getAttribute('data-gbdraw-feature-id'),
      rendered: element.getAttribute('data-gbdraw-rendered-feature-id'),
      biological: element.getAttribute('data-gbdraw-stable-feature-id'),
      record: element.getAttribute('data-gbdraw-record-id'),
      part: element.getAttribute('data-gbdraw-feature-part'),
      label: element.getAttribute('data-label-feature-id'),
      text: ['text','tspan','textPath','style'].includes(element.localName) ? element.textContent : null,
      attrs: attributes(element), ancestors
    };
  });
};

export const svgVisualSemantics = content => {
  if (!content) return null;
  const doc = new DOMParser().parseFromString(content, 'image/svg+xml');
  if (doc.querySelector('parsererror')) throw new Error('Invalid SVG in visual comparison');
  stripTransientPreviewState(doc.documentElement);
  return visualSemanticsFromRoot(doc.documentElement);
};

export const compareVisualSemantics = (expected, actual, { catalogIds = [] } = {}) => {
  const known = new Set(catalogIds);
  if (!expected || !actual) return expected === actual ? [] : [{ property: 'presence' }];
  const differences = [];
  for (let index = 0; index < Math.max(expected.length, actual.length); index++) {
    const before = expected[index], after = actual[index];
    if (!before || !after) { differences.push({ index, before, after }); continue; }
    const { label: oldBinding, ...oldVisual } = before;
    const { label: newBinding, ...newVisual } = after;
    // Only the established missing -> valid catalog binding enrichment is
    // nonvisual. Removal, replacement, or a dangling binding is a failure.
    const validEnrichment = !oldBinding && newBinding && known.has(newBinding);
    if (JSON.stringify(oldVisual) !== JSON.stringify(newVisual)
      || (oldBinding !== newBinding && !validEnrichment)) differences.push({ index, before, after });
  }
  return differences;
};

export const nonTargetFeatures = (semantics, targetIds = []) => {
  const targets = new Set(targetIds);
  return semantics.filter(node => node.feature && !targets.has(node.feature));
};
