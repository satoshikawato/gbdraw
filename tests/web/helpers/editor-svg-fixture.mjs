const selectorMatches = (element, selectorRaw) => {
  const selector = selectorRaw.trim();
  if (!selector) return false;
  if (selector.startsWith('#')) return element.id === selector.slice(1);
  if (selector === 'path' || selector === 'text' || selector === 'g') {
    return element.tagName === selector;
  }
  if (selector === 'g[data-legend-key]') {
    return element.tagName === 'g' && element.hasAttribute('data-legend-key');
  }
  const keyedLegend = selector.match(/^g\[data-legend-key="([^"]*)"\]$/);
  if (keyedLegend) {
    return element.tagName === 'g'
      && element.getAttribute('data-legend-key') === keyedLegend[1];
  }
  if (selector === '[transform]') return element.hasAttribute('transform');
  if (/^(?:path|polygon|rect)\[data-gbdraw-feature-id\]$/.test(selector)) {
    return element.tagName === selector.split('[', 1)[0]
      && element.hasAttribute('data-gbdraw-feature-id');
  }
  if (/^(?:path|polygon|rect)\[id\^="f"\]$/.test(selector)) {
    return element.tagName === selector.split('[', 1)[0] && element.id.startsWith('f');
  }
  if (selector === 'text[data-label-editable="true"]') {
    return element.tagName === 'text' && element.getAttribute('data-label-editable') === 'true';
  }
  if (selector === 'text[data-label-feature-id]') {
    return element.tagName === 'text' && element.hasAttribute('data-label-feature-id');
  }
  if (selector === 'text[data-label-source-text]') {
    return element.tagName === 'text' && element.hasAttribute('data-label-source-text');
  }
  if (selector.includes('data-gbdraw-role="comparison-legend"')) {
    if (element.getAttribute('data-gbdraw-role') !== 'comparison-legend') return false;
    const orientation = selector.match(/data-gbdraw-orientation="([^"]+)"/)?.[1];
    return !orientation || element.getAttribute('data-gbdraw-orientation') === orientation;
  }
  return false;
};

export class FakeSvgElement {
  constructor(tagName, attributes = {}, textContent = '') {
    this.tagName = tagName;
    this.attributes = new Map(
      Object.entries(attributes).map(([name, value]) => [name, String(value)])
    );
    this.children = [];
    this.parentElement = null;
    this.textContent = textContent;
    this.style = {};
  }

  get id() {
    return this.getAttribute('id') || '';
  }

  get parentNode() {
    return this.parentElement;
  }

  get childNodes() {
    return this.children;
  }

  get nextSibling() {
    if (!this.parentElement) return null;
    const index = this.parentElement.children.indexOf(this);
    return index >= 0 ? this.parentElement.children[index + 1] || null : null;
  }

  getAttribute(name) {
    return this.attributes.has(name) ? this.attributes.get(name) : null;
  }

  hasAttribute(name) {
    return this.attributes.has(name);
  }

  setAttribute(name, value) {
    this.attributes.set(name, String(value));
  }

  removeAttribute(name) {
    this.attributes.delete(name);
  }

  appendChild(child) {
    if (child.parentElement) {
      child.parentElement.children = child.parentElement.children.filter((item) => item !== child);
    }
    child.parentElement = this;
    this.children.push(child);
    return child;
  }

  insertBefore(child, reference = null) {
    if (child.parentElement) {
      child.parentElement.children = child.parentElement.children.filter((item) => item !== child);
    }
    child.parentElement = this;
    const index = reference ? this.children.indexOf(reference) : -1;
    if (index < 0) this.children.push(child);
    else this.children.splice(index, 0, child);
    return child;
  }

  replaceChild(replacement, current) {
    const index = this.children.indexOf(current);
    if (index < 0) throw new Error('replacement target is not a child');
    if (replacement.parentElement) {
      replacement.parentElement.children = replacement.parentElement.children.filter(
        (item) => item !== replacement
      );
    }
    current.parentElement = null;
    replacement.parentElement = this;
    this.children[index] = replacement;
    return current;
  }

  replaceChildren(...children) {
    this.children.forEach((child) => {
      child.parentElement = null;
    });
    this.children = [];
    children.forEach((child) => this.appendChild(child));
  }

  remove() {
    if (!this.parentElement) return;
    this.parentElement.children = this.parentElement.children.filter((item) => item !== this);
    this.parentElement = null;
  }

  cloneNode(deep = false) {
    const clone = new FakeSvgElement(
      this.tagName,
      Object.fromEntries(this.attributes),
      this.textContent
    );
    clone.style = { ...this.style };
    if (deep) this.children.forEach((child) => clone.appendChild(child.cloneNode(true)));
    return clone;
  }

  querySelectorAll(selector) {
    const selectors = String(selector).split(',');
    const matches = [];
    const visit = (node) => {
      node.children.forEach((child) => {
        if (selectors.some((candidate) => selectorMatches(child, candidate))) matches.push(child);
        visit(child);
      });
    };
    visit(this);
    return matches;
  }

  querySelector(selector) {
    return this.querySelectorAll(selector)[0] || null;
  }

  getElementById(id) {
    if (this.id === id) return this;
    let found = null;
    const visit = (node) => {
      node.children.forEach((child) => {
        if (found) return;
        if (child.id === id) found = child;
        else visit(child);
      });
    };
    visit(this);
    return found;
  }
}

export const appendFeature = (svg, id, attributes = {}) => svg.appendChild(
  new FakeSvgElement('path', {
    id,
    'data-gbdraw-feature-id': id,
    'data-gbdraw-feature-part': 'block',
    fill: '#111111',
    stroke: '#222222',
    'stroke-width': '1',
    ...attributes
  })
);

export const appendEditableLabel = (svg, featureId, sourceText) => svg.appendChild(
  new FakeSvgElement('text', {
    'data-label-editable': 'true',
    'data-label-feature-id': featureId,
    'data-label-source-text': sourceText
  }, sourceText)
);

export const appendFeatureLegend = (svg, caption, color = '#111111') => {
  const legend = svg.appendChild(new FakeSvgElement('g', { id: 'legend' }));
  const featureLegend = legend.appendChild(new FakeSvgElement('g', { id: 'feature_legend' }));
  const entry = featureLegend.appendChild(
    new FakeSvgElement('g', { 'data-legend-key': caption })
  );
  const swatch = entry.appendChild(new FakeSvgElement('path', {
    fill: color,
    stroke: '#222222',
    'stroke-width': '1'
  }));
  entry.appendChild(new FakeSvgElement('text', {}, caption));
  return { entry, featureLegend, legend, swatch };
};

export const serializeFakeSvg = (svg) => {
  const serialize = (element) => ({
    tag: element.tagName,
    attributes: Object.fromEntries([...element.attributes].sort(([left], [right]) => left.localeCompare(right))),
    text: element.textContent,
    children: element.children.map(serialize)
  });
  return JSON.stringify(serialize(svg));
};
