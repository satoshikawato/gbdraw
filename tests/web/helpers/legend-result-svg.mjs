// A mounted Result's feature Legend for node tests of the Legend owner: rows
// of Python's and editor rows (`data-legend-owner`), each a group with its
// swatch path and its text, placed by `translate`.
export class Node {
  constructor(tag, attributes = {}, children = [], text = '') {
    this.tagName = tag;
    this.localName = tag;
    this.attributes = new Map(Object.entries(attributes));
    this.children = [];
    this.parentElement = null;
    this.textContent = text;
    children.forEach((child) => { child.parentElement = this; this.children.push(child); });
  }

  get id() { return this.getAttribute('id') || ''; }
  getAttribute(name) { return this.attributes.has(name) ? this.attributes.get(name) : null; }
  hasAttribute(name) { return this.attributes.has(name); }
  setAttribute(name, value) { this.attributes.set(name, String(value)); }
  matches(selector) {
    if (selector === 'text' || selector === 'path') return this.tagName === selector;
    if (selector === 'g[data-legend-key]') return this.tagName === 'g' && this.hasAttribute('data-legend-key');
    if (selector.startsWith('#')) return this.id === selector.slice(1);
    return false;
  }
  querySelectorAll(selector) {
    return this.children.flatMap((child) => [...(child.matches(selector) ? [child] : []), ...child.querySelectorAll(selector)]);
  }
  querySelector(selector) { return this.querySelectorAll(selector)[0] || null; }
  getElementById(id) { return this.querySelector(`#${id}`); }
}
export const row = (key, fill, y, attributes = {}, swatch = {}) => new Node('g', { 'data-legend-key': key, transform: `translate(0,${y})`, ...attributes },
  [new Node('path', { fill, ...swatch }), new Node('text', {}, [], key)]);

/** A Result whose feature Legend holds `rows`. */
export const svgOf = (rows) => new Node('svg', {}, [new Node('g', { id: 'legend' }, [new Node('g', { id: 'feature_legend' }, rows)])]);
