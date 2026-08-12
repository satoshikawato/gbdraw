const parseAttributes = (source) => {
  const attributes = new Map();
  for (const match of String(source || '').matchAll(/([\w:-]+)\s*=\s*(["'])(.*?)\2/g)) {
    attributes.set(match[1], match[3]);
  }
  return attributes;
};

const selectorParts = (selector) => String(selector || '')
  .split(/\s*,\s*/)
  .map((part) => part.trim())
  .filter(Boolean);

const matchesSelectorPart = (element, selector) => {
  let source = String(selector || '').trim();
  if (!source) return false;
  const exclusions = Array.from(source.matchAll(/:not\(([^)]+)\)/g), (match) => match[1]);
  source = source.replace(/:not\([^)]+\)/g, '');

  const tag = source.match(/^([\w:-]+)/)?.[1] || '';
  if (tag && element.localName !== tag.toLowerCase()) return false;
  const id = source.match(/#([\w:-]+)/)?.[1] || '';
  if (id && element.id !== id) return false;
  const className = source.match(/\.([\w-]+)/)?.[1] || '';
  if (className) {
    const classes = String(element.getAttribute('class') || '').split(/\s+/).filter(Boolean);
    if (!classes.includes(className)) return false;
  }

  for (const match of source.matchAll(/\[([\w:-]+)(?:\s*(\^?=)\s*(["'])(.*?)\3)?\]/g)) {
    const [, name, operator, , expected = ''] = match;
    const actual = element.getAttribute(name);
    if (actual === null) return false;
    if (operator === '=' && actual !== expected) return false;
    if (operator === '^=' && !actual.startsWith(expected)) return false;
  }
  return exclusions.every((excluded) => !matchesSelectorPart(element, excluded));
};

class FakeSvgElement {
  constructor(source) {
    const tagName = String(source || '').match(/^<\s*([\w:-]+)/)?.[1] || '';
    this.tagName = tagName;
    this.localName = tagName.toLowerCase().replace(/^.*:/, '');
    this.attributes = parseAttributes(source);
    this.children = [];
    this.parentElement = null;
    this.style = {
      removeProperty: (name) => {
        const declarations = String(this.getAttribute('style') || '')
          .split(';')
          .map((entry) => entry.trim())
          .filter(Boolean)
          .filter((entry) => entry.split(':', 1)[0].trim() !== name);
        if (declarations.length > 0) this.setAttribute('style', declarations.join('; '));
        else this.removeAttribute('style');
      }
    };
  }

  get id() {
    return this.getAttribute('id') || '';
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

  querySelectorAll(selector) {
    const selectors = selectorParts(selector);
    const matches = [];
    const visit = (node) => {
      node.children.forEach((child) => {
        if (selectors.some((part) => matchesSelectorPart(child, part))) matches.push(child);
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
    return this.querySelector(`#${id}`);
  }

  getBBox() {
    const explicit = String(this.getAttribute('data-fake-bbox') || '')
      .trim()
      .split(/[\s,]+/)
      .map(Number);
    if (explicit.length === 4 && explicit.every(Number.isFinite)) {
      return { x: explicit[0], y: explicit[1], width: explicit[2], height: explicit[3] };
    }
    const viewBox = String(this.getAttribute('viewBox') || '')
      .trim()
      .split(/[\s,]+/)
      .map(Number);
    if (viewBox.length === 4 && viewBox.every(Number.isFinite)) {
      return { x: viewBox[0], y: viewBox[1], width: viewBox[2], height: viewBox[3] };
    }
    return { x: 0, y: 0, width: 100, height: 100 };
  }
}

class FakeSvgRoot extends FakeSvgElement {
  constructor(content) {
    super(String(content || '').match(/<svg\b[^>]*>/i)?.[0] || '<svg>');
    this.localName = /<svg\b/i.test(String(content || '')) ? 'svg' : 'html';
    this.content = String(content || '');
    if (this.localName === 'svg') this.#parseChildren();
  }

  cloneNode() {
    return this;
  }

  #parseChildren() {
    const stack = [];
    const tagPattern = /<\/?([\w:-]+)\b[^>]*>/g;
    for (const match of this.content.matchAll(tagPattern)) {
      const token = match[0];
      const tagName = match[1].toLowerCase().replace(/^.*:/, '');
      if (token.startsWith('</')) {
        if (stack.length > 0) stack.pop();
        continue;
      }
      if (tagName === 'svg' && stack.length === 0) {
        stack.push(this);
        continue;
      }
      const parent = stack[stack.length - 1];
      if (!parent) continue;
      const element = new FakeSvgElement(token);
      parent.appendChild(element);
      if (!/\/\s*>$/.test(token)) stack.push(element);
    }
  }
}

export const installFakeSvgDom = () => {
  globalThis.DOMParser = class {
    parseFromString(content) {
      const root = new FakeSvgRoot(content);
      return {
        documentElement: root,
        querySelector: (selector) => (selector === 'parsererror' && root.localName !== 'svg'
          ? { textContent: 'parse error' }
          : null)
      };
    }
  };
  globalThis.XMLSerializer = class {
    serializeToString(node) {
      return String(node?.content || '');
    }
  };
};
