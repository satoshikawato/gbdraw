import assert from 'node:assert/strict';
import { mkdir, mkdtemp, readFile, writeFile } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { pathToFileURL } from 'node:url';
import { join } from 'node:path';

globalThis.window = {
  Vue: {
    computed: (getter) => ({
      get value() {
        return getter();
      }
    })
  }
};

const repoRoot = process.cwd();
const sourceDir = join(repoRoot, 'gbdraw', 'web', 'js');
const tempDir = await mkdtemp(join(tmpdir(), 'gbdraw-feature-selection-'));
await writeFile(join(tempDir, 'package.json'), '{"type":"module"}\n', 'utf8');
await mkdir(join(tempDir, 'app'), { recursive: true });
await mkdir(join(tempDir, 'services'), { recursive: true });
await writeFile(
  join(tempDir, 'app', 'feature-selection.js'),
  await readFile(join(sourceDir, 'app', 'feature-selection.js'), 'utf8'),
  'utf8'
);
await writeFile(
  join(tempDir, 'services', 'svg-serialization.js'),
  await readFile(join(sourceDir, 'services', 'svg-serialization.js'), 'utf8'),
  'utf8'
);
const featureSelectionUrl = pathToFileURL(join(tempDir, 'app', 'feature-selection.js'));
const serializationUrl = pathToFileURL(join(tempDir, 'services', 'svg-serialization.js'));

const {
  getFeatureSelectionScope,
  normalizeFeatureSelectionId,
  stripFeatureSelectionClasses
} = await import(featureSelectionUrl);
const { stripTransientPreviewState } = await import(serializationUrl);

class FakeElement {
  constructor(className = '', style = '') {
    this.attributes = {};
    if (className) this.attributes.class = className;
    if (style) this.attributes.style = style;
    this.style = {
      removeProperty: (name) => {
        if (name !== 'cursor') return;
        const next = String(this.attributes.style || '')
          .split(';')
          .map((entry) => entry.trim())
          .filter((entry) => entry && !entry.toLowerCase().startsWith('cursor:'))
          .join('; ');
        if (next) {
          this.attributes.style = next;
        } else {
          delete this.attributes.style;
        }
      }
    };
    const element = this;
    this.classList = {
      remove: (token) => {
        const tokens = this.classTokens().filter((entry) => entry !== token);
        if (tokens.length) {
          this.attributes.class = tokens.join(' ');
        } else {
          delete this.attributes.class;
        }
      },
      get length() {
        return String(element.attributes.class || '').trim()
          ? String(element.attributes.class || '').trim().split(/\s+/).length
          : 0;
      }
    };
  }

  classTokens() {
    return String(this.attributes.class || '').trim().split(/\s+/).filter(Boolean);
  }

  getAttribute(name) {
    return Object.hasOwn(this.attributes, name) ? this.attributes[name] : null;
  }

  setAttribute(name, value) {
    this.attributes[name] = String(value);
  }

  removeAttribute(name) {
    delete this.attributes[name];
  }
}

const fakeSvg = (elements) => ({
  querySelectorAll: (selector) => {
    if (selector === '[style]') {
      return elements.filter((element) => element.getAttribute('style') !== null);
    }
    if (selector.includes(',')) {
      const selectors = selector.split(',').map((entry) => entry.trim());
      return elements.filter((element) =>
        selectors.some((entry) => entry.startsWith('.') && element.classTokens().includes(entry.slice(1)))
      );
    }
    if (selector.startsWith('.')) {
      const token = selector.slice(1);
      return elements.filter((element) => element.classTokens().includes(token));
    }
    return [];
  }
});

assert.equal(normalizeFeatureSelectionId('feature-1__part2'), 'feature-1');
assert.equal(
  getFeatureSelectionScope({
    fileIdx: 2,
    record_idx: 3,
    record_id: 'raw-record',
    displayRecordId: 'display-record'
  }),
  '2::3::display-record'
);

const selected = new FakeElement('gbdraw-feature-selected retained');
const anchor = new FakeElement('gbdraw-feature-selection-anchor');
const candidate = new FakeElement('gbdraw-feature-selection-candidate');
stripFeatureSelectionClasses(fakeSvg([selected, anchor, candidate]));
assert.equal(selected.getAttribute('class'), 'retained');
assert.equal(anchor.getAttribute('class'), null);
assert.equal(candidate.getAttribute('class'), null);

const searchMatch = new FakeElement('gbdraw-preview-feature-search-match gbdraw-feature-selected keep');
const cursorStyle = new FakeElement('', 'cursor: pointer; opacity: 0.8');
stripTransientPreviewState(fakeSvg([searchMatch, cursorStyle]));
assert.equal(searchMatch.getAttribute('class'), 'keep');
assert.equal(cursorStyle.getAttribute('style'), 'opacity: 0.8');

console.log('feature selection tests passed');
