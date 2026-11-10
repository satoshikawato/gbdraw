import assert from 'node:assert/strict';
import { cp, mkdtemp, writeFile } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { pathToFileURL } from 'node:url';
import { join } from 'node:path';

globalThis.window = {
  Vue: {
    computed: (getter) => ({
      get value() {
        return getter();
      }
    }),
    reactive: (value) => value
  }
};

const repoRoot = process.cwd();
const sourceDir = join(repoRoot, 'gbdraw', 'web', 'js');
const tempDir = await mkdtemp(join(tmpdir(), 'gbdraw-feature-selection-'));
// A copy of the module tree: feature-selection.js and the modules it imports (svg-serialization.js resolves a
// drawing's ID by its mode; services/region-axis.js orders a Shift+click range).
await cp(sourceDir, tempDir, { recursive: true });
await writeFile(join(tempDir, 'package.json'), '{"type":"module"}\n', 'utf8');
const featureSelectionUrl = pathToFileURL(join(tempDir, 'app', 'feature-selection.js'));
const serializationUrl = pathToFileURL(join(tempDir, 'services', 'svg-serialization.js'));

const {
  createFeatureSelection,
  getFeatureSelectionScope,
  normalizeFeatureSelectionId,
  stripFeatureSelectionClasses
} = await import(featureSelectionUrl);
const { stripTransientPreviewState } = await import(serializationUrl);
const { PREVIEW_MATCH_BUTTON_ATTRIBUTE, setClassToken } = await import(pathToFileURL(join(tempDir, 'services', 'preview-affordances.js')));

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

const fakeSvg = (elements, className = '') => {
  const root = new FakeElement(className);
  root.querySelectorAll = (selector) => {
    const attribute = /^\[([\w-]+)\]$/.exec(selector);
    if (attribute) return elements.filter((element) => element.getAttribute(attribute[1]) !== null);
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
  };
  return root;
};

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
const layoutTarget = new FakeElement('gbdraw-preview-layout-target keep-layout', 'cursor: help');
const cursorStyle = new FakeElement('', 'cursor: pointer; opacity: 0.8');
// OV-1842: the preview marks the match bands it makes buttons, and the label is
// translated (ja below), so saving finds them by the marker, not by the label.
/** @param {string} label */
const previewMatchButton = (label) => {
  const element = new FakeElement();
  for (const [name, value] of Object.entries({
    'data-gbdraw-match-id': 'comparison1_match1',
    role: 'button', tabindex: '0', 'aria-label': label, [PREVIEW_MATCH_BUTTON_ATTRIBUTE]: ''
  })) element.setAttribute(name, value);
  return element;
};
const pairwise = previewMatchButton('ペアワイズマッチ1');
// Hover and drag feedback that was removed leaves an empty inline style.
const emptiedStyle = new FakeElement();
emptiedStyle.setAttribute('style', '');
const previewRoot = fakeSvg(
  [searchMatch, cursorStyle, layoutTarget, pairwise, emptiedStyle],
  'gbdraw-preview-feature-search-results-active gbdraw-preview-feature-search-updating keep-root'
);
stripTransientPreviewState(previewRoot);
assert.equal(layoutTarget.getAttribute('class'), 'keep-layout');
assert.equal(layoutTarget.getAttribute('style'), null);
assert.equal(searchMatch.getAttribute('class'), 'keep');
assert.equal(cursorStyle.getAttribute('style'), 'opacity: 0.8');
assert.equal(emptiedStyle.getAttribute('style'), null);
assert.equal(previewRoot.getAttribute('class'), 'keep-root');
assert.equal(pairwise.getAttribute('data-gbdraw-match-id'), 'comparison1_match1');
for (const name of ['role', 'tabindex', 'aria-label', PREVIEW_MATCH_BUTTON_ATTRIBUTE]) assert.equal(pairwise.getAttribute(name), null);
// OV-3241: the one writer of a Preview class removes the attribute when no class is left. A browser's
// `classList.toggle` or `remove` of the last class leaves `class=""`, which Python never draws and the SVG download
// carried after a match popup marked and unmarked a match.
class DomClassElement {
  constructor(className = null) {
    this.attributes = className === null ? {} : { class: className };
    const element = this;
    const tokens = () => String(element.attributes.class ?? '').split(/\s+/).filter(Boolean);
    this.classList = {
      toggle(token, force) {
        const next = tokens().filter((entry) => entry !== token);
        if (force) next.push(token);
        if (force || tokens().includes(token)) element.attributes.class = next.join(' ');
        return Boolean(force);
      },
      get length() { return tokens().length; }
    };
  }

  getAttribute(name) { return Object.hasOwn(this.attributes, name) ? this.attributes[name] : null; }

  removeAttribute(name) { delete this.attributes[name]; }
}
const marked = new DomClassElement();
setClassToken(marked, 'gbdraw-match-selected', true);
assert.equal(marked.getAttribute('class'), 'gbdraw-match-selected');
setClassToken(marked, 'gbdraw-match-selected', false);
assert.equal(marked.getAttribute('class'), null);
const drawnClass = new DomClassElement('gbdraw-feature');
setClassToken(drawnClass, 'gbdraw-match-pending', true);
setClassToken(drawnClass, 'gbdraw-match-pending', false);
assert.equal(drawnClass.getAttribute('class'), 'gbdraw-feature');
setClassToken(marked, 'gbdraw-match-pending', false);
assert.equal(marked.getAttribute('class'), null);
const interactivePairwise = previewMatchButton('Pairwise match 1');
stripTransientPreviewState(fakeSvg([interactivePairwise]), { stripCursor: false });
assert.equal(interactivePairwise.getAttribute('role'), 'button');

// A click with Ctrl, Cmd or Shift on a feature: the selection owner's rule. Ctrl or Cmd alone toggles the feature;
// Shift selects the range from the anchor, replacing the selection, and adds it with Ctrl or Cmd. Four features of one
// record drawn left to right on a Linear Result.
{
  const ref = (value) => ({ value });
  const ids = ['f1', 'f2', 'f3', 'f4'];
  const features = ids.map((svg_id) => ({ svg_id, fileIdx: 0, record_idx: 0, record_id: 'REC' }));
  const elements = ids.map((id, index) => Object.assign(new FakeElement(), {
    getBoundingClientRect: () => ({ left: index * 20, right: index * 20 + 10, top: 0, bottom: 10 })
  }));
  elements.forEach((element, index) => element.setAttribute('data-gbdraw-feature-id', ids[index]));
  const svg = {
    querySelectorAll: (selector) => (selector.startsWith('.')
      ? elements.filter((element) => element.classTokens().includes(selector.slice(1)))
      : elements)
  };
  const state = {
    selectedFeatureIds: ref(new Set()), selectedFeatureAnchorId: ref(''), featureSelectionStatus: ref(''),
    featureSelectionSuppressNextClick: ref(false), featureSelectionDrag: { active: false, committed: false },
    extractedFeatures: ref(features), featuresBySvgId: ref(new Map(features.map((feature) => [feature.svg_id, feature]))),
    svgContainer: ref({ querySelector: () => svg }), clickedFeature: ref(null),
    clickedPairwiseMatch: ref(null), mode: ref('linear')
  };
  const selection = createFeatureSelection({ state });
  const click = (id, modifiers) => {
    selection.selectFeatureByClick(id, { ctrlKey: false, metaKey: false, shiftKey: false, ...modifiers });
    return [[...state.selectedFeatureIds.value].sort(), state.selectedFeatureAnchorId.value];
  };
  assert.deepEqual(click('f1', { ctrlKey: true }), [['f1'], 'f1'], 'Ctrl toggles a feature on');
  assert.deepEqual(click('f1', { ctrlKey: true }), [[], ''], 'Ctrl again toggles it off');
  assert.deepEqual(click('f1', { ctrlKey: true }), [['f1'], 'f1']);
  assert.deepEqual(click('f2', { metaKey: true }), [['f1', 'f2'], 'f2'], 'Cmd toggles a feature on');
  assert.deepEqual(click('f4', { shiftKey: true }), [['f2', 'f3', 'f4'], 'f2'], 'Shift replaces with the range');
  assert.deepEqual(click('f1', { shiftKey: true, ctrlKey: true }), [['f1', 'f2', 'f3', 'f4'], 'f2'],
    'Shift+Ctrl adds the range from the anchor');
  assert.deepEqual(click('f1', { shiftKey: true }), [['f1', 'f2'], 'f2']);
  assert.deepEqual(click('f3', { shiftKey: true, metaKey: true }), [['f1', 'f2', 'f3'], 'f2'],
    'Shift+Cmd adds the range from the anchor');
}

console.log('feature selection tests passed');
