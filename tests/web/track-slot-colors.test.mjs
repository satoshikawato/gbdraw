import assert from 'node:assert/strict';

const canvasContext = {
  _fillStyle: '#000000',
  get fillStyle() {
    return this._fillStyle;
  },
  set fillStyle(value) {
    const normalized = String(value || '').toLowerCase();
    if (/^#[0-9a-f]{6}$/.test(normalized)) this._fillStyle = normalized;
    if (normalized === 'tomato') this._fillStyle = '#ff6347';
  }
};
globalThis.document = {
  createElement: () => ({ getContext: () => canvasContext })
};

const {
  normalizeColorInputValue,
  resolveInheritedSkewSlotColor,
  resolveTrackSlotSkewColorValue
} = await import('../../gbdraw/web/js/app/track-slot-colors.js');

const vueRefPrototype = {};
Object.defineProperty(vueRefPrototype, 'value', {
  get() {
    return this._value;
  }
});
const vueRef = (value) => Object.assign(Object.create(vueRefPrototype), { _value: value });

const currentColors = vueRef({
  skew_high: '#112233',
  skew_low: '#445566'
});
const paletteDefinitions = vueRef({
  ocean: { skew_high: '#abcdef', skew_low: '#fedcba' },
  default: { skew_high: '#123456', skew_low: '#654321' }
});
const selectedPalette = vueRef('ocean');
const explicitSlot = {
  params: {
    positive_color: '#abc',
    negative_color: 'tomato'
  }
};
const untouchedParams = structuredClone(explicitSlot.params);

assert.equal(normalizeColorInputValue('#abc'), '#aabbcc');
assert.equal(normalizeColorInputValue('tomato'), '#ff6347');
assert.equal(normalizeColorInputValue('not-a-color'), null);
assert.equal(resolveTrackSlotSkewColorValue({
  slot: explicitSlot,
  key: 'positive_color',
  currentColors,
  paletteDefinitions,
  selectedPalette
}), '#aabbcc');
assert.equal(resolveTrackSlotSkewColorValue({
  slot: explicitSlot,
  key: 'negative_color',
  currentColors,
  paletteDefinitions,
  selectedPalette
}), '#ff6347');
assert.deepEqual(explicitSlot.params, untouchedParams);

const inheritedSlot = { params: {} };
assert.equal(resolveTrackSlotSkewColorValue({
  slot: inheritedSlot,
  key: 'positive_color',
  currentColors,
  paletteDefinitions,
  selectedPalette
}), '#112233');
assert.equal(resolveTrackSlotSkewColorValue({
  slot: inheritedSlot,
  key: 'negative_color',
  currentColors,
  paletteDefinitions,
  selectedPalette
}), '#445566');

delete explicitSlot.params.positive_color;
assert.equal(resolveTrackSlotSkewColorValue({
  slot: explicitSlot,
  key: 'positive_color',
  currentColors,
  paletteDefinitions,
  selectedPalette
}), '#112233');
assert.deepEqual(inheritedSlot.params, {});

assert.equal(resolveTrackSlotSkewColorValue({
  slot: { params: { positive_color: 'not-a-color' } },
  key: 'positive_color',
  currentColors,
  paletteDefinitions,
  selectedPalette
}), '#112233');
assert.equal(resolveInheritedSkewSlotColor({
  key: 'positive_color',
  currentColors: vueRef({}),
  paletteDefinitions,
  selectedPalette
}), '#abcdef');
assert.equal(resolveInheritedSkewSlotColor({
  key: 'negative_color',
  currentColors: {},
  paletteDefinitions: {},
  selectedPalette: 'missing'
}), '#ad72e3');
assert.equal(resolveInheritedSkewSlotColor({
  key: 'unknown',
  currentColors,
  paletteDefinitions,
  selectedPalette
}), '#777777');
