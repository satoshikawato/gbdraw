import assert from 'node:assert/strict';
import test from 'node:test';

const watchers = [];
globalThis.window = {
  Vue: {
    ref: (value) => ({ value }),
    reactive: (value) => value,
    computed: (getter) => ({ get value() { return getter(); } }),
    nextTick: async () => {},
    watch(source, callback) {
      watchers.push({ source, callback });
    }
  }
};

const { DefaultColorControl } = await import('../../gbdraw/web/js/components.js');

test('DefaultColorControl keeps picker input local until one accepted change', () => {
  const props = {
    modelValue: '#112233',
    fallback: '#cccccc',
    allowNone: true,
    ariaLabel: 'CDS feature color'
  };
  const emitted = [];
  const control = DefaultColorControl.setup(props, {
    emit(event, value) {
      emitted.push([event, value]);
      if (event === 'update:modelValue') props.modelValue = value;
    }
  });

  control.stageColor({ target: { value: '#445566' } });
  control.stageColor({ target: { value: '#778899' } });
  assert.equal(control.draftColor.value, '#778899');
  assert.deepEqual(emitted, [], 'input events must not submit semantic palette state');

  control.acceptColor({ target: { value: '#aabbcc' } });
  assert.deepEqual(emitted, [
    ['update:modelValue', '#aabbcc'],
    ['accept', '#aabbcc']
  ]);
});

test('DefaultColorControl accepts one explicit mode transition', () => {
  const props = {
    modelValue: '#112233',
    fallback: '#cccccc',
    allowNone: true,
    ariaLabel: 'CDS feature color'
  };
  const emitted = [];
  const control = DefaultColorControl.setup(props, {
    emit(event, value) {
      emitted.push([event, value]);
    }
  });

  control.updateMode({ target: { value: 'none' } });
  assert.deepEqual(emitted, [
    ['update:modelValue', 'none'],
    ['accept', 'none']
  ]);
});

test('DefaultColorControl resynchronizes its local draft after external acceptance', () => {
  const props = {
    modelValue: '#112233',
    fallback: '#cccccc',
    allowNone: true,
    ariaLabel: 'CDS feature color'
  };
  const control = DefaultColorControl.setup(props, { emit() {} });
  assert.equal(watchers.length >= 1, true);
  props.modelValue = '#abcdef';
  watchers.at(-1).callback();
  assert.equal(control.draftColor.value, '#abcdef');
});
