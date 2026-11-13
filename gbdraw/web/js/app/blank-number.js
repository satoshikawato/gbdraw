// @ts-check
// `v-blank-null`: a number field bound with `v-model` writes `null` when it is
// blank. Vue's `v-model.number` (and `v-model` on an `<input type="number">`)
// writes the empty string `''` for a cleared field, and also while the field
// holds text that is no number yet ("-", "1e"), so without this directive a
// blank field was held two ways, `''` and `null`, which its readers had to
// tell apart (OV-3351). This directive is the one owner of a blank number:
// every such field carries it next to its `v-model`
// (tests/web/blank-number.test.mjs), and no reader meets `''` from the app.
//
// Vue runs the `created` hooks of an element's directives in the order the
// template names them and adds `@input` listeners after them, so the listener
// below runs after `v-model`'s own and before the template's `@input`
// handler, which reads `null` already. Vue shows `null` as an empty field and
// leaves the field's text alone while it reads as the model's value, so the
// text being typed stays.

/** @typedef {(value: unknown) => void} ModelAssigner */
/** @typedef {{ props: Record<string, unknown> | null }} ModelVNode */

/** @type {WeakMap<HTMLInputElement, ModelAssigner | null>} */
const assigners = new WeakMap();

/**
 * The assignment `v-model` compiled for the field (`onUpdate:modelValue`);
 * null without a `v-model`, which tests/web/blank-number.test.mjs refuses in
 * every template.
 * @param {ModelVNode} vnode
 * @returns {ModelAssigner | null}
 */
const modelAssigner = (vnode) => {
  const update = vnode.props?.['onUpdate:modelValue'];
  /** @type {readonly unknown[]} */
  const updates = Array.isArray(update) ? update : [update];
  if (!updates.every((fn) => typeof fn === 'function')) return null;
  const fns = /** @type {readonly ModelAssigner[]} */ (updates);
  return (value) => { for (const fn of fns) fn(value); };
};

export const blankNull = {
  /**
   * @param {HTMLInputElement} field
   * @param {unknown} _binding
   * @param {ModelVNode} vnode
   */
  created(field, _binding, vnode) {
    assigners.set(field, modelAssigner(vnode));
    field.addEventListener('input', () => {
      if (field.value !== '') return;
      const assign = assigners.get(field);
      if (assign) assign(null);
    });
  },
  /**
   * The template may give the field another assignment on a render (a row of a list).
   * @param {HTMLInputElement} field
   * @param {unknown} _binding
   * @param {ModelVNode} vnode
   */
  beforeUpdate(field, _binding, vnode) {
    assigners.set(field, modelAssigner(vnode));
  },
  /** @param {HTMLInputElement} field */
  unmounted(field) {
    assigners.delete(field);
  }
};
