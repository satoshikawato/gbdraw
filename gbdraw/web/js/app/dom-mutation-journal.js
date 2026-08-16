export const createDomMutationJournal = () => {
  const undo = [];
  const rememberedAttributes = new WeakMap();
  const rememberedProperties = new WeakMap();
  const rememberedText = new WeakSet();
  let closed = false;
  const assertOpen = () => {
    if (closed) throw new Error('DOM mutation journal is already closed.');
  };
  const rememberAttribute = (element, name) => {
    let names = rememberedAttributes.get(element);
    if (!names) {
      names = new Set();
      rememberedAttributes.set(element, names);
    }
    if (names.has(name)) return;
    names.add(name);
    const previous = element?.getAttribute?.(name);
    const hadAttribute = typeof element?.hasAttribute === 'function'
      ? element.hasAttribute(name)
      : previous !== null && previous !== undefined;
    undo.push(() => {
      if (hadAttribute) element.setAttribute(name, previous);
      else element.removeAttribute(name);
    });
  };
  return {
    get changed() {
      return undo.length > 0;
    },
    setAttribute(element, name, value) {
      assertOpen();
      if (!element?.setAttribute) throw new Error('DOM edit target does not support attributes.');
      const normalized = String(value);
      if (element.getAttribute?.(name) === normalized) return false;
      rememberAttribute(element, name);
      element.setAttribute(name, normalized);
      return true;
    },
    removeAttribute(element, name) {
      assertOpen();
      if (!element?.removeAttribute || element.getAttribute?.(name) === null) return false;
      rememberAttribute(element, name);
      element.removeAttribute(name);
      return true;
    },
    setTextContent(element, value) {
      assertOpen();
      if (!element) throw new Error('DOM edit text target is unavailable.');
      const normalized = String(value ?? '');
      if (String(element.textContent ?? '') === normalized) return false;
      if (!rememberedText.has(element)) {
        rememberedText.add(element);
        const previous = element.textContent;
        undo.push(() => {
          element.textContent = previous;
        });
      }
      element.textContent = normalized;
      return true;
    },
    replaceChildren(element, ...children) {
      assertOpen();
      if (!element?.replaceChildren) throw new Error('DOM edit child target is unavailable.');
      const previous = Array.from(element.childNodes || []);
      undo.push(() => element.replaceChildren(...previous));
      element.replaceChildren(...children);
      return true;
    },
    appendChild(parent, child) {
      assertOpen();
      if (!parent?.appendChild || !child) throw new Error('DOM edit append target is unavailable.');
      const previousParent = child.parentNode || null;
      const previousNextSibling = child.nextSibling || null;
      undo.push(() => {
        if (previousParent) previousParent.insertBefore(child, previousNextSibling);
        else child.remove?.();
      });
      parent.appendChild(child);
      return true;
    },
    insertBefore(parent, child, reference = null) {
      assertOpen();
      if (!parent?.insertBefore || !child) throw new Error('DOM edit insertion target is unavailable.');
      const previousParent = child.parentNode || null;
      const previousNextSibling = child.nextSibling || null;
      if (previousParent === parent && previousNextSibling === reference) return false;
      undo.push(() => {
        if (previousParent) previousParent.insertBefore(child, previousNextSibling);
        else child.remove?.();
      });
      parent.insertBefore(child, reference);
      return true;
    },
    removeNode(element) {
      assertOpen();
      const parent = element?.parentNode || null;
      if (!parent) return false;
      const nextSibling = element.nextSibling || null;
      undo.push(() => parent.insertBefore(element, nextSibling));
      element.remove();
      return true;
    },
    replaceChild(parent, replacement, current) {
      assertOpen();
      if (!parent?.replaceChild || !replacement || !current) {
        throw new Error('DOM edit replacement target is unavailable.');
      }
      undo.push(() => {
        if (replacement.parentNode === parent) parent.replaceChild(current, replacement);
      });
      parent.replaceChild(replacement, current);
      return true;
    },
    setProperty(target, key, value) {
      assertOpen();
      if (!target || (typeof target !== 'object' && typeof target !== 'function')) {
        throw new Error('State edit target does not support properties.');
      }
      const hadProperty = key in target;
      if (hadProperty && Object.is(target[key], value)) return false;
      let keys = rememberedProperties.get(target);
      if (!keys) {
        keys = new Set();
        rememberedProperties.set(target, keys);
      }
      if (!keys.has(key)) {
        keys.add(key);
        const previous = target[key];
        undo.push(() => {
          if (hadProperty) target[key] = previous;
          else delete target[key];
        });
      }
      target[key] = value;
      return true;
    },
    deleteProperty(target, key) {
      assertOpen();
      if (!target || !Object.prototype.hasOwnProperty.call(target, key)) return false;
      let keys = rememberedProperties.get(target);
      if (!keys) {
        keys = new Set();
        rememberedProperties.set(target, keys);
      }
      if (!keys.has(key)) {
        keys.add(key);
        const previous = target[key];
        undo.push(() => {
          target[key] = previous;
        });
      }
      delete target[key];
      return true;
    },
    rollback() {
      if (closed) return;
      closed = true;
      for (let index = undo.length - 1; index >= 0; index -= 1) undo[index]();
      undo.length = 0;
    },
    commit() {
      if (closed) return;
      closed = true;
      undo.length = 0;
    }
  };
};
