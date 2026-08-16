const isObjectOwner = (value) => Boolean(
  value && (typeof value === 'object' || typeof value === 'function')
);

const isPlainObject = (value) => {
  if (!value || Object.prototype.toString.call(value) !== '[object Object]') return false;
  const prototype = Object.getPrototypeOf(value);
  return prototype === null || prototype === Object.prototype;
};

const nodePosition = (node) => {
  const parent = node?.parentNode || null;
  const siblings = parent
    ? Array.from(parent.childNodes || parent.children || [])
    : [];
  return {
    parent,
    nextSibling: node?.nextSibling || null,
    siblingIndex: siblings.indexOf(node)
  };
};

const restoreNodePosition = (node, position) => {
  if (!node) return;
  if (position.parent) position.parent.insertBefore(node, position.nextSibling);
  else node.remove?.();
};

const flattenRollbackError = (errors, error) => {
  if (error instanceof AggregateError && Array.isArray(error.errors)) {
    error.errors.forEach((nested) => flattenRollbackError(errors, nested));
    return;
  }
  errors.push(error instanceof Error ? error : new Error(String(error)));
};

const attemptRestorations = (restorations, message) => {
  const errors = [];
  restorations.forEach((restore) => {
    try {
      restore();
    } catch (error) {
      flattenRollbackError(errors, error);
    }
  });
  if (errors.length > 0) throw new AggregateError(errors, message);
};

export const createDomMutationJournal = () => {
  const undo = [];
  const rememberedAttributes = new WeakMap();
  const rememberedProperties = new WeakMap();
  const rememberedText = new WeakSet();
  const rememberedStateOwners = new WeakSet();
  let lifecycle = 'open';
  let rollbackErrors = null;
  let domChangeCount = 0;
  let stateChangeCount = 0;
  let capturedStateOwnerCount = 0;

  const assertOpen = () => {
    if (lifecycle === 'rolled-back') throw new Error('DOM mutation journal was already rolled back.');
    if (lifecycle === 'committed') throw new Error('DOM mutation journal was already committed.');
  };
  const pushUndo = (callback, { kind = 'dom', countsAsChange = true } = {}) => {
    undo.push(callback);
    if (!countsAsChange) return;
    if (kind === 'state') stateChangeCount += 1;
    else domChangeCount += 1;
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
    pushUndo(() => {
      if (hadAttribute) element.setAttribute(name, previous);
      else element.removeAttribute(name);
    });
  };
  const rememberProperty = (target, key, { countsAsChange = true } = {}) => {
    let keys = rememberedProperties.get(target);
    if (!keys) {
      keys = new Set();
      rememberedProperties.set(target, keys);
    }
    if (keys.has(key)) return false;
    keys.add(key);
    const hadProperty = key in target;
    const previous = target[key];
    pushUndo(() => {
      if (hadProperty) target[key] = previous;
      else delete target[key];
    }, { kind: 'state', countsAsChange });
    return true;
  };
  const captureStateGraph = (root, { deep = true } = {}) => {
    if (!isObjectOwner(root) || rememberedStateOwners.has(root)) return false;
    const records = [];
    const visit = (owner) => {
      if (!isObjectOwner(owner) || rememberedStateOwners.has(owner)) return;
      const supported = Array.isArray(owner)
        || owner instanceof Map
        || owner instanceof Set
        || isPlainObject(owner);
      if (!supported) return;
      rememberedStateOwners.add(owner);
      capturedStateOwnerCount += 1;
      if (Array.isArray(owner)) {
        const values = Array.from(owner);
        records.push({ kind: 'array', owner, values });
        if (deep) values.forEach(visit);
        return;
      }
      if (owner instanceof Map) {
        const entries = Array.from(owner.entries());
        records.push({ kind: 'map', owner, entries });
        if (deep) entries.forEach(([key, value]) => {
          visit(key);
          visit(value);
        });
        return;
      }
      if (owner instanceof Set) {
        const values = Array.from(owner.values());
        records.push({ kind: 'set', owner, values });
        if (deep) values.forEach(visit);
        return;
      }
      const entries = Object.keys(owner).map((key) => [key, owner[key]]);
      records.push({ kind: 'object', owner, entries });
      if (deep) entries.forEach(([, value]) => visit(value));
    };
    visit(root);
    if (records.length === 0) return false;
    pushUndo(() => {
      const errors = [];
      const attempt = (restore) => {
        try {
          restore();
        } catch (error) {
          flattenRollbackError(errors, error);
        }
      };
      for (let index = records.length - 1; index >= 0; index -= 1) {
        const record = records[index];
        if (record.kind === 'array') {
          attempt(() => {
            record.owner.splice(0, record.owner.length, ...record.values);
          });
          continue;
        }
        if (record.kind === 'map') {
          attempt(() => record.owner.clear());
          record.entries.forEach(([key, value]) => {
            attempt(() => record.owner.set(key, value));
          });
          continue;
        }
        if (record.kind === 'set') {
          attempt(() => record.owner.clear());
          record.values.forEach((value) => {
            attempt(() => record.owner.add(value));
          });
          continue;
        }
        const expectedKeys = new Set(record.entries.map(([key]) => key));
        Object.keys(record.owner).forEach((key) => {
          if (!expectedKeys.has(key)) attempt(() => { delete record.owner[key]; });
        });
        record.entries.forEach(([key, value]) => {
          attempt(() => { record.owner[key] = value; });
        });
      }
      if (errors.length > 0) throw new AggregateError(errors, 'State graph rollback failed.');
    }, { kind: 'state', countsAsChange: false });
    return true;
  };

  return {
    get changed() {
      return domChangeCount + stateChangeCount > 0;
    },
    get domChangeCount() {
      return domChangeCount;
    },
    get stateChangeCount() {
      return stateChangeCount;
    },
    get capturedStateOwnerCount() {
      return capturedStateOwnerCount;
    },
    captureProperty(target, key, { deep = false } = {}) {
      assertOpen();
      if (!isObjectOwner(target)) throw new Error('State edit target does not support properties.');
      const captured = rememberProperty(target, key, { countsAsChange: false });
      if (deep) captureStateGraph(target[key], { deep: true });
      return captured;
    },
    captureState(target, options = {}) {
      assertOpen();
      return captureStateGraph(target, options);
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
        pushUndo(() => {
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
      const previousSet = new Set(previous);
      const incoming = children
        .filter((child, index) => child && !previousSet.has(child) && children.indexOf(child) === index)
        .map((child) => ({ child, position: nodePosition(child) }));
      pushUndo(() => {
        const restorations = [() => element.replaceChildren(...previous)];
        const byOriginalParent = new Map();
        incoming.forEach((record) => {
          const owner = record.position.parent;
          if (!owner) {
            restorations.push(() => restoreNodePosition(record.child, record.position));
            return;
          }
          const records = byOriginalParent.get(owner) || [];
          records.push(record);
          byOriginalParent.set(owner, records);
        });
        byOriginalParent.forEach((records) => {
          // Restore later siblings first so an earlier node's captured
          // nextSibling is already back under the original parent.
          records
            .sort((left, right) => right.position.siblingIndex - left.position.siblingIndex)
            .forEach(({ child, position }) => {
              restorations.push(() => restoreNodePosition(child, position));
            });
        });
        attemptRestorations(restorations, 'Child replacement rollback failed.');
      });
      element.replaceChildren(...children);
      return true;
    },
    appendChild(parent, child) {
      assertOpen();
      if (!parent?.appendChild || !child) throw new Error('DOM edit append target is unavailable.');
      const previous = nodePosition(child);
      pushUndo(() => restoreNodePosition(child, previous));
      parent.appendChild(child);
      return true;
    },
    insertBefore(parent, child, reference = null) {
      assertOpen();
      if (!parent?.insertBefore || !child) throw new Error('DOM edit insertion target is unavailable.');
      const previous = nodePosition(child);
      if (previous.parent === parent && previous.nextSibling === reference) return false;
      pushUndo(() => restoreNodePosition(child, previous));
      parent.insertBefore(child, reference);
      return true;
    },
    removeNode(element) {
      assertOpen();
      const previous = nodePosition(element);
      if (!previous.parent) return false;
      pushUndo(() => restoreNodePosition(element, previous));
      element.remove();
      return true;
    },
    replaceChild(parent, replacement, current) {
      assertOpen();
      if (!parent?.replaceChild || !replacement || !current) {
        throw new Error('DOM edit replacement target is unavailable.');
      }
      const currentPosition = nodePosition(current);
      const replacementPosition = nodePosition(replacement);
      pushUndo(() => {
        attemptRestorations([
          () => restoreNodePosition(current, currentPosition),
          () => restoreNodePosition(replacement, replacementPosition)
        ], 'Node replacement rollback failed.');
      });
      parent.replaceChild(replacement, current);
      return true;
    },
    setProperty(target, key, value) {
      assertOpen();
      if (!isObjectOwner(target)) throw new Error('State edit target does not support properties.');
      const hadProperty = key in target;
      if (hadProperty && Object.is(target[key], value)) return false;
      rememberProperty(target, key);
      target[key] = value;
      return true;
    },
    deleteProperty(target, key) {
      assertOpen();
      if (!target || !Object.prototype.hasOwnProperty.call(target, key)) return false;
      rememberProperty(target, key);
      delete target[key];
      return true;
    },
    rollback() {
      if (lifecycle === 'rolled-back') return rollbackErrors;
      if (lifecycle === 'committed') throw new Error('DOM mutation journal was already committed.');
      lifecycle = 'rolled-back';
      const errors = [];
      for (let index = undo.length - 1; index >= 0; index -= 1) {
        try {
          undo[index]();
        } catch (error) {
          flattenRollbackError(errors, error);
        }
      }
      undo.length = 0;
      rollbackErrors = Object.freeze(errors);
      return rollbackErrors;
    },
    commit() {
      if (lifecycle === 'rolled-back') throw new Error('DOM mutation journal was already rolled back.');
      if (lifecycle === 'committed') return false;
      lifecycle = 'committed';
      undo.length = 0;
      return true;
    }
  };
};
