const TEXT_INPUT_TYPES = new Set([
  '',
  'date',
  'datetime-local',
  'email',
  'month',
  'number',
  'password',
  'search',
  'tel',
  'text',
  'time',
  'url',
  'week'
]);

export const isIgnoredTarget = (target) =>
  Boolean(
    target?.closest?.('[data-history-ignore], [data-history-managed], [data-history-scope="transient"]')
  );

const isEditableControl = (element) => {
  if (!element) return false;
  const tag = String(element.tagName || '').toLowerCase();
  if (tag === 'textarea' || tag === 'select') return true;
  if (tag !== 'input') return false;
  return TEXT_INPUT_TYPES.has(String(element.type || '').toLowerCase());
};

const controlLabel = (element) => {
  if (!element) return 'Edit';
  const tag = String(element.tagName || '').toLowerCase();
  if (tag === 'button') return 'Change setting';
  if (tag === 'select') return 'Change setting';
  const type = String(element.type || '').toLowerCase();
  if (type === 'file') return 'Change uploaded file';
  if (type === 'checkbox' || type === 'radio') return 'Change setting';
  if (type === 'color') return 'Change color';
  return 'Edit setting';
};

export const setupHistoryInputs = ({ root, history, nextTick }) => {
  const appRoot = root || document.getElementById('app');
  if (!appRoot || !history) return () => {};

  const txByElement = new WeakMap();
  const beginByElement = new WeakMap();
  let pendingCommit = null;

  const beginForElement = (element, source = 'input-adapter') => {
    if (!element || element.disabled || isIgnoredTarget(element)) return Promise.resolve(null);
    const pendingBegin = beginByElement.get(element);
    if (pendingBegin) return pendingBegin;
    const start = (async () => {
      // Focusout and the next focusin share one DOM event turn. Finish the old
      // capture before beginning the next control's transaction.
      if (pendingCommit) await pendingCommit;
      const existing = txByElement.get(element);
      if (existing && !existing.closed) return existing;
      // The control owns its transaction; another owner settles it first (R11).
      const tx = await history.begin(controlLabel(element), { source, owner: element });
      if (tx) txByElement.set(element, tx);
      return tx;
    })();
    beginByElement.set(element, start);
    void start.finally(() => beginByElement.delete(element));
    return start;
  };

  const commitElement = async (element) => {
    const pendingBegin = beginByElement.get(element);
    if (pendingBegin) await pendingBegin;
    const tx = txByElement.get(element);
    if (!tx) return;
    if (tx.deferAdapterCommit) return;
    txByElement.delete(element);
    if (tx.closed) return;
    const commit = (async () => {
      await nextTick();
      await history.commit(tx);
    })();
    pendingCommit = commit;
    try {
      await commit;
    } finally {
      if (pendingCommit === commit) pendingCommit = null;
    }
  };

  const findControl = (eventTarget) =>
    eventTarget?.closest?.('input, textarea, select, button, [contenteditable="true"], .upload-zone') || null;

  // R11: a discrete control begins in the capture phase of its committing
  // event (checkbox and radio: change; button: click), so a pointer, its label,
  // and the keyboard record the same one step. Text-like controls begin on focus.
  const onPointerDown = (event) => {
    const target = findControl(event.target);
    if (!target || isIgnoredTarget(target) || target.closest?.('button')) return;
    if (target.classList?.contains('upload-zone')) {
      const input = target.querySelector?.('input[type="file"]');
      if (input) void beginForElement(input);
      return;
    }
    const tag = String(target.tagName || '').toLowerCase();
    const type = String(target.type || '').toLowerCase();
    if (tag === 'select' || type === 'color' || type === 'file') {
      void beginForElement(target);
    }
  };

  const onFocusIn = (event) => {
    const target = findControl(event.target);
    if (!target || isIgnoredTarget(target)) return;
    if (isEditableControl(target) || target?.isContentEditable) {
      void beginForElement(target);
    }
  };

  const onKeyDown = (event) => {
    // A shortcut chord is not an edit of the focused control.
    if (event.ctrlKey || event.metaKey || event.altKey) return;
    const target = findControl(event.target);
    if (!target || isIgnoredTarget(target)) return;
    if (isEditableControl(target) || target?.isContentEditable) {
      void beginForElement(target);
    }
  };

  const onChangeCapture = (event) => {
    const target = findControl(event.target);
    if (!target || isIgnoredTarget(target)) return;
    const type = String(target.type || '').toLowerCase();
    if (type === 'checkbox' || type === 'radio') void beginForElement(target);
  };

  const onChange = (event) => {
    const target = findControl(event.target);
    if (!target || isIgnoredTarget(target)) return;
    if (!txByElement.has(target)) {
      void beginForElement(target).then(() => commitElement(target));
      return;
    }
    const tag = String(target.tagName || '').toLowerCase();
    if (tag === 'select' || tag === 'input' || tag === 'textarea') {
      void commitElement(target);
    }
  };

  const onFocusOut = (event) => {
    const target = findControl(event.target);
    if (!target || isIgnoredTarget(target)) return;
    if (isEditableControl(target) || target?.isContentEditable) {
      void commitElement(target);
    }
  };

  const onClick = (event) => {
    const button = event.target?.closest?.('button');
    if (!button || isIgnoredTarget(button)) return;
    // Commit after the click handlers and the state they update in this task.
    setTimeout(() => { void commitElement(button); }, 0);
  };

  const onClickCapture = (event) => {
    const button = event.target?.closest?.('button');
    if (!button || isIgnoredTarget(button)) return;
    void beginForElement(button);
  };

  appRoot.addEventListener('pointerdown', onPointerDown, true);
  appRoot.addEventListener('focusin', onFocusIn, true);
  appRoot.addEventListener('keydown', onKeyDown, true);
  appRoot.addEventListener('change', onChangeCapture, true);
  appRoot.addEventListener('change', onChange, false);
  appRoot.addEventListener('focusout', onFocusOut, false);
  appRoot.addEventListener('click', onClickCapture, true);
  appRoot.addEventListener('click', onClick, false);

  return () => {
    appRoot.removeEventListener('pointerdown', onPointerDown, true);
    appRoot.removeEventListener('focusin', onFocusIn, true);
    appRoot.removeEventListener('keydown', onKeyDown, true);
    appRoot.removeEventListener('change', onChangeCapture, true);
    appRoot.removeEventListener('change', onChange, false);
    appRoot.removeEventListener('focusout', onFocusOut, false);
    appRoot.removeEventListener('click', onClickCapture, true);
    appRoot.removeEventListener('click', onClick, false);
  };
};
